// ------------------------------------------------------------------
//   sam_gc_load_sa.c
//   Copyright (C) 2022-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "sections.h"
#include "codec.h"
#include "compressor.h"
#include "dispatcher.h"
#include "profiler.h"
#include "piz.h"
#include "reconstruct.h"
#include "zfile.h"
#include "qname.h"
#include "bit_array.h"
#include "writer.h"
#include "htscodecs/rANS_static4x16.h"
#include "libdeflate/libdeflate.h"

typedef struct {
    VBIType vblock_i;
    SAGroup first_grp_i; 
    uint32_t num_grps, num_alns;
    uint32_t seq_len;       // in bases
    uint32_t comp_qual_len; // size as of compressed in-memory in SA Groups. 
    uint32_t comp_cigars_len;
    uint32_t qname_len;
    uint64_t seq_start, qual_start, qname_start, first_aln_i;
} PlsgVbInfo;

#define plsg_info z_file->vb_info[1]
static uint32_t next_plsg_i = 0; // iterator for dispatching compute threads
static VBIType num_prim_vbs_loaded = 0; // number of PRIM VBs whose loading is complete
static Mutex num_prim_vbs_loaded_mutex = {};

static Mutex copy_qual_mutex = {}, copy_cigars_mutex = {};

#define vb_qual_buf   vb->z_data  // used for QUAL data being in-memory compressed in the loader compute thread
#define vb_cigars_buf vb->codec_bufs[1] // similar for cigar data

static Flags save_flag; 

//-------------------------------------------------------------------------
// PIZ: SA Group loader
// Method: decompress directly into memory, parallelizing with dispatcher and 3 threads per PRIM VB: SEQ, QUAL and Alignments
//-------------------------------------------------------------------------

ShowAln sam_show_sa_one_aln (const SAGroupType *g, const SAAlnType *a)
{
    ShowAln s;
    
    if (command==PIZ) {
        char cigar_info[64];
        if (a->cigar.piz.is_word) 
            sprintf (cigar_info, "word=%u", (WordIndex)a->cigar.piz.index);
        else
            sprintf (cigar_info, "in=%"PRIu64" l=%u comp=%u", (uint64_t)a->cigar.piz.index, ALN_CIGAR_LEN(a), (int)a->cigar.piz.comp_len);

        char rname_str[64]; // should be big enough...
        sprintf (rname_str, "\"%.48s\"(%u)", ctx_get_snip_by_word_index0 (ZCTX(OPTION_SA_RNAME), a->rname), a->rname);

        sprintf (s.s, "aln_i=%u.%u: sa_rname=%-13s pos=%-10u mapq=%-3u strand=%c nm=%-3u cigar[%u]=\"%s\"(%s)",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), rname_str,
                a->pos, a->mapq, "+-"[a->revcomp], a->nm,
                SA_CIGAR_DISPLAY_LEN, sam_piz_display_aln_cigar (a), cigar_info);
    }
    else 
        sprintf (s.s, "aln_i=%u.%u: cigar_sig=%-12.12s rname=%-4u pos=%-10u mapq=%-3u strand=%c  nm=%-3u",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), cigar_display_signature(a->cigar.signature).s, 
                a->rname, a->pos, a->mapq, "+-"[a->revcomp], a->nm);
    
    return s;
}

void sam_show_sa_one_grp (SAGroup grp_i)
{
    const SAGroupType *g = B(SAGroupType, z_file->sa_groups, grp_i);

    iprintf ("grp_i=%u: qname(i=%"PRIu64",l=%u,hash=%u)=%.*s seq=(i=%"PRIu64",l=%u) qual=(i=%"PRIu64",comp_l=%u)[%u]=\"%s\" strand=%c mul/fst/lst=%u,%u,%u num_alns=%u first_aln=%"PRIu64"%s\n",
             grp_i, (uint64_t)g->qname, g->qname_len, 
             QNAME_HASH (GRP_QNAME(g), g->qname_len, g->is_last), // the hash is part of the index in ZIP, and not part of the SAGroup struct. Its provided here for its usefulness in debugging
             g->qname_len, GRP_QNAME(g) , g->seq, g->seq_len, 
             g->qual, g->qual_comp_len, SA_QUAL_DISPLAY_LEN, sam_display_qual_from_SA_Group (g),  
             "+-"[g->revcomp], g->multi_segments, g->is_first, g->is_last,
             g->num_alns, (uint64_t)g->first_aln_i, g->first_grp_in_vb ? " FIRST_GRP_IN_VB" : "");

    for (uint64_t aln_i=g->first_aln_i; aln_i < g->first_aln_i + g->num_alns; aln_i++) 
        iprintf ("    %s\n", sam_show_sa_one_aln (g, B(SAAlnType, z_file->sa_alns, aln_i)).s);
}

void sam_show_sa (void)
{
    if (flag.show_sa >= 1) {
        ASSINP (flag.show_sa <= z_file->sa_groups.len, "--show-sa argument, for this file, should be between 0 and %"PRIu64, 
                z_file->sa_groups.len-1);
        
        sam_show_sa_one_grp (flag.show_sa - 1);
    }

    else {
        iprintf ("NUM_SA_GROUPS=%"PRIu64" NUM_SA_ALIGNMENTS=%"PRIu64"\n", z_file->sa_groups.len, z_file->sa_alns.len);

        for (SAGroup grp_i=0; grp_i < z_file->sa_groups.len; grp_i++) 
            sam_show_sa_one_grp (grp_i);
    }
}

// QNAME: reconstruct directly into z_file->sa_qname by overlaying txt_data
static inline void sam_load_groups_add_qname (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *vb_grps) 
{
    buf_set_overlayable (&z_file->sa_qnames);
    buf_overlay_partial (vb, &vb->txt_data, &z_file->sa_qnames, plsg->qname_start, "txt_data");
    
    vb->txt_data.len = 0;
    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {
        vb->line_i = grp_i; // needed for buddying to word
        vb->buddy_line_i = NO_LINE;

        vb_grps[grp_i].qname = plsg->qname_start + vb->txt_data.len;
        reconstruct_from_ctx (vb, SAM_QNAME, 0, true); // reconstructs into vb->txt_data, sets vb->buddy_line_i if SNIP_COPY_BUDDY
        vb_grps[grp_i].qname_len = vb->txt_data.len - (vb_grps[grp_i].qname - plsg->qname_start);

        // PRIMARY alignment with the QNAME has a buddy (but no info here about DEPN alignments!)
        // When reconstructing a primary VB, sam_piz_set_sa_grp will set the buddy if vb_grps[grp_i].prim_set_buddy is true
        vb_grps[grp_i].prim_set_buddy = piz_has_buddy; 
    }

    buf_free (vb->txt_data);
}

// FLAGS: get multi_segments, is_first, is_last from SAM_FLAGS. 
static inline void sam_load_groups_add_flags (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *vb_grps) 
{
    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {

        // get SAM_FLAG snip, which is a SPECIAL with the parameter being the flags without revcomp (see sam_seg_FLAG)
        STR0(snip);
        LOAD_SNIP(SAM_FLAG); 

        ASSERT0 (snip[0]==SNIP_SPECIAL && snip[1]==SAM_SPECIAL_SAGROUP, "expecting FLAG to be SAM_SPECIAL_SAGROUP");

        SamFlags sam_flags = { .value = atoi (&snip[2]) };
        vb_grps[grp_i].multi_segments = sam_flags.bits.multi_segments;
        vb_grps[grp_i].is_first       = sam_flags.bits.is_first;
        vb_grps[grp_i].is_last        = sam_flags.bits.is_last;
    }
}

// Grp loader compute thread: copy to z_file->sa_qual
static inline void sam_load_groups_move_comp_to_zfile (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *vb_grps, SAAlnType *vb_alns) 
{
    plsg->comp_qual_len   = vb_qual_buf.len; // update to actual (might not be the same as received from ZIP, if codec has changed)
    plsg->comp_cigars_len = vb_cigars_buf.len;

    // copy buffers in arbitrary order, based on mutex availability
    bool qual_done=false, cigars_done=false;
    uint64_t start_qual=0, start_cigars=0;

    while (!qual_done || !cigars_done) {

        bool achieved_something = false;
        
        // QUAL. note: we lock as we might realloc (therefore other threads can't access this Buffer concurrently)
        if (!qual_done && mutex_trylock (copy_qual_mutex)) {
            start_qual = z_file->sa_qual.len;
            buf_add_buf (evb, &z_file->sa_qual, &vb_qual_buf, uint8_t, NULL); // might grow the buffer
            qual_done = achieved_something = true;
            mutex_unlock (copy_qual_mutex);
        }

        // CIGARs
        if (!cigars_done && mutex_trylock (copy_cigars_mutex)) {
            start_cigars = z_file->sa_cigars.len;
            buf_add_buf (evb, &z_file->sa_cigars, &vb_cigars_buf, uint8_t, NULL); 
            cigars_done = achieved_something = true;
            mutex_unlock (copy_cigars_mutex);
        }

        if (!achieved_something) usleep (1000); // 1 ms
    }

    buf_free (vb_qual_buf);
    buf_free (vb_cigars_buf);

    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) 
        vb_grps[grp_i].qual += start_qual;

    for (SAGroup aln_i=0; aln_i < plsg->num_alns ; aln_i++) 
        if (!vb_alns[aln_i].cigar.piz.is_word)
            vb_alns[aln_i].cigar.piz.index += start_cigars; // update from index into vb_cigars_buf to z_file->sa_cigars 
}

static inline void sam_load_groups_add_seq (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *g, uint64_t z_seq_start)
{
    if (!CTX(SAM_SQBITMAP)->is_loaded) return; // sequence is skipped

    // reconstruct SEQ to vb->textual_seq, needed by sam_piz_prim_add_QUAL, and overlay txt_data on it
    buf_alloc (vb, &vb->textual_seq, 0, vb->seq_len, char, 0, "textual_seq");
    buf_set_overlayable (&vb->textual_seq);
    buf_overlay (vb, &vb->txt_data, &vb->textual_seq, "txt_data");

    reconstruct_from_ctx (vb, SAM_SQBITMAP, 0, true);
    vb->textual_seq.len = vb->txt_data.len;
    buf_free (vb->txt_data); // un-overlay

    ASSERT (vb->textual_seq.len == vb->seq_len, "Expecting textual_seq.len=%"PRIu64" == seq_len=%u", vb->textual_seq.len, vb->seq_len);

    // pack SEQ data into z_file->sa_seq
    BitArray *z_sa_seq = buf_get_bitarray (&z_file->sa_seq);
    sam_sa_native_to_acgt (vb, z_sa_seq, z_seq_start * 2, B1STc(vb->textual_seq), vb->seq_len, false, false, false); 
}

// compress QUAL data for in-memory storage, using fast rans codec.
static inline void sam_load_groups_add_qual (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *g)
{
    if (!CTX(SAM_QUAL)->is_loaded) return; // qual is skipped

    // to reconstruct into codec_bufs[0], we overlay txt_data on it, as the reconstruct machinery reconstructs to txt_data 
    buf_alloc (vb, &vb->codec_bufs[0], vb->seq_len, 0, char, 1, "codec_bufs[0]");
    buf_set_overlayable (&vb->codec_bufs[0]);
    buf_overlay (vb, &vb->txt_data, &vb->codec_bufs[0], "txt_data");

    // Reconstruct QUAL of one line into txt_data which is overlaid on codec_bufs[0]
    reconstruct_from_ctx (vb, SAM_QUAL, 0, true); // reconstructs into vb->txt_data
    vb->codec_bufs[0].len = vb->txt_data.len;
    g->no_qual = vb->qual_missing;
    buf_free (vb->txt_data); // un-overlay

    if (g->no_qual) return; // no quality data for this group - we're done

    ASSERT (vb->codec_bufs[0].len == vb->seq_len, "Expecting vb->txt_data.len=%"PRIu64" == vb->seq_len=%u", vb->txt_data.len, vb->seq_len);

    // Compress QUAL, or keep it as is - which ever is better
    uint32_t comp_len = rans_compress_bound_4x16 (g->seq_len, X_NOSZ); // maximum 
    buf_alloc (vb, &vb_qual_buf, comp_len, 0, uint8_t, 1, NULL); // likely already allocated in sam_load_groups_add_grps

    // uncompressed qual is is codec_bufs[0]; append compressed qual to vb_qual_buf
    ASSERT (rans_compress_to_4x16 (VB, B1ST8 (vb->codec_bufs[0]), g->seq_len, BAFT(uint8_t, vb_qual_buf), &comp_len, X_NOSZ) && comp_len,
            "Failed to compress PRIM qual of vb=%u grp_i=%u qual_len=%u", plsg->vblock_i, ZGRP_I(g), g->seq_len);

    g->qual = vb_qual_buf.len; // relative to the VB. we will update later to be relative to z_file->sa_qual.

    // add to buffer only if compression actually compresses
    if (comp_len < g->seq_len) { 
        g->qual_comp_len = comp_len;  
        vb_qual_buf.len += g->qual_comp_len;
    }
    else {
        buf_add_buf (vb, &vb_qual_buf, &vb->codec_bufs[0], char, NULL); // add uncompressed data instead of compressed
        g->qual_comp_len = 0; // we will use qual_comp_len==0 to mean "not compressed"
    }

    vb->codec_bufs[0].len = 0; // faster than buf_free
}

// populates the cigar data of one alignment
static void sam_load_groups_add_aln_cigar (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *g, SAAlnType *a, bool is_first_aln,
                                           bool is_all_the_same_LOOKUP, bool is_all_the_same_SQUANK,
                                           pSTRp(out_cigar)) // optional out
{
    ContextP ctx = CTX (OPTION_SA_CIGAR);

    STR(snip);
    WordIndex word_index = WORD_INDEX_NONE;
    
    // case: word_list.len==0 if all snips are in local, so dictionary was an all-the-same SNIP_LOOKUP, therefore dropped by ctx_drop_all_the_same, 
    if (ctx->word_list.len)
        word_index = ctx_get_next_snip (VB, ctx, false, pSTRa(snip)); 

    bool is_lookup = is_all_the_same_LOOKUP || *snip == SNIP_LOOKUP;
    bool is_squank = !is_lookup && (is_all_the_same_SQUANK || (snip_len > 2 && snip[0] == SNIP_SPECIAL && snip[1] == SAM_SPECIAL_SQUANK));

    // case: the CIGAR is in local of this vb we need to copy it (compressed) as the "SA Group loader" VB will be soon released. 
    if (is_lookup || is_squank) {

        if (is_lookup)
            LOAD_SNIP_FROM_LOCAL (ctx); // updates snip, snip_len

        else { // squank - temporarily reconstruct into txt_data - just for compressing. note: prim cigar is never segged as squank
            sam_piz_special_SQUANK (VB, ctx, &snip[2], snip_len-2, NULL/*reconstruct to vb->scratch*/, true);
            snip     = vb->scratch.data;
            snip_len = vb->scratch.len;
        }

        // note: we put our compressed cigars in vb_cigars_buf, to be copied to z_file->sa_cigars later.
        // this is because we can't assume our compressed length is the same as it was in ZIP, as codecs may change between genozip versions.
        uint32_t this_cigar_mem = rans_compress_bound_4x16 (snip_len, X_NOSZ); // initial allocation, we may grow it later if needed
        uint32_t all_cigar_mem_this_vb = is_first_aln ? (plsg->comp_cigars_len + rans_compress_bound_4x16 (snip_len, X_NOSZ)) : 0; // case first allocation: estimate the total memory for all cigars on alignments of this group, to save allocations
        buf_alloc (vb, &vb_cigars_buf, this_cigar_mem, all_cigar_mem_this_vb, char, 0, "scratch"); 

        uint32_t comp_len = snip_len; // initialize pessimistically- no compression

        if (snip_len > 30) { // no point even trying with shortish cigars
            comp_len = rans_compress_bound_4x16 (snip_len, X_NOSZ); // maximum 
            buf_alloc (vb, &vb_cigars_buf, comp_len, 0, uint8_t, 1, NULL);

            // append compressed cigars to vb_cigars_buf
            ASSERT (rans_compress_to_4x16 (VB, (uint8_t*)STRa(snip), BAFT(uint8_t, vb_cigars_buf), &comp_len, X_NOSZ) && comp_len,
                    "Failed to compress cigar of vb=%u grp_i=%u aln=%"PRId64" cigar_len=%u cigar=\"%.*s\"", 
                    plsg->vblock_i, ZGRP_I(g), ZALN_I(a), snip_len, STRf(snip));
        }
                
        // case: compression doesn't compress - abandon compression and store uncompressed - set comp_len to 0.
        if (comp_len >= snip_len) {
            memcpy (BAFT(uint8_t, vb_cigars_buf), snip, snip_len);
            comp_len = 0;
        }

        a->cigar.piz.is_word  = false; // cigar is in z_file->sa_cigars (long, possibly compressed)
        a->cigar.piz.index    = vb_cigars_buf.len;
        a->cigar.piz.comp_len = comp_len; // 0 if not compressed
        a->cigar.piz.len_lo   = snip_len & MAXB(ALN_CIGAR_LEN_BITS_LO);
        a->cigar.piz.len_hi   = (snip_len >> ALN_CIGAR_LEN_BITS_LO);

        vb_cigars_buf.len += comp_len ? comp_len : snip_len;
    }

    // case: the cigar is in dict - we therefore don't store - reconstruction will copy it from the dictionary
    else {
        a->cigar.piz.is_word = true; // cigar is in ZCTX(OPTION_SA_CIGAR).dict (short cigar)
        a->cigar.piz.index   = word_index;

        ASSERT (a->cigar.piz.index < ctx->word_list.len, "word_index=%d out of range [0,%d] snip_len=%u snip[0]=%u snip=\"%.*s", 
                (uint32_t)a->cigar.piz.index, ctx->word_list.len32-1, STRf(snip)[0], STRf(snip));
    }

    if (out_cigar) STRset (*out_cigar, snip);
    if (is_squank) buf_free (vb->scratch);
}

static void sam_load_groups_add_grp_cigars (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *g, SAAlnType *a, 
                                            bool is_all_the_same_LOOKUP, bool is_all_the_same_SQUANK,
                                            pSTRp(prim_cigar)) // out
{
    // primary alignment CIGAR 
    sam_load_groups_add_aln_cigar (vb, plsg, g, a, ZGRP_I(g)==0, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, STRa(prim_cigar));

    // calculate seq_len from primary CIGAR, and copy cigar to vb->textual_cigar
    sam_cigar_analyze (vb, STRa(*prim_cigar), false, &vb->seq_len);

    // populate SAAlnType.cigar with the CIGAR word index, for the non-primary alignments of this group
    for (uint8_t aln_i=1; aln_i < g->num_alns; aln_i++) 
        sam_load_groups_add_aln_cigar (vb, plsg, g, a + aln_i, false, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, 0, 0);
}

// SEQ - ACGT-pack uncompressed sequence directly to z_file->sa_seq
// Alns - populate CIGAR the cigar field of the alignments
// Grp  - populate seq, seq_len, num_alns
static inline void sam_load_groups_add_grps (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *vb_grps) 
{
    ASSERTNOTINUSE (vb_qual_buf);
    ASSERTNOTINUSE (vb_cigars_buf);

    // note: we put our compressed QUAL in vb_qual_buf, to be copied to z_file->sa_qual later.
    // this is because we can't assume our compressed length is the same as it was in ZIP, as codecs may change between genozip versions.
    uint32_t qual_comp_len_this_vb = plsg->comp_qual_len + rans_compress_bound_4x16 (vb_grps[0].seq_len, X_NOSZ); // initial allocation, we may grow it later if needed
    buf_alloc (vb, &vb_qual_buf, qual_comp_len_this_vb, 0, char, 0, "z_data"); 

    // get seq and seq_len from CIGAR
    uint32_t total_seq_len = 0;

    uint32_t sa_rname_len = CTX(OPTION_SA_RNAME)->word_list.len32;

    ContextP cigar_ctx = CTX(OPTION_SA_CIGAR);
    bool is_cigar_all_the_same_LOOKUP = !cigar_ctx->word_list.len && cigar_ctx->dict.len && *B1STc(cigar_ctx->dict)==SNIP_LOOKUP;
    bool is_cigar_all_the_same_SQUANK = !cigar_ctx->word_list.len && cigar_ctx->dict.len >= 3 && 
                                        *B1STc(cigar_ctx->dict)==SNIP_SPECIAL && *Bc(cigar_ctx->dict, 1)==SAM_SPECIAL_SQUANK;

    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {
        SAGroupType *g = &vb_grps[grp_i];

        ASSERT (g->first_aln_i >= plsg->first_aln_i && g->first_aln_i < z_file->sa_alns.len,
                "grp_i=%u z_grp_i=%u has first_aln_i=%"PRIu64" which is out of range [0,%"PRIu64"]", 
                grp_i, ZGRP_I(g), (uint64_t)g->first_aln_i, z_file->sa_alns.len-1);

        vb->line_i = grp_i;
        sam_reset_line (VB);

        SAAlnType *prim_aln = B(SAAlnType, z_file->sa_alns, g->first_aln_i); // first Aln of group

        // add all alignment cigars
        STR(prim_aln_cigar);
        if (CTX(OPTION_SA_CIGAR)->is_loaded) // not skipped
            sam_load_groups_add_grp_cigars (vb, plsg, g, prim_aln, is_cigar_all_the_same_LOOKUP, is_cigar_all_the_same_SQUANK, pSTRa(prim_aln_cigar));
        
        g->seq     = plsg->seq_start + total_seq_len; // in bases
        g->seq_len = vb->seq_len;   // in bases
        
        // populate vb with alignment details
        vb->last_int(SAM_POS) = prim_aln->pos;
        CTX(SAM_FLAG)->last_value.i = prim_aln->revcomp ? SAM_FLAG_REV_COMP : 0; // needed by codec_longr_reconstruct

        ASSERT (prim_aln->rname >= 0 && prim_aln->rname < sa_rname_len, 
                "rname=%d out of range: OPTION_SA_RNAME.word_len.len=%u. primary alignment of grp_i=%u: QNAME=(%"PRIu64",%u)=%.*s, CIGAR[12]=%*.*s",
                 prim_aln->rname, sa_rname_len, grp_i, (uint64_t)g->qname, (int)g->qname_len, (int)g->qname_len, Bc(z_file->sa_qnames, g->qname), MIN_(prim_aln_cigar_len,12),MIN_(prim_aln_cigar_len,12), prim_aln_cigar);

        ctx_get_snip_by_word_index (CTX(OPTION_SA_RNAME), prim_aln->rname, vb->chrom_name);
        
        vb->chrom_node_index = ctx_get_word_index_by_snip (VB, CTX(SAM_RNAME), STRa(vb->chrom_name)); // convert OPTION_SA_RNAME word_index to RNAME word_index
        ASSERT (vb->chrom_node_index != WORD_INDEX_NONE, "Cannot find rname=%u in context RNAME", prim_aln->rname);
        
        // add SEQ data to the group, packing it in-memory
        if (CTX(SAM_SQBITMAP)->is_loaded)
            sam_load_groups_add_seq (vb, plsg, g, plsg->seq_start + total_seq_len);

        // add QUAL data to the group, possibly compressing it in-memory
        if (CTX(SAM_QUAL)->is_loaded)
            sam_load_groups_add_qual (vb, plsg, g);

        total_seq_len += vb->seq_len;
        
        vb->textual_seq.len = vb->codec_bufs[0].len = 0;
    }

    buf_free (vb->textual_seq);
    buf_free (vb->codec_bufs[0]);
}

// Alns - populate RNAME, POS, MAPQ, STRAND and NM (CIGAR is added sam_load_groups_add_grps)
static inline void sam_load_groups_add_alns (VBlockSAMP vb, PlsgVbInfo *plsg, SAGroupType *vb_grps, SAAlnType *vb_alns) 
{
    ContextP sa_ctx    = LOADED_CTX(OPTION_SA_Z);
    ContextP rname_ctx = CTX(OPTION_SA_RNAME); // possibly not loaded

    // sanity check
    ARRAY (uint8_t, num_alns, sa_ctx->local); // this contains the number of alignments for each SA Group (inc. the primary alignment)
    ASSERT (num_alns_len == plsg->num_grps, "Mismatch in number of PRIM lines for prim vb=%u: expecting OPTION_SA_Z.local.len=%u == num_lines=%u OPTION_SA_Z.is_loaded=%s",
            plsg->vblock_i, (unsigned)num_alns_len, plsg->num_grps, TF(sa_ctx->is_loaded));

    for (uint32_t grp_i=0; grp_i < num_alns_len; grp_i++) {
        vb_grps[grp_i].first_aln_i = grp_i ? vb_grps[grp_i-1].first_aln_i + num_alns[grp_i-1] : plsg->first_aln_i;
        vb_grps[grp_i].num_alns    = num_alns[grp_i];
    }

    int32_t vb_num_alns = (vb_grps[plsg->num_grps-1].first_aln_i - plsg->first_aln_i) + num_alns[plsg->num_grps-1];
    ASSERT (vb_num_alns == plsg->num_alns, "Alignment count mismatch for vb=%u: expecting vb_num_alns=%d == plsg->num_alns=%u",
            plsg->vblock_i, vb_num_alns, plsg->num_alns);

    uint32_t sa_rname_len = rname_ctx->word_list.len32;

    for (uint32_t aln_i=0; aln_i < plsg->num_alns; aln_i++) {
        
        SAAlnType *a = &vb_alns[aln_i];

        // RNAME
        if (rname_ctx->is_loaded) { // RNAME not skipped
            reconstruct_from_ctx (vb, OPTION_SA_RNAME, 0, false);  // don't reconstruct, only set last_value to word_index
            a->rname = rname_ctx->last_value.i;         // WordIndex into OPTION_SA_RNAME, NOT RNAME!

            ASSERT (a->rname >= 0 && a->rname < sa_rname_len, 
                    "While loading aln_i=%u: rname=%u out of range: OPTION_SA_RNAME.word_len.len=%u", aln_i, a->rname, sa_rname_len);
        }

        // POS
        if (CTX(OPTION_SA_POS)->is_loaded) {
            reconstruct_from_ctx (vb, OPTION_SA_POS, 0, false);    // don't reconstruct, only set last_value to pos
            a->pos = CTX(OPTION_SA_POS)->last_value.i;
        }

        // STRAND - nodes initialized in sam_seg_0X_initialize to 0="-" 1="+"
        if (CTX(OPTION_SA_STRAND)->is_loaded) {
            reconstruct_from_ctx (vb, OPTION_SA_STRAND, 0, false); // don't reconstruct, only set last_value to word_index
            a->revcomp = !CTX(OPTION_SA_STRAND)->last_value.i;
        }
        
        // MAPQ
        if (CTX(OPTION_SA_MAPQ)->is_loaded) { // MAPQ is not skipped
            reconstruct_from_ctx (vb, OPTION_SA_MAPQ, 0, false);   // don't reconstruct, only set last_value to pos
            a->mapq = CTX(OPTION_SA_MAPQ)->last_value.i;
        }

        // NM
        if (CTX(OPTION_SA_NM)->is_loaded) { // NM is not skipped
            reconstruct_from_ctx (vb, OPTION_SA_NM, 0, false);     // don't reconstruct, only set last_value to pos
            a->nm = CTX(OPTION_SA_NM)->last_value.i;
        }
    }

    // group revcomp is the revcomp of its primary alignment
    for (uint32_t grp_i=0; grp_i < num_alns_len; grp_i++) 
        vb_grps[grp_i].revcomp = vb_alns[vb_grps[grp_i].first_aln_i - plsg->first_aln_i].revcomp;
}

// entry point of compute thread of SA Groups loading
static void sam_load_groups_add_one_prim_vb (VBlockP vb_)
{
    START_TIMER;
    VBlockSAMP vb = (VBlockSAMP)vb_;

    piz_uncompress_all_ctxs (VB);

    buf_free (vb->z_data); // we will use z_data for vb_qual_buf to avoid allocating a separate buffer

    PlsgVbInfo *plsg     = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    SAGroupType *vb_grps = B(SAGroupType, z_file->sa_groups, plsg->first_grp_i);
    SAAlnType *vb_alns   = B(SAAlnType, z_file->sa_alns, plsg->first_aln_i);

    vb_grps[0].first_grp_in_vb = true;

    sam_load_groups_add_alns (vb, plsg, vb_grps, vb_alns);
    if (CTX(SAM_QNAME)->is_loaded) sam_load_groups_add_qname (vb, plsg, vb_grps);
    sam_load_groups_add_flags (vb, plsg, vb_grps);
    sam_load_groups_add_grps (vb, plsg, vb_grps);
    sam_load_groups_move_comp_to_zfile (vb, plsg, vb_grps, vb_alns); // last, as it might block
    
    mutex_lock (num_prim_vbs_loaded_mutex);
    num_prim_vbs_loaded++;
    mutex_unlock (num_prim_vbs_loaded_mutex);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.

    COPY_TIMER_VB (evb, sam_load_groups_add_one_prim_vb);                          
}

// PIZ main thread - dispatches compute compute thread do read SA Groups from one PRIM VB. returns true if dispatched
bool sam_piz_dispatch_one_load_SA_Groups_vb (Dispatcher dispatcher)
{
    if (!VER(14) || flag.genocat_no_reconstruct) return false; // no need to load SA groups

    // check if we're done
    if (next_plsg_i == plsg_info.len) goto done_loading;

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, next_plsg_i++);

    VBlockP vb = dispatcher_generate_next_vb (dispatcher, plsg->vblock_i, SAM_COMP_PRIM);
    vb->comp_i              = SAM_COMP_PRIM;
    vb->preprocessing       = true; // we're done dispatching preprocessing VBs, however, preprocessing VB compute threads may still be running.
    vb->show_containers     = false; 

    piz_read_one_vb (vb, false);
    dispatcher_compute (dispatcher, sam_load_groups_add_one_prim_vb);

    return true;

done_loading:
    // restore globals
    flag = save_flag;  // also resets flag.preprocessing
    dispatcher_set_task_name (dispatcher, PIZ_TASK_NAME);

    return false;
}

// PIZ main thread: after joining a preprocessing VB
void sam_piz_after_preproc (VBlockP vb)
{
    z_file->num_preproc_vbs_joined++;

    // print SA after all groups are loaded
    if (flag.show_sa && sections_get_num_vbs (SAM_COMP_PRIM) == z_file->num_preproc_vbs_joined) {
        sam_show_sa();
        if (exe_type == EXE_GENOCAT) exit(0);
    }

    if (flag.show_vblocks) 
        iprintf ("LOADED_SA(id=%d) vb=%s\n", vb->id, VB_NAME);
}

// PIZ main thread: a callback of piz_after_global_area 
void sam_piz_load_SA_Groups (void)
{
    if (sections_get_num_comps() == 1 || // no PRIM/DEPN in this z_file
        flag.genocat_no_reconstruct) return; 

    uint32_t num_prim_vbs = sections_get_num_vbs (SAM_COMP_PRIM);
    if (!num_prim_vbs) return; // no actual PRIM lines (only DEPN)

    ARRAY_alloc (PlsgVbInfo, plsg, num_prim_vbs, false, plsg_info, evb, "z_file->plsg");

    uint64_t total_seq_len=0, total_qual_len=0, total_cigars_len=0, total_qname_len=0, total_alns=0; // total across all PRIM VBs of the file
    SAGroup total_grps=0;

    // surveys all PRIM VBs in the order the appear in the file (note: they needn't be consecutive in the file)
    Section vb_header_sec = NULL;
    for (uint32_t i=0; i < num_prim_vbs; i++) {
        
        sections_get_next_vb_of_comp_sec (SAM_COMP_PRIM, &vb_header_sec);

        ASSERT0 (vb_header_sec->st==SEC_VB_HEADER && vb_header_sec->comp_i==SAM_COMP_PRIM, "expecting a PRIM VB Header");

        SectionHeaderVbHeader header = zfile_read_section_header (evb, vb_header_sec->offset, vb_header_sec->vblock_i, SEC_VB_HEADER).vb_header;

        plsg[i] = (PlsgVbInfo){ .vblock_i        = vb_header_sec->vblock_i,
                                .first_grp_i     = BGEN32 (header.sam_prim_first_grp_i),
                                .first_aln_i     = total_alns,
                                .seq_start       = total_seq_len, // in bases
                                .qual_start      = total_qual_len,
                                .qname_start     = total_qname_len,
                                .num_grps        = vb_header_sec->num_lines,
                                .num_alns        = BGEN32 (header.sam_prim_num_alns),
                                .seq_len         = BGEN32 (header.sam_prim_seq_len),
                                .comp_qual_len   = BGEN32 (header.sam_prim_comp_qual_len), // this is just an estimate for PIZ, actual value will be updated after PIZ compresses
                                .qname_len       = BGEN32 (header.sam_prim_qname_len),
                                .comp_cigars_len = BGEN32 (header.sam_prim_comp_cigars_len) }; // this is just an estimate for PIZ, actual value will be updated after PIZ compresses };
        
        total_grps       += plsg[i].num_grps;
        total_alns       += plsg[i].num_alns;
        total_qual_len   += plsg[i].comp_qual_len; 
        total_cigars_len += plsg[i].comp_cigars_len;
        total_qname_len  += plsg[i].qname_len;
        total_seq_len    += (plsg[i].seq_len + 31) & ~(uint32_t)31;   // in bases: each VB's 2bit SEQ is word-aligned (=32 2bit bases), so VBs can be populated in parallel 
    }

    buf_alloc_exact_zero (evb, z_file->sa_groups, total_grps, SAGroupType, "z_file->sa_groups"); // also sets .len ; zero in case some fields not loaded
    buf_alloc_exact_zero (evb, z_file->sa_alns, total_alns, SAAlnType, "z_file->sa_alns");       
    buf_alloc_exact (evb, z_file->sa_qnames, total_qname_len, char, "z_file->sa_qnames");

    z_file->sa_seq.can_be_big = true; 
    buf_alloc_bitarr (evb, &z_file->sa_seq, total_seq_len * 2, "z_file->sa_seq");   // 2 bits per base

    // for in-memory compressed data, the size is an estimate, and we might grow it later, based on the compressed size
    // consumed on the ZIP side (it should be the same, but will not be if ZIP/PIZ codecs are not the same).
    // Note: we allocated 1M more for rANS memory requires, and 5% more than Seg reports just in case.
    z_file->sa_cigars.can_be_big = true; // suppress warning
    buf_alloc (evb, &z_file->sa_cigars, 0, (float)total_cigars_len * 1.05 + 1000000, uint8_t, 1, "z_file->sa_cigars");
    
    z_file->sa_qual.can_be_big = true; // suppress warning
    buf_alloc (evb, &z_file->sa_qual, 0, (float)total_qual_len * 1.05 + 1000000, uint8_t, 1, "z_file->sa_qual"); 

    mutex_initialize (copy_qual_mutex);
    mutex_initialize (copy_cigars_mutex);
    mutex_initialize (num_prim_vbs_loaded_mutex);

    next_plsg_i = num_prim_vbs_loaded = 0; // initialize for pizzing of THIS z_file

    save_flag = flag; // save in global variable

    flag.preprocessing = true; // we're currently dispatching compute threads for preprocessing (= loading SA Groups)
}

uint32_t sam_piz_get_plsg_i (VBIType vb_i)
{
    ASSERTNOTEMPTY (plsg_info);

    ARRAY (PlsgVbInfo, plsg, plsg_info);
    for (int i=0; i < plsg_len; i++)
        if (plsg[i].vblock_i == vb_i)
            return i;

    ABORT_R ("vb_i=%u is not a PRIM VB", vb_i);
}

// PIZ reconstruction compute thread
SAGroupType *sam_piz_get_prim_vb_first_sa_grp (VBlockSAMP vb)
{
    ASSERTNOTEMPTY (plsg_info);

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    return B(SAGroupType, z_file->sa_groups, plsg->first_grp_i);
}

bool sam_is_SA_Groups_loaded (void)
{
    uint32_t num_prim_vbs_loaded_now = __atomic_load_n (&num_prim_vbs_loaded, __ATOMIC_RELAXED);

    return num_prim_vbs_loaded_now == plsg_info.len;
}