// ------------------------------------------------------------------
//   sam_sag_load.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "compressor.h"
#include "dispatcher.h"
#include "zfile.h"
#include "qname.h"
#include "writer.h"
#include "huffman.h"
#include "htscodecs/arith_dynamic.h"

typedef struct { // one per PRIM VB
    VBIType vblock_i;
    SAGroup first_grp_i; 
    uint32_t num_grps, num_alns;
    uint64_t first_aln_i;
    uint64_t seq_start; 
    uint32_t seq_len;           // in bases
    uint32_t comp_qual_len;     // size as of compressed in-memory in Sag: huffman or arith compress as in ingest for files since 15.0.68, and arith-compressed for older files (length not know a-priori as compressed ingest using another codec) 
    uint64_t comp_qual_start;
    uint64_t solo_data_start;   // used by SAG_BY_SOLO 
    uint32_t solo_data_len;     // used by SAG_BY_SOLO: length of data for this Plsg in z_file->solo_data (each solo field is huffman-compressed for files since 15.0.68 and uncompressed for older files)
    uint32_t comp_cigars_len;   // used by SAG_BY_SA: length of data for this Plsg in z_file->sag_cigars (each cigar is huffman-compressed for files since 15.0.68 and uncompressed for older files)
    uint64_t comp_cigars_start; // used by SAG_BY_SA
    uint64_t comp_qnames_start;
    uint32_t comp_qnames_len;   // length of data for this Plsg in z_file->sag_qnames (each qname is huffman-compressed for files since 15.0.65 and uncompressed for older files)
} PlsgVbInfo;

#define plsg_info z_file->vb_info[1]
static uint32_t next_plsg_i = 0; // iterator for dispatching compute threads
static VBIType num_prim_vbs_loaded = 0; // number of PRIM VBs whose loading is complete

static Mutex copy_qual_mutex = {}, copy_cigars_mutex = {};

#define vb_qual_buf       vb->z_data        // used for QUAL data being in-memory compressed with ARITH (rather than huffman) in the sag_load compute thread
#define vb_qual_recon_buf vb->codec_bufs[0] // used for reconstructing a qual string when loading
#define vb_cigars_buf     vb->codec_bufs[1] // similar for cigar data (only for files up to 15.0.67)

static Flags save_flag; 

#define IS_PRECISE_QUAL   VER2(15,68) // condition for: precise in-memory QUAL length is known from sam_prim_comp_qual_len:     in files up to 15.0.67 this is now-obsolete rANS compression which we ignore and compress with arith instead, and since 15.0.68 either huffman or arith for which precise length is known
#define IS_PRECISE_CIGARS VER2(15,68) // condition for: precise in-memory CIGARS length is known from sam_prim_comp_cigars_len: in files up to 15.0.67 this is now-obsolete compression scheme, which we ignore and store uncompressed, and since 15.0.68 huffman for which precise length is known

//-------------------------------------------------------------------------
// PIZ: SA Group loader
// Method: decompress directly into memory, parallelizing with dispatcher and 3 threads per PRIM VB: SEQ, QUAL and Alignments
//-------------------------------------------------------------------------

ShowAln sam_show_sag_one_aln (const Sag *g, const SAAln *a)
{
    ShowAln s;
    
    if (IS_PIZ) {
        char cigar_info[64];
        if (a->cigar.piz.is_word) 
            snprintf (cigar_info, sizeof (cigar_info), "word=%u", (WordIndex)a->cigar.piz.index);
        else
            snprintf (cigar_info, sizeof (cigar_info), "in=%"PRIu64" l=%u", (uint64_t)a->cigar.piz.index, a->cigar.piz.len);

        char rname_str[64]; // should be big enough...
        snprintf (rname_str, sizeof (rname_str), "\"%.48s\"(%u)", ctx_get_snip_by_word_index0 (ZCTX(VER(15) ? SAM_RNAME/*alias since v15*/ : OPTION_SA_RNAME), a->rname), a->rname);

        snprintf (s.s, sizeof (s.s), "aln_i=%u.%u: sa_rname=%-13s pos=%-10u mapq=%-3u strand=%c nm=%-3u cigar[%u]=\"%s\"(%s)",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), rname_str,
                a->pos, a->mapq, "+-"[a->revcomp], a->nm,
                SA_CIGAR_DISPLAY_LEN, sam_piz_display_aln_cigar (a), cigar_info);
    }
    else 
        snprintf (s.s, sizeof (s.s), "aln_i=%u.%u: cigar_sig=%016"PRIx64" rname=%-4u pos=%-10u mapq=%-3u strand=%c  nm=%-3u",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), a->cigar.signature, a->rname, a->pos, a->mapq, "+-"[a->revcomp], a->nm);
    
    return s;
}

void sam_show_sag_one_grp (SAGroup grp_i)
{
    const Sag *g = B(Sag, z_file->sag_grps, grp_i);

    StrTextMegaLong extra = {};
    unsigned extra_len=0;

    if (IS_SAG_SA)
        SNPRINTF (extra, " first_aln=%"PRIu64, (uint64_t)g->first_aln_i);

    else if (IS_SAG_CC) {
        CCAln *aln = B(CCAln, z_file->sag_alns, grp_i);
        SNPRINTF (extra, " rname=%s pos=%u", 
                  (command==ZIP) ? ctx_get_z_snip_ex (ZCTX(SAM_RNAME), aln->rname, 0, 0) 
                                 : ctx_get_snip_by_word_index0 (ZCTX(SAM_RNAME), aln->rname), 
                  aln->pos);
    }

    else if (IS_SAG_SOLO) {
        SoloAln *aln = B(SoloAln, z_file->sag_alns, grp_i);
        SNPRINTF (extra, " aln_index=%"PRIu64" (uncomp_len,comp_len)=", ((uint64_t)aln->index_hi) << 32 | aln->index_lo);
        for (SoloTags solo=0; solo < NUM_SOLO_TAGS; solo++)
            if (aln->field_uncomp_len[solo])
                SNPRINTF (extra, " %s=(%u,%u)", ZCTX(solo_props[solo].did_i)->tag_name, aln->field_uncomp_len[solo], aln->field_comp_len[solo]);
    }

    char qname[g->qname_len];
    //xxx huffman_uncompress_or_copy (SAM_QNAME, GRP_QNAME(g), qname, g->qname_len);
    if (huffman_exists (SAM_QNAME)) 
        huffman_uncompress (SAM_QNAME, GRP_QNAME(g), qname, g->qname_len);
    else  
        memcpy (qname, GRP_QNAME(g), g->qname_len);

    iprintf ("grp_i=%u: qname(i=%"PRIu64",l=%u,hash=%u)=\"%.*s\" seq=(i=%"PRIu64",l=%u) qual=(i=%"PRIu64",comp_l=%u)[%u]=\"%s\" AS=%u strand=%c mul/fst/lst=%u,%u,%u %s=%u%s%s\n",
             grp_i, (uint64_t)g->qname, g->qname_len, 
             (uint32_t)qname_calc_hash (QNAME1, COMP_NONE, qname, g->qname_len, g->is_last, false, CRC32, NULL), // the hash is part of the index in ZIP, and not part of the SAGroup struct. Its provided here for its usefulness in debugging
             g->qname_len, qname, (uint64_t)g->seq, g->seq_len, 
             (uint64_t)g->qual, g->qual_comp_len, SA_QUAL_DISPLAY_LEN, sam_display_qual_from_sag (g).s,  
             (int)g->as, "+-"[g->revcomp], g->multi_segs, g->is_first, g->is_last,
             ((IS_SAG_NH || IS_SAG_CC || IS_SAG_SOLO) ? "NH" : "num_alns"), g->num_alns, 
             extra.s, g->first_grp_in_vb ? " FIRST_GRP_IN_VB" : "");

    if (IS_SAG_SA)
        for (uint64_t aln_i=g->first_aln_i; aln_i < g->first_aln_i + g->num_alns; aln_i++) 
            iprintf ("    %s\n", sam_show_sag_one_aln (g, B(SAAln, z_file->sag_alns, aln_i)).s);
}

void sam_show_sag (void)
{
    if (flag.show_sag >= 1) { // show a specific group
        ASSINP (flag.show_sag <= z_file->sag_grps.len, "--show-sag argument, for this file, should be between 0 and %"PRIu64, 
                z_file->sag_grps.len-1);
        
        sam_show_sag_one_grp (flag.show_sag - 1);
    }

    else {
        iprintf ("SAG_TYPE=%s NUM_SA_GROUPS=%"PRIu64" NUM_ALIGNMENTS=%"PRIu64"\n", 
                 sag_type_name (segconf.sag_type), z_file->sag_grps.len, z_file->sag_alns.count);

        for (SAGroup grp_i=0; grp_i < z_file->sag_grps.len; grp_i++) 
            sam_show_sag_one_grp (grp_i);
    }
}

static void reset_iterators (VBlockSAMP vb)
{
    CTX(SAM_AUX  )->iterator = (SnipIterator){}; 
    CTX(SAM_QNAME)->iterator = (SnipIterator){};
    CTX(SAM_BUDDY)->iterator = (SnipIterator){};
    CTX(SAM_BUDDY)->next_local = 0;
    ctx_unset_encountered (VB, CTX(SAM_AUX));
}

static inline void sam_load_groups_add_qname (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps) 
{
    START_TIMER;
    
    if (!CTX(SAM_QNAME)->is_loaded) return;

    ASSERTNOTINUSE (vb->txt_data);

    buf_alloc (vb, &vb->txt_data,  0, SAM_MAX_QNAME_LEN * plsg->num_grps, char, 0, "txt_data");

    decl_ctx (SAM_QNAME);
    bool is_huff = huffman_exists (SAM_QNAME);

    ContextP seq_len_ctx = segconf.seq_len_dict_id.num ? ECTX(segconf.seq_len_dict_id) : NULL;


    uint32_t max_comp_len = is_huff ? huffman_get_theoretical_max_comp_len (SAM_QNAME, SAM_MAX_QNAME_LEN) : 0;

    const uint8_t *start_qnames = B1ST8(z_file->sag_qnames);
    const uint8_t *first_qnames = B8(z_file->sag_qnames, plsg->comp_qnames_start); // first qname for this plsg
    const uint8_t *after_qnames = first_qnames + plsg->comp_qnames_len; // after qname data of this plsg
    uint8_t *next_qnames = (uint8_t *)first_qnames; 

    for (vb->line_i=0; vb->line_i < plsg->num_grps ; vb->line_i++) {
        sam_reset_line (VB);
        reconstruct_from_ctx (vb, SAM_BUDDY, 0, RECON_OFF); // set buddy (false = don't consume QNAME)
        
        rom uncomp_qname = BAFTtxt;
        ctx->mate_copied_exactly = false; 
        reconstruct_from_ctx (vb, SAM_QNAME, 0, RECON_ON);  // reconstructs into vb->txt_data, sets vb->buddy_line_i if SNIP_COPY_BUDDY
        uint32_t uncomp_qname_len = BAFTtxt - uncomp_qname;
        
        ASSPIZ (uncomp_qname_len <= SAM_MAX_QNAME_LEN, "unexpectedly, uncomp_qname_len=%u > SAM_MAX_QNAME_LEN=%u", uncomp_qname_len, SAM_MAX_QNAME_LEN);

        vb_grps[vb->line_i].qname_len = uncomp_qname_len; // always uncompressed length

        if (ctx->mate_copied_exactly)
            vb_grps[vb->line_i].qname = vb_grps[ctx->last_value.i].qname; // last_value is buddy_line_i     
    
        else {
            ASSERT (after_qnames > next_qnames, "PRIM/%u: plsg->comp_qnames_len=%u exceeded: (after-next)=%"PRId64, 
                    vb->plsg_i, plsg->comp_qnames_len, (uint64_t)(after_qnames - next_qnames));

            vb_grps[vb->line_i].qname = next_qnames - start_qnames;
        
            // since 15.0.65, a SEC_HUFFMAN section is available and we can compress, before that, we just copy
            //xxx next_qnames += huffman_compress_or_copy (SAM_QNAME, STRa(uncomp_qname), next_qnames, max_comp_len); // using huffman sent via SEC_HUFFMAN
            if (is_huff) {
                uint32_t comp_len = max_comp_len; // theoretical maximum 
                huffman_compress (SAM_QNAME, STRa(uncomp_qname), next_qnames, &comp_len); // using huffman sent via SEC_HUFFMAN
                next_qnames += comp_len;
            }

            // up to 15.0.64, store uncompressed
            else 
                next_qnames = mempcpy (next_qnames, uncomp_qname, uncomp_qname_len);
        }

        // if seq_len is carried by a QNAME item, set the last value here (needed for sam_cigar_special_CIGAR)
        if (seq_len_ctx) vb_grps[vb->line_i].seq_len = seq_len_ctx->last_value.i;
    }

    // in files up to 15.0.64, the uncompressed length didn't consider mates, so its longer than needed
    ASSERT (is_huff ? (after_qnames == next_qnames) : (after_qnames >= next_qnames), 
            "PRIM/%u: Expecting qnames for this prim vb to have length=%u, but its length=%"PRIu64,
            vb->plsg_i, plsg->comp_qnames_len, (uint64_t)(next_qnames - first_qnames));

    buf_destroy (vb->txt_data);

    reset_iterators (vb);

    COPY_TIMER (sam_load_groups_add_qname);
}

// FLAGS: get multi_segs, is_first, is_last from SAM_FLAGS. 
static inline void sam_load_groups_add_flags (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps) 
{
    START_TIMER;

    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {
        STR0(snip);
        LOAD_SNIP(SAM_FLAG); 

        // in SAG_SA, FLAG is segged as a SPECIAL
        if (IS_SAG_SA) {
            ASSERT0 (snip[0]==SNIP_SPECIAL && snip[1]==SAM_SPECIAL_pull_from_sag, "expecting FLAG to be SAM_SPECIAL_pull_from_sag");
            snip += 2;
        }

        SamFlags sam_flags = { .value = atoi (snip) };
        vb_grps[grp_i].multi_segs  = sam_flags.multi_segs;
        vb_grps[grp_i].is_first    = sam_flags.is_first;
        vb_grps[grp_i].is_last     = sam_flags.is_last;
        
        // case: non-SAG_SA - revcomp is segged in SAM_FLAG. Note: For SA-defined, it is always 0, and we retrieve it later from OPTION_SA_STRAND of the primary alignment
        if (!IS_SAG_SA)
            vb_grps[grp_i].revcomp = sam_flags.rev_comp; 
    }

    COPY_TIMER (sam_load_groups_add_flags);
}

// Grp loader compute thread: copy to z_file->sag_*
static inline void sam_load_groups_move_comp_to_zfile (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, SAAln *vb_alns) 
{
    START_TIMER;
    
    // copy buffers in arbitrary order, based on mutex availability. no need to copy if precise length was known and hence data was written directly.
    bool qual_done   = IS_PRECISE_QUAL;  
    bool cigars_done = IS_PRECISE_CIGARS || !IS_SAG_SA; 

    uint64_t start_qual=0, start_cigars=0;

    while (!qual_done || !cigars_done) {
        bool achieved_something = false;
        
        // note: we lock as we might realloc (therefore other threads can't access this Buffer concurrently)
        #define APPEND(x, vb_buf, max_index)                                    \
            if (!x##_done && mutex_trylock (copy_##x##_mutex)) {                \
                start_##x = z_file->sag_##x.len;                                \
                buf_append_buf (evb, &z_file->sag_##x, &vb_buf, uint8_t, NULL); \
                x##_done = achieved_something = true;                           \
                ASSERT (z_file->sag_##x.len <= max_index, "PRIM/%u: while loading SAGs: z_file->sag_"#x".len=%"PRIu64" > "#max_index"=%"PRIu64, \
                        plsg->vblock_i, z_file->sag_##x.len, max_index);        \
                mutex_unlock (copy_##x##_mutex);                                \
            }

        APPEND (qual,   vb_qual_buf,   MAX_SA_QUAL_INDEX);
        APPEND (cigars, vb_cigars_buf, MAX_SA_CIGAR_INDEX);

        // case: all muteces were locked by other threads - wait a bit
        if (!achieved_something) { 
            START_TIMER;
            usleep (1000); // 1 ms
            COPY_TIMER (sam_load_groups_move_comp_to_zfile_idle);
        }
    }

    // for buffers building in the VB first and then copied to z_file, update the indices
    if (!IS_PRECISE_QUAL) {
        buf_free (vb_qual_buf);
        for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) 
            vb_grps[grp_i].qual += start_qual;
    }

    if (!IS_PRECISE_CIGARS && IS_SAG_SA) {
        buf_free (vb_cigars_buf);
        for (SAGroup aln_i=0; aln_i < plsg->num_alns ; aln_i++) 
            if (!vb_alns[aln_i].cigar.piz.is_word)
                vb_alns[aln_i].cigar.piz.index += start_cigars; 
    }

    COPY_TIMER (sam_load_groups_move_comp_to_zfile);
}

static inline void sam_load_groups_add_seq (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, uint64_t z_seq_start)
{
    START_TIMER;

    if (!CTX(SAM_SQBITMAP)->is_loaded) return; // sequence is skipped

    // reconstruct SEQ to vb->textual_seq, needed by sam_piz_prim_add_QUAL
    buf_alloc (vb, &vb->textual_seq, 0, vb->seq_len, char, 0, "textual_seq");

    SWAP (vb->textual_seq, vb->txt_data); // careful as this doesn't modify buffer_list 
    reconstruct_from_ctx (vb, SAM_SQBITMAP, 0, RECON_ON);
    buf_verify (vb->txt_data, "SWAP-txt_data-textual_seq");
    SWAP (vb->textual_seq, vb->txt_data);

    vb->textual_seq_str = B1STc (vb->textual_seq);

    ASSERT (vb->textual_seq.len32 == vb->seq_len, "Expecting textual_seq.len=%u == seq_len=%u", vb->textual_seq.len32, vb->seq_len);

    // pack SEQ data into z_file->sag_seq
    BitsP z_sa_seq = (BitsP)&z_file->sag_seq;
    { START_TIMER; 
    sam_seq_pack (vb, z_sa_seq, z_seq_start * 2, B1STc(vb->textual_seq), vb->seq_len, false, false, HARD_FAIL); 
    COPY_TIMER (sam_load_groups_add_seq_pack); }

    COPY_TIMER (sam_load_groups_add_seq);
}

// PIZ: loads QUAL data of PRIM, and compresses it for in-memory storage, using fast huffman or arith compression.
static inline void sam_load_groups_add_qual (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g,
                                             bytes start_qual, bytes after_qual, uint8_t **next_qual)
{
    START_TIMER;
   
    buf_alloc (vb, &vb_qual_recon_buf, 0, vb->seq_len, char, 1, "codec_bufs[0]");

    // to reconstruct into codec_bufs[0], we exchange it with txt_data, as the reconstruct machinery reconstructs to txt_data 
    SWAP (vb_qual_recon_buf, vb->txt_data); // careful as this doesn't modify buffer_list 
    
    // reconstructs into vb->txt_data. note that in PRIM we segged the SPECIAL into QUALSA, not QUAL.
    // so this reconstructs as a LOOKUP from LT_CODEC / LT_BLOB, without going through sam_piz_special_QUAL.
    reconstruct_from_ctx (vb, SAM_QUAL, 0, RECON_ON); 
    buf_verify (vb->txt_data, "SWAP-txt_data-qual");
    SWAP (vb_qual_recon_buf, vb->txt_data);
    
    g->no_qual = (vb_qual_recon_buf.len32 == 1 && *B1STc(vb_qual_recon_buf) == '*');
    if (g->no_qual) goto done; // no quality data for this group - we're done

    ASSERT (vb_qual_recon_buf.len == vb->seq_len, "Expecting Ltxt=%u == vb->seq_len=%u", Ltxt, vb->seq_len);

    uint32_t comp_len = 0;
    
    // arith: possibly precise (since 15.0.68, in which case we move to next_qual) or not (in which case we keep in vb_qual_buf and move later)
    if (!huffman_exists (SAM_QUAL)) {
        comp_len = CTX(SAM_QUAL)->qual_longest_theoretical_comp_len; 
        buf_alloc (vb, &vb_qual_buf, comp_len, 0, uint8_t, 1, NULL); // likely already allocated in sam_load_groups_add_grps

        ASSERT (arith_compress_to (VB, B1ST8 (vb_qual_recon_buf), g->seq_len, BAFT8(vb_qual_buf), &comp_len, X_NOSZ) && comp_len,
                "Failed to arith_compress PRIM qual of vb=%u grp_i=%u qual_len=%u", plsg->vblock_i, ZGRP_I(g), g->seq_len);
    }

    // case: compress directly into z_file->sag_qual: huffman or arith
    if (IS_PRECISE_QUAL) {
        ASSERT (after_qual > *next_qual, "PRIM/%u: plsg->comp_qual_len=%u exceeded: (after-next)=%"PRId64, 
                vb->plsg_i, plsg->comp_qual_len, (uint64_t)(after_qual - *next_qual));

        // xxx comp_len = huffman_compress_or_copy (SAM_QUAL, B1STc (vb_qual_recon_buf), g->seq_len, *next_qual, 
        //            MIN_(after_qual - *next_qual, (uint64_t)CTX(SAM_QUAL)->qual_longest_theoretical_comp_len));

        if (huffman_exists (SAM_QUAL)) {
            comp_len = MIN_(after_qual - *next_qual, (uint64_t)CTX(SAM_QUAL)->qual_longest_theoretical_comp_len);
            huffman_compress (SAM_QUAL, B1STc (vb_qual_recon_buf), g->seq_len, *next_qual, &comp_len);
        }

        else    
            memcpy (*next_qual, B1ST8(vb_qual_buf), comp_len); // arith
        
        g->qual = *next_qual - start_qual;
        g->qual_comp_len = comp_len; 
        *next_qual += comp_len;  
        vb_qual_buf.len32 = 0;
    }
    
    // case: keep arith-compressed qual in vb_qual_buf - to be copied to z_file->sag_qual later
    else {
        g->qual = vb_qual_buf.len; // relative to the VB. we will update later to be relative to z_file->sag_qual.
        g->qual_comp_len = comp_len;  
        vb_qual_buf.len += comp_len;
    }

done:
    buf_free (vb_qual_recon_buf); 
    COPY_TIMER (sam_load_groups_add_qual);
}

// populates the cigar data of one alignment
static void sam_load_groups_add_aln_cigar (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, SAAln *a, bool is_first_aln,
                                           bool is_all_the_same_LOOKUP, bool is_all_the_same_SQUANK,
                                           bytes start_cigars, bytes after_cigars, uint8_t **next_cigars, 
                                           pSTRp(out_cigar)) // optional out
{
    START_TIMER;

    ContextP ctx = CTX (OPTION_SA_CIGAR);

    STR0(snip);
    WordIndex word_index = WORD_INDEX_NONE;
    
    // case: word_list.len==0 if all snips are in local, so dictionary was an all-the-same SNIP_LOOKUP, therefore dropped by ctx_drop_all_the_same, 
    if (ctx->word_list.len32)
        word_index = ctx_get_next_snip (VB, ctx, false, pSTRa(snip)); 

    bool is_lookup = is_all_the_same_LOOKUP || (snip_len && *snip == SNIP_LOOKUP);
    bool is_squank = !is_lookup && (is_all_the_same_SQUANK || (snip_len > 2 && snip[0] == SNIP_SPECIAL && snip[1] == SAM_SPECIAL_SQUANK));

    if (!start_cigars)
        buf_alloc (vb, &vb_cigars_buf, 0, 128 KB, char, 0, "vb->codec_bufs[1]"); // initial alloction
        
    // case: the CIGAR is in local of this vb we need to copy it (compressed) as the "SA Group loader" VB will be soon released. 
    if (is_lookup || is_squank) {

        if (is_lookup)
            ctx_get_next_snip_from_local (VB, ctx, pSTRa(snip)); 

        else { // squank - temporarily reconstruct into txt_data - just for compressing. note: prim cigar is never segged as squank
            sam_piz_special_SQUANK (VB, ctx, &snip[2], snip_len-2, NULL/*reconstruct to vb->scratch*/, true);
            snip     = vb->scratch.data;
            snip_len = vb->scratch.len;
        }

        // case: file starting 15.0.68: compress with huffman
        if (start_cigars) {
            ASSERT (after_cigars > *next_cigars, "PRIM/%u: plsg->comp_cigars_len=%u exceeded: (after-next)=%"PRId64, 
                    vb->plsg_i, plsg->comp_cigars_len, (uint64_t)(after_cigars - *next_cigars));

            a->cigar.piz.index = *next_cigars - start_cigars;
    
            uint32_t comp_len = MIN_(after_cigars - *next_cigars, 0xffffffffULL); // careful that it doesn't go beyond 32b

            huffman_compress (SAM_CIGAR, STRa(snip), *next_cigars, &comp_len);
            *next_cigars += comp_len;
        }

        // file up to 15.0.67: store uncompressed (in genozip up to 15.0.67 we use to compress with RANS CIGARs that are longer than 30 characters)
        // we put our compressed cigars in vb_cigars_buf, to be copied to z_file->sag_cigars later,
        // because we don't a-priori know the length of cigar data for each plsg
        else {
            a->cigar.piz.index = vb_cigars_buf.len;  
            buf_append (VB, vb_cigars_buf, char, snip, snip_len, NULL);
        }

        a->cigar.piz.is_word = false; // cigar is in z_file->sag_cigars (long, possibly compressed)
        a->cigar.piz.len     = snip_len;
    }

    // case: the cigar is in dict - we therefore don't store - reconstruction will copy it from the dictionary
    else {
        ASSERT (word_index >= 0 && word_index < ctx->word_list.len32, "word_index=%d of %s ∉ [0,%d] snip_len=%u snip[0]=%u snip=\"%.*s\"", 
                word_index, ctx->tag_name, ctx->word_list.len32-1, snip_len, snip_len ? (uint8_t)snip[0] : 0, STRf(snip));

        a->cigar.piz.is_word = true; // cigar is in ZCTX(OPTION_SA_CIGAR).dict (short cigar)
        a->cigar.piz.index   = word_index;
    }

    if (out_cigar) STRset (*out_cigar, snip);
    if (is_squank) buf_free (vb->scratch);
}

static void sam_load_groups_add_grp_cigars (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, SAAln *a, 
                                            bool is_all_the_same_LOOKUP, bool is_all_the_same_SQUANK,
                                            bytes start_cigars, bytes after_cigars, uint8_t **next_cigars,
                                            pSTRp(prim_cigar)) // out
{
    START_TIMER;

    // primary alignment CIGAR 
    sam_load_groups_add_aln_cigar (vb, plsg, g, a, ZGRP_I(g)==0, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, start_cigars, after_cigars, next_cigars, STRa(prim_cigar));

    // calculate seq_len from primary CIGAR, and copy cigar to vb->textual_cigar
    sam_cigar_analyze (vb, STRa(*prim_cigar), false, &vb->seq_len);

    // populate SAAln.cigar with the CIGAR word index, for the non-primary alignments of this group
    for (uint32_t aln_i=1; aln_i < g->num_alns; aln_i++) 
        sam_load_groups_add_aln_cigar (vb, plsg, g, a + aln_i, false, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, start_cigars, after_cigars, next_cigars, 0, 0);

    COPY_TIMER (sam_load_groups_add_grp_cigars);
}

static inline int32_t reconstruct_to_solo_aln (VBlockSAMP vb, Did did_i, qSTR8p(comp), pSTRp(prev), bool check_copy)
{
    decl_ctx (did_i);

    if (!ctx->is_loaded) {
        *comp_len = 0;
        return 0;
    }

    rom recon = BAFTtxt;
    reconstruct_from_ctx (VB, did_i, 0, RECON_ON); // note: we don't reset Ltxt after each field, because some fields rely on previous fields
    uint32_t recon_len = BAFTtxt - recon;

    // if identical, UR and UB, CR and CB, point to the same data in solo_data
    if (check_copy && str_issame (recon, *prev)) {
        *comp_len = 0;
        return recon_len; // note: no need to update prev, as its identical
    }

    // xxx *comp_len = huffman_compress_or_copy (did_i, STRa(recon), comp, *comp_len);

    if (huffman_exists (did_i))  // true since 15.0.68
        huffman_compress (did_i, STRa(recon), comp, comp_len);

    else { // files up to 15.0.67 - store solo uncompressed
        memcpy (comp, recon, recon_len);
        *comp_len = recon_len;
    }

    STRset (*prev, recon);

    return recon_len;
}

static inline void sam_load_groups_add_solo_data (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps)
{
    START_TIMER;

    SoloAln *solo_alns = B(SoloAln, z_file->sag_alns, plsg->first_grp_i);

    // get begining of solo data for this prim vb (adding up lengths of lens of solo data in prev VBs)
    const uint8_t *start_solo = B1ST8(z_file->solo_data);
    const uint8_t *first_solo = B8(z_file->solo_data, plsg->solo_data_start); // first solo data for this plsg
    const uint8_t *after_solo = first_solo + plsg->solo_data_len; // after solo data of this plsg
    uint8_t *next_solo = (uint8_t *)first_solo; 

    // for files since 15.0.68, we store only solos that appeared in segconf, i.e. have a huffman
    bool has_huff[NUM_SOLO_TAGS];
    for (int solo = 0; solo < NUM_SOLO_TAGS; solo++)
        has_huff[solo] = !VER2(15,68) || huffman_exists (solo_props[solo].did_i); 

    // note: multiple threads concurrently reconstruct directly into z_file->solo_data
    buf_alloc (vb, &vb->txt_data, 0, NUM_SOLO_TAGS * MAX_SOLO_UNCOMP_LEN * plsg->num_grps, char, 0, "txt_data"); // 1 MB is plenty for the Solo tags supported 

    for (vb->line_i=0; vb->line_i < plsg->num_grps ; vb->line_i++) {
        sam_reset_line (VB);

        reconstruct_from_ctx (VB, SAM_BUDDY, 0, RECON_ON); // set buddy 

        // get the index of the tag in this alignment's AUX container. we're really interested if they exist or not.
        ContainerPeekItem idxs[NUM_SOLO_TAGS] = SOLO_CON_PEEK_ITEMS;
        container_peek_get_idxs (VB, CTX(SAM_AUX), ARRAY_LEN(idxs), idxs, true);

        STR0(prev); // pointer into txt_data of previous tag reconstructed        
        SoloAln *solo_aln = &solo_alns[vb->line_i];

        ASSERT0 ((next_solo - start_solo) < 1 TB, "sag_alns for Solo is too big"); // this should never happen in PIZ as it was already tested in ZIP
        solo_aln->index_lo = (next_solo - start_solo) & 0xffffffff;
        solo_aln->index_hi = (next_solo - start_solo) >> 32; 

        for (SoloTags solo=0; solo < NUM_SOLO_TAGS; solo++)
            if (has_huff[solo] && idxs[solo].idx != -1) {
                ASSERT (after_solo > next_solo, "PRIM/%u: plsg->solo_data_len=%u exceeded: (after-next)=%"PRId64, 
                        vb->plsg_i, plsg->solo_data_len, (uint64_t)(after_solo - next_solo));

                uint32_t comp_len = MIN_(after_solo - next_solo, MAX_SOLO_COMP_LEN);

                int32_t uncomp_len = reconstruct_to_solo_aln (vb, solo_props[solo].did_i, next_solo, &comp_len, pSTRa(prev), solo_props[solo].maybe_same_as_prev);
                ASSERTNOTZERO (uncomp_len);
                
                solo_aln->field_uncomp_len[solo] = uncomp_len;

                if (comp_len) { // not a copy of the previous field
                    solo_aln->field_comp_len[solo] = comp_len;
                    next_solo += comp_len;
                }
                else {         // field is a copy of the previous field
                    solo_aln->field_comp_len[solo] = solo_aln->field_comp_len[solo-1];
                    solo_aln->field_comp_len[solo-1] = 0;     
                }
            }
    }

    ASSERT (after_solo == next_solo, "PRIM/%u: Expecting solo data for this prim vb to have length=%u, but its length=%"PRIu64,
            vb->plsg_i, plsg->solo_data_len, (uint64_t)(next_solo - first_solo));

    CTX(SAM_AUX)->iterator = (SnipIterator){}; // reset as it might be needed for AS
    buf_free (vb->txt_data); 

    reset_iterators (vb);

    COPY_TIMER (sam_load_groups_add_solo_data);
}

// SEQ - ACGT-pack uncompressed sequence directly to z_file->sag_seq
// Alns - populate CIGAR the cigar field of the alignments
// Grp  - populate seq, seq_len, num_alns
static inline void sam_load_groups_add_grps (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, CCAln *cc_alns) 
{
    START_TIMER;

    ASSERTNOTINUSE (vb_qual_buf);
    ASSERTNOTINUSE (vb_cigars_buf);

    // note: we put our compressed QUAL in vb_qual_buf, to be copied to z_file->sag_qual later.
    // this is because we can't assume our compressed length is the same as it was in ZIP, as codecs may change between genozip versions.
    uint32_t qual_comp_len_this_vb = plsg->comp_qual_len + (huffman_exists (SAM_QUAL) ? 0 : arith_compress_bound (vb_grps[0].seq_len, X_NOSZ)); // initial allocation, we may grow it later if needed
    buf_alloc (vb, &vb_qual_buf, qual_comp_len_this_vb, 0, char, 0, "z_data"); 

    // get seq and seq_len from CIGAR
    uint32_t total_seq_len = 0;

    uint32_t sa_rname_len = CTX(OPTION_SA_RNAME)->word_list.len32;

    ContextP cigar_ctx = CTX(OPTION_SA_CIGAR);
    bool is_cigar_all_the_same_LOOKUP = !cigar_ctx->word_list.len && cigar_ctx->dict.len && *B1STc(cigar_ctx->dict)==SNIP_LOOKUP;
    bool is_cigar_all_the_same_SQUANK = !cigar_ctx->word_list.len && cigar_ctx->dict.len >= 3 && 
                                        *B1STc(cigar_ctx->dict)==SNIP_SPECIAL && *Bc(cigar_ctx->dict, 1)==SAM_SPECIAL_SQUANK;

    if (CTX(OPTION_MC_Z)->is_loaded) // if we're going to reconstruct MC:Z into txt_data
        buf_alloc (vb, &vb->txt_data, 0, vb->recon_size, char, 0, "txt_data");    

    // if MC:Z attempts to copy an earlier CIGAR, it will get an empty snip bc CIGAR is not in txt_data. 
    // That's ok, bc this means this MC:Z will not be used to reconstruct its mate's CIGAR as its mate already appeared (and we don't use MC:Z for anything else during pre-processing)
    CTX(SAM_CIGAR)->empty_lookup_ok = true; 

    ContextP seq_len_ctx = segconf.seq_len_dict_id.num ? ECTX(segconf.seq_len_dict_id) : NULL;

    const uint8_t *start_cigars=0, *first_cigars=0, *after_cigars=0, *start_qual=0, *first_qual=0, *after_qual=0;
    uint8_t *next_cigars=0, *next_qual=0;

    // starting 15.0.68, we compress CIGARs directly to z_file->sag_cigars
    if (IS_SAG_SA && IS_PRECISE_CIGARS) { 
        start_cigars = B1ST8(z_file->sag_cigars);
        first_cigars = B8(z_file->sag_cigars, plsg->comp_cigars_start); // first cigar for this plsg
        after_cigars = first_cigars + plsg->comp_cigars_len; // after cigar data of this plsg
        next_cigars  = (uint8_t *)first_cigars; 
    }

    if (IS_PRECISE_QUAL) { 
        start_qual = B1ST8(z_file->sag_qual);
        first_qual = B8(z_file->sag_qual, plsg->comp_qual_start); // beginning of qual data for this plsg
        after_qual = first_qual + plsg->comp_qual_len; // after qual data of this plsg
        next_qual  = (uint8_t *)first_qual; 
    }

    // every alignment in the Primary VB represents a grp
    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {
        Sag *g = &vb_grps[grp_i];

        vb->line_i = grp_i;
        sam_reset_line (VB);

        reconstruct_from_ctx (VB, SAM_BUDDY, 0, RECON_ON); // set buddy (true = consume QNAME)

        // if this is an SA:Z-defined SAG, load primary alignment
        if (IS_SAG_SA) {

            ASSERT (g->first_aln_i >= plsg->first_aln_i && g->first_aln_i < z_file->sag_alns.len,
                    "grp_i=%u z_grp_i=%u has first_aln_i=%"PRIu64" ∉ [0,%"PRId64"]", 
                    grp_i, ZGRP_I(g), (uint64_t)g->first_aln_i, z_file->sag_alns.len-1);

            SAAln *prim_aln = B(SAAln, z_file->sag_alns, g->first_aln_i); // first Aln of group

            // add all alignment cigars
            STR(prim_aln_cigar);
            if (CTX(OPTION_SA_CIGAR)->is_loaded) // not skipped
                sam_load_groups_add_grp_cigars (vb, plsg, g, prim_aln, is_cigar_all_the_same_LOOKUP, is_cigar_all_the_same_SQUANK, 
                                                start_cigars, after_cigars, &next_cigars, 
                                                pSTRa(prim_aln_cigar));
                    
            // populate vb with alignment details
            vb->last_int(SAM_POS) = prim_aln->pos;
            CTX(SAM_FLAG)->last_value.i = prim_aln->revcomp ? SAM_FLAG_REV_COMP : 0; // needed by codec_longr_reconstruct

            ASSERT (prim_aln->rname >= 0 && prim_aln->rname < sa_rname_len, 
                    "rname=%d out of range: OPTION_SA_RNAME.word_len.len=%u. primary alignment of grp_i=%u: QNAME=(%"PRIu64",%u)=%.*s, CIGAR[12]=%*.*s",
                    prim_aln->rname, sa_rname_len, grp_i, (uint64_t)g->qname, (int)g->qname_len, (int)g->qname_len, Bc(z_file->sag_qnames, g->qname), MIN_(prim_aln_cigar_len,12),MIN_(prim_aln_cigar_len,12), prim_aln_cigar);

            // case: since v15, OPTION_SA_RNAME is an dict alias of SAM_RNAME
            if (VER(15))
                vb->chrom_node_index = prim_aln->rname; 
        
            // case: in v14, SAM_RNAME and OPTION_SA_RNAME are independent, so we lookup by name
            else {
                ctx_get_snip_by_word_index (CTX(OPTION_SA_RNAME), prim_aln->rname, vb->chrom_name);
                
                vb->chrom_node_index = ctx_get_word_index_by_snip (VB, CTX(SAM_RNAME), STRa(vb->chrom_name)); // convert OPTION_SA_RNAME word_index to RNAME word_index
                ASSERT (vb->chrom_node_index != WORD_INDEX_NONE, "Cannot find rname=%u in context RNAME", prim_aln->rname);
            }
        }

        // non-SAG_SA 
        else {
            // case: SAG_NH, SAG_CC, SAG_SOLO: we get num_alns from NH:i
            if (IS_SAG_NH || IS_SAG_CC || IS_SAG_SOLO) {
                reconstruct_from_ctx (VB, OPTION_NH_i, 0, RECON_OFF);
                g->num_alns = vb->last_int (OPTION_NH_i);
            }

            // needed to reconstruct SEQ
            CTX(SAM_FLAG)->last_value.i = g->revcomp ? SAM_FLAG_REV_COMP : 0; // needed by codec_longr_reconstruct - g->revcomp was populated in sam_load_groups_add_flags

            // set vb->chrom_* (needed to reconstruct SEQ)
            reconstruct_from_ctx (VB, SAM_RNAME, 0, RECON_OFF);
            vb->chrom_node_index = vb->last_int (SAM_RNAME);
            ctx_get_snip_by_word_index (CTX(SAM_RNAME), vb->chrom_node_index, vb->chrom_name);

            // set POS.last_value (which also requires RNEXT and PNEXT)
            reconstruct_from_ctx (VB, SAM_POS,   0, RECON_OFF);
            reconstruct_from_ctx (VB, SAM_RNEXT, 0, RECON_OFF); // needed by PNEXT...
            reconstruct_from_ctx (VB, SAM_PNEXT, 0, RECON_OFF); // store in history, in case POS of future-line mates needs it (if copy_buddy(PNEXT))
            
            // need to set before reconstructing CIGAR: if seq_len is carried by a QNAME item, set the last value here
            if (seq_len_ctx) ctx_set_last_value (VB, seq_len_ctx, (int64_t)g->seq_len);

            // analyze CIGAR (setting vb->seq_len, vb->ref_consumed etc)
            reconstruct_from_ctx (VB, SAM_CIGAR, 0, RECON_OFF);

            if (IS_SAG_CC) 
                cc_alns[grp_i] = (CCAln){ .rname = vb->chrom_node_index, 
                                          .pos   = vb->last_int(SAM_POS) };
        }

        ASSERT (plsg->seq_start + total_seq_len <= MAX_SA_SEQ_INDEX, "PRIM/%u: while loading SAGs: g->seq=%"PRIu64" > MAX_SA_SEQ_INDEX=%"PRIu64, 
                plsg->vblock_i, plsg->seq_start + total_seq_len, MAX_SA_SEQ_INDEX);

        ASSERT (vb->seq_len <= MAX_SA_SEQ_LEN, "PRIM/%u: while loading SAGs: g->seq_len=%u > MAX_SA_SEQ_LEN=%u", 
                plsg->vblock_i, vb->seq_len, MAX_SA_SEQ_LEN);

        g->seq     = plsg->seq_start + total_seq_len; // in bases
        g->seq_len = vb->seq_len;   // in bases

        // add SEQ data to the group, packing it in-memory
        if (CTX(SAM_SQBITMAP)->is_loaded)
            sam_load_groups_add_seq (vb, plsg, g, plsg->seq_start + total_seq_len);

        // add QUAL data to the group, possibly compressing it in-memory
        if (CTX(SAM_QUAL)->is_loaded) 
            sam_load_groups_add_qual (vb, plsg, g, start_qual, after_qual, &next_qual);
        else
            g->no_qual = true; // prevent attempting to reconstruct QUAL if its not loaded

        // get the index of the tag in this alignment's AUX container. we're really interested if they exist or not.
        ContainerPeekItem idxs[2] = { { _OPTION_AS_i, -1 }, { _OPTION_MC_Z, -1 } };
        container_peek_get_idxs (VB, CTX(SAM_AUX), ARRAY_LEN(idxs), idxs, true);

        // set AS.last_value, unless AS:i is skipped entirely
        if (segconf.sag_has_AS          && // depn lines' AS:i was segged against prim 
            CTX(OPTION_AS_i)->is_loaded && // we didn't skip AS:i due to subsetting command line options
            idxs[0].idx != -1) {           // this line has AS:i

            reconstruct_from_ctx (VB, OPTION_AS_i, 0, RECON_OFF);
            g->as = CAP_SA_AS (vb->last_int(OPTION_AS_i));  // capped at 0 and 65535
        }

        // if line contains MC:Z - reconstruct it into txt_data, as our mate's CIGAR might copy it
        // note: a better way to do this would be using reconstruct_to_history
        if (!IS_SAG_SA && CTX(OPTION_MC_Z)->is_loaded && idxs[1].idx != -1/*this line has MC:Z*/)  
            reconstruct_from_ctx (VB, OPTION_MC_Z, 0, RECON_ON);

        total_seq_len += vb->seq_len;
        
        vb->textual_seq.len = 0;
    }

    ASSERT (after_cigars == next_cigars || !start_cigars || !CTX(OPTION_SA_CIGAR)->is_loaded, 
            "PRIM/%u: Expecting CIGAR data for this prim vb to have compressed length=%u, but it is=%"PRIu64,
            vb->plsg_i, plsg->comp_cigars_len, (uint64_t)(next_cigars - first_cigars));

    ASSERT (after_qual == next_qual || !start_qual || !CTX(SAM_QUAL)->is_loaded, 
            "PRIM/%u: Expecting QUAL data for this prim vb to have compressed length=%u, but it is=%"PRIu64,
            vb->plsg_i, plsg->comp_qual_len, (uint64_t)(next_qual - first_qual));

    // update to compressed lengths actual if not precise
    if (!start_qual)   plsg->comp_qual_len   = vb_qual_buf.len;  
    if (!start_cigars) plsg->comp_cigars_len = vb_cigars_buf.len;

    buf_free (vb->textual_seq);
    buf_free (vb->txt_data);
    reset_iterators (vb);

    COPY_TIMER (sam_load_groups_add_grps);
}

// Alns - populate RNAME, POS, MAPQ, STRAND and NM (CIGAR is added sam_load_groups_add_grps)
static inline void sam_load_groups_add_SA_alns (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps) 
{
    START_TIMER;

    ContextP sa_ctx    = LOADED_CTX(OPTION_SA_Z);
    ContextP rname_ctx = CTX(OPTION_SA_RNAME); // possibly not loaded

    ASSERT (sa_ctx->local.len == plsg->num_grps, "PRIM/%u: Mismatch in number of PRIM lines: expecting OPTION_SA_Z.local.len=%"PRIu64" == num_lines=%u OPTION_SA_Z.is_loaded=%s",
            plsg->vblock_i, sa_ctx->local.len, plsg->num_grps, TF(sa_ctx->is_loaded));
    
    ASSERT (sa_ctx->ltype == LT_UINT8 || sa_ctx->ltype == LT_UINT16, "Unexpected OPTION_SA_Z.ltype=%s", lt_name (sa_ctx->ltype));

    for (uint32_t grp_i=0; grp_i < plsg->num_grps; grp_i++) {
        uint64_t first_aln_i       = grp_i ? vb_grps[grp_i-1].first_aln_i + vb_grps[grp_i-1].num_alns : plsg->first_aln_i;
        ASSERT (first_aln_i <= MAX_SA_GRP_ALNS, "PRIM/%u: while loading SAGs vb_grp_i=%u: first_aln_i%"PRIu64" > MAX_SA_GRP_ALNS=%"PRIu64, 
                plsg->vblock_i, grp_i, first_aln_i, MAX_SA_GRP_ALNS);

        vb_grps[grp_i].first_aln_i = first_aln_i;

        uint32_t grp_num_alns = (sa_ctx->ltype == LT_UINT8) ? *B8(sa_ctx->local, grp_i) : *B16(sa_ctx->local, grp_i);

        ASSERT (grp_num_alns <= (uint32_t)MAX_SA_NUM_ALNS, "PRIM/%u: while loading SAGs vb_grp_i=%u: grp_num_alns=%d > MAX_SA_NUM_ALNS=%u",
                plsg->vblock_i, grp_i, grp_num_alns, MAX_SA_NUM_ALNS);

        vb_grps[grp_i].num_alns = grp_num_alns;
    }

    int32_t vb_num_alns = (vb_grps[plsg->num_grps-1].first_aln_i - plsg->first_aln_i) + vb_grps[plsg->num_grps-1].num_alns;

    ASSERT (vb_num_alns == plsg->num_alns, "PRIM/%u: Alignment count mismatch: expecting vb_num_alns=%d == plsg->num_alns=%u",
            plsg->vblock_i, vb_num_alns, plsg->num_alns);

    uint32_t sa_rname_len = rname_ctx->word_list.len32;

    ARRAY (SAAln, alns, z_file->sag_alns);

    for (uint32_t grp_i=0; grp_i < plsg->num_grps; grp_i++) 
        for (uint32_t aln_i=vb_grps[grp_i].first_aln_i; aln_i < vb_grps[grp_i].first_aln_i + vb_grps[grp_i].num_alns; aln_i++) {
            SAAln *a = &alns[aln_i];

            // RNAME
            if (rname_ctx->is_loaded) { // RNAME not skipped
                reconstruct_from_ctx (vb, OPTION_SA_RNAME, 0, RECON_OFF);  // don't reconstruct, only set last_value to word_index
                a->rname = rname_ctx->last_value.i;         // WordIndex into OPTION_SA_RNAME, NOT RNAME!

                ASSERT (a->rname >= 0 && a->rname < sa_rname_len, 
                        "While loading aln_i=%u: rname=%u out of range: OPTION_SA_RNAME.word_len.len=%u", aln_i, a->rname, sa_rname_len);
            }

            // POS
            if (CTX(OPTION_SA_POS)->is_loaded) {
                reconstruct_from_ctx (vb, OPTION_SA_POS, 0, RECON_OFF);    // don't reconstruct, only set last_value to pos
                a->pos = CTX(OPTION_SA_POS)->last_value.i;
            }

            // STRAND - nodes initialized in sam_seg_0X_initialize to 0="-" 1="+"
            if (CTX(OPTION_SA_STRAND)->is_loaded) {
                reconstruct_from_ctx (vb, OPTION_SA_STRAND, 0, RECON_OFF); // don't reconstruct, only set last_value to word_index
                a->revcomp = !CTX(OPTION_SA_STRAND)->last_value.i;
            }
            
            // MAPQ
            if (CTX(OPTION_SA_MAPQ)->is_loaded) { // MAPQ is not skipped
                reconstruct_from_ctx (vb, OPTION_SA_MAPQ, 0, RECON_OFF);   // don't reconstruct, only set last_value to pos
                a->mapq = CTX(OPTION_SA_MAPQ)->last_value.i;
            }

            // NM
            if (CTX(OPTION_SA_NM)->is_loaded &&   // NM is not skipped
                !(segconf.SA_NM_by_CIGAR_X && aln_i == vb_grps[grp_i].first_aln_i)) { // in case of SA_NM_by_CIGAR_X, we get the NM of the primary alignment from the its CIGAR 
                reconstruct_from_ctx (vb, OPTION_SA_NM, 0, RECON_OFF);     // don't reconstruct, only set last_value to pos
                a->nm = CTX(OPTION_SA_NM)->last_value.i;
            }
        }

    // group revcomp is the revcomp of its primary alignment
    for (uint32_t grp_i=0; grp_i < plsg->num_grps; grp_i++) 
        vb_grps[grp_i].revcomp = alns[vb_grps[grp_i].first_aln_i].revcomp;

    COPY_TIMER (sam_load_groups_add_SA_alns);
}

// entry point of compute thread of sag loading
static void sam_load_groups_add_one_prim_vb (VBlockP vb_)
{
    START_TIMER;
    VBlockSAMP vb = (VBlockSAMP)vb_;

    piz_uncompress_all_ctxs (VB, PUR_SAM_LOAD_SAG);

    buf_free (vb->z_data); // we will use z_data for vb_qual_buf to avoid allocating a separate buffer

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    Sag *vb_grps     = B(Sag, z_file->sag_grps, plsg->first_grp_i);
    SAAln *vb_alns   = IS_SAG_SA ? B(SAAln, z_file->sag_alns, plsg->first_aln_i) : NULL;
    CCAln *cc_alns   = IS_SAG_CC ? B(CCAln, z_file->sag_alns, plsg->first_grp_i) : NULL;

    vb_grps[0].first_grp_in_vb = true;

    if (IS_SAG_SA)
        sam_load_groups_add_SA_alns (vb, plsg, vb_grps);

    else if (IS_SAG_SOLO)
        sam_load_groups_add_solo_data (vb, plsg, vb_grps);

    sam_load_groups_add_qname (vb, plsg, vb_grps);
    sam_load_groups_add_flags (vb, plsg, vb_grps);
    sam_load_groups_add_grps  (vb, plsg, vb_grps, cc_alns);

    // move compressed QNAME, QUAL, CIGAR to from vb to z_file 
    sam_load_groups_move_comp_to_zfile (vb, plsg, vb_grps, vb_alns); 

    __atomic_fetch_add (&num_prim_vbs_loaded, (int)1, __ATOMIC_ACQ_REL);

    vb_set_is_processed (VB); // tell dispatcher this thread is done and can be joined.

    COPY_TIMER_EVB (sam_load_groups_add_one_prim_vb);                          
}

// PIZ main thread - dispatches a compute thread do read SA Groups from one PRIM VB. returns true if dispatched
bool sam_piz_dispatch_one_load_sag_vb (Dispatcher dispatcher)
{
    if (!VER(14) || flag.genocat_no_reconstruct || flag.header_only
        || (flag.one_vb && IS_MAIN (sections_vb_header (flag.one_vb)))) { // --one-vb of a MAIN VB - no need for PRIM/DEPN
        flag.preprocessing = false;
        return false; // no need to load sags
    }

    // check if we're done
    if (next_plsg_i == plsg_info.len) goto done_loading;

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, next_plsg_i++);

    VBlockP vb = dispatcher_generate_next_vb (dispatcher, plsg->vblock_i, SAM_COMP_PRIM);
    vb->comp_i          = SAM_COMP_PRIM;
    vb->preprocessing   = true; // we're done dispatching preprocessing VBs, however, preprocessing VB compute threads may still be running.
    vb->show_containers = false; 

    piz_read_one_vb (vb, false); // only the needed contexts - filtered by sam_piz_is_skip_section
    dispatcher_compute (dispatcher, sam_load_groups_add_one_prim_vb);

    return true;

done_loading: // note: at this point, gencomp load VBs are still running in compute threads
    // restore globals
    flag = save_flag;  // also resets flag.preprocessing
    dispatcher_set_task_name (dispatcher, PIZ_TASK_NAME);

    return false;
}

// PIZ main thread: after joining a preprocessing VB
void sam_piz_after_preproc (VBlockP vb)
{
    z_file->num_preproc_vbs_joined++;

    // case: done loading - all groups are loaded
    if (sections_get_num_vbs (SAM_COMP_PRIM) == z_file->num_preproc_vbs_joined) {
        if (flag.show_sag) {
            sam_show_sag();
            if (is_genocat) exit_ok;
        }

        sam_gencomp_trim_memory();
    }

    if (flag_is_show_vblocks (PREPROCESSING_TASK_NAME) || flag_is_show_vblocks (PIZ_TASK_NAME))
        iprintf ("LOADED_SA(id=%d) vb=%s\n", vb->id, VB_NAME);
}

static void sam_sag_load_alloc_z (BufferP buf, bool is_precise, uint64_t precise_len, uint64_t estimated_len, MutexP mutex, rom buf_name)
{
    // if total length is known precisely, we allocate it, and each thread puts the data directly into its location
    if (is_precise) 
        buf_alloc_exact (evb, *buf, precise_len, uint8_t, buf_name);    
    
    // if total length is only an estimate, each thread will put the data into a VB buffer, and it will be added to the z_file buffer under mutex lock
    else { 
        buf_alloc (evb, buf, 0, estimated_len, uint8_t, 1, buf_name); 
        mutex_initialize (*mutex); 
    }
}

// PIZ main thread: a callback of piz_after_global_area 
void sam_piz_load_sags (void)
{
    next_plsg_i = num_prim_vbs_loaded = 0; // reset for new z_file

    if (flag.genocat_no_reconstruct) return; 

    uint32_t num_prim_vbs = sections_get_num_vbs (SAM_COMP_PRIM);
    if (!num_prim_vbs) return; // no actual PRIM lines (either not gencomp file, or only DEPN lines)

    ARRAY_alloc (PlsgVbInfo, plsg, num_prim_vbs, false, plsg_info, evb, "z_file->plsg");

    uint64_t total_seq_len=0, total_qual_len=0, total_cigars_len=0, total_solo_len=0, total_qname_len=0, total_alns=0; // total across all PRIM VBs of the file
    SAGroup total_grps=0;

    // surveys all PRIM VBs in the order of vb_i (i.e. order of creation) (note: this is not the order they appear in the file)
    Section vb_header_sec = NULL;
    for (uint32_t i=0; i < num_prim_vbs; i++) {
        
        sections_get_next_vb_header_sec (SAM_COMP_PRIM, &vb_header_sec);

        ASSERT0 (vb_header_sec->st==SEC_VB_HEADER && IS_PRIM(vb_header_sec), "expecting a PRIM VB Header");

        SectionHeaderVbHeader header = zfile_read_section_header (evb, vb_header_sec, SEC_VB_HEADER).vb_header;

        plsg[i] = (PlsgVbInfo){ .vblock_i          = vb_header_sec->vblock_i,
                                .first_grp_i       = BGEN32 (header.sam_prim_first_grp_i),
                                .first_aln_i       = total_alns,
                                .seq_start         = total_seq_len, // in bases
                                .seq_len           = BGEN32 (header.sam_prim_seq_len), // in bases
                                .num_grps          = vb_header_sec->num_lines,
                                .num_alns          = BGEN32 (header.sam_prim_num_sag_alns),
                                .comp_qual_start   = IS_PRECISE_QUAL ? total_qual_len : 0, 
                                .comp_qual_len     = IS_PRECISE_QUAL ? BGEN32 (header.sam_prim_comp_qual_len) : 0, // note: Value from files up to 15.0.57 is ignored, because we used rANS 
                                .comp_qnames_start = total_qname_len,
                                .comp_qnames_len   = BGEN32 (header.sam_prim_comp_qname_len) };  // huffman-compressed length since 15.0.65, uncompressed length up to 15.0.64 

        if (IS_SAG_SOLO) {
            // since 15.0.68 this is precise length of huffman-compressed Solo data (ZIP does ROUNDUP8 to full Bits words)
            // up to 15.0.67 this is precise length of uncompressed solo data (and load stores uncompressed too - since there is no huffman)
            plsg[i].solo_data_start = total_solo_len;
            plsg[i].solo_data_len = BGEN32 (header.sam_prim_solo_data_len); // note: round to word-align to avoid vbs(=threads) modifying the same 64b word in z_file->solo_data
            total_solo_len += ROUNDUP8 (plsg[i].solo_data_len);
        }

        else if (IS_SAG_SA && IS_PRECISE_CIGARS) {
            // since 15.0.68 this is precise length of huffman-compressed CIGAR data (ZIP does ROUNDUP8 to full Bits words)
            // up to 15.0.67 this is the length of a no-longer-used compression scheme, so we keep it at zero
            plsg[i].comp_cigars_start = total_cigars_len;
            plsg[i].comp_cigars_len = BGEN32 (header.sam_prim_comp_cigars_len);
            total_cigars_len += ROUNDUP8 (plsg[i].comp_cigars_len);
        }

        // note: for fields stored in Bits (huffman-compressed or seq-packed) - roundup to avoid having two threads access the same 64b-word
        total_seq_len    += ROUNDUP32 (plsg[i].seq_len);  // 32 bases = 8 bytes
        total_qname_len  += ROUNDUP8 (plsg[i].comp_qnames_len);   
        total_grps       += plsg[i].num_grps;
        total_alns       += plsg[i].num_alns;
        total_qual_len   += IS_PRECISE_QUAL ? ROUNDUP8 (plsg[i].comp_qual_len)   // since 15.0.58, comp_qual_len is an exact length of huffman / arith compression. 
                          : VER(15)         ? (BGEN32 (header.longest_seq_len) * plsg[i].num_grps / 3) // v15: estimate length: we use longest_seq_len which is available since v15
                          :                   (151 * plsg[i].num_grps / 3);      // v14: estimate length: we arbitrarily take seq_len=151
    }

    z_file->sag_grps.can_be_big   = z_file->sag_qnames.can_be_big = z_file->sag_alns.can_be_big = 
    z_file->solo_data.can_be_big  = z_file->sag_seq.can_be_big    = z_file->sag_qual.can_be_big = 
    z_file->sag_cigars.can_be_big = true; // suppress warning

    z_file->sag_alns.count = total_alns;

    buf_alloc_exact_zero (evb, z_file->sag_grps, total_grps, Sag, "z_file->sag_grps"); // also sets .len ; zero in case some fields not loaded
    
    // SEQ: we know the precise length     
    buf_alloc_bits (evb, &z_file->sag_seq, total_seq_len * 2, NOINIT, 0, "z_file->sag_seq");   // 2 bits per base

    // QUAL: we know the precise length only since 15.0.68
    sam_sag_load_alloc_z (&z_file->sag_qual, IS_PRECISE_QUAL, total_qual_len, total_qual_len, &copy_qual_mutex, "z_file->sag_qual");

    // note: since 15.0.65 for total_qname_len is total huffman-compressed length and we count
    // on the huffman alg in ZIP to be identical to PIZ. prior to that, it was total uncompress length
    buf_alloc_exact (evb, z_file->sag_qnames, total_qname_len, uint8_t, "z_file->sag_qnames");

    if (IS_SAG_SA) {
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_alns, SAAln, "z_file->sag_alns");       
        
        sam_sag_load_alloc_z (&z_file->sag_cigars, IS_PRECISE_CIGARS, total_cigars_len, 
                              MAX_(total_grps * 4, 1 MB), // very imprecise estimate
                              &copy_cigars_mutex, "z_file->sag_cigars");
    }

    else if (IS_SAG_CC)
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_grps, CCAln, "z_file->sag_alns");       

    // SOLO: we know the precise length     
    else if (IS_SAG_SOLO) { // we know the exact length we are going to store in memory: huffman-compressed since 15.0.68, not compressed for earlier files
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_grps, SoloAln, "z_file->sag_alns");       
        buf_alloc_exact (evb, z_file->solo_data, total_solo_len/*actually total_solo_len*/, char, "z_file->solo_data");  
    }

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

    ABORT ("vb_i=%u is not a PRIM VB", vb_i);
}

// PIZ reconstruction compute thread
Sag *sam_piz_get_prim_vb_first_sa_grp (VBlockSAMP vb)
{
    ASSERTNOTEMPTY (plsg_info);
    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    return B(Sag, z_file->sag_grps, plsg->first_grp_i);
}

bool sam_is_sag_loaded (void)
{
    uint32_t num_prim_vbs_loaded_now = __atomic_load_n (&num_prim_vbs_loaded, __ATOMIC_ACQUIRE);

    return num_prim_vbs_loaded_now == plsg_info.len;
}