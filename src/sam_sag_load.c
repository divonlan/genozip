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
#include "htscodecs/rANS_static4x16.h"

typedef struct {
    VBIType vblock_i;
    SAGroup first_grp_i; 
    uint32_t num_grps, num_alns;
    uint32_t seq_len;         // in bases
    uint32_t comp_qual_len;   // size as of compressed in-memory in SA Groups. 
    union {
    uint32_t comp_cigars_len; // used by SAG_BY_SA
    uint32_t solo_data_len;   // used by SAG_BY_SOLO 
    };
    uint32_t qname_len;
    uint64_t seq_start, qual_start, qname_start, first_aln_i;
} PlsgVbInfo;

#define plsg_info z_file->vb_info[1]
static uint32_t next_plsg_i = 0; // iterator for dispatching compute threads
static VBIType num_prim_vbs_loaded = 0; // number of PRIM VBs whose loading is complete

static Mutex copy_qual_mutex = {}, copy_cigars_mutex = {};

#define vb_qual_buf   vb->z_data  // used for QUAL data being in-memory compressed in the loader compute thread
#define vb_cigars_buf vb->codec_bufs[1] // similar for cigar data

static Flags save_flag; 

//-------------------------------------------------------------------------
// PIZ: SA Group loader
// Method: decompress directly into memory, parallelizing with dispatcher and 3 threads per PRIM VB: SEQ, QUAL and Alignments
//-------------------------------------------------------------------------

ShowAln sam_show_sag_one_aln (const Sag *g, const SAAln *a)
{
    ShowAln s;
    
    if (command==PIZ) {
        char cigar_info[64];
        if (a->cigar.piz.is_word) 
            snprintf (cigar_info, sizeof (cigar_info), "word=%u", (WordIndex)a->cigar.piz.index);
        else
            snprintf (cigar_info, sizeof (cigar_info), "in=%"PRIu64" l=%u comp=%u", (uint64_t)a->cigar.piz.index, ALN_CIGAR_LEN(a), (int)a->cigar.piz.comp_len);

        char rname_str[64]; // should be big enough...
        snprintf (rname_str, sizeof (rname_str), "\"%.48s\"(%u)", ctx_get_snip_by_word_index0 (ZCTX(VER(15) ? SAM_RNAME/*alias since v15*/ : OPTION_SA_RNAME), a->rname), a->rname);

        snprintf (s.s, sizeof (s.s), "aln_i=%u.%u: sa_rname=%-13s pos=%-10u mapq=%-3u strand=%c nm=%-3u cigar[%u]=\"%s\"(%s)",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), rname_str,
                a->pos, a->mapq, "+-"[a->revcomp], a->nm,
                SA_CIGAR_DISPLAY_LEN, sam_piz_display_aln_cigar (a), cigar_info);
    }
    else 
        snprintf (s.s, sizeof (s.s), "aln_i=%u.%u: cigar_sig=%-12.12s rname=%-4u pos=%-10u mapq=%-3u strand=%c  nm=%-3u",
                ZGRP_I(g), (unsigned)(ZALN_I(a) - g->first_aln_i), cigar_display_signature(a->cigar.signature).s, 
                a->rname, a->pos, a->mapq, "+-"[a->revcomp], a->nm);
    
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
        SoloAln *solo = B(SoloAln, z_file->sag_alns, grp_i);
        for (SoloTags tag_i=0; tag_i < NUM_SOLO_TAGS; tag_i++)
            if (solo->word[tag_i].len)
                SNPRINTF (extra, " %s=\"%.*s\"", 
                          ZCTX(solo_props[tag_i].did_i)->tag_name,
                          solo->word[tag_i].len, B(char, z_file->solo_data, solo->word[tag_i].index));
    }

    iprintf ("grp_i=%u: qname(i=%"PRIu64",l=%u,hash=%u)=%.*s seq=(i=%"PRIu64",l=%u) qual=(i=%"PRIu64",comp_l=%u)[%u]=\"%s\" AS=%u strand=%c mul/fst/lst=%u,%u,%u %s=%u%s%s\n",
             grp_i, (uint64_t)g->qname, g->qname_len, 
             qname_calc_hash (QNAME1, GRP_QNAME(g), g->qname_len, g->is_last, false, NULL), // the hash is part of the index in ZIP, and not part of the SAGroup struct. Its provided here for its usefulness in debugging
             g->qname_len, GRP_QNAME(g) , g->seq, g->seq_len, 
             (uint64_t)g->qual, g->qual_comp_len, SA_QUAL_DISPLAY_LEN, sam_display_qual_from_SA_Group (g),  
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

// QNAME: reconstruct directly into z_file->sa_qname by overlaying txt_data
static inline void sam_load_groups_add_qname (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps) 
{
    // note: multiple threads concurrently reconstruct directly into z_file->sag_qnames
    buf_set_shared (&z_file->sag_qnames);
    buf_overlay_partial (vb, &vb->txt_data, &z_file->sag_qnames, plsg->qname_start, "txt_data");
    
    ContextP seq_len_ctx = segconf.seq_len_dict_id.num ? ECTX(segconf.seq_len_dict_id) : NULL;

    Ltxt = 0;
    for (vb->line_i=0; vb->line_i < plsg->num_grps ; vb->line_i++) {
        sam_reset_line (VB);

        reconstruct_from_ctx (vb, SAM_BUDDY, 0, RECON_OFF); // set buddy (false = don't consume QNAME)

        vb_grps[vb->line_i].qname = plsg->qname_start + Ltxt;
        reconstruct_from_ctx (vb, SAM_QNAME, 0, RECON_ON); // reconstructs into vb->txt_data, sets vb->buddy_line_i if SNIP_COPY_BUDDY
        vb_grps[vb->line_i].qname_len = Ltxt - (vb_grps[vb->line_i].qname - plsg->qname_start); // 64 bit arithmetic

        // if seq_len is carried by a QNAME item, set the last value here (needed for sam_cigar_special_CIGAR)
        if (seq_len_ctx) vb_grps[vb->line_i].seq_len = seq_len_ctx->last_value.i;
    }

    buf_destroy (vb->txt_data);

    reset_iterators (vb);
}

// FLAGS: get multi_segs, is_first, is_last from SAM_FLAGS. 
static inline void sam_load_groups_add_flags (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps) 
{
    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) {
        STR0(snip);
        LOAD_SNIP(SAM_FLAG); 

        // in SAG_SA, FLAG is segged as a SPECIAL
        if (IS_SAG_SA) {
            ASSERT0 (snip[0]==SNIP_SPECIAL && snip[1]==SAM_SPECIAL_SAG, "expecting FLAG to be SAM_SPECIAL_SAG");
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
}

// Grp loader compute thread: copy to z_file->sag_qual
static inline void sam_load_groups_move_comp_to_zfile (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, SAAln *vb_alns) 
{
    plsg->comp_qual_len = vb_qual_buf.len; // update to actual (might not be the same as received from ZIP, if codec has changed)
    
    if (IS_SAG_SA)
        plsg->comp_cigars_len = vb_cigars_buf.len;

    // copy buffers in arbitrary order, based on mutex availability
    bool qual_done = false; 
    bool cigars_done = !IS_SAG_SA; // CIGARs only in SAG_BY_SA
    uint64_t start_qual=0, start_cigars=0;

    while (!qual_done || !cigars_done) {

        bool achieved_something = false;
        
        // QUAL. note: we lock as we might realloc (therefore other threads can't access this Buffer concurrently)
        if (!qual_done && mutex_trylock (copy_qual_mutex)) {
            start_qual = z_file->sag_qual.len;
            buf_add_buf (evb, &z_file->sag_qual, &vb_qual_buf, uint8_t, NULL); // might grow the buffer
            qual_done = achieved_something = true;
            mutex_unlock (copy_qual_mutex);
        }

        // CIGARs
        if (!cigars_done && mutex_trylock (copy_cigars_mutex)) {
            start_cigars = z_file->sag_cigars.len;
            buf_add_buf (evb, &z_file->sag_cigars, &vb_cigars_buf, uint8_t, NULL); 
            cigars_done = achieved_something = true;
            mutex_unlock (copy_cigars_mutex);
        }

        if (!achieved_something) usleep (1000); // 1 ms
    }

    buf_free (vb_qual_buf);
    buf_free (vb_cigars_buf);

    for (SAGroup grp_i=0; grp_i < plsg->num_grps ; grp_i++) 
        vb_grps[grp_i].qual += start_qual;

    if (IS_SAG_SA)
        for (SAGroup aln_i=0; aln_i < plsg->num_alns ; aln_i++) 
            if (!vb_alns[aln_i].cigar.piz.is_word)
                vb_alns[aln_i].cigar.piz.index += start_cigars; // update from index into vb_cigars_buf to z_file->sag_cigars 
}

static inline void sam_load_groups_add_seq (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, uint64_t z_seq_start)
{
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
    Bits *z_sa_seq = (BitsP)&z_file->sag_seq;
    sam_seq_pack (vb, z_sa_seq, z_seq_start * 2, B1STc(vb->textual_seq), vb->seq_len, false, false, HARD_FAIL); 
}

// PIZ: loads QUAL data of PRIM, and compresses it for in-memory storage, using fast rans codec.
static inline void sam_load_groups_add_qual (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g)
{
    if (!CTX(SAM_QUAL)->is_loaded) return; // qual is skipped

    buf_alloc (vb, &vb->codec_bufs[0], vb->seq_len, 0, char, 1, "codec_bufs[0]");

    // to reconstruct into codec_bufs[0], we exchange it with txt_data, as the reconstruct machinery reconstructs to txt_data 
    SWAP (vb->codec_bufs[0], vb->txt_data); // careful as this doesn't modify buffer_list 
    
    // reconstructs into vb->txt_data. note that in PRIM we segged the SPECIAL into QUALSA, not QUAL.
    // so this reconstructs as a LOOKUP from LT_CODEC / LT_BLOB, without going through sam_piz_special_QUAL.
    reconstruct_from_ctx (vb, SAM_QUAL, 0, RECON_ON); 
    buf_verify (vb->txt_data, "SWAP-txt_data-qual");
    SWAP (vb->codec_bufs[0], vb->txt_data);
    
    g->no_qual = (vb->codec_bufs[0].len32 == 1 && *B1STc(vb->codec_bufs[0]) == '*');
    if (g->no_qual) goto done; // no quality data for this group - we're done

    ASSERT (vb->codec_bufs[0].len == vb->seq_len, "Expecting Ltxt=%u == vb->seq_len=%u", Ltxt, vb->seq_len);

    // Compress QUAL, or keep it as is - which ever is better
    uint32_t comp_len = rans_compress_bound_4x16 (g->seq_len, X_NOSZ); // maximum 
    buf_alloc (vb, &vb_qual_buf, comp_len, 0, uint8_t, 1, NULL); // likely already allocated in sam_load_groups_add_grps

    // uncompressed qual is in codec_bufs[0]; append compressed qual to vb_qual_buf
    ASSERT (rans_compress_to_4x16 (VB, B1ST8 (vb->codec_bufs[0]), g->seq_len, BAFT8(vb_qual_buf), &comp_len, X_NOSZ) && comp_len,
            "Failed to compress PRIM qual of vb=%u grp_i=%u qual_len=%u", plsg->vblock_i, ZGRP_I(g), g->seq_len);

    g->qual = vb_qual_buf.len; // relative to the VB. we will update later to be relative to z_file->sag_qual.

    // add to buffer only if compression actually compresses
    if (comp_len < g->seq_len) { 
        g->qual_comp_len = comp_len;  
        vb_qual_buf.len += g->qual_comp_len;
    }
    else {
        buf_add_buf (vb, &vb_qual_buf, &vb->codec_bufs[0], char, NULL); // add uncompressed data instead of compressed
        g->qual_comp_len = 0; // we will use qual_comp_len==0 to mean "not compressed"
    }

done:
    buf_free (vb->codec_bufs[0]); 
}

// populates the cigar data of one alignment
static void sam_load_groups_add_aln_cigar (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *g, SAAln *a, bool is_first_aln,
                                           bool is_all_the_same_LOOKUP, bool is_all_the_same_SQUANK,
                                           pSTRp(out_cigar)) // optional out
{
    ContextP ctx = CTX (OPTION_SA_CIGAR);

    STR0(snip);
    WordIndex word_index = WORD_INDEX_NONE;
    
    // case: word_list.len==0 if all snips are in local, so dictionary was an all-the-same SNIP_LOOKUP, therefore dropped by ctx_drop_all_the_same, 
    if (ctx->word_list.len32)
        word_index = ctx_get_next_snip (VB, ctx, false, pSTRa(snip)); 

    bool is_lookup = is_all_the_same_LOOKUP || (snip_len && *snip == SNIP_LOOKUP);
    bool is_squank = !is_lookup && (is_all_the_same_SQUANK || (snip_len > 2 && snip[0] == SNIP_SPECIAL && snip[1] == SAM_SPECIAL_SQUANK));

    // case: the CIGAR is in local of this vb we need to copy it (compressed) as the "SA Group loader" VB will be soon released. 
    if (is_lookup || is_squank) {

        if (is_lookup)
            ctx_get_next_snip_from_local (VB, ctx, pSTRa(snip)); 

        else { // squank - temporarily reconstruct into txt_data - just for compressing. note: prim cigar is never segged as squank
            sam_piz_special_SQUANK (VB, ctx, &snip[2], snip_len-2, NULL/*reconstruct to vb->scratch*/, true);
            snip     = vb->scratch.data;
            snip_len = vb->scratch.len;
        }

        // note: we put our compressed cigars in vb_cigars_buf, to be copied to z_file->sag_cigars later.
        // this is because we can't assume our compressed length is the same as it was in ZIP, as codecs may change between genozip versions.
        uint32_t this_cigar_mem = rans_compress_bound_4x16 (snip_len, X_NOSZ); // initial allocation, we may grow it later if needed
        uint32_t all_cigar_mem_this_vb = is_first_aln ? (plsg->comp_cigars_len + rans_compress_bound_4x16 (snip_len, X_NOSZ)) : 0; // case first allocation: estimate the total memory for all cigars on alignments of this group, to save allocations
        buf_alloc (vb, &vb_cigars_buf, this_cigar_mem, all_cigar_mem_this_vb, char, 0, "scratch"); 

        uint32_t comp_len = snip_len; // initialize pessimistically- no compression

        if (snip_len > 30) { // no point even trying with shortish cigars
            comp_len = rans_compress_bound_4x16 (snip_len, X_NOSZ); // maximum 
            buf_alloc (vb, &vb_cigars_buf, comp_len, 0, uint8_t, 1, NULL);

            // append compressed cigars to vb_cigars_buf
            ASSERT (rans_compress_to_4x16 (VB, (uint8_t*)STRa(snip), BAFT8(vb_cigars_buf), &comp_len, X_NOSZ) && comp_len,
                    "Failed to compress cigar of vb=%u grp_i=%u aln=%"PRId64" cigar_len=%u cigar=\"%.*s\"", 
                    plsg->vblock_i, ZGRP_I(g), ZALN_I(a), snip_len, STRf(snip));
        }
                
        // case: compression doesn't compress - abandon compression and store uncompressed - set comp_len to 0.
        if (comp_len >= snip_len) {
            memcpy (BAFT8(vb_cigars_buf), snip, snip_len);
            comp_len = 0;
        }

        a->cigar.piz.is_word  = false; // cigar is in z_file->sag_cigars (long, possibly compressed)
        a->cigar.piz.index    = vb_cigars_buf.len;
        a->cigar.piz.comp_len = comp_len; // 0 if not compressed
        a->cigar.piz.len_lo   = snip_len & MAXB(ALN_CIGAR_LEN_BITS_LO);
        a->cigar.piz.len_hi   = (snip_len >> ALN_CIGAR_LEN_BITS_LO);

        vb_cigars_buf.len += comp_len ? comp_len : snip_len;
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
                                            pSTRp(prim_cigar)) // out
{
    // primary alignment CIGAR 
    sam_load_groups_add_aln_cigar (vb, plsg, g, a, ZGRP_I(g)==0, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, STRa(prim_cigar));

    // calculate seq_len from primary CIGAR, and copy cigar to vb->textual_cigar
    sam_cigar_analyze (vb, STRa(*prim_cigar), false, &vb->seq_len);

    // populate SAAln.cigar with the CIGAR word index, for the non-primary alignments of this group
    for (uint8_t aln_i=1; aln_i < g->num_alns; aln_i++) 
        sam_load_groups_add_aln_cigar (vb, plsg, g, a + aln_i, false, is_all_the_same_LOOKUP, is_all_the_same_SQUANK, 0, 0);
}

static inline ZWord reconstruct_to_solo_aln (VBlockSAMP vb, Did did_i, uint64_t solo_data_start, bool check_copy)
{
    decl_ctx (did_i);

    if (!ctx->is_loaded) return (ZWord){};

    uint32_t recon_start = Ltxt;
    reconstruct_from_ctx (VB, did_i, 0, RECON_ON); 
    uint32_t recon_len = Ltxt - recon_start;

    // if identical, UR and UB, CR and CB, point to the same data in solo_data
    if (check_copy && recon_start >= recon_len && !memcmp (Btxt (recon_start), Btxt (recon_start-recon_len), recon_len)) {
        recon_start -= recon_len;
        Ltxt -= recon_len;

        // update history if needed
        HistoryWord *hword = B(HistoryWord, ctx->history, vb->line_i);
        if (ctx->flags.store_per_line && hword->lookup == LookupTxtData)
            hword->index -= recon_len;
    }

    return (ZWord){ .index = solo_data_start + recon_start, .len = recon_len };
}

static inline void sam_load_groups_add_solo_data (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps)
{
    SoloAln *solo_alns = IS_SAG_SOLO ? B(SoloAln, z_file->sag_alns, plsg->first_grp_i) : NULL;

    // get beggining of solo data for this prim vb (adding up lengths of lens of solo data in prev VBs)
    uint64_t solo_data_start = 0;
    
    for (int i=0; i < vb->plsg_i; i++) 
        solo_data_start += B(PlsgVbInfo, plsg_info, i)->solo_data_len;

    // note: multiple threads concurrently reconstruct directly into z_file->solo_data
    buf_overlay_partial (vb, &vb->txt_data, &z_file->solo_data, solo_data_start, "txt_data");
    Ltxt = 0;

    for (vb->line_i=0; vb->line_i < plsg->num_grps ; vb->line_i++) {
        sam_reset_line (VB);
        reconstruct_from_ctx (VB, SAM_BUDDY, 0, RECON_ON); // set buddy (true = consume QNAME)

        // get the index of the tag in this alignment's AUX container. we're really interested if they exist or not.
        ContainerPeekItem idxs[NUM_SOLO_TAGS] = SOLO_CON_PEEK_ITEMS;
        container_peek_get_idxs (VB, CTX(SAM_AUX), ARRAY_LEN(idxs), idxs, true);

        for (SoloTags tag_i=0; tag_i < NUM_SOLO_TAGS; tag_i++)
            if (idxs[tag_i].idx != -1) 
                solo_alns[vb->line_i].word[tag_i] = reconstruct_to_solo_aln (vb, solo_props[tag_i].did_i, solo_data_start, solo_props[tag_i].maybe_same_as_prev);
    }

    CTX(SAM_AUX)->iterator = (SnipIterator){}; // reset as it might be needed for AS
    buf_destroy (vb->txt_data); // un-overlay

    reset_iterators (vb);
}

// SEQ - ACGT-pack uncompressed sequence directly to z_file->sag_seq
// Alns - populate CIGAR the cigar field of the alignments
// Grp  - populate seq, seq_len, num_alns
static inline void sam_load_groups_add_grps (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, CCAln *cc_alns) 
{
    ASSERTNOTINUSE (vb_qual_buf);
    ASSERTNOTINUSE (vb_cigars_buf);

    // note: we put our compressed QUAL in vb_qual_buf, to be copied to z_file->sag_qual later.
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

    if (CTX(OPTION_MC_Z)->is_loaded) // if we're going to reconstruct MC:Z into txt_data
        buf_alloc (vb, &vb->txt_data, 0, vb->recon_size, char, 0, "txt_data");

    // if MC:Z attempts to copy an earlier CIGAR, it will get an empty snip bc CIGAR is not in txt_data. 
    // That's ok, bc this means this MC:Z will not be used to reconstruct its mate's CIGAR as its mate already appeared (and we don't use MC:Z for anything else during pre-processing)
    CTX(SAM_CIGAR)->empty_lookup_ok = true; 

    ContextP seq_len_ctx = segconf.seq_len_dict_id.num ? ECTX(segconf.seq_len_dict_id) : NULL;

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
                sam_load_groups_add_grp_cigars (vb, plsg, g, prim_aln, is_cigar_all_the_same_LOOKUP, is_cigar_all_the_same_SQUANK, pSTRa(prim_aln_cigar));
                    
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

        g->seq     = plsg->seq_start + total_seq_len; // in bases
        g->seq_len = vb->seq_len;   // in bases

        // add SEQ data to the group, packing it in-memory
        if (CTX(SAM_SQBITMAP)->is_loaded)
            sam_load_groups_add_seq (vb, plsg, g, plsg->seq_start + total_seq_len);

        // add QUAL data to the group, possibly compressing it in-memory
        if (CTX(SAM_QUAL)->is_loaded) 
            sam_load_groups_add_qual (vb, plsg, g);
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
            g->as = CAP_SA_AS (vb->last_int(OPTION_AS_i));  // [0,255] capped at 0 and 255
        }

        // if line contains MC:Z - reconstruct it into txt_data, as our mate's CIGAR might copy it
        // note: a better way to do this would be using reconstruct_to_history
        if (!IS_SAG_SA && CTX(OPTION_MC_Z)->is_loaded && idxs[1].idx != -1/*this line has MC:Z*/)  
            reconstruct_from_ctx (VB, OPTION_MC_Z, 0, RECON_ON);

        total_seq_len += vb->seq_len;
        
        vb->textual_seq.len = vb->codec_bufs[0].len = 0;
    }

    buf_free (vb->textual_seq);
    buf_free (vb->codec_bufs[0]);
    buf_free (vb->txt_data);
    reset_iterators (vb);
}

// Alns - populate RNAME, POS, MAPQ, STRAND and NM (CIGAR is added sam_load_groups_add_grps)
static inline void sam_load_groups_add_SA_alns (VBlockSAMP vb, PlsgVbInfo *plsg, Sag *vb_grps, SAAln *vb_alns) 
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
        
        SAAln *a = &vb_alns[aln_i];

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
        if (CTX(OPTION_SA_NM)->is_loaded) { // NM is not skipped
            reconstruct_from_ctx (vb, OPTION_SA_NM, 0, RECON_OFF);     // don't reconstruct, only set last_value to pos
            a->nm = CTX(OPTION_SA_NM)->last_value.i;
        }
    }

    // group revcomp is the revcomp of its primary alignment
    for (uint32_t grp_i=0; grp_i < num_alns_len; grp_i++) 
        vb_grps[grp_i].revcomp = vb_alns[vb_grps[grp_i].first_aln_i - plsg->first_aln_i].revcomp;
}

// entry point of compute thread of sag loading
static void sam_load_groups_add_one_prim_vb (VBlockP vb_)
{
    START_TIMER;
    VBlockSAMP vb = (VBlockSAMP)vb_;

    piz_uncompress_all_ctxs (VB);

    buf_free (vb->z_data); // we will use z_data for vb_qual_buf to avoid allocating a separate buffer

    PlsgVbInfo *plsg = B(PlsgVbInfo, plsg_info, vb->plsg_i);
    Sag *vb_grps     = B(Sag, z_file->sag_grps, plsg->first_grp_i);
    SAAln *vb_alns   = IS_SAG_SA ? B(SAAln, z_file->sag_alns, plsg->first_aln_i) : NULL;
    CCAln *cc_alns   = IS_SAG_CC ? B(CCAln, z_file->sag_alns, plsg->first_grp_i) : NULL;

    vb_grps[0].first_grp_in_vb = true;

    if (IS_SAG_SA)
        sam_load_groups_add_SA_alns (vb, plsg, vb_grps, vb_alns);

    else if (IS_SAG_SOLO)
        sam_load_groups_add_solo_data (vb, plsg, vb_grps);

    if (CTX(SAM_QNAME)->is_loaded) sam_load_groups_add_qname (vb, plsg, vb_grps);
    sam_load_groups_add_flags (vb, plsg, vb_grps);
    sam_load_groups_add_grps (vb, plsg, vb_grps, cc_alns);
    sam_load_groups_move_comp_to_zfile (vb, plsg, vb_grps, vb_alns); // last, as it might block

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
    if (flag.show_sag && sections_get_num_vbs (SAM_COMP_PRIM) == z_file->num_preproc_vbs_joined) {
        sam_show_sag();
        if (is_genocat) exit_ok;
    }

    if (flag_is_show_vblocks (PREPROCESSING_TASK_NAME) || flag_is_show_vblocks (PIZ_TASK_NAME))
        iprintf ("LOADED_SA(id=%d) vb=%s\n", vb->id, VB_NAME);
}

// PIZ main thread: a callback of piz_after_global_area 
void sam_piz_load_sags (void)
{
    next_plsg_i = num_prim_vbs_loaded = 0; // reset for new z_file

    if (flag.genocat_no_reconstruct) return; 

    uint32_t num_prim_vbs = sections_get_num_vbs (SAM_COMP_PRIM);
    if (!num_prim_vbs) return; // no actual PRIM lines (either not gencomp file, or only DEPN lines)

    ARRAY_alloc (PlsgVbInfo, plsg, num_prim_vbs, false, plsg_info, evb, "z_file->plsg");

    uint64_t total_seq_len=0, total_qual_len=0, total_cigars_len=0, total_qname_len=0, total_alns=0; // total across all PRIM VBs of the file
    SAGroup total_grps=0;

    // surveys all PRIM VBs in the order the appear in the file (note: they needn't be consecutive in the file)
    Section vb_header_sec = NULL;
    for (uint32_t i=0; i < num_prim_vbs; i++) {
        
        sections_get_next_vb_header_sec (SAM_COMP_PRIM, &vb_header_sec);

        ASSERT0 (vb_header_sec->st==SEC_VB_HEADER && IS_PRIM(vb_header_sec), "expecting a PRIM VB Header");

        SectionHeaderVbHeader header = zfile_read_section_header (evb, vb_header_sec, SEC_VB_HEADER).vb_header;

        plsg[i] = (PlsgVbInfo){ .vblock_i        = vb_header_sec->vblock_i,
                                .first_grp_i     = BGEN32 (header.sam_prim_first_grp_i),
                                .first_aln_i     = total_alns,
                                .seq_start       = total_seq_len, // in bases
                                .qual_start      = total_qual_len,
                                .qname_start     = total_qname_len,
                                .num_grps        = vb_header_sec->num_lines,
                                .num_alns        = BGEN32 (header.sam_prim_num_sag_alns),
                                .seq_len         = BGEN32 (header.sam_prim_seq_len),
                                .comp_qual_len   = BGEN32 (header.sam_prim_comp_qual_len), // this is just an estimate for PIZ, actual value will be updated after PIZ compresses
                                .qname_len       = BGEN32 (header.sam_prim_qname_len),
                                // note: both the plsg and the header fields are unions with solo_data_len
                                .comp_cigars_len = BGEN32 (header.sam_prim_comp_cigars_len) }; // this is just an estimate for PIZ, actual value will be updated after PIZ compresses 
        
        total_grps       += plsg[i].num_grps;
        total_alns       += plsg[i].num_alns;
        total_qual_len   += plsg[i].comp_qual_len; 
        total_cigars_len += plsg[i].comp_cigars_len; // this is also an alias with solo_data_len
        total_qname_len  += plsg[i].qname_len + 1;   // +1 to leave 1 byte gap, to avoid terminal \0 of the final QNAME of one VB overwriting the first QNAME of the next
        total_seq_len    += (plsg[i].seq_len + 31) & ~(uint32_t)31;   // in bases: each VB's 2bit SEQ is word-aligned (=32 2bit bases), so VBs can be populated in parallel 
    }

    z_file->sag_grps.can_be_big   = z_file->sag_qnames.can_be_big = z_file->sag_alns.can_be_big = 
    z_file->solo_data.can_be_big  = z_file->sag_seq.can_be_big    = z_file->sag_qual.can_be_big = 
    z_file->sag_cigars.can_be_big = true; // suppress warning
    
    buf_alloc_exact_zero (evb, z_file->sag_grps, total_grps, Sag, "z_file->sag_grps"); // also sets .len ; zero in case some fields not loaded
    buf_alloc_exact (evb, z_file->sag_qnames, total_qname_len, char, "z_file->sag_qnames");

    z_file->sag_alns.count = total_alns;

    if (IS_SAG_SA)
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_alns, SAAln, "z_file->sag_alns");       

    else if (IS_SAG_CC)
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_grps, CCAln, "z_file->sag_alns");       

    else if (IS_SAG_SOLO) {
        buf_alloc_exact_zero (evb, z_file->sag_alns, total_grps, SoloAln, "z_file->sag_alns");       
        buf_alloc_exact (evb, z_file->solo_data, total_cigars_len + 65536, char, "z_file->solo_data");  // add some extra as it is created vs reconstruction that sometimes assumes extra bytes available (for peek, or for reconstructing and then abandoning being the same as previous field)
        buf_set_shared (&z_file->solo_data); //  we're going to overlay vb->txt_data on it     
    }
    
    buf_alloc_bits (evb, &z_file->sag_seq, total_seq_len * 2, NOINIT, 0, "z_file->sag_seq");   // 2 bits per base

    // for in-memory compressed data, the size is an estimate, and we might grow it later, based on the compressed size
    // consumed on the ZIP side (it should be the same, but will not be if ZIP/PIZ codecs are not the same).
    // Note: we allocated 1M more for rANS memory requires, and 5% more than Seg reports just in case.
    if (IS_SAG_SA) {
        buf_alloc (evb, &z_file->sag_cigars, 0, (float)total_cigars_len * 1.05 + 1000000, uint8_t, 1, "z_file->sag_cigars");
    
        mutex_initialize (copy_cigars_mutex);
    }

    buf_alloc (evb, &z_file->sag_qual, 0, (float)total_qual_len * 1.05 + 1000000, uint8_t, 1, "z_file->sag_qual"); 

    mutex_initialize (copy_qual_mutex);

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