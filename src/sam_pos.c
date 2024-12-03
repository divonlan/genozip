// ------------------------------------------------------------------
//   sam_pos.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// a module for handling POS and PNEXT

#include "sam_private.h"
#include "random_access.h"

PosType32 sam_seg_POS (VBlockSAMP vb, ZipDataLineSAMP dl, WordIndex prev_line_chrom, unsigned add_bytes)
{
    ZipDataLineSAMP mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1
    PosType32 pos = dl->POS;
    PosType32 prev_line_pos = vb->line_i ? (dl-1)->POS : 0;

    ContextP pos_ctx = CTX(SAM_POS);

    bool do_mux = IS_MAIN(vb) && segconf.is_paired && !segconf.sam_is_unmapped; // for simplicity. To do: also for prim/depn components
    int channel_i = sam_has_mate?1 : sam_has_prim?2 : 0;
    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, SAM_POS, (MultiplexerP)&vb->mux_POS, channel_i) 
                                  : pos_ctx;

    if (segconf.sam_is_unmapped && !dl->POS)
        seg_by_ctx (VB, "0", 1, pos_ctx, add_bytes); // expecting all-the-same

    // case: DEPN or PRIM line.
    // Note: in DEPN, pos already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    else if (sam_seg_has_sag_by_SA (vb)) {
        sam_seg_against_sa_group (vb, channel_ctx, add_bytes);
        ctx_set_last_value (VB, pos_ctx, (int64_t)pos);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (IS_PRIM(vb)) {
            seg_pos_field (VB, OPTION_SA_POS, OPTION_SA_POS, 0, 0, 0, 0, dl->POS, 0);

            // count POS field contribution to OPTION_SA_POS, so sam_stats_reallocate can allocate the z_data between POS and SA:Z
            CTX(OPTION_SA_POS)->counts.count += add_bytes; 
        }
    }
    
    // case: seg against mate's PNEXT
    else if (sam_has_mate && mate_dl->PNEXT == pos) 
        seg_by_ctx (VB, STRa(copy_mate_PNEXT_snip), channel_ctx, add_bytes); // copy POS from earlier-line mate PNEXT

    // case: seg against POS in the predicted alignment in prim line SA:Z 
    else if (sam_has_prim && sam_seg_is_item_predicted_by_prim_SA (vb, SA_POS, pos)) 
        seg_special1 (VB, SAM_SPECIAL_COPY_PRIM, '0'+SA_POS, channel_ctx, add_bytes); 

    else if (!segconf.is_sorted) 
        seg_integer (VB, channel_ctx, pos, true, add_bytes);

    else 
        pos = seg_pos_field (VB, channel_ctx->did_i, SAM_POS, 0, 0, 0, 0, pos, add_bytes);

    ctx_set_last_value (VB, pos_ctx, (int64_t)pos);

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_POS.snip), SAM_POS, 0); // de-multiplexer

    random_access_update_pos (VB, SAM_POS);

    if (segconf_running) segconf_test_sorted (VB, prev_line_chrom, pos, prev_line_pos);

    return pos;
}

static inline int sam_PNEXT_get_mux_channel (VBlockSAMP vb, bool rnext_is_equal)
{
    return sam_has_mate?0 : sam_has_prim?1 : rnext_is_equal?2 : 3;
}

void sam_seg_PNEXT (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(pnext_str)/* option 1 */, PosType32 pnext/* option 2 */, unsigned add_bytes)
{
    if (pnext_str) 
        ASSSEG (str_get_int_range32 (STRa(pnext_str), 0, MAX_POS_SAM, &pnext),
                "PNEXT=\"%.*s\" ∉ [0,%d] (pnext_str=%p)", pnext_str_len, pnext_str, (int)MAX_POS_SAM, pnext_str);

    if (pnext && segconf_running) 
        segconf.has[SAM_PNEXT]++; // "has" means we found evidence of non-zero PNEXT

    if (segconf.has[SAM_PNEXT]) {
        int channel_i = sam_PNEXT_get_mux_channel (vb, vb->RNEXT_is_equal);
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, SAM_PNEXT, (MultiplexerP)&vb->mux_PNEXT, channel_i);

        // case: copy mate's POS
        if (channel_i==0 && DATA_LINE (vb->mate_line_i)->POS == pnext) 
            seg_by_ctx (VB, STRa(copy_mate_POS_snip), channel_ctx, add_bytes); // copy PNEXT from earlier-line mate POS

        // case: copy prim lines PNEXT
        else if (channel_i==1 && DATA_LINE (vb->saggy_line_i)->PNEXT == pnext) 
            seg_by_ctx (VB, STRa(copy_saggy_PNEXT_snip), channel_ctx, add_bytes); // copy from PNEXT to the channel ctx

        else 
            pnext = seg_pos_field (VB, channel_ctx->did_i, SAM_POS, 0, 0, 0, 0, pnext, add_bytes);

        seg_by_did (VB, STRa(vb->mux_PNEXT.snip), SAM_PNEXT, 0);
    }

    // expecting PNEXT to be 0, but can accommodate if its not
    else {
        if (!pnext) 
            seg_by_did (VB, "0", 1, SAM_PNEXT, add_bytes); // this is expected to be all-the-same
        else
            seg_integer_as_snip_do (VB, CTX(SAM_PNEXT), pnext, add_bytes); // this is not expected to happen usually. note: segged as text to avoid making local and int buffer, which would prevent singletons in the expected case
    }
    
    ctx_set_last_value (VB, CTX(SAM_PNEXT), (int64_t)pnext);
    dl->PNEXT = pnext;
}

// v14: De-multiplex PNEXT
SPECIAL_RECONSTRUCTOR (sam_piz_special_PNEXT)
{
    // compare RNAME and RNEXT without assuming that their word_index's are the same (they might not be in headerless SAM)
    ASSPIZ0 (ctx_has_value_in_line_(vb, CTX(SAM_RNAME)), "RNAME has no value in line");
    ASSPIZ0 (ctx_has_value_in_line_(vb, CTX(SAM_RNEXT)), "RNEXT has no value in line");

    STR(RNAME); ctx_get_snip_by_word_index (CTX(SAM_RNAME), CTX(SAM_RNAME)->last_value.i, RNAME);
    STR(RNEXT); ctx_get_snip_by_word_index (CTX(SAM_RNEXT), CTX(SAM_RNEXT)->last_value.i, RNEXT);
     
    bool rnext_is_equal = IS_EQUAL_SIGN (RNEXT) || str_issame (RNAME, RNEXT);
    
    int channel_i = sam_PNEXT_get_mux_channel (VB_SAM, rnext_is_equal);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// 12.0.41 up to v13: For collated files, in lines where PNEXT equal the previous line's POS.
SPECIAL_RECONSTRUCTOR (sam_piz_special_PNEXT_IS_PREV_POS_old)
{
    new_value->i = CTX(SAM_POS)->last_value.i - CTX(SAM_POS)->last_delta;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return true; // new value
}
