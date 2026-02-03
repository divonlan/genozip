// ------------------------------------------------------------------
//   sam_tlen.c
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"

//---------
// SEG
//---------

static SamTlenType sam_seg_predict_TLEN (VBlockSAMP vb, ZipDataLineSAMP dl, bool is_rname_rnext_same)
{
    PosType32 pnext_pos_delta = dl->PNEXT - dl->POS;
    SamTlenType prediction;

    if (dl->FLAG.supplementary || dl->FLAG.next_unmapped) 
        prediction = 0;
    
    else if (!dl->FLAG.multi_segs)
        prediction = vb->ref_consumed;

    else if (!is_rname_rnext_same) 
        prediction = 0;
    
    else if (!IS_PRIM(vb)               && // doesn't yet work in PRIM (bug 1213)
             vb->line_i && (dl->FLAG.duplicate || (dl-1)->FLAG.duplicate) && 
             dl->POS   == (dl-1)->POS   && // fail fast
             dl->PNEXT == (dl-1)->PNEXT &&
             dl->RNAME == (dl-1)->RNAME &&
             dl->RNEXT == (dl-1)->RNEXT)
        prediction = CTX(SAM_TLEN)->last_value.i;

    else if (!dl->FLAG.rev_comp && dl->FLAG.next_rev_comp) {
        
        uint32_t approx_mate_ref_consumed;
        
        // note: we only apply this logic if has[MC], so to save PIZ the need to lookup AUX and MC unnessarily in files that don't have MCs
        // note: until 15.0.75, if has[MC] but not encountered in this line, TLEN was segged as integer (prediction not used)
        if (segconf.has[OPTION_MC_Z] && ctx_encountered_in_line(VB, OPTION_MC_Z)) // TLEN is segged after AUX in bam/sam_seg_txt_line
            approx_mate_ref_consumed = CTX(OPTION_MC_Z)->last_value.i;

        // if this line has a mate that appeared earlier in the same VB, get the mate's ref_consumed. 
        else if (sam_has_mate) 
            approx_mate_ref_consumed = DATA_LINE (vb->mate_line_i)->ref_consumed; // most of the time we're lucky and this is our mate

        // if no MC and no mate, approximate mate's ref_consumed as this line's ref_consumed
        else 
            approx_mate_ref_consumed = vb->ref_consumed;

        prediction = pnext_pos_delta + approx_mate_ref_consumed;
    }

    else if (dl->FLAG.rev_comp && !dl->FLAG.next_rev_comp) 
        prediction = pnext_pos_delta - vb->ref_consumed + 2 * !dl->FLAG.is_aligned;

    else 
        prediction = pnext_pos_delta - vb->ref_consumed;

    // xcons modifications to prediction
    if (segconf.sam_has_xcons) { // 15.0.76
        if (dl->FLAG.multi_segs && !dl->FLAG.is_aligned && prediction < 0)
            prediction -= 2; 

        if (dl->FLAG.rev_comp && dl->FLAG.next_rev_comp)
            prediction = -prediction;

        if (dl->PNEXT == 0)
            prediction = 0;
    }
        
    return prediction;
}

// TLEN - 4 cases: 
// 1. case: a non-zero value that is the negative of its mate in a sorted file - a COPY_MATE_TLEN special
// 2. case: a non-zero value that is the negative of the previous line (usually a mate in a collated file) - a SNIP_DELTA & "-" (= value negation)
// 3. case: tlen>0 and pnext_pos_delta>0 and seq_len>0 tlen is stored as SNIP_SPECIAL & tlen-pnext_pos_delta-seq_len
// 4. otherwise: stored as is
void sam_seg_TLEN (VBlockSAMP vb, ZipDataLineSAMP dl, 
                   STRp(tlen), SamTlenType tlen_value, // option 1 and 2
                   bool is_rname_rnext_same)
{
    Context *ctx = CTX(SAM_TLEN);
    unsigned add_bytes = IS_BAM_ZIP ? sizeof (uint32_t) : tlen_len + 1;

    if (tlen) { // get tlen_value
        ASSSEG0 (tlen_len, "empty TLEN");

        bool is_int = str_get_int_range32 (STRa(tlen), MIN_TLEN, MAX_TLEN, &tlen_value); // note: tlen_value remains 0 if not a valid integer
        ASSSEG (is_int, "expecting TLEN to be an integer [%d,%d], but found \"%.*s\"", MIN_TLEN, MAX_TLEN, STRf(tlen));
    }

    if (segconf_running) {
        if (tlen_value) segconf.has_TLEN_non_zero = true;
    }

    SamTlenType prediction = sam_seg_predict_TLEN (vb, dl, is_rname_rnext_same);

    if (flag.show_tlen_pred && !segconf_running)
        printf ("%s\t%.*s\ttlen=%d %s predicted=%d\tFLAG=%s RNAME≡RNEXT=%s POS=%u PNEXT=%u RefConsumed=%u mate_RefConsumed=%u\n", 
                VB_NAME, STRfQNAME, (int)tlen_value, (prediction==tlen_value ? "=" : "≠"), (int)prediction, 
                sam_dis_FLAG(dl->FLAG).s, VX(is_rname_rnext_same), (int)dl->POS, (int)dl->PNEXT, vb->ref_consumed, (int)CTX(OPTION_MC_Z)->last_value.i);

    if (segconf.has_TLEN_non_zero && ABS (tlen_value - prediction) <= 7) {
        SNIPi4 (SNIP_SPECIAL, SAM_SPECIAL_TLEN, 
                '0' + (segconf.has[OPTION_MC_Z] > 0), 
                '0' + segconf.sam_has_xcons, 
                tlen_value - prediction);
        seg_by_ctx (VB, STRa(snip), ctx, add_bytes);

        if (tlen_value == prediction) vb->num_tlen_pred++;  // for stats
    }

    else if (segconf.has_TLEN_non_zero && !segconf_running) // note: only in these cases ltype=LT_DYN_INT is set in sam_seg_initialize
        seg_integer (VB, ctx, tlen_value, true, add_bytes);

    else
        seg_integer_as_snip_do (VB, ctx, tlen_value, add_bytes); // likely all 0, so all-the-same

    ctx_set_last_value (VB, ctx, (int64_t)tlen_value);
}

//---------
// PIZ
//---------

static SamTlenType sam_piz_predict_TLEN (VBlockSAMP vb, bool has_mc, bool has_xcons)
{
    SamTlenType prediction;
    STRlast (last_rname, SAM_RNAME);
    STRlast (last_rnext, SAM_RNEXT);
    PosType32 pnext_pos_delta = (PosType32)CTX(SAM_PNEXT)->last_value.i - (PosType32)CTX(SAM_POS)->last_value.i;

    if (last_flags.supplementary || last_flags.next_unmapped) 
        prediction = 0;

    else if (!last_flags.multi_segs) 
        prediction = vb->ref_consumed;

    else if (!OUT_DT(SAM) && *(int32_t*)last_rname != *(int32_t*)last_rnext) 
        prediction = 0;

    else if (OUT_DT(SAM) && !IS_EQUAL_SIGN (last_rnext) && !str_issame (last_rname, last_rnext)) 
        prediction = 0;

    else if (VER2(15,76) && !IS_PRIM(vb) && vb->line_i && 
        (last_flags.duplicate || CTX(SAM_FLAG)->prev_flags.duplicate) && 
        CTX(SAM_POS)  ->last_value.i == CTX(SAM_POS)  ->prev_pos &&
        CTX(SAM_PNEXT)->last_value.i == CTX(SAM_PNEXT)->prev_pos &&
        CTX(SAM_RNAME)->last_value.i == CTX(SAM_RNAME)->prev_wi  &&
        CTX(SAM_RNEXT)->last_value.i == CTX(SAM_RNEXT)->prev_wi)

        prediction = CTX(SAM_TLEN)->last_value.i;

    else if (!last_flags.rev_comp && last_flags.next_rev_comp) {
        STR(MC);

        if (has_mc &&
            // note: until 15.0.75, if has_mc this line was guaranteed to have MC. since 15.0.76, prediction is used even if has_mc but line doesn't have MC 
            ({ sam_piz_peek_OPTION (vb, CTX(OPTION_MC_Z), pSTRa(MC), &has_mc); has_mc; }) ) { // peeking confirmed MC exists

            unsigned MC_ref_consumed = sam_cigar_get_MC_ref_consumed (STRa(MC));
            prediction = pnext_pos_delta + MC_ref_consumed;
        }
                
        else if (sam_has_mate) {
            uint32_t mate_ref_consumed = B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->mate_line_i)->ref_consumed;

            prediction = pnext_pos_delta + mate_ref_consumed;
        }

        else
            prediction = pnext_pos_delta + vb->ref_consumed;
    }

    else if (last_flags.rev_comp && !last_flags.next_rev_comp) 
        prediction = pnext_pos_delta - vb->ref_consumed + 2 * !last_flags.is_aligned;

    else 
        prediction = pnext_pos_delta - vb->ref_consumed;

    // xcons modifications to prediction
    if (has_xcons) { // 15.0.76
        if (last_flags.multi_segs && !last_flags.is_aligned && prediction < 0)
            prediction -= 2; 

        if (last_flags.rev_comp && last_flags.next_rev_comp)
            prediction = -prediction;

        if (CTX(SAM_PNEXT)->last_value.i == 0)
            prediction = 0;
    }

    return prediction;
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_TLEN)
{
    ASSPIZ0 (snip_len, "snip_len=0");

    int base_snip_len = VER2(15,76) ? 2 : 1;
    bool has_mc = snip[0] - '0';
    bool has_xcons = (base_snip_len == 2) && (snip[1] - '0'); 

    int32_t tlen_by_calc = (snip_len > base_snip_len) ? atoi (&snip[base_snip_len]) : 0;

    new_value->i = tlen_by_calc + sam_piz_predict_TLEN (VB_SAM, has_mc, has_xcons);

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return true; // new value
}

// Used for files compressed with Genozip up to 12.0.42
SPECIAL_RECONSTRUCTOR (sam_piz_special_TLEN_old)
{
    ASSPIZ0 (snip_len, "snip_len=0");

    int32_t tlen_by_calc = atoi (snip);
    int32_t tlen_val = tlen_by_calc + CTX(SAM_PNEXT)->last_delta + vb->seq_len;

    new_value->i = tlen_val;

    if (reconstruct) RECONSTRUCT_INT (tlen_val);

    return true; // new value
}

// Used for files compressed with Genozip 12.0.41 and 12.0.42
SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_MATE_TLEN_old)
{
    ASSPIZ0 (VB_SAM->mate_line_i >= 0, "No mate line is set for the current line");

    new_value->i = -*B(int64_t, ctx->history, VB_SAM->mate_line_i); // minus the buddy
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return true; // new value
}

// place value in correct location in alignment
TRANSLATOR_FUNC (sam_piz_sam2bam_TLEN)
{
    BAMAlignmentFixedP alignment = (BAMAlignmentFixedP)Btxt (vb->line_start);
    alignment->tlen = LTEN32 ((int32_t)ctx->last_value.i);
    return 0;
}
