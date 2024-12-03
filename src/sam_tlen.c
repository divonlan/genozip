// ------------------------------------------------------------------
//   sam_tlen.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"

//---------
// SEG
//---------

static inline bool sam_seg_predict_TLEN (VBlockSAMP vb, ZipDataLineSAMP dl, bool is_rname_rnext_same,
                                         SamTlenType *predicted_tlen)
{
    PosType32 pnext_pos_delta = dl->PNEXT - dl->POS;

    if (dl->FLAG.supplementary || dl->FLAG.next_unmapped) 
        *predicted_tlen = 0;
    
    else if (!dl->FLAG.multi_segs)
        *predicted_tlen = vb->ref_consumed;

    else if (!is_rname_rnext_same) 
        *predicted_tlen = 0;
    
    else if (!dl->FLAG.rev_comp && dl->FLAG.next_rev_comp) {
        
        uint32_t approx_mate_ref_consumed;

        if (!segconf_running && segconf.has[OPTION_MC_Z]) {
            if (!ctx_encountered_in_line(VB, OPTION_MC_Z)) return false; // in a has[OPTION_MC_Z] file, if we use this formula referring to MC, we need to assure PIZ that MC exists on this line, as it cannot easily check 
            
            approx_mate_ref_consumed = CTX(OPTION_MC_Z)->last_value.i;
        }

        // if this line has a mate, get the mate's ref_consumed. Note: we hope this is our mate, but this is not guaranteed -
        // it is just an earlier read with the same QNAME - could be a supplamentary alignment
        else if (sam_has_mate) 
            approx_mate_ref_consumed = DATA_LINE (vb->mate_line_i)->ref_consumed; // most of the time we're lucky and this is our mate

        // if no MC and no mate, approximate mate's ref_consumed as this line's ref_consumed
        else 
            approx_mate_ref_consumed = vb->ref_consumed;

        *predicted_tlen = pnext_pos_delta + approx_mate_ref_consumed;
    }

    else if (dl->FLAG.rev_comp && !dl->FLAG.next_rev_comp) 
        *predicted_tlen = pnext_pos_delta - vb->ref_consumed + 2 * !dl->FLAG.is_aligned;

    else 
        *predicted_tlen = pnext_pos_delta - vb->ref_consumed;

    return true;
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

    if (segconf_running && tlen_value) segconf.has_TLEN_non_zero = true;

    SamTlenType predicted_tlen;
    if (segconf.has_TLEN_non_zero && sam_seg_predict_TLEN (vb, dl, is_rname_rnext_same, &predicted_tlen)
        && ABS (tlen_value - predicted_tlen) <= 7) {

        // if (predicted_tlen != tlen_value)
        //     printf ("WRONG FLAG=%x tlen=%d expected=%d : line_i=%u POS=%u PNEXT=%u RefConsumed=%u mate_RefConsumed=%u sup=%u next_unmapped=%u is_first=%u rev=%u next_rev=%u aligned=%u\n", 
        //             dl->FLAG.value, (int)tlen_value, (int)predicted_tlen,
        //             vb->line_i, (int)dl->POS, (int)dl->PNEXT, vb->ref_consumed, (int)CTX(OPTION_MC_Z)->last_value.i,
        //             dl->FLAG.supplementary, dl->FLAG.next_unmapped, dl->FLAG.is_first, dl->FLAG.rev_comp, dl->FLAG.next_rev_comp, dl->FLAG.is_aligned);

        SNIPi3 (SNIP_SPECIAL, SAM_SPECIAL_TLEN, '0' + (segconf.has[OPTION_MC_Z] > 0), tlen_value - predicted_tlen);
        seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
    }

    else if (segconf.has_TLEN_non_zero && !segconf_running) // note: only in these cases ltype=LT_DYN_INT is set in sam_seg_initialize
        seg_integer (VB, ctx, tlen_value, true, add_bytes);

    else
        seg_integer_as_snip_do (VB, ctx, tlen_value, add_bytes); // likely all 0, so all-the-same

    ctx->last_value.i = tlen_value;
}

//---------
// PIZ
//---------

static inline PosType32 sam_piz_predict_TLEN (VBlockSAMP vb, bool has_mc)
{
    if (last_flags.supplementary || last_flags.next_unmapped) return 0;

    if (!last_flags.multi_segs) return vb->ref_consumed;

    STRlast (last_rname, SAM_RNAME);
    STRlast (last_rnext, SAM_RNEXT);
         
    if (!OUT_DT(SAM) && *(int32_t*)last_rname != *(int32_t*)last_rnext) return 0;

    if (OUT_DT(SAM) && !IS_EQUAL_SIGN (last_rnext) && !str_issame (last_rname, last_rnext)) return 0;

    PosType32 pnext_pos_delta = (PosType32)CTX(SAM_PNEXT)->last_value.i - (PosType32)CTX(SAM_POS)->last_value.i;

    if (!last_flags.rev_comp && last_flags.next_rev_comp) {

        if (has_mc) {
            STR(MC);
            reconstruct_peek (VB, CTX(OPTION_MC_Z), pSTRa(MC));
            unsigned MC_ref_consumed = sam_cigar_get_MC_ref_consumed (STRa(MC));
            return pnext_pos_delta + MC_ref_consumed;
        }
                
        else if (sam_has_mate) {
            uint32_t mate_ref_consumed = B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->mate_line_i)->ref_consumed;

            return pnext_pos_delta + mate_ref_consumed;
        }

        else
            return pnext_pos_delta + vb->ref_consumed;
    }

    if (last_flags.rev_comp && !last_flags.next_rev_comp) 
        return pnext_pos_delta - vb->ref_consumed + 2 * !last_flags.is_aligned;

    else 
        return pnext_pos_delta - vb->ref_consumed;
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_TLEN)
{
    ASSPIZ0 (snip_len, "snip_len=0");

    bool has_mc = snip[0] - '0';
    int32_t tlen_by_calc = (snip_len > 1) ? atoi (snip+1) : 0;

    new_value->i = tlen_by_calc + sam_piz_predict_TLEN (VB_SAM, has_mc);

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
