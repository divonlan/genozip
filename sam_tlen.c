// ------------------------------------------------------------------
//   sam_tlen.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "codec.h"

//---------
// SEG
//---------

static inline bool sam_seg_predict_TLEN (VBlockSAM *vb, ZipDataLineSAM *dl, bool is_rname_rnext_same,
                                          int64_t *predicted_tlen)
{
    int64_t pnext_pos_delta = dl->PNEXT - dl->POS;

    //printf ("RNEXT=%.*s\n", vb->last_txt_len(SAM_RNAME), last_txt(vb, SAM_RNAME));
    if (dl->FLAG.bits.supplementary || dl->FLAG.bits.next_unmapped) 
        *predicted_tlen = 0;
    
    else if (!dl->FLAG.bits.multi_segments)
        *predicted_tlen = vb->ref_consumed;

    else if (!is_rname_rnext_same) 
        *predicted_tlen = 0;
    
    else if (!dl->FLAG.bits.rev_comp && dl->FLAG.bits.next_rev_comp) {
        
        uint32_t approx_mate_ref_consumed;

        if (!segconf.running && segconf.has[OPTION_MC_Z]) {
            if (!ctx_encountered_in_line(VB, OPTION_MC_Z)) return false; // in a has[OPTION_MC_Z] file, if we use this formula referring to MC, we need to assure PIZ that MC exists on this line, as it cannot easily check 
            
            approx_mate_ref_consumed = CTX(OPTION_MC_Z)->last_value.i;
        }

        // if this line has a buddy, get the buddy's ref_consumed. Note: we hope this is our mate, but this is not guaranteed -
        // it is just an earlier read with the same QNAME - could be a supplamentary alignment
        else if (vb->buddy_line_i != -1) 
            approx_mate_ref_consumed = DATA_LINE (vb->buddy_line_i)->ref_consumed; // most of the time we're lucky and this is our mate

        // if no MC and no buddy, approximate buddy's ref_consumed as this line's ref_consumed
        else 
            approx_mate_ref_consumed = vb->ref_consumed;

        *predicted_tlen = pnext_pos_delta + approx_mate_ref_consumed;
    }

    else if (dl->FLAG.bits.rev_comp && !dl->FLAG.bits.next_rev_comp) 
        *predicted_tlen = pnext_pos_delta - vb->ref_consumed + 2 * !dl->FLAG.bits.is_aligned;

    else 
        *predicted_tlen = pnext_pos_delta - vb->ref_consumed;

    return true;
}

// TLEN - 4 cases: 
// 1. case: a non-zero value that is the negative of its mate in a sorted file - a COPY_BUDDY_TLEN special
// 2. case: a non-zero value that is the negative of the previous line (usually a mate in a collated file) - a SNIP_DELTA & "-" (= value negation)
// 3. case: tlen>0 and pnext_pos_delta>0 and seq_len>0 tlen is stored as SNIP_SPECIAL & tlen-pnext_pos_delta-seq_len
// 4. otherwise: stored as is
void sam_seg_TLEN (VBlockSAM *vb, ZipDataLineSAM *dl, 
                   STRp(tlen), int64_t tlen_value, // option 1 and 2
                   bool is_rname_rnext_same)
{
    Context *ctx = CTX(SAM_TLEN);
    unsigned add_bytes = IS_BAM ? sizeof (uint32_t) : tlen_len + 1;

    if (tlen) { // get tlen_value
        ASSSEG0 (tlen_len, tlen, "empty TLEN");

        bool is_int = str_get_int (tlen, tlen_len, &tlen_value); // note: tlen_value remains 0 if not a valid integer
        ASSSEG (is_int, tlen, "expecting TLEN to be an integer, but found \"%.*s\"", tlen_len, tlen);
    }

    if (segconf.running && tlen_value) segconf.has_TLEN_non_zero = true;

    int64_t predicted_tlen;
    if (segconf.has_TLEN_non_zero && sam_seg_predict_TLEN (vb, dl, is_rname_rnext_same, &predicted_tlen)) {

/*        if (predicted_tlen != tlen_value)
            printf ("WRONG FLAG=%x tlen=%d expected=%d : line_i=%u POS=%u PNEXT=%u RefConsumed=%u mate_RefConsumed=%u sup=%u next_unmapped=%u is_first=%u rev=%u next_rev=%u aligned=%u\n", 
                    dl->FLAG.value, (int)tlen_value, (int)predicted_tlen,
                    (int)vb->line_i, (int)dl->POS, (int)dl->PNEXT, vb->ref_consumed, (int)CTX(OPTION_MC_Z)->last_value.i,
                    dl->FLAG.bits.supplementary, dl->FLAG.bits.next_unmapped, dl->FLAG.bits.is_first, dl->FLAG.bits.rev_comp, dl->FLAG.bits.next_rev_comp, dl->FLAG.bits.is_aligned);
*/
        char tlen_by_calc[30] = { SNIP_SPECIAL, SAM_SPECIAL_TLEN, '0' + segconf.has[OPTION_MC_Z] };
        unsigned tlen_by_calc_len = 
            (tlen_value == predicted_tlen) ? 0 : str_int (tlen_value - predicted_tlen, &tlen_by_calc[3]);

        seg_by_ctx (VB, tlen_by_calc, tlen_by_calc_len + 3, ctx, add_bytes);
    }
    else if (tlen)
        seg_by_ctx (VB, STRa(tlen), ctx, add_bytes);

    else
        seg_integer_as_text_do (VB, ctx, tlen_value, add_bytes);

    ctx->last_value.i = tlen_value;
}

//---------
// PIZ
//---------

static inline int64_t sam_piz_predict_TLEN (VBlockSAM *vb, bool has_mc)
{
    SamFlags sam_flag = { .value = CTX(SAM_FLAG)->last_value.i };

    if (sam_flag.bits.supplementary || sam_flag.bits.next_unmapped) return 0;

    if (!sam_flag.bits.multi_segments) return vb->ref_consumed;

    const char *last_rname  = last_txt(vb, SAM_RNAME);
    const char *last_rnext  = last_txt(vb, SAM_RNEXT);
    unsigned last_rname_len = vb->last_txt_len(SAM_RNAME);
    unsigned last_rnext_len = vb->last_txt_len(SAM_RNEXT);
    
    if (flag.out_dt == DT_BAM && *(int32_t*)last_rname != *(int32_t*)last_rnext) return 0;

    if (flag.out_dt == DT_SAM && !(last_rnext_len==1 && *last_rnext=='=') &&
        (last_rname_len != last_rnext_len || memcmp (last_rname, last_rnext, last_rnext_len))) return 0;

    int64_t pnext_pos_delta = CTX(SAM_PNEXT)->last_value.i - CTX(SAM_POS)->last_value.i;

    if (!sam_flag.bits.rev_comp && sam_flag.bits.next_rev_comp) {

        if (has_mc) {
            STR(MC);
            reconstruct_peek (VB, CTX(OPTION_MC_Z), pSTRa(MC));
            unsigned MC_ref_consumed = sam_cigar_get_MC_ref_consumed (STRa(MC));
            return pnext_pos_delta + MC_ref_consumed;
        }
                
        else if (vb->buddy_line_i != NO_BUDDY) {
            uint32_t approx_mate_ref_consumed = *ENT (uint32_t, CTX(SAM_CIGAR)->piz_ctx_specific_buf, vb->buddy_line_i);

            return pnext_pos_delta + approx_mate_ref_consumed;
        }

        else
            return pnext_pos_delta + vb->ref_consumed;
    }

    if (sam_flag.bits.rev_comp && !sam_flag.bits.next_rev_comp) 
        return pnext_pos_delta - vb->ref_consumed + 2 * !sam_flag.bits.is_aligned;

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
SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_BUDDY_TLEN_old)
{
    ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line");

    new_value->i = -*ENT (int64_t, ctx->history, vb->buddy_line_i); // minus the buddy
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return true; // new value
}

// place value in correct location in alignment
TRANSLATOR_FUNC (sam_piz_sam2bam_TLEN)
{
    BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
    alignment->tlen = LTEN32 (ctx->last_value.i);
    return 0;
}
