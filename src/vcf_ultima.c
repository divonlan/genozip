// ------------------------------------------------------------------
//   vcf_ultima.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "reconstruct.h"
#include "reference.h"

// not 100% but close enough
static PosType64 X_LM_RM_predict_pos (VBlockVCFP vb, ConstRangeP range, PosType64 line_pos, bool is_rm, unsigned lm_seq_len)
{
    #define VT(x) (vb->var_types[0] == VT_##x)

    if ((VT(DEL) || VT(INS))) {
        int n_reps = 0;
        if (is_rm) {
            rom payload     = (vb->var_types[0] == VT_DEL)  ? &vb->main_ref[1] : &vb->main_alt[1];
            int payload_len = ((vb->var_types[0] == VT_DEL) ? vb->main_ref_len : vb->main_alt_len) - 1;

            PosType64 pos = line_pos + 1;
            int payload_i=0;
            decl_acgt_decode;
            while (REFp(pos) == payload[payload_i]) {
                pos++; 
                payload_i=(payload_i + 1) % payload_len;
                if (!payload_i) n_reps++; // counts full repeats
            }
            
            if ((VT(DEL) && !str_is_monochar (STRa(payload))) || (VT(INS) && !n_reps))
                pos -= payload_i; // exclude partial payload matches
            
            else if (VT(INS) && n_reps && payload_i) 
                pos = (pos - payload_i) + 1;

            return pos;
        }
        else
            return line_pos - (lm_seq_len - 1);
    }

    else if (VT(SNP)) {
        if (is_rm) 
            return line_pos + 1;
        else
            return line_pos - lm_seq_len;
    }

    return 0; // fail
} 

void vcf_seg_INFO_X_LM_RM (VBlockVCFP vb, ContextP ctx, STRp(seq))
{
    decl_acgt_decode;
    RefLock lock = REFLOCK_NONE;

    // to do: handle multi-ALT cases, eg "ATTTC,TATTT". Each item corresponds to one of the ALTs
    if (!flag.reference || !str_is_ACGT (STRa(seq), NULL)) 
        goto fallback;

    PosType64 line_pos = DATA_LINE(vb->line_i)->pos[0];

    RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), line_pos, seq_len, WORD_INDEX_NONE, 
                                      (IS_REF_EXT_STORE ? &lock : NULL));
    
    if (!range || line_pos < range->first_pos || line_pos + seq_len - 1 > range->last_pos)
        goto fallback;

    bool is_rm = (ctx->did_i == INFO_X_RM);
    PosType64 pos = X_LM_RM_predict_pos (vb, range, line_pos, is_rm, seq_len);

    // check again, with pos
    if (!pos || pos < range->first_pos || pos + seq_len - 1 > range->last_pos)
        goto fallback;

    uint32_t index_within_range = pos - range->first_pos;

    // find delta
    int delta = 999;
    for (int delta_i=0; delta_i < 200; delta_i++) {
        int candidate_delta = DEINTERLACE(typeof(delta_i), delta_i); // 0, -1, 1, -2, 2...
        
        bool found = true; // optimistic
        for (int i=0; i < seq_len; i++)
            if (seq[i] != REFp (pos + i + candidate_delta)) {
                found = false;
                break;
            }
        
        if (found) {
            delta = candidate_delta;
            break;
        }
    }
    if (delta == 999) 
        goto fallback;
    
    SNIPi2_2 (SNIP_SPECIAL, VCF_SPECIAL_X_LM_RM, seq_len, delta);
    seg_by_ctx (VB, STRa(snip), ctx, seq_len);

    if (IS_REF_EXT_STORE)
        bits_set_region (&range->is_set, index_within_range, seq_len);

    ref_unlock (gref, &lock); 
    return;

fallback:
    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
    seg_by_ctx (VB, STRa(seq), ctx, seq_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_X_LM_RM)
{    
    if (!reconstruct) return NO_NEW_VALUE;

    str_split_ints (snip, snip_len, 2, ',', item, true);
    int seq_len = items[0];
    int delta   = items[1];

    ConstRangeP range = ref_piz_get_range (vb, gref, HARD_FAIL);

    PosType64 pos = delta + X_LM_RM_predict_pos (VB_VCF, range, CTX(VCF_POS)->last_value.i, (ctx->did_i == INFO_X_RM), seq_len);
    
    char *next = BAFTtxt;
    decl_acgt_decode;
    
    for (uint32_t i=0; i < seq_len; i++)
        *next++ = REFp (pos++);

    Ltxt += seq_len;

    return NO_NEW_VALUE;
}
