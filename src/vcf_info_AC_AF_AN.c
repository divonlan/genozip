// ------------------------------------------------------------------
//   vcf_ac_af_an.c
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"

// -------
// INFO/AC
// -------

static bool vcf_seg_INFO_AC_item (VBlockP vb_, ContextP arr_ctx, STRp(ac_item), uint32_t alt_i)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    bool is_mle = (arr_ctx->did_i == INFO_MLEAC || arr_ctx->st_did_i == INFO_MLEAC); // MLEAC M0LEAC etc
    
    double af;
    int64_t ac;
    if (str_get_int (STRa(ac_item), &ac) &&
        ( (!is_mle && str_item_i_float (STRa(BII(AF)->value), ',', alt_i, &af)) ||
          ( is_mle && str_item_i_float (STRa(BII(MLEAF)->value), ',', alt_i, &af)))) {
        
        int delta = ac - (int64_t)round (af * CTX(INFO_AN)->last_value.i); // last_value was set in vcf_seg_INFO_AC

        if (delta) { 
            SNIPi3 (SNIP_SPECIAL, VCF_SPECIAL_INFO_AC, is_mle ? 'M' : 'D', delta);
            seg_by_ctx (VB, STRa(snip), arr_ctx, ac_item_len); 
        }
        else
            seg_special1 (VB, VCF_SPECIAL_INFO_AC, is_mle ? 'M' : 'D', arr_ctx, ac_item_len);
    }

    else
        seg_by_ctx (VB, STRa(ac_item), arr_ctx, ac_item_len);
        
    return true;
}

void vcf_seg_INFO_AC (VBlockVCFP vb, ContextP ctx, STRp(ac_str))
{
    SEGCONF_RECORD_WIDTH (INFO_AC, ac_str_len);

    int64_t ac, ac_hom, ac_het, ac_hemi;

    // case: expecting AC ≅ AN * AF for each alt allele (might not be exact due to rounding errors, esp if AF is a very small fraction)
    if (has(AN) && has(AF) && str_get_int (STRa(BII(AN)->value), &CTX(INFO_AN)->last_value.i)) 
        seg_array_by_callback (VB, ctx, STRa(ac_str), ',', vcf_seg_INFO_AC_item, VCF_SPECIAL_N_ALTS, N_ALTS, ac_str_len); // since 15.0.36, we seg AC items as an array

    // GIAB: AC = AC_Hom + AC_Het + AC_Hemo
    else if (has(AC_Hom) && has(AC_Het) && has(AC_Hemi) &&
             str_get_int (STRa(ac_str), &ac) &&
             str_get_int (STRa(BII(AC_Hom)->value),  &ac_hom) && 
             str_get_int (STRa(BII(AC_Het)->value),  &ac_het) && 
             str_get_int (STRa(BII(AC_Hemi)->value), &ac_hemi) &&
             ac == ac_hom + ac_het + ac_hemi) { 
             
        seg_special1 (VB, VCF_SPECIAL_INFO_AC, '1', ctx, ac_str_len);
    }

    // fallback
    else 
        seg_by_ctx (VB, STRa(ac_str), ctx, ac_str_len);
}

void vcf_seg_INFO_MLEAC (VBlockVCFP vb, ContextP ctx, STRp(mleac))
{
    SEGCONF_RECORD_WIDTH (INFO_MLEAC, mleac_len);

    // case: expecting MLEAC ≅ AN * MLEAF for each alt allele (might not be exact due to rounding errors, esp if AF is a very small fraction)
    if (has(AN) && has(MLEAF) && str_get_int (STRa(BII(AN)->value), &CTX(INFO_AN)->last_value.i)) 
        seg_array_by_callback (VB, ctx, STRa(mleac), ',', vcf_seg_INFO_AC_item, VCF_SPECIAL_N_ALTS, N_ALTS, mleac_len); 

    // fallback
    else 
        seg_by_ctx (VB, STRa(mleac), ctx, mleac_len);
}

// reconstruct: AC = AN * AF 
SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_AC)
{
    // note: update last_value too, so its available to vcf_piz_luft_A_AN, which is called becore last_value is updated

    // case: called for item in AC / MLEAC array, possibly with a delta (since 15.0.36)
    if (snip_len && (*snip == 'D' || *snip == 'M')) {
        STR(af_str);
        reconstruct_peek (vb, CTX((*snip == 'D') ? INFO_AF : INFO_MLEAF), pSTRa(af_str));

        uint32_t alt_i = current_con.repeat;
        
        double af;
        ASSPIZ (str_item_i_float (STRa(af_str), ',', alt_i, &af), "Failed get AF of alt_i=%u", alt_i);

        int32_t delta = (snip_len > 1) ? atoi (&snip[1]) : 0;

        ctx->last_value.i = new_value->i = delta + (int64_t)round (reconstruct_peek (vb, CTX(INFO_AN), 0, 0).i * af);
    }

    // Backward compatability note: In files v6->11, snip has 2 bytes for AN, AF which mean: '0'=appears after AC, '1'=appears before AC. We ignore them.
    // used in files up to 15.0.35
    else if (!snip_len || !VER(15))
        ctx->last_value.i = new_value->i = (int64_t)round (reconstruct_peek (vb, CTX(INFO_AN), 0, 0).i * reconstruct_peek(vb, CTX(INFO_AF), 0, 0).f);

    // case: GIAB
    else if (*snip == '1') 
        ctx->last_value.i = new_value->i = 
            reconstruct_peek (vb, CTX(INFO_AC_Hom),  0, 0).i +
            reconstruct_peek (vb, CTX(INFO_AC_Het),  0, 0).i +
            reconstruct_peek (vb, CTX(INFO_AC_Hemi), 0, 0).i;

    else 
        ABORT_PIZ ("unrecognized snip '%c'(%u). %s", *snip, (uint8_t)*snip, genozip_update_msg());

    if (reconstruct) RECONSTRUCT_INT (new_value->i); 

    return HAS_NEW_VALUE;
}

static bool vcf_seg_INFO_MLEAF_item (VBlockP vb_, ContextP arr_ctx, STRp(mleaf_item), uint32_t alt_i)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    STR(af_item);
    str_item_i (STRa(BII(AF)->value), ',', alt_i, pSTRa(af_item));

    if (str_issame (af_item, mleaf_item))
        seg_special0 (VB, VCF_SPECIAL_INFO_MLEAF, arr_ctx, mleaf_item_len);

    else
        seg_by_ctx (VB, STRa(mleaf_item), arr_ctx, mleaf_item_len);
        
    return true;
}

void vcf_seg_INFO_MLEAF (VBlockVCFP vb, ContextP ctx, STRp(mleaf))
{
    // case: in many cases MLEAF == AF for any alt allele 
    if (has(AF)) 
        seg_array_by_callback (VB, ctx, STRa(mleaf), ',', vcf_seg_INFO_MLEAF_item, VCF_SPECIAL_N_ALTS, N_ALTS, mleaf_len); 

    // fallback
    else 
        seg_by_ctx (VB, STRa(mleaf), ctx, mleaf_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_MLEAF)
{
    STR(af_str);
    reconstruct_peek (vb, CTX(INFO_AF), pSTRa(af_str));

    uint32_t alt_i = current_con.repeat;
    
    STR(af_item);
    ASSPIZ (str_item_i (STRa(af_str), ',', alt_i, pSTRa(af_item)), "Failed get AF of alt_i=%u", alt_i);

    RECONSTRUCT_str (af_item);

    return NO_NEW_VALUE;
}

