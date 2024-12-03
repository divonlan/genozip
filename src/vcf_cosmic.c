// ------------------------------------------------------------------
//   vcf_cosmic.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

void vcf_cosmic_seg_initialize (VBlockVCFP vb)
{
}

// eg: "COSN8882887", "COSN30315184"
void vcf_seg_INFO_LEGACY_ID (VBlockVCFP vb, ContextP ctx, STRp(lid))
{
    seg_id_field (VB, ctx, STRa(lid), false, lid_len);
}

void vcf_seg_INFO_SO_TERM (VBlockVCFP vb, ContextP ctx, STRp(st))
{
    if ((vb->REF_len == 1 && vb->ALT_len == 1 && str_issame_(st, st_len, "SNV"      , 3)) ||
        (vb->REF_len == 1 && vb->ALT_len >  1 && str_issame_(st, st_len, "insertion", 9)) ||
        (vb->REF_len >  1 && vb->ALT_len == 1 && str_issame_(st, st_len, "deletion" , 8)))
        seg_special0 (VB, VCF_SPECIAL_SO_TERM, ctx, st_len);

    else
        seg_by_ctx (VB, STRa(st), ctx, st_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_SO_TERM)
{
    if (VT0(SNP))      RECONSTRUCT ("SNV"      , 3);
    else if (VT0(INS)) RECONSTRUCT ("insertion", 9);
    else                RECONSTRUCT ("deletion" , 8);
    
    return NO_NEW_VALUE;
}
