// ------------------------------------------------------------------
//   vcf_1000G.c : 1000 Genome Project fields
//   Copyright (C) 2024-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

sSTRl(half_AN_snip, 32);

void vcf_1000G_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_snip_special_other_int (VCF_SPECIAL_DIVIDE_BY, _INFO_AN, half_AN_snip, 2);
    }
}

void vcf_1000G_seg_initialize (VBlockVCFP vb)
{
    ctx_set_same_line (VB, INFO_MAF, INFO_NS, DID_EOL);
}

void vcf_seg_INFO_MAF (VBlockVCFP vb, ContextP ctx, STRp(maf))
{
    if (has(AF) && str_issame_(STRa(maf), STRa(BII(AF)->value)))
        seg_by_ctx (VB, STRa(copy_INFO_AF_snip), ctx, maf_len);
    
    else
        seg_by_ctx (VB, STRa(maf), ctx, maf_len);
}

void vcf_seg_INFO_NS (VBlockVCFP vb, ContextP ctx, STRp(ns_str))
{
    // TO DO: a better method would be to postpone NS to after samples and actually count samples with data
    int64_t an, ns;

    if (has(AN) && str_get_int (STRa(BII(AN)->value), &an) && 
        str_get_int (STRa(ns_str), &ns) && ns == an / 2)

        seg_by_ctx (VB, STRa(half_AN_snip), ctx, ns_str_len);
    
    else
        seg_integer_or_not (VB, ctx, STRa(ns_str), ns_str_len);
}
