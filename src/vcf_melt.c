// ------------------------------------------------------------------
//   vcf_melt.c
//   Copyright (C) 2024-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

sSTRl(END_minus_SVLEN_snip, 32);
sSTRl(copy_END_snip, 32);

void vcf_me_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_minus_snip (VCF, _VCF_POS, _INFO_SVLEN, END_minus_SVLEN_snip);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_POS, true, 0, copy_END_snip); // END is an alias of POS
    }
}

void vcf_seg_melt_ADJLEFT (VBlockVCFP vb, ContextP ctx, STRp(adjleft_str))
{
    int64_t svlen = CTX(INFO_SVLEN)->last_value.i;
    int64_t end   = CTX(VCF_POS)->last_value.i; // END is an alias of POS
    int64_t adjleft;

    if (ctx_encountered_in_line (VB, INFO_SVLEN)  &&
        ctx_encountered_in_line (VB, INFO_END)    &&
        str_get_int (STRa(adjleft_str), &adjleft) &&
        adjleft == end - svlen)

        seg_by_ctx (VB, STRa(END_minus_SVLEN_snip), ctx, adjleft_str_len);
    
    else if (str_issame_(STRa(adjleft_str), "0", 1))
        seg_by_ctx (VB, "0", 1, ctx, 1); // seg as snip

    else
        seg_integer_or_not (VB, ctx, STRa(adjleft_str), adjleft_str_len);
}

void vcf_seg_melt_ADJRIGHT (VBlockVCFP vb, ContextP ctx, STRp(adjright_str))
{
    int64_t end = CTX(VCF_POS)->last_value.i; // END is an alias of POS
    int64_t adjright;

    if (ctx_encountered_in_line (VB, INFO_END)    &&
        str_get_int (STRa(adjright_str), &adjright) &&
        adjright == end)

        seg_by_ctx (VB, STRa(copy_END_snip), ctx, adjright_str_len);

    else if (str_issame_(STRa(adjright_str), "0", 1))
        seg_by_ctx (VB, "0", 1, ctx, 1); // seg as snip

    else
        seg_integer_or_not (VB, ctx, STRa(adjright_str), adjright_str_len);
}
