// ------------------------------------------------------------------
//   vcf_cosmic.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"

void vcf_cosmic_seg_initialize (VBlockVCFP vb)
{
    seg_id_field_init (CTX(INFO_LEGACY_ID));
}

void vcf_seg_INFO_LEGACY_ID (VBlockVCFP vb, ContextP ctx, STRp(lid))
{
    seg_id_field_do (VB, ctx, STRa(lid));
}

void vcf_seg_INFO_SO_TERM (VBlockVCFP vb, ContextP ctx, STRp(st))
{
    if ((vb->main_ref_len == 1 && vb->main_alt_len == 1 && str_issame_(st, st_len, "SNV"      , 3)) ||
        (vb->main_ref_len == 1 && vb->main_alt_len >  1 && str_issame_(st, st_len, "insertion", 9)) ||
        (vb->main_ref_len >  1 && vb->main_alt_len == 1 && str_issame_(st, st_len, "deletion" , 8)))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_SO_TERM }, 2, ctx, st_len);

    else
        seg_by_ctx (VB, STRa(st), ctx, st_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_SO_TERM)
{
    STRlast (refalt, VCF_REFALT);

    if (refalt_len == 3)        RECONSTRUCT ("SNV"      , 3);
    else if (refalt[1] == '\t') RECONSTRUCT ("insertion", 9);
    else                        RECONSTRUCT ("deletion" , 8);
    
    return NO_NEW_VALUE;
}
