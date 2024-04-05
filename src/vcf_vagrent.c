// ------------------------------------------------------------------
//   vcf_vagrent.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"

sSTRl(copy_VD_snip, 30);

void vcf_vagrent_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _INFO_VD, false, 0, copy_VD_snip);
    }
}

static void vcf_seg_VD_VW (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    // examples of values:
    // ENSG00000197049|ENST00000358533|downstream
    // SAMD11|CCDS2.2|intronic
    // ISG15|CCDS6.1|r.44_45insa|c.?|p.?
    // CDK11B|CCDS44042.1|r.452_453insaaagaa|c.372_373insAAAGAA|p.E126_R127insKE

    // TO DO: improve this 
    seg_array (VB, ctx, ctx->did_i, STRa(value), '|', 0, false, STORE_NONE, _INFO_VDVW_ARR, value_len);
}

// ##INFO=<ID=VD,Number=1,Type=String,Description="Vagrent Default Annotation">
void vcf_seg_INFO_VD (VBlockVCFP vb, ContextP ctx, STRp(vd))
{
    vcf_seg_VD_VW (vb, ctx, STRa(vd));

    set_last_txt (INFO_VD, vd);
}

// ##INFO=<ID=VW,Number=1,Type=String,Description="Vagrent Most Deleterious Annotation">
void vcf_seg_INFO_VW (VBlockVCFP vb, ContextP ctx, STRp(vw))
{
    if (ctx_encountered_in_line (VB, INFO_VD) && str_issame_(STRa(vw), STRlst(INFO_VD)))
        seg_by_ctx (VB, STRa(copy_VD_snip), ctx, vw_len);

    else
        vcf_seg_VD_VW (vb, ctx, STRa(vw));
}
