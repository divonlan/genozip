// ------------------------------------------------------------------
//   vcf_gwas.c
//   Copyright (C) 2022-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

// ---------------
// ZIP
// ---------------

void vcf_gwas_zip_initialize (void)
{
}

void vcf_gwas_seg_initialize (VBlockVCFP vb)
{
    // even though sparse, ES, ES, LP and AF are rarely unique - we're slightly better off without singletons
    ctx_set_no_stons (VB, INFO_AF, FORMAT_ES, FORMAT_SE, FORMAT_LP, DID_EOL); 
}

//--------------------------------------------------------------------------------------------------------------
// VCF-GWAS: <ID=ID,Number=1,Type=String,Description="Study variant identifier">
// example: rs1885866
//--------------------------------------------------------------------------------------------------------------
void vcf_gwas_seg_FORMAT_ID (VBlockVCFP vb, ContextP ctx, STRp(id))
{
    if (ctx_encountered_in_line (VB, VCF_ID) && 
        CTX(VCF_ID)->last_txt.index != INVALID_LAST_TXT_INDEX &&
        str_issame_(STRa(id), STRtxt(CTX(VCF_ID)->last_txt))) 
       
        seg_by_ctx (VB, STRa(copy_VCF_ID_snip), ctx, id_len);

    else 
        seg_id_field (VB, ctx, STRa(id), false, id_len);
}

// ---------------
// PIZ
// ---------------
