// ------------------------------------------------------------------
//   vcf_dbsnp.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

sSTRl(delta_ID_snip, 32);
sSTRl(delta_POS_snip, 32);

void vcf_dbsnp_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_ID,  true, 0, delta_ID_snip);   // copies numeric part of ID, expecting delta to be 0
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_POS, true, 0, delta_POS_snip);
    }
}

void vcf_dbsnp_seg_initialize (VBlockVCFP vb)
{
    seg_mux_init (vb, INFO_VC, VCF_SPECIAL_MUX_BY_VARTYPE, true, VC);

    ctx_set_no_stons (VB, VCF_QUAL, INFO_FREQ, DID_EOL);
}

// ##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
// *might* be the same as the numeric value of the ID
void vcf_seg_INFO_RS (VBlockVCFP vb, ContextP ctx, STRp(rs))
{
    // case: eg ID=rs3844233 RS=3844233. We use a SNIP_OTHER_DELTA with delta=0 to copy last_value from ID.
    // (can't use SNIP_COPY bc it would copy the entire txt "rs3844233")
    int64_t rs_value;
    if (ctx_has_value_in_line_(VB, CTX(VCF_ID)) &&
        str_get_int (STRa(rs), &rs_value) &&
        rs_value == CTX(VCF_ID)->last_value.i) {

        seg_by_ctx (VB, STRa(delta_ID_snip), ctx, rs_len);
    }
    
    else 
        seg_integer_or_not (VB, ctx, STRa(rs), rs_len);
}

// ##INFO=<ID=RSPOS,Number=1,Type=Integer,Description="Chr position reported in dbSNP">
// *might* be the same as POS
void vcf_seg_INFO_RSPOS (VBlockVCFP vb, ContextP ctx, STRp(rspos))
{
    // case: eg ID=rs3844233 RS=3844233. We use a SNIP_OTHER_DELTA with delta=0 to copy last_value from ID.
    // (can't use SNIP_COPY bc it would copy the entire txt "rs3844233")
    int64_t rspos_value;
    if (ctx_has_value_in_line_(VB, CTX(VCF_POS)) &&
        str_get_int (STRa(rspos), &rspos_value)) {

        if (rspos_value == CTX(VCF_POS)->last_value.i) // shortcut for most common case
            seg_by_ctx (VB, STRa(delta_POS_snip), ctx, rspos_len);

        else {
            STRlic(snip, 32);
            seg_prepare_snip_other (SNIP_OTHER_DELTA, _VCF_POS, true, rspos_value - CTX(VCF_POS)->last_value.i, snip);
            seg_by_ctx (VB, STRa(snip), ctx, rspos_len);            
        }
    }
    
    else 
        seg_integer_or_not (VB, ctx, STRa(rspos), rspos_len);
}

typedef enum { VARTYPE_SNP=0, 
               VARTYPE_NONSNP_1REF=1,    // non-SNP, single-base REF - usually an insertion
               VARTYPE_MREF=2 } VarType; // multi-base REF - usually a deletion

static VarType get_vartype (VBlockVCFP vb)
{
    if (vb->REF_len > 1) return VARTYPE_MREF;

    else if (vb->ALT_len == 1 || vb->ALT_len == str_count_char (STRa(vb->ALT), ',') * 2 + 1) return VARTYPE_SNP;

    else return VARTYPE_NONSNP_1REF;
}

// <ID=VC,Number=1,Type=String,Description="Variation Class">
// Possible values include: SNV, INS, INDEL...
void vcf_seg_INFO_VC (VBlockVCFP vb, ContextP ctx, STRp(vc))
{
    if (!segconf_running) {
        VarType vartype = get_vartype (vb);
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, INFO_VC, (MultiplexerP)&vb->mux_VC, vartype);

        seg_by_ctx (VB, STRa(vc), channel_ctx, vc_len);
        seg_by_ctx (VB, STRa(vb->mux_VC.snip), ctx, 0);
    }

    else
        seg_by_ctx (VB, STRa(vc), ctx, vc_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_VARTYPE)
{
    VarType vartype = get_vartype (VB_VCF);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), vartype, new_value, reconstruct);
}
