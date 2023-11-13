// ------------------------------------------------------------------
//   sam_isaac.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "vcf_private.h"
#include "reconstruct.h"

void vcf_isaac_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, FORMAT_GQ, FORMAT_GQX, INFO_REFREP, DID_EOL); 

    // GQX stuff
    seg_mux_init (VB, CTX(FORMAT_GQX), 3, VCF_SPECIAL_MUX_GQX, false, (MultiplexerP)&vb->mux_GQX); 
    seg_by_did (VB, STRa(vb->mux_GQX.snip), FORMAT_GQX, 0); // all the same

    CTX(INFO_REFREP)->ltype = LT_DYN_INT;
}

void vcf_seg_FORMAT_GQX (VBlockVCFP vb, ContextP ctx, STRp(gqx))
{
    Allele *ht_data = B(Allele, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line + vb->ploidy * vb->sample_i);

    int channel_i = (!ctx_encountered (VB, FORMAT_GT) || *ht_data == '.') ? 0
                  : !ctx_encountered (VB, FORMAT_GQ)                      ? 1
                  :                                                         2;

    ContextP chan_ctx = seg_mux_get_channel_ctx (VB, FORMAT_GQX, (MultiplexerP)&vb->mux_GQX, channel_i);
    
    if (channel_i == 2) 
        seg_delta_vs_other (VB, chan_ctx, CTX(FORMAT_GQ), STRa(gqx));

    else
        seg_by_ctx (VB, STRa(gqx), chan_ctx, gqx_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_GQX)
{
    rom gt = last_txt (VB, FORMAT_GT);

    int channel_i = (!ctx_encountered (VB, FORMAT_GT) || *gt == '.') ? 0
                  : !ctx_encountered (VB, FORMAT_GQ)                 ? 1
                  :                                                    2;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs longer than 20 bases are not reported.
// example: REF=AC ALT=A RU=C
void vcf_seg_INFO_RU (VBlockVCFP vb, ContextP ctx, STRp(ru))
{
    if (*vb->main_ref == *vb->main_alt && vb->line_coords==DC_PRIMARY && 
        ((vb->main_alt_len > 1 && vb->main_ref_len == 1 && str_issame_(STRa(ru), vb->main_alt+1, vb->main_alt_len-1)) ||
         (vb->main_ref_len > 1 && vb->main_alt_len == 1 && str_issame_(STRa(ru), vb->main_ref+1, vb->main_ref_len-1))))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_RU }, 2, ctx, ru_len);

    else
        seg_by_ctx (VB, STRa(ru), ctx, ru_len); // possibly multiple ALT. TO DO: support multiple ALT
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_RU)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (reconstruct) {
        if (vb->main_ref_len > 1) 
            RECONSTRUCT (vb->main_ref+1, vb->main_ref_len-1);
        else
            RECONSTRUCT (vb->main_alt+1, vb->main_alt_len-1);
    }

    return NO_NEW_VALUE;
}

void vcf_seg_INFO_IDREP (VBlockVCFP vb, ContextP ctx, STRp(idrep_str))
{
    ContextP refrep_ctx = CTX(INFO_REFREP);
    int64_t refrep = refrep_ctx->last_value.i;
    int64_t idrep;

    if (ctx_has_value_in_line_(VB, refrep_ctx) && str_get_int (STRa(idrep_str), &idrep)) {
        SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_IDREP, (vb->main_alt_len > vb->main_ref_len ? (idrep - refrep) : (refrep - idrep)));
        seg_by_ctx (VB, STRa(snip), ctx, idrep_str_len);
    }
    else 
        seg_by_ctx (VB, STRa(idrep_str), ctx, idrep_str_len); // possibly multiple ALT. TO DO: support multiple ALT
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_IDREP)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    int64_t delta = atoi (snip);
    int64_t refrep = CTX(INFO_REFREP)->last_value.i;
    new_value->i = (vb->main_alt_len > vb->main_ref_len) ? (refrep + delta) : (refrep - delta) ;

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// INFO - mux by filter
int vcf_isaac_info_channel_i (VBlockP vb)
{
    STRlast (filter, VCF_FILTER);

    #define FILT(f) str_issame_(STRa(filter), f, STRLEN(f)) 
    return FILT("PASS") || FILT("LowGQX") || FILT("HighDPFRatio") || FILT("LowGQX;HighDPFRatio");
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_ISAAC_FILTER)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), vcf_isaac_info_channel_i (vb), new_value, reconstruct);
}
