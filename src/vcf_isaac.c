// ------------------------------------------------------------------
//   vcf_isaac.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "vcf_private.h"

void vcf_isaac_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, FORMAT_GQ, FORMAT_GQX, INFO_REFREP, DID_EOL); 

    // GQX stuff
    seg_mux_init (vb, FORMAT_GQX, VCF_SPECIAL_MUX_GQX, false, GQX); 
    seg_by_did (VB, STRa(vb->mux_GQX.snip), FORMAT_GQX, 0); // all the same

    ctx_set_dyn_int (VB, INFO_REFREP, DID_EOL);

    CTX(INFO_EVS)->con_rep_special = VCF_SPECIAL_N_ALTS;
}

void vcf_seg_FORMAT_GQX (VBlockVCFP vb, ContextP ctx, STRp(gqx))
{        
    Allele *ht_data = this_sample_GT (vb);
    
    int channel_i = (!ctx_encountered (VB, FORMAT_GT) || *ht_data == '.') ? 0
                  : !ctx_encountered (VB, FORMAT_GQ)                      ? 1
                  :                                                         2;

    ContextP chan_ctx = seg_mux_get_channel_ctx (VB, FORMAT_GQX, (MultiplexerP)&vb->mux_GQX, channel_i);
    
    if (channel_i == 2) 
        seg_delta_vs_other_localS (VB, chan_ctx, CTX(FORMAT_GQ), STRa(gqx), -1);

    else
        seg_integer_or_not (VB, chan_ctx, STRa(gqx), gqx_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_GQX)
{
    rom gt = last_txt (VB, FORMAT_GT);

    int channel_i = (!ctx_encountered (VB, FORMAT_GT) || *gt == '.') ? 0
                  : !ctx_encountered (VB, FORMAT_GQ)                 ? 1
                  :                                                    2;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

void vcf_seg_INFO_IDREP (VBlockVCFP vb, ContextP ctx, STRp(idrep_str))
{
    ContextP refrep_ctx = CTX(INFO_REFREP);
    int64_t refrep = refrep_ctx->last_value.i;
    int64_t idrep;

    if (ctx_has_value_in_line_(VB, refrep_ctx) && str_get_int (STRa(idrep_str), &idrep)) {
        SNIPi2 (SNIP_SPECIAL, VCF_SPECIAL_IDREP, (vb->ALT_len > vb->REF_len ? (idrep - refrep) : (refrep - idrep)));
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
    new_value->i = (vb->ALT_len > vb->REF_len) ? (refrep + delta) : (refrep - delta) ;

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

static bool vcf_seg_GMAF_allele_cb (VBlockP vb_, ContextP ctx, STRp(value), uint32_t unused_rep)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    bool has_refminor = ctx_encountered_in_line (VB, INFO_RefMinor);
    bool is_ref = str_is_1char (vb->REF, value[0]);
    bool is_alt = vb->alts[0][0] == value[0] && vb->alt_lens[0] == 1;

    if ((has_refminor && is_ref) || (!has_refminor && is_alt)) // usually case
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_GMAF_allele, '0' }, 3, ctx, value_len);
    
    else if ((!has_refminor && is_ref) || (has_refminor && is_alt)) // reverse case
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_GMAF_allele, '1' }, 3, ctx, value_len);

    else // can happen, for example, if ALT='.'
        seg_by_ctx (VB, STRa(value), ctx, value_len);

    return true; // segged successfully
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_GMAF_allele)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    
    // check INFO container prefixes for RefMinor (a valueless item, so no dict_id)
    #define prfx     vb->con_stack[vb->con_stack_len-2].prefixes // -1 is the GMAF container, -2 is the INFO container 
    #define prfx_len vb->con_stack[vb->con_stack_len-2].prefixes_len

    SAFE_NULT (prfx);
    bool has_refminor = !!strstr (prfx, CON_PX_SEP_ "RefMinor" CON_PX_SEP_);
    SAFE_RESTORE;

    // case: REF
    if ((has_refminor && *snip == '0') || (!has_refminor && *snip == '1'))
        RECONSTRUCT1 (*vb->REF);

    // case: ALT
    else
        RECONSTRUCT1 (*vb->alts[0]);

    return NO_NEW_VALUE;

    #undef prfx
}

static void vcf_GMAF_AF_prediction (STRp (af1000g_str), int pred_i, qSTRp (prediction))
{
    SAFE_NULT(af1000g_str);
    double af1000g = atof (af1000g_str);
    SAFE_RESTORE;

    *prediction_len = snprintf (prediction, 16, "%.4g", pred_i ? af1000g : (1 - af1000g));
}

static bool vcf_seg_GMAF_AF_cb (VBlockP vb_, ContextP ctx, STRp(value), uint32_t unused_rep)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    STRl (prediction,16);

    if (!has(AF1000G)) goto fallback;

    for (int pred_i=0; pred_i <= 1; pred_i++) {
        vcf_GMAF_AF_prediction (STRa(BII(AF1000G)->value), pred_i, qSTRa(prediction));

        if (str_issame (value, prediction)) {
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_GMAF_AF, '0' + pred_i }, 3, ctx, value_len);
            return true;
        }
    }

fallback:
    seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_SIMPLE, value_len);
    
    return true; // segged successfully
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_GMAF_AF)
{
    int pred_i = *snip - '0';

    STR(af1000g_str);
    reconstruct_peek (vb, CTX(INFO_AF1000G), pSTRa(af1000g_str));

    STRl (prediction, 16);
    vcf_GMAF_AF_prediction (STRa(af1000g_str), pred_i, qSTRa(prediction));

    RECONSTRUCT_str (prediction);

    return NO_NEW_VALUE;
}

// <ID=GMAF,Number=A,Type=String,Description="Global minor allele frequency (GMAF); technically, the frequency of the second most frequent allele.  Format: GlobalMinorAllele|AlleleFreqGlobalMinor">
// case 1: REF=A ALT=G GMAF=G|0.4822 
// case 2: REF=C ALT=. RefMinor;GMAF=C|0.04812
void vcf_seg_INFO_GMAF (VBlockVCFP vb, ContextP ctx, STRp(gmaf))
{        
    MediumContainer con = { .repeats   = 1,
                            .nitems_lo = 2,
                            .items[0]  = { .dict_id.num = DICT_ID_MAKE1_5("G0MAF"), .separator = "|" },
                            .items[1]  = { .dict_id.num = DICT_ID_MAKE1_5("G1MAF")   } };
    
    seg_struct (VB, ctx, con, STRa(gmaf), (SegCallback[]){ vcf_seg_GMAF_allele_cb, vcf_seg_GMAF_AF_cb }, gmaf_len, true);
}


void vcf_seg_INFO_EVS (VBlockVCFP vb, ContextP ctx, STRp(evs))
{        
    MediumContainer con = { .nitems_lo = 3,
                            .repsep    = ",",
                            .drop_final_repsep = true,
                            .items[0]  = { .dict_id.num = DICT_ID_MAKE1_7("E0VS_AF"),  .separator = "|" },
                            .items[1]  = { .dict_id.num = DICT_ID_MAKE1_8("E1VS_COV"), .separator = "|" },
                            .items[2]  = { .dict_id.num = DICT_ID_MAKE1_8("E2VS_CNT")                   } };

    seg_array_of_struct_ (VB, ctx, con, 0, 0, STRa(evs), 
                          (SegCallback[]){ seg_add_to_local_string_cb, seg_integer_or_not_cb, seg_integer_or_not_cb },
                          VCF_SPECIAL_N_ALTS, vb->n_alts, 0, evs_len);
}

