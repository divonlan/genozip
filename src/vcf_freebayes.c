// ------------------------------------------------------------------
//   vcf_freebayes.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "zip_dyn_int.h"

sSTRl(dp_minus_ro_snip,32);

void vcf_freebayes_zip_initialize (void)
{
    seg_prepare_minus_snip (VCF, _FORMAT_DP, _FORMAT_RO, dp_minus_ro_snip);
}

void vcf_freebayes_seg_initialize (VBlockVCFP vb)
{
    if (segconf.FMT_RO_AO_method == RO_AO_by_AD) {
        CTX(FORMAT_RO)->other_ctx = ctx_get_ctx (vb, make_array_item_dict_id (_FORMAT_AD, 0));
        CTX(FORMAT_AO)->other_ctx = ctx_get_ctx (vb, make_array_item_dict_id (_FORMAT_AD, 1));

        ctx_set_store (VB, STORE_INT, CTX(FORMAT_RO)->other_ctx->did_i, CTX(FORMAT_AO)->other_ctx->did_i, DID_EOL);
    }

    ctx_set_store (VB, STORE_INT, FORMAT_RO, FORMAT_AO, INFO_DPB, DID_EOL);
}

void vcf_seg_FORMAT_RO_AO (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    if (segconf.FMT_RO_AO_method == RO_AO_by_AD) {
        if (ctx_has_value (VB, ctx->other_ctx->did_i)) 
            seg_delta_vs_other_dictS (VB, ctx, ctx->other_ctx, STRa(value), -1);
        else 
            seg_by_ctx (VB, STRa(value), ctx, value_len); // usually "."
    }    

    else if (segconf.FMT_RO_AO_method == RO_AO_by_DP) {
        int64_t ao;
        // predicting: AO + RO = DP
        if ((ctx->did_i == FORMAT_AO) && ctx_has_value (VB, FORMAT_DP) && ctx_has_value (VB, FORMAT_RO) &&
            str_get_int (STRa(value), &ao) && ao == CTX(FORMAT_DP)->last_value.i - CTX(FORMAT_RO)->last_value.i) {
            
            seg_by_ctx (VB, STRa(dp_minus_ro_snip), ctx, value_len);
        }

        else 
            seg_integer_or_not (VB, ctx, STRa(value), value_len);
    } 
    
    set_last_txt (ctx->did_i, value);
}

void vcf_seg_FORMAT_QR_QA (VBlockVCFP vb, ContextP ctx, STRp(value_str))
{
    ContextP other_ctx = (ctx->did_i == FORMAT_QA) ? CTX(FORMAT_AO) : CTX(FORMAT_RO);

    int64_t value;

    // value is '.' as predicted
    if (str_is_1char (value_str, '.') && 
        other_ctx->last_txt.len==1 && *last_txtx (vb, other_ctx) == '.') {
        
        seg_special0 (VB, VCF_SPECIAL_QR_QA, ctx, value_str_len); 
        return;
    }

    else if (!str_get_int (STRa(value_str), &value)) fallback:
        seg_by_ctx (VB, STRa(value_str), ctx, value_str_len);

    else {
        if (segconf_running) {
            if (ctx_has_value (VB, other_ctx->did_i) && other_ctx->last_value.i) {
                segconf.Q_to_O += (float)value / (float)other_ctx->last_value.i;
                segconf.n_Q_to_O++;
            }
            goto fallback;
        }

        if (!ctx_has_value (VB, other_ctx->did_i) || !segconf.Q_to_O) goto fallback;

        int64_t prediction = round ((float)other_ctx->last_value.i * (float)segconf.Q_to_O);
        int64_t delta = value - prediction;

        seg_special0 (VB, VCF_SPECIAL_QR_QA, ctx, value_str_len);
        dyn_int_append (VB, ctx, delta, 0);
    }
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_QR_QA)
{
    ContextP other_ctx = (ctx->did_i == FORMAT_QA) ? CTX(FORMAT_AO) : CTX(FORMAT_RO);

    if (other_ctx->last_txt.len==1 && *last_txtx (vb, other_ctx) == '.') {
        RECONSTRUCT1 ('.');
        return NO_NEW_VALUE;
    }
    int64_t prediction = round ((float)other_ctx->last_value.i * (float)segconf.Q_to_O);

    int64_t delta = reconstruct_from_local_int (VB, ctx, 0, RECON_OFF);

    new_value->i = prediction + delta;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);  

    return HAS_NEW_VALUE;
}

void vcf_seg_INFO_DPB (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    decl_ctx(INFO_DPB);
    STRlast (dpb_str, INFO_DPB);
 
    SEGCONF_RECORD_WIDTH (INFO_DPB, dpb_str_len);

    int64_t dpb;    
    if (ctx_has_value (VB, INFO_DP) && str_get_int (STRa(dpb_str), &dpb)) {
        int64_t delta = dpb - CTX(INFO_DP)->last_value.i; 

        // prepare snip: VCF_SPECIAL_DEFER followed by snip to be reconstructed after samples
        SNIP(32);
        snip[0] = SNIP_SPECIAL;
        snip[1] = VCF_SPECIAL_DEFER;
        snip_len -= 2;
        seg_prepare_snip_other_do (SNIP_OTHER_DELTA, _INFO_DP, true, delta, 0, snip+2, &snip_len);
    
        seg_by_ctx (VB, snip, snip_len+2, ctx, dpb_str_len);
        // seg_delta_vs_other_dictN (VB, ctx, CTX(INFO_DP), dpb, -1, dpb_str_len); // store in dict, because float values go to local
    }

    else
        seg_add_to_local_string (VB, ctx, STRa(dpb_str), LOOKUP_SIMPLE, dpb_str_len); // usually floats, but possibly also other values, eg if DP is invalid
}

