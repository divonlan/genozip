// ------------------------------------------------------------------
//   vcf_info_QD.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "strings.h"
#include "dict_id.h"
#include "reconstruct.h"

//------------
// Seg
//------------

// add up DP's of samples with GT!=0/0, for consumption by INFO/QD predictor
void vcf_seg_sum_DP_for_QD (VBlockVCFP vb, int64_t value)
{
    if (ctx_encountered_in_line (VB, INFO_QD) && 
        ctx_has_value (VB, FORMAT_GT) && CTX(FORMAT_GT)->last_value.i >= 1 /*dosage*/)
        
        CTX(INFO_QD)->qd.sum_dp_with_dosage += value;
}

// we need a two-decimal-digit number with no trailing zeros
static bool vcf_seg_QD_verify_prediction (ContextP ctx, double qual_value, double sum_dp, STRp(qd), int add)
{
    STRl(qd_pred, 32);

    bool has_2_decimals = (qd_len >= 3 && qd[qd_len-3] == '.');
    bool has_trailing_zeros = (has_2_decimals && qd[qd_len-1] == '0');

    if (ctx->is_initialized && ((!ctx->flags.trailing_zero_ok && has_trailing_zeros) || // has trailing zero despite policy saying no
                                (ctx->flags.trailing_zero_ok && !has_2_decimals)))      // not exactly 2 decimals despite mandating trailing zeros
        return false;

    int pred = 100.0 * qual_value / sum_dp + add + 0.5; // note: if (100.0 * qual_value / sum_dp) decimal is exactly 0.xx5, it is rounded inconsistently in the data - some times up and sometimes down
    qd_pred_len = sprintf (qd_pred, "%.2f", (double)pred / 100.0); // with trailing zeros if number is round

    bool pred_has_trailing_zeros = (qd_pred[qd_pred_len-1] == '0');

    // remove trailing zeros if we are obliged to
    if (ctx->is_initialized && !ctx->flags.trailing_zero_ok && pred_has_trailing_zeros)
        qd_pred_len -= (qd_pred[qd_pred_len-2] != '0') ? 1 : 3; // if both are 0, remove decimal point too

    // case: prediction successful. if round number - trailing zeros according to policy, or if no policy - with zeros
    if (str_issame (qd, qd_pred)) {
        if (!ctx->is_initialized && pred_has_trailing_zeros) { // we encounted the first example of trailing zeros - set the policy
            ctx->is_initialized = true;
            ctx->flags.trailing_zero_ok = true;
        }
        return true;
    }

    // case: prediction failed. if round number - trailing zeros according to policy, or if no policy - with zeros 
    else {
        // case: we don't have a policy yet try again after removing trailing zeros
        if (!ctx->is_initialized && !has_2_decimals && pred_has_trailing_zeros) {
            qd_pred_len -= (qd_pred[qd_pred_len-2] != '0') ? 1 : 3; // if both are 0, remove decimal point too

            // we encounted the first example of "no trailing zeros" - set the policy
            if (str_issame (qd, qd_pred)) {
                ctx->is_initialized = true;
                ctx->flags.trailing_zero_ok = false;
                return true; 
            }
        }

        return false; // prediction failed
    }
}

static QdPredType vcf_seg_is_QD_predictable (VBlockVCFP vb, ContextP ctx, STRp(qd))
{
    STRlast (qual, VCF_QUAL);

    double qual_value, qd_value;
    if (!str_get_float (STRa(qual), &qual_value, 0, 0) || !str_get_float (STRa(qd), &qd_value, 0, 0) || qd_value <= 0.0) 
        return QD_PRED_NONE;
    
    int ratio = (int)(qual_value / qd_value + 0.5);

    bool has_info_dp = ctx_has_value_in_line_(vb, CTX(INFO_DP));
    int info_dp = CTX(INFO_DP)->last_value.i;

    // if single sample and we have INFO/DP, we use INFO/DP instead of the sample DPs
    if (has_info_dp && info_dp >= ratio-1 && info_dp <= ratio+1) {
        if (vcf_seg_QD_verify_prediction (ctx, qual_value, info_dp, STRa(qd), 0))  return QD_PRED_INFO_DP;
        if (vcf_seg_QD_verify_prediction (ctx, qual_value, info_dp, STRa(qd), 1))  return QD_PRED_INFO_DP_P001;
        if (vcf_seg_QD_verify_prediction (ctx, qual_value, info_dp, STRa(qd), -1)) return QD_PRED_INFO_DP_M001;
    }

    // prediction based on sum of FORMAT/DP, excluding samples with 0/0 or ./.
    if ((vcf_num_samples + !has_info_dp >= 2) && !LO_IS_OK_SWITCH (last_ostatus) &&  // can't use this in case of a REF⇆ALT switch, because GT changes)
         ctx->qd.sum_dp_with_dosage >= ratio-1 && ctx->qd.sum_dp_with_dosage <= ratio+1) {
        if (vcf_seg_QD_verify_prediction (ctx, qual_value, ctx->qd.sum_dp_with_dosage, STRa(qd), 0 )) return QD_PRED_SUM_DP;
        if (vcf_seg_QD_verify_prediction (ctx, qual_value, ctx->qd.sum_dp_with_dosage, STRa(qd), 1 )) return QD_PRED_SUM_DP_P001;
        if (vcf_seg_QD_verify_prediction (ctx, qual_value, ctx->qd.sum_dp_with_dosage, STRa(qd), -1)) return QD_PRED_SUM_DP_M001;
    }
    
    // case: prediction failed
    return QD_PRED_NONE;
}

// called after all samples have been segged, and sum_dp_with_dosage is set
void vcf_seg_INFO_QD (VBlockP vb)
{
    decl_ctx (INFO_QD);
    STRlast (qd, INFO_QD);
    
    SEGCONF_RECORD_WIDTH (QD, qd_len);

    // case: we can't generate a prediction or prediction is wrong - seg normally
    QdPredType pd;
    if (!(pd = vcf_seg_is_QD_predictable (VB_VCF, ctx, STRa(qd)))) 
        seg_by_ctx (vb, STRa(qd), ctx, qd_len);

    else
        seg_by_ctx (vb, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_QD, '0' + pd }, 3, ctx, qd_len);
}

//------------
// Piz
//------------

// called from toplevel callback
void vcf_piz_insert_INFO_QD (VBlockVCFP vb)
{
    decl_ctx (INFO_QD);
    if (!ctx->recon_insertion) return;

    QdPredType type = ctx->qd.pred_type;

    double qual_value = CTX(VCF_QUAL)->last_value.f;
    STRl(qd, 32); // prediction

    uint32_t sum_dp = (type == QD_PRED_INFO_DP || type == QD_PRED_INFO_DP_P001 || type == QD_PRED_INFO_DP_M001) 
                        ? CTX(INFO_DP)->last_value.i : ctx->qd.sum_dp_with_dosage;

    int add = (type == QD_PRED_INFO_DP_P001 || type == QD_PRED_SUM_DP_P001) ? 1
            : (type == QD_PRED_INFO_DP_M001 || type == QD_PRED_SUM_DP_M001) ? -1
            :                                                                 0;


    int pred = 100.0 * qual_value / sum_dp + 0.5 + add; 
    qd_len = sprintf (qd, "%.2f", (double)pred / 100.0); // with trailing zeros if number is round

    // remove trailing zeros if we are obliged to
    if (qd[qd_len-1] == '0' && !ctx->flags.trailing_zero_ok)
        qd_len -= (qd[qd_len-2] != '0') ? 1 : 3; // if both are 0, remove decimal point too

    vcf_piz_insert_field (vb, ctx, STRa(qd), segconf.wid_QD.width);
}

void vcf_piz_sum_DP_for_QD (VBlockP vb, STRp(recon))
{
    int64_t dp;
    if (vcf_piz_GT_get_last_dosage (vb) >= 1 && str_get_int (STRa(recon), &dp))
        CTX(INFO_QD)->qd.sum_dp_with_dosage += dp;
}

// just store pred_type - actual reconstruction is deferred to vcf_piz_reconstruct_QD
SPECIAL_RECONSTRUCTOR (vcf_piz_special_QD)
{
    ctx->qd.pred_type = (uint32_t)(snip[0] - '0');

    ASSPIZ (ctx->qd.pred_type < NUM_QD_PRED_TYPES, 
            "Unknown pred_type=%d. %s", ctx->qd.pred_type, genozip_update_msg());

    vcf_piz_defer_to_after_samples (QD);

    return NO_NEW_VALUE;
}
