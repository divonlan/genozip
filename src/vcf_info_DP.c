// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "zip_dyn_int.h"

//--------
// INFO/DP
// -------

static void vcf_seg_INFO_DP_by_FORMAT_DP (VBlockP vb); // forward
static void vcf_seg_INFO_DP_by_BaseCounts (VBlockP vb);

// return true if caller still needs to seg 
void vcf_seg_INFO_DP (VBlockVCFP vb, ContextP ctx, STRp(dp_str))
{
    SEGCONF_RECORD_WIDTH (INFO_DP, dp_str_len);

    // used in: vcf_seg_one_sample (for 1-sample files), vcf_seg_INFO_DP_by_FORMAT_DP (multi sample files)
    int64_t dp;
    bool has_value = str_get_int (STRa(dp_str), &dp);

    if (!has_value || segconf.INFO_DP_method == INFO_DP_DEFAULT || segconf_running) 
        seg_integer_or_not (VB, ctx, STRa(dp_str), dp_str_len);

    else if (segconf.INFO_DP_method == BY_BaseCounts) 
        vb_add_to_deferred_q (VB, ctx, vcf_seg_INFO_DP_by_BaseCounts, vb->idx_DP, INFO_BaseCounts);

    // defer segging to vcf_seg_INFO_DP_by_FORMAT_DP called after samples are done
    else if (segconf.INFO_DP_method == BY_FORMAT_DP)  
        vb_add_to_deferred_q (VB, ctx, vcf_seg_INFO_DP_by_FORMAT_DP, vb->idx_DP, DID_NONE);

    else
        ABOSEG ("Unknown method INFO_DP_method=%d", segconf.INFO_DP_method);

    if (has_value) 
        ctx_set_last_value (VB, ctx, dp);
}

static void vcf_seg_INFO_DP_by_FORMAT_DP (VBlockP vb)
{
    decl_ctx (INFO_DP);

    int value_len = str_int_len (ctx->last_value.i);

    // note: INFO/DP >= sum(FORMAT/DP) as the per-sample value is filtered, see: https://gatk.broadinstitute.org/hc/en-us/articles/360036891012-DepthPerSampleHC
    // note: up to 15.0.35, we had the value_len before the \t
    SNIPi3 (SNIP_SPECIAL, VCF_SPECIAL_deferred_DP, '\t', ctx->last_value.i - ctx->dp.sum_format_dp);
    seg_by_ctx (VB, STRa(snip), ctx, value_len); 
}

static void vcf_seg_INFO_DP_by_BaseCounts (VBlockP vb)
{
    ContextP ctx_dp = CTX(INFO_DP);
    ContextP ctx_basecounts = CTX(INFO_BaseCounts);
    unsigned add_bytes = ctx_dp->last_txt.len;

    if (ctx_has_value_in_line_(vb, ctx_basecounts) && ctx_has_value_in_line_(vb, ctx_dp)) {
        if (ctx_basecounts->last_value.i == ctx_dp->last_value.i) // expected: big majority of cases
            seg_by_ctx (VB, STRa(copy_BaseCounts_sum), ctx_dp, add_bytes);
        
        else {
            int64_t delta = ctx_dp->last_value.i - ctx_basecounts->last_value.i; 

            if (!ctx_dp->snip_cache.len32) {
                buf_alloc_exact (vb, ctx_dp->snip_cache, 48, uint8_t, "snip_cache");
                *Bc(ctx_dp->snip_cache, 0) = SNIP_SPECIAL;
                *Bc(ctx_dp->snip_cache, 1) = VCF_SPECIAL_deferred_DP;

                uint32_t snip2_len = ctx_dp->snip_cache.len32 - 2;
                seg_prepare_snip_other_do (SNIP_OTHER_DELTA, ctx_basecounts->dict_id, true, 0, '$', Bc(ctx_dp->snip_cache, 2), &snip2_len);
                ctx_dp->snip_cache.len32 = snip2_len + 2;
            }

            seg_by_ctx (VB, STRb(ctx_dp->snip_cache), ctx_dp, add_bytes);

            dyn_int_append (vb, ctx_dp, delta, 0); 
        }
    }

    else
        seg_integer_or_not (VB, ctx_dp, STRtxt(ctx_dp->last_txt), add_bytes);
}

// initialize reconstructing INFO/DP by sum(FORMAT/DP) - save space in txt_data, and initialize delta
SPECIAL_RECONSTRUCTOR (vcf_piz_special_deferred_DP)
{
    if (!flag.drop_genotypes && !flag.gt_only && !flag.samples) {
        if (segconf.INFO_DP_method == BY_FORMAT_DP) {
            int64_t sum_format_dp;
            str_item_i_int (STRa(snip), '\t', 1, &sum_format_dp); // note: up to 15.0.35, items[0] was the length of the integer to be inserted. we ignore it now.
            ctx->dp.sum_format_dp = sum_format_dp; // initialize with delta
        }
        vcf_piz_defer (ctx);

        return NO_NEW_VALUE; // we don't have the value yet - it will be set in vcf_piz_insert_INFO_DP
    }
    else {
        if (reconstruct) 
            RECONSTRUCT ("-1", 2); // bc we can't calculate INFO/DP in these cases bc we need FORMAT/DP of all samples
    
        new_value->i = -1;
        return HAS_NEW_VALUE;
    }
}

// finalize reconstructing INFO/DP by sum(FORMAT/DP) - called after reconstructing all samples
void vcf_piz_insert_INFO_DP (VBlockVCFP vb)
{
    decl_ctx (INFO_DP);
    
    if (IS_RECON_INSERTION(ctx)) {
        STRl(info_dp,16);

        if (segconf.INFO_DP_method == BY_FORMAT_DP) {
            info_dp_len = str_int_ex (ctx->dp.sum_format_dp, info_dp, false);
            ctx_set_last_value (VB, ctx, (int64_t)ctx->dp.sum_format_dp); // consumed by eg vcf_piz_insert_INFO_QD
        }
        
        else { // BY_BaseCounts
            rom recon = BAFTtxt;
            STR(snip);
            ctx_get_snip_by_word_index (ctx, ctx->last_wi, snip);
            reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, snip+2, snip_len-2, true, __FUNCTION__, __LINE__);

            info_dp_len = BAFTtxt - recon;
            memcpy (info_dp, recon, info_dp_len);
            Ltxt -= info_dp_len;

            int64_t info_dp_value;
            ASSPIZ (str_get_int (STRa(info_dp), &info_dp_value), "bad INFO_DP=\"%.*s\"", STRf(info_dp));

            ctx_set_last_value (VB, ctx, info_dp_value);
        } 

        vcf_piz_insert_field (vb, ctx, STRa(info_dp));
    }
}

// used starting v13.0.5, replaced in v14 with a new vcf_piz_special_deferred_DP
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DP_by_DP_v13)
{
    str_split (snip, snip_len, 2, '\t', item, 2);

    int num_dps_this_line   = atoi (items[0]);
    int64_t value_minus_sum = atoi (items[1]);

    ContextP format_dp_ctx = CTX(FORMAT_DP);

    int64_t sum=0;

    ASSPIZ (format_dp_ctx->next_local + num_dps_this_line <= format_dp_ctx->local.len, "Not enough data in FORMAT/DP.local to reconstructed INFO/DP: next_local=%u local.len=%u but needed num_dps_this_line=%u",
            format_dp_ctx->next_local, format_dp_ctx->local.len32, num_dps_this_line);
            
    uint32_t invalid = lt_desc[format_dp_ctx->ltype].max_int; // represents '.'
    for (int i=0; i < num_dps_this_line; i++) {
        uint32_t format_dp = (format_dp_ctx->ltype == LT_UINT8)  ? (uint32_t)*B8 ( format_dp_ctx->local, format_dp_ctx->next_local + i)
                           : (format_dp_ctx->ltype == LT_UINT16) ? (uint32_t)*B16 (format_dp_ctx->local, format_dp_ctx->next_local + i)
                           : /* LT_UINT32 */                       (uint32_t)*B32 (format_dp_ctx->local, format_dp_ctx->next_local + i);

        if (format_dp != invalid) sum += format_dp; 
    }

    new_value->i = value_minus_sum + sum;

    RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

static bool vcf_seg_INFO_DP4_delta (VBlockP vb, ContextP ctx, STRp(value), uint32_t unused_rep)
{
    if (ctx_encountered_in_line (vb, (ctx-1)->did_i)) {
        seg_delta_vs_other_localS (VB, ctx, ctx-1, STRa(value), -1);
        return true; // segged successfully
    }
    else
        return false;
}

// <ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
// Expecting first two values to be roughly similar, as the two last bases roughly similar
void vcf_seg_INFO_DP4 (VBlockVCFP vb, ContextP ctx, STRp(dp4))
{
    static const MediumContainer container_DP4 = {
        .repeats      = 1, 
        .nitems_lo    = 4, 
        .items        = { { .dict_id.num = _INFO_DP4_RF, .separator = "," }, 
                          { .dict_id.num = _INFO_DP4_RR, .separator = "," }, 
                          { .dict_id.num = _INFO_DP4_AF, .separator = "," }, 
                          { .dict_id.num = _INFO_DP4_AR                   } } };

    SegCallback callbacks[4] = { 0, vcf_seg_INFO_DP4_delta, 0, vcf_seg_INFO_DP4_delta }; 

    seg_struct (VB, ctx, container_DP4, STRa(dp4), callbacks, dp4_len, true);
}
