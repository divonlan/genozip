// ------------------------------------------------------------------
//   vcf_gatk.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent Pending.
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

#define adjustment ctx->last_delta

static SmallContainer RAW_MQandDP_con = {
    .nitems_lo = 2, 
    .repeats   = 1, 
    .items     = { { .dict_id={ _INFO_RAW_MQandDP_MQ }, .separator = {','} },
                   { .dict_id={ _INFO_RAW_MQandDP_DP }                     } }
};

sSTRl(RAW_MQandDP_snip, 64 + 2 * 16);     

void vcf_gatk_zip_initialize (void)
{
    container_prepare_snip ((ContainerP)&RAW_MQandDP_con, 0, 0, qSTRa(RAW_MQandDP_snip));
}

void vcf_gatk_seg_initialize (VBlockVCFP vb)
{
    ctx_set_ltype (VB, LT_STRING, INFO_TLOD, INFO_NALOD, INFO_NLOD, DID_EOL);

    ctx_set_dyn_int (VB, INFO_RPA, DID_EOL);
}

// ------------------
// INFO/RAW_MQandDP
// ------------------

// ##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
// comma-seperated two numbers: RAW_MQandDP=720000,200: 1. sum of squared MQ values and 2. total reads over variant genotypes (note: INFO/MQ is sqrt(#1/#2))
void vcf_seg_INFO_RAW_MQandDP (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    seg_by_container (VB, ctx, (ContainerP)&RAW_MQandDP_con, STRa(value), STRa(RAW_MQandDP_snip), NULL, true, value_len);
}

// -----------------
// INFO/RU  : Tandem repeat unit (bases)
// INFO/RPA : Number of times tandem repeat unit is repeated, for each allele (including reference)
// Observed examples:
// REF=TA ALT=T RU=A RPA=9,8
// REF=G ALT=GT RU=T RPA=9,10
// REF=C ALT=CCCCT RU=CCCT RPA=3,4
// REF=CTTTTTT ALT=C,CT,CTT,CTTT,CTTTT,CTTTTT RU=T RPA=22,16,17,18,19,20,21
// REF=CTTAT ALT=C RU=TTAT RPA=2,1
// -----------------

static bool is_repeat (STRp(seq), int ru_len)
{
    for (int rep_i=1; rep_i < seq_len / ru_len; rep_i++)
        if (memcmp (seq, &seq[ru_len * rep_i], ru_len)) return false;

    return true;
}

static void predict_RU (VBlockVCFP vb)
{
    TxtWord *predicted_RU = &CTX(INFO_RU)->predicted_RU;

    if (predicted_RU->index) return; // already predicted

    rom seq = vb->REF + 1; // payload only (without anchor base)
    int seq_len = vb->REF_len - 1;

    // set seq to longest allele (REF or one of the ALTs)
    for (int alt_i=0; alt_i < vb->n_alts; alt_i++)
        if (vb->alt_lens[alt_i] > seq_len + 1) {
            seq     = vb->alts[alt_i]     + 1; // payload only (without anchor base)
            seq_len = vb->alt_lens[alt_i] - 1; 
        }

    // shortcut common cases. note: seq_len=0 happens for SNPs - these are not expected to have RU/RPA fields, but we support for completeness
    if (seq_len <= 1)
        *predicted_RU = TXTWORD (seq);
        
    else if (str_is_monochar (STRa(seq))) {
        seq_len = 1;
        *predicted_RU = TXTWORD (seq);
    }

    else {
        *predicted_RU = TXTWORD (seq); // initialize: the whole payload is a single "repeat" (modified if we find a shorter repeat)

        for (int ru_len=2; ru_len <= seq_len / 2; ru_len++)
            if ((seq_len % ru_len)==0 && is_repeat (STRa(seq), ru_len)) {
                seq_len = ru_len;
                *predicted_RU = TXTWORD (seq);
                break;
            }
    }
}

// Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs longer than 20 bases are not reported.
// Appears in GATK, Isaac Variant Caller / starling
void vcf_seg_INFO_RU (VBlockVCFP vb, ContextP ctx, STRp(ru))
{
    predict_RU (vb);

    txtSTR (predicted_RU, ctx->predicted_RU);

    if (str_issame (ru, predicted_RU))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_RU, '1'/*version*/ }, 3, ctx, ru_len);

    else
        seg_by_ctx (VB, STRa(ru), ctx, ru_len); 
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_RU)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (reconstruct) {
        if (!snip_len) { // 15.0.13 to 15.0.40
            if (vb->REF_len > 1) 
                RECONSTRUCT (vb->REF+1, vb->REF_len-1);
            else
                RECONSTRUCT (vb->ALT+1, vb->ALT_len-1);
        }

        else { // since 15.0.41
            predict_RU (vb);

            txtSTR (predicted_RU, ctx->predicted_RU);
            RECONSTRUCT (predicted_RU, predicted_RU_len);        
        }
    }
    return NO_NEW_VALUE;
}

// ZIP/PIZ
static int get_n_repeats (STRp(seq), STRp(ru))
{
    seq++; seq_len--; // skip anchor base

    if (!seq_len) return 0; // trivially true

    if (seq_len % ru_len) return -1;  // fail: length is not a multiple of predicted_RU_len

    for (int rep_i=0; rep_i < seq_len / ru_len; rep_i++)
        if (memcmp (&seq[ru_len * rep_i], ru, ru_len)) return -1; // fail: not a repeat of RU

    return seq_len / ru_len;
}

void vcf_seg_INFO_RPA (VBlockVCFP vb, ContextP ctx, STRp(rpa_str))
{
    predict_RU (vb);
    txtSTR (ru, CTX(INFO_RU)->predicted_RU);

    str_split_ints (rpa_str, rpa_str_len, vb->n_alts + 1, ',', rpa, true);
    if (!n_rpas) goto fallback;

    int visible_reps = get_n_repeats (STRa(vb->REF), STRa(ru));
    if (visible_reps == -1) goto fallback;

    int delta = rpas[0] - visible_reps;

    for (int alt_i=0; alt_i < vb->n_alts; alt_i++)
        if (get_n_repeats (STRi(vb->alt, alt_i), STRa(ru)) + delta != rpas[alt_i + 1])
            goto fallback;

    seg_integer (VB, ctx, delta, false, 0); // add delta to local
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_RPA }, 2, ctx, rpa_str_len);
    return;

fallback:
    seg_by_ctx (VB, STRa(rpa_str), ctx, rpa_str_len);
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_RPA)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    predict_RU (vb);
    txtSTR (ru, CTX(INFO_RU)->predicted_RU);

    int delta = reconstruct_from_local_int (VB, ctx, 0, false);

    if (reconstruct) {
        RECONSTRUCT_INT (get_n_repeats (STRa(vb->REF), STRa(ru)) + delta);

        for (int alt_i=0; alt_i < vb->n_alts; alt_i++) {
            RECONSTRUCT1 (',');
            RECONSTRUCT_INT (get_n_repeats (STRi(vb->alt, alt_i), STRa(ru)) + delta);
        }
    }

    return NO_NEW_VALUE;
}

// ---------------
// INFO/BaseCounts
// ---------------

// ##INFO=<ID=genozip BugP.vcf -ft ,Number=4,Type=Integer,Description="Counts of each base">
// Sorts BaseCounts vector with REF bases first followed by ALT bases, as they are expected to have the highest values
bool vcf_seg_INFO_BaseCounts (VBlockVCFP vb, ContextP ctx_basecounts, STRp(value)) // returns true if caller still needs to seg 
{
    if (vb->REF_len != 1 || vb->ALT_len != 1) 
        return true; // not a bi-allelic SNP or line is a luft line without easy access to REFALT - caller should seg

    char *str = (char *)value;
    int64_t sum = 0;

    uint32_t counts[4], sorted_counts[4] = {}; // corresponds to A, C, G, T

    SAFE_NUL (&value[value_len]);
    for (unsigned i=0; i < 4; i++) {
        counts[i] = strtoul (str, &str, 10);
        str++; // skip comma separator
        sum += counts[i];
    }
    SAFE_RESTORE;

    if (str - value != value_len + 1 /* +1 due to final str++ */) return true; // invalid BaseCounts data - caller should seg

    unsigned ref_i = acgt_encode[(int)*vb->REF];
    unsigned alt_i = acgt_encode[(int)*vb->ALT];

    bool used[4] = {};
    sorted_counts[0] = counts[ref_i]; // first - the count of the REF base
    sorted_counts[1] = counts[alt_i]; // second - the count of the ALT base
    used[ref_i] = used[alt_i] = true;

    // finally - the other two cases in their original order (usually these are 0)
    for (unsigned sc_i=2; sc_i <= 3; sc_i++)
        for (unsigned c_i=0; c_i <= 3; c_i++)
            if (!used[c_i]) { // found a non-zero count
                sorted_counts[sc_i] = counts[c_i];
                used[c_i] = true;
                break;
            }

    char snip[2 + value_len + 1]; // +1 for \0
    sprintf (snip, "%c%c%u,%u,%u,%u", SNIP_SPECIAL, VCF_SPECIAL_BaseCounts, 
             sorted_counts[0], sorted_counts[1], sorted_counts[2], sorted_counts[3]);

    seg_by_ctx (VB, snip, value_len+2, ctx_basecounts, value_len); 
    
    ctx_basecounts->flags.store = STORE_INT;
    ctx_set_last_value (VB, ctx_basecounts, sum);

    return false; // we already segged - caller needn't seg
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_INFO_BaseCounts)
{
    uint32_t counts[4], sorted_counts[4] = {}; // counts of A, C, G, T

    new_value->i = 0;
    char *str = (char *)snip;

    for (unsigned i=0; i < 4; i++) {
        sorted_counts[i] = strtoul (str, &str, 10);
        str++; // skip comma separator
        new_value->i += sorted_counts[i];
    }

    if (!reconstruct) goto done; // just return the new value

    ASSVCF (str - snip == snip_len + 1, "expecting (str-snip)=%d == (snip_len+1)=%u", (int)(str - snip), snip_len+1);

    unsigned ref_i = acgt_encode[(uint8_t)VB_VCF->REF[0]];
    unsigned alt_i = acgt_encode[(uint8_t)VB_VCF->alts[0][0]];
    
    counts[ref_i] = sorted_counts[0];
    counts[alt_i] = sorted_counts[1];
    
    unsigned sc_i=2;
    for (unsigned i=0; i <= 3; i++)
        if (ref_i != i && alt_i != i) counts[i] = sorted_counts[sc_i++];

    bufprintf (vb, &vb->txt_data, "%u,%u,%u,%u", counts[0], counts[1], counts[2], counts[3]);

done:
    return HAS_NEW_VALUE;
}

// ---------------
// INFO/SF
// ---------------

static void vcf_seg_INFO_SF_seg (VBlockP vb); // forward

// INFO/SF contains a comma-seperated list of the 0-based index of the samples that are NOT '.'
// example: "SF=0,1,15,22,40,51,59,78,88,89,90,112,124,140,147,155,156,164,168,183,189,197,211,215,216,217,222,239,244,256,269,270,277,281,290,291,299,323,338,340,348" 
// Algorithm: SF is segged either as an as-is string, or as a SPECIAL that includes the index of all the non-'.' samples. 
// if sf.SF_by_GT=YES, we use SNIP_SPECIAL and we validate the correctness during vcf_seg_FORMAT_GT -
// if it is wrong we set sf.SF_by_GT=NO. The assumption is that normally, it is either true for all lines or false.
bool vcf_seg_INFO_SF_init (VBlockVCFP vb, ContextP ctx, STRp(sf))
{
    SEGCONF_RECORD_WIDTH (SF, sf_len);

    switch (ctx->sf.SF_by_GT) {
        case no: 
            return true; // "special" is suppressed - caller should go ahead and seg normally

        case unknown: 
            ctx->sf.SF_by_GT = yes; // first call to this function, after finding that we have an INFO/SF field - set and fall through

        case yes: 
            seg_set_last_txt (VB, ctx, STRa(sf));

            vb_add_to_deferred_q (VB, ctx, vcf_seg_INFO_SF_seg, vb->idx_SF);

            ctx->sf.next = 0; // relative to last_txt
            adjustment = 0;      
            
            // snip being contructed 
            buf_alloc (vb, &ctx->deferred_snip, 0, sf_len + 20, char, 2, "contexts->deferred_snip"); // initial value - we will increase if needed
            BNXTc (ctx->deferred_snip) = SNIP_SPECIAL;
            BNXTc (ctx->deferred_snip) = VCF_SPECIAL_SF;

            return false; // caller should not seg as we already did

        default:
            ABOSEG ("invalid sf.SF_by_GT=%d", ctx->sf.SF_by_GT);
    }
}

// verify next number on the list of the SF field is sample_i (called from vcf_seg_FORMAT_GT)
void vcf_seg_INFO_SF_one_sample (VBlockVCFP vb)
{
    decl_ctx(INFO_SF);
    STRlast (sf, INFO_SF);

    // case: no more SF values left to compare - we ignore this sample
    while (ctx->sf.next < sf_len) {

        buf_alloc (vb, &ctx->deferred_snip, 10, 0, char, 2, "contexts->deferred_snip");

        rom sf_one_value = &sf[ctx->sf.next]; 
        char *after;
        int32_t value = strtol (sf_one_value, &after, 10); // terminated with , ; \t (not \n bc we know there are samples)

        int32_t adjusted_sample_i = (int32_t)(vb->sample_i + adjustment); // adjustment is the number of values in SF that are not in samples

        // case: badly formatted SF field
        if (*after != ',' && *after != ';' && *after != '\t') {
            CTX(INFO_SF)->sf.SF_by_GT = no; // failed - turn off for the remainder of this vb
            break;
        }
        
        // case: value exists in SF and samples
        else if (value == adjusted_sample_i) {
            BNXTc (ctx->deferred_snip) = ',';
            ctx->sf.next = after - sf + 1; // +1 to skip comma
            break;
        }

        // case: value in SF file doesn't appear in samples - keep the value in the snip
        else if (value < adjusted_sample_i) {
            ctx->deferred_snip.len += str_int (value, BAFTc (ctx->deferred_snip));
            BNXTc (ctx->deferred_snip) = ',';
            adjustment++;
            ctx->sf.next = after - sf + 1; // +1 to skip comma
            // continue and read the next value
        }

        // case: value in SF is larger than current sample - don't advance iterator - perhaps future sample will cover it
        else { // value > adjusted_sample_i
            BNXTc (ctx->deferred_snip) = '~'; // skipped sample
            break; 
        }
    }
}

static void vcf_seg_INFO_SF_seg (VBlockP vb_)
{   
    VBlockVCFP vb = (VBlockVCFP)vb_;
    decl_ctx(INFO_SF);
    STRlast (sf, INFO_SF);

    // case: SF data remains after all samples - copy it
    int32_t remaining_len = (int32_t)ctx->last_txt.len - (int32_t)ctx->sf.next; // -1 if all done, because we skipped a non-existing comma
    if (remaining_len > 0) {
        buf_add_more (VB, &ctx->deferred_snip, &sf[ctx->sf.next], remaining_len, "contexts->deferred_snip");
        BNXTc (ctx->deferred_snip) = ','; // buf_add_more allocates one character extra
    }

    if (CTX(INFO_SF)->sf.SF_by_GT == yes) 
        seg_by_ctx (VB, STRb(ctx->deferred_snip), ctx, sf_len);
    
    else if (CTX(INFO_SF)->sf.SF_by_GT == no)
        seg_by_ctx (VB, STRa(sf), ctx, sf_len);

    buf_free (ctx->deferred_snip);
}

#define sample_i ctx->last_value.i
#define snip_i   ctx->deferred_snip.param

// leave space for reconstructing SF - actual reconstruction will be in vcf_piz_container_cb
SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_INFO_SF)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (reconstruct) {
        if (!flag.drop_genotypes && !flag.samples) {
            adjustment = sample_i = snip_i = 0;

            // temporary place for SF
            buf_alloc (vb, &ctx->insertion, 0, 5 * vcf_header_get_num_samples(), char, 1, "contexts->insertion"); // initial estimate, we may further grow it later
            ctx->insertion.len = 0;

            // save snip for later (note: the SNIP_SPECIAL+code are already removed)
            buf_add_moreS (vb, &ctx->deferred_snip, snip, "contexts->deferred_snip");

            vcf_piz_defer_to_after_samples (SF);
        }
        else
            RECONSTRUCT ("N/A", 3);
    }

    return NO_NEW_VALUE;
}

// While reconstructing the GT fields of the samples - calculate the INFO/SF field
void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vb, unsigned rep, char *recon, int32_t recon_len)
{
    decl_ctx (INFO_SF);

    if (rep != 0 || !IS_RECON_INSERTION(ctx)) return; // we only look at the first ht in a sample, and only if its not '.'/'%'

    if (*recon == '.' || *recon == '%') { // . can be written as % in vcf_seg_FORMAT_GT
        sample_i++;
        return;
    }

    ARRAY (const char, snip, ctx->deferred_snip);

    while (snip_i < snip_len) {

        buf_alloc (vb, &ctx->insertion, 12, 0, char, 2, "insertion"); // sufficient for int32 + ','

        int32_t adjusted_sample_i = (int32_t)(sample_i + adjustment); // last_delta is the number of values in SF that are not in samples

        // case: we derive the value from sample_i
        if (snip[snip_i] == ',') {
            snip_i++;

            buf_add_int_as_text (&ctx->insertion, adjusted_sample_i);

            if (snip_i < snip_len) // add comma if not done yet
                BNXTc (ctx->insertion) = ',';

            break;
        }

        // case: sample was not used to derive any value
        if (snip[snip_i] == '~') {
            snip_i++;
            break;
        }

        // case: value quoted verbatim - output and continue searching for , or ~ for this sample
        else {
            char *comma = strchr (&snip[snip_i], ',');
            unsigned len = comma - &snip[snip_i];
            bool add_comma = (snip_i + len + 1 < snip_len); // add comma if not last

            buf_add (&ctx->insertion, &snip[snip_i], len + add_comma);
            snip_i += len + 1;

            adjustment++;
        }
    }

    sample_i++;
}

// Upon completing the line - insert the calculated INFO/SF field to its place
void vcf_piz_insert_INFO_SF (VBlockVCFP vb)
{
    decl_ctx(INFO_SF);

    if (!IS_RECON_INSERTION(ctx)) return;

    ARRAY (const char, snip, ctx->deferred_snip);

    // if there are some items remaining in the snip (values that don't appear in samples) - copy them
    if (snip_i < snip_len) {
        unsigned remaining_len = snip_len - snip_i - 1; // all except the final comma
        buf_add_more ((VBlockP)vb, &ctx->insertion, &snip[snip_i], remaining_len, "contexts->insertion"); 
    }

    // make room for the SF txt and copy it to its final location
    vcf_piz_insert_field (vb, ctx, STRb(ctx->insertion), segconf.wid_SF.width);

    buf_free (ctx->deferred_snip);
    buf_free (ctx->insertion);
}

// ---------------
// INFO/QD
// ---------------

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
    if ((vcf_num_samples + !has_info_dp >= 2) &&  
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

// called from toplevel callback
void vcf_piz_insert_INFO_QD (VBlockVCFP vb)
{
    decl_ctx (INFO_QD);
    if (!IS_RECON_INSERTION(ctx)) return;

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

// just store pred_type - actual reconstruction is deferred to vcf_piz_insert_INFO_QD
SPECIAL_RECONSTRUCTOR (vcf_piz_special_QD)
{
    ctx->qd.pred_type = (uint32_t)(snip[0] - '0');

    ASSPIZ (ctx->qd.pred_type < NUM_QD_PRED_TYPES, 
            "Unknown pred_type=%d. %s", ctx->qd.pred_type, genozip_update_msg());

    vcf_piz_defer_to_after_samples (QD);

    return NO_NEW_VALUE;
}

// -------------------
// INFO/AS_SB_TABLE
// -------------------

void vcf_seg_INFO_AS_SB_TABLE (VBlockP vb)
{
    decl_ctx (INFO_AS_SB_TABLE);

    STRlast (as_sb_table, INFO_AS_SB_TABLE);
    SAFE_NULT (as_sb_table);

    if (VB_VCF->n_alts != 1) goto fallback;

    SEGCONF_RECORD_WIDTH (AS_SB_TABLE, as_sb_table_len);

    char seps[4] = {',', '|', ',', '\0'}; // structure expected if bi-allelic

    // extract each item, from eg: 62,42|18,12
    char *next = (char *)as_sb_table;
    for (int i=0; i < 4; i++) {
        uint32_t val = strtoul (next, &next, 10);
        if (*next != seps[i]) goto fallback; // check if format is as expected for bi-alleic AS_SB_TABLE
    
        // Limitation: sum_sb is only 16bit. Will likely fail on files with many samples or extremely high coverage.
        if (val != CTX(FORMAT_SB)->sum_sb[i]) goto fallback; // check if preditable from ΣSB
        next++; // skip separator
    }

    if (next-1 != as_sb_table + as_sb_table_len) goto fallback; // verify that entire string was consumed

    seg_by_ctx (vb, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_DEFER }, 2, ctx, as_sb_table_len);
    goto done;

fallback:
    seg_by_ctx (vb, STRa(as_sb_table), ctx, as_sb_table_len);

done:
    SAFE_RESTORE;
}

void vcf_piz_sum_SB_for_AS_SB_TABLE (VBlockP vb, STRp(recon))
{
    str_split_ints (recon, recon_len, 4, ',', sb, true);
    
    if (n_sbs == 4) 
        for (int i=0; i < 4; i++)
            CTX(FORMAT_SB)->sum_sb[i] += sbs[i];
}

void vcf_piz_insert_INFO_AS_SB_TABLE (VBlockVCFP vb)
{
    decl_ctx (INFO_AS_SB_TABLE);
    
    if (!IS_RECON_INSERTION(ctx)) return;
    
    uint16_t *sb = CTX(FORMAT_SB)->sum_sb;

    char as_sb_table[64];
    int as_sb_table_len = sprintf (as_sb_table, "%u,%u|%u,%u", sb[0], sb[1], sb[2], sb[3]);

    vcf_piz_insert_field (vb, ctx, STRa(as_sb_table), segconf.wid_AS_SB_TABLE.width);
}

// defer reconstruction after reconstruction of FORMAT/SB - happens in vcf_piz_insert_INFO_AS_SB_TABLE
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DEFER)
{
    switch (ctx->did_i) {
        case INFO_AS_SB_TABLE : vcf_piz_defer_to_after_samples (AS_SB_TABLE); break;
        case VCF_ID           : vcf_piz_defer_to_later (ID); break;
        default               : ABORT_PIZ ("unsupported ctx=%s", ctx->tag_name);
    }

    return NO_NEW_VALUE;
}
