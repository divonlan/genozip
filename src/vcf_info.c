// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "piz.h"
#include "optimize.h"
#include "file.h"
#include "dict_id.h"
#include "codec.h"
#include "reconstruct.h"
#include "stats.h"
#include "zip_dyn_int.h"

#define info_items CTX(VCF_INFO)->info_items

STRl(copy_VCF_POS_snip, 16);
sSTRl(copy_BaseCounts_sum, 32);

// called after reading VCF header, before segconf
void vcf_info_zip_initialize (void) 
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _VCF_POS, false, 0, copy_VCF_POS_snip);
        
        // copy_BaseCounts_sum is a DELTA-0 of BaseCounts (=copy) prefixes with SPECIAL
        copy_BaseCounts_sum[0] = SNIP_SPECIAL;
        copy_BaseCounts_sum[1] = VCF_SPECIAL_deferred_DP;
        copy_BaseCounts_sum_len -= 2;
        seg_prepare_snip_other_do (SNIP_OTHER_DELTA, _INFO_BaseCounts, true, 0, 0, copy_BaseCounts_sum + 2, &copy_BaseCounts_sum_len);   
        copy_BaseCounts_sum_len += 2;
    }
}

void vcf_info_seg_initialize (VBlockVCFP vb) 
{
    #define T(cond, did_i) ((cond) ? (did_i) : DID_NONE)

    ctx_set_store (VB, STORE_INT, INFO_AN, INFO_AC, INFO_ADP, INFO_DP, INFO_MLEAC, 
                   INFO_DP4_RF, INFO_DP4_RR, INFO_DP4_AF, INFO_DP4_AR, 
                   INFO_AC_Hom, INFO_AC_Het, INFO_AC_Hemi,
                   DID_EOL);

    CTX(INFO_AF)->flags.store = STORE_FLOAT;
    // xxx (is this really needed for --indels-only?) CTX(INFO_SVTYPE)-> flags.store = STORE_INDEX; // since v13 - consumed by vcf_refalt_piz_is_variant_indel

    CTX(INFO_SF)->sf.SF_by_GT = unknown;
    
    CTX(INFO_AGE_HISTOGRAM_HET)->nothing_char = CTX(INFO_AGE_HISTOGRAM_HOM)->nothing_char = '.';
    
    ctx_set_dyn_int (VB, INFO_SVLEN, INFO_DP4_RF, INFO_DP4_AF,
                     T(segconf.INFO_DP_method == INFO_DP_DEFAULT, INFO_DP),
                     DID_EOL);
    
    if (segconf.INFO_DP_method == BY_BaseCounts)
        ctx_set_same_line (VB, INFO_DP, DID_EOL);

    if (segconf.has[INFO_CLNHGVS]) vcf_seg_hgvs_consolidate_stats (vb, INFO_CLNHGVS);
    if (segconf.has[INFO_HGVSG])   vcf_seg_hgvs_consolidate_stats (vb, INFO_HGVSG);
    if (segconf.has[INFO_ANN])     vcf_seg_hgvs_consolidate_stats (vb, INFO_ANN); // subfield HGVS_c

    #undef T
}

//--------
// INFO/DP
// -------

static void vcf_seg_INFO_DP_by_FORMAT_DP (VBlockP vb); // forward
static void vcf_seg_INFO_DP_by_BaseCounts (VBlockP vb);

// return true if caller still needs to seg 
static void vcf_seg_INFO_DP (VBlockVCFP vb, ContextP ctx, STRp(dp_str))
{
    SEGCONF_RECORD_WIDTH (DP, dp_str_len);

    // used in: vcf_seg_one_sample (for 1-sample files), vcf_seg_INFO_DP_by_FORMAT_DP (multi sample files)
    int64_t dp;
    bool has_value = str_get_int (STRa(dp_str), &dp);

    if (!has_value || segconf.INFO_DP_method == INFO_DP_DEFAULT || segconf.running) 
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
        ctx->dp.is_deferred = true;             // DP needs to be inserted by vcf_piz_insert_INFO_DP

        if (segconf.INFO_DP_method == BY_FORMAT_DP) {
            int64_t sum_format_dp;
            str_item_i_int (STRa(snip), '\t', 1, &sum_format_dp); // note: up to 15.0.35, items[0] was the length of the integer to be inserted. we ignore it now.
            ctx->dp.sum_format_dp = sum_format_dp; // initialize with delta
        }

        vcf_piz_defer_to_after_samples (DP);

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

        vcf_piz_insert_field (vb, ctx, STRa(info_dp), segconf.wid_DP.width);
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
static void vcf_seg_INFO_DP4 (VBlockVCFP vb, ContextP ctx, STRp(dp4))
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

// -------
// INFO/AA
// -------

// return 0 if the allele equals main REF, the alt number if it equals one of the ALTs, or -1 if none or -2 if '.'
static int vcf_INFO_ALLELE_get_allele (VBlockVCFP vb, STRp (value), bool *lowercase)
{
    // check for '.'
    if (IS_PERIOD (value)) return -1;

    if ((*lowercase = str_islower (STRa(value)))) // since 15.0.51
        str_toupper_(value, (char*)value, value_len); // temporarily uppercase

    // check if its equal REF
    int res = -1; // default - not any of the ref or alts

    if (str_issame (value, vb->REF)) 
        res = 0;

    // check if its equal one of the ALTs
    else
        for (int alt_i=0; alt_i < vb->n_alts; alt_i++)
            if (str_issame_(STRa(value), STRi(vb->alt, alt_i)))
                res = alt_i + 1;

    if (*lowercase)
        str_tolower_(value, (char*)value, value_len); // restore

    return res;
}

// checks if value is identifcal to the REF or one of the ALT alleles, and if so segs a SPECIAL snip
// Used for INFO/AA, INFO/CSQ/Allele, INFO/ANN/Allele. Any field using this should have the VCF2VCF_ALLELE translator set in vcf_lo_luft_trans_id.
bool vcf_seg_INFO_allele (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat)  
{
    VBlockVCFP vb = (VBlockVCFP)vb_;
    bool lowercase;

    int allele = vcf_INFO_ALLELE_get_allele (vb, STRa(value), &lowercase);

    // case: this is one of the alleles in REF/ALT - our special alg will just copy from that allele
    if (allele >= 0) {
        char snip[] = { SNIP_SPECIAL, VCF_SPECIAL_ALLELE, '0', '0' + allele /* ASCII 48...147 */, 'L' };
        seg_by_ctx (VB, snip, 4 + lowercase, ctx, value_len);
    }

    // case: a unique allele and no xstrand - we just leave it as-is
    else
        seg_by_ctx (VB, STRa(value), ctx, value_len); 

    return true; // segged successfully
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_ALLELE)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    int allele = snip[1] - '0'; // note: snip[0] was used by DVCF for coordinate, up to 15.0.41
    bool lowercase = snip_len >= 3 && snip[2] == 'L';
    char *recon = BAFTtxt;

    if (allele == 0)
        RECONSTRUCT (vb->REF, vb->REF_len);
    
    else {
        ASSPIZ (allele <= vb->n_alts, "allele=%d but there are only %d ALTs", allele, vb->n_alts);
        RECONSTRUCT (vb->alts[allele-1], vb->alt_lens[allele-1]);
    }

    if (lowercase) // since 15.0.51
        str_tolower_(recon, recon, BAFTtxt - recon); 

    return NO_NEW_VALUE;
}

// -------
// INFO/VT
// -------

// Variant type
void vcf_seg_INFO_VT (VBlockVCFP vb, ContextP ctx, STRp(vt_str))
{
    VariantType vt = vb->var_types[0];
    #define V(s,s_len) str_issame_(STRa(vt_str), s, s_len)

    if (segconf.running && segconf.INFO_VT_type == INFO_VT_UNKNOWN) {
        if (segconf.vcf_is_callmom)
            segconf.INFO_VT_type = INFO_VT_CALLMOM;

        else if (V("Sub", 3) || V("Del", 3) || V("Ins", 3) || V("Complex", 7))
            segconf.INFO_VT_type = INFO_VT_VAGrENT;
        
        else if (V("SNP", 3) || V("INDEL", 5) || V("SV", 2))
            segconf.INFO_VT_type = INFO_VT_1KG;
    }

    if ((segconf.INFO_VT_type == INFO_VT_VAGrENT && 
          ((vt == VT_SNP && V("Sub", 3)) || 
           (vt == VT_DEL && V("Del", 3)) || 
           (vt == VT_INS && V("Ins", 3)))) 
    ||
        (segconf.INFO_VT_type == INFO_VT_1KG &&
          ((has(SVTYPE)  && V("SV", 2))  ||
           (vt == VT_SNP && V("SNP", 3)) || 
           ((vt == VT_INS || vt == VT_DEL) && V("INDEL", 5)))))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_VT, '0' + segconf.INFO_VT_type }, 3, ctx, vt_str_len);

    else // fallback
        seg_by_ctx (VB, STRa(vt_str), ctx, vt_str_len);
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_VT)
{    
    VBlockVCFP vb = (VBlockVCFP)vb_;
    if (!reconstruct) goto done;

    switch (snip[0] - '0') {
        case INFO_VT_VAGrENT:
            switch (vb->var_types[0]) {
                case VT_SNP: RECONSTRUCT ("Sub", 3); break;
                case VT_DEL: RECONSTRUCT ("Del", 3); break;
                case VT_INS: RECONSTRUCT ("Ins", 3); break;
                default    : ABORT ("Unexpected for INFO_VT_VAGrENT: vt=%u", vb->var_types[0]);
            };
            break;

        case INFO_VT_1KG:
            if (curr_container_has (VB, _INFO_SVTYPE))
                RECONSTRUCT ("SV", 2); 

            else switch (vb->var_types[0]) {
                case VT_SNP: RECONSTRUCT ("SNP", 3); break;
                case VT_DEL: 
                case VT_INS: RECONSTRUCT ("INDEL", 5); break;
                default    : ABORT ("Unexpected for INFO_VT_1KG: vt=%u", vb->var_types[0]);
            };
            break;

        default: 
            ABORT ("Invalid INFO_VT_type=%d", snip[0] - '0');
    }

done:
    return NO_NEW_VALUE;
}    

static void vcf_seg_info_one_subfield (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    unsigned modified_len = value_len + 20;
    char modified[modified_len]; // used for 1. fields that are optimized 2. fields translated luft->primary. A_1 transformed 4.321e-03->0.004321
        
    // note: since we use modified for both optimization and luft_back - we currently don't support
    // subfields having both translators and optimization. This can be fixed if needed.

    ctx->line_is_luft_trans = false; // initialize
    
    #define ADJUST_FOR_MODIFIED ({                                  \
        int32_t growth = (int32_t)modified_len - (int32_t)value_len;\
        if (growth) vb->recon_size += growth;                       \
        STRset (value, modified); })                                       

    // many fields, for example: ##INFO=<ID=AN_amr_male,Number=1,Type=Integer,Description="Total number of alleles in male samples of Latino ancestry">
    #define T(i,c) (ctx->tag_name[i] == (c))
    if ((T(0,'A') && T(1,'N') && (T(2,'_') || T(2,'-'))) ||
        ctx->dict_id.num == _INFO_AC_Hom || ctx->dict_id.num == _INFO_AC_Het || ctx->dict_id.num == _INFO_AC_Hemi || // single integer, even if n_alt > 1
        !memcmp (ctx->tag_name, "nhomalt", 7)) 
        seg_integer_or_not (VB, ctx, STRa(value), value_len);

    else if ((T(0,'A') && T(1,'C') && (T(2,'_') || T(2,'-'))) ||
             (T(0,'H') && T(1,'o') && T(2,'m') && T(3,'_')))

        vcf_seg_array_of_N_ALTS_numbers (vb, ctx, STRa(value), STORE_INT);

    else if ((T(0,'A') && T(1,'F') && (T(2,'_') || T(2,'-'))))
        vcf_seg_array_of_N_ALTS_numbers (vb, ctx, STRa(value), STORE_FLOAT);

    #undef T

    else switch (ctx->dict_id.num) {
        #define CALL(f) ({ (f); break; })
        #define CALL_IF(cond,f)  if (cond) { (f); break; } else goto standard_seg 
        #define CALL_IF0(cond,f) if (cond) { (f); break; } else 
        #define CALL_WITH_FALLBACK(f) if (f(vb, ctx, STRa(value))) { seg_by_ctx (VB, STRa(value), ctx, value_len); } break
        #define STORE_AND_SEG(store_type) ({ seg_set_last_txt_store_value (VB, ctx, STRa(value), store_type); seg_by_ctx (VB, STRa(value), ctx, value_len); break; })
        #define DEFER(f,seg_after_did_i) ({ vb_add_to_deferred_q (VB, ctx, vcf_seg_INFO_##f, vb->idx_##f, seg_after_did_i); break; })

        // ---------------------------------------
        // Fields defined in the VCF specification
        // ---------------------------------------
        case _INFO_AC:              CALL (vcf_seg_INFO_AC (vb, ctx, STRa(value))); 
        case _INFO_AF:              CALL (vcf_seg_array_of_N_ALTS_numbers (vb, ctx, STRa(value), STORE_FLOAT));
        case _INFO_AA:              CALL_IF (!segconf.vcf_is_cosmic, vcf_seg_INFO_allele (VB, ctx, STRa(value), 0)); // note: in COSMIC, INFO/AA is something entirely different
        case _INFO_DP:              CALL (vcf_seg_INFO_DP (vb, ctx, STRa(value)));
        case _INFO_END:             CALL (vcf_seg_INFO_END (vb, ctx, STRa(value))); // note: END is an alias of POS - they share the same delta stream - the next POS will be a delta vs this END)

        // ---------------------------------------
        // Population AF fields
        // ---------------------------------------
        case _INFO_AMR_AF:
        case _INFO_AFR_AF:
        case _INFO_EUR_AF:
        case _INFO_ASN_AF:
        case _INFO_SAS_AF:
        case _INFO_EAS_AF:
        case _INFO_AF1000G:
        case _INFO_GNOMAD_AF:       CALL (seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        
        // ---------------------------------------
        // GATK fields
        // ---------------------------------------
        case _INFO_RAW_MQandDP:     CALL_IF (segconf.has[INFO_RAW_MQandDP], vcf_seg_INFO_RAW_MQandDP (vb, ctx, STRa(value)));
        case _INFO_BaseCounts:      DEFER(BaseCounts, DID_NONE);
        case _INFO_SF:              CALL_WITH_FALLBACK (vcf_seg_INFO_SF_init); // Source File
        case _INFO_MLEAC:           CALL (vcf_seg_INFO_MLEAC (vb, ctx, STRa(value)));
        case _INFO_MLEAF:           CALL (vcf_seg_INFO_MLEAF (vb, ctx, STRa(value)));
        case _INFO_QD:              CALL_IF (segconf.has[INFO_QD], DEFER (QD, INFO_DP)); // deferred seg to after samples - then call vcf_seg_INFO_QD
        case _INFO_RU:              CALL (vcf_seg_INFO_RU (vb, ctx, STRa(value)));
        case _INFO_RPA:             CALL (vcf_seg_INFO_RPA (vb, ctx, STRa(value)));
        case _INFO_MFRL:            
        case _INFO_MBQ:
        case _INFO_MMQ:             CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_INT, DICT_ID_NONE, value_len));
        case _INFO_HaplotypeScore:
        case _INFO_R2_5P_bias:
        case _INFO_BaseQRankSum:
        case _INFO_ReadPosRankSum:
        case _INFO_MQRankSum:
        case _INFO_ClippingRankSum:
        case _INFO_SOR:
        case _INFO_FS:
        case _INFO_ExcessHet:
        case _INFO_InbreedingCoeff:
        case _INFO_NALOD:
        case _INFO_NLOD:
        case _INFO_TLOD:            CALL (seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        case _INFO_GERMQ:
        case _INFO_CONTQ:
        case _INFO_SEQQ:
        case _INFO_STRANDQ:
        case _INFO_STRQ:
        case _INFO_ECNT:            CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_AS_SB_TABLE:     CALL_IF (segconf.AS_SB_TABLE_by_SB, DEFER(AS_SB_TABLE, DID_NONE));
        
        case _INFO_VQSLOD: // Optimize VQSLOD 
            if (flag.optimize_VQSLOD && optimize_float_2_sig_dig (STRa(value), 0, qSTRa(modified))) 
                ADJUST_FOR_MODIFIED;
            CALL (seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));

        // ---------------------------------------
        // VEP fields
        // ---------------------------------------
        case _INFO_CSQ:
        case _INFO_vep:             CALL_IF (segconf.vcf_is_vep, vcf_seg_INFO_CSQ (vb, ctx, STRa(value)));
        case _INFO_AGE_HISTOGRAM_HET:
                                    CALL (seg_integer_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, VCF_SPECIAL_N_ALTS, vb->n_alts, value_len));
        case _INFO_AGE_HISTOGRAM_HOM: 
                                    CALL (seg_integer_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, VCF_SPECIAL_N_ALTS, vb->n_alts, value_len));
        case _INFO_DP_HIST:
        case _INFO_GQ_HIST:         CALL (seg_integer_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, VCF_SPECIAL_N_ALLELES, vb->n_alts + 1, value_len));
        case _INFO_DP4:             CALL (vcf_seg_INFO_DP4 (vb, ctx, STRa(value)));
        case _INFO_MC:              CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_CLNDN:           CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), '|', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_CLNID:           CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_NCC:             CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_CLNHGVS: // ClinVar & dbSNP
        case _INFO_HGVSG:   // COSMIC & Mastermind
                                    CALL_IF0 (segconf.vcf_is_mastermind, vcf_seg_mastermind_HGVSG (vb, ctx, STRa(value)))
                                    CALL (vcf_seg_INFO_HGVS (VB, ctx, STRa(value), 0)); 

        // ##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
        // example: CPIC:0b3ac4db1d8e6e08a87b6942|CPIC:647d4339d5c1ddb78daff52f|CPIC:9968ce1c4d35811e7175cd29|CPIC:PA166160951|CPIC:c6c73562e2b9e4ebceb0b8bc
        // I tried seg_array_of_struct - it is worse than simple seg

        case _INFO_ALLELEID:        CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_RSID:            CALL (seg_id_field (VB, ctx, STRa(value), false, value_len));

        // bcftools csq
        case _INFO_BCSQ:            CALL (seg_array (VB, ctx, INFO_BCSQ, STRa(value), ',', '|', false, STORE_NONE, DICT_ID_NONE, value_len));

        // ---------------------------------------
        // SnpEff fields
        // ---------------------------------------
        case _INFO_ANN:             CALL (vcf_seg_INFO_ANN (vb, ctx, STRa(value)));
        case _INFO_EFF:             CALL (vcf_seg_INFO_EFF (vb, ctx, STRa(value)));

        // ---------------------------------------
        // ICGC
        // ---------------------------------------
        case _INFO_mutation:        CALL_IF (segconf.vcf_is_icgc, vcf_seg_INFO_mutation (vb, ctx, STRa(value)));
        case _INFO_CONSEQUENCE:     CALL_IF (segconf.vcf_is_icgc, seg_array (VB, ctx, INFO_CONSEQUENCE, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_OCCURRENCE:      CALL_IF (segconf.vcf_is_icgc, seg_array (VB, ctx, INFO_OCCURRENCE,  STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        // ---------------------------------------
        // COSMIC
        // ---------------------------------------
        case _INFO_LEGACY_ID:       CALL_IF (segconf.vcf_is_cosmic, vcf_seg_INFO_LEGACY_ID   (vb, ctx, STRa(value)));
        case _INFO_SO_TERM:         CALL_IF (segconf.vcf_is_cosmic, vcf_seg_INFO_SO_TERM     (vb, ctx, STRa(value)));

        // ---------------------------------------
        // Mastermind
        // ---------------------------------------
        case _INFO_MMID3:           CALL_IF (segconf.vcf_is_mastermind, vcf_seg_INFO_MMID3   (vb, ctx, STRa(value)));
        case _INFO_MMURI3:          CALL_IF (segconf.vcf_is_mastermind, vcf_seg_INFO_MMURI3  (vb, ctx, STRa(value)));
        case _INFO_MMURI:           CALL_IF (segconf.vcf_is_mastermind, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        case _INFO_GENE:            CALL_IF (segconf.vcf_is_mastermind, STORE_AND_SEG (STORE_NONE)); // consumed by vcf_seg_INFO_MMID3

        // ---------------------------------------
        // Illumina genotyping
        // ---------------------------------------
        case _INFO_PROBE_A:         CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_PROBE_A      (vb, ctx, STRa(value)));
        case _INFO_PROBE_B:         CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_PROBE_B      (vb, ctx, STRa(value)));
        case _INFO_ALLELE_A:        CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ALLELE_A     (vb, ctx, STRa(value)));
        case _INFO_ALLELE_B:        CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ALLELE_B     (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_CHR:    CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_CHR (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_POS:    CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_POS (vb, ctx, STRa(value)));
        case _INFO_ILLUMINA_STRAND: CALL_IF (segconf.vcf_illum_gtyping, vcf_seg_ILLUMINA_STRAND (vb, ctx, STRa(value)));
        case _INFO_refSNP:          CALL_IF (segconf.vcf_illum_gtyping, seg_id_field         (VB, ctx, STRa(value), false, value_len));

        // ---------------------------------------
        // dbSNP
        // ---------------------------------------
        case _INFO_dbSNPBuildID:    CALL_IF (segconf.vcf_is_dbSNP, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_RS:              CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_RS (vb, ctx, STRa(value)));
        case _INFO_RSPOS:           CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_RSPOS (vb, ctx, STRa(value)));
        case _INFO_GENEINFO:        CALL_IF (segconf.vcf_is_dbSNP, seg_array (VB, ctx, INFO_GENEINFO, STRa(value), '|', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_VC:              CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_VC (vb, ctx, STRa(value)));
        case _INFO_FREQ:            CALL_IF (segconf.vcf_is_dbSNP, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        // case _INFO_TOPMED: // better leave as simple snip as the items are allele frequencies which are correleted

        // ---------------------------------------
        // dbNSFP
        // ---------------------------------------
        case _INFO_Polyphen2_HDIV_score : 
        case _INFO_PUniprot_aapos :
        case _INFO_SiPhy_29way_pi : CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));

        case _INFO_VEST3_score    :
        case _INFO_FATHMM_score   : CALL (seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));

        // ---------------------------------------
        // gnomAD
        // ---------------------------------------
        case _INFO_age_hist_het_bin_freq:
        //case _INFO_age_hist_hom_bin_freq: // same dict_id as _INFO_age_hist_het_bin_freq
        case _INFO_gq_hist_alt_bin_freq:
        //case _INFO_gq_hist_all_bin_freq:
        case _INFO_dp_hist_alt_bin_freq:
        //case _INFO_dp_hist_all_bin_freq:
        case _INFO_ab_hist_alt_bin_freq:
                                    CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), '|', 0, false, STORE_INT, DICT_ID_NONE, value_len));

        case _INFO_Genes:           CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_VRS_Allele_IDs:  CALL (vcf_seg_INFO_VRS_Allele_IDs (vb, ctx, STRa(value)));
        case _INFO_VRS_Starts:      CALL (vcf_seg_INFO_VRS_Starts (vb, ctx, STRa(value)));
        case _INFO_VRS_Ends:        CALL (vcf_seg_INFO_VRS_Ends (vb, ctx, STRa(value)));
        case _INFO_VRS_States:      CALL (vcf_seg_INFO_VRS_States (vb, ctx, STRa(value)));
        case _INFO_AS_QD:           CALL (vcf_seg_INFO_AS_QD (vb, ctx, STRa(value)));
        case _INFO_AS_SOR:          CALL (vcf_seg_INFO_AS_SOR (vb, ctx, STRa(value)));
        case _INFO_AS_MQ:           CALL (vcf_seg_INFO_AS_MQ (vb, ctx, STRa(value)));
        case _INFO_AS_MQRankSum:    CALL (vcf_seg_INFO_AS_MQRankSum (vb, ctx, STRa(value)));
        case _INFO_AS_FS:           CALL (vcf_seg_INFO_AS_FS (vb, ctx, STRa(value)));
        case _INFO_AS_QUALapprox:   CALL (vcf_seg_INFO_AS_QUALapprox (vb, ctx, STRa(value)));
        case _INFO_AS_VarDP:        CALL (vcf_seg_INFO_AS_VarDP (vb, ctx, STRa(value)));
        case _INFO_AS_ReadPosRankSum: CALL (vcf_seg_INFO_AS_ReadPosRankSum (vb, ctx, STRa(value)));

        // ---------------------------------------
        // VAGrENT
        // ---------------------------------------
        case _INFO_VD:              CALL_IF (segconf.vcf_is_vagrent, vcf_seg_INFO_VD (vb, ctx, STRa(value)));
        case _INFO_VW:              CALL_IF (segconf.vcf_is_vagrent, vcf_seg_INFO_VW (vb, ctx, STRa(value)));

        // ---------------------------------------
        // Illumina IsaacVariantCaller / starling
        // ---------------------------------------
        case _INFO_REFREP:          CALL_IF (segconf.vcf_is_isaac, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_IDREP:           CALL_IF (segconf.vcf_is_isaac, vcf_seg_INFO_IDREP (vb, ctx, STRa(value)));
        case _INFO_CSQT:            CALL_IF (segconf.vcf_is_isaac, seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_cosmic:          CALL_IF (segconf.vcf_is_isaac, seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_NONE, DICT_ID_NONE, value_len));
        case _INFO_phyloP:          CALL_IF (segconf.vcf_is_isaac, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        case _INFO_SNVHPOL:         CALL_IF (segconf.vcf_is_isaac, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_GMAF:            CALL_IF (segconf.vcf_is_isaac, vcf_seg_FORMAT_GMAF (vb, ctx, STRa(value)));

        // ---------------------------------------
        // Illumina DRAGEN fields
        // ---------------------------------------
        case _INFO_REFLEN:          CALL (vcf_seg_INFO_REFLEN (vb, ctx, STRa(value)));

        // ---------------------------------------
        // Ultima Genomics 
        // ---------------------------------------
        case _INFO_X_LM:
        case _INFO_X_RM:            CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_X_LM_RM (vb, ctx, STRa(value)));
        case _INFO_X_IL:            CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_X_IL (vb, ctx, STRa(value)));
        case _INFO_X_IC:            CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_X_IC (vb, ctx, STRa(value)));
        case _INFO_X_HIN:           CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_X_HIN (vb, ctx, STRa(value)));
        case _INFO_X_HIL:           CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_X_HIL (vb, ctx, STRa(value)));
        case _INFO_X_GCC:           
        case _INFO_HAPDOM:
        case _INFO_TREE_SCORE:      CALL_IF (segconf.vcf_is_ultima, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        case _INFO_VARIANT_TYPE:    CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_VARIANT_TYPE (vb, ctx, STRa(value)));
        case _INFO_HAPCOMP:
        case _INFO_ASSEMBLED_HAPS:  CALL_IF (segconf.vcf_is_ultima, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_FILTERED_HAPS:   CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_FILTERED_HAPS (vb, ctx, STRa(value)));
        
        // ---------------------------------------
        // Platypus
        // ---------------------------------------
        case _INFO_TR:
        case _INFO_TC:
        case _INFO_TCF:
        case _INFO_PP:
        case _INFO_MGOF:            CALL_IF (segconf.vcf_is_platypus, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_SC:              CALL_IF (segconf.vcf_is_platypus, vcf_seg_playpus_INFO_SC (vb, ctx, STRa(value)));
        case _INFO_HP:              CALL_IF (segconf.vcf_is_platypus, vcf_seg_playpus_INFO_HP (vb, ctx, STRa(value)));
        case _INFO_WS:
        case _INFO_WE:              CALL_IF (segconf.vcf_is_platypus, vcf_seg_playpus_INFO_WS_WE (vb, ctx, STRa(value)));
        case _INFO_TCR:             CALL_IF (segconf.vcf_is_platypus, vcf_seg_playpus_INFO_TCR (vb, ctx, STRa(value)));

        // ---------------------------------------
        // Structural variants
        // ---------------------------------------
        case _INFO_CIPOS:           CALL (vcf_seg_INFO_CIPOS (vb, ctx, STRa(value)));
        case _INFO_CIEND:           CALL (vcf_seg_INFO_CIEND (vb, ctx, STRa(value)));
        case _INFO_SVLEN:           CALL (vcf_seg_INFO_SVLEN (vb, ctx, STRa(value)));
        case _INFO_SVTYPE:          CALL (vcf_seg_SVTYPE (vb, ctx, STRa(value)));
        case _INFO_HOMSEQ:          CALL (vcf_seg_HOMSEQ (vb, ctx, STRa(value)));
        case _INFO_MEINFO:          CALL (seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len));
        case _INFO_MATEID:          CALL_IF0(segconf.vcf_is_svaba, vcf_seg_svaba_MATEID (vb, ctx, STRa(value)))
                                    CALL_IF (segconf.vcf_is_pbsv,  vcf_seg_pbsv_MATEID (vb, ctx, STRa(value)));
        case _INFO_MAPQ:            CALL_IF (segconf.vcf_is_svaba, vcf_seg_svaba_MAPQ (vb, ctx, STRa(value)));
        case _INFO_SPAN:            CALL_IF (segconf.vcf_is_svaba, vcf_seg_svaba_SPAN (vb, ctx, STRa(value)));
        case _INFO_SCTG:            CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_SCTG,           TW_SCTG,      false, value_len));
        case _INFO_INSERTION:       CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_INSERTION,      TW_INSERTION, false, value_len));
        case _INFO_EVDNC:           CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_EVDNC,          TW_EVDNC,     false, value_len));
        case _INFO_NUMPARTS:        CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_NUMPARTS,       TW_NUMPARTS,  false, value_len));
        case _INFO_NM:              CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_NM,             TW_MATENM,    false, value_len));
        case _INFO_MATENM:          CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_MATENM,         TW_NM,        false, value_len));
        case _INFO_MATEMAPQ:        CALL_IF (segconf.vcf_is_svaba, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_HOMSEQ,         TW_HOMSEQ,    false, value_len));

        case _INFO_HOMLEN:          CALL_IF (segconf.vcf_is_sv    && has(HOMSEQ),      vcf_seg_LEN_OF (vb, ctx, STRa(value), vb->idx_HOMSEQ,      STRa(homlen_snip)));
        case _INFO_DUPHOMLEN:       CALL_IF (segconf.vcf_is_manta && has(DUPHOMSEQ),   vcf_seg_LEN_OF (vb, ctx, STRa(value), vb->idx_DUPHOMSEQ,   STRa(duphomlen_snip)));
        case _INFO_SVINSLEN:        CALL_IF (segconf.vcf_is_manta && has(SVINSSEQ),    vcf_seg_LEN_OF (vb, ctx, STRa(value), vb->idx_SVINSSEQ,    STRa(svinslen_snip)));
        case _INFO_DUPSVINSLEN:     CALL_IF (segconf.vcf_is_manta && has(DUPSVINSSEQ), vcf_seg_LEN_OF (vb, ctx, STRa(value), vb->idx_DUPSVINSSEQ, STRa(dupsvinslen_snip)));
        case _INFO_CIGAR:           CALL_IF (segconf.vcf_is_manta, vcf_seg_manta_CIGAR (vb, ctx, STRa(value)));
        case _INFO_BND_DEPTH:       CALL_IF (segconf.vcf_is_manta, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_BND_DEPTH, TW_MATE_BND_DEPTH, false, value_len));
        case _INFO_MATE_BND_DEPTH:  CALL_IF (segconf.vcf_is_manta, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_MATE_BND_DEPTH, TW_BND_DEPTH, false, value_len));

        case _INFO_VT:              CALL (vcf_seg_INFO_VT (vb, ctx, STRa(value)));

        default: standard_seg:
            seg_by_ctx (VB, STRa(value), ctx, value_len);
            
            if (ctx->flags.store == STORE_INT) {
                int64_t val;
                if (str_get_int (STRa(value), &val))
                    ctx_set_last_value (VB, ctx, val);
            }
    }

    seg_set_last_txt (VB, ctx, STRa(value));
}

static SORTER (sort_by_subfield_name)
{ 
    InfoItem *ina = (InfoItem *)a;
    InfoItem *inb = (InfoItem *)b;
    
    // END comes first (as eg vcf_INFO_SVLEN_prediction depends on POS.last_delta)
    if (str_issame_(STRa(ina->name), "END=", 4)) return -1;
    if (str_issame_(STRa(inb->name), "END=", 4)) return 1;
    
    return strncmp (ina->name, inb->name, MIN_(ina->name_len, inb->name_len));
}

void vcf_seg_info_subfields (VBlockVCFP vb, STRp(info))
{
    info_items.len = 0; // reset from previous line

    // case: INFO field is '.' (empty) 
    if (IS_PERIOD (info) && !segconf.vcf_is_isaac) { // note: in Isaac, it slightly better to mux the "."
        seg_by_did (VB, ".", 1, VCF_INFO, 2); // + 1 for \t or \n
        return;
    }

    // parse the info string
    str_split (info, info_len, MAX_FIELDS-2, ';', pair, false); // -2 - leave room for LUFT + PRIM
    ASSVCF (n_pairs, "Too many INFO subfields, Genozip supports up to %u", MAX_FIELDS-2);

    buf_alloc (vb, &info_items, 0, n_pairs + 2, InfoItem, CTX_GROWTH, "info_items");

    // pass 1: initialize info items + get indices of AC
    for (unsigned i=0; i < n_pairs; i++) {
        rom equal_sign = memchr (pairs[i], '=', pair_lens[i]);
        unsigned name_len = (unsigned)(equal_sign - pairs[i]); // nonsense if no equal sign
        unsigned tag_name_len = equal_sign ? name_len : pair_lens[i];

        InfoItem ii = { .name_len  = equal_sign ? name_len + 1                : pair_lens[i], // including the '=' if there is one
                        .value     = equal_sign ? equal_sign + 1              : NULL,
                        .value_len = equal_sign ? pair_lens[i] - name_len - 1 : 0  };
        memcpy (ii.name, pairs[i], ii.name_len); // note: we make a copy the name, because vcf_seg_FORMAT_GT might overwrite the INFO field
        
        // create context if it doesn't already exist (also verifies tag is not too long)
        DictId dict_id = dict_id_make (pairs[i], tag_name_len, DTYPE_1);
        ii.ctx = ctx_get_ctx_tag (vb, dict_id, pairs[i], tag_name_len); // create if it doesn't already exist
        
        if (segconf.running) segconf.has[ii.ctx->did_i]++;

        #define X(x) case INFO_##x : vb->idx_##x = info_items.len32; break
        switch (ii.ctx->did_i) {
            X(AN); X(AF); X(AC); X(MLEAC); X(MLEAF); X(AC_Hom); X(AC_Het); X(AC_Hemi); X(DP); X(QD); X(SF);
            X(AS_SB_TABLE); X(SVINSSEQ) ; X(SVTYPE); X(HOMSEQ); X(END) ; X(SVLEN); X(CIPOS); X(LEFT_SVINSSEQ);
            X(BaseCounts);
            default: {}
        }
        #undef X

        BNXT (InfoItem, info_items) = ii;
    }

    ARRAY (InfoItem, ii, info_items);

    // pass 2: seg all subfields except AC
    for (unsigned i=0; i < ii_len; i++)
        if (ii[i].value) 
            vcf_seg_info_one_subfield (vb, ii[i].ctx, STRa(ii[i].value));
        else
            ctx_set_encountered (VB, ii[i].ctx);
}

// Seg INFO fields that were deferred to after all samples are segged
void vcf_seg_finalize_INFO_fields (VBlockVCFP vb)
{
    if (!info_items.len) return; // no INFO items on this line (except if dual-coords - we will add them in a sec)

    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true }; 
 
    ARRAY (InfoItem, ii, info_items);

    con_set_nitems (con, ii_len);

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag.optimize_sort && ii_len > 1) 
        qsort (ii, ii_len, sizeof(InfoItem), sort_by_subfield_name);

    char prefixes[CONTAINER_MAX_PREFIXES_LEN];  // these are the Container prefixes
    prefixes[0] = prefixes[1] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Container 
    uint32_t total_names_len=0;
    for (unsigned i=0; i < ii_len; i++) {
        // Set the Container item and find (or create) a context for this name
        con.items[i] = (ContainerItem){ .dict_id   = !ii[i].value                        ? DICT_ID_NONE 
                                                   : ii[i].ctx->dict_id.num == _INFO_END ? (DictId)_VCF_POS
                                                   :                                       ii[i].ctx->dict_id,
                                        .separator = { ';' } }; 
            
        // add to the prefixes
        ASSVCF (prefixes_len + ii[i].name_len + 1 <= CONTAINER_MAX_PREFIXES_LEN, 
                "INFO contains tag names that, combined (including the '='), exceed the maximum of %u characters", CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii[i].name, ii[i].name_len);
        prefixes_len += ii[i].name_len;
        prefixes[prefixes_len++] = CON_PX_SEP;

        total_names_len += ii[i].name_len + 1; // +1 for ; \t or \n separator
    }

    // seg deferred fields 
    if (flag.debug_seg) vb_display_deferred_q (VB, __FUNCTION__);

    for (uint8_t i=0; i < vb->deferred_q_len; i++) 
        vb->deferred_q[i].seg (VB);

    // case INFO is muxed: multiplex by has_RGQ or FILTER in Isaac
    if (!segconf.running && (segconf.has[FORMAT_RGQ] || segconf.vcf_is_isaac)) {
        ContextP channel_ctx = 
            seg_mux_get_channel_ctx (VB, VCF_INFO, (MultiplexerP)&vb->mux_INFO, (segconf.has[FORMAT_RGQ] ? CTX(FORMAT_RGQ)->line_has_RGQ : vcf_isaac_info_channel_i (VB)));
        
        seg_by_did (VB, STRa(vb->mux_INFO.snip), VCF_INFO, 0);

        container_seg (vb, channel_ctx, &con, prefixes, prefixes_len, total_names_len /* names inc. = and separator */);
    }

    // case: INFO not muxed
    else 
        container_seg (vb, CTX(VCF_INFO), &con, prefixes, prefixes_len, total_names_len /* names inc. = and separator */);
}

