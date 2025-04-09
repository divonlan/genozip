// ------------------------------------------------------------------
//   vcf_info.c
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "stats.h"
#include "zip_dyn_int.h"

STRl(copy_VCF_POS_snip, 16);
STRl(copy_VCF_ID_snip, 16);
STRl(copy_INFO_AF_snip, 16);
STRl(copy_BaseCounts_sum, 32);

// called after reading VCF header, before segconf
void vcf_info_zip_initialize (void) 
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _VCF_POS, false, 0, copy_VCF_POS_snip);
        seg_prepare_snip_other (SNIP_COPY, _VCF_ID, false, 0, copy_VCF_ID_snip);
        seg_prepare_snip_other (SNIP_COPY, _INFO_AF, false, 0, copy_INFO_AF_snip);
        
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
    ctx_set_store (VB, STORE_INT, INFO_AN, INFO_AC, INFO_ADP, INFO_DP, INFO_MLEAC, 
                   INFO_DP4_RF, INFO_DP4_RR, INFO_DP4_AF, INFO_DP4_AR, 
                   INFO_AC_Hom, INFO_AC_Het, INFO_AC_Hemi,
                   INFO_MA,
                   DID_EOL);

    ctx_set_store (VB, STORE_FLOAT, INFO_AF, DID_EOL);

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
        for_alt2
            if (str_issame_(STRa(value), STRa(alt->alt)))
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

    // case: an allele that is not one of REF or ALT - we just leave it as-is
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
        ASSPIZ (allele <= N_ALTS, "allele=%d but there are only %d ALTs", allele, N_ALTS);
        RECONSTRUCT_str (ALTi(allele-1)->alt);
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
    VariantType vt = ALTi(0)->var_type;
    #define V(s,s_len) str_issame_(STRa(vt_str), s, s_len)

    if (segconf_running && segconf.INFO_VT_type == INFO_VT_UNKNOWN) {
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
           ((vt == VT_INS || vt == VT_DEL) && V("INDEL", 5))))
    ||
        (segconf.INFO_VT_type == INFO_VT_CALLMOM &&
           ((vt == VT_SNP && V("S", 1))  ||
            ((vt == VT_INS || vt == VT_DEL) && V("I", 1)) ||
            ((vt != VT_SNP && vt != VT_INS && vt != VT_DEL) && V("M", 1)))))

        seg_special1 (VB, VCF_SPECIAL_VT, '0' + segconf.INFO_VT_type, ctx, vt_str_len);

    else // fallback
        seg_by_ctx (VB, STRa(vt_str), ctx, vt_str_len);
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_VT)
{    
    VBlockVCFP vb = (VBlockVCFP)vb_;
    if (!reconstruct) goto done;

    VariantType vt = ALTi(0)->var_type;

    switch (snip[0] - '0') {
        case INFO_VT_VAGrENT:
            switch (vt) {
                case VT_SNP: RECONSTRUCT ("Sub", 3); break;
                case VT_DEL: RECONSTRUCT ("Del", 3); break;
                case VT_INS: RECONSTRUCT ("Ins", 3); break;
                default    : ABORT ("Unexpected for INFO_VT_VAGrENT: vt=%u", vt);
            };
            break;

        case INFO_VT_1KG:
            if (curr_container_has (VB, _INFO_SVTYPE))
                RECONSTRUCT ("SV", 2); 

            else switch (vt) {
                case VT_SNP: RECONSTRUCT ("SNP", 3); break;
                case VT_DEL: 
                case VT_INS: RECONSTRUCT ("INDEL", 5); break;
                default    : ABORT ("Unexpected for INFO_VT_1KG: vt=%u", vt);
            };
            break;

        case INFO_VT_CALLMOM:
            switch (vt) {
                case VT_SNP: RECONSTRUCT1 ('S'); break;
                case VT_DEL:
                case VT_INS: RECONSTRUCT1 ('I'); break;
                default:     RECONSTRUCT1 ('M'); break;
            }
            break;

        default: 
            ABORT ("Invalid INFO_VT_type=%d", snip[0] - '0');
    }

done:
    return NO_NEW_VALUE;
}    

// Minor Allele - expected to be identical to REF or one of the ALTs
bool vcf_seg_INFO_MA (VBlockVCFP vb, ContextP ctx, STRp(ma))
{
    for_alt2
        if (str_issame(ma, alt->alt)) {
            char snip[12] = { SNIP_SPECIAL, VCF_SPECIAL_COPY_REForALT };
            int snip_len = 2 + str_int (alt_i + 1, &snip[2]);
            seg_by_ctx (VB, STRa(snip), ctx, ma_len);

            ctx_set_last_value (VB, ctx, (int64_t)(alt_i + 1)); // to be used by for segging AF against MAF
            return false;
        }

    if (str_issame(ma, vb->REF)) {
        seg_special1 (VB, VCF_SPECIAL_COPY_REForALT, '0', ctx, ma_len);

        ctx_set_last_value (VB, ctx, (int64_t)0);
        return false;
    }

    return true; // caller should seg
}

void vcf_seg_string (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    if (flag.optimize)
        seg_by_ctx (VB, STRa(value), ctx, value_len);
    
    else
        seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len);
}

static void vcf_seg_info_one_subfield (VBlockVCFP vb, ContextP ctx, STRp(value))
{            
    // many fields, for example: ##INFO=<ID=AN_amr_male,Number=1,Type=Integer,Description="Total number of alleles in male samples of Latino ancestry">
    #define N(i,c) (ctx->tag_name[i] == (c))
    if ((N(0,'A') && N(1,'N') && (N(2,'_') || N(2,'-'))) ||
        ctx->dict_id.num == _INFO_AC_Hom || ctx->dict_id.num == _INFO_AC_Het || ctx->dict_id.num == _INFO_AC_Hemi || // single integer, even if n_alt > 1
        !memcmp (ctx->tag_name, "nhomalt", 7)) 
        seg_integer_or_not (VB, ctx, STRa(value), value_len);

    else if ((N(0,'A') && N(1,'C') && (N(2,'_') || N(2,'-'))) ||
             (N(0,'H') && N(1,'o') && N(2,'m') && N(3,'_')))

        vcf_seg_array_of_N_ALTS_numbers (vb, ctx, STRa(value), STORE_INT);

    else if ((N(0,'A') && N(1,'F') && (N(2,'_') || N(2,'-'))))
        vcf_seg_array_of_N_ALTS_numbers (vb, ctx, STRa(value), STORE_FLOAT);

    #undef N

    else switch (ctx->dict_id.num) {
        #define CALL(f) ({ (f); break; })
        #define CALL_IF(cond,f)  if (cond) { (f); break; } else goto standard_seg 
        #define CALL_IF0(cond,f) if (cond) { (f); break; } else 
        #define CALL_WITH_FALLBACK(f) if (f(vb, ctx, STRa(value))) { seg_by_ctx (VB, STRa(value), ctx, value_len); } break
        #define STORE_AND_SEG(store_type) ({ seg_set_last_txt_store_value (VB, ctx, STRa(value), store_type); seg_by_ctx (VB, STRa(value), ctx, value_len); break; })
        #define DEFER(f,seg_after_did_i) ({ vb_add_to_deferred_q (VB, ctx, vcf_seg_INFO_##f, vb->idx_##f, seg_after_did_i); break; })
        // important: when adding DEFER, also update vcf_seg_copy_one_sample

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
        case _INFO_GNOMAD_AF:       CALL (vcf_seg_string (vb, ctx, STRa(value)));
        
        // ---------------------------------------
        // GATK fields
        // ---------------------------------------
        case _INFO_RAW_MQandDP:     CALL_IF (segconf.has[INFO_RAW_MQandDP], vcf_seg_INFO_RAW_MQandDP (vb, ctx, STRa(value)));
        case _INFO_BaseCounts:      DEFER(BaseCounts, DID_NONE); // depends on FORMAT_AD
        case _INFO_SF:              CALL_WITH_FALLBACK (vcf_seg_INFO_SF_init); // Source File
        case _INFO_MLEAC:           CALL (vcf_seg_INFO_MLEAC (vb, ctx, STRa(value)));
        case _INFO_MLEAF:           CALL (vcf_seg_INFO_MLEAF (vb, ctx, STRa(value)));
        case _INFO_QD:              CALL_IF (segconf.has[INFO_QD], DEFER (QD, INFO_DP)); // depends on VCF_QUAL and INFO_DP (that may depend on FORMAT_DP)
        case _INFO_RU:              CALL (vcf_seg_INFO_RU (vb, ctx, STRa(value)));
        case _INFO_RPA:             CALL (vcf_seg_INFO_RPA (vb, ctx, STRa(value)));
        case _INFO_MFRL:            
        case _INFO_MBQ:
        case _INFO_MMQ:             CALL (seg_array (VB, ctx, ctx->did_i, STRa(value), ',', 0, false, STORE_INT, DICT_ID_NONE, value_len));
        case _INFO_VDB:
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
        case _INFO_TLOD:            CALL (vcf_seg_string (vb, ctx, STRa(value)));
        
        case _INFO_VQSLOD:          CALL_IF0 (segconf.vcf_is_dragen, seg_textual_float (VB, ctx, STRa(value), value_len))
                                    CALL (vcf_seg_string (vb, ctx, STRa(value)));
        case _INFO_GERMQ:
        case _INFO_CONTQ:
        case _INFO_SEQQ:
        case _INFO_STRANDQ:
        case _INFO_STRQ:
        case _INFO_ECNT:            CALL (seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_AS_SB_TABLE:     CALL_IF (segconf.AS_SB_TABLE_by_SB, DEFER(AS_SB_TABLE, DID_NONE)); // depends on FORMAT_SB
        
        // ---------------------------------------
        // VEP fields
        // ---------------------------------------
        case _INFO_CSQ:
        case _INFO_vep:             CALL_IF (segconf.vcf_is_vep, vcf_seg_INFO_CSQ (vb, ctx, STRa(value)));
        case _INFO_AGE_HISTOGRAM_HET:
                                    CALL (seg_integer_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, VCF_SPECIAL_N_ALTS, N_ALTS, value_len));
        case _INFO_AGE_HISTOGRAM_HOM: 
                                    CALL (seg_integer_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, VCF_SPECIAL_N_ALTS, N_ALTS, value_len));
        case _INFO_DP_HIST:
        case _INFO_GQ_HIST:         CALL (seg_integer_matrix (VB, ctx, ctx->did_i, STRa(value), ',', '|', false, VCF_SPECIAL_N_ALLELES, N_ALTS + 1, value_len));
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
        case _INFO_MMURI:           CALL_IF (segconf.vcf_is_mastermind, seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_NONE, value_len)); // value is a URL
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
        case _INFO_CAF:             CALL_IF (segconf.vcf_is_dbSNP, vcf_seg_INFO_CAF (vb, ctx, STRa(value)));
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
        case _INFO_phyloP:          CALL_IF (segconf.vcf_is_isaac, vcf_seg_string (vb, ctx, STRa(value)));
        case _INFO_GMAF:            CALL_IF (segconf.vcf_is_isaac, vcf_seg_INFO_GMAF (vb, ctx, STRa(value)));
        case _INFO_EVS:             CALL_IF (segconf.vcf_is_isaac, vcf_seg_INFO_EVS (vb, ctx, STRa(value)));
        case _INFO_SNVHPOL:         CALL_IF (segconf.vcf_is_isaac, vcf_seg_INFO_SNVHPOL (vb, ctx, STRa(value)));

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
        case _INFO_TREE_SCORE:      CALL_IF (segconf.vcf_is_ultima, vcf_seg_string (vb, ctx, STRa(value)));
        case _INFO_VARIANT_TYPE:    CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_VARIANT_TYPE (vb, ctx, STRa(value)));
        case _INFO_HAPCOMP:
        case _INFO_ASSEMBLED_HAPS:  CALL_IF (segconf.vcf_is_ultima, seg_integer_or_not (VB, ctx, STRa(value), value_len));
        case _INFO_FILTERED_HAPS:   CALL_IF (segconf.vcf_is_ultima, vcf_seg_INFO_FILTERED_HAPS (vb, ctx, STRa(value)));
        
        // ---------------------------------------
        // freebayes
        // ---------------------------------------
        case _INFO_DPB:             CALL_IF (segconf.vcf_is_freebayes, DEFER(DPB, DID_NONE)); // depends on INFO_DP (that may depend on FORMAT_DP)

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
        // GIAB
        // ---------------------------------------
        case _INFO_callsetnames:    goto standard_seg; // better than seg_array
        case _INFO_callsets:        CALL_IF (has(callsetnames),  seg_by_ARRAY_LEN_OF (VB, ctx, STRa(value), STRa(BII(callsetnames)->value),  STRa(callsets_snip)));
        case _INFO_datasets:        CALL_IF (has(datasetnames),  seg_by_ARRAY_LEN_OF (VB, ctx, STRa(value), STRa(BII(datasetnames)->value),  STRa(datasets_snip)));
        case _INFO_platforms:       CALL_IF (has(platformnames), seg_by_ARRAY_LEN_OF (VB, ctx, STRa(value), STRa(BII(platformnames)->value), STRa(platforms_snip)));

        // ---------------------------------------
        // 1000 Genomes
        // ---------------------------------------
        case _INFO_ID:              CALL_IF (str_issame_(STRa(value), STRlst(VCF_ID)), seg_by_ctx (VB, STRa(copy_VCF_ID_snip), ctx, value_len));
        case _INFO_MAF:             CALL (vcf_seg_INFO_MAF (vb, ctx, STRa(value)));
        case _INFO_NS:              CALL (vcf_seg_INFO_NS (vb, ctx, STRa(value)));

        // ---------------------------------------
        // Structural variants
        // ---------------------------------------
        case _INFO_CIPOS:           CALL (vcf_seg_INFO_CIPOS (vb, ctx, STRa(value)));
        case _INFO_CIEND:           CALL (vcf_seg_INFO_CIEND (vb, ctx, STRa(value)));
        case _INFO_SVLEN:           CALL (vcf_seg_INFO_SVLEN (vb, ctx, STRa(value)));
        case _INFO_SVTYPE:          CALL (vcf_seg_SVTYPE (vb, ctx, STRa(value)));
        case _INFO_HOMSEQ:          CALL (vcf_seg_HOMSEQ (vb, ctx, STRa(value)));
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

        case _INFO_HOMLEN:          CALL_IF (segconf.vcf_is_sv    && has(HOMSEQ),      seg_LEN_OF (VB, ctx, STRa(value), BII(HOMSEQ)->value_len,      STRa(homlen_snip)));
        case _INFO_DUPHOMLEN:       CALL_IF (segconf.vcf_is_manta && has(DUPHOMSEQ),   seg_LEN_OF (VB, ctx, STRa(value), BII(DUPHOMSEQ)->value_len,   STRa(duphomlen_snip)));
        case _INFO_SVINSLEN:        CALL_IF (segconf.vcf_is_manta && has(SVINSSEQ),    seg_LEN_OF (VB, ctx, STRa(value), BII(SVINSSEQ)->value_len,    STRa(svinslen_snip)));
        case _INFO_DUPSVINSLEN:     CALL_IF (segconf.vcf_is_manta && has(DUPSVINSSEQ), seg_LEN_OF (VB, ctx, STRa(value), BII(DUPSVINSSEQ)->value_len, STRa(dupsvinslen_snip)));
        case _INFO_CIGAR:           CALL_IF (segconf.vcf_is_manta, vcf_seg_manta_CIGAR (vb, ctx, STRa(value)));
        case _INFO_BND_DEPTH:       CALL_IF (segconf.vcf_is_manta, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_BND_DEPTH, TW_MATE_BND_DEPTH, false, value_len));
        case _INFO_MATE_BND_DEPTH:  CALL_IF (segconf.vcf_is_manta, vcf_seg_sv_copy_mate (vb, ctx, STRa(value), TW_MATE_BND_DEPTH, TW_BND_DEPTH, false, value_len));

        // ---------------------------------------
        // Mobile elements
        // ---------------------------------------
        case _INFO_MEINFO:          CALL (vcf_seg_MEINFO (vb, ctx, STRa(value)));
        case _INFO_ADJLEFT:         CALL_IF (segconf.vcf_is_melt, vcf_seg_melt_ADJLEFT (vb, ctx, STRa(value)));
        case _INFO_ADJRIGHT:        CALL_IF (segconf.vcf_is_melt, vcf_seg_melt_ADJRIGHT (vb, ctx, STRa(value)));
        case _INFO_INTERNAL:        CALL_IF (segconf.vcf_is_melt, vcf_seg_melt_INTERNAL (vb, ctx, STRa(value)));
        case _INFO_DIFF:            CALL_IF (segconf.vcf_is_melt, vcf_seg_melt_DIFF (vb, ctx, STRa(value)));
        case _INFO_LP:
        case _INFO_RP:
        case _INFO_SR:
        case _INFO_ASSESS:          CALL_IF (segconf.vcf_is_melt, seg_integer_or_not (VB, ctx, STRa(value), value_len));

        case _INFO_VT:              CALL (vcf_seg_INFO_VT (vb, ctx, STRa(value)));

        // ---------------------------------------
        // cncb
        // ---------------------------------------
        case _INFO_MA:              CALL_WITH_FALLBACK (vcf_seg_INFO_MA);

        default: standard_seg:
            vcf_seg_field_fallback (vb, ctx, STRa(value));
            
            if (ctx->flags.store == STORE_INT) {
                int64_t val;
                if (str_get_int (STRa(value), &val))
                    ctx_set_last_value (VB, ctx, val);
            }
    }

    seg_set_last_txt (VB, ctx, STRa(value));
}

void vcf_seg_field_fallback (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    #define ISNUM(x)  (ctx->header_info.vcf.Number == (x))
    #define ISTYPE(x) (ctx->header_info.vcf.Type == VCF_##x)
    
    if (ISTYPE(Integer) && ISNUM(1)) 
        seg_integer_or_not (VB, ctx, STRa(value), value_len);

    else if (ISNUM(NUMBER_A)) 
        seg_array_(VB, ctx, ctx->did_i, STRa(value), ',', 0, false, 
                   ISTYPE(Integer) ? STORE_INT : STORE_NONE/*goes by ctx->seg_to_local*/, 
                   DICT_ID_NONE, VCF_SPECIAL_N_ALTS, N_ALTS, value_len);

    else if (ISTYPE(Float) && ISNUM(1) && ctx->seg_to_local) 
        seg_add_to_local_string (VB, ctx, STRa(value), LOOKUP_SIMPLE, value_len);

    else
        // note: in case R/G, >1: values are usally correlated, so we're better off with segging to dict
        seg_by_ctx (VB, STRa(value), ctx, value_len);
}

void vcf_parse_info_subfields (VBlockVCFP vb, STRp(info))
{
    // parse the info string
    str_split (info, info_len, MAX_FIELDS, ';', pair, false); 
    ASSVCF (n_pairs, "Too many INFO subfields, Genozip supports up to %u", MAX_FIELDS);

    buf_alloc (vb, &ii_buf, 0, n_pairs + 2, InfoItem, CTX_GROWTH, "info_items");

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
        
        if (segconf_running) segconf.has[ii.ctx->did_i]++;

        #define X(x) case INFO_##x : vb->idx_##x = ii_buf.len32; break
        switch (ii.ctx->did_i) {
            X(AN); X(AF); X(AC); X(MLEAC); X(MLEAF); X(AC_Hom); X(AC_Het); X(AC_Hemi); X(DP); X(QD); X(SF);
            X(AS_SB_TABLE); X(SVINSSEQ) ; X(SVTYPE); X(SVLEN); X(HOMSEQ); X(END) ; X(CIPOS); X(LEFT_SVINSSEQ);
            X(BaseCounts); X(platformnames); X(datasetnames); X(callsetnames); X(AF1000G); X(RefMinor);
            X(DPB); X(COMMON);
            default: {}
        }
        #undef X

        BNXT (InfoItem, ii_buf) = ii;
    }
}

void vcf_seg_info_subfields (VBlockVCFP vb, STRp(info))
{
    START_TIMER;

    // case: INFO field is '.' (empty) 
    if (IS_PERIOD (info) && !segconf.vcf_is_isaac) { // note: in Isaac, it slightly better to mux the "."
        seg_by_did (VB, ".", 1, VCF_INFO, 2); // + 1 for \t or \n
        return;
    }

    vcf_parse_info_subfields (vb, STRa(info));

    for_buf (InfoItem, ii, ii_buf)
        if (ii->value) 
            vcf_seg_info_one_subfield (vb, ii->ctx, STRa(ii->value));
        else
            ctx_set_encountered (VB, ii->ctx);

    set_last_txt (VCF_INFO, info);

    COPY_TIMER (vcf_seg_info_subfields);
}

// Seg INFO fields that were deferred to after all samples are segged
void vcf_seg_finalize_INFO_fields (VBlockVCFP vb)
{
    if (!ii_buf.len) return; // no INFO items on this line (except if dual-coords - we will add them in a sec)
    
    START_TIMER;
    decl_ctx (VCF_INFO);

    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true }; 
 
    con_set_nitems (con, ii_buf.len32);

    char prefixes[CONTAINER_MAX_PREFIXES_LEN];  // these are the Container prefixes
    prefixes[0] = prefixes[1] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Container 
    uint32_t total_names_len=0;
    for_buf2 (InfoItem, ii, i, ii_buf) {
        // Set the Container item and find (or create) a context for this name
        con.items[i] = (ContainerItem){ .dict_id   = !ii->value                        ? DICT_ID_NONE 
                                                   : ii->ctx->dict_id.num == _INFO_END ? (DictId)_VCF_POS
                                                   :                                     ii->ctx->dict_id,
                                        .separator = { ';' } }; 
            
        // add to the prefixes
        ASSVCF (prefixes_len + ii->name_len + 1 <= CONTAINER_MAX_PREFIXES_LEN, 
                "INFO contains tag names that, combined (including the '='), exceed the maximum of %u characters", CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii->name, ii->name_len);
        prefixes_len += ii->name_len;
        prefixes[prefixes_len++] = CON_PX_SEP;

        total_names_len += ii->name_len + 1; // +1 for ; \t or \n separator
    }

    // seg deferred fields 
    if (flag.debug_seg) vb_display_deferred_q (VB, __FUNCTION__);

    for (uint8_t i=0; i < vb->deferred_q_len; i++) // deferred_q_len=0 if is_copied
        vb->deferred_q[i].seg (VB);

    // case INFO is muxed: multiplex by has_RGQ or FILTER in Isaac
    ContextP channel_ctx;
    if (!segconf_running && (segconf.has[FORMAT_RGQ] || segconf.vcf_is_isaac)) {
        channel_ctx = seg_mux_get_channel_ctx (VB, VCF_INFO, (MultiplexerP)&vb->mux_INFO, (segconf.has[FORMAT_RGQ] ? CTX(FORMAT_RGQ)->line_has_RGQ : vcf_isaac_info_channel_i (VB)));
        seg_by_ctx (VB, STRa(vb->mux_INFO.snip), ctx, 0);
    }

    // case: INFO not muxed
    else 
        channel_ctx = ctx;

    container_seg (vb, channel_ctx, &con, prefixes, prefixes_len, total_names_len /* names inc. = and separator */);
    
    COPY_TIMER (vcf_seg_finalize_INFO_fields);
}
