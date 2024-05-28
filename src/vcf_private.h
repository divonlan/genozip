// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "vblock.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"

#define VCF_MAGIC "##fileformat=VCF"
#define BCF_MAGIC "BCF"

#define VCF_FIELD_NAMES "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
#define VCF_FIELD_NAMES_LONG VCF_FIELD_NAMES "\tFORMAT"

// note: whenever the same field is used (eg VCF_QUAL), it has to have the same tw
typedef enum          { TW_CHROM,  TW_MATE_CHROM,  TW_POS,  TW_MATE_POS,  TW_QUAL,  TW_SAMPLES,  TW_EVDNC,   TW_NUMPARTS,   TW_HOMSEQ,   TW_INSERTION,   TW_SCTG,   TW_NM,   TW_MATENM,   TW_MAPQ,   TW_MATEMAPQ, NUM_SVABA_TWs } MateCopiesSvaba;
#define SVABA_TW_DIDS { VCF_CHROM, VCF_MATE_CHROM, VCF_POS, VCF_MATE_POS, VCF_QUAL, VCF_SAMPLES, INFO_EVDNC, INFO_NUMPARTS, INFO_HOMSEQ, INFO_INSERTION, INFO_SCTG, INFO_NM, INFO_MATENM, INFO_MAPQ, INFO_MATEMAPQ              }

typedef enum          { TW_CHROM0, TW_MATE_CHROM0, TW_POS0, TW_MATE_POS0, TW_QUAL0, TW_SAMPLES0, TW_FILTER,  TW_BND_DEPTH,   TW_MATE_BND_DEPTH,  NUM_MANTA_TWs } MateCopiesManta;
#define MANTA_TW_DIDS { VCF_CHROM, VCF_MATE_CHROM, VCF_POS, VCF_MATE_POS, VCF_QUAL, VCF_SAMPLES, VCF_FILTER, INFO_BND_DEPTH, INFO_MATE_BND_DEPTH               }

typedef enum          { TW_CHROM1, TW_MATE_CHROM1, TW_POS1, TW_MATE_POS1, TW_unused1, TW_unused2, TW_FILTER0, NUM_PBSV_TWs } MateCopiesPbsv;
#define PBSV_TW_DIDS  { VCF_CHROM, VCF_MATE_CHROM, VCF_POS, VCF_MATE_POS, DID_NONE,   DID_NONE,   VCF_FILTER,              } 

#define NUM_TWs NUM_SVABA_TWs // the maximal one

typedef struct {
    WordIndex chrom;        // Seg: enter as node_index ; Merge: convert to word_index
    WordIndex format_node_i;// the node_index into contexts[VCF_FORMAT].nodes and also format_mapper_buf that applies to this line. Data on the fields is in format_mapper_buf[dl.format_node_i]
    PosType32 pos;          // 
    TxtWord BND_id;         // BND variants: a number as close as possible to unique of a BND event
    TxtWord tw[NUM_TWs];    // used by vcf_seg_sv_copy_mate 
} ZipDataLineVCF;

#define DATA_LINE(i) B(ZipDataLineVCF, vb->lines, i)

typedef enum { VCF_v_UNKNOWN, VCF_v4_1, VCF_v4_2, VCF_v4_3, VCF_v4_4, VCF_v4_5 } VcfVersion;

#define ZIP_MAX_PLOIDY_FOR_MUX 4 // ZIP only. In PIZ use z_file->max_ploidy_for_mux
#define ZIP_NUM_DOSAGES_FOR_MUX (ZIP_MAX_PLOIDY_FOR_MUX+1) // ZIP only: 0 to ZIP_MAX_PLOIDY_FOR_MUX
typedef MULTIPLEXER(ZIP_NUM_DOSAGES_FOR_MUX) DosageMultiplexer, *DosageMultiplexerP;

#define VCF_MAX_ARRAY_ITEMS SMALL_CON_NITEMS

typedef packed_enum { 
    VT_UNKNOWN, 
    VT_SNP,         // A G 
    VT_DEL,         // AGG A   - must be left anchored
    VT_INS,         // A AGG   - must be left anchored
    VT_SNP_LONG,    // AGG CGG 
    VT_DEL_LONG,    // AGG AG
    VT_INS_LONG,    // AG AGG 
    VT_SUBST_DEL,   // ATATGTG ATG - left anchored, part of remaining bases after deletion are substituted 
    VT_SUBST_INS,   // ATG ATATGTG - left anchored, part of ref is substituted, and then insertion 
    VT_SUBST,       // ATG ACG - left anchored
    VT_UPSTRM_DEL,  // A * - upstream deletion
    VT_BND,         // eg ]11:69541170]T or G[11:69485520[
    VT_SYM_DEL,     // <DEL> or <DEL:ME>
    VT_SYM_INS,     // <INS> or <INS:ME> or <INS:ME:ALU>, <INS:ME:LINE1>, <INS:ME:SVA>, <INS:UNK>
    VT_SYM_DUP,     // <DUP> or <DUP:TANDEM> Duplication
    VT_SYM_INV,     // <INV> Inversion
    VT_SYM_CNV,     // <CNV> Copy Number Polymorphism
    VT_SYM_CPX,     // <CPX> Complex SV - see https://gatk.broadinstitute.org/hc/en-us/articles/5334587352219-How-to-interpret-SV-VCFs
    VT_SYM_CTX,     // <CTX> Reciprocal chromosomal translocation - see https://gatk.broadinstitute.org/hc/en-us/articles/5334587352219-How-to-interpret-SV-VCFs
    VT_SYM_BND,     // <BND> Translocation - see https://gatk.broadinstitute.org/hc/en-us/articles/5334587352219-How-to-interpret-SV-VCFs
    VT_NO_ALT,      // A .
    NUM_VTs
    // Symbolic variants observed:
    // 1000 Genome conventions: https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
    // ##ALT=<ID=DEL,Description="Deletion">
    // ##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
    // ##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">
    // ##ALT=<ID=DUP,Description="Duplication">
    // ##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
    // ##ALT=<ID=INS,Description="Insertion of novel sequence">
    // ##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
    // ##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
    // ##ALT=<ID=INV,Description="Inversion">
    // ##ALT=<ID=CNV,Description="Copy number variable region">

    // GATK:
    // ALT=<ID=BND,Description="Translocation">
    // ALT=<ID=CPX,Description="Complex SV">
    // ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">
    // ALT=<ID=DEL,Description="Deletion">
    // ALT=<ID=DUP,Description="Duplication">
    // ALT=<ID=INS,Description="Insertion">
    // ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">
    // ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">
    // ALT=<ID=INS:ME:SVA,Description="SVA element insertion">
    // ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">

    // https://forgemia.inra.fr/genotoul-bioinfo/jflow-toolshed/blob/b1c02354bad4e40cb429edb90356d7516e81a8eb/cnv_detection/lib/svreader/resources/template.vcf
    // ##ALT=<ID=IDP,Description="Interspersed Duplication">
    // ##ALT=<ID=ITX,Description="Intra-chromosomal translocation">

} VariantType;
#define ALT0(x) (VB_VCF->var_types[0] == VT_##x)

#define SVTYPE_BY_VT {                                              \
    [VT_BND]="BND",                                                 \
    [VT_DEL]="DEL",     [VT_SYM_DEL]="DEL", [VT_SUBST_DEL]="DEL",   \
    [VT_INS]="INS",     [VT_SYM_INS]="INS", [VT_SUBST_INS]="INS",   \
    [VT_SYM_INV]="INV", [VT_SYM_DUP]="DUP", [VT_SYM_BND]="BND",     \
    [VT_SYM_CPX]="CPX", [VT_SYM_CTX]="CTX", [VT_SYM_CNV]="CNV", /* note: in pbsv this is "cnv" */ \
    /* [VT_UPSTRM_SEL]=? not encountered yet */                     \
}

typedef struct VBlockVCF {
    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    Ploidy ploidy;           // ZIP
    VcfVersion vcf_version;
    uint64_t first_line;     // ZIP: used for --add_line_numbers  
    
    // This line's REFALT
    rom REF, ALT;            // Seg: pointer into txt_data of REF in main field, set by vcf_seg_txt_line
    uint32_t REF_len, ALT_len;
    int8_t n_alts;           // generated by vb_parse_ALT: 0 means not parsed yet, -1 means parse fails (eg bc too many alts)
    rom alts[MAX_ALLELES-1];
    uint32_t alt_lens[MAX_ALLELES-1];
    VariantType var_types[MAX_ALLELES-1]; // of main ALT[] fields

    // BND type 0, eg: REF=A ALT="AAACTCCT[hs37d5:33588521["
    // BND type 1: eg: REF=A ALT="AAACTCCT]hs37d5:33588521]"
    // BND type 2: REF=G ALT="[hs37d5:35428323[TAAGAGCCGCTGGCTGGCTGTCCGGGCAGGCCTCCTGGCTGCACCTGCCACAGTGCACAGGCTGACTGAGGTGCACG"
    // BND type 3: REF=G ALT="]hs37d5:35428323]TAAGAGCCGCTGGCTGGCTGTCCGGGCAGGCCTCCTGGCTGCACCTGCCACAGTGCACAGGCTGACTGAGGTGCACG"
    uint8_t BND_type;        // type of BND variant (0 to 3) - valid only if var_types[0]==VT_BND 
    
    TxtWord BND_INS;         // insertion (if any) included in the BND variant

    // structural variants
    STR(mate_chrom_name);    // ZIP/PIZ: mate's chrom in case of VT_BND (valid only if n_alts>0 and var_tyeps[0]==VT_BND)
    PosType32 mate_pos;      // ZIP/PIZ: mate's pos in case of VT_BND
    LineIType mate_line_i;   // Seg/PIZ: the mate of this line.
    uint32_t mate_line_count;// for stats
    Multiplexer2 mux_mate[NUM_TWs]; // BND variants: mux by mate 
    Multiplexer6 mux_pbsv_I0D;
    Multiplexer3 mux_pbsv_I1D;

    // FORMAT/AD
    int64_t ad_values[VCF_MAX_ARRAY_ITEMS];
    
    #define first_idx idx_AN        // ZIP: INFO fields indices within INFO
    // IMPORTANT: when adding, add to X() in vcf_seg_info_subfields
    int16_t idx_AN, idx_AC, idx_AF, idx_MLEAC, idx_MLEAF, idx_AC_Hom, idx_AC_Het, idx_AC_Hemi, idx_QD, idx_DP, idx_SF, 
            idx_AS_SB_TABLE, idx_END, idx_SVLEN, idx_CIPOS, idx_BaseCounts,
            idx_SVTYPE, idx_HOMSEQ, idx_DUPHOMSEQ, idx_SVINSSEQ, idx_DUPSVINSSEQ, idx_LEFT_SVINSSEQ,
            idx_platformnames, idx_datasetnames, idx_callsetnames, idx_AF1000G; 

    #define has(f)   (vb->idx_##f != -1)
    #define after_idx mux_PLn

    // Multiplexers
    #define first_mux mux_PLn
    DosageMultiplexer mux_PLn, mux_GL, mux_GP, mux_PRI, mux_DS, mux_PP, mux_PVAL, mux_FREQ, mux_RD, 
                      mux_VAF, mux_AD[2], mux_ADALL[2];
    
    MULTIPLEXER(51 * ZIP_NUM_DOSAGES_FOR_MUX) mux_PLy; // TODO: 60 would be better than 51 as it was up to 15.0.35, but mux is currently limited to 256 channels
    MULTIPLEXER(7 * ZIP_NUM_DOSAGES_FOR_MUX) mux_GQ;
    #define MAX_DP_FOR_MUX 60       
    MULTIPLEXER(MAX_DP_FOR_MUX) mux_RGQ;   

    Multiplexer2 mux_POS;           // GVCF: multiplex by whether this field is END or POS
    Multiplexer2 mux_QUAL, mux_INFO;// GVCF: multiplex by has_RGQ
    Multiplexer2 mux_FORMAT_DP;     // channel=1 by AD or SDP, channel=0 transposed if AD/SDP not available           
    Multiplexer2 mux_IGT, mux_IPS;  // multiplex by (sample_i>0)
    Multiplexer3 mux_VC;            // multiplex dbSNP's INFO/VC by VARTYPE
    Multiplexer3 mux_GQX;           // multiplex Isaac's FORMAT/GQX
    Multiplexer3 mux_BAF, mux_X, mux_Y; // Illumina genotyping: by adjusted dosage 

    #define after_mux PL_mux_by_DP

    thool PL_mux_by_DP;
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;
#define VB_VCF ((VBlockVCFP)vb)

typedef ContextP ContextPBlock[MAX_FIELDS];

// VCF standard keys
#define HK_GENOZIP_CMD   "##genozip_command="
#define TAG_SOURCE       Source=\""GENOZIP_URL"\""
#define KH_INFO          "##INFO=<ID="
#define HK_CONTIG        "##contig="

extern uint32_t vcf_num_samples; // ZIP
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag.samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples
extern VcfVersion vcf_header_get_version (void);

#define BII(x) B(InfoItem, vb->contexts[VCF_INFO].info_items, vb->idx_##x)

extern void vcf_seg_array_of_N_ALTS_numbers (VBlockVCFP vb, ContextP ctx, STRp(value), StoreType type);

// POS stuff
extern void vcf_seg_pos (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(pos_str));
extern void vcf_seg_INFO_END (VBlockVCFP vb, ContextP end_ctx, STRp(end_str));

// QUAL stuff
extern void vcf_seg_QUAL (VBlockVCFP vb, STRp(qual));
extern void vcf_segconf_finalize_QUAL (VBlockVCFP vb);
extern void vcf_piz_insert_QUAL_by_GP (VBlockVCFP vb);

// ID stuff
extern void vcf_piz_insert_ID_is_variant (VBlockVCFP vb);

// AC / AF / AN
extern void vcf_seg_INFO_AC (VBlockVCFP vb, ContextP ac_ctx, STRp(ac_str));
extern void vcf_seg_INFO_AN (VBlockVCFP vb);
extern void vcf_seg_INFO_MLEAC (VBlockVCFP vb, ContextP ac_ctx, STRp(ac_str));
extern void vcf_seg_INFO_MLEAF (VBlockVCFP vb, ContextP ctx, STRp(mleaf));
extern void vcf_piz_insert_INFO_AN (VBlockVCFP vb);

// 1000G stuff
extern void vcf_1000G_zip_initialize (void);
extern void vcf_1000G_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_MAF (VBlockVCFP vb, ContextP ctx, STRp(maf));
extern void vcf_seg_INFO_NS (VBlockVCFP vb, ContextP ctx, STRp(ns_str));

// Samples stuff
extern void vcf_seg_FORMAT (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(fmt));
extern void vcf_samples_zip_initialize (void);
extern void vcf_samples_seg_initialize (VBlockVCFP vb);
extern void vcf_samples_seg_finalize (VBlockVCFP vb);

extern rom vcf_seg_samples (VBlockVCFP vb, ZipDataLineVCF *dl, int32_t len, char *next_field, bool *has_13);
extern int vcf_seg_get_mux_channel_i (VBlockVCFP vb);
extern int vcf_piz_get_mux_channel_i (VBlockP vb);
extern ContextP vcf_seg_FORMAT_mux_by_dosage (VBlockVCFP vb, ContextP ctx, STRp(cell), const DosageMultiplexer *mux);
extern void vcf_seg_FORMAT_mux_by_dosagexDP (VBlockVCFP vb, ContextP ctx, STRp(cell), void *mux_p);

// FORMAT/GT stuff
extern WordIndex vcf_seg_FORMAT_GT (VBlockVCFP vb, ContextP ctx, ZipDataLineVCF *dl, STRp(cell), bool has_ps, bool has_null_dp);
extern void vcf_seg_FORMAT_GT_finalize_line (VBlockVCFP vb, uint32_t line_n_samples);
extern void vcf_piz_FORMAT_GT_rewrite_predicted_phase (VBlockP vb, char *recon, uint32_t recon_len);
extern void vcf_piz_GT_cb_null_GT_if_null_DP (VBlockP vb , char *recon);
extern int vcf_piz_GT_get_last_dosage (VBlockP vb);

static inline Allele *this_sample_GT (VBlockVCFP vb) {
    ContextP ht_ctx = CTX(FORMAT_GT_HT);
    return B(Allele, ht_ctx->local, ht_ctx->HT_n_lines * ht_ctx->ht_per_line + vb->ploidy * vb->sample_i);
}

// GIAB trio stuff
extern void vcf_giab_zip_initialize (void);
extern void vcf_giab_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_FORMAT_IGT (VBlockVCFP vb, ContextP ctx, STRp(igt));
extern void vcf_seg_FORMAT_IPS (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP ctx, STRp(ips));
extern void vcf_seg_ADALL_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values);
eSTRl(datasets_snip); eSTRl(callsets_snip); eSTRl(platforms_snip);

#define IS_TRIVAL_FORMAT_SUBFIELD ((!recon_len || (recon_len==1 && *recon=='.')) && dict_id_is_vcf_format_sf (ctx->dict_id))
extern void vcf_FORMAT_PL_decide (VBlockVCFP vb);
extern void vcf_FORMAT_PL_after_vbs (void);

// FORMAT/PS and FORMAT/PID stuff
extern void vcf_samples_zip_initialize_PS_PID (void);
extern void vcf_samples_seg_initialize_LOOKBACK (VBlockVCFP vb);
extern void vcf_samples_seg_finalize_PS_PID (VBlockVCFP vb);
extern void vcf_seg_FORMAT_PS_PID (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP ctx, STRp(value));
extern void vcf_seg_FORMAT_PS_PID_missing_value (VBlockVCFP vb, ContextP ctx, rom end_of_sample);
extern void vcf_samples_seg_initialize_PS_PID (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_piz_ps_pid_lookback_insert (VBlockP vb, Did did_i, STRp(recon)); 
extern void vcf_piz_ps_pid_lookback_shift (VBlockP vb, STRp(insert));

// FORMAT/GQ
extern void vcf_segconf_finalize_GQ (VBlockVCFP vb);
extern void vcf_seg_FORMAT_GQ (VBlockVCFP vb);

// GATK stuff
extern void vcf_gatk_zip_initialize (void);
extern void vcf_gatk_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_RAW_MQandDP (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_INFO_BaseCounts (VBlockP vb);
extern void vcf_seg_INFO_AS_SB_TABLE (VBlockP vb);
extern void vcf_piz_sum_SB_for_AS_SB_TABLE (VBlockP vb, STRp(recon));
extern void vcf_piz_insert_INFO_AS_SB_TABLE (VBlockVCFP vb);
extern void vcf_seg_INFO_RU (VBlockVCFP vb, ContextP ctx, STRp(ru));
extern void vcf_seg_INFO_RPA (VBlockVCFP vb, ContextP ctx, STRp(rpa_str));

// INFO/SF
extern bool vcf_seg_INFO_SF_init (VBlockVCFP vb, ContextP sf_ctx, STRp(value));
extern void vcf_seg_INFO_SF_one_sample (VBlockVCFP vb);
extern void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *recon, int32_t recon_len);
extern void vcf_piz_insert_INFO_SF (VBlockVCFP vcf_vb);

// INFO/QD stuff
typedef packed_enum { QD_PRED_NONE, QD_PRED_INFO_DP, QD_PRED_INFO_DP_P001, QD_PRED_INFO_DP_M001, 
                      QD_PRED_SUM_DP,  QD_PRED_SUM_DP_P001,  QD_PRED_SUM_DP_M001, NUM_QD_PRED_TYPES } QdPredType;
extern void vcf_seg_sum_DP_for_QD (VBlockVCFP vb, int64_t value);
extern void vcf_seg_INFO_QD (VBlockP vb);
extern void vcf_piz_sum_DP_for_QD (VBlockP vb, STRp(recon));
extern void vcf_piz_insert_INFO_QD (VBlockVCFP vb);
extern void vcf_piz_insert_INFO_BaseCounts_by_AD (VBlockVCFP vb);

// INFO stuff

typedef struct { char name[MAX_TAG_LEN]; // not nul-terminated, including '=' if there is one
                 rom value; 
                 unsigned name_len, value_len; 
                 ContextP ctx; } InfoItem, *InfoItemP;

extern void vcf_info_zip_initialize (void);
extern void vcf_info_seg_initialize (VBlockVCFP vb);
extern void vcf_piz_insert_INFO_DP (VBlockVCFP vb);

extern void vcf_seg_info_subfields (VBlockVCFP vb, STRp(info));
extern void vcf_seg_finalize_INFO_fields (VBlockVCFP vb);
extern bool vcf_seg_INFO_allele (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat);

// REFALT stuff
extern void vcf_refalt_zip_initialize (void);
extern void vcf_refalt_seg_initialize (VBlockVCFP vb);
extern void vcf_refalt_seg_REF_ALT (VBlockVCFP vb, STRp(ref), STRp(alt));
typedef enum { EQUALS_NEITHER, EQUALS_REF, EQUALS_ALT, EQUALS_MISSING } RefAltEquals;
extern bool vcf_refalt_piz_is_variant_snp (VBlockVCFP vb);
extern bool vcf_refalt_piz_is_variant_indel (VBlockVCFP vb);
extern void vcf_piz_refalt_parse (VBlockVCFP vb);
extern void vb_parse_ALT (VBlockVCFP vb);

// GVCF stuff
extern bool vcf_piz_line_has_RGQ (VBlockVCFP vb);

// Illumina genotyping stuff
extern void vcf_illum_gtyping_zip_initialize (void);
extern void vcf_illum_gtyping_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_PROBE_A (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_PROBE_B (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_ILLUMINA_CHR (VBlockVCFP vb, ContextP ctx, STRp(chr));
extern void vcf_seg_ILLUMINA_POS (VBlockVCFP vb, ContextP ctx, STRp(pos));
extern void vcf_seg_ILLUMINA_STRAND (VBlockVCFP vb, ContextP ctx, STRp(strand));
extern void vcf_seg_ALLELE_A (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_ALLELE_B (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_mux_by_adjusted_dosage (VBlockVCFP vb, ContextP ctx, STRp(baf), const Multiplexer3 *mux);

// dbSNP
extern void vcf_dbsnp_zip_initialize (void);
extern void vcf_dbsnp_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_RS (VBlockVCFP vb, ContextP ctx, STRp(rs));
extern void vcf_seg_INFO_RSPOS (VBlockVCFP vb, ContextP ctx, STRp(rspos));
extern void vcf_seg_INFO_VC (VBlockVCFP vb, ContextP ctx, STRp(vc));

// GWAS-VCF stuff
extern void vcf_gwas_zip_initialize (void);
extern void vcf_gwas_seg_initialize (VBlockVCFP vb);
extern void vcf_gwas_seg_FORMAT_ID (VBlockVCFP vb, ContextP ctx, STRp(id));

// VAGrENT stuff
extern void vcf_vagrent_zip_initialize (void);
extern void vcf_seg_INFO_VD (VBlockVCFP vb, ContextP ctx, STRp(vd));
extern void vcf_seg_INFO_VW (VBlockVCFP vb, ContextP ctx, STRp(vw));

// VEP stuff
extern uint64_t _INFO_Allele;
extern void vcf_vep_zip_initialize (void);
extern void vcf_vep_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_CSQ (VBlockVCFP vb, ContextP ctx, STRp(value));

// SnpEff stuff
extern void vcf_seg_INFO_ANN (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_INFO_EFF (VBlockVCFP vb, ContextP ctx, STRp(value));

// COSMIC stuff
extern void vcf_cosmic_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_LEGACY_ID (VBlockVCFP vb, ContextP ctx, STRp(lid));
extern void vcf_seg_INFO_SO_TERM (VBlockVCFP vb, ContextP ctx, STRp(st));

// ClinVar stuff
extern void vcf_seg_hgvs_consolidate_stats (VBlockVCFP vb, Did parent);
extern bool vcf_seg_INFO_HGVS (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat);

// Mastermind stuff
extern void vcf_mastermind_zip_initialize (void);
extern void vcf_mastermind_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_mastermind_HGVSG (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_INFO_MMID3 (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_INFO_MMURI3 (VBlockVCFP vb, ContextP ctx, STRp(value));

// gnomAD stuff
extern void vcf_gnomad_zip_initialize (void);
extern void vcf_gnomad_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_VRS_Starts (VBlockVCFP vb, ContextP ctx, STRp(arr));
extern void vcf_seg_INFO_VRS_Ends (VBlockVCFP vb, ContextP ctx, STRp(arr));
extern void vcf_seg_INFO_VRS_States (VBlockVCFP vb, ContextP ctx, STRp(arr));
extern void vcf_seg_INFO_VRS_Allele_IDs (VBlockVCFP vb, ContextP ctx, STRp(ids));
extern void vcf_seg_INFO_AS_QD (VBlockVCFP vb, ContextP ctx, STRp(as_qd));
extern void vcf_seg_INFO_AS_SOR (VBlockVCFP vb, ContextP ctx, STRp(as_sor));
extern void vcf_seg_INFO_AS_MQ (VBlockVCFP vb, ContextP ctx, STRp(as_mq));
extern void vcf_seg_INFO_AS_MQRankSum (VBlockVCFP vb, ContextP ctx, STRp(as_mqranksum));
extern void vcf_seg_INFO_AS_FS (VBlockVCFP vb, ContextP ctx, STRp(as_fs));
extern void vcf_seg_INFO_AS_VarDP (VBlockVCFP vb, ContextP ctx, STRp(as_vardp));
extern void vcf_seg_INFO_AS_QUALapprox (VBlockVCFP vb, ContextP ctx, STRp(as_qualapprox));
extern void vcf_seg_INFO_AS_ReadPosRankSum (VBlockVCFP vb, ContextP ctx, STRp(as_readposranksum));

// ICGC stuff
extern void vcf_seg_INFO_mutation (VBlockVCFP vb, ContextP ctx, STRp(mut));

// ISAAC stuff
extern void vcf_isaac_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_FORMAT_GQX (VBlockVCFP vb, ContextP ctx, STRp(gqx));
extern void vcf_seg_INFO_GMAF (VBlockVCFP vb, ContextP ctx, STRp(gmaf));
extern void vcf_seg_INFO_EVS (VBlockVCFP vb, ContextP ctx, STRp(evs));
extern void vcf_seg_INFO_IDREP (VBlockVCFP vb, ContextP ctx, STRp(idrep));
extern int vcf_isaac_info_channel_i (VBlockP vb);

// structural variants stuff
extern rom svtype_by_vt[];
extern void vcf_sv_zip_initialize (Did *tw_dids, int num_tw_dids);
extern void vcf_sv_seg_initialize (VBlockVCFP vb, Did *tw_dids, int num_tw_dids);
extern void vcf_seg_SVTYPE (VBlockVCFP vb, ContextP ctx, STRp(svtype));
extern void vcf_seg_INFO_SVLEN (VBlockVCFP vb, ContextP ctx, STRp(svlen_str));
extern void vcf_seg_INFO_REFLEN (VBlockVCFP vb, ContextP ctx, STRp(reflen_str));
extern void vcf_seg_INFO_CIPOS (VBlockVCFP vb, ContextP ctx, STRp(cipos));
extern void vcf_seg_INFO_CIEND (VBlockVCFP vb, ContextP ctx, STRp(ciend));
extern void vcf_seg_HOMSEQ (VBlockVCFP vb, ContextP ctx, STRp(homseq));
extern void vcf_seg_BND_mate (VBlockVCFP vb, STRp(id), STRp(mate_id), uint64_t hash);
extern ContextP vcf_seg_sv_SAMPLES (VBlockVCFP vb, rom samples, uint32_t remaining_txt_len, ContextP *ctxs, uint32_t n_ctxs);
extern ContextP vcf_seg_sv_copy_mate (VBlockVCFP vb, ContextP ctx, STRp(value), int tw, int her_tw, bool seg_only_if_mated, unsigned add_bytes);

// manta stuff
eSTRl(homlen_snip); eSTRl(duphomlen_snip); eSTRl(svinslen_snip); eSTRl(dupsvinslen_snip);

extern void vcf_manta_zip_initialize (void);
extern void vcf_manta_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_manta_ID (VBlockVCFP vb, STRp(id));
extern void vcf_seg_manta_CIGAR (VBlockVCFP vb, ContextP ctx, STRp(cigar));
extern void vcf_seg_LEN_OF (VBlockVCFP vb, ContextP ctx, STRp(len_str), int16_t idx, STRp(special_snip));

// SvABA stuff
#define vcf_has_mate (VB_VCF->mate_line_i != NO_LINE)
extern void vcf_svaba_zip_initialize (void);
extern void vcf_svaba_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_svaba_ID (VBlockVCFP vb, STRp(id));
extern void vcf_seg_svaba_MATEID (VBlockVCFP vb, ContextP ctx, STRp(mate_id));
extern void vcf_seg_svaba_MAPQ (VBlockVCFP vb, ContextP ctx, STRp(mapq));
extern void vcf_seg_svaba_SPAN (VBlockVCFP vb, ContextP ctx, STRp(span_str));

// PBSV stuff
extern void vcf_pbsv_zip_initialize (void);
extern void vcf_pbsv_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_pbsv_ID (VBlockVCFP vb, STRp(id));
extern void vcf_piz_insert_pbsv_ID (VBlockVCFP vb);
extern void vcf_seg_pbsv_MATEID (VBlockVCFP vb, ContextP ctx, STRp(mate_id));

// Mobile Elements stuff
extern void vcf_me_zip_initialize (void);
extern void vcf_me_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_MEINFO (VBlockVCFP vb, ContextP ctx, STRp(meinfo));
extern void vcf_seg_melt_ADJLEFT (VBlockVCFP vb, ContextP ctx, STRp(adjleft_str));
extern void vcf_seg_melt_ADJRIGHT (VBlockVCFP vb, ContextP ctx, STRp(adjright_str));
extern void vcf_seg_melt_INTERNAL (VBlockVCFP vb, ContextP ctx, STRp(internal));
extern void vcf_seg_melt_DIFF (VBlockVCFP vb, ContextP ctx, STRp(diff));

// Ultima Genomics stuff
extern void vcf_ultima_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_INFO_X_LM_RM (VBlockVCFP vb, ContextP ctx, STRp(seq));
extern void vcf_seg_INFO_X_IL (VBlockVCFP vb, ContextP ctx, STRp(il_str));
extern void vcf_seg_INFO_X_IC (VBlockVCFP vb, ContextP ctx, STRp(ic_str));
extern void vcf_seg_INFO_X_HIN (VBlockVCFP vb, ContextP ctx, STRp(hin_str));
extern void vcf_seg_INFO_X_HIL (VBlockVCFP vb, ContextP ctx, STRp(hil_str));
extern void vcf_seg_INFO_FILTERED_HAPS (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_INFO_VARIANT_TYPE (VBlockVCFP vb, ContextP ctx, STRp(vt));
extern void vcf_seg_INFO_X_GCC (VBlockVCFP vb, ContextP ctx, STRp(x));

// Platypus stuff
extern void vcf_platypus_zip_initialize (void);
extern void vcf_platypus_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_playpus_INFO_SC (VBlockVCFP vb, ContextP ctx, STRp(seq));
extern void vcf_seg_playpus_INFO_TCR (VBlockVCFP vb, ContextP ctx, STRp(tcr_str));
extern void vcf_seg_playpus_INFO_WS_WE (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_playpus_INFO_HP (VBlockVCFP vb, ContextP ctx, STRp(hp_str));    
extern void vcf_seg_platypus_FORMAT_GOF (VBlockVCFP vb, ContextP ctx, STRp(gof_str));

eSTRl(snip_copy_af);
eSTRl(copy_VCF_POS_snip);
eSTRl(copy_VCF_ID_snip);
eSTRl(copy_INFO_AF_snip);

#define VCF_ERR_PREFIX { progress_newline(); fprintf (stderr, "Error %s:%u in variant %.*s:%"PRId64": ", __FUNCLINE, STRf(vb->chrom_name), vb->last_int (VCF_POS)); }
#define ASSVCF(condition, format, ...) ({ if (!(condition)) { VCF_ERR_PREFIX; fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }})
#define ASSVCF0(condition, msg)        ASSVCF ((condition), msg "%s", "")
#define WARNVCF(format, ...)           ({ if (!flag.quiet)  { VCF_ERR_PREFIX; fprintf (stderr, format "\n", __VA_ARGS__); } })

// inserting INFO fields after all Samples have been reconstructed
extern void vcf_piz_insert_field (VBlockVCFP vb, ContextP ctx, STRp(value), int chars_reserved);

#define vcf_piz_defer_to_later(x) ({    \
    ctx->special_res = SPEC_RES_DEFERRED;       \
    if (reconstruct) {                          \
        ctx->recon_insertion = vb->line_i + 1;  \
        Ltxt += segconf.wid_##x.width; /* our best guess - minimize memory moves during vcf_piz_insert_field */ \
    }                                           \
})

#define vcf_piz_defer_to_after_samples(x) vcf_piz_defer_to_later(x)

#define IS_RECON_INSERTION(ctx) ((ctx)->recon_insertion == vb->line_i + 1)