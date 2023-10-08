// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "vblock.h"
#include "vcf.h"
#include "website.h"
#include "seg.h"

#define VCF_FIELD_NAMES "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
#define VCF_FIELD_NAMES_LONG VCF_FIELD_NAMES "\tFORMAT"

typedef struct {
    WordIndex chrom[2];          // Seg: enter as node_index ; Merge: convert to word_index
    PosType32 pos[2];            // arrays of [2] - { primary-coord, luft-coord } 
    PosType32 end_delta;         // Delta of INFO/END vs POS (same in both coordinates) - used in case chrom and pos are the same
    uint32_t tie_breaker;        // tie-breaker in case chrom, pos and end are the same
    
    bool has_haplotype_data : 1; // FORMAT field contains GT
    bool has_genotype_data  : 1; // FORMAT field contains subfields other than GT

    WordIndex format_node_i;     // the node_index into contexts[VCF_FORMAT].nodes and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_node_i]
} ZipDataLineVCF;
#define DATA_LINE(i) B(ZipDataLineVCF, vb->lines, i)

typedef enum { VCF_v_UNKNOWN, VCF_v4_1, VCF_v4_2, VCF_v4_3, VCF_v4_4, VCF_v4_5 } VcfVersion;

typedef MULTIPLEXER(4) DosageMultiplexer;

#define VCF_MAX_ARRAY_ITEMS SMALL_CON_NITEMS

#define NUM_COORDS 4
typedef enum __attribute__ ((__packed__)) { DC_NONE, DC_PRIMARY, DC_LUFT, DC_BOTH } Coords; // 2 bits, part of the file format, goes into FlagsTxtHeader, FlagsVbHeader

#define OTHER_COORDS(c) ((c)==DC_PRIMARY ? DC_LUFT : DC_PRIMARY)
#define SEL(prim,luft) ((vb->line_coords == DC_PRIMARY) ? (prim) : (luft))

typedef struct VBlockVCF {
    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    Ploidy ploidy;           // ZIP only
    VcfVersion vcf_version;
    uint64_t first_line;     // ZIP: used for --add_line_numbers  
    
    // used for segging INFO
    Buffer info_items;       // Seg: INFO items of the line being segged

    rom main_ref;            // Seg: pointer into txt_data of REF in main field, set by vcf_seg_txt_line
    rom main_alt;

    unsigned main_ref_len, main_alt_len;

    // GVCF stuff
    enum { RGQ_UNKNOWN=-1, RGQ_HASNT=false, RGQ_HAS=true } line_has_RGQ;

    // INFO/SF stuff
    enum { USE_SF_UNKNOWN, USE_SF_YES, USE_SF_NO } use_special_sf;
    Buffer sf_txt, sf_snip; // INFO/SF data as it appears in the snip being constructed
        
    // FORMAT/AD
    int64_t ad_values[VCF_MAX_ARRAY_ITEMS];
    
    // FORMAT stuff 
    Buffer format_mapper_buf;       // ZIP only: an array of type Container - one entry per entry in CTX(VCF_FORMAT)->nodes   
    Buffer format_contexts;         // ZIP only: an array of format_mapper_buf.len of ContextBlock
    Buffer last_format;             // ZIP only: cache previous line's FORMAT string

    // Multiplexers
    #define first_mux mux_PLn
    DosageMultiplexer mux_PLn, mux_GL, mux_GP, mux_PRI, mux_DS, mux_PP, mux_PVAL, mux_FREQ, mux_RD, 
                      mux_BAF, mux_X, mux_Y, mux_VAF,
                      mux_AD[2], mux_ADALL[2];

    PLMuxByDP PL_mux_by_DP;
    
    #define MAX_DP_FOR_MUX 60
    MULTIPLEXER(1 + MAX_DP_FOR_MUX * 3) mux_PLy;
    MULTIPLEXER(1 + 7 * 3) mux_GQ;
    MULTIPLEXER(MAX_DP_FOR_MUX) mux_RGQ;   

    MULTIPLEXER(2) mux_QUAL, mux_INFO; // multiplex by has_RGQ (in GVCF)
    MULTIPLEXER(2) mux_IGT, mux_IPS;   // multiplex by (sample_i>0)
    MULTIPLEXER(3) mux_VC;             // multiplex dbSNP's INFO/VC by VARTYPE
    MULTIPLEXER(3) mux_GQX;            // multiplex Isaac's FORMAT/GQX

    #define after_mux hapmat_helper_index_buf
    // used by CODEC_HAPM (for VCF haplotype matrix) 
    Buffer hapmat_helper_index_buf; // ZIP: used by codec_hapmat_count_alt_alleles 
    Buffer hapmat_columns_data;     // used by codec_hapmat_piz_get_one_line 
    Buffer hapmat_column_of_zeros;  // used by codec_hapmat_piz_calculate_columns   
    Buffer hapmat_one_array;        // one line or column 

    // DVCF stuff
    bool sort;                      // ZIP: true if this VB will be sorted
    Coords vb_coords;               // ZIP: DC_PRIMARY, DC_LUFT or DC_BOTH
                                    // PIZ: DC_PRIMARY or DC_LUFT - influenced by FlagsVbHeader.coords and flag.luft 
    bool is_rejects_vb;             // PIZ/ZIP: this is a VB of rejects variants for header ##primary_only/##luft_only
    Coords line_coords;             // Seg: coords of current line - DC_PRIMARY or DC_LUFT
    uint32_t pos_aln_i;             // ZIP: chain alignment of POS (used to compare to that of END)
    int32_t recon_size_luft;        // ZIP only: expected reconstruction size if this VB is reconstructed in LUFT coords inc. as ##luft_only in a DC_LUFT rejects VB) 

    Buffer tags;                    // Seg: used for FORMAT and INFO tag renaming.
    Buffer rejects_report;          // human readable report about rejects
    char new_ref;                   // SNP: new REF that is neither REF nor ALT; left-anchored INDEL with XSTRAND: new anchor
    bool is_del_sv;                 // is ALT == "<DEL>"
    bool prev_line_rejected;        // ZIP only
    int32_t reject_bytes;           // ZIP of a Luft file: number of bytes of reject data in this VB (data originating from ##primary_only/##luft_only) 
    bool is_unsorted[2];            // ZIP: line order of this VB[primary, luft] is unsorted 
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;
#define VB_VCF ((VBlockVCFP)vb)

typedef ContextP ContextPBlock[MAX_FIELDS];

// Liftover stuff

// The file format contains the string values defined in DVCF_STATUS_NAMES. The numeric values or order can be changed.
typedef enum {
    // Algorithm interim statuses not outputed to the genozip file
    LO_NA = -1,                     // not a Liftover item - never written z_file
    LO_UNKNOWN = 0,                 // Status of this Liftover line not known yet - never written z_file

    // OK statuses - this are internal to genozip and z_file.oStatus - not expressed in the txt file
    LO_OK,                          // internal progress step in ZIP - never written z_file
    LO_OK_REF_SAME_SNP,             // REF and ALT are the same for a SNP
    LO_OK_REF_SAME_SNP_REV,         // REF and ALT are the same for a SNP - reverse complemented
    LO_OK_REF_SAME_SNP_IUPAC,       // REF considered unchanged as it matches a IUPAC \"base\" in the Luft reference
    LO_OK_REF_SAME_INDEL,           // REF and ALT are the same for an INDEL, and confirmed not to be a REF⇆ALT switch
    LO_OK_REF_SAME_NDNI_REV,        // Same, non-Ins non-Del left-aligned Indel, reverse complemented
    LO_OK_REF_SAME_DEL_REV,         // Same deletion - reverse complemented
    LO_OK_REF_SAME_INS_REV,         // Same insertion - reverse complemented
#define LO_OK_REF_SAME_INDEL_LAST LO_OK_REF_SAME_INS_REV
    LO_OK_REF_SAME_NLA,             // REF and ALT are the same for a non-left-aligned, not SNP variant
    LO_OK_REF_SAME_SV,              // REF and ALT are the same for a variant with a symbolic ALT allele
    LO_OK_REF_ALT_SWITCH_SNP,       // REF and ALT are switched, and INFO and FORMAT subfields updated
    LO_OK_REF_ALT_SWITCH_INDEL,     // REF and ALT are switched, and INFO and FORMAT subfields updated
    LO_OK_REF_ALT_SWITCH_DEL_TO_INS,// Deletion in Primary became incorporated in the Luft reference
    LO_OK_REF_ALT_SWITCH_INDEL_RPTS,// Switched number of payload repeats in reference
    LO_OK_REF_ALT_SWITCH_INDEL_FLANKING,// REF bases the same, but switch called based on flanking regions
    LO_OK_REF_ALT_SWITCH_NDNI,      // Left-anchored non=Ins non-Del INDEL
    LO_OK_REF_ALT_SWITCH_INDEL_WITH_GAP,// Deletion with payload in chain file gap
#define LO_OK_REF_ALT_SWITCH_INDEL_LAST LO_OK_REF_ALT_SWITCH_INDEL_WITH_GAP
    LO_OK_REF_ALT_SWITCH_NLA,       // REF and ALT are switched, and INFO and FORMAT subfields updated
    LO_OK_REF_NEW_SNP,              // REF in a SNP was replaced with a new REF - can only happen if AF=1

    // Rejection reasons - MUST BE AFTER LO_REJECTED (see LO_IS_REJECTED) 
    LO_REJECTED,                    // generic rejection status

    // Reasons due to data (any tool should reject)
    LO_CHROM_NOT_IN_PRIM_REF,       // Primary reference doesn't contain CHROM
    LO_CHROM_NOT_IN_CHAIN,          // chain file doesn't contain a qName which has CHROM  
#define LO_NO_MAPPING_IN_CHAIN_FIRST LO_NO_MAPPING_IN_CHAIN_REF
    LO_NO_MAPPING_IN_CHAIN_REF,     // chain file doesn't contain a mapping covering REF
    LO_NO_MAPPING_IN_CHAIN_ANCHOR,  // New left-anchor base (after reverse-complementing) is before beginning of the chromosome
    LO_NO_MAPPING_REF_SPLIT,        // REF is not fully within a single alignment in the chain file
#define LO_NO_MAPPING_IN_CHAIN_LAST LO_NO_MAPPING_REF_SPLIT
    LO_REF_MISMATCHES_REFERENCE,    // REF different than reference 

    // Reasons due to Genozip limitations (more sophisticated tools might lift)
    LO_REF_MULTIALT_SWITCH_SNP,     // REF changes for a multi-allelic SNP, or REF change would make a bi-allelic into a tri-allelic SNP
    LO_NEW_ALLELE_SNP,              // The Luft reference represents an allele that is neither REF or ALT
    LO_REF_MULTIALT_SWITCH_INDEL,   // REF changes for a multi-allelic INDEL, or REF change would make a bi-allelic into a tri-allelic INDEL
    LO_NEW_ALLELE_DEL_REF_CHANGED_MISSING,  // REF changed in a Deletion variant that has a "*" ALT
    LO_NEW_ALLELE_DEL_REF_CHANGED,  // REF changed in Deletion variant, but not REF⇆ALT switch (i.e. Deletion not integrated into new reference)
    LO_NEW_ALLELE_DEL_SAME_REF,     // REF bases match, but this is a new Deletion allele based on context
    LO_NEW_ALLELE_INS_REF_CHANGED,  // REF changed in Insertion variant
    LO_NEW_ALLELE_INS_SAME_REF,     // REF bases match, but this is a new Insertion allele based on context
    LO_NEW_ALLELE_INDEL_NO_SWITCH,  // REF switched with one of the ALTs, but flanking regions differ
    LO_NEW_ALLELE_NDNI,             // REF is a new allele a left-anchored non-Del non-Ins indel

    LO_XSTRAND_NLA,                 // Used in v12.0.0-12.0.37: Genozip limiation: A complex indel variant mapped to the reverse strand
    LO_NEW_ALLELE_NLA,              // The Luft reference represents an allele that is neither REF or ALT
    LO_NEW_ALLELE_SV,               // Genozip limiation: The Luft reference is different than REF for a variant with a symbolic ALT allele
    LO_XSTRAND_SV,                  // Genozip limiation: A variant with a symbolic ALT allele mapped to the reverse strand
    LO_COMPLEX_REARRANGEMENTS,      // Genozip limiation: Variant contains "complex rearrangements" (see VCF spec)

    // Reasons when genozipping a DVCF file (without --chain)
    LO_ADDED_VARIANT,               // Variant was added to file after it was already lifted over
    LO_UNSUPPORTED_REFALT,          // INFO/DVCF fields indicate an unsupported REF/ALT configuration - perhaps other tool, perhaps later version of Genozip

    LO_INFO,                        // An error cross-rending an INFO subfield (including INFO/END)
    LO_FORMAT,                      // An error cross-rending an FORMAT subfield
    
    LO_ALTS_NOT_SAME_LEN_INS_REV,   // A multi-allelic insertion with xstrand contains ALTs of different length, eg: "C CT,T". Re-anchoring whould make each ALT have a different POS        
    NUM_LO_STATUSES
} LiftOverStatus;

// NOTE: The rejection strings should not change or vcf_lo_seg_INFO_REJX won't identify oStatus in rejected lines when parsing old dual coordinate files (it will fallback to LO_REJECTED)
// It *IS OK* to add more statuses, change their order or change their numeric values, it is *NOT OK* to modify the names as they are part of the file format
extern rom dvcf_status_names[NUM_LO_STATUSES];
#define DVCF_STATUS_NAMES { /* for display esthetics - max 25 characters - note: these strings are defined in the dual-coordinate specification */\
    "UNKNOWN", \
    "OK", "OkRefSameSNP", "OkRefSameSNPRev", "OkRefSameSNPIupac", \
    "OkRefSameIndel", "OkRefSameNDNIRev", "OkRefSameDelRev", "OkRefSameInsRev", \
    "OkRefSameNotLeftAnc", "OkRefSameStructVariant", "OkRefAltSwitchSNP", \
    "OkRefAltSwitchIndel", "OkRefAltSwitchDelToIns", "OkRefAltSwitchIndelRpts", "OkRefAltSwitchIndelFlank", "OkRefAltSwitchNDNI", "OkRefAltSwitchWithGap", \
    "OkRefAltSwitchNotLeftAnc", "OkNewRefSNP", \
    "Rejected", "ChromNotInPrimReference", "ChromNotInChainFile", \
    "RefNotMappedInChain", "NewAnchorNotInChrom", "RefSplitInChain", \
    "RefMismatchesReference", "RefMultiAltSwitchSNP", "RefNewAlleleSNP", "RefMultiAltSwitchIndel", \
    "RefNewAlleleDelRefChgHas*", "RefNewAlleleDelRefChanged", "RefNewAlleleDelSameRef", \
    "RefNewAlleleInsRefChanged", "RefNewAlleleInsSameRef", "RefNewAlleleIndelNoSwitch", "RefNewAlleleNDNI",\
    "XstrandNotLeftAnc", "RefNewAlleleNotLeftAnc","RefNewAlleleSV", "XstrandSV", "ComplexRearrangements", \
    "AddedVariant", "UnsupportedRefAlt", \
    "INFO", "FORMAT", "AltsNotSameLenInsRev" \
}

#define LO_IS_REJECTED(ost) ((ost) >= LO_REJECTED) // note: this condition works also for unrecognized reject strings (that an have index >= NUM_LO_STATUSES)
#define LO_IS_NO_MAPPING(ost) ((ost)>=LO_NO_MAPPING_IN_CHAIN_FIRST && (ost)<=LO_NO_MAPPING_IN_CHAIN_LAST)
#define LO_IS_OK(ost) ((ost) >= LO_OK && (ost) < LO_REJECTED)
#define LO_IS_OK_SAME_INDEL(ost)           ((ost)>=LO_OK_REF_SAME_INDEL       && (ost)<=LO_OK_REF_SAME_INDEL_LAST)
#define LO_IS_OK_REF_ALT_SWITCH_INDEL(ost) ((ost)>=LO_OK_REF_ALT_SWITCH_INDEL && (ost)<=LO_OK_REF_ALT_SWITCH_INDEL_LAST)
#define LO_IS_OK_SNP(ost) ((ost)==LO_OK_REF_SAME_SNP || (ost)==LO_OK_REF_ALT_SWITCH_SNP || (ost)==LO_OK_REF_NEW_SNP)
#define LO_IS_OK_INDEL(ost) (LO_IS_OK_SAME_INDEL(ost) || LO_IS_OK_REF_ALT_SWITCH_INDEL(ost))
#define LO_IS_OK_COMPLEX(ost) ((ost)==LO_OK_REF_SAME_NLA || (ost)==LO_OK_REF_ALT_SWITCH_NLA)
#define LO_IS_OK_SWITCH(ost) ((ost)==LO_OK_REF_ALT_SWITCH_SNP || LO_IS_OK_REF_ALT_SWITCH_INDEL(ost) || (ost)==LO_OK_REF_ALT_SWITCH_NLA)
#define LO_IS_OK_SAME(ost) ((ost)==LO_OK_REF_SAME_SNP || LO_IS_OK_SAME_INDEL(ost) || (ost)==LO_OK_REF_SAME_NLA)
#define last_ostatus (ctx_has_value_in_line_(vb, CTX(VCF_oSTATUS)) ? (LiftOverStatus)vb->last_index (VCF_oSTATUS) : LO_UNKNOWN)
#define last_ostatus_name_piz (last_ostatus < CTX(VCF_oSTATUS)->word_list.len ? ctx_get_words_snip (CTX(VCF_oSTATUS), last_ostatus) : "Invalid oStatus")

// fields and their order of INFO/LIFTOVER and INFO/LIFTBACK containers (part of the file format, see vcf_lo_zip_initialize)
typedef enum { IL_CHROM, IL_POS, IL_REF, IL_XSTRAND, NUM_IL_FIELDS } InfoLiftFields; 

// Liftover - header keys
#define HK_GENOZIP_CMD  "##genozip_command="
#define HK_CHAIN        "##chain="
#define HK_DC           "##dual_coordinates"
#define HK_ORIGINAL_REF "##original_reference="
#define HK_DC_PRIMARY  HK_DC"=PRIMARY"
#define HK_DC_LUFT     HK_DC"=LUFT"
#define HK_RENDALG            "RendAlg"
#define HK_RENDALG_ATTR       HK_RENDALG"="
#define HK_RENAME_REFALT_ATTR "RenameRefalt="
#define HK_RENAME_STRAND_ATTR "RenameStrand="
#define HK_RENAME_TLAFER_ATTR "RenameTlafer="
#define HK_RENAME_ALWAYS_ATTR "RenameAlways="
#define TAG_SOURCE      "Source=\""GENOZIP_URL"\""
#define KH_INFO         "##INFO=<ID="
#define KH_INFO_LUFT    KH_INFO INFO_LUFT_NAME ",Number=4,Type=String,Description=\"Info for rendering variant in LUFT coords. See " WEBSITE_DVCF "\","TAG_SOURCE",Version=\"%s\"," HK_RENDALG_ATTR "\"NONE\">"
#define KH_INFO_PRIM    KH_INFO INFO_PRIM_NAME ",Number=4,Type=String,Description=\"Info for rendering variant in PRIMARY coords\","TAG_SOURCE",Version=\"%s\"," HK_RENDALG_ATTR "\"NONE\">"
#define KH_INFO_LREJ    KH_INFO INFO_LREJ_NAME ",Number=1,Type=String,Description=\"Reason variant was rejected for LUFT coords\","TAG_SOURCE",Version=\"%s\"," HK_RENDALG_ATTR "\"NONE\">"
#define KH_INFO_PREJ    KH_INFO INFO_PREJ_NAME ",Number=1,Type=String,Description=\"Reason variant was rejected for PRIMARY coords\","TAG_SOURCE",Version=\"%s\"," HK_RENDALG_ATTR "\"NONE\">"
#define KH_INFO_oSTATUS KH_INFO "oSTATUS,Number=1,Type=String,Description=\"Lift status\","TAG_SOURCE",Version=\"%s\"," HK_RENDALG_ATTR "\"NONE\">"

// VCF standard keys
#define HK_CONTIG      "##contig="

// header keys that appear only in a Primary VCF file
#define HK_LUFT_CONTIG "##luft_contig="     
#define HK_LUFT_REF    "##luft_reference="
#define HK_LUFT_ONLY   "##luft_only="

// header keys that appear only in Luft VCF file
#define HK_PRIM_CONTIG "##primary_contig="
#define HK_PRIM_REF    "##primary_reference="
#define HK_PRIM_ONLY   "##primary_only="

// Header stuff
typedef enum { RA_REFALT=0, RA_STRAND, RA_TLAFER, RA_ALWAYS, NUM_RENAME_ATTRS } RenameAttr;
extern rom vcf_header_rename_attrs[NUM_RENAME_ATTRS];
extern const unsigned vcf_header_rename_attr_lens[NUM_RENAME_ATTRS];

extern uint32_t vcf_num_samples; // ZIP
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag.samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples
extern VcfVersion vcf_header_get_version (void);

// Samples stuff
extern void vcf_seg_FORMAT (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(fmt));
extern void vcf_samples_zip_initialize (void);
extern void vcf_samples_seg_initialize (VBlockVCFP vb);
extern void vcf_samples_seg_finalize (VBlockVCFP vb);

extern rom vcf_seg_samples (VBlockVCFP vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13);
extern int vcf_seg_get_mux_channel_i (VBlockVCFP vb, bool fail_if_dvcf_refalt_switch);
extern int vcf_piz_get_mux_channel_i (VBlockP vb);
extern ContextP vcf_seg_FORMAT_mux_by_dosage (VBlockVCFP vb, ContextP ctx, STRp(cell), const DosageMultiplexer *mux);

eSTRl(af_snip);

// FORMAT/GT stuff
extern WordIndex vcf_seg_FORMAT_GT (VBlockVCFP vb, ContextP ctx, ZipDataLineVCF *dl, STRp(cell), bool has_ps, bool has_null_dp);
extern void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCFP vb);
extern void vcf_piz_FORMAT_GT_rewrite_predicted_phase (VBlockP vb, char *recon, uint32_t recon_len);
extern void vcf_piz_GT_cb_null_GT_if_null_DP (VBlockP vb , char *recon);
extern int vcf_piz_GT_get_last_dosage (VBlockP vb);

// GIAB trio stuff
extern void vcf_giab_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_FORMAT_IGT (VBlockVCFP vb, ContextP ctx, STRp(igt));
extern void vcf_seg_FORMAT_IPS (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP ctx, STRp(ips));
extern void vcf_seg_ADALL_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values);

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

// INFO/SF
extern bool vcf_seg_INFO_SF_init (VBlockVCFP vb, ContextP sf_ctx, STRp(value));
extern void vcf_seg_INFO_SF_seg (VBlockVCFP vb);
extern void vcf_seg_INFO_SF_one_sample (VBlockVCFP vb);
extern void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *recon, int32_t recon_len);
extern int vcf_piz_TOPLEVEL_cb_insert_INFO_SF (VBlockVCFP vcf_vb);

// INFO/QD stuff
typedef enum { QD_PRED_NONE, QD_PRED_INFO_DP, QD_PRED_INFO_DP_P001, QD_PRED_INFO_DP_M001, 
                             QD_PRED_SUM_DP,  QD_PRED_SUM_DP_P001,  QD_PRED_SUM_DP_M001, NUM_QD_PRED_TYPES } QdPredType;
extern void vcf_seg_sum_DP_for_QD (VBlockVCFP vb, int64_t value);
extern void vcf_seg_INFO_QD (VBlockVCFP vb);
extern void vcf_piz_sum_DP_for_QD (VBlockP vb, STRp(recon));
extern void vcf_piz_insert_QD (VBlockVCFP vb);

// INFO stuff

typedef struct { char name[MAX_TAG_LEN]; // not nul-terminated, including '=' if there is one
                 rom value; 
                 unsigned name_len, value_len; 
                 ContextP ctx; } InfoItem;

extern void vcf_info_zip_initialize (void);
extern void vcf_info_seg_initialize (VBlockVCFP vb);
extern void vcf_piz_finalize_DP_by_DP (VBlockVCFP vb);

extern void vcf_seg_info_subfields (VBlockVCFP vb, STRp(info));
extern void vcf_finalize_seg_info (VBlockVCFP vb);
extern bool vcf_seg_INFO_allele (VBlockP vb_, ContextP ctx, STRp(value), uint32_t repeat);

// Refalt stuff
extern void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, STRp(ref), STRp(alt));
extern void vcf_refalt_seg_other_REFALT (VBlockVCFP vb, Did did_i, LiftOverStatus ostatus, bool is_xstrand, unsigned add_bytes);
extern LiftOverStatus vcf_refalt_lift (VBlockVCFP vb, const ZipDataLineVCF *dl, bool xstrand, WordIndex luft_ref_index, bool *is_left_anchored);
typedef enum { EQUALS_NEITHER, EQUALS_REF, EQUALS_ALT, EQUALS_MISSING } RefAltEquals;
RefAltEquals vcf_refalt_oref_equals_ref_or_alt (char oref, char ref, STRp(alt), bool is_xstrand);
extern bool vcf_refalt_piz_is_variant_snp (VBlockVCFP vb);
extern bool vcf_refalt_piz_is_variant_indel (VBlockVCFP vb);
extern void vcf_refalt_seg_convert_to_primary (VBlockVCFP vb, LiftOverStatus ostatus);
extern void vcf_piz_refalt_parse (VBlockVCFP vb, STRp(refalt));

// GVCF stuff
extern bool vcf_piz_line_has_RGQ (VBlockVCFP vb);

// Illumina genotyping stuff
extern void vcf_illum_gtyping_initialize (VBlockVCFP vb);
extern void vcf_seg_PROBE_A (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_PROBE_B (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_ILLUMINA_CHR (VBlockVCFP vb, ContextP ctx, STRp(chr));
extern void vcf_seg_ILLUMINA_POS (VBlockVCFP vb, ContextP ctx, STRp(pos));
extern void vcf_seg_ILLUMINA_STRAND (VBlockVCFP vb, ContextP ctx, STRp(strand));
extern void vcf_seg_ALLELE_A (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_ALLELE_B (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_seg_mux_by_adjusted_dosage (VBlockVCFP vb, ContextP ctx, STRp(baf), const DosageMultiplexer *mux);

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
extern void vcf_vep_zip_initialize (rom spec, rom buf_1st);
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

// ICGC stuff
extern void vcf_seg_INFO_mutation (VBlockVCFP vb, ContextP ctx, STRp(mut));

// ISAAC stuff
extern void vcf_isaac_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_FORMAT_GQX (VBlockVCFP vb, ContextP ctx, STRp(gqx));
extern void vcf_seg_INFO_RU (VBlockVCFP vb, ContextP ctx, STRp(ru));
extern void vcf_seg_INFO_IDREP (VBlockVCFP vb, ContextP ctx, STRp(idrep));

// manta stuff
extern void vcf_manta_seg_initialize (VBlockVCFP vb);
extern void vcf_seg_manta_ID (VBlockVCFP vb, STRp(id));

// Tags stuff

// tag sources - ordered from least authorative to most 
typedef enum { TAG_NO_SRC, TAG_GENOZIP, TAG_HEADER_DST, TAG_HEADER, TAG_CMDLINE_DST, TAG_CMDLINE } VcfTagSource;
#define TAG_SOURCE_NAMES { "NoSrc", "Genozip", "HeaderDst", "Header", "CmdLineDst", "CmdLine" }

#define MAX_NUMBER_LEN 8            // maximum length of Number attribute 
#define MAX_TYPE_LEN 12
#define MAX_RENDALG_LEN (MAX_TAG_LEN+16)
typedef struct { 
    char tag_name[MAX_TAG_LEN];     // this can also be refered to as dests[RA_NAME] since MAX_TAG_LEN is word-aligned
    char dests[NUM_RENAME_ATTRS][MAX_TAG_LEN];
    unsigned tag_name_len, dest_lens[NUM_RENAME_ATTRS];
    DictIdType dtype;
    STRl(number, MAX_NUMBER_LEN);   // Number attribute  (only used in --chain)
    STRl(type, MAX_TYPE_LEN);       // Type attribute    (only used in --chain)
    STRl(rendalg, MAX_RENDALG_LEN); // RendAlg attribute (only used in --chain)
    VcfTagSource source;            // where this tag originated
} Tag;

extern void vcf_tags_populate_tags_from_command_line (void);
extern void vcf_tags_add_tag (VBlockVCFP vb, ContextP ctx, DictIdType dtype, STRp(tag_name));
extern unsigned vcf_tags_rename (VBlockVCFP vb, unsigned num_tags, const ContextPBlock ctxs, rom sf_names[], const unsigned sf_name_lens[], const InfoItem *info_items, char *renamed);
extern void vcf_tags_finalize_tags_from_vcf_header (void);
extern bool vcf_tags_add_attr_from_header (DictIdType dtype, STRp(tag_name), RenameAttr attr, STRp(number), STRp (type), STRp (rendalg), pSTRp(dest), bool recursive);
extern Tag *vcf_tags_get_next_missing_tag (Tag *tag);

// Liftover Zip
extern void vcf_lo_zip_initialize (void);
extern TranslatorId vcf_lo_luft_trans_id (DictId dict_id, char number);
extern void vcf_lo_seg_generate_INFO_DVCF (VBlockVCFP vb, ZipDataLineVCF *dl);
extern void vcf_lo_set_rollback_point (VBlockVCFP vb);
extern void vcf_lo_seg_rollback_and_reject (VBlockVCFP vb, LiftOverStatus ostatus, Context *ctx);
extern LiftOverStatus vcf_lo_get_liftover_coords (VBlockVCFP vb, PosType32 pos, WordIndex *dst_contig_index, PosType32 *dst_1pos, bool *xstrand, uint32_t *aln_i); // out
extern void vcf_lo_seg_INFO_LUFT_and_PRIM (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_lo_seg_INFO_REJX (VBlockVCFP vb, ContextP ctx, STRp(value));
extern bool vcf_lo_seg_cross_render_to_primary (VBlockVCFP vb, ContextP ctx, STRp (this_value), qSTRp (primary_snip), bool validate_only);

// Line sorter
typedef struct {
    WordIndex chrom_wi; 
    uint32_t tie_breaker;
    PosType32 start_pos, end_pos;
} LineCmpInfo; 
extern bool vcf_is_sorting (CompIType comp_i);
extern int vcf_linesort_cmp (LineCmpInfo a, LineCmpInfo b);
extern void vcf_linesort_merge_vb (VBlockP vb);

#define vcf_set_ostatus(ostatus) ctx_set_last_value (VB, CTX(VCF_oSTATUS), (int64_t)(ostatus))

// Liftover Piz
extern void vcf_lo_piz_TOPLEVEL_cb_filter_line (VBlockVCFP vb);

#define VCF_ERR_PREFIX { progress_newline(); fprintf (stderr, "Error %s:%u in variant %s=%.*s %s=%"PRId64": ", __FUNCLINE, (VB_VCF->line_coords == DC_PRIMARY ? "CHROM" : "oCHROM"), vb->chrom_name_len, vb->chrom_name, (VB_VCF->line_coords == DC_PRIMARY ? "POS" : "oPOS"), vb->last_int (VB_VCF->line_coords == DC_PRIMARY ? VCF_POS : VCF_oPOS)); }
#define ASSVCF(condition, format, ...) ({ if (!(condition)) { VCF_ERR_PREFIX; fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }})
#define ASSVCF0(condition, msg)        ASSVCF ((condition), msg "%s", "")
#define WARNVCF(format, ...)           ({ if (!flag.quiet)  { VCF_ERR_PREFIX; fprintf (stderr, format "\n", __VA_ARGS__); } })

#define REJECT(ostatus, reason, ...)              ({ if (!vb->comp_i) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%d\t%.*s%s\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN_(100, vb->main_ref_len), vb->main_ref, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); })
#define REJECT_MAPPING(reason)                    ({ if (!vb->comp_i) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%d\t.\t"      reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0]); return; })
#define REJECT_SUBFIELD(ostatus, ctx, reason,...) ({ if (!vb->comp_i) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%d\t.\t"      reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], __VA_ARGS__); \
                                                       vcf_lo_seg_rollback_and_reject (vb, (ostatus), (ctx)); })
// misc
extern void vcf_piz_insert_field (VBlockVCFP vb, Did did, STRp(value));
