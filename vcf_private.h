// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "vblock.h"
#include "vcf.h"
#include "website.h"
#include "seg.h"
#include "container.h"

#define VCF_FIELD_NAMES "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
#define VCF_FIELD_NAMES_LONG VCF_FIELD_NAMES "\tFORMAT"

// IMPORTANT: if changing fields in DataLine, also update vb_release_vb
typedef struct {
    VBLOCK_COMMON_LINES_ZIP
    
    bool has_haplotype_data : 1; // FORMAT field contains GT
    bool has_genotype_data  : 1; // FORMAT field contains subfields other than GT

    WordIndex format_node_i; // the node_index into contexts[VCF_FORMAT].nodes and also format_mapper_buf that applies to this line. Data on the fields is in vb->format_mapper_buf[dl.format_node_i]
} ZipDataLineVCF;
#define DATA_LINE(i) ENT (ZipDataLineVCF, vb->lines, i)

typedef enum { VCF_v_UNKNOWN, VCF_v4_1, VCF_v4_2, VCF_v4_3, VCF_v4_4, VCF_v4_5 } VcfVersion;

typedef MULTIPLEXER(4) DosageMultiplexer;

#define VCF_MAX_ARRAY_ITEMS SMALL_CON_NITEMS

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
typedef struct VBlockVCF {

    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    
    uint16_t ploidy;                // ZIP only
    VcfVersion vcf_version;

    // used for segging FORMAT/GT
    uint32_t gt_prev_ploidy;
    char gt_prev_phase;
    
    // used for segging INFO
    Buffer info_items;              // Seg: INFO items of the line being segged

    const char *main_refalt;        // used by vcf_refalt_lift and vcf_seg_INFO_BaseCounts, set by vcf_seg_txt_line
    unsigned main_ref_len, main_alt_len;

    // INFO/SF stuff
    enum { USE_SF_UNKNOWN, USE_SF_YES, USE_SF_NO } use_special_sf;
    Buffer sf_txt, sf_snip; // INFO/SF data as it appears in the snip being constructed

    // INFO/END
    uint64_t last_end_line_i;       // PIZ: last line on which INFO/END was encountered
    
    // FORMAT/DP
    int64_t sum_dp_this_line;       // ZIP: sum of values of FORMAT/DP in samples of this line - '.' counts as 0.
    uint32_t num_dps_this_line;     // ZIP: possibley less that vcf_num_samples, in case of missing samples or missing DP fields in some of the samples
    
    // FORMAT/AD
    int64_t ad_values[VCF_MAX_ARRAY_ITEMS];

    // FORMAT/PS
    bool PS_encountered_last_line;
    
    // FORMAT stuff 
    Buffer format_mapper_buf;       // ZIP only: an array of type Container - one entry per entry in CTX(VCF_FORMAT)->nodes   
    Buffer format_contexts;         // ZIP only: an array of format_mapper_buf.len of ContextBlock
    Buffer last_format;             // ZIP only: cache previous line's FORMAT string

    // Multiplexers
    DosageMultiplexer mux_PLn, mux_GL, mux_GP, mux_PRI, mux_DS, mux_PP, mux_PVAL, mux_FREQ, mux_RD,
                      mux_AD[2], mux_ADALL[2];

    PLMuxByDP PL_mux_by_DP;
    
    MULTIPLEXER(1 + 50 * 3) mux_PLy; // num_DP x 3 + 1
    MULTIPLEXER(1 + 7 * 3) mux_GQ;

    // used by CODEC_HAPM (for VCF haplotype matrix) 
    Buffer hapmat_helper_index_buf; // ZIP: used by codec_hapmat_count_alt_alleles 
    Buffer hapmat_columns_data;     // used by codec_hapmat_piz_get_one_line 
    Buffer hapmat_column_of_zeros;  // used by codec_hapmat_piz_calculate_columns   
    Buffer hapmat_one_array;        // one line or column 

    // DVCF stuff
    Buffer tags;                    // Seg: used for FORMAT and INFO tag renaming.
    Buffer rejects_report;          // human readable report about rejects
    char new_ref;                   // SNP: new REF that is neither REF nor ALT; left-anchored INDEL with XSTRAND: new anchor
    bool is_del_sv;                 // is ALT == "<DEL>"
    bool prev_line_rejected;        // ZIP only
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
    LO_OK_REF_SAME_INDEL,           // REF and ALT are the same for an INDEL, and confirmed not to be a REF<>ALT switch
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
    LO_NEW_ALLELE_DEL_REF_CHANGED,  // REF changed in Deletion variant, but not REF<>ALT switch (i.e. Deletion not integrated into new reference)
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
    
    NUM_LO_STATUSES
} LiftOverStatus;

// NOTE: The rejection strings should not change or vcf_lo_seg_INFO_REJX won't identify oStatus in rejected lines when parsing old dual coordinate files (it will fallback to LO_REJECTED)
// It *IS OK* to add more statuses, change their order or change their numeric values, it is *NOT OK* to modify the names as they are part of the file format
extern const char *dvcf_status_names[NUM_LO_STATUSES];
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
    "RefNewAllelInsRefChanged", "RefNewAlleleInsSameRef", "RefNewAlleleIndelNoSwitch", "RefNewAlleleNDNI",\
    "XstrandNotLeftAnc", "RefNewAlleleNotLeftAnc","RefNewAlleleSV", "XstrandSV", "ComplexRearrangements", \
    "AddedVariant", "UnsupportedRefAlt", \
    "INFO", "FORMAT" \
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
#define KH_INFO_LUFT    KH_INFO INFO_LUFT_NAME ",Number=4,Type=String,Description=\"Info for rendering variant in LUFT coords. See " WEBSITE_COORDS "\","TAG_SOURCE",Version=\"%s\"," HK_RENDALG_ATTR "\"NONE\">"
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
extern const char *vcf_header_rename_attrs[NUM_RENAME_ATTRS];
extern const unsigned vcf_header_rename_attr_lens[NUM_RENAME_ATTRS];

extern uint32_t vcf_num_samples; // ZIP
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag.samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples
extern VcfVersion vcf_header_get_version (void);

// Samples stuff
extern void vcf_samples_zip_initialize (void);
extern void vcf_samples_seg_initialize (VBlockVCFP vb);
extern void vcf_samples_seg_finalize (VBlockVCFP vb);

extern const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13, const char *backup_luft_samples, uint32_t backup_luft_samples_len);
extern void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCF *vb);
extern void vcf_piz_FORMAT_GT_rewrite_predicted_phase (VBlockP vb, char *recon, uint32_t recon_len);

// FORMAT/GT stuff
extern WordIndex vcf_seg_FORMAT_GT (VBlockVCFP vb, ContextP ctx, ZipDataLineVCF *dl, STRp(cell), bool has_ps, bool has_null_dp);
extern void vcf_piz_GT_cb_null_GT_if_null_DP (VBlockP vb , char *recon);

#define IS_TRIVAL_FORMAT_SUBFIELD ((!recon_len || (recon_len==1 && *recon=='.')) && dict_id_is_vcf_format_sf (ctx->dict_id))
extern void vcf_FORMAT_PL_decide (VBlockVCFP vb);
extern void vcf_FORMAT_PL_after_vbs (void);

// INFO stuff

typedef struct { char name[MAX_TAG_LEN]; // not nul-terminated, including '=' if there is one
                 const char *value; 
                 unsigned name_len, value_len; 
                 ContextP ctx; } InfoItem;

extern void vcf_info_zip_initialize (void);
extern void vcf_info_seg_initialize (VBlockVCFP vb);
extern void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *recon, int32_t recon_len);
extern void vcf_piz_TOPLEVEL_cb_insert_INFO_SF (VBlockVCFP vcf_vb);
extern void vcf_seg_INFO_SF_one_sample (VBlockVCF *vb);
extern void vcf_seg_info_subfields (VBlockVCF *vb, const char *info_str, unsigned info_len);
extern void vcf_finalize_seg_info (VBlockVCF *vb);

// Refalt stuff
extern void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, STRp(ref), STRp(alt));
extern void vcf_refalt_seg_other_REFALT (VBlockVCFP vb, DidIType did_i, LiftOverStatus ostatus, bool is_xstrand, unsigned add_bytes);
extern LiftOverStatus vcf_refalt_lift (VBlockVCFP vb, const ZipDataLineVCF *dl, bool xstrand, WordIndex luft_ref_index, bool *is_left_anchored);
typedef enum { EQUALS_NEITHER, EQUALS_REF, EQUALS_ALT, EQUALS_MISSING } RefAltEquals;
RefAltEquals vcf_refalt_oref_equals_ref_or_alt (char oref, char ref, STRp(alt), bool is_xstrand);
extern bool vcf_refalt_piz_is_variant_snp (VBlockP vb);
extern bool vcf_refalt_piz_is_variant_indel (VBlockP vb);
extern void vcf_refalt_seg_convert_to_primary (VBlockVCFP vb, LiftOverStatus ostatus);

// Tags stuff

// tag sources - ordered from least authorative to most 
typedef enum { TAG_NO_SRC, TAG_GENOZIP, TAG_HEADER_DST, TAG_HEADER, TAG_CMDLINE_DST, TAG_CMDLINE } VcfTagSource;
#define TAG_SOURCE_NAMES { "NoSrc", "Genozip", "HeaderDst", "Header", "CmdLineDst", "CmdLine" }

#define MAX_NUMBER_LEN 8 // maximum length of Number attribute 
#define MAX_TYPE_LEN 12
#define MAX_RENDALG_LEN (MAX_TAG_LEN+16)
typedef struct { 
    char tag_name[MAX_TAG_LEN]; // this can also be refered to as dests[RA_NAME] since MAX_TAG_LEN is word-aligned
    char dests[NUM_RENAME_ATTRS][MAX_TAG_LEN];
    unsigned tag_name_len, dest_lens[NUM_RENAME_ATTRS];
    DictIdType dtype;
    char number[MAX_NUMBER_LEN];   unsigned number_len;  // Number attribute  (only used in --chain)
    char type[MAX_TYPE_LEN];       unsigned type_len;    // Type attribute    (only used in --chain)
    char rendalg[MAX_RENDALG_LEN]; unsigned rendalg_len; // RendAlg attribute (only used in --chain)
    VcfTagSource source; // where did this tag originate
} Tag;

extern void vcf_tags_populate_tags_from_command_line (void);
extern void vcf_tags_add_tag (VBlockVCFP vb, ContextP ctx, DictIdType dtype, STRp(tag_name));
extern unsigned vcf_tags_rename (VBlockVCFP vb, unsigned num_tags, const ContextPBlock ctxs, const char *sf_names[], const unsigned sf_name_lens[], const InfoItem *info_items, char *renamed);
extern void vcf_tags_finalize_tags_from_vcf_header (void);
extern bool vcf_tags_add_attr_from_header (DictIdType dtype, STRp(tag_name), RenameAttr attr, STRp(number), STRp (type), STRp (rendalg), pSTRp(dest), bool recursive);
extern Tag *vcf_tags_get_next_missing_tag (Tag *tag);

// Liftover Zip
extern void vcf_lo_zip_initialize (void);
extern void vcf_lo_append_rejects_file (VBlockP vb, Coords coord);
extern void vcf_lo_seg_generate_INFO_DVCF (VBlockVCFP vb, ZipDataLineVCF *dl);
extern void vcf_lo_set_rollback_point (VBlockVCFP vb);
extern void vcf_lo_seg_rollback_and_reject (VBlockVCFP vb, LiftOverStatus ostatus, Context *ctx);
extern LiftOverStatus vcf_lo_get_liftover_coords (VBlockVCFP vb, PosType pos, WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i); // out
extern void vcf_lo_seg_INFO_LUFT_and_PRIM (VBlockVCFP vb, ContextP ctx, STRp(value));
extern void vcf_lo_seg_INFO_REJX (VBlockVCFP vb, ContextP ctx, STRp(value));
extern bool vcf_lo_seg_cross_render_to_primary (VBlockVCFP vb, ContextP ctx, STRp (this_value), char *modified_snip, unsigned *modified_snip_len);

#define vcf_set_ostatus(ostatus) ctx_set_last_value (VB, CTX(VCF_oSTATUS), (int64_t)(ostatus))

// Liftover Piz
extern void vcf_lo_piz_TOPLEVEL_cb_filter_line (VBlockP vb);

#define VCF_ERR_PREFIX { progress_newline(); fprintf (stderr, "Error %s:%u in variant %s=%.*s %s=%"PRId64": ", __FUNCTION__, __LINE__, (vb->line_coords == DC_PRIMARY ? "CHROM" : "oCHROM"), vb->chrom_name_len, vb->chrom_name, (vb->line_coords == DC_PRIMARY ? "POS" : "oPOS"), vb->last_int (vb->line_coords == DC_PRIMARY ? VCF_POS : VCF_oPOS)); }
#define ASSVCF(condition, format, ...) do { if (!(condition)) { VCF_ERR_PREFIX; fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSVCF0(condition, msg)        ASSVCF ((condition), msg "%s", "")
#define WARNVCF(format, ...)           do { if (!flag.quiet)  { VCF_ERR_PREFIX; fprintf (stderr, format "\n", __VA_ARGS__); } } while(0)

#define REJECT(ostatus, reason, ...)              do { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%"PRId64"\t%.*s%s\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN_(100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } while(0)
#define REJECT_MAPPING(reason)                    do { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%"PRId64 "\t.\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0]); return; } while(0)
#define REJECT_SUBFIELD(ostatus, ctx, reason,...) do { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%"PRId64"\t.\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], __VA_ARGS__); \
                                                       vcf_lo_seg_rollback_and_reject (vb, (ostatus), (ctx)); } while(0)
#define LIFTOK(ostatus, reason, ...) do { if (flag.show_lift && !flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%"PRId64"\t%.*s%s\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN_(100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } while(0)
#define LIFTOKEXT(ostatus, reason, ...) do { if (flag.show_lift && !flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%"PRId64"\t%.*s%s\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN_(100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } while(0)
#define LIFTOK0(ostatus, reason) LIFTOK(ostatus, reason "%s", "")
#define REJECTIF(condition, ostatus, reason, ...) do { if (condition) { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s\t%.*s\t%"PRId64"\t%.*s%s\t" reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN_(100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } } while(0)
#define REJECTIF0(condition, ostatus, reason)      do { if (condition) REJECT (ostatus, reason "%s", ""); } while (0)
