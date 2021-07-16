// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef VCF_PRIVATE_INCLUDED
#define VCF_PRIVATE_INCLUDED

#include "vblock.h"
#include "vcf.h"
#include "website.h"

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

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
typedef struct VBlockVCF {

    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    
    uint16_t ploidy;           // ZIP only
    VcfVersion vcf_version;

    uint32_t sample_i;         // ZIP: current sample in line (0-based) being segmented 

    // used for segging FORMAT/GT
    Context *gt_ctx;
    uint32_t gt_prev_ploidy;
    char gt_prev_phase;
    
    // used for segging INFO
    Buffer info_items;              // Seg: INFO items of the line being segged

    const char *main_refalt;        // used by vcf_refalt_lift and vcf_seg_INFO_BaseCounts, set by vcf_seg_txt_line
    unsigned main_ref_len, main_alt_len;

    // INFO/SF stuff
    Context *sf_ctx;
    enum { USE_SF_UNKNOWN, USE_SF_YES, USE_SF_NO } use_special_sf;
    Buffer sf_txt, sf_snip; // INFO/SF data as it appears in the snip being constructed

    // INFO/END
    uint32_t last_end_line_i;       // PIZ: last line on which INFO/END was encountered
    
    // FORMAT/AD
    #define MAX_ARG_ARRAY_ITEMS 36  // same as MAX_COMPOUND_COMPONENTS
    int64_t ad_values[MAX_ARG_ARRAY_ITEMS];
    Context *adall_ctx;             // save FORMAT/ADALL context, to avoid searching for ADALL if it does not exist, since it will conflict in the map with AD

    // dictionaries stuff 
    Buffer format_mapper_buf;       // ZIP only: an array of type Container - one entry per entry in CTX(VCF_FORMAT)->nodes   
    Buffer format_contexts;         // ZIP only: an array of MAX_FIELDS * format_mapper_buf.len of ContextP

    // used by CODEC_HAPM (for VCF haplotype matrix) 
    Context *hapmat_index_ctx; 
    Buffer hapmat_helper_index_buf; // ZIP: used by codec_hapmat_count_alt_alleles 
    Buffer hapmat_columns_data;     // used by codec_hapmat_piz_get_one_line 
    Buffer hapmat_column_of_zeros;  // used by codec_hapmat_piz_calculate_columns   
    Buffer hapmat_one_array;        // one line or column 

    // DVCF stuff
    Buffer rejects_report;          // human readable report about rejects
    char new_ref;                   // new REF that is neither REF nor ALT
    bool is_del_sv;                 // is ALT == "<DEL>"
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;

// Liftover stuff

// this is part of the file format (in the oSTATUS context) and CANNOT BE CHANGED (but can be extended)
typedef enum {
    // Algorithm interim statuses not outputed to the genozip file
    LO_NA                        = -1, // not a Liftover item - never written z_file
    LO_UNKNOWN                   = 0,  // Status of this Liftover line not known yet - never written z_file

    // OK statuses - this are internal to genozip and z_file.oStatus - not expressed in the txt file
    LO_OK                        = 1,  // internal progress step in ZIP - never written z_file
    LO_OK_REF_SAME_SNP           = 2,  // REF and ALT are the same for a SNP
    LO_OK_REF_SAME_INDEL         = 3,  // REF and ALT are the same for an INDEL, and confirmed not to be a REF<>ALT switch
    LO_OK_REF_SAME_SV            = 4,  // REF and ALT are the same for a variant with a symbolic ALT allele
    LO_OK_REF_ALT_SWITCH_SNP     = 5,  // REF and ALT are switched, and INFO and FORMAT subfields updated
    LO_OK_REF_ALT_SWITCH_INDEL   = 6,  // REF and ALT are switched, and INFO and FORMAT subfields updated
    LO_OK_REF_NEW_SNP            = 7,  // REF in a SNP was replaced with a new REF - can only happen if AF=1

    // Rejection reasons - MUST BE AFTER LO_REJECTED (see LO_IS_REJECTED) 
    LO_REJECTED                  = 8,  // generic rejection status

    // Reasons due to data (any tool should reject)
    LO_CHROM_NOT_IN_PRIM_REF     = 9,  // Primary reference doesn't contain CHROM
    LO_CHROM_NOT_IN_CHAIN        = 10, // chain file doesn't contain a qName which has CHROM  
    LO_NO_MAPPING_IN_CHAIN       = 11, // chain file doesn't contain a mapping for the primary POS
    LO_REF_MISMATCHES_REFERENCE  = 12, // REF different than reference 

    // Reasons due to Genozip limitations (more sophisticated tools might lift)
    LO_REF_MULTIALT_SWITCH_SNP   = 13, // REF changes for a multi-allelic SNP, or REF change would make a bi-allelic into a tri-allelic SNP
    LO_NEW_ALLELE_SNP            = 14, // The Luft reference represents an allele that is neither REF or ALT
    LO_REF_MULTIALT_SWITCH_INDEL = 15, // REF changes for a multi-allelic INDEL, or REF change would make a bi-allelic into a tri-allelic INDEL
    LO_NEW_ALLELE_INDEL          = 16, // The Luft reference represents an allele that is neither REF or ALT
    LO_NEW_ALLELE_SV             = 17, // Genozip limiation: The Luft reference is different than REF for a variant with a symbolic ALT allele
    LO_XSTRAND_SV                = 18, // Genozip limiation: A variant with a symbolic ALT allele mapped to the reverse strand
    LO_COMPLEX_REARRANGEMENTS    = 19, // Genozip limiation: Variant contains VCF 4.2 "complex rearrangements"
    LO_NOT_LEFT_ANCHORED         = 20, // Genozip limiation: We require non-SNP variants to be left anchored (i.e. A AG, not G AG)

    // Reasons when genozipping a DVCF file (without --chain)
    LO_ADDED_VARIANT             = 21, // Variant was added to file after it was already lifted over
    LO_UNSUPPORTED_REFALT        = 22, // INFO/DVCF fields indicate an unsupported REF/ALT configuration - perhaps other tool, perhaps later version of Genozip

    LO_INFO                      = 23, // An error cross-rending an INFO subfield (including INFO/END)
    LO_FORMAT                    = 24, // An error cross-rending an FORMAT subfield
    
    NUM_LO_STATUSES
} LiftOverStatus;

// NOTE: The rejection strings should not change or vcf_lo_seg_INFO_REJX won't identify oStatus in rejected lines when parsing old dual coordinate files (it will fallback to LO_REJECTED)
// It *IS OK* to add more statuses, change their order or change their numeric values.
extern const char *dvcf_status_names[];
#define DVCF_STATUS_NAMES { /* for display esthetics - max 25 characters - note: these strings are defined in the dual-coordinate specification */\
    "UNKNOWN", \
    "OK", "OkRefSameSNP", "OkRefSameIndel", "OkRefSameStructVariant", "OkRefAltSwitchSNP", "OkRefAltSwitchIndel", "OkNewRefSNP", \
    "Rejected", "ChromNotInPrimReference", "ChromNotInChainFile", "NoMappingInChainFile", "REFMismatchesReference", \
    "RefMultiAltSwitchSNP", "RefNewAlleleSNP", \
    "RefMultiAltSwitchIndel", "RefNewAlleleIndel", \
    "RefNewAlleleSV", "XstrandSV", "ComplexRearrangements", "NotLeftAnchored", \
    "AddedVariant", "UnsupportedRefAlt", \
    "INFO", "FORMAT" \
}

#define LO_IS_REJECTED(ost) ((ost) >= LO_REJECTED) // note: this condition works also for unrecognized reject strings (that an have index >= NUM_LO_STATUSES)
#define LO_IS_OK(ost) ((ost) >= LO_OK && (ost) < LO_REJECTED)
#define LO_IS_OK_SNP(ost) ((ost)==LO_OK_REF_SAME_SNP || (ost)==LO_OK_REF_ALT_SWITCH_SNP || (ost)==LO_OK_REF_NEW_SNP)
#define LO_IS_OK_INDEL(ost) ((ost)==LO_OK_REF_SAME_INDEL || (ost)==LO_OK_REF_ALT_SWITCH_INDEL)
#define LO_IS_OK_SWITCH(ost) ((ost)==LO_OK_REF_ALT_SWITCH_SNP || (ost)==LO_OK_REF_ALT_SWITCH_INDEL)

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
#define HK_RENDERALG_ATTR "RendAlg"
#define KH_INFO_LUFT    "##INFO=<ID=" INFO_LUFT ",Number=4,Type=String,Description=\"Info for rendering variant in LUFT coords. See " WEBSITE_COORDS "\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_PRIM    "##INFO=<ID=" INFO_PRIM ",Number=4,Type=String,Description=\"Info for rendering variant in PRIMARY coords\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_LREJ    "##INFO=<ID=" INFO_LREJ ",Number=1,Type=String,Description=\"Reason variant was rejected for LUFT coords\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_PREJ    "##INFO=<ID=" INFO_PREJ ",Number=1,Type=String,Description=\"Reason variant was rejected for PRIMARY coords\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_oSTATUS "##INFO=<ID=oSTATUS,Number=1,Type=String,Description=\"Lift status\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"

// header keys that appear only in a Primary VCF file
#define HK_LUFT_CONTIG "##luft_contig="     
#define HK_LUFT_REF    "##luft_reference="
#define HK_LUFT_ONLY   "##luft_only="

// header keys that appear only in Luft VCF file
#define HK_PRIM_CONTIG "##primary_contig="
#define HK_PRIM_REF    "##primary_reference="
#define HK_PRIM_ONLY   "##primary_only="

// Header stuff
extern uint32_t vcf_num_samples; // ZIP
extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag.samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples

#define MAX_VCF_ID_LEN 100 // including terminating nul
extern const char *vcf_header_get_VCF_ID_by_dict_id (DictId dict_id, bool must_exist);
extern VcfVersion vcf_header_get_version (void);

// Samples stuff
extern void vcf_seg_samples_initialize (void);
extern const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13, const char *backup_luft_samples, uint32_t backup_luft_samples_len);
extern void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCF *vb);
#define IS_TRIVAL_FORMAT_SUBFIELD ((!recon_len || (recon_len==1 && *recon=='.')) && dict_id_is_vcf_format_sf (ctx->dict_id))

// INFO/SF stuff
extern void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *recon, int32_t recon_len);
extern void vcf_piz_TOPLEVEL_cb_insert_INFO_SF (VBlockVCFP vcf_vb);
extern void vcf_seg_INFO_SF_one_sample (VBlockVCF *vb);
extern void vcf_seg_info_subfields (VBlockVCF *vb, const char *info_str, unsigned info_len);
extern void vcf_finalize_seg_info (VBlockVCF *vb);

// Refalt stuff
extern void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len);
extern void vcf_refalt_seg_other_REFALT (VBlockVCFP vb, DidIType did_i, LiftOverStatus ostatus, bool is_xstrand, unsigned add_bytes);
extern LiftOverStatus vcf_refalt_lift (VBlockVCFP vb, const ZipDataLineVCF *dl, bool xstrand);
typedef enum { EQUALS_NEITHER, EQUALS_REF, EQUALS_ALT, EQUALS_MISSING } RefAltEquals;
RefAltEquals vcf_refalt_oref_equals_ref_or_alt (char oref, char ref, const char *alt, unsigned alt_len, bool is_xstrand);
extern bool vcf_refalt_piz_is_variant_snp (VBlockP vb);
extern bool vcf_refalt_piz_is_variant_indel (VBlockP vb);

// Liftover Zip
extern void vcf_lo_zip_initialize (void);
extern void vcf_lo_append_rejects_file (VBlockP vb, Coords coord);
extern void vcf_lo_seg_generate_INFO_DVCF (VBlockVCFP vb, ZipDataLineVCF *dl);
extern void vcf_lo_set_rollback_point (VBlockVCFP vb);
extern void vcf_lo_seg_rollback_and_reject (VBlockVCFP vb, LiftOverStatus ostatus, Context *ctx);
extern LiftOverStatus vcf_lo_get_liftover_coords (VBlockVCFP vb, PosType pos, WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i); // out
extern void vcf_lo_seg_INFO_LUFT_and_PRIM (VBlockVCFP vb, DictId dict_id, const char *value, int value_len);
extern void vcf_lo_seg_INFO_REJX (VBlockVCFP vb, DictId dict_id, const char *value, int value_len);
extern bool vcf_lo_seg_cross_render_to_primary (VBlockVCFP vb, ContextP ctx, const char *this_value, unsigned this_value_len, char *modified_snip, unsigned *modified_snip_len);

#define vcf_set_ostatus(ostatus) ctx_set_last_value ((VBlockP)(vb), &(vb)->contexts[VCF_oSTATUS], (int64_t)(ostatus))

// Liftover Piz
extern void vcf_lo_piz_TOPLEVEL_cb_filter_line (VBlockP vb);

#define VCF_ERR_PREFIX { progress_newline; fprintf (stderr, "Error %s:%u in variant %s=%.*s %s=%"PRId64": ", __FUNCTION__, __LINE__, (vb->line_coords == DC_PRIMARY ? "CHROM" : "oCHROM"), vb->chrom_name_len, vb->chrom_name, (vb->line_coords == DC_PRIMARY ? "POS" : "oPOS"), vb->last_int (vb->line_coords == DC_PRIMARY ? VCF_POS : VCF_oPOS)); }
#define ASSVCF(condition, format, ...) do { if (!(condition)) { VCF_ERR_PREFIX; fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSVCF0(condition, msg)        ASSVCF ((condition), msg "%s", "")
#define WARNVCF(format, ...)           do { if (!flag.quiet)  { VCF_ERR_PREFIX; fprintf (stderr, format "\n", __VA_ARGS__); } } while(0)

#define REJECT(ostatus, reason, ...)              do { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s CHROM=%.*s POS=%"PRId64" REF=%.*s%s " reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN (100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } while(0)
#define REJECT_MAPPING(reason)                    do { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s CHROM=%.*s POS=%"PRId64 " " reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0]); return; } while(0)
#define REJECT_SUBFIELD(ostatus, ctx, reason,...) do { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s CHROM=%.*s POS=%"PRId64" " reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], __VA_ARGS__); \
                                                       vcf_lo_seg_rollback_and_reject (vb, (ostatus), (ctx)); } while(0)
#define LIFTOK(ostatus, reason, ...) do { if (flag.show_lift && !flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s CHROM=%.*s POS=%"PRId64" REF=%.*s%s " reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN (100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } while(0)
#define LIFTOK0(ostatus, reason) LIFTOK(ostatus, reason "%s", "")
#define REJECTIF(condition, ostatus, reason, ...) do { if (condition) { if (!flag.rejects_coord) bufprintf (vb, &vb->rejects_report, "%s CHROM=%.*s POS=%"PRId64" REF=%.*s%s " reason "\n", dvcf_status_names[ostatus], vb->chrom_name_len, vb->chrom_name, DATA_LINE (vb->line_i)->pos[0], MIN (100, vb->main_ref_len), vb->main_refalt, vb->main_ref_len > 100 ? "..." :"", __VA_ARGS__); return (ostatus); } } while(0)
#define REJECTIF0(condition, ostatus, reason)     do { if (condition) REJECT (ostatus, reason "%s", ""); } while (0)

#endif

