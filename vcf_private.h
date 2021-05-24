// ------------------------------------------------------------------
//   vcf_private.h
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
typedef struct VBlockVCF {

    VBLOCK_COMMON_FIELDS

    // charactaristics of the data
    uint16_t ploidy;           // ZIP only
            
    // used for segging FORMAT/GT
    Context *gt_ctx;
    uint32_t gt_prev_ploidy;
    char gt_prev_phase;
    
    // used for segging INFO
    Buffer info_items;              // Seg: INFO items of the line being segged

    const char *main_refalt;        // used by vcf_refalt_check_oref and vcf_seg_INFO_BaseCounts, set by vcf_seg_txt_line
    unsigned main_ref_len, main_alt_len;

    // INFO/SF stuff
    Context *sf_ctx;
    enum { USE_SF_UNKNOWN, USE_SF_YES, USE_SF_NO } use_special_sf;
    Buffer sf_txt, sf_snip; // INFO/SF data as it appears in the snip being constructed

    // INFO/END
    uint32_t last_end_line_i;       // PIZ: last line on which INFO/END was encountered
    
    // FORMAT/AD
    #define MAX_ARG_ARRAY_ITEMS 36 // same as MAX_COMPOUND_COMPONENTS
    int64_t ad_values[MAX_ARG_ARRAY_ITEMS];

    // dictionaries stuff 
    Buffer format_mapper_buf;       // ZIP only: an array of type Container - one entry per entry in vb->contexts[VCF_FORMAT].nodes   
    Buffer format_contexts;         // ZIP only: an array of MAX_SUBFIELDS * format_mapper_buf.len of ContextP

    // used by CODEC_HAPM (for VCF haplotype matrix) 
    Context *hapmat_index_ctx; 
    Buffer hapmat_helper_index_buf; // ZIP: used by codec_hapmat_count_alt_alleles 
    Buffer hapmat_columns_data;     // used by codec_hapmat_piz_get_one_line 
    Buffer hapmat_column_of_zeros;  // used by codec_hapmat_piz_calculate_columns   
    Buffer hapmat_one_array;        // one line or column 
} VBlockVCF;

typedef VBlockVCF *VBlockVCFP;

// Liftover stuff

// this is part of the file format (in the oSTATUS context) and CANNOT BE CHANGED (but can be extended)
typedef enum {
    // Algorithm interim statuses not outputed to the genozip file
    LO_NA               = -1, // not a Liftover item - never written z_file
    LO_UNKNOWN          = 0,  // Status of this Liftover line not known yet - never written z_file

    // OK statuses - this are internal to genozip and z_file.oStatus - not expressed in the txt file
    LO_OK               = 1,  // internal progress step in ZIP - never written z_file
    LO_OK_REF_SAME      = 2,
    LO_OK_REF_ALT_SWTCH = 3,
    LO_OK_FFU4          = 4,
    LO_OK_FFU5          = 5,
    LO_OK_FFU6          = 6,
    LO_OK_FFU7          = 7,
    LO_OK_FFU8          = 8,

    // Rejection reasons defined in the "Dual coordinate VCF" spec 
    LO_REJECTED         = 9,  // generic rejection status
    LO_NO_CHROM         = 10, // chain file doesn't contain a qName which is the primary chrom as 
    LO_NO_MAPPING       = 11, // chain file doesn't contain a mapping for the primary POS
    LO_REF_TOO_LONG     = 12, // Genozip only supports --chain (even without change) for REF of length 1
    LO_REF_LONG_CHANGE  = 13, // Genozip only supports REF change, for REF of length 1
    LO_REF_LONG_XSTRAND = 14, // Genozip can only reverse complement an REF of length 1 
    LO_ALT_LONG_XSTRAND = 15, // Genozip can only reverse complement an ALT of length 1 
    LO_ALT_LONG_SWITCH  = 16, // Genozip can only switch REF<>ALT for an ALT of length 1
    LO_REF_CHNG_NOT_ALT = 17, // Genozip can only liftover a changed REF if its switched with the ALT
    LO_ADDED_VARIANT    = 18, // Variant was added to file after it was already lifted over

    // these must be in the same order as VCF_TRANSLATORS, starting from VCF2VCF_G
    LO_INFO             = 19,
    LO_FORMAT           = 20,
    
    NUM_LO_STATUSES
} LiftOverStatus;

#define LO_IS_REJECTED(ost) ((ost) >= LO_REJECTED) // note: this condition works also for unrecognized reject strings (that an have index >= NUM_LO_STATUSES)
#define LO_IS_OK(ost) ((ost) >= LO_OK && (ost) < LO_REJECTED)

#define last_ostatus (ctx_has_value_in_line_(vb, &vb->contexts[VCF_oSTATUS]) ? vb->last_index (VCF_oSTATUS) : LO_UNKNOWN)
#define last_ostatus_name_piz (last_ostatus < vb->contexts[VCF_oSTATUS].word_list.len ? ctx_get_words_snip (&vb->contexts[VCF_oSTATUS], last_ostatus) : "Invalid oStatus")

// fields and their order of INFO/LIFTOVER and INFO/LIFTBACK containers (part of the file format, see vcf_lo_zip_initialize)
typedef enum { IL_CHROM, IL_POS, IL_REF, IL_XSTRAND, NUM_IL_FIELDS } InfoLiftFields; 

// NOTE: The rejection strings should not change or vcf_lo_seg_INFO_REJX won't identify oStatus in rejected lines when parsing old dual coordinate files (it will fallback to LO_REJECTED)
extern const char *dvcf_status_names[];
#define DVCF_STATUS_NAMES { /* for display esthetics - max 14 characters - note: these strings are defined in the dual-coordinate specification */\
    "UNKNOWN", \
    "OK", "OKRefSame", "OkRefAltSwitch", "OKffu4", "OKffu5", "OKffu6", "OKffu7", "OKffu8", \
    "Rejected", "NoChrom", "NoMapping", \
    "RefTooLong", "RefLongChange", "RefLongXstrand", "AltLongXstrand", "AltLongSwitch", "RefChngeNotAlt", \
    "AddedVariant", "INFO", "FORMAT" \
}

// Liftover - header keys
#define HK_GENOZIP_CMD "##genozip_command="
#define HK_CHAIN       "##chain="
#define HK_DC          "##dual_coordinates"
#define HK_DC_PRIMARY  HK_DC"=PRIMARY"
#define HK_DC_LUFT     HK_DC"=LUFT"
#define HK_RENDERALG_ATTR "RenderAlg"
#define KH_INFO_LUFT   "##INFO=<ID=" INFO_LUFT ",Number=4,Type=String,Description=\"Info for rendering variant in LUFT coords. See " WEBSITE_COORDS "\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_PRIM   "##INFO=<ID=" INFO_PRIM ",Number=4,Type=String,Description=\"Info for rendering variant in PRIMARY coords\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_LREJ   "##INFO=<ID=" INFO_LREJ ",Number=1,Type=String,Description=\"Reason variant was rejected for LUFT coords\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"
#define KH_INFO_PREJ   "##INFO=<ID=" INFO_PREJ ",Number=1,Type=String,Description=\"Reason variant was rejected for PRIMARY coords\",Source=\"genozip\",Version=\"%s\"," HK_RENDERALG_ATTR "=NONE>"

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

// Samples stuff
extern void vcf_seg_samples_initialize (void);
extern const char *vcf_seg_samples (VBlockVCF *vb, ZipDataLineVCF *dl, int32_t *len, char *next_field, bool *has_13, const char *backup_luft_samples, uint32_t backup_luft_samples_len);
extern void vcf_seg_FORMAT_GT_complete_missing_lines (VBlockVCF *vb);

// INFO/SF stuff
extern void vcf_piz_GT_cb_calc_INFO_SF (VBlockVCFP vcf_vb, unsigned rep, char *recon, int32_t recon_len);
extern void vcf_piz_TOPLEVEL_cb_insert_INFO_SF (VBlockVCFP vcf_vb);
extern void vcf_seg_INFO_SF_one_sample (VBlockVCF *vb, unsigned sample_i);
extern void vcf_seg_info_subfields (VBlockVCF *vb, const char *info_str, unsigned info_len);
extern void vcf_finalize_seg_info (VBlockVCF *vb);

// Refalt stuff
extern void vcf_refalt_seg_main_ref_alt (VBlockVCFP vb, const char *ref, unsigned ref_len, const char *alt, unsigned alt_len);
extern LiftOverStatus vcf_refalt_check_oref (VBlockVCFP vb, const ZipDataLineVCF *dl, bool xstrand, unsigned *oref_len);
typedef enum { EQUALS_NEITHER, EQUALS_REF2, EQUALS_ALT, EQUALS_MISSING } RefAltEquals;
RefAltEquals vcf_refalt_ref_equals_ref2_or_alt (char ref, char ref2, const char *alt, unsigned alt_len, bool is_xstrand);

// Liftover Zip
extern void vcf_lo_zip_initialize (void);
extern void vcf_lo_append_rejects_file (VBlockP vb, Coords coord);
extern void vcf_lo_seg_generate_INFO_DVCF (VBlockVCFP vb, ZipDataLineVCF *dl);
extern void vcf_lo_set_rollback_point (VBlockVCFP vb);
extern void vcf_lo_seg_rollback_and_reject (VBlockVCFP vb, LiftOverStatus ostatus, Context *ctx);
extern LiftOverStatus vcf_lo_get_liftover_coords (VBlockVCFP vb, WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i); // out
extern void vcf_lo_seg_INFO_LUFT_and_PRIM (VBlockVCFP vb, DictId dict_id, const char *value, int value_len);
extern void vcf_lo_seg_INFO_REJX (VBlockVCFP vb, DictId dict_id, const char *value, int value_len);
extern bool vcf_lo_seg_cross_render_to_primary (VBlockVCFP vb, ContextP ctx, const char *this_value, unsigned this_value_len, char *modified_snip, unsigned *modified_snip_len);

#define vcf_set_ostatus(ostatus) ctx_set_last_value ((VBlockP)(vb), &(vb)->contexts[VCF_oSTATUS], (int64_t)(ostatus))

// Liftover Piz
extern void vcf_lo_piz_TOPLEVEL_cb_filter_line (VBlockP vb);

#define ASSVCF(condition, format, ...) do { if (!(condition)) { progress_newline; fprintf (stderr, "Error %s:%u in variant CHROM=%.*s POS=%"PRId64": ",  __FUNCTION__, __LINE__, vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS)); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSVCF0(condition, msg)        do { if (!(condition)) { progress_newline; fprintf (stderr, "Error %s:%u in variant CHROM=%.*s POS=%"PRId64": %s\n",  __FUNCTION__, __LINE__, vb->chrom_name_len, vb->chrom_name, vb->last_int (VCF_POS), (msg)); exit_on_error(true); }} while(0)

#endif

