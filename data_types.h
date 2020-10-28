// ------------------------------------------------------------------
//   header.h
//   Copyrigh8t (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DATA_TYPES_INCLUDED
#define DATA_TYPES_INCLUDED

#include "genozip.h"
#include "md5.h"
#include "sections.h"

#include "vcf.h"
#include "sam.h"
#include "fastq.h"
#include "fasta.h"
#include "gff3.h"
#include "me23.h"
#include "reference.h"
#include "txtfile.h"

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
typedef enum { DT_NONE=-1, // used in the code logic, never written to the file
               DT_REF=0, DT_VCF=1, DT_SAM=2, DT_FASTQ=3, DT_FASTA=4, DT_GFF3=5, DT_ME23=6, // these values go into SectionHeaderGenozipHeader.data_type
               DT_BAM=7, NUM_DATATYPES 
             } DataType; 

typedef bool (*PizSpecialCtxHandler)(VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len, LastValueTypeP new_value);

typedef struct DataTypeProperties {
    
    // ZIP properties and functions
    const char *name;
    enum {NO_RA, RA} has_random_access;
    unsigned line_height; // how many actual txt file lines are in one seg line (seg lines are counted in lines.len)
    unsigned (*sizeof_vb)(void);
    unsigned (*sizeof_zip_dataline)(void);

    // TXT file properties
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    char txt_header_1st_char;  // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)
    int32_t (*is_header_done) (void);  // header length if header read is complete, 0 if not complete yet + sets lines.len
    uint32_t (*unconsumed) (VBlockP);  // called by I/O thread called by txtfile_read_vblock to get the length of uncosumed txt beyond the vblock txt
    bool (*zip_inspect_txt_header) (BufferP txt_header); // called by I/O thread to verify the txt header. returns false if this txt file should be skipped

    // ZIP callbacks
    void (*zip_initialize)(void);      // called by I/O thread when after the txt header is read
    void (*zip_read_one_vb)(VBlockP);  // called by I/O thread after reading txt of one vb into vb->txt_data
    void (*seg_initialize)(VBlockP);   // called by Compute thread at the beginning of Seg
    const char *(*seg_txt_line)(VBlockP, const char *field_start_line, bool *has_13);  // Called by Compute thread to Seg one line
    void (*seg_finalize)(VBlockP);     // called by Compute thread at the end of Seg
    void (*compress)(VBlockP);

    // PIZ callbacks
    void (*piz_initialize)(void); // called after reconstructing the txt header and before compute threads
    bool (*piz_read_one_vb)(VBlockP, SectionListEntryP);
    bool (*is_skip_secetion)(VBlockP, SectionType, DictId);
    void (*reconstruct_seq)(VBlockP, ContextP, const char *, unsigned);
    bool (*container_filter)(VBlockP, DictId, ConstContainerP, unsigned rep, int item); // returns true if rep, item should be reconstructed. if item=-1, this applies to the entire rep  

    unsigned num_special;
    PizSpecialCtxHandler special[10];

    // VBlock functions
    void (*release_vb)(VBlockP);
    void (*destroy_vb)(VBlockP);
    void (*cleanup_memory)(VBlockP);

    // misc properties and functions
    const char *show_stats_line_name; // the header displayed in --show-stats
    const char *stat_dict_types[3]; // the dictionary type displayed in --show-stats
} DataTypeProperties;

#define usz(type) ((unsigned)sizeof(type))
#define DATA_TYPE_PROPERTIES { \
/*    name       has_ra ht sizeof_vb     sizeof_zip_dataline  txt_headr 1st  is_header_done       unconsumed        zip_inspect_txt_header, zip_initialize        zip_read_one_vb        seg_initialize        seg_txt_line        seg_finalize,       compress                  piz_initialize         piz_read_one_vb        is_skip_secetion           reconstruct_seq            container_filter num_special        special        release_vb          destroy_vb           cleanup_memory          show_sections_line stat_dict_types                 */ \
    { "REFERENCE", RA,  1, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  NULL,                fasta_unconsumed, NULL,                   ref_make_ref_init,    NULL,                  fasta_seg_initialize, fasta_seg_txt_line, NULL,               ref_make_create_range,    NULL,                  NULL,                  NULL,                      NULL,                      NULL,             0,                 {},            fast_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "VCF",     RA,    1, vcf_vb_size,  vcf_vb_zip_dl_size,  HDR_MUST, '#', NULL,                NULL,             vcf_inspect_txt_header, NULL,                 NULL,                  vcf_seg_initialize,   vcf_seg_txt_line,   vcf_seg_finalize,   NULL,                     NULL,                  NULL,                  vcf_piz_is_skip_section,   NULL,                      vcf_piz_filter,   NUM_VCF_SPECIAL,   VCF_SPECIAL,   vcf_vb_release_vb,  vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "SAM",     RA,    1, sam_vb_size,  sam_vb_zip_dl_size,  HDR_OK,   '@', NULL,                NULL,             sam_inspect_txt_header, sam_zip_initialize,   NULL,                  sam_seg_initialize,   sam_seg_txt_line,   sam_seg_finalize,   NULL,                     NULL,                  NULL,                  sam_piz_is_skip_section,   sam_piz_reconstruct_seq,   NULL,             NUM_SAM_SPECIAL,   SAM_SPECIAL,   sam_vb_release_vb,  sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",   NO_RA, 4, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  NULL,                fastq_unconsumed, NULL,                   fastq_zip_initialize, fastq_zip_read_one_vb, fastq_seg_initialize, fastq_seg_txt_line, fastq_seg_finalize, NULL,                     NULL,                  fastq_piz_read_one_vb, fastq_piz_is_skip_section, fastq_piz_reconstruct_seq, fastq_piz_filter, 0,                 {},            fast_vb_release_vb, NULL,                NULL,                   "Entries",         { "FIELD", "DESC",   "ERROR!" } }, \
    { "FASTA",   RA,    1, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  NULL,                fasta_unconsumed, NULL,                   NULL,                 NULL,                  fasta_seg_initialize, fasta_seg_txt_line, fasta_seg_finalize, NULL,                     fasta_piz_initialize,  fasta_piz_read_one_vb, fasta_piz_is_skip_section, NULL,                      NULL,             NUM_FASTA_SPECIAL, FASTA_SPECIAL, fast_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "GVF",     RA,    1, 0,            0,                   HDR_OK,   '#', NULL,                NULL,             NULL,                   NULL,                 NULL,                  gff3_seg_initialize,  gff3_seg_txt_line,  gff3_seg_finalize,  NULL,                     NULL,                  NULL,                  NULL,                      NULL,                      NULL,             0,                 {},            NULL,               NULL,                NULL,                   "Sequences",       { "FIELD", "ATTRS",  "ITEMS"  } }, \
    { "23ANDME", RA,    1, 0,            0,                   HDR_OK,   '#', NULL,                NULL,             NULL,                   NULL,                 NULL,                  me23_seg_initialize,  me23_seg_txt_line,  me23_seg_finalize,  NULL,                     NULL,                  NULL,                  NULL,                      NULL,                      NULL,             0,                 {},            NULL,               NULL,                NULL,                   "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }, \
    { "BAM",     RA,    0, sam_vb_size,  sam_vb_zip_dl_size,  HDR_MUST, -1,  bam_is_header_done,  bam_unconsumed,   NULL,                   sam_zip_initialize,   NULL,                  bam_seg_initialize,   bam_seg_txt_line,   sam_seg_finalize,   NULL,                     NULL,                  NULL,                  NULL,                      NULL,                      NULL,             NUM_SAM_SPECIAL,   SAM_SPECIAL,   sam_vb_release_vb,  sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
}  
#define DATA_TYPE_FUNCTIONS_DEFAULT /* only applicable to (some) functions */ \
    { "DEFAULT", 0,     0, 0,            0,                   0,        0,   def_is_header_done,  def_unconsumed,   0,                      0,                    0,                     0,                    0,                  0,                  0,                        0,                     0,                     0,                         0,                         0,                0,                 {},            0,                  0,                   0,                      "",                { }                             }

extern DataTypeProperties dt_props[NUM_DATATYPES], dt_props_def;

#define DT_(src,prop) (dt_props[src->data_type].prop)
#define DTP(prop)  DT_((vb),prop)
#define DTPZ(prop) DT_(z_file,prop) 
#define DTPT(prop) DT_(txt_file,prop)

#define ASSERT_DT_FUNC(src,prop) /* use this to assert the a function exists if it is mandatory */ \
    ASSERT (DT_(src,prop) || dt_props_def.prop, "Error in %s:%u: undefined callback function " #src "->" #prop " for %s", \
            __FUNCTION__, __LINE__, dt_name(src->data_type))

// expession evaluates to def_value if neither the src nor the default function exists
#define DT_FUNC_OPTIONAL(src,prop,def_value) (!DT_(src,prop) && !dt_props_def.prop) ? def_value : (DT_(src,prop) ? DT_(src,prop) : dt_props_def.prop) // caller adds function arguments here

// expession evaluates to 0 if neither the src nor the default function exists
// !!WARNING!!: to use DT_FUNC in an expression it MUST be enclosed in parathesis: (DT_FUNC(vb,func)(args))
#define DT_FUNC(src,prop) DT_FUNC_OPTIONAL (src, prop, 0)

// Fields - the CHROM field MUST be the first field (because of mtf_copy_reference_contig_to_chrom_ctx)
typedef enum { REF_CONTIG, NUM_REF_FIELDS } RefFields;
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_SAMPLES, VCF_EOL, VCF_TOPLEVEL, NUM_VCF_FIELDS } VcfFields;
typedef enum { SAM_RNAME, SAM_QNAME, SAM_FLAG, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_RNEXT, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL, SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, SAM_QUAL, SAM_DOMQRUNS, SAM_EOL, BAM_BIN, SAM_TOPLEVEL, SAM_TOPLEVEL_BIN, NUM_SAM_FIELDS } SamFields;
typedef enum { FASTQ_CONTIG /* copied from reference */, FASTQ_DESC, FASTQ_E1L, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND, FASTQ_E2L, FASTQ_QUAL, FASTQ_DOMQRUNS, FASTQ_TOPLEVEL, NUM_FASTQ_FIELDS } FastqFields;
typedef enum { FASTA_CONTIG, FASTA_LINEMETA, FASTA_EOL, FASTA_DESC, FASTA_COMMENT, FASTA_SQBITMAP, FASTA_NONREF, FASTA_NONREF_X, FASTA_GPOS, FASTA_STRAND, FASTA_TOPLEVEL, NUM_FASTA_FIELDS } FastaFields;
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS, GFF3_EOL, GFF3_TOPLEVEL, NUM_GFF3_FIELDS } Gff3Fields;
typedef enum { ME23_CHROM, ME23_POS, ME23_ID, ME23_GENOTYPE, ME23_EOL, ME23_TOPLEVEL, NUM_ME23_FIELDS } Me23Fields;  

#define MAX_NUM_FIELDS_PER_DATA_TYPE MAX ((int) NUM_REF_FIELDS,    \
                                     MAX ((int) NUM_VCF_FIELDS,    \
                                     MAX ((int) NUM_SAM_FIELDS,    \
                                     MAX ((int) NUM_FASTQ_FIELDS,  \
                                     MAX ((int) NUM_FASTA_FIELDS,  \
                                     MAX ((int) NUM_GFF3_FIELDS,   \
                                          (int) NUM_ME23_FIELDS     ))))))

#define MAX_DICTS (MAX_SUBFIELDS*2 + MAX_NUM_FIELDS_PER_DATA_TYPE)  
//#if MAX_DICTS > 253 // 254 is for future use and 255 is DID_I_NONE
//#error "MAX_DICTS cannot go beyond 253"
//#endif

typedef struct DataTypeFields {
    unsigned num_fields;
    int pos, info, nonref, qual, eol; // the fields, or -1 if this data type doesn't have them
    char *names[MAX_NUM_FIELDS_PER_DATA_TYPE]; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
} DataTypeFields;

#define CHROM (DidIType)0 // chrom is always the first field

#define TOPLEVEL "TOPLEVEL"
#define TOPLEVEL_BIN "TOPLVBIN"
#define DATA_TYPE_FIELDS { \
/* num_fields        pos         info        nonref        qual        eol        names (including extend fields) - max 8 characters - 2 first chars must be unique within each data type (for dict_id_to_did_i_map) */ \
  {NUM_REF_FIELDS,   -1,         -1,         -1,           -1,         -1,        { "CONTIG", }, }, \
  {NUM_VCF_FIELDS,   VCF_POS,    VCF_INFO,   -1,           -1,         VCF_EOL,   { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLES", "EOL", TOPLEVEL } }, \
  {NUM_SAM_FIELDS,   SAM_POS,    -1,         SAM_NONREF,   SAM_QUAL,   SAM_EOL,   { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "QUAL", "DOMQRUNS", "EOL", "BIN", TOPLEVEL, TOPLEVEL_BIN } }, \
  {NUM_FASTQ_FIELDS, -1,         -1,         FASTQ_NONREF, FASTQ_QUAL, FASTQ_E1L, { "CONTIG", "DESC", "E1L", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "E2L", "QUAL", "DOMQRUNS", TOPLEVEL } }, \
  {NUM_FASTA_FIELDS, -1,         -1,         FASTA_NONREF, -1,         FASTA_EOL, { "CONTIG", "LINEMETA", "EOL", "DESC", "COMMENT", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", TOPLEVEL } }, \
  {NUM_GFF3_FIELDS,  GFF3_START, GFF3_ATTRS, -1,           -1,         GFF3_EOL,  { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRS", "EOL", TOPLEVEL } }, \
  {NUM_ME23_FIELDS,  ME23_POS,   -1,         -1,           -1,         ME23_EOL,  { "CHROM", "POS", "ID", "GENOTYPE", "EOL", TOPLEVEL } }, \
  {NUM_SAM_FIELDS,   SAM_POS,    -1,         SAM_NONREF,   SAM_QUAL,   SAM_EOL,   { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "QUAL", "DOMQRUNS", "EOL", "BIN", TOPLEVEL, TOPLEVEL_BIN } }, \
}
extern DataTypeFields dt_fields[NUM_DATATYPES];
#define DTF(prop)  (dt_fields[vb->      data_type].prop)
#define DTFZ(prop) (dt_fields[z_file->  data_type].prop)

// list of ctx who's local data is compressed via a callback function
#define LOCAL_GET_LINE_CALLBACKS {  \
    VCF_LOCAL_GET_LINE_CALLBACKS    \
    GFF3_LOCAL_GET_LINE_CALLBACKS   \
    SAM_LOCAL_GET_LINE_CALLBACKS    \
    BAM_LOCAL_GET_LINE_CALLBACKS    \
    FASTQ_LOCAL_GET_LINE_CALLBACKS  \
    FASTA_LOCAL_GET_LINE_CALLBACKS  \
}

// aliases - these are used only in PIZ, as a way for multiple dict_id's to get access to the same data storage
// on the ZIP side, the data is just placed directly in the primary ctx
#define DICT_ID_ALIASES { \
    VCF_DICT_ID_ALIASES   \
    SAM_DICT_ID_ALIASES   \
    BAM_DICT_ID_ALIASES   \
    FASTQ_DICT_ID_ALIASES \
    FASTA_DICT_ID_ALIASES \
    GFF3_DICT_ID_ALIASES  \
    ME23_DICT_ID_ALIASES  \
}

// possible transformations on Container items in case Container is reconstructed
// with "transform" (eg for binary files)
typedef enum __attribute__ ((__packed__)) { // 1 byte
    TRS_NONE=0,  // keep as is

    // Tranform textual integer to Little Endian binary
    TRS_LTEN_U8=1, TRS_LTEN_I8=2, TRS_LTEN_U16=3, TRS_LTEN_I16=4, TRS_LTEN_U32=5, TRS_LTEN_I32=6,

    // SAM-to-BAM tranformations
    TRS_SAM_RNAME    =10, // reconstructs the b250 index or -1 if "*"
    TRS_SAM_POS      =11, // textual 1-based POS to Little Endian U32 0-based POS. 
    TRS_SAM_SEQ      =12, // textual SEQ to BAM-format SEQ
    TRS_SAM_QUAL     =13, // textual QUAL to BAM-format QUAL
    TRS_SAM_OPTIONAL =14, // transform prefixes in Optional Container from SAM to BAM format
    TRS_SAM_FLOAT    =15, // SAM_SPECIAL_FLOAT snip to Little Endian 32bit float
} ContainerItemTransform;

extern const char *dt_name (DataType data_type);

#endif
