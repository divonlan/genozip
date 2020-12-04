// ------------------------------------------------------------------
//   header.h
//   Copyrigh8t (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DATA_TYPES_INCLUDED
#define DATA_TYPES_INCLUDED

#include "genozip.h"
#include "digest.h"
#include "sections.h"
#include "flags.h"
#include "reference.h"
#include "txtfile.h"

// data type files
#include "vcf.h"
#include "sam.h"
#include "fastq.h"
#include "fasta.h"
#include "gff3.h"
#include "me23.h"
#include "generic.h"

#define MAX_NUM_SPECIAL MAX ((int)NUM_VCF_SPECIAL, \
                        MAX ((int)NUM_SAM_SPECIAL, \
                             (int)NUM_FASTA_SPECIAL))

typedef SPECIAL_RECONSTRUCTOR ((*PizSpecialReconstructor));
typedef TRANSLATOR_FUNC ((*PizTranslator));
#define MAX_NUM_TRANSLATORS NUM_SAM_TRANS

typedef struct DataTypeProperties {
    
    // ZIP properties and functions
    const char *name;
    bool is_binary; 
    DataType bin_type;                 // the binary equivalent of this textual file - exists for every data type that might have genozip_header.txt_is_bin=true
    enum {NO_RA, RA} has_random_access;
    unsigned line_height;              // how many actual txt file lines are in one seg line (seg lines are counted in lines.len)
    unsigned (*sizeof_vb)(void);
    unsigned (*sizeof_zip_dataline)(void);

    // TXT file properties
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    char txt_header_1st_char;          // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)
    int32_t (*is_header_done) (void);  // header length if header read is complete, -1 if not complete yet + sets lines.len
    int32_t (*unconsumed) (VBlockP, uint32_t first_i, int32_t *i);  // called by I/O thread called by txtfile_get_unconsumed_to_pass_up to get the length of unconsumed txt to pass to next vb. returns -1 if first_i is too high and it needs more data.
    bool (*inspect_txt_header) (BufferP txt_header); // called by I/O thread to verify the txt header. returns false if this txt file should be skipped

    // ZIP callbacks
    void (*zip_initialize)(void);      // called by I/O thread when after the txt header is read
    void (*zip_finalize)(void);        // called by I/O thread after each txt file compressing is done
    void (*zip_read_one_vb)(VBlockP);  // called by I/O thread after reading txt of one vb into vb->txt_data
    void (*seg_initialize)(VBlockP);   // called by Compute thread at the beginning of Seg
    const char *(*seg_txt_line)(VBlockP, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13);  // Called by Compute thread to Seg one line
    void (*seg_finalize)(VBlockP);     // called by Compute thread at the end of Seg
    void (*compress)(VBlockP);

    // PIZ callbacks
    void (*piz_initialize)(void); // called after reconstructing the txt header and before compute threads
    void (*piz_finalize)(void);   // called by I/O thread after each z_file reconstruction is done
    bool (*piz_read_one_vb)(VBlockP, ConstSectionListEntryP);
    bool (*is_skip_secetion)(VBlockP, SectionType, DictId);
    void (*reconstruct_seq)(VBlockP, ContextP, const char *, unsigned);
    CONTAINER_FILTER_FUNC ((*container_filter)); // called for containers as defined in the container

    // Special reconstruction functions corresponding to SNIP_SPECIAL snips
    unsigned num_special;
    PizSpecialReconstructor special[MAX_NUM_SPECIAL];

    // Container item translators, invoked after reconstruction, for translating from one data type to another
    unsigned num_translators;
    PizTranslator translator[MAX_NUM_TRANSLATORS];
    
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
/*    name         is_bin bin_type has_ra ht sizeof_vb     sizeof_zip_dataline  txt_headr 1st  is_header_done       unconsumed        inspect_txt_header,     zip_initialize        zip_finalize      zip_read_one_vb        seg_initialize        seg_txt_line        seg_finalize,       compress                  piz_initialize         piz_finalize         piz_read_one_vb        is_skip_secetion           reconstruct_seq            container_filter       num_special        special        num_trans        translators        release_vb          destroy_vb           cleanup_memory          show_sections_line stat_dict_types                 */ \
    { "REFERENCE", false, DT_NONE, RA,    1, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  NULL,                fasta_unconsumed, NULL,                   ref_make_ref_init,    NULL,             NULL,                  fasta_seg_initialize, fasta_seg_txt_line, NULL,               ref_make_create_range,    NULL,                  NULL,                NULL,                  NULL,                      NULL,                      NULL,                  0,                 {},            0,               {},                fast_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "VCF",       false, DT_NONE, RA,    1, vcf_vb_size,  vcf_vb_zip_dl_size,  HDR_MUST, '#', NULL,                NULL,             vcf_inspect_txt_header, NULL,                 NULL,             NULL,                  vcf_seg_initialize,   vcf_seg_txt_line,   vcf_seg_finalize,   NULL,                     NULL,                  NULL,                NULL,                  vcf_piz_is_skip_section,   NULL,                      vcf_piz_filter,        NUM_VCF_SPECIAL,   VCF_SPECIAL,   0,               {},                vcf_vb_release_vb,  vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "SAM",       false, DT_BAM,  RA,    1, sam_vb_size,  sam_vb_zip_dl_size,  HDR_OK,   '@', NULL,                NULL,             sam_header_inspect,     NULL,                 sam_header_finalize, NULL,               sam_seg_initialize,   sam_seg_txt_line,   sam_seg_finalize,   NULL,                     NULL,                  sam_header_finalize, NULL,                  sam_piz_is_skip_section,   sam_piz_reconstruct_seq,   sam_piz_sam2fq_filter, NUM_SAM_SPECIAL,   SAM_SPECIAL,   NUM_SAM_TRANS,   SAM_TRANSLATORS,   sam_vb_release_vb,  sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",     false, DT_NONE, NO_RA, 4, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  NULL,                fastq_unconsumed, NULL,                   fastq_zip_initialize, NULL,             fastq_zip_read_one_vb, fastq_seg_initialize, fastq_seg_txt_line, fastq_seg_finalize, NULL,                     NULL,                  NULL,                fastq_piz_read_one_vb, fastq_piz_is_skip_section, fastq_piz_reconstruct_seq, fastq_piz_filter,      0,                 {},            0,               {},                fast_vb_release_vb, NULL,                NULL,                   "Entries",         { "FIELD", "DESC",   "ERROR!" } }, \
    { "FASTA",     false, DT_NONE, RA,    1, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  NULL,                fasta_unconsumed, NULL,                   NULL,                 NULL,             NULL,                  fasta_seg_initialize, fasta_seg_txt_line, fasta_seg_finalize, NULL,                     fasta_piz_initialize,  NULL,                fasta_piz_read_one_vb, fasta_piz_is_skip_section, NULL,                      fasta_piz_filter,      NUM_FASTA_SPECIAL, FASTA_SPECIAL, 0,               {},                fast_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "GVF",       false, DT_NONE, RA,    1, 0,            0,                   HDR_OK,   '#', NULL,                NULL,             NULL,                   NULL,                 NULL,             NULL,                  gff3_seg_initialize,  gff3_seg_txt_line,  gff3_seg_finalize,  NULL,                     NULL,                  NULL,                NULL,                  NULL,                      NULL,                      NULL,                  0,                 {},            0,               {},                NULL,               NULL,                NULL,                   "Sequences",       { "FIELD", "ATTRS",  "ITEMS"  } }, \
    { "23ANDME",   false, DT_NONE, RA,    1, 0,            0,                   HDR_MUST, '#', NULL,                NULL,             me23_header_inspect,    NULL,                 NULL,             NULL,                  me23_seg_initialize,  me23_seg_txt_line,  me23_seg_finalize,  NULL,                     NULL,                  NULL,                NULL,                  NULL,                      NULL,                      NULL,                  0,                 {},            NUM_ME23_TRANS,  ME23_TRANSLATORS,  NULL,               NULL,                NULL,                   "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }, \
    { "BAM",       true,  DT_NONE, RA,    0, sam_vb_size,  sam_vb_zip_dl_size,  HDR_MUST, -1,  bam_is_header_done,  bam_unconsumed,   sam_header_inspect,     NULL,                 sam_header_finalize, NULL,               bam_seg_initialize,   bam_seg_txt_line,   sam_seg_finalize,   NULL,                     NULL,                  sam_header_finalize, NULL,                  NULL,                      NULL,                      sam_piz_sam2fq_filter, NUM_SAM_SPECIAL,   SAM_SPECIAL,   NUM_SAM_TRANS,   SAM_TRANSLATORS,   sam_vb_release_vb,  sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "BCF",       false, DT_NONE, RA,    1, vcf_vb_size,  vcf_vb_zip_dl_size,  HDR_MUST, '#', NULL,                NULL,             vcf_inspect_txt_header, NULL,                 NULL,             NULL,                  vcf_seg_initialize,   vcf_seg_txt_line,   vcf_seg_finalize,   NULL,                     NULL,                  NULL,                NULL,                  vcf_piz_is_skip_section,   NULL,                      vcf_piz_filter,        NUM_VCF_SPECIAL,   VCF_SPECIAL,   0,               {},                vcf_vb_release_vb,  vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "GENERIC",   true,  DT_GENERIC, NO_RA, 0, 0,         0,                   HDR_NONE, -1,  NULL,                generic_unconsumed, NULL,                 NULL,                 NULL,             NULL,                  NULL,                 NULL,               generic_seg_finalize,NULL,                    NULL,                  NULL,                NULL,                  NULL,                      NULL,                      NULL,                  NUM_GNRIC_SPECIAL, GNRIC_SPECIAL, 0,               {},                NULL,               NULL,                NULL,                   "N/A",             { "FIELD", "ERROR!", "ERROR!" } }, \
    { "PHYLIP",    false, DT_NONE, RA,    1, fast_vb_size, fast_vb_zip_dl_size, HDR_MUST, -1,  NULL,                fasta_unconsumed, NULL,                   NULL,                 NULL,             NULL,                  fasta_seg_initialize, fasta_seg_txt_line, fasta_seg_finalize, NULL,                     fasta_piz_initialize,  NULL,                fasta_piz_read_one_vb, fasta_piz_is_skip_section, NULL,                      fasta_piz_filter,      NUM_FASTA_SPECIAL, FASTA_SPECIAL, 0,               {}, fast_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
}  
#define DATA_TYPE_FUNCTIONS_DEFAULT /* only applicable to (some) functions */ \
    { "DEFAULT",   false, DT_NONE, NO_RA, 0, 0,            0,                   0,        0,   def_is_header_done,  def_unconsumed,   0,                      0,                    0,                0,                     0,                    0,                  0,                  0,                        0,                     NULL,                0,                     0,                         0,                         0,                     0,                 {},            0,               {},                0,                  0,                   0,                      "",                { }                             }

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

// Fields - the CHROM field MUST be the first field (because of ctx_copy_ref_contigs_to_zf)
typedef enum { REF_CONTIG, NUM_REF_FIELDS } RefFields;
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_SAMPLES, VCF_EOL, VCF_TOPLEVEL, NUM_VCF_FIELDS } VcfFields;
typedef enum { SAM_RNAME, SAM_QNAME, SAM_FLAG, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_RNEXT, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL, 
               SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, 
               SAM_QUAL, SAM_DOMQRUNS, SAM_EOL, SAM_BAM_BIN,
               SAM_TOPLEVEL, SAM_TOP2BAM, SAM_TOP2FQ, 
               SAM_E2_Z, SAM_2NONREF, SAM_N2ONREFX, SAM_2GPOS, SAM_S2TRAND, // E2 data - we put them in primary fields bc they need to be sequential - merging VBs dictionaries doesn't necessarily make sequentials ones in ZF
               SAM_U2_Z, SAM_D2OMQRUN,                                      // U2 data - same reason
               NUM_SAM_FIELDS } SamFields;
typedef enum { FASTQ_CONTIG /* copied from reference */, FASTQ_DESC, FASTQ_E1L, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND, FASTQ_E2L, FASTQ_QUAL, FASTQ_DOMQRUNS, FASTQ_TOPLEVEL, NUM_FASTQ_FIELDS } FastqFields;
typedef enum { FASTA_CONTIG, FASTA_LINEMETA, FASTA_EOL, FASTA_DESC, FASTA_COMMENT, FASTA_SQBITMAP, FASTA_NONREF, FASTA_NONREF_X, FASTA_GPOS, FASTA_STRAND, FASTA_TOPLEVEL, NUM_FASTA_FIELDS } FastaFields;
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS, GFF3_EOL, GFF3_TOPLEVEL, NUM_GFF3_FIELDS } Gff3Fields;
typedef enum { ME23_CHROM, ME23_POS, ME23_ID, ME23_GENOTYPE, ME23_EOL, ME23_TOPLEVEL, ME23_TOP2VCF, NUM_ME23_FIELDS } Me23Fields;  
typedef enum { GNRIC_DATA, GNRIC_TOPLEVEL, NUM_GNRIC_FIELDS } GenericFields;

#define MAX_NUM_FIELDS_PER_DATA_TYPE MAX ((int) NUM_REF_FIELDS,    \
                                     MAX ((int) NUM_VCF_FIELDS,    \
                                     MAX ((int) NUM_SAM_FIELDS,    \
                                     MAX ((int) NUM_FASTQ_FIELDS,  \
                                     MAX ((int) NUM_FASTA_FIELDS,  \
                                     MAX ((int) NUM_GFF3_FIELDS,   \
                                     MAX ((int) NUM_ME23_FIELDS,   \
                                          (int) NUM_GNRIC_FIELDS )))))))

#define MAX_DICTS (MAX_SUBFIELDS*2 + MAX_NUM_FIELDS_PER_DATA_TYPE)  
//#if MAX_DICTS > 253 // 254 is for future use and 255 is DID_I_NONE
//#error "MAX_DICTS cannot go beyond 253"
//#endif

typedef struct DataTypeFields {
    unsigned num_fields;
    int pos, info, nonref, eol, toplevel; // the fields, or -1 if this data type doesn't have them
    char *names[MAX_NUM_FIELDS_PER_DATA_TYPE]; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
} DataTypeFields;

#define CHROM (DidIType)0 // chrom is always the first field

#define TOPLEVEL "TOPLEVEL"
#define DATA_TYPE_FIELDS { \
/* num_fields        pos         info        nonref        eol        toplevel names (including extend fields) - max 8 characters - 2 first chars must be unique within each data type (for dict_id_to_did_i_map) */ \
  {NUM_REF_FIELDS,   -1,         -1,         -1,           1,         -1,             { "CONTIG", }, }, \
  {NUM_VCF_FIELDS,   VCF_POS,    VCF_INFO,   -1,           VCF_EOL,   VCF_TOPLEVEL,   { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLES", "EOL", TOPLEVEL } }, \
  {NUM_SAM_FIELDS,   SAM_POS,    -1,         SAM_NONREF,   SAM_EOL,   SAM_TOPLEVEL,   { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "QUAL", "DOMQRUNS", "EOL", "BAM_BIN", TOPLEVEL, "TOP2BAM", "TOP2FQ", "E2:Z", "2NONREF", "N2ONREFX", "2GPOS", "S2TRAND", "U2:Z", "D2OMQRUN" } }, \
  {NUM_FASTQ_FIELDS, -1,         -1,         FASTQ_NONREF, FASTQ_E1L, FASTQ_TOPLEVEL, { "CONTIG", "DESC", "E1L", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "E2L", "QUAL", "DOMQRUNS", TOPLEVEL } }, \
  {NUM_FASTA_FIELDS, -1,         -1,         FASTA_NONREF, FASTA_EOL, FASTA_TOPLEVEL, { "CONTIG", "LINEMETA", "EOL", "DESC", "COMMENT", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", TOPLEVEL } }, \
  {NUM_GFF3_FIELDS,  GFF3_START, GFF3_ATTRS, -1,           GFF3_EOL,  GFF3_TOPLEVEL,  { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRS", "EOL", TOPLEVEL } }, \
  {NUM_ME23_FIELDS,  ME23_POS,   -1,         -1,           ME23_EOL,  ME23_TOPLEVEL,  { "CHROM", "POS", "ID", "GENOTYPE", "EOL", TOPLEVEL, "TOP2VCF" } }, \
  {NUM_SAM_FIELDS,   SAM_POS,    -1,         SAM_NONREF,   SAM_EOL,   SAM_TOP2BAM,    { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "QUAL", "DOMQRUNS", "EOL", "BAM_BIN", TOPLEVEL, "TOP2BAM", "TOP2FQ", "E2:Z", "2NONREF", "N2ONREFX", "2GPOS", "S2TRAND", "U2:Z", "D2OMQRUN" } }, \
  {NUM_VCF_FIELDS,   VCF_POS,    VCF_INFO,   -1,           VCF_EOL,   VCF_TOPLEVEL,   { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLES", "EOL", TOPLEVEL } }, \
  {NUM_GNRIC_FIELDS, -1,        -1,          -1,           -1,        GNRIC_TOPLEVEL, { "DATA", TOPLEVEL } }, \
  {NUM_FASTA_FIELDS, -1,         -1,         FASTA_NONREF, FASTA_EOL, FASTA_TOPLEVEL, { "CONTIG", "LINEMETA", "EOL", "DESC", "COMMENT", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", TOPLEVEL } }, \
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

typedef struct DtTranslation { 
    DataType src_z_non_bin_dt; 
    uint8_t src_z_is_binary;  // same type as z_file.flags
    DataType dst_txt_dt;
    DidIType toplevel;
    float factor; 
    TXTHEADER_TRANSLATOR ((*txtheader_translator));
    bool trans_containers;    // whether to invoke translators in containers
    bool is_src_dt;           // this translation dst is actually the source file dt that was compressed (of which we have the MD5)
} DtTranslation;

#define TRANSLATIONS { \
    /*                          non-bin-dt binary dst_txt_dt toplevel     factor  txtheader_transl.    trans_con src_dt*/ \
    /* SAM to BAM          */ { DT_SAM,    false, DT_BAM,    SAM_TOP2BAM,    2,   txtheader_sam2bam,   true,     false }, /* BAM is expected to be smaller, but in edge cases numbers "1\t" -> uint32 and QUAL="*"->0xff X seq_len */ \
    /* SAM to FASTQ        */ { DT_SAM,    false, DT_FASTQ,  SAM_TOP2FQ,     1,   txtheader_sam2fq,    true,     false }, /* sizes of SEQ, QUAL, QNAME the same */ \
    /* Binary SAM to SAM   */ { DT_SAM,    true,  DT_SAM,    SAM_TOPLEVEL,   3,   txtheader_bam2sam,   false,    false }, /* BAM stores sequences in 2x and numbers in 1-2.5x */ \
    /* Binary SAM to BAM   */ { DT_SAM,    true,  DT_BAM,    SAM_TOP2BAM,    2,   NULL,                true,     true  }, /* BAM is expected to be smaller, but in edge cases numbers "1\t" -> uint32 */ \
    /* Binary BAM to FASTQ */ { DT_SAM,    true,  DT_FASTQ,  SAM_TOP2FQ,     1.5, txtheader_sam2fq,    true,     false }, /* sizes of QNAME, QUAL the same, SEQ x2 */ \
    /* 23andMe to VCF      */ { DT_ME23,   false, DT_VCF,    ME23_TOP2VCF,   4,   txtheader_me232vcf,  true,     false }, \
    /* FASTA to Phylip     */ { DT_FASTA,  false, DT_PHYLIP, FASTA_TOPLEVEL, 4,   txtheader_fa2phylip, true,     false }, \
}

extern const DtTranslation dt_get_translation (void);
extern DataType dt_get_txt_dt (DataType dt);
extern const char *dt_name (DataType data_type);

#endif
