// ------------------------------------------------------------------
//   data_types.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef DATA_TYPES_INCLUDED
#define DATA_TYPES_INCLUDED

#include "dict_id.h"
#include "genozip.h"
#include "digest.h"
#include "sections.h"
#include "flags.h"
#include "reference.h"
#include "txtfile.h"
#include "container.h"
#include "dict_id_gen.h"

// data type files
#include "vcf.h"
#include "sam.h"
#include "fastq.h"
#include "fasta.h"
#include "gff3.h"
#include "me23.h"
#include "generic.h"
#include "phylip.h"
#include "chain.h"
#include "kraken.h"

#ifndef MAX // can't use MAX_ in a global variable
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

#define MAX_NUM_SPECIAL MAX ((int)NUM_GNRIC_SPECIAL, \
                        MAX ((int)NUM_VCF_SPECIAL,   \
                        MAX ((int)NUM_SAM_SPECIAL,   \
                        MAX ((int)NUM_CHAIN_SPECIAL, \
                             (int)NUM_FASTA_SPECIAL))))

typedef TRANSLATOR_FUNC ((*PizTranslator));
#define MAX_NUM_TRANSLATORS MAX (NUM_VCF_TRANS,      \
                            MAX (NUM_SAM_TRANS,      \
                            MAX (NUM_ME23_TRANS,     \
                            MAX (NUM_PHYLIP_TRANS,   \
                                 NUM_FASTA_TRANS))))

typedef struct DataTypeProperties {
    
    // ZIP properties and functions
    const char *name;
    bool is_binary; 
    DataType bin_type;                 // the binary equivalent of this textual file - exists for every data type that might have genozip_header.txt_is_bin=true
    unsigned line_height;              // how many actual txt file lines are in one seg line (seg lines are counted in lines.len). drop_lines in container_reconstruct_do needs to match the maximum.
    unsigned (*sizeof_vb)(DataType dt);
    unsigned (*sizeof_zip_dataline)(void);

    // TXT HEDEAR stuff
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    const char *hdr_contigs;        // format of contigs in the txtheader
    char txt_header_1st_char;          // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)
    int32_t (*is_header_done) (bool is_eof);  // ZIP: header length if header read is complete, -1 if not complete yet + sets lines.len
    int32_t (*unconsumed) (VBlockP, uint32_t first_i, int32_t *i);  // called by main thread called by txtfile_get_unconsumed_to_pass_up to get the length of unconsumed txt to pass to next vb. returns -1 if first_i is too high and it needs more data.
    bool (*inspect_txt_header) (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags); // called by main thread to verify the txt header. returns false if this txt file should be skipped

    // ZIP callbacks
    void (*zip_initialize)(void);       // called by main thread when after the txt header is read
    void (*zip_finalize)(void);         // called by main thread after each txt file compressing is done
    void (*zip_read_one_vb)(VBlockP);   // called by main thread after reading txt of one vb into vb->txt_data
    void (*zip_after_compute)(VBlockP); // called by main thread after completing compute thread of VB
    bool (*zip_dts_flag)(void);         // called to set FlagsGenozipHeader.dt_specific
    void (*seg_initialize)(VBlockP);    // called by Compute thread at the beginning of Seg
    const char *(*seg_txt_line)(VBlockP, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13);  // Called by Compute thread to Seg one line
    bool (*seg_is_small)(ConstVBlockP, DictId); // returns true if dict_id is known to have less than ~ 30K values (i.e. its context needs only a small hash table)
    void (*seg_finalize)(VBlockP);      // called by Compute thread at the end of Seg
    void (*compress)(VBlockP);

    // PIZ callbacks
    void (*piz_header_init)(void);// called at the beginning of reconstruction of every txt header that is reconstructed
    bool (*piz_initialize)(void); // called at the beginning of each output txt_file (after the global sections have been read already)
    void (*piz_finalize)(void);   // called by main thread after each z_file reconstruction is done
    bool (*piz_read_one_vb)(VBlockP, Section); // called by main thread after all sections of a VB have been read, before dispatching the compute thread
    bool (*is_skip_section)(VBlockP, SectionType, DictId);
    void (*reconstruct_seq)(VBlockP, ContextP, const char *, unsigned);
    CONTAINER_FILTER_FUNC ((*container_filter)); // called for containers as defined in the container
    CONTAINER_CALLBACK ((*container_cb)); // callback called after of every repeat in containers with .callback=true

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
    const char *show_stats_line_name; // the header displayed in --stats
    const char *dtype_names[4];       // the dictionary type names
} DataTypeProperties;

#define usz(type) ((unsigned)sizeof(type))
#define DATA_TYPE_PROPERTIES { \
/*    name         is_bin bin_type    ht sizeof_vb      sizeof_zip_dataline   txt_headr hdr_contigs     1st  is_header_done       unconsumed        inspect_txt_header,     zip_initialize          zip_finalize         zip_read_one_vb        zip_after_compute          zip_dts_flag          seg_initialize          seg_txt_line          seg_is_small,         seg_finalize,         compress                  piz_header_init      piz_initialize          piz_finalize         piz_read_one_vb          is_skip_section             reconstruct_seq            container_filter       container_cb              num_special        special        num_trans        translators        release_vb           destroy_vb           cleanup_memory          show_sections_line dtype_names                 */ \
    { "REFERENCE", false, DT_NONE,    1, fasta_vb_size, fasta_vb_zip_dl_size, HDR_NONE, NULL,           -1,  NULL,                fasta_unconsumed, NULL,                   ref_make_ref_init,      ref_make_finalize,   NULL,                  ref_make_after_compute,    NULL,                 ref_make_seg_initialize,fasta_seg_txt_line,   fasta_seg_is_small,   NULL,                 ref_make_create_range,    NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            0,               {},                fasta_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "VCF",       false, DT_NONE,    1, vcf_vb_size,   vcf_vb_zip_dl_size,   HDR_MUST, VCF_CONTIG_FMT, '#', NULL,                NULL,             vcf_inspect_txt_header, vcf_zip_initialize,     NULL,                vcf_zip_read_one_vb,   vcf_zip_after_compute,     NULL,                 vcf_seg_initialize,     vcf_seg_txt_line,     vcf_seg_is_small,     vcf_seg_finalize,     NULL,                     vcf_header_piz_init, NULL,                   NULL,                vcf_piz_read_one_vb,     vcf_piz_is_skip_section,    NULL,                      vcf_piz_filter,        vcf_piz_container_cb,     NUM_VCF_SPECIAL,   VCF_SPECIAL,   NUM_VCF_TRANS,   VCF_TRANSLATORS,   vcf_vb_release_vb,   vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT", "BOTH" } }, \
    { "SAM",       false, DT_BAM,     1, sam_vb_size,   sam_vb_zip_dl_size,   HDR_OK,   SAM_CONTIG_FMT, '@', NULL,                NULL,             sam_header_inspect,     sam_zip_initialize,     sam_header_finalize, NULL,                  NULL,                      sam_zip_dts_flag,     sam_seg_initialize,     sam_seg_txt_line,     sam_seg_is_small,     sam_seg_finalize,     NULL,                     NULL,                NULL,                   sam_header_finalize, sam_piz_read_one_vb,     sam_piz_is_skip_section,    sam_reconstruct_seq,       sam_piz_filter,        sam_piz_container_cb,     NUM_SAM_SPECIAL,   SAM_SPECIAL,   NUM_SAM_TRANS,   SAM_TRANSLATORS,   sam_vb_release_vb,   sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",     false, DT_NONE,    4, fastq_vb_size, fastq_vb_zip_dl_size, HDR_NONE, NULL,           -1,  NULL,                fastq_unconsumed, NULL,                   fastq_zip_initialize,   NULL,                fastq_zip_read_one_vb, NULL,                      fastq_zip_dts_flag,   fastq_seg_initialize,   fastq_seg_txt_line,   fastq_seg_is_small,   fastq_seg_finalize,   NULL,                     NULL,                NULL,                   NULL,                fastq_piz_read_one_vb,   fastq_piz_is_skip_section,  fastq_reconstruct_seq,     fastq_piz_filter,      fastq_piz_container_cb,   0,                 {},            0,               {},                fastq_vb_release_vb, fastq_vb_destroy_vb, NULL,                   "Sequences",       { "FIELD", "DESC",   "ERROR!" } }, \
    { "FASTA",     false, DT_NONE,    1, fasta_vb_size, fasta_vb_zip_dl_size, HDR_NONE, NULL,           -1,  NULL,                fasta_unconsumed, NULL,                   fasta_zip_initialize,   NULL,                NULL,                  NULL,                      NULL,                 fasta_seg_initialize,   fasta_seg_txt_line,   fasta_seg_is_small,   fasta_seg_finalize,   NULL,                     NULL,                fasta_piz_initialize,   NULL,                fasta_piz_read_one_vb,   fasta_piz_is_skip_section,  NULL,                      fasta_piz_filter,      NULL,                     NUM_FASTA_SPECIAL, FASTA_SPECIAL, NUM_FASTA_TRANS, FASTA_TRANSLATORS, fasta_vb_release_vb, fasta_vb_destroy_vb, NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "GFF3",      false, DT_NONE,    1, 0,             0,                    HDR_OK,   NULL,           '#', NULL,                NULL,             NULL,                   NULL,                   NULL,                NULL,                  NULL,                      NULL,                 gff3_seg_initialize,    gff3_seg_txt_line,    gff3_seg_is_small,    gff3_seg_finalize,    NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      gff3_piz_filter,       NULL,                     0,                 {},            0,               {},                NULL,                NULL,                NULL,                   "Sequences",       { "FIELD", "ATTRS",  "ENST"  } }, \
    { "23ANDME",   false, DT_NONE,    1, 0,             0,                    HDR_MUST, NULL,           '#', NULL,                NULL,             me23_header_inspect,    NULL,                   NULL,                NULL,                  NULL,                      NULL,                 me23_seg_initialize,    me23_seg_txt_line,    me23_seg_is_small,    me23_seg_finalize,    NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            NUM_ME23_TRANS,  ME23_TRANSLATORS,  NULL,                NULL,                NULL,                   "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }, \
    { "BAM",       true,  DT_NONE,    0, sam_vb_size,   sam_vb_zip_dl_size,   HDR_MUST, SAM_CONTIG_FMT, -1,  bam_is_header_done,  bam_unconsumed,   sam_header_inspect,     sam_zip_initialize,     sam_header_finalize, NULL,                  NULL,                      sam_zip_dts_flag,     bam_seg_initialize,     bam_seg_txt_line,     sam_seg_is_small,     sam_seg_finalize,     NULL,                     NULL,                NULL,                   sam_header_finalize, sam_piz_read_one_vb,     NULL,                       NULL,                      sam_piz_filter,        sam_piz_container_cb,     NUM_SAM_SPECIAL,   SAM_SPECIAL,   NUM_SAM_TRANS,   SAM_TRANSLATORS,   sam_vb_release_vb,   sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "BCF",       false, DT_NONE,    1, vcf_vb_size,   vcf_vb_zip_dl_size,   HDR_MUST, VCF_CONTIG_FMT, '#', NULL,                NULL,             vcf_inspect_txt_header, vcf_zip_initialize,     NULL,                vcf_zip_read_one_vb,   vcf_zip_after_compute,     NULL,                 vcf_seg_initialize,     vcf_seg_txt_line,     vcf_seg_is_small,     vcf_seg_finalize,     NULL,                     vcf_header_piz_init, NULL,                   NULL,                vcf_piz_read_one_vb,     vcf_piz_is_skip_section,    NULL,                      vcf_piz_filter,        vcf_piz_container_cb,     NUM_VCF_SPECIAL,   VCF_SPECIAL,   NUM_VCF_TRANS,   VCF_TRANSLATORS,   vcf_vb_release_vb,   vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT", "BOTH" } }, \
    { "GENERIC",   true,  DT_GENERIC, 0, 0,             0,                    HDR_NONE, NULL,           -1,  NULL,                generic_unconsumed, NULL,                 NULL,                   NULL,                NULL,                  NULL,                      NULL,                 NULL,                   NULL,                 generic_seg_is_small, generic_seg_finalize, NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     NUM_GNRIC_SPECIAL, GNRIC_SPECIAL, 0,               {},                NULL,                NULL,                NULL,                   "N/A",             { "FIELD", "ERROR!", "ERROR!" } }, \
    { "PHYLIP",    false, DT_NONE,    1, 0,             phy_vb_zip_dl_size,   HDR_MUST, NULL,           -1,  phy_is_header_done,  NULL,             phy_header_inspect,     NULL,                   NULL,                NULL,                  NULL,                      NULL,                 phy_seg_initialize,     phy_seg_txt_line,     phy_seg_is_small,     phy_seg_finalize,     NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            NUM_PHYLIP_TRANS,PHYLIP_TRANSLATORS,NULL,                NULL,                NULL,                   "Sequences",       { "FIELD", "ERROR!", "ERROR!" } }, \
    { "CHAIN",     false, DT_NONE,    1, chain_vb_size, 0,                    HDR_NONE, NULL,           -1,  NULL,                chain_unconsumed, NULL,                   chain_zip_initialize,   NULL,                NULL,                  NULL,                      chain_zip_dts_flag,   chain_seg_initialize,   chain_seg_txt_line,   chain_seg_is_small,   chain_seg_finalize,   NULL,                     NULL,                chain_piz_initialize,   NULL,                NULL,                    NULL,                       NULL,                      chain_piz_filter,      NULL,                     NUM_CHAIN_SPECIAL, CHAIN_SPECIAL, 0,               {},                chain_vb_release_vb, NULL,                NULL,                   "Alignment sets",  { "FIELD", "ERROR!", "ERROR!" } }, \
    { "KRAKEN",    false, DT_NONE,    1, 0,             0,                    HDR_NONE, NULL,           -1,  NULL,                NULL,             NULL,                   kraken_zip_initialize,  NULL,                NULL,                  kraken_zip_after_compute,  NULL,                 kraken_seg_initialize,  kraken_seg_txt_line,  kraken_seg_is_small,  kraken_seg_finalize,  NULL,                     NULL,                kraken_piz_initialize,  NULL,                NULL,                    kraken_piz_is_skip_section, NULL,                      NULL,                  kraken_piz_container_cb,  0,                 {},            0,               {},                NULL,                NULL,                NULL,                   "Sequences",       { "FIELD", "QNAME",  "SEQLEN" } }, \
}  
#define DATA_TYPE_FUNCTIONS_DEFAULT /* only applicable to (some) functions */ \
    { "DEFAULT",   false, DT_NONE,    0, def_vb_size,   0,                    0,        0,              0,   def_is_header_done,  def_unconsumed,   0,                      0,                      0,                   0,                     0,                         0,                    0,                      0,                    0,                    0,                    0,                        0,                   0,                      NULL,                0,                       0,                          0,                         container_no_filter,   0,                        0,                 {},            0,               {},                0,                   0,                   0,                      "",                {                             } }

extern DataTypeProperties dt_props[NUM_DATATYPES], dt_props_def;

#define DT_(src,prop) (dt_props[src->data_type].prop)
#define DTP(prop)  DT_((vb),prop)
#define DTPZ(prop) DT_(z_file,prop) 
#define DTPT(prop) DT_(txt_file,prop)

#define ASSERT_DT_FUNC(src,prop) /* use this to assert the a function exists if it is mandatory */ \
    ASSERT (DT_(src,prop) || dt_props_def.prop, "undefined callback function " #src "->" #prop " for %s", dt_name(src->data_type))

// expession evaluates to def_value if neither the src nor the default function exists
#define DT_FUNC_OPTIONAL(src,prop,def_value) (!DT_(src,prop) && !dt_props_def.prop) ? def_value : (DT_(src,prop) ? DT_(src,prop) : dt_props_def.prop) // caller adds function arguments here

// expession evaluates to 0 if neither the src nor the default function exists
// !!WARNING!!: to use DT_FUNC in an expression it MUST be enclosed in parathesis: (DT_FUNC(vb,func)(args))
#define DT_FUNC(src,prop) DT_FUNC_OPTIONAL (src, prop, 0)

#define MAX_DICTS 2048 // 11 bits -matches CONTAINER_FIELDS.nitems_lo+nitems_hi

typedef struct DataTypeFields {
    const unsigned num_fields;
    const DidIType pos, nonref, prim_chrom, luft_chrom, luft_pos; // the fields, or DID_I_NONE if this data type doesn't have them
    DictId eol, toplevel;
    struct { 
        DictId dict_id; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
        STR (tag_name);
    } predefined[MAX_NUM_FIELDS_PER_DATA_TYPE]; // consumed by ctx_initialize_predefined_ctxs 
} DataTypeFields;
#define CHROM (DidIType)0 // chrom is always the first field

#define TOPLEVEL "TOPLEVEL"

#define DATA_TYPE_FIELDS { \
/* num_fields        pos              nonref        prim_chrom,     luft_chrom,     luft_pos,        eol             toplevel             predefined        */ \
  {NUM_REF_FIELDS,   DID_I_NONE,      DID_I_NONE,   CHROM,          DID_I_NONE,     DID_I_NONE,      {0},            {0},                 REF_PREDEFINED    }, \
  {NUM_VCF_FIELDS,   VCF_POS,         DID_I_NONE,   CHROM,          VCF_oCHROM,     VCF_oPOS,        {_VCF_EOL},     {_VCF_TOPLEVEL},     VCF_PREDEFINED    }, \
  {NUM_SAM_FIELDS,   SAM_POS,         SAM_NONREF,   CHROM,          DID_I_NONE,     DID_I_NONE,      {_SAM_EOL},     {_SAM_TOPLEVEL},     SAM_PREDEFINED    }, \
  {NUM_FASTQ_FIELDS, DID_I_NONE,      FASTQ_NONREF, CHROM,          DID_I_NONE,     DID_I_NONE,      {_FASTQ_E1L},   {_FASTQ_TOPLEVEL},   FASTQ_PREDEFINED  }, \
  {NUM_FASTA_FIELDS, DID_I_NONE,      FASTA_NONREF, CHROM,          DID_I_NONE,     DID_I_NONE,      {_FASTA_EOL},   {_FASTA_TOPLEVEL},   FASTA_PREDEFINED  }, \
  {NUM_GFF3_FIELDS,  GFF3_START,      DID_I_NONE,   CHROM,          DID_I_NONE,     DID_I_NONE,      {_GFF3_EOL},    {_GFF3_TOPLEVEL},    GFF3_PREDEFINED   }, \
  {NUM_ME23_FIELDS,  ME23_POS,        DID_I_NONE,   CHROM,          DID_I_NONE,     DID_I_NONE,      {_ME23_EOL},    {_ME23_TOPLEVEL},    ME23_PREDEFINED   }, \
  {NUM_SAM_FIELDS,   SAM_POS,         SAM_NONREF,   CHROM,          DID_I_NONE,     DID_I_NONE,      {_SAM_EOL},     {_SAM_TOP2BAM},      SAM_PREDEFINED    }, \
  {NUM_VCF_FIELDS,   VCF_POS,         DID_I_NONE,   CHROM,          VCF_oCHROM,     VCF_oPOS,        {_VCF_EOL},     {_VCF_TOPLEVEL},     VCF_PREDEFINED    }, \
  {NUM_GNRIC_FIELDS, DID_I_NONE,      DID_I_NONE,   DID_I_NONE,     DID_I_NONE,     DID_I_NONE,      {0},            {_GNRIC_TOPLEVEL},   GNRIC_PREDEFINED  }, \
  {NUM_PHY_FIELDS,   DID_I_NONE,      DID_I_NONE,   DID_I_NONE,     DID_I_NONE,     DID_I_NONE,      {_PHY_EOL},     {_PHY_TOPLEVEL},     PHY_PREDEFINED    }, \
  {NUM_CHAIN_FIELDS, CHAIN_STARTPRIM, DID_I_NONE,   CHAIN_NAMEPRIM, CHAIN_NAMELUFT, CHAIN_STARTLUFT, {_CHAIN_EOL},   {_CHAIN_TOPLEVEL},   CHAIN_PREDEFINED  }, \
  {NUM_KRAKEN_FIELDS,DID_I_NONE,      DID_I_NONE,   DID_I_NONE,     DID_I_NONE,     DID_I_NONE,      {_KRAKEN_EOL},  {_KRAKEN_TOPLEVEL},  KRAKEN_PREDEFINED }, \
}

extern DataTypeFields dt_fields[NUM_DATATYPES];
#define DTF(prop)  (dt_fields[vb->      data_type].prop)
#define DTFZ(prop) (dt_fields[z_file->  data_type].prop)
#define DTFT(prop) (dt_fields[txt_file->data_type].prop)

// list of ctx who's local data is compressed via a callback function
#define LOCAL_GET_LINE_CALLBACKS {  \
    VCF_LOCAL_GET_LINE_CALLBACKS    \
    GFF3_LOCAL_GET_LINE_CALLBACKS   \
    SAM_LOCAL_GET_LINE_CALLBACKS    \
    BAM_LOCAL_GET_LINE_CALLBACKS    \
    FASTQ_LOCAL_GET_LINE_CALLBACKS  \
    FASTA_LOCAL_GET_LINE_CALLBACKS  \
    PHY_LOCAL_GET_LINE_CALLBACKS    \
}

// aliases - these are used only in PIZ, as a way for multiple dict_id's to get access to the same data storage
// on the ZIP side, the data is just placed directly in the destination ctx
#define DICT_ID_ALIASES { \
    VCF_DICT_ID_ALIASES   \
    GFF3_DICT_ID_ALIASES  \
}

typedef struct DtTranslation { 
    DataType src_z_non_bin_dt; 
    uint8_t src_z_is_binary;  // same type as z_file.flags
    DataType dst_txt_dt;
    DictId toplevel;
    float factor;             // specifies the size of vb.txt_data allocated in PIZ vs the size recorded in SectionHeaderVbHeader.recon_size.
                              // It is fully allocated in advance and cannot be extended.
    TXTHEADER_TRANSLATOR ((*txtheader_translator));
    bool trans_containers;    // whether to invoke translators in containers
    bool is_src_dt;           // this translation dst is actually the source file dt that was compressed (of which we have the MD5)
    bool (*is_translation)(VBlockP vb); // if non-NULL - called to determine if this VB should be translated
} DtTranslation;

#define TRANSLATIONS { \
    /*                           non-bin-dt  binary dst_txt_dt toplevel        factor  txtheader_transl.        trans_con. src_dt is_trans. */ \
    /* SAM to BAM           */ { DT_SAM,     false, DT_BAM,    { _SAM_TOP2BAM },       2,   sam_header_sam2bam,  true,     false, NULL }, /* BAM is expected to be smaller, but in edge cases numbers "1\t" -> uint32 and QUAL="*"->0xff X seq_len */ \
    /* SAM to FASTQ         */ { DT_SAM,     false, DT_FASTQ,  { _SAM_TOP2FQ },        1,   txtheader_sam2fq,    true,     false, NULL }, /* sizes of SEQ, QUAL, QNAME the same */ \
    /* SAM to FASTQ extend. */ { DT_SAM,     false, DT_FASTQ,  { _SAM_TOP2FQEX },      1.1, txtheader_sam2fq,    true,     false, NULL }, /* very slightly bigger due to added + */ \
    /* Binary SAM to SAM    */ { DT_SAM,     true,  DT_SAM,    { _SAM_TOPLEVEL },      3,   sam_header_bam2sam,  false,    false, NULL }, /* BAM stores sequences in 2x and numbers in 1-2.5x */ \
    /* Binary SAM to BAM    */ { DT_SAM,     true,  DT_BAM,    { _SAM_TOP2BAM },       2,   NULL,                true,     true,  NULL }, /* BAM is expected to be smaller, but in edge cases numbers "1\t" -> uint32 */ \
    /* Binary BAM to FASTQ  */ { DT_SAM,     true,  DT_FASTQ,  { _SAM_TOP2FQ },        1.5, txtheader_sam2fq,    true,     false, NULL }, /* sizes of QNAME, QUAL the same, SEQ x2 */ \
    /* 23andMe to VCF       */ { DT_ME23,    false, DT_VCF,    { _ME23_TOP2VCF },      4,   txtheader_me232vcf,  true,     false, NULL }, \
    /* FASTA to Phylip      */ { DT_FASTA,   false, DT_PHYLIP, { _FASTA_TOPLEVEL },    1.1, txtheader_fa2phy,    true,     false, NULL }, \
    /* Phylip to FASTA      */ { DT_PHYLIP,  false, DT_FASTA,  { _PHY_TOP2FASTA },     1.1, txtheader_phy2fa,    true,     false, NULL }, \
    /* VCF: genocat --luft  */ { DT_VCF,     false, DT_VCF,    { _VCF_TOPLUFT },       1,   NULL,                true,     false, vcf_vb_is_luft }, /* Length of --luft is exactly recon_size_luft*/ \
    /* genocat --taxid      */ { DT_KRAKEN,  false, DT_NONE,   { _KRAKEN_TOP2TAXID },  1,   NULL,                true,     false, kraken_is_translation }, \
}

extern const DtTranslation dt_get_translation (VBlockP vb);
extern DataType dt_get_txt_dt (DataType dt);
extern const char *dt_name (DataType data_type);

#endif
