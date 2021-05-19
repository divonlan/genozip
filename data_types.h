// ------------------------------------------------------------------
//   data_types.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DATA_TYPES_INCLUDED
#define DATA_TYPES_INCLUDED

#include "genozip.h"
#include "digest.h"
#include "sections.h"
#include "flags.h"
#include "reference.h"
#include "txtheader.h"
#include "txtfile.h"
#include "container.h"

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

#define MAX_NUM_SPECIAL MAX ((int)NUM_GNRIC_SPECIAL, \
                        MAX ((int)NUM_VCF_SPECIAL,   \
                        MAX ((int)NUM_SAM_SPECIAL,   \
                        MAX ((int)NUM_CHAIN_SPECIAL, \
                             (int)NUM_FASTA_SPECIAL))))

typedef TRANSLATOR_FUNC ((*PizTranslator));
#define MAX_NUM_TRANSLATORS MAX (NUM_VCF_TRANS,    \
                            MAX (NUM_SAM_TRANS,    \
                            MAX (NUM_ME23_TRANS,   \
                            MAX (NUM_PHYLIP_TRANS, \
                                 NUM_FASTA_TRANS))))

typedef struct DataTypeProperties {
    
    // ZIP properties and functions
    const char *name;
    bool is_binary; 
    DataType bin_type;                 // the binary equivalent of this textual file - exists for every data type that might have genozip_header.txt_is_bin=true
    enum {NO_RA, RA} has_random_access;
    unsigned line_height;              // how many actual txt file lines are in one seg line (seg lines are counted in lines.len). drop_lines in container_reconstruct_do needs to match the maximum.
    unsigned (*sizeof_vb)(void);
    unsigned (*sizeof_zip_dataline)(void);

    // TXT HEDEAR stuff
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    const char *header_contigs;        // format of contigs in the txtheader
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
    const char *show_stats_line_name; // the header displayed in --show-stats
    const char *stat_dict_types[3]; // the dictionary type displayed in --show-stats
} DataTypeProperties;

#define usz(type) ((unsigned)sizeof(type))
#define DATA_TYPE_PROPERTIES { \
/*    name         is_bin bin_type has_ra ht sizeof_vb      sizeof_zip_dataline   txt_headr header_contigs  1st  is_header_done       unconsumed        inspect_txt_header,     zip_initialize          zip_finalize        zip_read_one_vb        zip_after_compute          zip_dts_flag          seg_initialize          seg_txt_line          seg_is_small,         seg_finalize,         compress                  piz_header_init      piz_initialize          piz_finalize         piz_read_one_vb          is_skip_section             reconstruct_seq            container_filter       container_cb              num_special        special        num_trans        translators        release_vb           destroy_vb           cleanup_memory          show_sections_line stat_dict_types                 */ \
    { "REFERENCE", false, DT_NONE, RA,    1, fasta_vb_size, fasta_vb_zip_dl_size, HDR_NONE, NULL,           -1,  NULL,                fasta_unconsumed, NULL,                   ref_make_ref_init,      NULL,               NULL,                  NULL,                      NULL,                 fasta_seg_initialize,   fasta_seg_txt_line,   fasta_seg_is_small,   NULL,                 ref_make_create_range,    NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            0,               {},                fasta_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "VCF",       false, DT_NONE, RA,    1, vcf_vb_size,   vcf_vb_zip_dl_size,   HDR_MUST, VCF_CONTIG_FMT, '#', NULL,                NULL,             vcf_inspect_txt_header, vcf_zip_initialize,     txtheader_finalize, NULL,                  vcf_zip_after_compute,     NULL,                 vcf_seg_initialize,     vcf_seg_txt_line,     vcf_seg_is_small,     vcf_seg_finalize,     NULL,                     vcf_header_piz_init, NULL,                   txtheader_finalize,  vcf_piz_read_one_vb,     vcf_piz_is_skip_section,    NULL,                      vcf_piz_filter,        vcf_piz_container_cb,     NUM_VCF_SPECIAL,   VCF_SPECIAL,   NUM_VCF_TRANS,   VCF_TRANSLATORS,   vcf_vb_release_vb,   vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "SAM",       false, DT_BAM,  RA,    1, sam_vb_size,   sam_vb_zip_dl_size,   HDR_OK,   SAM_CONTIG_FMT, '@', NULL,                NULL,             sam_header_inspect,     sam_zip_initialize,     txtheader_finalize, NULL,                  NULL,                      sam_zip_dts_flag,     sam_seg_initialize,     sam_seg_txt_line,     sam_seg_is_small,     sam_seg_finalize,     NULL,                     NULL,                NULL,                   txtheader_finalize,  NULL,                    sam_piz_is_skip_section,    sam_reconstruct_seq,       NULL,                  sam_piz_container_cb,     NUM_SAM_SPECIAL,   SAM_SPECIAL,   NUM_SAM_TRANS,   SAM_TRANSLATORS,   sam_vb_release_vb,   sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",     false, DT_NONE, NO_RA, 4, fastq_vb_size, fastq_vb_zip_dl_size, HDR_NONE, NULL,           -1,  NULL,                fastq_unconsumed, NULL,                   fastq_zip_initialize,   NULL,               fastq_zip_read_one_vb, NULL,                      fastq_zip_dts_flag,   fastq_seg_initialize,   fastq_seg_txt_line,   fastq_seg_is_small,   fastq_seg_finalize,   NULL,                     NULL,                NULL,                   NULL,                fastq_piz_read_one_vb,   fastq_piz_is_skip_section,  fastq_reconstruct_seq,     fastq_piz_filter,      fastq_piz_container_cb,   0,                 {},            0,               {},                fastq_vb_release_vb, fastq_vb_destroy_vb, NULL,                   "Sequences",       { "FIELD", "DESC",   "ERROR!" } }, \
    { "FASTA",     false, DT_NONE, RA,    1, fasta_vb_size, fasta_vb_zip_dl_size, HDR_NONE, NULL,           -1,  NULL,                fasta_unconsumed, NULL,                   NULL,                   NULL,               NULL,                  NULL,                      NULL,                 fasta_seg_initialize,   fasta_seg_txt_line,   fasta_seg_is_small,   fasta_seg_finalize,   NULL,                     NULL,                fasta_piz_initialize,   NULL,                fasta_piz_read_one_vb,   fasta_piz_is_skip_section,  NULL,                      fasta_piz_filter,      NULL,                     NUM_FASTA_SPECIAL, FASTA_SPECIAL, NUM_FASTA_TRANS, FASTA_TRANSLATORS, fasta_vb_release_vb, fasta_vb_destroy_vb, NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "GVF",       false, DT_NONE, RA,    1, 0,             0,                    HDR_OK,   NULL,           '#', NULL,                NULL,             NULL,                   NULL,                   NULL,               NULL,                  NULL,                      NULL,                 gff3_seg_initialize,    gff3_seg_txt_line,    gff3_seg_is_small,    gff3_seg_finalize,    NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            0,               {},                NULL,                NULL,                NULL,                   "Sequences",       { "FIELD", "ATTRS",  "ITEMS"  } }, \
    { "23ANDME",   false, DT_NONE, RA,    1, 0,             0,                    HDR_MUST, NULL,           '#', NULL,                NULL,             me23_header_inspect,    NULL,                   NULL,               NULL,                  NULL,                      NULL,                 me23_seg_initialize,    me23_seg_txt_line,    me23_seg_is_small,    me23_seg_finalize,    NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            NUM_ME23_TRANS,  ME23_TRANSLATORS,  NULL,                NULL,                NULL,                   "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }, \
    { "BAM",       true,  DT_NONE, RA,    0, sam_vb_size,   sam_vb_zip_dl_size,   HDR_MUST, SAM_CONTIG_FMT, -1,  bam_is_header_done,  bam_unconsumed,   sam_header_inspect,     sam_zip_initialize,     txtheader_finalize, NULL,                  NULL,                      sam_zip_dts_flag,     bam_seg_initialize,     bam_seg_txt_line,     sam_seg_is_small,     sam_seg_finalize,     NULL,                     NULL,                NULL,                   txtheader_finalize,  NULL,                    NULL,                       NULL,                      NULL          ,        sam_piz_container_cb,     NUM_SAM_SPECIAL,   SAM_SPECIAL,   NUM_SAM_TRANS,   SAM_TRANSLATORS,   sam_vb_release_vb,   sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "BCF",       false, DT_NONE, RA,    1, vcf_vb_size,   vcf_vb_zip_dl_size,   HDR_MUST, VCF_CONTIG_FMT, '#', NULL,                NULL,             vcf_inspect_txt_header, vcf_zip_initialize,     txtheader_finalize, NULL,                  vcf_zip_after_compute,     NULL,                 vcf_seg_initialize,     vcf_seg_txt_line,     vcf_seg_is_small,     vcf_seg_finalize,     NULL,                     vcf_header_piz_init, NULL,                   txtheader_finalize,  vcf_piz_read_one_vb,     vcf_piz_is_skip_section,    NULL,                      vcf_piz_filter,        vcf_piz_container_cb,     NUM_VCF_SPECIAL,   VCF_SPECIAL,   NUM_VCF_TRANS,   VCF_TRANSLATORS,   vcf_vb_release_vb,   vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "GENERIC",   true,  DT_GENERIC, NO_RA, 0, 0,          0,                    HDR_NONE, NULL,           -1,  NULL,                generic_unconsumed, NULL,                 NULL,                   NULL,               NULL,                  NULL,                      NULL,                 NULL,                   NULL,                 generic_seg_is_small, generic_seg_finalize, NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     NUM_GNRIC_SPECIAL, GNRIC_SPECIAL, 0,               {},                NULL,                NULL,                NULL,                   "N/A",             { "FIELD", "ERROR!", "ERROR!" } }, \
    { "PHYLIP",    false, DT_NONE, NO_RA, 1, 0,             phy_vb_zip_dl_size,   HDR_MUST, NULL,           -1,  phy_is_header_done,  NULL,             phy_header_inspect,     NULL,                   NULL,               NULL,                  NULL,                      NULL,                 phy_seg_initialize,     phy_seg_txt_line,     phy_seg_is_small,     phy_seg_finalize,     NULL,                     NULL,                NULL,                   NULL,                NULL,                    NULL,                       NULL,                      NULL,                  NULL,                     0,                 {},            NUM_PHYLIP_TRANS,PHYLIP_TRANSLATORS,NULL,                NULL,                NULL,                   "Sequences",       { "FIELD", "ERROR!", "ERROR!" } }, \
    { "CHAIN",     false, DT_NONE, NO_RA, 1, 0,             0,                    HDR_NONE, NULL,           -1,  NULL,                chain_unconsumed, NULL,                   chain_zip_initialize,   NULL,               NULL,                  NULL,                      NULL,                 chain_seg_initialize,   chain_seg_txt_line,   chain_seg_is_small,   chain_seg_finalize,   NULL,                     NULL,                chain_piz_initialize,   NULL,                NULL,                    NULL,                       NULL,                      chain_piz_filter,      NULL,                     NUM_CHAIN_SPECIAL, CHAIN_SPECIAL, 0,               {},                NULL,                NULL,                NULL,                   "Alignment sets",  { "FIELD", "ERROR!", "ERROR!" } }, \
    { "KRAKEN",    false, DT_NONE, NO_RA, 1, 0,             0,                    HDR_NONE, NULL,           -1,  NULL,                NULL,             NULL,                   kraken_zip_initialize,  NULL,               NULL,                  kraken_zip_after_compute,  NULL,                 kraken_seg_initialize,  kraken_seg_txt_line,  kraken_seg_is_small,  kraken_seg_finalize,  NULL,                     NULL,                kraken_piz_initialize,  NULL,                NULL,                    kraken_piz_is_skip_section, NULL,                      NULL,                  kraken_piz_container_cb,  0,                 {},            0,               {},                NULL,                NULL,                NULL,                   "Sequences",       { "FIELD", "QNAME",  "ERROR!" } }, \
}  
#define DATA_TYPE_FUNCTIONS_DEFAULT /* only applicable to (some) functions */ \
    { "DEFAULT",   false, DT_NONE, NO_RA, 0, 0,             0,                    0,        0,              0,   def_is_header_done,  def_unconsumed,   0,                      0,                      0,                  0,                     0,                         0,                    0,                      0,                    0,                    0,                    0,                        0,                   0,                      NULL,                0,                       0,                          0,                         container_no_filter,   0,                        0,                 {},            0,               {},                0,                   0,                   0,                      "",                {                             } }

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

// Fields - the CHROM field MUST be the first field (because of ctx_build_zf_ctx_from_contigs)
typedef enum { REF_CONTIG, NUM_REF_FIELDS } RefFields;
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_SAMPLES, VCF_EOL, VCF_TOPLEVEL, 
               VCF_oCHROM, VCF_oPOS, VCF_oREFALT, VCF_oXSTRAND, VCF_COORDS, VCF_oSTATUS, VCF_COPYPOS, VCF_LIFT_REF, VCF_COPYSTAT, VCF_TOPLUFT, // Liftover data - must appear in same order in any data type that has it
               NUM_VCF_FIELDS } VcfFields;
typedef enum { SAM_RNAME, SAM_QNAME, SAM_FLAG, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_RNEXT, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL, 
               SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, 
               SAM_QUAL, SAM_DOMQRUNS, SAM_EOL, SAM_BAM_BIN,
               SAM_TOPLEVEL, SAM_TOP2BAM, SAM_TOP2FQ, 
               SAM_E2_Z, SAM_2NONREF, SAM_N2ONREFX, SAM_2GPOS, SAM_S2TRAND, // E2 data - we put them in primary fields bc they need to be sequential - merging VBs dictionaries doesn't necessarily make sequentials ones in ZF
               SAM_U2_Z, SAM_D2OMQRUN,                                      // U2 data - same reason
               SAM_TAXID,
               NUM_SAM_FIELDS } SamFields;
typedef enum { FASTQ_CONTIG /* copied from reference */, FASTQ_DESC, FASTQ_E1L, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND, FASTQ_E2L, FASTQ_QUAL, FASTQ_DOMQRUNS, FASTQ_TOPLEVEL, FASTQ_TAXID, NUM_FASTQ_FIELDS } FastqFields;
typedef enum { FASTA_CONTIG, FASTA_LINEMETA, FASTA_EOL, FASTA_DESC, FASTA_COMMENT, FASTA_NONREF, FASTA_NONREF_X, FASTA_TOPLEVEL, FASTA_TAXID, NUM_FASTA_FIELDS } FastaFields;
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS, GFF3_EOL, GFF3_TOPLEVEL, NUM_GFF3_FIELDS } Gff3Fields;
typedef enum { ME23_CHROM, ME23_POS, ME23_ID, ME23_GENOTYPE, ME23_EOL, ME23_TOPLEVEL, ME23_TOP2VCF, NUM_ME23_FIELDS } Me23Fields;  
typedef enum { GNRIC_DATA, GNRIC_TOPLEVEL, NUM_GNRIC_FIELDS } GenericFields;
typedef enum { PHY_ID, PHY_SEQ, PHY_EOL, PHY_TOPLEVEL, PHY_TOP2FASTA, NUM_PHY_FIELDS } PhyFields;
typedef enum { CHAIN_NAMEDST, CHAIN_STRNDDST, CHAIN_STARTDST, CHAIN_ENDDST, CHAIN_SIZEDST,
               CHAIN_NAMESRC, CHAIN_STRNDSRC, CHAIN_STARTSRC, CHAIN_ENDSRC, CHAIN_SIZESRC,
               CHAIN_CHAIN, CHAIN_SCORE, CHAIN_ID, CHAIN_SET, 
               CHAIN_SIZE, CHAIN_GAPS, CHAIN_EOL, CHAIN_TOPLEVEL, CHAIN_SEP, NUM_CHAIN_FIELDS } ChainFields;
typedef enum { KRAKEN_CU, KRAKEN_QNAME, KRAKEN_TAXID, KRAKEN_SEQLEN, KRAKEN_KMERS, KRAKEN_KMERTAX, KRAKEN_KMERLEN, 
               KRAKEN_EOL, KRAKEN_TOPLEVEL, KRAKEN_TOP2TAXID, NUM_KRAKEN_FIELDS } KrakenFields;

#define MAX_NUM_FIELDS_PER_DATA_TYPE MAX ((int) NUM_REF_FIELDS,     \
                                     MAX ((int) NUM_VCF_FIELDS,     \
                                     MAX ((int) NUM_SAM_FIELDS,     \
                                     MAX ((int) NUM_FASTQ_FIELDS,   \
                                     MAX ((int) NUM_FASTA_FIELDS,   \
                                     MAX ((int) NUM_GFF3_FIELDS,    \
                                     MAX ((int) NUM_ME23_FIELDS,    \
                                     MAX ((int) NUM_GNRIC_FIELDS,   \
                                     MAX ((int) NUM_CHAIN_FIELDS,   \
                                     MAX ((int) NUM_KRAKEN_FIELDS, \
                                          (int) NUM_PHY_FIELDS ))))))))))

#define MAX_DICTS 2048

typedef struct DataTypeFields {
    unsigned num_fields;
    DidIType pos, test_regions, nonref, ochrom, eol, toplevel; // the fields, or DID_I_NONE if this data type doesn't have them
    char *names[MAX_NUM_FIELDS_PER_DATA_TYPE]; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
} DataTypeFields;

#define CHROM (DidIType)0 // chrom is always the first field

#define TOPLEVEL "TOPLEVEL"
#define DATA_TYPE_FIELDS { \
/* num_fields        pos         test_regions    nonref        ochrom,     eol          toplevel names (including extend fields) - max 8 characters - 2 first chars must be unique within each data type (for dict_id_to_did_i_map) */ \
  {NUM_REF_FIELDS,   DID_I_NONE, DID_I_NONE,     DID_I_NONE,   DID_I_NONE, DID_I_NONE,  DID_I_NONE,       { "CONTIG", }, }, \
  {NUM_VCF_FIELDS,   VCF_POS,    VCF_POS,        DID_I_NONE,   VCF_oCHROM, VCF_EOL,     VCF_TOPLEVEL,     { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLES", "EOL", TOPLEVEL, "oCHROM", "oPOS", "oREFALT", "oXSTRAND", "COORDS", "o$TATUS", "C0PYPOS", "LIFT_REF", "CoPYSTAT", "ToPLUFT", } }, \
  {NUM_SAM_FIELDS,   SAM_POS,    SAM_POS,        SAM_NONREF,   DID_I_NONE, SAM_EOL,     SAM_TOPLEVEL,     { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "QUAL", "DOMQRUNS", "EOL", "BAM_BIN", TOPLEVEL, "TOP2BAM", "TOP2FQ", "E2:Z", "2NONREF", "N2ONREFX", "2GPOS", "S2TRAND", "U2:Z", "D2OMQRUN", "TAXID" } }, \
  {NUM_FASTQ_FIELDS, DID_I_NONE, DID_I_NONE,     FASTQ_NONREF, DID_I_NONE, FASTQ_E1L,   FASTQ_TOPLEVEL,   { "CONTIG", "DESC", "E1L", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "E2L", "QUAL", "DOMQRUNS", TOPLEVEL, "TAXID" } }, \
  {NUM_FASTA_FIELDS, DID_I_NONE, FASTA_LINEMETA, FASTA_NONREF, DID_I_NONE, FASTA_EOL,   FASTA_TOPLEVEL,   { "CONTIG", "LINEMETA", "EOL", "DESC", "COMMENT", "NONREF", "NONREF_X", TOPLEVEL, "TAXID" } }, \
  {NUM_GFF3_FIELDS,  GFF3_START, GFF3_START,     DID_I_NONE,   DID_I_NONE, GFF3_EOL,    GFF3_TOPLEVEL,    { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRS", "EOL", TOPLEVEL } }, \
  {NUM_ME23_FIELDS,  ME23_POS,   ME23_POS,       DID_I_NONE,   DID_I_NONE, ME23_EOL,    ME23_TOPLEVEL,    { "CHROM", "POS", "ID", "GENOTYPE", "EOL", TOPLEVEL, "TOP2VCF" } }, \
  {NUM_SAM_FIELDS,   SAM_POS,    SAM_POS,        SAM_NONREF,   DID_I_NONE, SAM_EOL,     SAM_TOP2BAM,      { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SQBITMAP", "NONREF", "NONREF_X", "GPOS", "STRAND", "QUAL", "DOMQRUNS", "EOL", "BAM_BIN", TOPLEVEL, "TOP2BAM", "TOP2FQ", "E2:Z", "2NONREF", "N2ONREFX", "2GPOS", "S2TRAND", "U2:Z", "D2OMQRUN", "TAXID" } }, \
  {NUM_VCF_FIELDS,   VCF_POS,    VCF_POS,        DID_I_NONE,   VCF_oCHROM, VCF_EOL,     VCF_TOPLEVEL,     { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLES", "EOL", TOPLEVEL, "oCHROM", "oPOS", "oREFALT", "oXSTRAND", "COORDS", "o$TATUS", "C0PYPOS", "LIFT_REF", "CoPYSTAT", "ToPLUFT", } }, \
  {NUM_GNRIC_FIELDS, DID_I_NONE, DID_I_NONE,     DID_I_NONE,   DID_I_NONE, DID_I_NONE,  GNRIC_TOPLEVEL,   { "DATA", TOPLEVEL } }, \
  {NUM_PHY_FIELDS,   DID_I_NONE, DID_I_NONE,     DID_I_NONE,   DID_I_NONE, PHY_EOL,     PHY_TOPLEVEL,     { "ID", "SEQ", "EOL", TOPLEVEL, "TOP2FA" } }, \
  {NUM_CHAIN_FIELDS, DID_I_NONE, CHAIN_NAMEDST,  DID_I_NONE,   DID_I_NONE, CHAIN_EOL,   CHAIN_TOPLEVEL,   { "NaMEDST", "SrANDDST", "StARTDST", "EnDDST", "SiZEDST", "NAMESRC", "SRANDSRC", "STARTSRC", "ENDSRC", "SIZESRC", "CHAIN", "SCORE", "ID", "SET", "SIZE", "GAPS", "EOL", TOPLEVEL, "SEP" } }, /* unique first 2 letters */ \
  {NUM_KRAKEN_FIELDS,DID_I_NONE, DID_I_NONE,     DID_I_NONE,   DID_I_NONE, KRAKEN_EOL,  KRAKEN_TOPLEVEL,  { "CU", "QNAME", "TAXID", "SEQLEN", "KMERS", "KMERTAX", "KMERLEN", "EOL", TOPLEVEL, "TOP2HASH" } }, \
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
    float factor;             // specifies the size of vb.txt_data allocated in PIZ vs the size recorded in SectionHeaderVbHeader.vb_data_size.
                              // It is fully allocated in advance and cannot be extended.
    TXTHEADER_TRANSLATOR ((*txtheader_translator));
    bool trans_containers;    // whether to invoke translators in containers
    bool is_src_dt;           // this translation dst is actually the source file dt that was compressed (of which we have the MD5)
    bool (*is_translation)(VBlockP vb); // if non-NULL - called to determine if this VB should be translated
} DtTranslation;

#define TRANSLATIONS { \
    /*                           non-bin-dt  binary dst_txt_dt toplevel        factor  txtheader_transl.   trans_con. src_dt is_trans. */ \
    /* SAM to BAM           */ { DT_SAM,     false, DT_BAM,    SAM_TOP2BAM,       2,   txtheader_sam2bam,   true,     false, NULL }, /* BAM is expected to be smaller, but in edge cases numbers "1\t" -> uint32 and QUAL="*"->0xff X seq_len */ \
    /* SAM to FASTQ         */ { DT_SAM,     false, DT_FASTQ,  SAM_TOP2FQ,        1,   txtheader_sam2fq,    true,     false, NULL }, /* sizes of SEQ, QUAL, QNAME the same */ \
    /* Binary SAM to SAM    */ { DT_SAM,     true,  DT_SAM,    SAM_TOPLEVEL,      3,   txtheader_bam2sam,   false,    false, NULL }, /* BAM stores sequences in 2x and numbers in 1-2.5x */ \
    /* Binary SAM to BAM    */ { DT_SAM,     true,  DT_BAM,    SAM_TOP2BAM,       2,   NULL,                true,     true,  NULL }, /* BAM is expected to be smaller, but in edge cases numbers "1\t" -> uint32 */ \
    /* Binary BAM to FASTQ  */ { DT_SAM,     true,  DT_FASTQ,  SAM_TOP2FQ,        1.5, txtheader_sam2fq,    true,     false, NULL }, /* sizes of QNAME, QUAL the same, SEQ x2 */ \
    /* 23andMe to VCF       */ { DT_ME23,    false, DT_VCF,    ME23_TOP2VCF,      4,   txtheader_me232vcf,  true,     false, NULL }, \
    /* FASTA to Phylip      */ { DT_FASTA,   false, DT_PHYLIP, FASTA_TOPLEVEL,    1.1, txtheader_fa2phy,    true,     false, NULL }, \
    /* Phylip to FASTA      */ { DT_PHYLIP,  false, DT_FASTA,  PHY_TOP2FASTA,     1.1, txtheader_phy2fa,    true,     false, NULL }, \
    /* VCF: genocat --luft  */ { DT_VCF,     false, DT_VCF,    VCF_TOPLUFT,       1.3, NULL,                true,     false, vcf_vb_is_luft }, /* Length of --luft is may vary slightly due to different CHROM, POS, AC, END lengths */ \
    /* genocat --taxid      */ { DT_KRAKEN,  false, DT_NONE,   KRAKEN_TOP2TAXID,  1,   NULL,                true,     false, kraken_is_translation }, \
}

extern const DtTranslation dt_get_translation (VBlockP vb);
extern DataType dt_get_txt_dt (DataType dt);
extern const char *dt_name (DataType data_type);

#endif
