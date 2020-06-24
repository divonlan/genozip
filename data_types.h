// ------------------------------------------------------------------
//   header.h
//   Copyrigh8t (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DATA_TYPES_INCLUDED
#define DATA_TYPES_INCLUDED

#include "genozip.h"
#include "md5.h"
#include "section_types.h"

#include "vcf.h"
#include "sam.h"
#include "fastq.h"
#include "fasta.h"
#include "gff3.h"
#include "me23.h"

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
#define NUM_DATATYPES 6
typedef enum { DT_NONE=-1, // these values are used in the code logic, they are never written to the file
               DT_VCF=0, DT_SAM=1, DT_FASTQ=2, DT_FASTA=3, DT_GFF3=4, DT_ME23=5 } DataType; // these values go into SectionHeaderGenozipHeader.data_type

typedef void (*PizSpecialCtxHandler)(VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len);

typedef struct DataTypeProperties {
    // ZIP properties and functions
    const char *name;
    enum {NO_RA, RA} has_random_access;
    unsigned (*sizeof_vb)(void);
    unsigned (*sizeof_zip_dataline)(void);
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    char txt_header_1st_char;  // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)
    void (*seg_initialize)(VBlockP);
    const char *(*seg_txt_line)(VBlockP, const char *field_start_line, bool *has_13);
    void (*compress)(VBlockP);
    void (*update_header)(VBlockP, uint32_t vcf_first_line_i);

    // PIZ functions
    bool (*read_one_vb)(VBlockP, SectionListEntryP);
    void (*uncompress)(VBlockP);
    bool (*is_skip_secetion)(VBlockP, SectionType, DictId);
    unsigned num_special;
    PizSpecialCtxHandler special[10];

    // VBlock functions
    void (*release_vb)(VBlockP);
    void (*destroy_vb)(VBlockP);
    void (*cleanup_memory)(VBlockP);

    // misc properties and functions
    const char *show_sections_line_name; // the header displayed in --show-sections
    const char *stat_dict_types[3]; // the dictionary type displayed in --show-sections
} DataTypeProperties;

#define usz(type) ((unsigned)sizeof(type))
#define DATA_TYPE_PROPERTIES { \
/*    name       has_ra sizeof_vb     sizeof_zip_dataline  txt_headr 1st  seg_initialize        seg_txt_line        compress                  update_header                          read_one_vb           uncompress                is_skip_secetion           num_special        special        release_vb          destroy_vb           cleanup_memory          show_sections_line stat_dict_types                 */ \
    { "VCF",     RA,    vcf_vb_size,  vcf_vb_zip_dl_size,  HDR_MUST, '#', vcf_seg_initialize,   vcf_seg_txt_line,   vcf_zip_compress_one_vb,  vcf_zfile_update_compressed_vb_header, vcf_piz_read_one_vb,  vcf_piz_uncompress_vb,    vcf_piz_is_skip_section,   NUM_VCF_SPECIAL,   VCF_SPECIAL,   vcf_vb_release_vb,  vcf_vb_destroy_vb,   vcf_vb_cleanup_memory,  "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "SAM",     RA,    sam_vb_size,  sam_vb_zip_dl_size,  HDR_OK,   '@', sam_seg_initialize,   sam_seg_txt_line,   NULL,                     zfile_update_compressed_vb_header,     NULL,                 sam_piz_reconstruct_vb,   NULL,                      NUM_SAM_SPECIAL,   SAM_SPECIAL,   sam_vb_release_vb,  sam_vb_destroy_vb,   NULL,                   "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",   NO_RA, fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  fastq_seg_initialize, fastq_seg_txt_line, NULL,                     zfile_update_compressed_vb_header,     fast_piz_read_one_vb, fastq_piz_reconstruct_vb, fastq_piz_is_skip_section, 0,                 {},            fast_vb_release_vb, NULL,                NULL,                   "Entries",         { "FIELD", "DESC",   "ERROR!" } }, \
    { "FASTA",   RA,    fast_vb_size, fast_vb_zip_dl_size, HDR_NONE, -1,  fasta_seg_initialize, fasta_seg_txt_line, NULL,                     zfile_update_compressed_vb_header,     fast_piz_read_one_vb, fasta_piz_reconstruct_vb, fasta_piz_is_skip_section, NUM_FASTA_SPECIAL, FASTA_SPECIAL, fast_vb_release_vb, NULL,                NULL,                   "Lines",           { "FIELD", "DESC",   "ERROR!" } }, \
    { "GVF",     RA,    0,            0,                   HDR_OK,   '#', gff3_seg_initialize,  gff3_seg_txt_line,  NULL,                     zfile_update_compressed_vb_header,     NULL,                 gff3_piz_reconstruct_vb,  NULL,                      0,                 {},            NULL,               NULL,                NULL,                   "Sequences",       { "FIELD", "ATTRS",  "ITEMS"  } }, \
    { "23ANDME", RA,    0,            0,                   HDR_OK,   '#', me23_seg_initialize,  me23_seg_txt_line,  NULL,                     zfile_update_compressed_vb_header,     NULL,                 me23_piz_reconstruct_vb,  NULL,                      0,                 {},            NULL,               NULL,                NULL,                   "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }  \
}
extern DataTypeProperties dt_props[NUM_DATATYPES];
#define DTP(prop)  (dt_props[(vb)->    data_type].prop)
#define DTPZ(prop) (dt_props[z_file->  data_type].prop)
#define DTPT(prop) (dt_props[txt_file->data_type].prop)

// Fields - the CHROM field, if there is one, MUST be the first field (because of mtf_copy_reference_contig_to_chrom_ctx)
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_GT, VCF_EOL, NUM_VCF_FIELDS } VcfFields;
typedef enum { SAM_RNAME, SAM_QNAME, SAM_FLAG, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_RNEXT, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL, SAM_SEQ, SAM_QUAL, SAM_EOL, NUM_SAM_FIELDS } SamFields;
typedef enum { FASTQ_DESC, FASTQ_E1L, FASTQ_SEQ, FASTQ_E2L, FASTQ_PLUS, FASTQ_E3L, FASTQ_QUAL, FASTQ_E4L, NUM_FASTQ_FIELDS } FastqFields;
typedef enum { FASTA_CONTIG, FASTA_LINEMETA, FASTA_EOL, NUM_FASTA_FIELDS } FastaFields;
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS, GFF3_EOL, NUM_GFF3_FIELDS } Gff3Fields;
typedef enum { ME23_CHROM, ME23_POS, ME23_ID, ME23_GENOTYPE, ME23_EOL, NUM_ME23_FIELDS } Me23Fields;  

#define MAX_NUM_FIELDS_PER_DATA_TYPE NUM_SAM_FIELDS // CAREFUL! for now its the biggest, change if it changes

typedef struct DataTypeFields {
    unsigned num_fields;
    int chrom, pos, info, eol; // the fields, or -1 if this data type doesn't have them
    char *names[MAX_NUM_FIELDS_PER_DATA_TYPE]; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
} DataTypeFields;

#define DATA_TYPE_FIELDS { \
/* num_fields        chrom         pos         info        eol        names (including extend fields)                                                                   */ \
  {NUM_VCF_FIELDS,   VCF_CHROM,    VCF_POS,    VCF_INFO,   VCF_EOL,   { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GT", "EOL" },                    }, \
  {NUM_SAM_FIELDS,   SAM_RNAME,    SAM_POS,    -1,         SAM_EOL,   { "RNAME", "QNAME", "FLAG", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "OPTIONAL", "SEQ", "QUAL", "EOL" }, }, \
  {NUM_FASTQ_FIELDS, -1,           -1,         -1,         FASTQ_E1L, { "DESC", "E1L", "SEQ", "E2L", "PLUS", "E3L", "QUAL", "E4L" },                                                            }, \
  {NUM_FASTA_FIELDS, FASTA_CONTIG, -1,         -1,         FASTA_EOL, { "CONTIG", "LINEMETA", "EOL" },                                                         }, \
  {NUM_GFF3_FIELDS,  GFF3_SEQID,   GFF3_START, GFF3_ATTRS, GFF3_EOL,  { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRS", "EOL" },               }, \
  {NUM_ME23_FIELDS,  ME23_CHROM,   ME23_POS,   -1,         ME23_EOL,  { "CHROM", "POS", "ID", "GENOTYPE", "EOL" },                                                    }, \
}
extern DataTypeFields dt_fields[NUM_DATATYPES];
#define DTF(prop)  (dt_fields[vb->      data_type].prop)
#define DTFZ(prop) (dt_fields[z_file->  data_type].prop)
#define DTFT(prop) (dt_fields[txt_file->data_type].prop)

extern const char *dt_name (DataType data_type);

#endif