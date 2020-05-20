// ------------------------------------------------------------------
//   header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DATA_TYPES_INCLUDED
#define DATA_TYPES_INCLUDED

#include "genozip.h"
#include "md5.h"
#include "section_types.h"

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
#define NUM_DATATYPES 6
typedef enum { DT_VCF_V1=-2, DT_NONE=-1, // these values are used in the code logic, they are never written to the file
               DT_VCF=0, DT_SAM=1, 
               DT_FASTQ=2, DT_FASTA=3, DT_GFF3=4,
               DT_ME23=5 } DataType; // these values go into SectionHeaderGenozipHeader.data_type

typedef struct DataTypeProperties {
    // ZIP properties and functions
    const char *name;
    enum {NO_RA, RA} has_random_access;
    unsigned sizeof_vb, sizeof_zip_dataline;
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    char txt_header_1st_char;  // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)
    void (*seg_initialize)(VBlockP);
    const char *(*seg_data_line)(VBlockP, const char *field_start_line);
    void (*compress)(VBlockP);
    void (*update_header)(VBlockP, uint32_t vcf_first_line_i);

    // PIZ functions
    bool (*read_one_vb)(VBlockP, SectionListEntryP);
    void (*uncompress)(VBlockP);

    // VBlock functions
    void (*release_vb)(VBlockP);
    void (*destroy_vb)(VBlockP);
    void (*initialize_vb)(VBlockP);
    void (*cleanup_memory)(VBlockP);

    // misc properties and functions
    const char *show_sections_line_name; // the header displayed in --show-sections
    const char *stat_dict_types[3]; // the dictionary type displayed in --show-sections
} DataTypeProperties;

#define usz(type) ((unsigned)sizeof(type))
#define DATA_TYPE_PROPERTIES { \
    { "VCF",     RA,    usz(VBlockVCF),  usz(ZipDataLineVCF),  HDR_MUST, '#', seg_vcf_initialize,   seg_vcf_data_line,   zip_vcf_compress_one_vb,  zfile_vcf_update_compressed_vb_header, piz_vcf_read_one_vb,  piz_vcf_uncompress_vb,    vb_vcf_release_vb,  vb_vcf_destroy_vb,  NULL, vb_vcf_cleanup_memory, "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "SAM",     RA,    usz(VBlockSAM),  usz(ZipDataLineSAM),  HDR_OK,   '@', NULL,                 seg_sam_data_line,   NULL,                     zfile_update_compressed_vb_header,     NULL,                 piz_sam_reconstruct_vb,   vb_sam_release_vb,  vb_sam_destroy_vb,  vb_sam_initialize_vb, NULL,  "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",   NO_RA, usz(VBlockFAST), usz(ZipDataLineFAST), HDR_NONE, -1,  NULL,                 seg_fastq_data_line, NULL,                     zfile_update_compressed_vb_header,     piz_fast_read_one_vb, piz_fastq_reconstruct_vb, vb_fast_release_vb, NULL,               NULL,                 NULL,  "Entries",         { "FIELD", "ERROR!", "DESC"   } }, \
    { "FASTA",   NO_RA, usz(VBlockFAST), usz(ZipDataLineFAST), HDR_NONE, -1,  seg_fasta_initialize, seg_fasta_data_line, NULL,                     zfile_update_compressed_vb_header,     piz_fast_read_one_vb, piz_fasta_reconstruct_vb, vb_fast_release_vb, NULL,               NULL,                 NULL,  "Lines",           { "FIELD", "ERROR!", "DESC"   } }, \
    { "GVF",     RA,    usz(VBlockGFF3), usz(ZipDataLineGFF3), HDR_OK,   '#', seg_gff3_initialize,  seg_gff3_data_line,  NULL,                     zfile_update_compressed_vb_header,     piz_gff3_read_one_vb, piz_gff3_reconstruct_vb,  vb_gff3_release_vb, vb_gff3_destroy_vb, NULL,                 NULL,  "Sequences",       { "FIELD", "ATTRS",  "ERROR!" } }, \
    { "23ANDME", RA,    usz(VBlock),     0,                    HDR_OK,   '#', NULL,                 seg_me23_data_line,  NULL,                     zfile_update_compressed_vb_header,     NULL,                 piz_me23_reconstruct_vb,  NULL,               NULL,               NULL,                 NULL,  "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }  \
}
extern DataTypeProperties dt_props[NUM_DATATYPES];
#define DTP(prop)  (dt_props[(vb)->    data_type].prop)
#define DTPZ(prop) (dt_props[z_file->  data_type].prop)
#define DTPT(prop) (dt_props[txt_file->data_type].prop)

// VCF fields: CHROM up to the FORMAT field - excluding the samples. Note: we treat REF and ALT and the tab between them as a 
// single field as they are correlated so compress better together
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT, VCF_GT, NUM_VCF_FIELDS } VcfFields;
typedef enum { SAM_QNAME, SAM_FLAG, SAM_RNAME, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL, SAM_SEQ, SAM_QUAL, NUM_SAM_FIELDS } SamFields;
typedef enum { FAST_DESC,  FAST_LINEMETA,  FAST_SEQ,  FASTQ_QUAL,    NUM_FASTQ_FIELDS } FastqFields;
typedef enum { FASTA_DESC, FASTA_LINEMETA, FASTA_SEQ, FASTA_COMMENT, NUM_FASTA_FIELDS } FastaFields;
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS, NUM_GFF3_FIELDS } Gff3Fields;
typedef enum { ME23_CHROM, ME23_POS, ME23_ID, ME23_GENOTYPE, ME23_HAS13, NUM_ME23_FIELDS } Me23Fields;  

#define MAX_NUM_FIELDS_PER_DATA_TYPE NUM_SAM_FIELDS // CAREFUL! for now its the biggest, change if it changes
#if MAX_NUM_FIELDS_PER_DATA_TYPE + 2*MAX_SUBFIELDS > MAX_DICTS 
#error "MAX_NUM_FIELDS_PER_DATA_TYPE too large"
#endif

typedef struct DataTypeFields {
    unsigned num_fields;
    int chrom, pos, info; // the fields, or -1 if this data type doesn't have them
    char *names[MAX_NUM_FIELDS_PER_DATA_TYPE]; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
    uint64_t *local_by_lzma[MAX_NUM_FIELDS_PER_DATA_TYPE];
    uint64_t *dont_move_to_local[MAX_NUM_FIELDS_PER_DATA_TYPE]; // dict_ids that don't get moved from dict to local, even if everything is a singleton
} DataTypeFields;

#define DATA_TYPE_FIELDS { \
/* num_fields        chrom       pos         info        names (including extend fields)                                                                   local_by_lzma, dont_move_to_local                                                         */ \
  {NUM_VCF_FIELDS,   VCF_CHROM,  VCF_POS,    VCF_INFO,   { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GT" },                    { &dict_id_fields[VCF_POS], &dict_id_fields[VCF_ID], 0 }, { &dict_id_fields[VCF_CHROM], &dict_id_fields[VCF_POS], &dict_id_fields[VCF_ID], &dict_id_INFO_END, 0 } }, \
  {NUM_SAM_FIELDS,   SAM_RNAME,  SAM_POS,    -1,         { "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "PNEXT", "TLEN", "OPTIONAL", "SEQ", "QUAL" }, { &dict_id_fields[SAM_POS], &dict_id_fields[SAM_PNEXT], &dict_id_fields[SAM_SEQ], &dict_id_OPTION_BD, &dict_id_OPTION_BI, 0 }, { &dict_id_fields[SAM_RNAME], &dict_id_fields[SAM_POS], &dict_id_fields[SAM_PNEXT], &dict_id_OPTION_mc, 0}  }, \
  {NUM_FASTQ_FIELDS, 1,         -1,          -1,         { "DESC", "LINEMETA", "SEQ", "QUAL" },                                                            { &dict_id_fields[FAST_SEQ], 0 },  { &dict_id_fields[FAST_LINEMETA],  0 }  },   \
  {NUM_FASTA_FIELDS, 1,         -1,          -1,         { "DESC", "LINEMETA", "SEQ", "COMMENT" },                                                         { &dict_id_fields[FASTA_SEQ], 0 }, { &dict_id_fields[FASTA_LINEMETA], 0 }  }, \
  {NUM_GFF3_FIELDS,  GFF3_SEQID, GFF3_START, GFF3_ATTRS, { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRS" },               { &dict_id_fields[GFF3_START], &dict_id_ATTR_Dbxref, &dict_id_ENSTid, &dict_id_ATTR_Variant_seq, 0}, {&dict_id_fields[GFF3_SEQID], &dict_id_fields[GFF3_START], &dict_id_fields[GFF3_END], &dict_id_fields[GFF3_ATTRS], &dict_id_ATTR_Dbxref, 0}  }, \
  {NUM_ME23_FIELDS,  ME23_CHROM, ME23_POS,   -1,         { "CHROM", "POS", "ID", "GENOTYPE", "HAS13" },                                                    { &dict_id_fields[ME23_ID], 0}, {0} }, \
}
extern DataTypeFields dt_fields[NUM_DATATYPES];
#define DTF(prop)  (dt_fields[vb->      data_type].prop)
#define DTFZ(prop) (dt_fields[z_file->  data_type].prop)
#define DTFT(prop) (dt_fields[txt_file->data_type].prop)

// list of ctx who's local data is compressed via a callback function
#define LOCAL_COMP_CALLBACKS { \
  { DT_SAM,   &dict_id_OPTION_BD,          zip_sam_get_start_len_line_i_bd    }, \
  { DT_SAM,   &dict_id_OPTION_BI,          zip_sam_get_start_len_line_i_bi    }, \
  { DT_SAM,   &dict_id_fields[SAM_SEQ],    zip_sam_get_start_len_line_i_seq   }, \
  { DT_SAM,   &dict_id_fields[SAM_QUAL],   zip_sam_get_start_len_line_i_qual  }, \
  { DT_FASTQ, &dict_id_fields[FASTQ_QUAL], zip_fast_get_start_len_line_i_qual }, \
  { DT_FASTQ, &dict_id_fields[FAST_SEQ],   zip_fast_get_start_len_line_i_seq  }, \
  { DT_FASTA, &dict_id_fields[FASTA_SEQ],  zip_fast_get_start_len_line_i_seq  }, \
}

extern const char *dt_name (DataType data_type);

#endif