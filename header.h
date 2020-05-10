// ------------------------------------------------------------------
//   header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HEADER_INCLUDED
#define HEADER_INCLUDED

#include "genozip.h"
#include "md5.h"
#include "section_types.h"

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
#define NUM_DATATYPES 6
typedef enum { DT_VCF_V1=-2, DT_NONE=-1, // these values are used in the code logic, they are never written to the file
               DT_VCF=0, DT_SAM=1, 
               DT_FASTQ=2, DT_FASTA=3, DT_GFF3=4,
               DT_ME23=5 } DataType; // these values go into SectionHeaderGenozipHeader.data_type

typedef void (*ComputeFunc)(VBlockP);

typedef struct DataTypeProperties {
    // ZIP properties and functions
    const char *name;
    enum {NO_RA, RA} has_random_access;
    unsigned sizeof_vb, sizeof_zip_dataline;
    enum {HDR_NONE, HDR_OK, HDR_MUST} txt_header_required;
    char txt_header_1st_char;  // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)
    void (*seg_initialize)(VBlockP);
    const char *(*seg_data_line)(VBlockP, const char *field_start_line);
    ComputeFunc compress;
    void (*update_header)(VBlockP, uint32_t vcf_first_line_i);

    // PIZ functions
    void (*read_one_vb)(VBlockP);
    ComputeFunc uncompress;  

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
    { "VCF",     RA,    usz(VBlockVCF),  usz(ZipDataLineVCF),  HDR_MUST, '#', seg_vcf_initialize,   seg_vcf_data_line,   zip_vcf_compress_one_vb,  zfile_vcf_update_compressed_vb_header, piz_vcf_read_one_vb,  piz_vcf_uncompress_one_vb,  vb_vcf_release_vb,  vb_vcf_destroy_vb,  NULL, vb_vcf_cleanup_memory, "Variants",        { "FIELD", "INFO",   "FORMAT" } }, \
    { "SAM",     RA,    usz(VBlockSAM),  usz(ZipDataLineSAM),  HDR_OK,   '@', seg_sam_initialize,   seg_sam_data_line,   zip_sam_compress_one_vb,  zfile_update_compressed_vb_header,     piz_sam_read_one_vb,  piz_sam_uncompress_one_vb,  vb_sam_release_vb,  vb_sam_destroy_vb,  vb_sam_initialize_vb, NULL,  "Alignment lines", { "FIELD", "QNAME",  "OPTION" } }, \
    { "FASTQ",   NO_RA, usz(VBlockFAST), usz(ZipDataLineFAST), HDR_NONE, -1,  NULL,                 seg_fastq_data_line, zip_fast_compress_one_vb, zfile_update_compressed_vb_header,     piz_fast_read_one_vb, piz_fast_uncompress_one_vb, vb_fast_release_vb, vb_fast_destroy_vb, NULL,                 NULL,  "Entries",         { "FIELD", "ERROR!", "DESC"   } }, \
    { "FASTA",   NO_RA, usz(VBlockFAST), usz(ZipDataLineFAST), HDR_NONE, -1,  seg_fasta_initialize, seg_fasta_data_line, zip_fast_compress_one_vb, zfile_update_compressed_vb_header,     piz_fast_read_one_vb, piz_fast_uncompress_one_vb, vb_fast_release_vb, vb_fast_destroy_vb, NULL,                 NULL,  "Lines",           { "FIELD", "ERROR!", "DESC"   } }, \
    { "GVF",     RA,    usz(VBlockGFF3), usz(ZipDataLineGFF3), HDR_OK,   '#', seg_gff3_initialize,  seg_gff3_data_line,  zip_gff3_compress_one_vb, zfile_update_compressed_vb_header,     piz_gff3_read_one_vb, piz_gff3_uncompress_one_vb, vb_gff3_release_vb, vb_gff3_destroy_vb, NULL,                 NULL,  "Sequences",       { "FIELD", "ATTRS",  "ERROR!" } }, \
    { "23ANDME", RA,    usz(VBlockME23), 0,                    HDR_OK,   '#', seg_me23_initialize,  seg_me23_data_line,  zip_me23_compress_one_vb, zfile_update_compressed_vb_header,     piz_me23_read_one_vb, piz_me23_uncompress_one_vb, vb_me23_release_vb, vb_me23_destroy_vb, NULL,                 NULL,  "SNPs",            { "FIELD", "ERROR!", "ERROR!" } }  \
}
extern DataTypeProperties dt_props[NUM_DATATYPES];
#define DTP(prop)  (dt_props[(vb)->    data_type].prop)
#define DTPZ(prop) (dt_props[z_file->  data_type].prop)
#define DTPT(prop) (dt_props[txt_file->data_type].prop)

// VCF related global parameters - set before any thread is created, and never change
extern uint32_t global_vcf_num_samples, global_vcf_num_displayed_samples;

// VCF fields: CHROM up to the FORMAT field - excluding the samples. Note: we treat REF and ALT and the tab between them as a 
// single field as they are correlated so compress better together
#define NUM_VCF_FIELDS 8 
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT } VcfFields;

// Note on SAM fields: 
// - the QNAME field is broken into its :-separated components
// - the POS and PNEXT are stored in vb->pos_data and compressed with lzma
// - RNEXT is stored in the RNAME dictionary
// - OPTIONAL fields are stored as a template and then each subfield has its own dictionary.
//   the template looks like eg: "CT:Z:NM:i:" for two subfields CT and NM
#define NUM_SAM_FIELDS 9
typedef enum { SAM_QNAME, SAM_FLAG, SAM_RNAME, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL } SamFields;

// FASTQ/FASTA fields
#define NUM_FAST_FIELDS 2
typedef enum { FAST_DESC, FAST_LINEMETA } FastqFields;

#define NUM_GFF3_FIELDS 9 // https://m.ensembl.org/info/website/upload/gff3.html
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GFF3_ATTRS } Gff3Fields;

// 23ANDME fields
#define NUM_ME23_FIELDS 3
typedef enum { ME23_CHROM, ME23_POS, ME23_ID } Me23Fields; // same order as VCF

#define MAX_NUM_FIELDS_PER_DATA_TYPE 9 // maximum between NUM_*_FIELDS
#if MAX_NUM_FIELDS_PER_DATA_TYPE + 2*MAX_SUBFIELDS > MAX_DICTS 
#error "MAX_NUM_FIELDS_PER_DATA_TYPE too large"
#endif

typedef struct DataTypeFields {
    unsigned num_fields;
    int chrom; // used for random access data ; -1 if DATATYPE_HAS_RANDOM_ACCESS is false
    int info;
    SectionType first_dict_sec;
    SectionType info_dict_sec;
    SectionType info_sf_dict_sec;
    SectionType chrom_dict_sec; // used by --regions for subsetting
    char *names[MAX_NUM_FIELDS_PER_DATA_TYPE]; // these names go into the dictionary names on disk. to preserve backward compatibility, they should not be changed (names are not longer than 8=DICT_ID_LEN as the code assumes it)
} DataTypeFields;

#define DATA_TYPE_FIELDS { \
  {NUM_VCF_FIELDS,  VCF_CHROM,   VCF_INFO,  SEC_CHROM_DICT,      SEC_VCF_INFO_DICT,   SEC_VCF_INFO_SF_DICT,   SEC_CHROM_DICT,      { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT" } }, \
  {NUM_SAM_FIELDS,  SAM_RNAME,  -1,         SEC_SAM_QNAME_DICT, -1,                   -1,                     SEC_SAM_RNAME_DICT,  { "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "PNEXT", "TLEN", "OPTIONAL" } }, \
  {NUM_FAST_FIELDS, -1,         -1,         SEC_FAST_DESC_DICT, -1,                   -1,                     -1,                  { "DESC", "LINEMETA" } }, \
  {NUM_FAST_FIELDS, -1,         -1,         SEC_FAST_DESC_DICT, -1,                   -1,                     -1,                  { "DESC", "LINEMETA" } }, \
  {NUM_GFF3_FIELDS, GFF3_SEQID, GFF3_ATTRS, SEC_GFF3_SEQID_DICT, SEC_GFF3_ATTRS_DICT, SEC_GFF3_ATTRS_SF_DICT, SEC_GFF3_SEQID_DICT, { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTRS" } }, \
  {NUM_ME23_FIELDS, ME23_CHROM, -1,         SEC_CHROM_DICT,     -1,                   -1,                     SEC_CHROM_DICT,      { "CHROM", "POS", "ID" } }, \
}
extern DataTypeFields dt_fields[NUM_DATATYPES];
#define DTF(prop)  (dt_fields[vb->      data_type].prop)
#define DTFZ(prop) (dt_fields[z_file->  data_type].prop)
#define DTFT(prop) (dt_fields[txt_file->data_type].prop)

extern void header_initialize(void);
extern bool header_txt_to_genozip (uint32_t *vcf_line_i);
extern bool header_genozip_to_txt (Md5Hash *digest);

extern const char *dt_name (DataType data_type);

// v1 compatibility (VCF only)
extern bool v1_header_genozip_to_vcf (Md5Hash *digest);
extern bool v1_vcf_header_get_vcf_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,uint64_t *num_items_concat,
                                          Md5Hash  *md5_hash_concat, char *created, unsigned created_len);


#endif
