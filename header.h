// ------------------------------------------------------------------
//   header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HEADER_INCLUDED
#define HEADER_INCLUDED

#include "genozip.h"
#include "md5.h"

// IMPORTANT: DATATYPES GO INTO THE FILE FORMAT - THEY CANNOT BE CHANGED
#define NUM_DATATYPES 6
typedef enum { DT_VCF_V1=-2, DT_NONE=-1, // these values are used in the code logic, they are never written to the file
               DT_VCF=0, DT_SAM=1, 
               DT_FASTQ=2, DT_FASTA=3, DT_GFF3=4,
               DT_ME23=5 } DataType; // these values go into SectionHeaderGenozipHeader.data_type
#define DATATYPE_NAMES { "VCF", "SAM", "FASTQ", "FASTA", "GVF" /* actually GFF3, but for now we can call it GVF */, \
                         "23ANDME" } // index in array matches values in DataType

#define DATATYPE_HAS_RANDOM_ACCESS { true, true, false, false, true, true }

typedef void (*ComputeFunc)(VBlockP);
#define COMPRESS_FUNC_BY_DT { zip_vcf_compress_one_vb, zip_sam_compress_one_vb,  \
                              zip_fast_compress_one_vb, zip_fast_compress_one_vb, zip_gff3_compress_one_vb, zip_me23_compress_one_vb }

#define UNCOMPRESS_FUNC_BY_DT { piz_vcf_uncompress_one_vb, piz_sam_uncompress_one_vb, \
                                piz_fast_uncompress_one_vb, piz_fast_uncompress_one_vb, \
                                piz_gff3_uncompress_one_vb, piz_me23_uncompress_one_vb }

typedef void (*UpdateHeaderFunc) (VBlockP vb, uint32_t vcf_first_line_i);
#define UPDATE_HEADER_FUNC_BY_DT { zfile_vcf_update_compressed_vb_header,     \
                                   zfile_update_compressed_vb_header, \
                                   zfile_update_compressed_vb_header, \
                                   zfile_update_compressed_vb_header, \
                                   zfile_update_compressed_vb_header, \
                                   zfile_update_compressed_vb_header  }         

typedef void (*IOFunc) (VBlockP vb);
#define READ_ONE_VB_FUNC_BY_DT { zfile_vcf_read_one_vb,  zfile_sam_read_one_vb,   \
                                 zfile_fast_read_one_vb, zfile_fast_read_one_vb, \
                                 zfile_gff3_read_one_vb, zfile_me23_read_one_vb  }

// by data type - related to the header of the txt file of each data type
#define TXT_HEADER_IS_ALLOWED      { true, true, false, false, true, true } // is it possible to have a header in this data_type
#define TXT_HEADER_IS_REQUIRED     { true, false, false, false, false, false } // should we error if the header is missing
#define TXT_HEADER_LINE_FIRST_CHAR { '#', '@', -1, -1, '#', '#' }; // first character in each line in the text file header (-1 if TXT_HEADER_IS_ALLOWED is false)

// by data type - the header displayed in --show-sections
#define STAT_SHOW_SECTIONS_LINE_NAME { "Variants", "Alignment lines", "Entries", "Lines", "Sequences", "SNPs" }

// by data type - the dictionary type displayed in --show-sections
#define STAT_DICT_TYPES { { "FIELD", "INFO",   "FORMAT" }, \
                          { "FIELD", "QNAME",  "OPTION" }, \
                          { "FIELD", "ERROR!", "DESC"   }, \
                          { "FIELD", "ERROR!", "DESC"   }, \
                          { "FIELD", "ATTRS",  "ERROR!" }, \
                          { "FIELD", "ERROR!", "ERROR!" } };

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

#define NUM_GFF3_FIELDS 10 // https://m.ensembl.org/info/website/upload/gff3.html
typedef enum { GFF3_SEQID, GFF3_SOURCE, GFF3_TYPE, GFF3_START, GFF3_END, GFF3_SCORE, GFF3_STRAND, GFF3_PHASE, GVF_SEQ, GFF3_ATTRS } Gff3Fields;

// 23ANDME fields
#define NUM_ME23_FIELDS 3
typedef enum { ME23_CHROM, ME23_POS, ME23_ID } Me23Fields; // same order as VCF

#define MAX_NUM_FIELDS_PER_DATA_TYPE 10 // maximum between NUM_*_FIELDS

#define DATATYPE_LAST_FIELD { NUM_VCF_FIELDS-1, NUM_SAM_FIELDS-1, NUM_FAST_FIELDS-1, NUM_FAST_FIELDS-1, NUM_GFF3_FIELDS-1, NUM_ME23_FIELDS-1 }
extern const unsigned datatype_last_field[NUM_DATATYPES];

#define CHROM_DID_I_BY_DT   { VCF_CHROM, SAM_RNAME, -1, -1, GFF3_SEQID, ME23_CHROM } // -1 if DATATYPE_HAS_RANDOM_ACCESS is false
extern const uint8_t chrom_did_i_by_dt[NUM_DATATYPES];  // used for random access data

#define INFO_DID_I_BY_DT { VCF_INFO, -1, -1, -1, GFF3_ATTRS, -1 }

#define FIRST_FIELD_DICT_SECTION { SEC_CHROM_DICT, SEC_SAM_QNAME_DICT, \
                                   SEC_FAST_DESC_DICT, SEC_FAST_DESC_DICT, SEC_GFF3_SEQID_DICT, SEC_CHROM_DICT }

// by data type - the chrom fields used by --regions for subsetting
#define CHROM_FIELD_DICT_SECTION { SEC_CHROM_DICT, SEC_SAM_RNAME_DICT, -1, -1, SEC_GFF3_SEQID_DICT, SEC_CHROM_DICT }
#define INFO_FIELD_DICT_SECTION  { SEC_VCF_INFO_DICT, -1, -1, -1, SEC_GFF3_ATTRS_DICT, -1 }
#define INFO_SF_DICT_SECTION     { SEC_VCF_INFO_SF_DICT, -1, -1, -1, SEC_GFF3_ATTRS_SF_DICT, -1 }

#define FIELD_NAMES /* max 8 chars per name */ \
    { { "CHROM", "POS", "ID", "REF+ALT", "QUAL", "FILTER", "INFO", "FORMAT" },\
      { "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "PNEXT", "TLEN", "OPTIONAL" },\
      { "DESC", "LINEMETA" },\
      { "DESC", "LINEMETA" },\
      { "SEQID", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "SEQ", "ATTRS" },\
      { "CHROM", "POS", "ID" }\
    };
extern const char *field_names[NUM_DATATYPES][MAX_NUM_FIELDS_PER_DATA_TYPE];

extern void header_initialize(void);
extern bool header_txt_to_genozip (uint32_t *vcf_line_i);
extern bool header_genozip_to_txt (Md5Hash *digest);

extern const char *dt_name (DataType data_type);

// v1 compatibility (VCF only)
extern bool v1_header_genozip_to_vcf (Md5Hash *digest);
extern bool v1_vcf_header_get_vcf_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,uint64_t *num_items_concat,
                                          Md5Hash  *md5_hash_concat, char *created, unsigned created_len);


#endif
