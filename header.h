// ------------------------------------------------------------------
//   header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HEADER_INCLUDED
#define HEADER_INCLUDED

#include "genozip.h"
#include "md5.h"

// data types genozip can compress
#define NUM_DATATYPES 2
typedef enum { DATA_TYPE_VCF_V1=-2, DATA_TYPE_NONE=-1, DATA_TYPE_VCF=0, DATA_TYPE_SAM=1 } DataType; // these values go into SectionHeaderGenozipHeader.data_type
#define DATATYPE_NAMES { "VCF", "SAM" } // index in array matches values in DataType
extern const unsigned datatype_last_field[NUM_DATATYPES];
const extern unsigned chrom_did_i_by_data_type[NUM_DATATYPES];

// VCF related global parameters - set before any thread is created, and never change
extern uint32_t global_vcf_num_samples, global_vcf_num_displayed_samples;

// VCF fields: CHROM up to the FORMAT field - excluding the samples. Note: we treat REF and ALT and the tab between them as a 
// single field as they are correlated so compress better together
#define NUM_VCF_FIELDS 8 
typedef enum { VCF_CHROM, VCF_POS, VCF_ID, VCF_REFALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT } VcfFields;
extern const char *vcf_field_names[NUM_VCF_FIELDS];

// Note on SAM fields: 
// - the QNAME field is broken into its :-separated components
// - the POS and PNEXT are stored in vb->pos_data and compressed with lzma
// - RNEXT is stored in the RNAME dictionary
// - OPTIONAL fields are stored as a template and then each subfield has its own dictionary.
//   the template looks like eg: "CT:Z:NM:i:" for two subfields CT and NM
#define NUM_SAM_FIELDS 9
typedef enum { SAM_QNAME, SAM_FLAG, SAM_RNAME, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_PNEXT, SAM_TLEN, SAM_OPTIONAL } SamFields;
extern const char *sam_field_names[NUM_SAM_FIELDS];

extern void header_initialize(void);
extern bool header_txt_to_genozip (uint32_t *vcf_line_i);
extern bool header_genozip_to_txt (Md5Hash *digest);

// v1 compatibility (VCF only)
extern bool v1_header_genozip_to_vcf (Md5Hash *digest);
extern bool v1_vcf_header_get_vcf_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,uint64_t *num_items_concat,
                                          Md5Hash  *md5_hash_concat, char *created, unsigned created_len);

#endif
