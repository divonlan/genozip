// ------------------------------------------------------------------
//   vcf_header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCF_HEADER_INCLUDED
#define VCF_HEADER_INCLUDED

#include "genozip.h"
#include "md5.h"
#include "sections.h"

// CHROM up to the FORMAT field - excluding the samples. Note: we treat REF and ALT and the tab between them as a 
// single field as they are correlated so compress better together
# define NUM_VCF_B250S 8 
typedef enum { CHROM, POS, ID, REFALT, QUAL, FILTER, INFO, FORMAT } VcfFields;
extern const char *vcf_field_names[];

// reads VCF header and writes its compressed form to the GENOZIP file. returns num_samples.
extern bool vcf_header_vcf_to_genozip (uint32_t *vcf_line_i);
extern bool vcf_header_genozip_to_vcf (Md5Hash *digest);

extern void vcf_header_initialize(void);

// v1 compatability
extern bool v1_vcf_header_genozip_to_vcf (Md5Hash *digest);
extern bool v1_vcf_header_get_vcf_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,uint64_t *num_items_concat,
                                          Md5Hash  *md5_hash_concat, char *created, unsigned created_len);

#endif
