// ------------------------------------------------------------------
//   vcf_header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCF_HEADER_INCLUDED
#define VCF_HEADER_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "md5.h"
#include "sections.h"

// CHROM up to the FORMAT field - excluding the samples. Note: we treat REF and ALT and the tab between them as a 
// single field as they are correlated so compress better together
# define NUM_VCF_B250S 8 
typedef enum { CHROM, POS, ID, REFALT, QUAL, FILTER, INFO, FORMAT } VcfFields;
extern const char *vcf_field_names[];

// reads VCF header and writes its compressed form to the GENOZIP file. returns num_samples.
extern bool vcf_header_vcf_to_genozip (VariantBlockP vb, unsigned *line_i, Buffer **first_data_line);
extern bool vcf_header_genozip_to_vcf (VariantBlockP vb, Md5Hash *digest /* out */);

// v1 compatability
extern bool v1_vcf_header_genozip_to_vcf (VariantBlockP vb, Md5Hash *digest);

#endif
