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

// reads VCF header and writes its compressed form to the GENOZIP file. returns num_samples.
extern bool vcf_header_vcf_to_genozip (VariantBlockP vb, unsigned *line_i, Buffer **first_data_line);
extern bool vcf_header_genozip_to_vcf (VariantBlockP vb, Md5Hash *digest /* out */);
extern bool vcf_header_get_vcf_header (FileP z_file, SectionHeaderVCFHeader *vcf_header_header, bool *encrypted);

#endif
