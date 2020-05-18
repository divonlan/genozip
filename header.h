// ------------------------------------------------------------------
//   header.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HEADER_INCLUDED
#define HEADER_INCLUDED

#include "genozip.h"
#include "md5.h"
#include "section_types.h"

extern void header_initialize(void);
extern bool header_txt_to_genozip (uint32_t *vcf_line_i);
extern bool header_genozip_to_txt (Md5Hash *digest);

// v1 compatibility (VCF only)
extern bool v1_header_genozip_to_vcf (Md5Hash *digest);
extern bool v1_vcf_header_get_vcf_header (uint64_t *uncompressed_data_size, uint32_t *num_samples,uint64_t *num_items_concat,
                                          Md5Hash  *md5_hash_concat, char *created, unsigned created_len);


extern uint32_t last_txt_header_len; // for stats

// VCF related global parameters - set before any thread is created, and never change
extern uint32_t global_vcf_num_samples, global_vcf_num_displayed_samples;

#endif
