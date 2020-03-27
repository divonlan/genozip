// ------------------------------------------------------------------
//   samples.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SAMPLES_INCLUDED
#define SAMPLES_INCLUDED

extern void samples_add (const char *samples_str);
extern void samples_digest_vcf_header (Buffer *vcf_header_buf);

extern char *vcf_samples_is_included;
#define samples_am_i_included(sample_i) (!flag_samples || ((bool)(vcf_samples_is_included[sample_i]))) // macro for speed - this is called in the critical loop of reconstructing samples

extern bool samples_is_sb_included (uint32_t num_samples_per_block, uint32_t sb_i);

#endif

