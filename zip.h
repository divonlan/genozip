// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

extern void zip_dispatcher (const char *vcf_basename, unsigned max_threads, bool is_last_file);

extern void zip_vcf_set_global_samples_per_block (const char *num_samples_str);

extern void zip_generate_b250_section (VBlockP vb, MtfContextP ctx);

extern void zip_vcf_compress_one_vb (VBlockP vb_);
extern void zip_sam_compress_one_vb (VBlockP vb_);

#endif
