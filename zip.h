// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

extern void zip_dispatcher (const char *vcf_basename, unsigned max_threads, bool is_last_file);

extern void zip_set_global_samples_per_block (const char *num_samples_str);

#endif
