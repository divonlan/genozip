// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

// tradeoff: larger is better compression, but in some cases might be slower retrieval speed
#define SAMPLES_PER_BLOCK      1024  // default unless changed with --sblock
#define SAMPLES_PER_BLOCK_STR "1024" // used in help text -- needs to be identical to SAMPLES_PER_BLOCK

extern void zip_dispatcher (const char *vcf_basename, unsigned max_threads, bool is_last_file);

extern void zip_set_global_samples_per_block (const char *num_samples_str);

#endif
