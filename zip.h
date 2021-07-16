// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef ZIP_INCLUDED
#define ZIP_INCLUDED

#include "genozip.h"

#define DYN_SET_MEM_TASK "dynamically_set_memory"

extern void zip_one_file (const char *vcf_basename, bool *is_last_file, bool z_closes_after_me);

#endif
