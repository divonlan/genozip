// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

extern void zip_one_file (const char *vcf_basename, bool *is_last_file, bool z_closes_after_me);
extern void zip_compress_all_contexts_b250 (VBlockP vb);
