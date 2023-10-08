// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern void zip_one_file (rom vcf_basename, bool is_last_user_txt_file);
extern void zip_compress_all_contexts_b250 (VBlockP vb);
extern void zip_init_vb (VBlockP vb);
extern bool zip_is_input_exhausted (void);