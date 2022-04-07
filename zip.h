// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

extern void zip_one_file (rom vcf_basename, bool is_last_user_txt_file);
extern void zip_compress_all_contexts_b250 (VBlockP vb);
extern void zip_init_vb (VBlockP vb);
