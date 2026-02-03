// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern void zip_one_file (rom vcf_basename, bool is_last_user_txt_file);
extern void zip_init_vb (VBlockP vb);
extern bool zip_is_input_exhausted (void);
extern uint64_t zip_get_target_progress (void);

extern LocalGetLineCB *zip_get_local_data_callback (DataType dt, ContextP ctx);
extern void zip_set_no_stons_if_callback (VBlockP vb);
