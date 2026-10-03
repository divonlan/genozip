// ------------------------------------------------------------------
//   biopsy.h
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

#define BIOPSY_Z_FILE_NAME "biopsy.genozip" // R1 data when --biopsy + --pair

extern void biopsy_init (rom optarg);
extern void biopsy_take (VBlockP vb);
extern bool biopsy_is_done (void);
extern void biopsy_data_is_exhausted (void);
extern void biopsy_finalize (void);
extern void biopsy_compress (void);

extern void biopsy_bytes_init (rom optarg);
extern noreturn void biopsy_bytes (rom filename);

