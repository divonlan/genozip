// ------------------------------------------------------------------
//   biopsy.h
//   Copyright (C) 2021-2025 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

extern void biopsy_init (rom optarg);
extern void biopsy_take (VBlockP vb);
extern bool biopsy_is_done (void);
extern void biopsy_data_is_exhausted (void);
extern void biopsy_finalize (void);
