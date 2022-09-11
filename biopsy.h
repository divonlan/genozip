// ------------------------------------------------------------------
//   biopsy.h
//   Copyright (C) 2021-2022 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

extern void biopsy_init (rom optarg);
extern void biopsy_take (VBlockP vb);

