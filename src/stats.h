// ------------------------------------------------------------------
//   stats.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "buffer.h"

extern void stats_generate (void); // ZIP
extern void stats_display (void);  // PIZ and ZIP
extern void stats_read_and_display (void);     // PIZ
extern void stats_add_txt_name (rom fn);

extern Buffer stats_programs;
