// ------------------------------------------------------------------
//   zriter.h
//   Copyright (C) 2023-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

extern void zriter_write (BufferP buf, BufferP section_list, int64_t offset_in_z_file, bool background);
extern void zriter_wait_for_bg_writing (void);
extern void zriter_flush (void);
