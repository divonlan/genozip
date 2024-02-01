// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"

extern void progress_new_component (rom component_name, rom message, bool test_mode);
extern void progress_update (rom task, uint64_t sofar, uint64_t total, bool done);
extern void progress_finalize_component (rom status);
extern void progress_finalize_component_time (rom status, Digest md5);
extern void progress_finalize_component_time_ratio (rom me, float ratio, Digest md5);
extern void progress_finalize_component_time_ratio_better (rom me, float ratio, rom better_than, float ratio_than, Digest md5);
extern void progress_concatenated_md5 (rom me, Digest md5);
extern void progress_erase (void);
