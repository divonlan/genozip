// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "digest.h"

extern char *progress_new_component (rom component_name, rom message, int test_mode);
extern void progress_update (char **prefix, uint64_t sofar, uint64_t total, bool done);
extern void progress_update_status (char **prefix, rom status);
extern void progress_finalize_component (rom status);
extern void progress_finalize_component_time (rom status, Digest md5);
extern void progress_finalize_component_time_ratio (rom me, float ratio, Digest md5);
extern void progress_finalize_component_time_ratio_better (rom me, float ratio, rom better_than, float ratio_than, Digest md5);
extern void progress_concatenated_md5 (rom me, Digest md5);
