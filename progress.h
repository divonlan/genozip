// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "digest.h"

extern char *progress_new_component (const char *component_name, const char *message, int test_mode);
extern void progress_update (char **prefix, uint64_t sofar, uint64_t total, bool done);
extern void progress_update_status (char **prefix, const char *status);
extern void progress_finalize_component (const char *status);
extern void progress_finalize_component_time (const char *status, Digest md5);
extern void progress_finalize_component_time_ratio (const char *me, float ratio, Digest md5);
extern void progress_finalize_component_time_ratio_better (const char *me, float ratio, const char *better_than, float ratio_than, Digest md5);
extern void progress_concatenated_md5 (const char *me, Digest md5);
