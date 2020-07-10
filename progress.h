// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PROGRESS_INCLUDED
#define PROGRESS_INCLUDED

#include "genozip.h"
#include "md5.h"

extern void progress_new_component (const char *component_name, const char *status, int test_mode);
extern void progress_update (uint64_t sofar, uint64_t total, bool done);
extern void progress_udpate_status (const char *status);
extern void progress_finalize_component (const char *status);
extern void progress_finalize_component_time (const char *status, Md5Hash md5);
extern void progress_finalize_component_time_ratio (const char *me, double ratio, Md5Hash md5);
extern void progress_finalize_component_time_ratio_better (const char *me, double ratio, const char *better_than, double ratio_than, Md5Hash md5);
extern void progress_concatenated_md5 (const char *me, Md5Hash md5);

#endif