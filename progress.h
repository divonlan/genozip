// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PROGRESS_INCLUDED
#define PROGRESS_INCLUDED

#include "genozip.h"

extern void progress_new_component (const char *component_name, bool txt_file_size_unknown, int test_mode);

extern void progress_update (uint64_t sofar, uint64_t total, bool done);
extern void progress_udpate_status (const char *status);
extern void progress_finalize_component (const char *status);
extern void progress_finalize_component_time (const char *status);
extern void progress_finalize_component_time_ratio (const char *me, double ratio);
extern void progress_finalize_component_time_ratio_better (const char *me, double ratio, const char *better_than, double ratio_than);

#endif