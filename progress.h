// ------------------------------------------------------------------
//   zip.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PROGRESS_INCLUDED
#define PROGRESS_INCLUDED

#include "genozip.h"

extern void progress_new_component (const char *component_name, bool txt_file_size_unknown, int test_mode);

extern void progress_update (uint64_t sofar, uint64_t total, bool done);

extern const char *progress_ellapsed_time (bool ever);

#endif