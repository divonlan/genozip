// ------------------------------------------------------------------
//   ploptimize.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PLOPTIMIZE_INCLUDED
#define PLOPTIMIZE_INCLUDED

#define PL_MAX_SNIP_LEN 300

#include <stdbool.h>

extern bool pl_optimize (const char *snip, unsigned len, char *updated_snip, unsigned *updated_len);

#endif
