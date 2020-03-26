// ------------------------------------------------------------------
//   visual_c_misc_funcs.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifdef _MSC_VER

#ifndef VISUAL_C_MISC_FUNCS_INCLUDED
#define VISUAL_C_MISC_FUNCS_INCLUDED

#include "../genozip.h"

extern void usleep(uint64_t usec);

extern double log2 (double n);

extern int my_round (double n);

#define PRIx64       "I64x"
#define PRId64       "I64d"
#define PRIu64       "I64u"

#endif

#endif