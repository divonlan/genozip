// ------------------------------------------------------------------
//   bzlib_mod.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//
// This is a hacky modification of bzlib.c that is separated out so that I can link with the standard bzlib

#include "genozip.h"
#include <bzlib.h>

// note that this struct is not aligned to 32/64 bit. we need to trust that the bzlib compiler
// is similar and used with similar options to ours...
typedef struct {
    void *a;
    char  b[5000];
    int32_t c;
    uint8_t d;
    bz_stream strm;
} bzFile;

// this should go into bzlib.h
BZ_EXTERN unsigned long long BZ_API(BZ2_bzoffset) (
      BZFILE* b
   );

// this should go into bzlib.c
unsigned long long BZ_API(BZ2_bzoffset) (BZFILE* b)
{
   return  (((unsigned long long)((bzFile*)b)->strm.total_in_hi32) << 32) |
            ((unsigned long long)((bzFile*)b)->strm.total_in_lo32);
}
