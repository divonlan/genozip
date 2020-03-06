// ------------------------------------------------------------------
//   bzlib_mod.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//
// This is a hacky modification of bzlib.c that is separated out so that I can link with the standard bzlib

#include "genozip.h"
#include <bzlib.h>

// note that this struct is not aligned to 32/64 bit. we need to trust that the bzlib compiler
// and its options produce a similar alignment to ours...
typedef struct {
    void *a;
    char  b[5000];
    int32_t c;
    uint8_t d;
    bz_stream strm;
} bzFile;

const char *BZ2_errstr (int err)
{
    switch (err) {
        case BZ_OK:               return "BZ_OK";
        case BZ_RUN_OK:           return "BZ_RUN_OK";
        case BZ_FLUSH_OK:         return "BZ_FLUSH_OK";
        case BZ_FINISH_OK:        return "BZ_FINISH_OK";
        case BZ_STREAM_END:       return "BZ_STREAM_END";
        case BZ_SEQUENCE_ERROR:   return "BZ_SEQUENCE_ERROR";
        case BZ_PARAM_ERROR:      return "BZ_PARAM_ERROR";
        case BZ_MEM_ERROR:        return "BZ_MEM_ERROR";
        case BZ_DATA_ERROR:       return "BZ_DATA_ERROR";
        case BZ_DATA_ERROR_MAGIC: return "BZ_DATA_ERROR_MAGIC";
        case BZ_IO_ERROR:         return "BZ_IO_ERROR";
        case BZ_UNEXPECTED_EOF:   return "BZ_UNEXPECTED_EOF";
        case BZ_OUTBUFF_FULL:     return "BZ_OUTBUFF_FULL";
        case BZ_CONFIG_ERROR:     return "BZ_CONFIG_ERROR";
        default:                  return "Unknown BZ2 error";
    }
}

// this should go into bzlib.c
unsigned long long BZ2_bzoffset (BZFILE* b)
{
   return  (((unsigned long long)((bzFile*)b)->strm.total_in_hi32) << 32) |
            ((unsigned long long)((bzFile*)b)->strm.total_in_lo32);
}
