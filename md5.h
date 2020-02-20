// ------------------------------------------------------------------
//   md5.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef MD5_INCLUDED
#define MD5_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#else
#include "compatability/visual_c_stdint.h"
#endif

// Md5Hash must be packed as it appears in a Section in the Genozip file format (will only be meaningful on CPUs with more than 128 bit though...)
#pragma pack(push, 1) 

typedef union { 
    uint8_t  bytes[16]; 
    uint32_t words[4];
    uint64_t ulls[2];
} Md5Hash;

#pragma pack(pop)

typedef struct {
    uint32_t     lo, hi;
    uint32_t     a, b, c, d;
    union {
        uint8_t  bytes[64];
        uint32_t words[16];
    } buffer;
} Md5Context;

extern void md5_do (const void *data, unsigned len, Md5Hash *digest);
extern void md5_update (Md5Context *ctx, const void *data, unsigned len, bool initialize);
extern void md5_finalize (Md5Context *ctx, Md5Hash *digest);
extern const char *md5_display (const Md5Hash *digest, bool prefix_space);

#define md5_is_equal(digest1,digest2) ((digest1).ulls[0] == (digest2).ulls[0] && (digest1).ulls[1] == (digest2).ulls[1])
#define md5_is_zero(digest) (!(digest).ulls[0] && !(digest).ulls[1])
#define md5_set_zero(digest) { (digest)->ulls[0] = (digest)->ulls[1] = 0; }

#endif
