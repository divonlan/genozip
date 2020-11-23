// ------------------------------------------------------------------
//   md5.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef MD5_INCLUDED
#define MD5_INCLUDED

#include "genozip.h"

// Md5Hash must be packed as it appears in a Section in the Genozip file format (will only be meaningful on CPUs with more than 128 bit though...)
#pragma pack(1) 

typedef union { 
    uint8_t  bytes[16]; 
    uint32_t words[4];
    uint64_t ulls[2];
} Md5Hash;
#define MD5HASH_NONE (Md5Hash){ .ulls = { 0, 0 } }

#pragma pack()

typedef struct {
    uint32_t     lo, hi;
    uint32_t     a, b, c, d;
    union {
        uint8_t  bytes[64];
        uint32_t words[16];
    } buffer;
    bool initialized;
} Md5Context;
#define MD5CONTEXT_NONE (Md5Context){ .lo=0, .hi=0, .buffer.words = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, .initialized=false }

extern Md5Hash md5_do (const void *data, uint32_t len);
extern void md5_update (Md5Context *ctx, ConstBufferP buf);
extern Md5Hash md5_finalize (Md5Context *ctx);

typedef struct { char s[34]; } MD5Display;
extern MD5Display md5_display (Md5Hash digest);
extern Md5Hash md5_snapshot (const Md5Context *ctx);
extern void md5_display_ctx (const Md5Context *ctx); // for debugging

#define md5_is_equal(digest1,digest2) ((digest1).ulls[0] == (digest2).ulls[0] && (digest1).ulls[1] == (digest2).ulls[1])
#define md5_is_zero(digest) (!(digest).ulls[0] && !(digest).ulls[1])

#endif
