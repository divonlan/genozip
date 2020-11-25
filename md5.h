// ------------------------------------------------------------------
//   md5.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef MD5_INCLUDED
#define MD5_INCLUDED

#include "genozip.h"
#include "digest.h"

extern void md5_initialize (Md5Context *ctx);
extern Digest md5_finalize (Md5Context *ctx);
extern Digest md5_do (const void *data, uint32_t len);
extern void md5_update (Md5Context *ctx, const void *data, uint32_t len);
extern void md5_display_ctx (const Md5Context *ctx); // for debugging

#endif
