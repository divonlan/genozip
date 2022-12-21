// ------------------------------------------------------------------
//   md5.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"

extern void md5_initialize (Md5Context *ctx);
extern Digest md5_finalize (Md5Context *ctx);
extern Digest md5_do (const void *data, uint32_t len);
extern void md5_update (Md5Context *ctx, const void *data, uint32_t len);
extern void md5_display_ctx (const Md5Context *ctx); // for debugging
