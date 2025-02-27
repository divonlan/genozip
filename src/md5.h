// ------------------------------------------------------------------
//   md5.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"

extern void md5_initialize (Md5StateP ctx);
extern Digest md5_finalize (Md5StateP ctx);
extern Digest md5_do (const void *data, uint32_t len);
extern void md5_update (Md5StateP ctx, const void *data, uint32_t len);
extern void md5_display_state (const Md5State *ctx); // for debugging
extern Digest md5_read (const char str[32]);