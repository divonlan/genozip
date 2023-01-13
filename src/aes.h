// ------------------------------------------------------------------
//   aes.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define AES_BLOCKLEN       16   // aes encrypts blocks of 128 bits
#define AES_KEYLEN         32   // 256 bit

extern void aes_initialize (VBlockP vb, bytes key);
extern void aes_xcrypt_buffer (VBlockP vb, uint8_t *data, uint32_t length);
extern char *aes_display_key (bytes key);
extern char *aes_display_data (bytes data, unsigned data_len);

