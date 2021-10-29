// ------------------------------------------------------------------
//   aes.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#define AES_BLOCKLEN       16   // aes encrypts blocks of 128 bits
#define AES_KEYLEN         32   // 256 bit

extern void aes_initialize (VBlockP vb, const uint8_t *key);
extern void aes_xcrypt_buffer (VBlockP vb, uint8_t *data, uint32_t length);
extern char *aes_display_key (const uint8_t* key);
extern char *aes_display_data (const uint8_t *data, unsigned data_len);

