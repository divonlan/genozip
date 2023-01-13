// ------------------------------------------------------------------
//   base64.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern unsigned base64_encode (bytes in, unsigned in_len, char *b64_str);
extern uint32_t base64_decode (rom b64_str, unsigned *b64_str_len /* in / out */, uint8_t *data);

#define base64_size(plain_size) ((unsigned)(((plain_size) + 2) / 3) * 4)
#define base64_sizeof(type_or_variable) base64_size(sizeof(type_or_variable))
