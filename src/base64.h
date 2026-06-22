// ------------------------------------------------------------------
//   base64.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define BASE_DECODE_OVERFLOW 4
extern unsigned base64_encode (STR8𐤐(in), char *restrict b64_str);
extern uint32_t base64_decode (STR𐤐(b64_str), STR8c𐤐(out));
extern DictId base64_decode_dict_id (rom b64_str);
extern uint32_t base64_get_container_nitems (rom b64_str);

#define base64_size(plain_size) ((unsigned)(((plain_size) + 2) / 3) * 4)
#define base64_sizeof(type_or_variable) base64_size(sizeof(type_or_variable))
