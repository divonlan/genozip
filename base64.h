// ------------------------------------------------------------------
//   base64.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef BASE64_INCLUDED
#define BASE64_INCLUDED

#include "genozip.h"

extern unsigned base64_encode (const uint8_t *in, unsigned in_len, char *b64_str);
extern void base64_decode (const char *b64_str, unsigned len, uint8_t *out, unsigned *out_len /* in / out */);

#define base64_size(plain_size) ((unsigned)(((plain_size) + 2) / 3) * 4)
#define base64_sizeof(type_or_variable) base64_size(sizeof(type_or_variable))
#endif