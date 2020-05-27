// ------------------------------------------------------------------
//   base64.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//
//   partially based on:

/*
 * Base64 encoding/decoding (RFC1341)
 * Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
 *
 * This software may be distributed under the terms of the BSD license.
 * See README for more details.
 */
#include "base64.h"

static const uint8_t base64_table[65] = 
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// returns length of encoded (which is at most base64_sizeof)
// out must be allocated base64_sizeof bytes
unsigned base64_encode (const uint8_t *in, unsigned in_len, char *b64_str)
{
    ASSERT0 (in, "Error in base64_encode: in is NULL");

	char *pos;
	const uint8_t *end;

	end = in + in_len;
	pos = b64_str;
	while (end - in >= 3) {
		*pos++ = base64_table[in[0] >> 2];
		*pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
		*pos++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
		*pos++ = base64_table[in[2] & 0x3f];
		in += 3;
	}

	if (end - in) {
		*pos++ = base64_table[in[0] >> 2];
		if (end - in == 1) {
			*pos++ = base64_table[(in[0] & 0x03) << 4];
			*pos++ = '=';
		} else {
			*pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
			*pos++ = base64_table[ (in[1] & 0x0f) << 2];
		}
		*pos++ = '=';
	}

    return pos - b64_str;
}


void base64_decode (const char *b64_str, unsigned len, uint8_t *out, unsigned *out_len /* in/out */)
{
	uint8_t dtable[256], block[4];

	memset (dtable, 0x80, 256);
	for (unsigned i=0; i < 64; i++)
		dtable[base64_table[i]] = (uint8_t)i;
	dtable['='] = 0;

	unsigned count=0, pad=0;
	for (unsigned i=0; i < len; i++)
		if (dtable[(unsigned)b64_str[i]] != 0x80)
			count++;

	ASSERT (count && !(count % 4), "Error in base64_decode: bad base64 - expecting it to be a string with length divisable by 4 but its length is %u: %.*s",
            len, len, b64_str);

	uint8_t *pos = out;
	count = 0;
	for (unsigned i=0; i < len; i++) {
		uint8_t tmp = dtable[(unsigned)b64_str[i]];
		if (tmp == 0x80) continue;

		if (b64_str[i] == '=') pad++;
		block[count] = tmp;
		count++;
		if (count == 4) {
			             *pos++ = (block[0] << 2) | (block[1] >> 4);
			if (!pad)    *pos++ = (block[1] << 4) | (block[2] >> 2);
			if (pad < 2) *pos++ = (block[2] << 6) |  block[3];
			count = 0;
		}

		ASSERT (*out_len >= (unsigned)(pos - out), "Error in base64_decode: 'out' is too small for decoding. out_len=%u b64_str=%.*s", *out_len, len, b64_str);
	}

	*out_len = pos - out;
}