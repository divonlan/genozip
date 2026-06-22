// ------------------------------------------------------------------
//   base64.c
//
// derived from: https://github.com/launchdarkly/c-client-sdk/blob/master/base64.c
//
// Original unmodified license statement:
//  * Base64 encoding/decoding (RFC1341)
//  * Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
//  *
//  * This software may be distributed under the terms of the BSD license.
//

#include "base64.h"

static alignas(64) const char encode[64] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static alignas(64) const uint8_t decode[256] = {
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 0-15
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 16-31
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 62  , 0x80, 0x80, 0x80, 63  , // ASCII 32-47
	52  , 53  , 54  , 55  , 56  , 57  , 58  , 59  , 60  , 61  , 0x80, 0x80, 0x80, 0   , 0x80, 0x80, // ASCII 48-63 (outside of encode '=' = 0) 
    0x80, 0   , 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , // ASCII 64-79
	15  , 16  , 17  , 18  , 19  , 20  , 21  , 22  , 23  , 24  , 25  , 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 80-95
	0x80, 26  , 27  , 28  , 29  , 30  , 31  , 32  , 33  , 34  , 35  , 36  , 37  , 38  , 39  , 40  , // ASCII 96-111
	41  , 42  , 43  , 44  , 45  , 46  , 47  , 48  , 49  , 50  , 51  , 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 112-127
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 128-255
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80 };

#include "vblock.h"

// returns length of encoded (which is at most base64_sizeof)
// data must be allocated base64_sizeof bytes
unsigned base64_encode (STR8𐤐(data), char *restrict b64)
{
    ASSERTNOTNULL (data);

	bytes end = data + data_len;
	char *next = b64;
	while (end - data >= 3) {
		// get 3 bytes
        uint32_t bits = (data[0] << 16) | (data[1] << 8) | data[2];

		// translate 3 bytes to 4 base64 characaters (note: don't use a union of uint8_t[4] and uint32_t !!! that would force this variable to be in memory rather than a register)
		uint8_t out0 = encode[(bits >> 18) & 0x3F];
        uint8_t out1 = encode[(bits >> 12) & 0x3F];
        uint8_t out2 = encode[(bits >> 6 ) & 0x3F];
        uint8_t out3 = encode[(bits >> 0 ) & 0x3F];
        
        // Write all 4 bytes to memory at once
        PUT_UINT32 (next, out0 | (out1 << 8) | (out2 << 16) | (out3 << 24));
        next += 4;
        data += 3;		
	}

	if (end - data == 1) {
		uint8_t out0 = encode[data[0] >> 2];
        uint8_t out1 = encode[(data[0] & 0x03) << 4];
		uint8_t out2 = '=';
		uint8_t out3 = '=';

        PUT_UINT32 (next, out0 | (out1 << 8) | (out2 << 16) | (out3 << 24));
        next += 4;
	}

	else if (end - data == 2) {
		uint8_t out0 = encode[data[0] >> 2];
        uint8_t out1 = encode[((data[0] & 0x03) << 4) | (data[1] >> 4)];
		uint8_t out2 = encode[ (data[1] & 0x0f) << 2];
		uint8_t out3 = '=';

        PUT_UINT32 (next, out0 | (out1 << 8) | (out2 << 16) | (out3 << 24));
        next += 4;
	}	

    return next - b64;
}

// returns length of decoded data. b64_len is updated to the length of b64 string that was consumed.
// IMPORTANT: data allocated must be BASE_DECODE_OVERFLOW=4 bytes longer than expected data (overflow space)
uint32_t base64_decode (STR𐤐(b64), STR8c𐤐(out)) // out_len==1 if asking to read to end of snip (a non-b64 char like \0 or \t will terminate the b64 string) 
{
	if (out_len == -1) // read entire b64
		out_len = (b64_len / 4) * 3 
		        - (b64[b64_len - 2] == '=')
		        - (b64[b64_len - 1] == '=');

	uint8_t *after_out = &out[out_len];
	uint32_t save_overflow = GET_UINT32 (&out[out_len]); // we might overwrite these. actually only 3 bytes, but its faster to read/write 4 bytes
		
	for (; out < after_out; out += 3, b64 += 4) {
		uint8_t out0 = (decode[(int)b64[0]] << 2) | (decode[(int)b64[1]] >> 4);
        uint8_t out1 = (decode[(int)b64[1]] << 4) | (decode[(int)b64[2]] >> 2);
		uint8_t out2 = (decode[(int)b64[2]] << 6) |  decode[(int)b64[3]];

        PUT_UINT32 (out, out0 | (out1 << 8) | (out2 << 16)); // overwrites one byte because of writing 4-byte (instead of 3), and up to 2 additional bytes because we might need only 1 or 2 out of the 3 bytes
	}

	PUT_UINT32 (after_out, save_overflow); // restore

	return out_len;
}

// fast decoding of the common case of a dict_id 
DictId base64_decode_dict_id (rom b64)
{
	return (DictId){ .id = {
		[0] = (decode[(int)b64[0]] << 2) | (decode[(int)b64[1]] >> 4),
		[1] = (decode[(int)b64[1]] << 4) | (decode[(int)b64[2]] >> 2),
		[2] = (decode[(int)b64[2]] << 6) | (decode[(int)b64[3]] >> 0),
		[3] = (decode[(int)b64[4]] << 2) | (decode[(int)b64[5]] >> 4),
		[4] = (decode[(int)b64[5]] << 4) | (decode[(int)b64[6]] >> 2),
		[5] = (decode[(int)b64[6]] << 6) | (decode[(int)b64[7]] >> 0),
		[6] = (decode[(int)b64[8]] << 2) | (decode[(int)b64[9]] >> 4),
		[7] = (decode[(int)b64[9]] << 4) | (decode[(int)b64[10]]>> 2)
	}};
}

uint32_t base64_get_container_nitems (rom b64)
{
	uint8_t nitems_hi = ((decode[(int)b64[0]] << 2) | (decode[(int)b64[1]] >> 4)) & (CONTAINER_MAX_DICTS >> 8);
	uint8_t nitems_lo = (decode[(int)b64[5]] << 4) | (decode[(int)b64[6]] >> 2);

	return ((uint32_t)nitems_hi << 8) | (uint32_t)nitems_lo;
}
