// ------------------------------------------------------------------
//   base250.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef BASE250_INCLUDED
#define BASE250_INCLUDED

#include "genozip.h"

// the number is stored as a 4-byte Big Endian (MSB first), unless it is one of the below, in which case 
// only 1 byte is used
#define BASE250_EMPTY_SF   250 // empty string
#define BASE250_MISSING_SF 251 // container item missing, remove preceding separator
#define BASE250_ONE_UP     252 // value is one higher than previous value. used in B250_ENC_8
#define BASE250_MOST_FREQ0 253 // this translates to 0,1,2 representing the most frequent values (according to vb_i=1 sorting). used in B250_ENC_16
#define BASE250_MOST_FREQ1 254
#define BASE250_MOST_FREQ2 255

#define MAX_BASE250_NUMERALS 4
typedef struct {
    WordIndex n;                                 // the number being encoded
    union {
        uint8_t numerals[MAX_BASE250_NUMERALS];  // 4-byte big endian except if a 1 byte value 250-255
        int8_t shortcode;
        uint32_t bgen;                           // big endian
    } encoded;
} Base250;

extern Base250 base250_encode (WordIndex n);
extern WordIndex base250_decode (const uint8_t **str_p); // decodes and advances str_p
#define base250_len(data) ((*(data) < 250) ? 4 : 1)
#define base250_copy(dst, b250) memcpy (dst, (b250).encoded.numerals, base250_len((b250).encoded.numerals))

#endif