// ------------------------------------------------------------------
//   base250.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef BASE250_INCLUDED
#define BASE250_INCLUDED

#include "genozip.h"

// values 0 to 249 are used as numerals in base-250. 
// The remaining 6 values below are control characters, and can only appear in numerals[0].
#define BASE250_EMPTY_SF   250 // subfield declared in FORMAT is empty, terminating : present
#define BASE250_MISSING_SF 251 // subfield declared in FORMAT is missing at end of cell, no :
#define BASE250_2_NUMERALS 253 // this number has 2 numerals, starting from numerals[1]
#define BASE250_3_NUMERALS 254 // this number has 3 numerals
#define BASE250_4_NUMERALS 255 // this number has 4 numerals
typedef struct {
    uint8_t numerals[5];       // number in base-250 (i.e. 0-249), digit[0] is least significant digit. If num_numerals=2,3 or 4 then numerals[1] is 250,251 or 252
    uint8_t num_numerals;      // legal values - 1,2,3,4
    uint8_t unused[2];         // padding to 64 bit
} Base250;

extern Base250 base250_encode (uint32_t n);
extern uint32_t base250_decode (const uint8_t *str);

#endif