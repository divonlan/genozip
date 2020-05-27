// ------------------------------------------------------------------
//   base250.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "base250.h"
#include "endianness.h"
#include "move_to_front.h"

// Used by ZIP only
Base250 base250_encode (uint32_t n) // number to encode
{
    static const uint32_t MAX_B250 = 250UL*256UL*256UL*256UL-1;

    ASSERT (n <= MAX_B250, "Error: number too large for base250. n=%u, MAX_B250=%u", n, MAX_B250);

    // get numberals in base 250 (i.e. each numeral is 0 to 249) - least-signifcant-first order (little endian)
    Base250 result;
    result.n = n;
    if (n <= 2) {
        result.encoded.bgen = 0;
        result.encoded.numerals[0] = BASE250_MOST_FREQ0 + n;
    }
    else result.encoded.bgen = BGEN32 (n);

    return result;
}

uint32_t base250_decode (const uint8_t **str)
{
    ASSERT0 (*str, "Error in base250_decode: *str is NULL");

    switch ((*str)[0]) {
        case BASE250_MOST_FREQ0: (*str)++; return 0;
        case BASE250_MOST_FREQ1: (*str)++; return 1;
        case BASE250_MOST_FREQ2: (*str)++; return 2;
        case BASE250_ONE_UP:     (*str)++; return WORD_INDEX_ONE_UP;
        case BASE250_EMPTY_SF:   (*str)++; return WORD_INDEX_EMPTY_SF;
        case BASE250_MISSING_SF: (*str)++; return WORD_INDEX_MISSING_SF;
        default : 
            *str += 4;
            return (*str)[-4] * 16777216 + (*str)[-3] * 65536 + (*str)[-2] * 256 + (*str)[-1]; // careful not to use BGEN32 as string might not be aligned to word boundary
    }
}

#define V1_BASE250
#include "vcf_v1.c"
