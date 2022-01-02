// ------------------------------------------------------------------
//   base250.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "base250.h"
#include "endianness.h"
#include "context.h"

// Used by ZIP only
Base250 base250_encode (WordIndex n) // number to encode
{
    // note: the actual maximum that the format can carry is 250*256*256*256-1 but we restrict it to MAX_WORD_INDEX
    ASSERT (n >= 0 && n <= MAX_WORD_INDEX, "n=%d is out of range 0-%u", n, MAX_WORD_INDEX);

    Base250 result;
    result.n = n;
    if (n <= 2) {
        result.encoded.bgen = 0;
        result.encoded.numerals[0] = BASE250_MOST_FREQ0 + n;
    }
    else result.encoded.bgen = BGEN32 (n);

    return result;
}

WordIndex base250_decode (const uint8_t **str, bool advance, const char *ctx_name)
{
    #define ADVANCE(n) if (advance) *str += n

    ASSERT (*str, "*str is NULL in ctx=%s", ctx_name);

    switch ((*str)[0]) {
        case BASE250_MOST_FREQ0: ADVANCE(1); return 0;
        case BASE250_MOST_FREQ1: ADVANCE(1); return 1;
        case BASE250_MOST_FREQ2: ADVANCE(1); return 2;
        case BASE250_ONE_UP:     ADVANCE(1); return WORD_INDEX_ONE_UP;
        case BASE250_EMPTY_SF:   ADVANCE(1); return WORD_INDEX_EMPTY;
        case BASE250_MISSING_SF: ADVANCE(1); return WORD_INDEX_MISSING;
        default : {
            WordIndex value = (*str)[0] * 16777216 + (*str)[1] * 65536 + (*str)[2] * 256 + (*str)[3]; // careful not to use BGEN32 as string might not be aligned to word boundary
            ADVANCE(4);
            return value;
        }
    }
}
