// ------------------------------------------------------------------
//   base250.c
//   Copyright (C) 2019-2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

Base250 base250_encode (uint32_t n)
{
    static const uint32_t MAX_B250 = 250UL*250UL*250UL*250UL-1;

    ASSERT (n <= MAX_B250, "Error: number too large for base250. n=%u, MAX_B250=%u", n, MAX_B250);

    // least-signifcant first order
    Base250 result;
    result.num_numerals = 1;
    
    for (unsigned i=1; i <= 4; i++) {
        result.numerals[i] = n % 250;
        n /= 250;
        if (n) 
            result.num_numerals++;
        else 
            break;
    }

    if (result.num_numerals > 1) 
        result.numerals[0] = (BASE250_2_NUMERALS-2) + result.num_numerals++; // 253, 254 or 255 for 2,3 or 4
    else
        result.numerals[0] = result.numerals[1];

    return result;
}

uint32_t base250_decode (const uint8_t *str)
{
    if (str[0] < BASE250_2_NUMERALS) return str[0]; // single numeral or control character

    uint32_t result = 0;
    uint32_t factor = 1;
    unsigned num_numerals = str[0] - BASE250_2_NUMERALS + 2; // 2, 3 or 4
    for (unsigned i=1; i <= num_numerals; i++) {
        result += str[i] * factor;
        factor *= 250;
    }

    return result;
}