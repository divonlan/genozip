// ------------------------------------------------------------------
//   base250.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "base250.h"
#include "move_to_front.h"

Base250 base250_encode (uint32_t n)
{
    static const uint32_t MAX_B250 = 250UL*250UL*250UL*250UL-1;

    ASSERT (n <= MAX_B250, "Error: number too large for base250. n=%u, MAX_B250=%u", n, MAX_B250);

    // least-signifcant-first order (little endian)
    Base250 result;
    result.num_numerals = 1;
    result.n = n;

    for (unsigned i=1; i <= 4; i++) {
        result.numerals[i] = n % 250;
        n /= 250;
        if (n) 
            result.num_numerals++;
        else 
            break;
    }
/* OLD ENCODING
    if (result.num_numerals > 1) 
        result.numerals[0] = (BASE250_2_NUMERALS-2) + result.num_numerals++; // 253, 254 or 255 for 2,3 or 4
    else
        result.numerals[0] = result.numerals[1];
*/
//// NEW ENCODING
    if (result.num_numerals >= 3) 
        result.numerals[0] = (BASE250_2_NUMERALS-2) + result.num_numerals++; // 253, 254 or 255 for 2,3 or 4
    else if (result.num_numerals == 2) {
        result.num_numerals = 2;
        result.numerals[0] = result.numerals[1];
        result.numerals[1] = result.numerals[2];
    }
    else { // == 1
        if (n==0) { // #1 most frequent snip
            result.numerals[0] = 253;
        }
        else {
            result.num_numerals = 2;
            result.numerals[0] = result.numerals[1];
            result.numerals[1] = 0;
        }
    }
////    

    return result;
}
/*
static inline uint32_t base250_decode_do (const uint8_t *str, unsigned *num_numerals)
{
    if (str[0] < 250) {
        *num_numerals = 1;
        return str[0]; // single numeral 
    }
    else if (str[0] == BASE250_EMPTY_SF) {
        *num_numerals = 1;
        return WORD_INDEX_EMPTY_SF;
    }

    else if (str[0] == BASE250_MISSING_SF) {
        *num_numerals = 1;
        return WORD_INDEX_MISSING_SF;
    }

    else if (str[0] == BASE250_ONE_UP) {
        *num_numerals = 1;
        return WORD_INDEX_ONE_UP;
    }

    uint32_t result = 0;
    uint32_t factor = 1;
    *num_numerals = str[0] - BASE250_2_NUMERALS + 2; // 2, 3 or 4
    for (unsigned i=1; i <= *num_numerals; i++) {
        result += str[i] * factor;
        factor *= 250;
    }

    return result;
}
*/


static inline uint32_t base250_decode_do (const uint8_t *str, unsigned *num_numerals)
{
    if (str[0] == 253) {
        *num_numerals = 1;
        return 0; // most frequent snip
    }
    if (str[0] < 250) {
        *num_numerals = 2;
        return str[0] + 250 * str[1]; // single numeral 
    }
    else if (str[0] == BASE250_EMPTY_SF) {
        *num_numerals = 1;
        return WORD_INDEX_EMPTY_SF;
    }

    else if (str[0] == BASE250_MISSING_SF) {
        *num_numerals = 1;
        return WORD_INDEX_MISSING_SF;
    }

    else if (str[0] == BASE250_ONE_UP) {
        *num_numerals = 1;
        return WORD_INDEX_ONE_UP;
    }

    uint32_t result = 0;
    uint32_t factor = 1;
    *num_numerals = str[0] - BASE250_2_NUMERALS + 2; // 2, 3 or 4
    for (unsigned i=1; i <= *num_numerals; i++) {
        result += str[i] * factor;
        factor *= 250;
    }

    return result;
}


// decode 250 and also advance pointer to next item
uint32_t base250_decode (const uint8_t **str_p)
{
    unsigned num_numerals;
    uint32_t result = base250_decode_do (*str_p, &num_numerals);

    *str_p += num_numerals + (num_numerals > 1);
    
    return result;
}


unsigned base250_len (const uint8_t *data)
{
        if (*data <= 249) return 2;
        else if (*data >= 250 && *data <= 253) return 1;
        else return *data - 250; // 4 or 5
}
