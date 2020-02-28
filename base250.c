// ------------------------------------------------------------------
//   base250.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "base250.h"
#include "move_to_front.h"

// Used by ZIP only: encodes in both 8 and 16 bits, as this gets stored in z_file and then used by 
// different VBs to construct b250 sections, which might use different encodings even for the same
// dictionary (e.g. if first vb still has < 250 words, it will choose 8 bits, but if the next vb
// adds words to exceed 250, it will select the 16 bit encoding)
Base250 base250_encode (uint32_t n)
{
    static const uint32_t MAX_B250 = 250UL*250UL*250UL*250UL-1;

    ASSERT (n <= MAX_B250, "Error: number too large for base250. n=%u, MAX_B250=%u", n, MAX_B250);

    // get numberals in base 250 (i.e. each numeral is 0 to 249) - least-signifcant-first order (little endian)
    Base250 result;
    result.num_numerals[B250_ENC_8] = 1;
    result.n = n;

    for (unsigned i=1; i <= 4; i++) {
        result.numerals[B250_ENC_8][i] = result.numerals[B250_ENC_16][i] = n % 250;
        n /= 250;
        if (n) 
            result.num_numerals[B250_ENC_8]++;
        else 
            break;
    }
    result.num_numerals[B250_ENC_16] = result.num_numerals[B250_ENC_8];

    // calculate the 8bit encoding
    if (result.num_numerals[B250_ENC_8] > 1) 
        result.numerals[B250_ENC_8][B250_ENC_8] = (BASE250_2_NUMERALS-2) + result.num_numerals[B250_ENC_8]++; // 253, 254 or 255 for 2,3 or 4
    else
        result.numerals[B250_ENC_8][B250_ENC_8] = result.numerals[B250_ENC_8][B250_ENC_16];


    // calculate the 16b encoding
    if (result.num_numerals[B250_ENC_16] >= 3) 
        result.numerals[B250_ENC_16][B250_ENC_8] = (BASE250_2_NUMERALS-2) + result.num_numerals[B250_ENC_16]++; // 253, 254 or 255 for 2,3 or 4

    else if (result.num_numerals[B250_ENC_16] == 2) {
        result.numerals[B250_ENC_16][B250_ENC_8] = result.numerals[B250_ENC_16][B250_ENC_16];
        result.numerals[B250_ENC_16][B250_ENC_16] = result.numerals[B250_ENC_16][2];
    }
    else { // == 1
        if (result.n==0)  // #1 most frequent snip
            result.numerals[B250_ENC_16][B250_ENC_8] = 253;
        
        else {
            result.num_numerals[B250_ENC_16] = 2;
            result.numerals[B250_ENC_16][B250_ENC_8] = result.numerals[B250_ENC_16][B250_ENC_16];
            result.numerals[B250_ENC_16][B250_ENC_16] = 0;
        }
    }
/*
    printf ("n=%u numerals=%u: ", result.n, result.num_numerals);
    for (unsigned i=0; i<result.num_numerals; i++) 
        printf ("%u ", result.numerals[i]);
    printf ("\n");   
*/
    return result;
}

uint32_t base250_decode (const uint8_t **str, Base250Encoding encoding)
{
    ASSERT ((encoding == B250_ENC_16) || 
            (encoding == B250_ENC_8)  ||
            (encoding == B250_ENC_NONE && (*str)[B250_ENC_8] == BASE250_MISSING_SF), // if line format has no non-GT subfields 
            "Error: invalid encoding=%d", encoding);

    uint32_t ret, bytes_consumed=1;

    if (encoding == B250_ENC_16 && (*str)[B250_ENC_8] == BASE250_MOST_FREQ) // note BASE250_MOST_FREQ is the same value as BASE250_2_NUMERALS
        ret = 0; // most frequent snip
    
    else switch ((*str)[B250_ENC_8]) {
        case BASE250_ONE_UP:     ret = WORD_INDEX_ONE_UP     ; break;
        case BASE250_EMPTY_SF:   ret = WORD_INDEX_EMPTY_SF   ; break;
        case BASE250_MISSING_SF: ret = WORD_INDEX_MISSING_SF ; break;
        case BASE250_2_NUMERALS: bytes_consumed=3; ret = (*str)[B250_ENC_16] + 250 * (*str)[2]; break;  // only happens in B250_ENC_8, because BASE250_MOST_FREQ (also 253) was handled earlier
        case BASE250_3_NUMERALS: bytes_consumed=4; ret = (*str)[B250_ENC_16] + 250 * (*str)[2] + 62500 * (*str)[3]; break;
        case BASE250_4_NUMERALS: bytes_consumed=5; ret = (*str)[B250_ENC_16] + 250 * (*str)[2] + 62500 * (*str)[3] + 15625000 * (*str)[4]; break;
        default: // 0 to 249
            if (encoding == B250_ENC_16) { // B250_ENC_16 - default is two numerals
                bytes_consumed = 2;
                ret = (*str)[B250_ENC_8] + 250 * (*str)[B250_ENC_16];  
            }
            else ret = (*str)[B250_ENC_8]; // B250_ENC_8 - default is a single numeral
    }

    *str += bytes_consumed;

    return ret;
}


unsigned base250_len (const uint8_t *data, Base250Encoding encoding)
{
    switch (encoding) {
        case B250_ENC_16:
            if (*data <= 249) return 2;
            else if (*data >= 250 && *data <= 253) return 1;
            else return *data - 250; // 4 or 5

        case B250_ENC_8:
            if (*data <= 252) return 1;
            else return *data - 250; // 3 or 4 or 5

        case B250_ENC_NONE: // happens when called for a line that has no non-GT subfields
            ASSERT0 (*data == BASE250_MISSING_SF, "Error: data has no encoding, but a byte other than BASE250_MISSING_SF was found");
            return 1;
    }
    
    ABORT ("Invalid encoding type: %d", encoding);
    return 0; // never reaches here - avoid compiler warning
}
 