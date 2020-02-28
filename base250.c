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

    uint8_t numerals[4];
    unsigned num_numerals = 1;
    for (unsigned i=0; i < 4; i++) {
        if (n==0) numerals[i] = 0;
        else {
            numerals[i] = n % 250;
            n /= 250;
            if (n) num_numerals++;
        }
    }

    result.num_numerals[B250_ENC_24] = result.num_numerals[B250_ENC_16] = result.num_numerals[B250_ENC_8];

    // calculate the 8bit encoding
    if (num_numerals > 1) {
        result.num_numerals[B250_ENC_8] = num_numerals + 1;
        result.numerals[B250_ENC_8][0] = (BASE250_2_NUMERALS-2) + num_numerals; // 253, 254 or 255 for 2,3 or 4
        result.numerals[B250_ENC_8][1] = numerals[0];
        result.numerals[B250_ENC_8][2] = numerals[1];
        result.numerals[B250_ENC_8][3] = numerals[2];
        result.numerals[B250_ENC_8][4] = numerals[3];
    }
    else {
        result.num_numerals[B250_ENC_8] = num_numerals;
        result.numerals[B250_ENC_8][0] = numerals[0];
    }

    // calculate the 16b encoding
    switch (num_numerals) {
        case 1:
            if (result.n==0) { // #1 most frequent snip
                result.numerals[B250_ENC_16][0]  = 253;
                result.num_numerals[B250_ENC_16] = 1;
                break;
            }
            // else for through
        case 2:
            result.num_numerals[B250_ENC_16] = 2;
            result.numerals[B250_ENC_16][0]  = numerals[0];
            result.numerals[B250_ENC_16][1]  = numerals[1];
            break;

        case 3:
        case 4:
            result.num_numerals[B250_ENC_16] = num_numerals + 1;
            result.numerals[B250_ENC_16][0]  = (BASE250_2_NUMERALS-2) + num_numerals;
            result.numerals[B250_ENC_16][1]  = numerals[0];
            result.numerals[B250_ENC_16][2]  = numerals[1];
            result.numerals[B250_ENC_16][3]  = numerals[2];
            result.numerals[B250_ENC_16][4]  = numerals[3];
            break;
    }

    // calculate the 24b encoding
    switch (num_numerals) {
        case 1:
            if (result.n==0) { // #1 most frequent snip
                result.numerals[B250_ENC_24][0]  = 253;
                result.num_numerals[B250_ENC_24] = 1;
                break;
            }
            // else for through
        case 2:
        case 3:
            result.num_numerals[B250_ENC_24] = 3;
            result.numerals[B250_ENC_24][0]  = numerals[0];
            result.numerals[B250_ENC_24][1]  = numerals[1];
            result.numerals[B250_ENC_24][2]  = numerals[2];
            break;

        case 4:
            result.num_numerals[B250_ENC_24] = 5;
            result.numerals[B250_ENC_24][0]  = BASE250_4_NUMERALS;
            result.numerals[B250_ENC_24][1]  = numerals[0];
            result.numerals[B250_ENC_24][2]  = numerals[1];
            result.numerals[B250_ENC_24][3]  = numerals[2];
            result.numerals[B250_ENC_24][4]  = numerals[3];
            break;
    }

    return result;
}

uint32_t base250_decode (const uint8_t **str, Base250Encoding encoding)
{
    ASSERT ((encoding >= B250_ENC_8 && encoding <= B250_ENC_24)  ||
            (encoding == B250_ENC_NONE && (*str)[0] == BASE250_MISSING_SF), // if line format has no non-GT subfields 
            "Error: invalid encoding=%d", encoding);

    uint32_t ret, bytes_consumed=1;

    if ((encoding == B250_ENC_16 || encoding == B250_ENC_24) && (*str)[0] == BASE250_MOST_FREQ) // note BASE250_MOST_FREQ is the same value as BASE250_2_NUMERALS
        ret = 0; // most frequent snip
    
    else switch ((*str)[0]) {
        case BASE250_ONE_UP:     ret = WORD_INDEX_ONE_UP     ; break;
        case BASE250_EMPTY_SF:   ret = WORD_INDEX_EMPTY_SF   ; break;
        case BASE250_MISSING_SF: ret = WORD_INDEX_MISSING_SF ; break;
        case BASE250_2_NUMERALS: bytes_consumed=3; ret = (*str)[1] + 250 * (*str)[2]; break;  // only happens in B250_ENC_8, because BASE250_MOST_FREQ (also 253) was handled earlier
        case BASE250_3_NUMERALS: bytes_consumed=4; ret = (*str)[1] + 250 * (*str)[2] + 62500 * (*str)[3]; break;
        case BASE250_4_NUMERALS: bytes_consumed=5; ret = (*str)[1] + 250 * (*str)[2] + 62500 * (*str)[3] + 15625000 * (*str)[4]; break;
        default: // 0 to 249
            if (encoding == B250_ENC_24) { // B250_ENC_24 - default is 3 numerals
                bytes_consumed = 3;
                ret = (*str)[0] + 250 * (*str)[1] + 62500 * (*str)[2];  
            }
            else if (encoding == B250_ENC_16) { // B250_ENC_16 - default is 2 numerals
                bytes_consumed = 2;
                ret = (*str)[0] + 250 * (*str)[1];  
            }
            else ret = (*str)[0]; // B250_ENC_8 - default is a single numeral
    }

    *str += bytes_consumed;

    return ret;
}


unsigned base250_len (const uint8_t *data, Base250Encoding encoding)
{
    switch (encoding) {
        case B250_ENC_24:
            if (*data <= 249) return 3;
            else if (*data >= 250 && *data <= 253) return 1;
            else return 5; // this is 255

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
 
const char *enc_name (Base250Encoding encoding)
{
    static const char *names[] = {"none", "8 bit", "16 bit", "24 bit"};
    ASSERT (encoding >= B250_ENC_NONE && encoding <= B250_ENC_24, "Error: invalid encoding %d", encoding);
    return names[encoding+1];
}

void base250_unit_test()
{
    static unsigned n[4] = {200, 20000, 100000, 100000000};

    for (Base250Encoding e=B250_ENC_8; e <= B250_ENC_24; e++)
        for (unsigned i=0; i<4; i++) {
            Base250 b = base250_encode (n[i]);
            uint8_t r[10];
            const uint8_t *p = r;
            memcpy (r, b.numerals[e], b.num_numerals[e]);
            unsigned d = base250_decode (&p, e);
            printf ("Decoded %u in %s: result=%u", n[i], enc_name(e), d);
            printf (" moved=%u bytes\n", (unsigned)(p-r));
            ASSERT0 (d == n[i], "Error: failed!");
        }
}