// ------------------------------------------------------------------
//   ploptimize.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimize.h"
#include "strings.h"

// optimize numbers in the range (-99.5,99.5) to 2 significant digits
bool optimize_float_2_sig_dig (rom snip, unsigned len, float cap_value_at /* 0 if no cap */,
                               char *optimized_snip, unsigned *optimized_snip_len)
{
    if (!IS_DIGIT(snip[0]) && snip[0] != '.' && snip[0] != '-') return false; // not a number

    bool negative = (snip[0] == '-');

    // temporarily nul-terminate string and get number
    SAFE_NUL (&snip[len]);
    float fp = atof (snip);
    SAFE_RESTORE;

    // cap value if requested
    if (cap_value_at && fp > cap_value_at) fp = cap_value_at;

    if (negative) fp = -fp; // fp is always positive

    if (fp >= 99.49999999) return false; // numbers must be in the range (-9.5,9.5) for this optimization (add epsilon to account for floating point rounding)

    char *writer = optimized_snip;

    // effecient outputing of two significant digits - a lot faster that sprintf
    #define NUM_EXPS 8
    #define MAX_NUM_LEN (1 /* hyphen */ + (NUM_EXPS-1) /* prefix */ + 2 /* digits */)

    static const float exps[NUM_EXPS]    = { 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001 };
    static const float mult_by[NUM_EXPS] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
    static rom prefix = "0.0000000000000000000";
    unsigned e=0; for (; e < NUM_EXPS; e++)
        if (fp >= exps[e]) {
            int twodigits = round (fp * mult_by[e]); // eg 4.31->43 ; 4.39->44 ; 0.0451->45
            if (twodigits == 100) { // cannot happen with e=0 because we restrict the integer to be up to 8 in the condition above
                e--;
                twodigits = 10;
            }
            if (negative) *(writer++) = '-';

            if (e >= 2) 
                writer = mempcpy (writer, prefix, e);

            *(writer++) = twodigits / 10 + '0';
            unsigned second_digit = twodigits % 10;
            if (e==1 && second_digit) *(writer++) = '.';
            if (!e || second_digit) *(writer++) = twodigits % 10 + '0'; // trailing 0: we write 2 (not 2.0) and 0.2 (not 0.20)
            break;
        }

    if (e == NUM_EXPS) *(writer++) = '0'; // rounding a very small positive or negative number to 0

    *optimized_snip_len = writer - optimized_snip;
    
    //printf ("snip:%.*s optimized:%.*s\n", len, snip, *optimized_snip_len, optimized_snip);
    return true;
}

bool optimize_vector_2_sig_dig (rom snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len /* in / out */)
{
    char *writer = optimized_snip;
    unsigned digit_i=0;
    for (unsigned i=0; i <= len; i++) { 

        if (snip[i] == ',' || i == len) { // end of number

            // optimize might actually increase the length in edge cases, e.g. -.1 -> -0.1, so we
            // make sure we have enough room for another number
            if ((writer - optimized_snip) + MAX_NUM_LEN > *optimized_snip_len) return false;

            unsigned one_number_len;
            bool ret = optimize_float_2_sig_dig (&snip[i-digit_i], digit_i, 0, writer, &one_number_len);
            if (!ret) return false;
            writer += one_number_len;

            if (i < len) *(writer++) = ',';
            digit_i=0;
        }
        else digit_i++;
    }

    *optimized_snip_len = writer - optimized_snip;
    return true;
}

// change the quality scores to be in a small number of bins, similar to Illumina: https://sapac.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf
// This is assuming standard Sanger format of Phred scores between 0 and 93 encoded in ASCII 33->126
// See here: https://pythonhosted.org/OBITools/fastq.html
// 0,1,2 - unchanged
// 2–9   -> 6
// 10–19 -> 15
// 20–24 -> 22
// 25–29 -> 27
// 30–34 -> 33
// 35–39 -> 37
// ≥ 40  -> 40
void optimize_phred_quality_string (char *str, unsigned len)
{
#   define P(x) ((x)+33)
    static uint8_t phred_mapper[256] = {
        // non Phread ASCII 0-32 - unchanged
        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32, 
        
        // Phread values 0-93 -> ASCII 33-126
        P(0),  P(1),  P(6),  P(6),  P(6),  // Phred 0-4: same as Illumina
        P(6),  P(6),  P(6),  P(6),  P(6),  // Phred 5-9: same as Illumina
        P(15), P(15), P(15), P(15), P(15), // Phred 10-14: same as Illumina
        P(15), P(15), P(15), P(15), P(15), // Phred 14-19: same as Illumina
        P(22), P(22), P(22), P(22), P(22), // Phred 20-24: same as Illumina
        P(27), P(27), P(27), P(27), P(27), // Phred 25-29: same as Illumina
        P(33), P(33), P(33), P(33), P(33), // Phred 30-34: same as Illumina
        P(37), P(37), P(37), P(37), P(37), // Phred 35-39: same as Illumina
        P(42), P(42), P(42), P(42), P(42), // Phred 40-44
        P(47), P(47), P(47), P(47), P(47), // Phred 45-49
        P(52), P(52), P(52), P(52), P(52), // Phred 50-54
        P(57), P(57), P(57), P(57), P(57), // Phred 59-59
        P(62), P(62), P(62), P(62), P(62), // Phred 60-64
        P(67), P(67), P(67), P(67), P(67), // Phred 65-69
        P(72), P(72), P(72), P(72), P(72), // Phred 70-74
        P(77), P(77), P(77), P(77), P(77), // Phred 74-79
        P(82), P(82), P(82), P(82), P(82), // Phred 80-84
        P(87), P(87), P(87), P(87), P(87), // Phred 85-89
        P(91), P(91), P(91),               // Phred 90-92
        P(93),                             // Phred 93: keep maximum value as is (common in Pac Bio)

        // non Phread ASCII 127-255 - unchanged
        127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,
        157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,
        187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,
        217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,
        247,248,249,250,251,252,253,254,255
    };

    // do the mapping
    for (; len; str++, len--) *str = (char)phred_mapper[(uint8_t)*str];
}

