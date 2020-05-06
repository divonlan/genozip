// ------------------------------------------------------------------
//   ploptimize.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimize.h"
#include "dict_id.h"
#include "header.h"
#include "strings.h"

// optimize numbers in the range (-99.5,99.5) to 2 significant digits
static inline bool optimize_float_2_sig_dig (const char *snip, unsigned len, double cap_value_at /* 0 if no cap */,
                                             char *optimized_snip, unsigned *optimized_snip_len)
{
    if (!IS_DIGIT(snip[0]) && snip[0] != '.' && snip[0] != '-') return false; // not a number

    bool negative = (snip[0] == '-');

    // temporarily null-terminate string and get number
    char save = snip[len];
    ((char*)snip)[len] = 0;
    double fp = atof (snip);
    ((char*)snip)[len] = save;

    // cap value if requested
    if (cap_value_at && fp > cap_value_at) fp = cap_value_at;

    if (negative) fp = -fp; // fp is always positive

    if (fp >= 99.49999999) return false; // numbers must be in the range (-9.5,9.5) for this optimization (add epsilon to account for floating point rounding)

    char *writer = optimized_snip;

    // effecient outputing of two significant digits - a lot faster that sprintf
    #define NUM_EXPS 8
    #define MAX_NUM_LEN (1 /* hyphen */ + (NUM_EXPS-1) /* prefix */ + 2 /* digits */)

    static const double exps[NUM_EXPS]    = { 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001 };
    static const double mult_by[NUM_EXPS] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
    static const char *prefix = "0.0000000000000000000";
    unsigned e=0; for (; e < NUM_EXPS; e++)
        if (fp >= exps[e]) {
            int twodigits = round (fp * mult_by[e]); // eg 4.31->43 ; 4.39->44 ; 0.0451->45
            if (twodigits == 100) { // cannot happen with e=0 because we restrict the integer to be up to 8 in the condition above
                e--;
                twodigits = 10;
            }
            if (negative) *(writer++) = '-';

            if (e >= 2) {
                memcpy (writer, prefix, e);
                writer += e;
            }

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

static inline bool optimize_vector_2_sig_dig (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (len > OPTIMIZE_MAX_SNIP_LEN) return false; // too long - we can't optimize - return unchanged

    char *writer = optimized_snip;
    unsigned digit_i=0;
    for (unsigned i=0; i <= len; i++) { 

        if (snip[i] == ',' || i == len) { // end of number

            // optimize might actually increase the length in edge cases, e.g. -.1 -> -0.1, so we
            // make sure we have enough room for another number
            if ((writer - optimized_snip) + MAX_NUM_LEN > OPTIMIZE_MAX_SNIP_LEN) return false;

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

static inline bool optimize_vcf_pl (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (len > OPTIMIZE_MAX_SNIP_LEN) return false; // too long - we can't optimize - return unchanged

    char *writer = optimized_snip;
    unsigned digit_i=0;
    for (unsigned i=0; i <= len; i++) { 
        if (snip[i] == ',' || i == len) { // end of number
            if (digit_i == 1) 
                *(writer++) = snip[i-1];
                
            else if (digit_i > 2 || (digit_i == 2 && snip[i-2] >= '6')) { // optimize - phred score of 60 or more (= 1 in 1 ^ 10^-10) is changed to 99
                *(writer++) = '6';
                *(writer++) = '0';
            }
            else if (digit_i == 2) { 
                *(writer++) = snip[i-2];
                *(writer++) = snip[i-1];
            }
            else return false; // digit_i==0

            if (i < len) *(writer++) = ',';
            digit_i=0;
        }
        else if (IS_DIGIT (snip[i]))
            digit_i++;
        
        else return false; // another character
    }
    
    *optimized_snip_len = writer - optimized_snip;
    return true;
}

// we separate to too almost identical functions optimize_vcf_format, optimize_vcf_info to gain a bit of performace -
// less conditions - as the caller knows if he is format or info
bool optimize_vcf_format (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (dict_id.num == dict_id_FORMAT_GL) return optimize_vector_2_sig_dig (snip, len, optimized_snip, optimized_snip_len);
    if (dict_id.num == dict_id_FORMAT_GP) return optimize_vector_2_sig_dig (snip, len, optimized_snip, optimized_snip_len);
    if (dict_id.num == dict_id_FORMAT_PL) return optimize_vcf_pl (snip, len, optimized_snip, optimized_snip_len);
    
    ABORT ("Error in optimize: unsupport dict %s", err_dict_id (dict_id));
    return 0; // never reaches here, avoid compiler warning
}

bool optimize_vcf_info (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (dict_id.num == dict_id_INFO_VQSLOD) 
        return optimize_float_2_sig_dig (snip, len, 0, optimized_snip, optimized_snip_len);
    
    ABORT ("Error in optimize: unsupport dict %s", err_dict_id (dict_id));
    return 0; // never reaches here, avoid compiler warning
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
    static unsigned char phred_mapper[256];
    static bool phread_mapper_initialized = false;

#   define PHRED(i) ((i)-33)

    // one time initialization
    if (!phread_mapper_initialized) {
        for (unsigned i=0; i < 256; i++) {
            if (i <= 34) phred_mapper[i]=i;  // ASCII 0->34 unchanged
            else if (PHRED(i) >= 2  && PHRED(i) <= 9 ) phred_mapper[i] = 33 + 6;
            else if (PHRED(i) >= 10 && PHRED(i) <= 19) phred_mapper[i] = 33 + 15;
            else if (PHRED(i) >= 20 && PHRED(i) <= 24) phred_mapper[i] = 33 + 22;
            else if (PHRED(i) >= 25 && PHRED(i) <= 29) phred_mapper[i] = 33 + 27;
            else if (PHRED(i) >= 30 && PHRED(i) <= 34) phred_mapper[i] = 33 + 33;
            else if (PHRED(i) >= 35 && PHRED(i) <= 39) phred_mapper[i] = 33 + 37;
            else if (PHRED(i) >= 40                  ) phred_mapper[i] = 33 + 40;
        }
        phread_mapper_initialized = true;
    }

    // do the mapping
    for (; len; str++, len--) *str = (char)phred_mapper[(uint8_t)*str];
}

