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

static inline bool optimize_pl (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (len > OPTIMIZE_MAX_SNIP_LEN) goto fail; // too long - we can't optimize - return unchanged

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
            else goto fail; // digit_i==0

            if (i < len) *(writer++) = ',';
            digit_i=0;
        }
        else if (snip[i] >= '0' && snip[i] <= '9')
            digit_i++;
        
        else goto fail;// another character
    }
    
    *optimized_snip_len = writer - optimized_snip;
    return true;

fail:
    *optimized_snip_len = len;
    return false;
}

static inline bool optimize_gl (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (len > OPTIMIZE_MAX_SNIP_LEN) goto fail; // too long - we can't optimize - return unchanged

    char *writer = optimized_snip;
    unsigned digit_i=0;
    for (unsigned i=0; i <= len; i++) { 

        if (snip[i] == ',' || i == len) { // end of number

            // temporarily null-terminate string and get number
            char save = snip[i];
            ((char*)snip)[i] = 0;
            double fp = atof (&snip[i-digit_i]);
            ((char*)snip)[i] = save;
            
            if (fp > 0 || fp <= -9) goto fail; // GL numbers must be in the range (-9,0]

            // effecient outputing of two significant digits - a lot faster that sprintf
            #define NUM_EXPS 7
            static const double exps[NUM_EXPS]    = { -1.0, -0.1, -0.01, -0.001, -0.0001, -0.00001, -0.000001 };
            static const double mult_by[NUM_EXPS] = { 10,   100,  1000,  10000,  100000,  1000000,  10000000  };
            static const char *prefix = "-0.0000000000000000000";
            unsigned e=0; for (; e < NUM_EXPS; e++)
                if (fp <= exps[e]) {
                    int twodigits = -round (fp * mult_by[e]); // eg -4.31->43 -4.39->44 -0.0451->45
                    if (twodigits == 100) { // cannot happen with e=0 because we restrict the integer to be up to 8 in the condition above
                        e--;
                        twodigits = 10;
                    }
                    memcpy (writer, prefix, e+1 + (e>=1));
                    writer += e+1 + (e>=1);
                    *(writer++) = twodigits / 10 + '0';
                    if (!e) *(writer++) = '.';
                    unsigned second_digit = twodigits % 10;
                    if (second_digit || !e) *(writer++) = twodigits % 10 + '0'; // trailing 0: we write 2.0 but 0.2 (not 0.20)
                    break;
                }

            if (e == NUM_EXPS) {
                memcpy (writer, "-0.0", 4);
                writer += 4;
            }

            if (i < len) *(writer++) = ',';
            digit_i=0;
        }
        else digit_i++;
    }

    *optimized_snip_len = writer - optimized_snip;
    return true;

fail:
    *optimized_snip_len = len;
    return false;
}

bool optimize (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (dict_id.num == dict_id_GL) return optimize_gl (snip, len, optimized_snip, optimized_snip_len);
    if (dict_id.num == dict_id_PL) return optimize_pl (snip, len, optimized_snip, optimized_snip_len);
    
    ABORT ("Error in optimize: unsupport dict %s", dict_id_printable (dict_id).id);
    return 0; // never reaches here, avoid compiler warning
}
