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
        else if (snip[i] >= '0' && snip[i] <= '9')
            digit_i++;
        
        else // another character
            return false;
    }
    
    *optimized_snip_len = writer - optimized_snip;
    return true;
}

static inline bool optimize_gl (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (len > OPTIMIZE_MAX_SNIP_LEN) return false; // too long - we can't optimize - return unchanged

    char *writer = optimized_snip;
    unsigned digit_i=0;
    for (unsigned i=0; i <= len; i++) { 

        if (snip[i] == ',' || i == len) { // end of number

            // temporarily null-terminate string and get number
            char save = snip[i];
            ((char*)snip)[i] = 0;
            double fp = atof (&snip[i-digit_i]);
            ((char*)snip)[i] = save;
            
            if (fp > 0) return false; // GL numbers must be <= 0

            #define NUM_EXPS 7
            static const double exps[NUM_EXPS]    = { -1.0, -0.1, -0.01, -0.001, -0.0001, -0.00001, -0.000001 };
            static const double mult_by[NUM_EXPS] = { 10,   100,  1000,  10000,  100000,  1000000,  10000000  };
            static const char *prefix = "-0.0000000000000000000";
            unsigned e=0; for (; e < NUM_EXPS; e++)
                if (fp <= exps[e]) {
                    int twodigits = -round (fp * mult_by[e]); // eg -4.31->43 -4.39->44 -0.0451->45
                    memcpy (writer, prefix, e+1 + (e>=1));
                    writer += e+1 + (e>=1);
                    *(writer++) = twodigits / 10 + '0';
                    if (!e) *(writer++) = '.';
                    *(writer++) = twodigits % 10 + '0';
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
}
/*
            // we its a simple number, do it the fast way without any external functions
            if (!e_fomrat) {
                char *d1=NULL; // first significant digit
                bool dec_point=false;
                bool found_integer
                for (unsigned j=i-digit_i; j < i; j++) {
                    // case: we've not reached a significant digit yet
                    if ((!d1 && snip[j] == '0')
                        *(writer++) = snip[j];

                    else if (snip[j] == '-') {
                        if (j) return false; // - can only appear as the first character
                        *(writer++) = snip[j];
                    } 
                    
                    else (snip[j] == '.') {
                        if (dec_point) return false;  // only one decimal point allowed
                        dec_point = true;
                        *(writer++) = snip[j];
                    }

                    else if (snip[j] < '0' || snip[j] > '9')
                        return false; // not a number - can't optimize

                    else if (!dec_point && snip[j] > )
                    // case: first significant digit
                    else if (!d1) { 
                        d1 = writer;
                        *(writer++) = snip[j];
                    }

                    // case second significant digit - no rounding needed
                    else if (j+1 == i || snip[j+1] <= '4') { 
                        *(writer++) = snip[j];
                        break; // we have 2 significant digits - we're done
                    }

                    // case second significant digit - need to round up, but only the second digit
                    else if (snip[j] < '9') {
                        *(writer++) = snip[j] + 1;
                        break; // we have 2 significant digits - we're done
                    }
                    
                    // case second significant digit - need to round up first digit too, but it is not 9
                    else if (*d1 < '9') {
                        *(writer++) = '0'; // second digit rounds fom 9 to 0
                        *d1++; // first digit increments
                        break; // we have 2 significant digits - we're done
                    }

                    // case - rounding, 9.9 -> 10.0
                    else if (snip[j-1] == '.')
                        *(writer-2) = '1';
                        *(writer-1) = '0';
                        
                        i                            
                        }
                        if ()
                    }
                }
*/


bool optimize (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (dict_id.num == dict_id_GL) return optimize_gl (snip, len, optimized_snip, optimized_snip_len);
    if (dict_id.num == dict_id_PL) return optimize_pl (snip, len, optimized_snip, optimized_snip_len);
    
    ABORT ("Error in optimize: unsupport dict %s", dict_id_printable (dict_id).id);
    return 0; // never reaches here, avoid compiler warning
}
