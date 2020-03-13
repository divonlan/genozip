// ------------------------------------------------------------------
//   ploptimize.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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
            double fp = atof (&snip[i-digit_i]);
            if (fp > 0) return false; // GL numbers must be <= 0

            #define NUM_EXPS 7
            static const double exps[NUM_EXPS] = { -1.0, -0.1, -0.01, -0.001, -0.0001, -0.00001, -0.000001 };
            unsigned e=0; for (; e < sizeof(exps)/sizeof(exps[0]); e++)
                if (fp <= exps[e]) {      // eg -1, -4.3212
                    sprintf (writer, "%.*f", e+1, fp);
                    writer += strlen (writer);
                    break;
                }
            if (e == NUM_EXPS) {
                strcpy (writer, "-0.0");
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

bool optimize (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len)
{
    if (dict_id.num == dict_id_GL) return optimize_gl (snip, len, optimized_snip, optimized_snip_len);
    if (dict_id.num == dict_id_PL) return optimize_pl (snip, len, optimized_snip, optimized_snip_len);
    
    ABORT ("Error in optimize: unsupport dict %s", dict_id_printable (dict_id).id);
    return 0; // never reaches here, avoid compiler warning
}
