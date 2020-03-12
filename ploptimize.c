// ------------------------------------------------------------------
//   ploptimize.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <string.h>
#include "ploptimize.h"

bool pl_optimize (const char *snip, unsigned len, char *updated_snip, unsigned *updated_len)
{
    if (len > PL_MAX_SNIP_LEN) return false; // too long - we can't optimize - return unchanged

    char *writer = updated_snip;
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
    
    *updated_len = writer - updated_snip;
    return true;
}