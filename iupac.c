// ------------------------------------------------------------------
//   iupac.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "flags.h"
#include "iupac.h"

void iupac_set (const char *optarg)
{
    bool neg = optarg[0] == '^';
    flag.iupac = neg ? IUP_NEGATIVE : IUP_POSITIVE;

    for (const char *c = &optarg[neg]; *c; c++)
        flag.iupac_mask[(int)*c] = true;
        
    // missing sequences (SEQ=*) (appear only in SAM/BAM) - always included in positive --iupac and excluded in negative
    flag.iupac_mask[(int)'*'] = true; 
}

void iupac_show (void)
{
    if (flag.iupac) {
        iprintf ("iupac=%d: ", flag.iupac);
        for (unsigned i=0; i < 256; i++)
            if (flag.iupac_mask[i]) iprintf ("%c", i);
        iprint0 ("\n");
    }
    else
        iprint0 ("iupac=false\n");
}

bool iupac_is_included (const char *seq, unsigned seq_len)
{    
    bool caught=false;
    for (unsigned i=0; i < seq_len; i++)
        if (!flag.iupac_mask[(int)seq[i]]) {
            caught=true;
            break;
        }

    if (flag.iupac == IUP_POSITIVE) return !caught;
    else                            return  caught;
}
