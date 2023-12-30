// ------------------------------------------------------------------
//   bases_filter.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "genozip.h"
#include "flags.h"
#include "bases_filter.h"
#include "flags.h"

uint8_t iupac_ascii_mask[256] = {}; // for --bases filter - '1' for every ASCII included in a positive or negative bases
uint8_t iupac_bam_mask[16]    = {}; // each entry corresponds to: =ACMGRSVTWYHKDBN (defined page 16: https://samtools.github.io/hts-specs/SAMv1.pdf)

void iupac_set (rom optarg)
{
    bool neg = optarg[0] == '^';
    flag.bases = neg ? IUP_NEGATIVE : IUP_POSITIVE;

    // ascii mask
    for (rom c = &optarg[neg]; *c; c++)
        iupac_ascii_mask[(int)*c] = true;
        
    // missing sequences (SEQ=*) (appear only in SAM/BAM) - always included in positive --bases and excluded in negative
    iupac_ascii_mask[(int)'*'] = true; 

    // bam mask
    for (rom c = &optarg[neg]; *c; c++) {
        // the characters "=ACMGRSVTWYHKDBN" are mapped to BAM 0->15, in this matrix we add 0x80 as a validity bit. All other characters are 0x00 - invalid
        static const uint8_t sam2bam_seq_map[256] = { ['=']=0x80, ['A']=0x81, ['C']=0x82, ['M']=0x83, ['G']=0x84, ['R']=0x85, ['S']=0x86, ['V']=0x87, 
                                                      ['T']=0x88, ['W']=0x89, ['Y']=0x8a, ['H']=0x8b, ['K']=0x8c, ['D']=0x8d, ['B']=0x8e, ['N']=0x8f };
        uint8_t bam_value = sam2bam_seq_map[(int)*c]; 
        ASSINP (bam_value & 0x80, "Invalid --bases argument \"%s\": each charcters should be a valid IUPAC character, one of \"=ACMGRSVTWYHKDBN\". See: https://www.bioinformatics.org/sms/iupac.html",
                optarg);
        iupac_bam_mask[(int)(bam_value & 0xf)] = 1;
    }
}

void iupac_show (void)
{
    if (flag.bases) {
        iprintf ("bases=%d: ", flag.bases);
        for (unsigned i=0; i < 256; i++)
            if (iupac_ascii_mask[i]) iprintf ("%c", i);
        iprint0 ("\n");
    }
    else
        iprint0 ("bases=false\n");
}

bool iupac_is_included_ascii (STRp(seq))
{   
    bool caught=false;
    for (unsigned i=0; i < seq_len; i++)
        if (!iupac_ascii_mask[(int)seq[i]]) {
            caught=true;
            break;
        }

    if (flag.bases == IUP_POSITIVE) return !caught;
    else                            return  caught;
}

bool iupac_is_included_bam (STRp(seq))
{    
    if (!seq_len) 
        return (flag.bases == IUP_POSITIVE) ? true : false;

    bool caught=false;
    for (unsigned i=0; i < seq_len; i++) {
        uint8_t value = (!(i%2)) ? ((seq[i/2] & 0xf0) >> 4)
                                 : seq[i/2] & 0x0f;
        if (!iupac_bam_mask[value]) {
            caught=true;
            break;
        }
    }

    if (flag.bases == IUP_POSITIVE) return !caught;
    else                            return  caught;
}
