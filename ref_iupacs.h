// ------------------------------------------------------------------
//   ref_iupacs.h
//   Copyright (C) 2021-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REF_IUPACS_INCLUDED
#define REF_IUPACS_INCLUDED

#include "genozip.h"

// make-reference side
extern void ref_iupacs_compress (void);
extern void ref_iupacs_after_compute (VBlockP vb);

extern void ref_iupacs_add_do (VBlockP vb, uint64_t idx, char iupac);
static inline void ref_iupacs_add (VBlockP vb, uint64_t idx, char base)
{
    // true for multi-base iupac codes, except N: http://www.bioinformatics.org/sms/iupac.html
    static const char base2iupac[256] = { ['U']='U', ['R']='R', ['Y']='Y', ['S']='S', ['W']='W', ['K']='K', ['M']='M', ['B']='B', ['D']='D', ['H']='H', ['V']='V',
                                          ['u']='U', ['r']='R', ['y']='Y', ['s']='S', ['w']='W', ['k']='K', ['m']='M', ['b']='B', ['d']='D', ['h']='H', ['v']='V' };                                        
    if (base2iupac[(int)base]) ref_iupacs_add_do (vb, idx, base2iupac[(int)base]);
}

#define IUPAC_IS_INCLUDED(ref_base,vcf_base) hxcgcb

// using with --chain side
extern void ref_iupacs_load (Reference ref);

#define ref_iupacs_is_included(vb, luft_range, opos, vcf_base) \
    (((luft_range) == (vb)->iupacs_last_range && (opos) > (vb)->iupacs_last_opos && (opos) < (vb)->iupacs_next_opos) ? false /* quick negative */ \
     : ref_iupacs_is_included_do ((vb), (luft_range), (opos), (vcf_base)))
extern bool ref_iupacs_is_included_do (VBlockP vb, const Range *luft_range, PosType opos, char vcf_base);

#endif
