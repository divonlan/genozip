// ------------------------------------------------------------------
//   ref_iupacs.h
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

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

#define ref_iupacs_is_included(ref, vb, range, pos, vcf_base) \
    (((range) == (vb)->iupacs_last_range[(ref)==gref] && (pos) > (vb)->iupacs_last_pos[(ref)==gref] && (pos) < (vb)->iupacs_next_pos[(ref)==gref]) ? false /* quick negative */ \
     : ref_iupacs_is_included_do ((ref), (vb), (range), (pos), (vcf_base)))
extern bool ref_iupacs_is_included_do (Reference ref, VBlockP vb, const Range *range, PosType64 pos, char vcf_base);
extern char ref_iupacs_get (Reference ref, const Range *r, PosType64 pos, bool reverse, PosType64 *next_pos);
