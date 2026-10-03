// ------------------------------------------------------------------
//   ref_iupacs.h
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

extern const char base2iupac[256];

// make-reference side
extern void ref_iupacs_compress (void);
extern void ref_iupacs_after_compute (VBlockP vb);

extern void ref_iupacs_add_do (VBlockP vb, uint64_t idx, char iupac);
static inline void ref_iupacs_add (VBlockP vb, uint64_t idx, char base)
{
    if (base2iupac[(uint8_t)base]) ref_iupacs_add_do (vb, idx, base2iupac[(uint8_t)base]);
}

#define IUPAC_IS_INCLUDED(ref_base,vcf_base) hxcgcb

// using with --chain side
extern void ref_iupacs_load (void);

#define ref_iupacs_is_included(vb, range, pos, vcf_base) \
    (((range) == (vb)->iupacs_last_range && (pos) > (vb)->iupacs_last_pos && (pos) < (vb)->iupacs_next_pos) ? false /* quick negative */ \
     : ref_iupacs_is_included_do ((vb), (range), (pos), (vcf_base)))
extern bool ref_iupacs_is_included_do (VBlockP vb, const Range *range, PosType64 pos, char vcf_base);
extern char ref_iupacs_get (const Range *r, PosType64 pos, bool reverse, PosType64 *next_pos);
