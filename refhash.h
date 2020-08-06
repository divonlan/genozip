// ------------------------------------------------------------------
//   refhash.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REFHASH_INCLUDED
#define REFHASH_INCLUDED

#include "genozip.h"
#include "reference.h"

extern void refhash_initialize (void);
extern void refhash_free (void);

// make-reference stuff
extern void refhash_compress_refhash (void);
extern void refhash_calc_one_range (const Range *r, const Range *next_r);

// stuff for loading and using refhash when ZIPping a fastq or fasta file
extern void refhash_load(void);
extern void refhash_load_standalone (void);

extern bool refhash_best_match (VBlockP vb, const char *seq, const uint32_t seq_len, PosType *start_gpos, bool *is_forward, bool *is_all_ref);

// globals
extern const char complement[256];

#endif
