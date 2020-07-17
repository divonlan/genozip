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

#define REFHASH_NOMATCH -10000000000000000LL
extern int64_t refhash_best_match (const char *seq, int64_t seq_len);

///// testing
int64_t refhash_get_final_match_len (const char *seq, int64_t gpos, int64_t fwd_len);

#endif
