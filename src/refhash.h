// ------------------------------------------------------------------
//   refhash.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "digest.h"
#include "buf_struct.h"

extern Buffer refhash_buf;

// make-reference stuff
extern void refhash_make_initialize (void);
extern void refhash_make_refhash (void);
extern Digest refhash_get_digest (void);

// stuff for loading and using refhash when ZIPping a fastq or sam/bam file
extern void refhash_set_ref_file_info (Digest digest, uint8_t ref_bases_per_hash, uint8_t ref_bits_per_hash_out, uint8_t ref_gpos_bytes, MakeRefSize make_ref_size);
extern uint64_t refhash_get_refhash_size (void);
extern void refhash_load (void);
extern void refhash_load_standalone (void);
extern void refhash_destroy (void);
extern bool refhash_exists (void);

typedef enum { REFHASH_COUNT_INSTANCES, REFHASH_OCCUPY } RefhashCalcType;
extern void refhash_calc_one_range (VBlockP vb, RefhashCalcType calc_type);

extern int gpos_bytes; // 4 or 5 - since 15.0.81
