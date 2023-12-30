// ------------------------------------------------------------------
//   refhash.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "reference.h"

#define HOOK 'G'
#define HOOK_REV 'C' // complement of HOOK

#define MAX_ALIGNER_GPOS ((PosType64)0xfffffffe)
#define NO_GPOS  ((PosType64)0xffffffff)

extern unsigned num_layers;
extern bool bits_per_hash_is_odd;  // true bits_per_hash is odd
extern uint32_t nukes_per_hash;    // = layer_bits[0] / 2
extern uint32_t layer_bitmask[64]; // 1s in the layer_bits[] LSbs
extern uint32_t **refhashs;
extern Buffer refhash_buf;

// make-reference stuff
extern void refhash_initialize_for_make (void);
extern void refhash_compress_refhash (void);
extern void refhash_calc_one_range (VBlockP vb, ConstRangeP r, ConstRangeP next_r);
extern Digest refhash_get_digest (void);

// stuff for loading and using refhash when ZIPping a fastq or sam/bam file
extern uint64_t refhash_get_refhash_size (void);
extern void refhash_load (Reference ref);
extern void refhash_load_standalone (void);
extern void refhash_destroy (void);
extern void refhash_set_digest (Digest digest);
extern void refhash_update_layers (int64_t delta_bytes);
extern ConstBitsP refhash_get_emoneg (void);
