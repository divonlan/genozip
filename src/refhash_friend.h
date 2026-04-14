// ------------------------------------------------------------------
//   refhash_friend.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

// shared between refhash_make.c refhash_load.c and aligner.c

#include "genozip.h"
#include "refhash.h"

#define HOOK 'G'
#define encoded_HOOK 2 // see acgt_decoder
#define HOOK_REV 'C'   // complement of HOOK

#define NO_HASH_ENT32    0xffffffffULL
#define NO_HASH_ENT40    0xffffffffffULL
#define MAX_ALIGNER_GPOS (NO_HASH_ENT40 - 1)

extern bool refhash_is_flat; // true iff reference file was created with 15.0.81 or later
extern int bases_per_hash; 
extern int bits_per_hash_in; // bases_per_hash * 2
extern int bits_per_hash_out;
extern Digest refhash_digest;

// FLAT (new) refhash: since 15.0.81
extern int gpos_bytes;
extern int hash_shift;
extern MakeRefSize make_ref_size;

#define ref_hash_len (1ULL << bits_per_hash_out)  // number of hash enrtries. each entry is uint32_t or uint40_t depending on reference length

#define fibonacci_hash(kmer) (((kmer) * 11400714819323198485ULL) >> hash_shift)

// LAYERED (old) refhash : up to 15.0.80
#define NUM_LAYERS 4 
#define LAYER_START ((uint64_t[NUM_LAYERS]){ 0, 256 MB, 384 MB, 448 MB }) // start of layer within refhash_buf (units of uint32_t: layer lengths are: 256, 128, 64, 32 Mega entries)