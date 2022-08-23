// ------------------------------------------------------------------
//   refhash.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "reference.h"

#define HOOK 'G'
#define HOOK_REV 'C' // complement of HOOK

#define MAX_ALIGNER_GPOS ((PosType)0xfffffffe)
#define NO_GPOS  ((PosType)0xffffffff)

extern unsigned num_layers;
extern bool bits_per_hash_is_odd;  // true bits_per_hash is odd
extern uint32_t bits_per_hash;     // = layer_bits[0]
extern uint32_t nukes_per_hash;    // = layer_bits[0] / 2
extern uint32_t layer_bits[64];    // number of bits in each layer - layer_bits[0] is the base (widest) layer
extern uint32_t layer_size[64];    // size (in bytes) of each layer - layer_size[0] is the base (biggest) layer
extern uint32_t layer_bitmask[64]; // 1s in the layer_bits[] LSbs
extern uint32_t **refhashs;

extern void refhash_initialize (bool *dispatcher_invoked);
extern void refhash_destroy (bool destroy_only_if_not_mmap);
extern void refhash_verify_before_exit (void);

// make-reference stuff
extern void refhash_compress_refhash (void);
extern void refhash_calc_one_range (const Range *r, const Range *next_r);

// stuff for loading and using refhash when ZIPping a fastq or fasta file
extern void refhash_load_standalone (void);
extern void refhash_create_cache_join (bool free_mem);
extern void refhash_remove_cache (void);
extern bool refhash_has_refhash (void);

// globals
extern const char complement[256];
