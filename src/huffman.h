// ------------------------------------------------------------------
//   huffman.h
//   Copyright (C) 2023-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

// compression
extern HuffmanP huffman_initialize (void);
extern void huffman_destroy (HuffmanP *h);
extern int  huffman_chew_one_sample (HuffmanP h, bytes sample, uint32_t sample_len);
extern void huffman_produce_compressor (HuffmanP h);
extern bool huffman_is_produced (HuffmanP h);
extern void huffman_compress (HuffmanP h, bytes uncomp, uint32_t uncomp_len, uint8_t *comp, uint32_t *comp_len);
extern int huffman_decompress (HuffmanP h, bytes comp, uint8_t *uncomp, uint32_t uncomp_len);

static inline int huffman_comp_len_required_allocation (int uncomp_len) { return uncomp_len * 8; }
