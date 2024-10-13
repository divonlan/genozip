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
extern void huffman_start_chewing (Did did_i, STRp(master), uint8_t max_prefix_len);
extern void huffman_chew_one_sample (Did did_i, STRp(sample), bool skip_if_same_as_prev);
extern void huffman_produce_compressor (Did did_i, const bool mask[256]);
extern void huffman_get_master (Did did_i, pSTRp(master));
extern void huffman_compress (Did did_i, STRp (uncomp), qSTR8p(comp));
extern uint32_t huffman_compress_or_copy (Did did_i, STRp(uncomp), STR8c(comp));
extern uint32_t huffman_get_theoretical_max_comp_len (Did did_i, uint32_t uncomp_len);
extern uint32_t huffman_compress_len (Did did_i, STRp(uncomp));
extern int huffman_uncompress (Did did_i, bytes comp, STRc(uncomp));
extern int huffman_uncompress_or_copy (Did did_i, bytes comp, STRc(uncomp));
extern int RECONSTRUCT_huffman_or_copy (VBlockP vb, Did did_i, uint32_t uncomp_len, bytes comp);
extern int huffman_uncompress_len (Did did_i, bytes comp, uint32_t uncomp_len);
extern bool huffman_issame (Did did_i, bytes comp, uint32_t uncomp_len, STRp(uncomp_compare_to));
extern void huffman_compress_section (Did did_i);
extern void huffman_piz_read_all (void);

#define huffman_exists(did_i) (ZCTX(did_i)->huffman.len == 1)

