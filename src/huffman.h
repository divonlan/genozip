// ------------------------------------------------------------------
//   huffman.h
//   Copyright (C) 2023-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "context.h"

// compression
extern void huffman_start_chewing (Did did_i, STRp(master), uint8_t max_prefix_len, int n_channels);
extern void huffman_chew_one_sample (Did did_i, STRp(sample), bool skip_if_same_as_prev);
extern uint32_t huffman_get_num_chewed (Did did_i);
typedef bool HuffmanMask[256];
extern void huffman_produce_compressor (Did did_i, const HuffmanMask *mask);
extern void huffman_get_master (Did did_i, pSTRp(master));
extern uint32_t huffman_compress (VBlockP vb, Did did_i, STRp (uncomp), qSTR8p(comp));
extern uint32_t huffman_get_theoretical_max_comp_len (Did did_i, uint32_t uncomp_len);
extern uint32_t huffman_compressed_len (Did did_i, STRp(uncomp));

extern int huffman_uncompress (Did did_i, bytes comp, STRc(uncomp));
extern int RECONSTRUCT_huffman (VBlockP vb, Did did_i, uint32_t uncomp_len, bytes comp);
extern int huffman_uncompress_len (Did did_i, bytes comp, uint32_t uncomp_len);
extern bool huffman_issame (Did did_i, bytes comp, uint32_t uncomp_len, STRp(uncomp_compare_to));
extern void huffman_compress_section (Did did_i);
extern void huffman_piz_read_all (void);
extern void huffman_piz_backcomp_produce_qual (STRp(qual));

typedef enum { HUFF_NOT_PRODUCED, HUFF_PRODUCED_BY_ZIP, HUFF_PRODUCED_BY_PIZ } HuffProducedStatus; // values of huffman.param
#define huffman_exists(did_i) (ZCTX(did_i)->huffman.param != HUFF_NOT_PRODUCED)

// NICO stuff
extern void nico_start_chewing (Did did_i);
extern void nico_chew_one_cigar (Did did_i, BamCigarOpP cigar, uint32_t n_cigar_op);
extern void nico_chew_one_textual_cigar (VBlockP vb, Did did_i, STRp(textual_cigar));
extern void nico_produce_compressor (Did did_i, bool is_deep_cigar);
extern uint32_t nico_max_comp_len (Did did_i, uint32_t n_cigar_op);
extern uint32_t nico_compress_cigar (VBlockP vb, Did did_i, BamCigarOpP cigar, uint32_t n_cigar_op, STR8c(comp));
extern uint32_t nico_compress_textual_cigar (VBlockP vb, Did did_i, STRp(textual_cigar), STR8c(comp));
extern uint32_t nico_compressed_len_textual_cigar (VBlockP vb, Did did_i, STRp(textual_cigar));
extern uint32_t nico_uncompress_cigar (VBlockP vb, Did did_i, bytes comp, BufferP cigar, rom buf_name);
extern uint32_t nico_uncompress_textual_cigar (Did did_i, bytes comp, BufferP textual_cigar, uint32_t textual_len, bool do_htos);
