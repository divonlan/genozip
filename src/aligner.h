// ------------------------------------------------------------------
//   aligner.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef enum { MAPPING_NO_MAPPING, MAPPING_ALIGNED, MAPPING_PERFECT } MappingType;

extern MappingType aligner_seg_seq (VBlockP vb, ContextP bitmap_ctx, STRp(seq), bool no_bitmap_if_perfect, bool is_pair_2, PosType pair_gpos, bool pair_is_forward);
extern void aligner_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, uint32_t seq_len, bool is_pair_2, bool is_perfect_alignment, ReconType reconstruct);

extern Bits aligner_seq_to_bitmap (rom seq, uint64_t seq_len, uint64_t *bitmap_words, bool *seq_is_all_actg);
