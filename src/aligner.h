// ------------------------------------------------------------------
//   aligner.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define declare_seq_contexts ContextP __attribute__((unused))   \
    bitmap_ctx = CTX(FASTQ_SQBITMAP),                           \
    nonref_ctx = CTX(FASTQ_NONREF),                             \
    gpos_ctx   = CTX(FASTQ_GPOS),                               \
    gpos_d_ctx = CTX(FASTQ_GPOS_DELTA),                         \
    strand_ctx = CTX(FASTQ_STRAND),                             \
    seqmis_ctx = CTX(FASTQ_SEQMIS_A)

typedef enum { MAPPING_NO_MAPPING, MAPPING_ALIGNED, MAPPING_PERFECT } MappingType;

extern MappingType aligner_seg_seq (VBlockP vb, STRp(seq), bool no_bitmap_if_perfect, bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward);
extern void aligner_reconstruct_seq (VBlockP vb, uint32_t seq_len, bool is_pair_2, bool is_perfect_alignment, ReconType reconstruct, char *first_mismatch_base, uint32_t *first_mismatch_offset, uint32_t *num_mismatches);
