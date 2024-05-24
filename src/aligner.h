// ------------------------------------------------------------------
//   aligner.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#define declare_seq_contexts ContextP __attribute__((unused)) \
    bitmap_ctx    = CTX(SAM_SQBITMAP),                        \
    nonref_ctx    = CTX(SAM_NONREF),                          \
    gpos_ctx      = CTX(SAM_GPOS),                            \
    gpos_d_ctx    = CTX(SAM_GPOS_DELTA),                      \
    gpos_r2_ctx   = CTX(SAM_GPOS_R2),                         \
    strand_ctx    = CTX(SAM_STRAND),                          \
    strand_r2_ctx = CTX(SAM_STRAND_R2),                       \
    seqmis_ctx    = CTX(SAM_SEQMIS_A),                        \
    seqins_ctx    = CTX(SAM_SEQINS_A)


typedef enum { MAPPING_NO_MAPPING, MAPPING_ALIGNED, MAPPING_PERFECT } MappingType;

extern MappingType aligner_seg_seq (VBlockP vb, STRp(seq), bool is_pair_2, PosType64 pair_gpos, bool pair_is_forward);
extern void aligner_reconstruct_seq (VBlockP vb, uint32_t seq_len, bool is_pair_2, bool is_perfect_alignment, ReconType reconstruct, char *first_mismatch_base, uint32_t *first_mismatch_offset, uint32_t *num_mismatches);
