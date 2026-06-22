// ------------------------------------------------------------------
//   aligner.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

// note: pointers are restricted if vb is restricted
#define declare_seq_contexts ContextP           \
    UNUSED bitmap_ctx    = CTX(SAM_SQBITMAP),   \
    UNUSED nonref_ctx    = CTX(SAM_NONREF),     \
    UNUSED gpos_ctx      = CTX(SAM_GPOS),       \
    UNUSED gpos_Δ_ctx    = CTX(SAM_GPOS_DELTA), \
    UNUSED gpos_r2_ctx   = CTX(SAM_GPOS_R2),    \
    UNUSED gpos_gap_ctx  = CTX(SAM_GPOS_GAP),   \
    UNUSED junction_ctx  = CTX(SAM_JUNCTION),   \
    UNUSED strand_ctx    = CTX(SAM_STRAND),     \
    UNUSED strand_r2_ctx = CTX(SAM_STRAND_R2),  \
    UNUSED seqmis_ctx    = CTX(SAM_SEQMIS_A),   \
    UNUSED seqins_ctx    = CTX(SAM_SEQINS_A)


typedef enum { MAPPING_NO_MAPPING, MAPPING_ALIGNED, MAPPING_SPLICED, MAPPING_PERFECT, MAPPING_PERFECT_SPLICED } MappingType;

extern MappingType aligner_seg_seq (VBlockP vb, STRp(seq), bool am_i_R2, PosType64 gpos_R1, bool is_forward_R1);
extern void aligner_seg_gpos_and_fwd (VBlockP vb, uint32_t seq_len, PosType64 gpos, PosType64 gpos2, bool is_forward, uint32_t junction, bool am_i_R2, PosType64 gpos_R1, bool is_forward_R1, PosType64𐤐 G1);

extern void aligner_reconstruct_seq (VBlockP vb, uint32_t seq_len, bool am_i_R2, bool is_spliced_seg2, bool is_perfect_alignment, ReconType reconstruct, int mismatches_len, char *mismatch_base, uint32_t *mismatch_offset, uint32_t *num_mismatches);
extern void aligner_recon_get_gpos_and_fwd (VBlockP vb, bool am_i_R2, bool spliced_2nd_segment, PosType64𐤐 gpos, bool𐤐 is_forward);
