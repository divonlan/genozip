// ------------------------------------------------------------------
//   bai.h
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

typedef uint16_t BaiBinType;

extern bool bai_is_native_indexing_supported (DataType dt);
extern void bai_set_show_bai (rom optarg);
extern void bai_initialize (rom txt_header, uint32_t num_header_contigs, PosType64 len_of_longest_contig);
extern void bai_calculate_one_vb (VBlockP vb);
extern void bai_update_offsets_one_bb (VBlockP vb, uint32_t bb_i);
extern void bai_finalize_one_vb (VBlockP vb);
extern void bai_write (void);
extern void bai_show (rom filename);
extern BaiBinType bai_reg2bin (PosType32 first_pos, PosType32 last_pos);
