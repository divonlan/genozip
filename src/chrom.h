// ------------------------------------------------------------------
//   chrom.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "file.h"

extern void chrom_2ref_load (void);

// seg
extern WordIndex chrom_seg_ex (VBlockP vb, Did did_i, STRp(chrom), PosType64 LN, int add_bytes, bool *is_new);
static inline WordIndex chrom_seg (VBlockP vb, STRp(chrom)) { return chrom_seg_ex (vb, CHROM, STRa(chrom), 0, chrom_len+1, NULL); }
static inline WordIndex chrom_seg_by_did_i (VBlockP vb, Did did_i, STRp(chrom), unsigned add_bytes) 
    { return chrom_seg_ex (vb, did_i, STRa(chrom), 0, add_bytes, NULL); }
extern bool chrom_seg_cb (VBlockP vb, ContextP ctx, STRp (chrom), uint32_t repeat);
extern WordIndex chrom_seg_no_b250 (VBlockP vb, STRp(chrom), bool *is_new);

// sorter
extern void chrom_index_by_name (Did chrom_did_i); // ZIP and PIZ
extern WordIndex chrom_get_by_name (STRp (chrom_name)); // ZIP and PIZ

// chrom2ref
extern void chrom_calculate_ref2chrom (uint64_t num_ref_contigs);
extern void chrom_2ref_compress (void);
extern WordIndex chrom_2ref_seg_get (ConstVBlockP vb, WordIndex chrom_index);
#define chrom_2ref_seg_is_needed(did_i) ((did_i) == CHROM && IS_REF_CHROM2REF)

extern void chrom_2ref_load (void);
static inline WordIndex chrom_2ref_piz_get (WordIndex chrom_index) 
    { return ZCTX(CHROM)->chrom2ref_map.len ? *B(WordIndex, ZCTX(CHROM)->chrom2ref_map, chrom_index) : chrom_index; }

extern void chrom_finalize (void);