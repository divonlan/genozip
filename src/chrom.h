// ------------------------------------------------------------------
//   chrom.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"
#include "reference.h"
#include "file.h"

extern void chrom_2ref_load (Reference ref);

// seg
extern WordIndex chrom_seg_ex (VBlockP vb, Did did_i, STRp(chrom), PosType64 LN, bool *is_alt_out, int add_bytes, bool recon_changes_if_match, bool *is_new);
static inline WordIndex chrom_seg (VBlockP vb, STRp(chrom)) { return chrom_seg_ex (vb, CHROM, STRa(chrom), 0, 0, chrom_len+1, true, NULL); }
static inline WordIndex chrom_seg_by_did_i (VBlockP vb, Did did_i, STRp(chrom), unsigned add_bytes) 
    { return chrom_seg_ex (vb, did_i, STRa(chrom), 0, 0, add_bytes, true, NULL); }
extern bool chrom_seg_cb (VBlockP vb, ContextP ctx, STRp (chrom), uint32_t repeat);
extern WordIndex chrom_seg_no_b250 (VBlockP vb, STRp(chrom), bool *is_new);

// sorter
extern void chrom_index_by_name (Did chrom_did_i); // ZIP and PIZ
extern WordIndex chrom_get_by_name (STRp (chrom_name)); // ZIP and PIZ

// chrom2ref
extern void chrom_calculate_ref2chrom (uint64_t num_ref_contigs);
extern void chrom_2ref_compress (Reference ref);
extern WordIndex chrom_2ref_seg_get (Reference ref, ConstVBlockP vb, WordIndex chrom_index);
#define chrom_2ref_seg_is_needed(did_i) ((did_i) == CHROM && (flag.reference & REF_ZIP_CHROM2REF))

extern void chrom_2ref_load (Reference ref);
static inline WordIndex chrom_2ref_piz_get (WordIndex chrom_index) 
    { return z_file->chrom2ref_map.len ? *B(WordIndex, z_file->chrom2ref_map, chrom_index) : chrom_index; }

