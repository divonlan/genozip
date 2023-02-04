// ------------------------------------------------------------------
//   fastq_private.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "fastq.h"
#include "vblock.h"

typedef struct {
    TxtWord seq;
    uint32_t qual_index;         // start within vb->txt_data (qual_len==seq.len)
    bool dont_compress_QUAL : 1;
    bool monobase           : 1; // sequence is entirely of the same base
} ZipDataLineFASTQ;

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockFASTQ {
    VBLOCK_COMMON_FIELDS

    // pairing stuff - used if we are the 2nd file in the pair 
    uint32_t pair_vb_i;      // the equivalent vb_i in the first file, or 0 if this is the first file
    uint32_t pair_num_lines; // number of lines in the equivalent vb in the first file
    char *optimized_desc;    // base of desc in flag.optimize_DESC 
    uint32_t optimized_desc_len;
    uint64_t first_line;     // ZIP: used for optimize_DESC  
} VBlockFASTQ;

typedef VBlockFASTQ *VBlockFASTQP;

#define VB_FASTQ ((VBlockFASTQP)vb)

#define DATA_LINE(i) B(ZipDataLineFASTQ, vb->lines, i)

// SEQ
extern void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool deep);

// QUAL
extern void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, uint32_t qual_len, bool deep);

// Deep stuff
extern void fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(desc), STRp(seq), STRp(qual), bool *deep_desc, bool *deep_seq, bool *deep_qual);
extern void fastq_deep_seg_finalize_segconf (uint32_t n_lines);
extern void fastq_deep_seg_initialize (VBlockP vb);
extern void fastq_deep_show_entries_stats (void);
