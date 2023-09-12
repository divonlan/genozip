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
#include "file.h"

typedef struct {
    TxtWord seq;
    TxtWord qual;                // start within vb->txt_data (qual.len==seq.len except if condensed by codec_homp_compress)
    uint32_t sam_seq_len;        // Deep: seq_len of matching sam alignment
    bool dont_compress_QUAL : 1; // true in case --deep and fully copied from SAM
    bool monochar           : 1; // sequence is entirely of the same character (eg: NNNNN)
} ZipDataLineFASTQ;

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockFASTQ {
    VBLOCK_COMMON_FIELDS

    // pairing stuff - used if we are the 2nd file in the pair 
    uint32_t pair_vb_i;          // ZIP/PIZ: in R2: the equivalent vb_i in the R1 (vb_i >= 1), or 0 if this is R1
    uint32_t pair_num_lines;     // R2: number of reads (FASTQ lines) in the equivalent vb in the R1
    uint32_t pair_txt_data_len;  // ZIP R2: populated if flag.debug

    STRw (optimized_desc);       // base of desc in flag.optimize_DESC 
    uint64_t first_line;         // ZIP: used for optimize_DESC  

    bool has_extra;              // ZIP: a VB-private copy of segconf.has_extra

    #define FASTQ_NUM_TOP_LEVEL_FIELDS 16
    bool item_filter[FASTQ_NUM_TOP_LEVEL_FIELDS];  // PIZ: item filter - true if to reconstrut, false if not. index is item in toplevel container. length must be >= nitems_lo of all genozip verions

    // stats
    uint32_t deep_stats[NUM_DEEP_STATS];  // ZIP: stats collection regarding Deep

    Multiplexer2 mux_SULTIMA;
} VBlockFASTQ;

typedef VBlockFASTQ *VBlockFASTQP;

#define VB_FASTQ ((VBlockFASTQP)vb)

#define DATA_LINE(i) B(ZipDataLineFASTQ, vb->lines, i)

// DESC
extern void fastq_seg_QNAME (VBlockFASTQP vb, STRp(qname), uint32_t line1_len, bool deep, uint32_t uncanonical_suffix_len);
extern void fastq_seg_DESC (VBlockFASTQP vb, STRp(desc), bool deep_qname2, uint32_t uncanonical_suffix_len);
extern void fastq_seg_LINE3 (VBlockFASTQP vb, STRp(qline3), STRp(qname), STRp(desc));
extern void fastq_segconf_analyze_DESC (VBlockFASTQP vb, STRp(desc));
extern bool fastq_is_line3_copy_of_line1 (STRp(qname), STRp(line3), uint32_t desc_len);

// SEQ
extern void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool deep);
extern bool fastq_piz_R1_test_aligned (VBlockFASTQP vb);

// QUAL
extern void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual), bool deep);

// Deep stuff
extern void fastq_deep_zip_initialize (void);
extern void fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(qname2), STRp(seq), STRp(qual), bool *deep_qname, bool *deep_seq, bool *deep_qual, uint32_t *uncanonical_suffix_len);
extern void fastq_deep_seg_finalize_segconf (uint32_t n_lines);
extern void fastq_deep_seg_initialize (VBlockFASTQP vb);
extern void fastq_deep_seg_QNAME (VBlockFASTQP vb, Did did_i, STRp(qname), uint32_t uncanonical_suffix_len, unsigned add_bytes);
extern void fastq_deep_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), ContextP bitmap_ctx, ContextP nonref_ctx);
extern void fastq_deep_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, ContextP qual_ctx, uint32_t qual_len);
extern void fastq_deep_zip_after_compute (VBlockFASTQP vb);
extern void fastq_deep_zip_show_entries_stats (void);
extern void fastq_deep_piz_wait_for_deep_data (void);
