// ------------------------------------------------------------------
//   fastq_private.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "file.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "deep.h"

#define DTYPE_QNAME     DTYPE_1
#define DTYPE_FASTQ_AUX DTYPE_2

#define MAX_DESC_FIELDS (MAX_FIELDS-100)

typedef struct {
    TxtWord seq;
    TxtWord qual;               // start within vb->txt_data (qual.len==seq.len except if condensed by codec_homp_compress)
    uint32_t sam_seq_len;       // Deep: seq_len of matching sam alignment
    uint32_t sam_seq_offset;    // Deep: start of matching SAM SEQ within FASTQ SEQ
    bool dont_compress_QUAL;    // true in case of Deep and fully copied from SAM
    bool deep_qual;             // QUAL is to be fully or partially copied from Deep
    bool monochar;              // sequence is entirely of the same character (eg: NNNNN)
    bool qual_is_trims_only;    // Deep: fastq_zip_qual has modified the QUAL data to be trims only 
} ZipDataLineFASTQ;

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockFASTQ {
    VBLOCK_COMMON_FIELDS

    // --------- current line - reset before every line by fastq_reset_line() ----------
    #define first_fastq_vb_zero_per_line piz_deep_flags
    PizDeepSeqFlags piz_deep_flags;// PIZ Deep
    uint32_t deep_seq_len;       // PIZ Deep
    uint32_t sam_seq_offset;     // PIZ Deep: offset of start of SEQ / QUAL copied from SAM with the FASTQ SEQ / QUAL

    #define after_fastq_vb_zero_per_line R1_vb_i
    // --------- END OF current line ----------

    // data used for segging R2
    uint32_t R1_vb_i;            // ZIP/PIZ R2: the equivalent vb_i in the R1 (vb_i >= 1)
    STR (R1_last_qname);         // ZIP R2: pointer into z_file->R1_last_qname: last reversed, nul-terminated canonical qname of the correspondining R1 VB
    uint32_t R1_num_lines;       // ZIP R2: number of reads (FASTQ lines) in the corresponding R1 vb
    int32_t R2_lowest_read;      // ZIP: in fastq_unconsumed: lowest/highest index in txt_data of reads that we tested against R1_last_qname
    int32_t R2_highest_read; 

    uint64_t first_line;         // ZIP R2: used for optimizing QNAME  

    bool has_extra;              // ZIP: a VB-private copy of segconf.has_extra

    #define FASTQ_NUM_TOP_LEVEL_FIELDS 16
    bool item_filter[FASTQ_NUM_TOP_LEVEL_FIELDS];  // PIZ: item filter - true if to reconstrut, false if not. index is item in toplevel container. length must be >= nitems_lo of all genozip verions

    // stats
    uint32_t deep_stats[NUM_DEEP_STATS_ZIP];  // ZIP: stats collection regarding Deep

    Multiplexer2 mux_ultima_c;
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

// SAUX
extern bool fastq_segconf_analyze_saux (VBlockFASTQP vb, STRp(saux));
extern void fastq_seg_saux (VBlockFASTQP vb, STRp(saux));

// Agilent stuff
extern void agilent_seg_initialize (VBlockP vb);
extern void agilent_seg_RX (VBlockP vb, ContextP ctx, STRp(rx), unsigned add_bytes); // RX and QX are also in sam_private.h.
extern void agilent_seg_QX (VBlockP vb, ContextP ctx, STRp(qx), unsigned add_bytes);

// SEQ
extern void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool deep);
extern bool fastq_piz_R1_test_aligned (VBlockFASTQP vb);

// QUAL
extern void fastq_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qual));

// Deep stuff
extern void fastq_deep_zip_initialize (void);
extern bool fastq_seg_deep (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(qname), STRp(qname2), STRp(seq), STRp(qual), uint32_t *uncanonical_suffix_len);
extern void fastq_deep_seg_finalize_segconf (uint32_t n_lines);
extern void fastq_deep_seg_initialize (VBlockFASTQP vb);
extern void fastq_deep_seg_QNAME (VBlockFASTQP vb, Did did_i, STRp(qname), uint32_t uncanonical_suffix_len, unsigned add_bytes);
extern void fastq_deep_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), ContextP bitmap_ctx, ContextP nonref_ctx);
extern void fastq_deep_seg_QUAL (VBlockFASTQP vb, ZipDataLineFASTQ *dl, ContextP qual_ctx, uint32_t qual_len);
extern void fastq_deep_zip_after_compute (VBlockFASTQP vb);
extern void fastq_deep_zip_finalize (void);
extern void fastq_deep_piz_wait_for_deep_data (void);
