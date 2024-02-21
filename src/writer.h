// ------------------------------------------------------------------
//   writer.h
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

// Reconstruction plan 
extern bool writer_create_plan (void);

// Writer thread
extern bool writer_handover_data (VBlockP *vb_p);
extern bool writer_handover_txtheader (VBlockP *txt_header_vb_p);
extern void writer_finish_writing (bool is_last_txt_file);

// BGZF threads
extern uint32_t writer_get_max_bgzf_threads (void);

// VBlock and Txt Header properties
extern bool writer_am_i_pair_2 (VBIType vb_i, VBIType *pair_1_vb_i);
extern bool writer_does_txtheader_need_write (CompIType comp_i);
extern bool writer_does_txtheader_need_recon (CompIType comp_i);
extern bool writer_does_vb_need_recon (VBIType vb_i);
extern bool writer_does_vb_need_write (VBIType vb_i);
extern BitsP writer_get_is_dropped (VBIType vb_i);
extern void writer_destroy_is_vb_info (void);
extern bool writer_get_fasta_contig_grepped_out (VBIType vb_i);
extern void writer_set_fasta_contig_grepped_out (VBIType vb_i);
extern int64_t writer_get_txt_line_i (VBlockP vb, LineIType line_in_vb);
extern void writer_set_num_txtheader_lines (CompIType comp_i, uint32_t num_txtheader_lines);

extern VBlockP wvb;