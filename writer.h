// ------------------------------------------------------------------
//   writer.h
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// Reconstruction plan 
extern void writer_create_plan (void);

// Writer thread
extern void writer_handover_data (VBlockP *vb_p);
extern void writer_handover_txtheader (VBlockP *txt_header_vb_p);
extern void writer_finish_writing (bool is_last_txt_file);

// VBlock and Txt Header properties
extern unsigned writer_get_pair (VBIType vb_i, VBIType *pair_vb_i);
extern bool writer_does_txtheader_need_write (Section sec);
extern bool writer_does_txtheader_need_recon (Section sec);
extern bool writer_does_vb_need_recon (VBIType vb_i);
extern BitsP writer_get_is_dropped (VBIType vb_i);
extern bool writer_is_bgzf_by_compute_thread (VBIType vb_i);
extern bool writer_get_fasta_contig_grepped_out (VBIType vb_i);
extern void writer_set_fasta_contig_grepped_out (VBIType vb_i);
extern uint64_t writer_get_txt_line_i (VBlockP vb);
extern void writer_set_num_txtheader_lines (CompIType comp_i, uint32_t num_txtheader_lines);
