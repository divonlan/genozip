// ------------------------------------------------------------------
//   writer_private.h
//   Copyright (C) 2021-2025 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "vblock.h"
#include "writer.h"

typedef struct {
    Mutex wait_for_data;   // initialized locked, unlocked when data is moved to txt_data and is_loaded is set
    VBlockP vb;            // data handed over main -> compute -> main -> writer
    VBIType vblock_i;      // this is the same as the index of this entry, and also of vb->vblock_i (after vb is loaded). for convenience.
    VBIType pair_vb_i;     // vb_i of pair vb, or 0 if there isn't any. if pair_vb_i > vblock_i then we are pair_1
    uint32_t num_lines;    // number of lines in this VB - read from section list (up to V14: from VB header). SAM MAIN: excluding gencomp lines
    CompIType comp_i;      // comp_i of this VB 
    bool is_loaded;        // has data moving to txt_data completed
    bool needs_recon;      // VB/TXT_HEADER should be reconstructed, but writing it depends on needs_write
    bool needs_write;      // VB/TXT_HEADER is written to output file, and hence is in recon_plan
    bool encountered;      // used by writer_create_txt_section_list
    bool fasta_contig_grepped_out; // FASTA: a VB that starts with SEQ data, inherits this from the last contig of the previous VB 
    bool vb_plan_expanded; // PLAN_VB_PLAN was expanded during creation of z_file->recon_plan, before creating reading_list, so no need to expand it again when updating the reading list
    bool range_vb_in_use;  // writer main loop: true if this vb contains PLAN_RANGE/PLAN_END_OF_VB items and is currently being written
    PairType pair;         // TXT_HEADER: used for FASTQ and SAM deep: may be NOT_PAIRED, PAIR_R1, PAIR_R2
    Buffer is_dropped_buf; // a Bits with a bit set is the line is marked for dropping
    BitsP is_dropped;      // pointer into is_dropped_buf

    // data only set for SAM MAIN VBs
    uint32_t num_gc_vbs;   // number of PRIM/DEPN VBs integrated into this MAIN VB 
    uint32_t num_gc_lines[2]; // SAM MAIN VB: number of PRIM/DEPN lines integrated into this MAIN VB
    VBIType first_gc_vb[2];// [0]=PRIM [1]=DEPN: the first PRIM/DEPN line in this MAIN VB, corresponds to this PRIM/DEPN VB and line_i within it
    uint32_t first_gc_vb_line_i[2];

    // data only set for SAM PRIM/DEPN VBs
    VBIType ending_vb;     // highest vb_i of any MAIN vb that consumes lines from this gc VB
} VbInfo; // used both for VBs and txt headers

#define VBINFO(vb_i) ({ ASSERTINRANGX ((vb_i), 1, z_file->num_vbs); B(VbInfo, z_file->vb_info[0], (vb_i)); })
#define VBINFO_NUM(v) BNUM (z_file->vb_info[0], (v))

extern void writer_add_trivial_plan (CompIType comp_i, PlanFlavor flavor);
extern void writer_add_one_vb_plan_prim_or_depn (VBIType vblock_i);
