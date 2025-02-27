// ------------------------------------------------------------------
//   fasta_private.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "vblock.h"
#include "contigs.h"

typedef struct VBlockFASTA {    
    VBLOCK_COMMON_FIELDS

    bool contig_grepped_out;
    // note: last_line is initialized to FASTA_LINE_SEQ (=0) so that a ; line as the first line of the VB is interpreted as a description, not a comment
    enum { FASTA_LINE_SEQ, FASTA_LINE_DESC, FASTA_LINE_COMMENT } last_line; // ZIP & PIZ (FASTA only, not REF)

    bool vb_has_no_newline;         // ZIP: VB ends contains part of a sequence, not ending in a newline, which is continued in the next VB
    bool ra_initialized;            // ZIP: RA was initialized for this VB
    uint32_t lines_this_contig;     // ZIP

    // --make-reference
    bool has_contig_metadata;       // used by make-reference
    ContigMetadata contig_metadata; 
} VBlockFASTA, *VBlockFASTAP;

#define VB_FASTA ((VBlockFASTAP)vb)
