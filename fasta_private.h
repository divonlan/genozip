// ------------------------------------------------------------------
//   fasta_private.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FASTA_PRIVATE_INCLUDED
#define FASTA_PRIVATE_INCLUDED

#include "genozip.h"
#include "vblock.h"
#include "fasta.h"

// IMPORTANT: if changing fields in VBlockFASTA, also update fasta_vb_release_vb 
typedef struct VBlockFASTA {    
    VBLOCK_COMMON_FIELDS

    bool contig_grepped_out;
    // note: last_line is initialized to FASTA_LINE_SEQ (=0) so that a ; line as the first line of the VB is interpreted as a description, not a comment
    enum { FASTA_LINE_SEQ, FASTA_LINE_DESC, FASTA_LINE_COMMENT } last_line; // ZIP & PIZ (FASTA only, not REF)

    uint32_t lines_this_contig;     // ZIP

    // caching of seq line snips
    uint32_t std_line_len;          // ZIP: determined by first non-first-line seq line in VB 
    WordIndex std_line_node_index;  // ZIP: node index for non-first lines, with same length as std_line_len

    // --make-reference
    bool has_contig_metadata;       // used by make-reference
    ContigMetadata contig_metadata; 
} VBlockFASTA;

#endif


