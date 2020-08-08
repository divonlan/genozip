// ------------------------------------------------------------------
//   fast_private.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FAST_PRIVATE_INCLUDED
#define FAST_PRIVATE_INCLUDED

#include "vblock.h"
#include "fastq.h"
#include "fasta.h"

typedef struct {
    uint32_t seq_data_start, qual_data_start; // start within vb->txt_data
    uint32_t seq_len;                         // length within vb->txt_data (in case of FASTQ, this length is also applies to quality, per FASTQ spec)
} ZipDataLineFAST;

// IMPORTANT: if changing fields in VBlockFAST, also update vb_fast_release_vb 
typedef struct VBlockFAST {     // for FASTA and FASTQ
    VBLOCK_COMMON_FIELDS
    SubfieldMapper desc_mapper; // FASTA and FASTQ - ZIP & PIZ

    // FASTQ - pairing stuff - used if we are the 2nd file in the pair 
    uint32_t pair_vb_i;      // the equivalent vb_i in the first file, or 0 if this is the first file
    uint32_t pair_num_lines; // number of lines in the equivalent vb in the first file
    char *optimized_desc;    // base of desc in flag_optimize_DESC 
    uint32_t optimized_desc_len;
    
    // FASTA stuff
    bool contig_grepped_out;
    // note: last_line is initialized to FASTA_LINE_SEQ (=0) so that a ; line as the first line of the VB is interpreted as a description, not a comment
    enum { FASTA_LINE_SEQ, FASTA_LINE_DESC, FASTA_LINE_COMMENT } last_line; // FASTA ZIP
} VBlockFAST;

#define DATA_LINE(i) ENT (ZipDataLineFAST, vb->lines, i)

extern bool fast_piz_test_grep (VBlockFAST *vb);
extern bool fasta_initialize_contig_grepped_out (VBlockFAST *vb, bool does_vb_have_any_desc, bool last_desc_in_this_vb_matches_grep);
extern bool fasta_is_grepped_out_due_to_regions (VBlockFAST *vb, const char *line_start);

extern Structured structured_DESC;

#define dict_id_is_fast_desc_sf dict_id_is_type_1
#define dict_id_fast_desc_sf dict_id_type_1

#endif
