// ------------------------------------------------------------------
//   zip_fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "vblock.h"
#include "buffer.h"
#include "zfile.h"
#include "header.h"
#include "zip.h"
#include "seg.h"
#include "txtfile.h"
#include "optimize.h"

#define DATA_LINE(i) ENT (ZipDataLineFAST, vb->lines, i)

// get lengths of sequence data (SEQ+E2) and quality data (QUAL+U2)
static void zip_fast_get_seq_len (VBlockFAST *vb, uint32_t *seq_len)
{
    *seq_len = 0;

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) 
        *seq_len  += DATA_LINE (vb_line_i)->seq_len;
}

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
static void zip_fast_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                               char **line_seq_data, uint32_t *line_seq_len,  // out 
                                               char **unused_data,  uint32_t *unused_len)
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);
    *line_seq_data = ENT (char, vb->txt_data, dl->seq_data_start);
    *line_seq_len  = dl->seq_len;
    *unused_data   = NULL;
    *unused_len    = 0;
}   

// callback function for compress to get data of one line (called by comp_compress_bzlib)
static void zip_fast_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                                char **line_qual_data, uint32_t *line_qual_len, // out
                                                char **unused_data,   uint32_t *unused_len) 
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->seq_len;
    *unused_data    = NULL;
    *unused_len     = 0;

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag_optimize_QUAL) optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

// this function receives all lines of a vblock and processes them
// in memory to the compressed format. This thread then terminates, and the I/O thread writes the output.
void zip_fast_compress_one_vb (VBlockP vb_)
{ 
    VBlockFAST *vb = (VBlockFAST *)vb_;

    // DESC subfields
    zip_generate_and_compress_subfields (vb_, &vb->desc_mapper);

    // SEQ
    uint32_t seq_len;
    zip_fast_get_seq_len (vb, &seq_len);
    zfile_compress_section_data_alg (vb_, SEC_SEQ_DATA,  NULL, zip_fast_get_start_len_line_i_seq,  seq_len, COMP_LZMA);

    // FASTQ only: QUAL
    if (vb->data_type == DT_FASTQ) 
        zfile_compress_section_data_alg (vb_, SEC_QUAL_DATA, NULL, zip_fast_get_start_len_line_i_qual, seq_len, COMP_BZ2);

    // FASTA obly: COMMENT
    else 
        zfile_compress_section_data (vb_, SEC_FASTA_COMMENT_DATA, &vb->comment_data);
  }
