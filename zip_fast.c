// ------------------------------------------------------------------
//   zip_fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "vblock.h"
#include "zip.h"
#include "optimize.h"

#define DATA_LINE(i) ENT (ZipDataLineFAST, vb->lines, i)

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
void zip_fast_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
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
void zip_fast_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
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
