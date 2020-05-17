// ------------------------------------------------------------------
//   zip_gff3.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "vblock.h"
#include "buffer.h"
#include "zfile.h"
#include "zip.h"

// this function receives all lines of a vblock and processes them
// in memory to the compressed format. This thread then terminates, and the I/O thread writes the output.
void zip_gff3_compress_one_vb (VBlockP vb_)
{ 
    VBlockGFF3 *vb = (VBlockGFF3 *)vb_;

    // optional data sections
    COMPRESS_DATA_SECTION (SEC_RANDOM_POS_DATA, random_pos_data, uint32_t, COMP_LZMA, true);
    COMPRESS_DATA_SECTION (SEC_SEQ_DATA, seq_data, char, COMP_LZMA, true);
    COMPRESS_DATA_SECTION (SEC_NUMERIC_ID_DATA, dbxref_numeric_data, char, COMP_LZMA, true);
    COMPRESS_DATA_SECTION (SEC_ENST_DATA, enst_data, char, COMP_LZMA, true);
}
