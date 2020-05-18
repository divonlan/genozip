// ------------------------------------------------------------------
//   zip_sam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "vblock.h"
#include "zip.h"
#include "optimize.h"

#define DATA_LINE(i) ENT (ZipDataLineSAM, vb->lines, i)

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
void zip_sam_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                       char **line_seq_data, uint32_t *line_seq_len,  // out 
                                       char **line_e2_data,  uint32_t *line_e2_len)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    *line_seq_data = ENT (char, vb->txt_data, dl->seq_data_start);
    *line_seq_len  = dl->seq_data_len;
    *line_e2_data  = dl->e2_data_start ? ENT (char, vb->txt_data, dl->e2_data_start) : NULL;
    *line_e2_len   = dl->e2_data_len;

    // if SEQ is just "*" (i.e. unavailable) replace it by " " for consistency with QUAL (= simplify PIZ code)
    if (dl->seq_data_len == 1 && (*line_seq_data)[0] == '*') 
        *line_seq_data = " "; // pointer to static string
}   

// callback function for compress to get data of one line (called by comp_compress_bzlib)
void zip_sam_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                        char **line_qual_data, uint32_t *line_qual_len, // out
                                        char **line_u2_data,   uint32_t *line_u2_len) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->qual_data_len;
    *line_u2_data   = dl->u2_data_start ? ENT (char, vb->txt_data, dl->u2_data_start) : NULL;
    *line_u2_len    = dl->u2_data_len;

    // if QUAL is just "*" (i.e. unavailable) replace it by " " because '*' is a legal PHRED quality value that will confuse PIZ
    if (dl->qual_data_len == 1 && (*line_qual_data)[0] == '*') 
        *line_qual_data = " "; // pointer to static string

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    else if (flag_optimize_QUAL) {
        optimize_phred_quality_string (*line_qual_data, *line_qual_len);
        if (*line_u2_data) optimize_phred_quality_string (*line_u2_data, *line_u2_len);
    }
}

// callback function for compress to get data of one line
void zip_sam_get_start_len_line_i_bd (VBlock *vb, uint32_t vb_line_i, 
                                      char **line_bd_data, uint32_t *line_bd_len,  // out 
                                      char **unused1,  uint32_t *unused2)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_bd_data = dl->bd_data_len ? ENT (char, vb->txt_data, dl->bd_data_start) : NULL;
    *line_bd_len  = dl->bd_data_len;
    *unused1 = NULL;
    *unused2 = 0;
}   

// callback function for compress to get BI data of one line
// if BD data is available for this line too - we output the character-wise delta as they are correlated
// we expect the length of BI and BD to be the same, the length of the sequence. this is enforced by seg_sam_optional_field.
void zip_sam_get_start_len_line_i_bi (VBlock *vb, uint32_t vb_line_i, 
                                      char **line_bi_data, uint32_t *line_bi_len,  // out 
                                      char **unused1,  uint32_t *unused2)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    if (dl->bi_data_len && dl->bd_data_len) {

        ASSERT (dl->bi_data_len == dl->bd_data_len, "Error: expecting dl->bi_data_len=%u to be equal to dl->bd_data_len=%u",
                dl->bi_data_len, dl->bd_data_len);

        char *bi       = ENT (char, vb->txt_data, dl->bi_data_start);
        const char *bd = ENT (char, vb->txt_data, dl->bd_data_start);

        // calculate character-wise delta
        for (unsigned i=0; i < dl->bi_data_len; i++) *(bi++) -= *(bd++);
    }

    *line_bi_data = dl->bi_data_len ? ENT (char, vb->txt_data, dl->bi_data_start) : NULL;
    *line_bi_len  = dl->bi_data_len;
    *unused1 = NULL;
    *unused2 = 0;
}   
