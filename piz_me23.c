// ------------------------------------------------------------------
//   piz_me23.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.

#include "vblock.h"
#include "piz.h"
#include "profiler.h"
#include "endianness.h"
#include "zfile.h"
#include "buffer.h"
#include "header.h"
#include "move_to_front.h"
#include "regions.h"
#include "file.h"

static void piz_me23_reconstruct_vb (VBlockME23 *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        uint32_t txt_line_i = vb->first_line + vb_line_i;
        bool doesnt_have_13=true; // note: for simplicity, we don't recover \r in case of --strip

        DECLARE_SNIP;
        IFNOTSTRIP(".",1) { RECONSTRUCT_ID (ME23_ID, &vb->id_numeric_data, &vb->next_id_numeric, &doesnt_have_13, true); }

        uint32_t chrom_word_index = RECONSTRUCT_FROM_DICT (ME23_CHROM, true);
        RECONSTRUCT_FROM_DICT_POS (ME23_POS, vb->last_pos, true, NULL, true); // reconstruct from delta
        RECONSTRUCT_FROM_TABLESS_BUF (vb->genotype_data, vb->next_genotype, 2, false, "HT_DATA");

        // remove the extra * added for ploidy=1 genotypes
        if (*LASTENT(char, vb->txt_data) == '*') 
            vb->txt_data.len--; 

        // add the end-of-line
        RECONSTRUCT (doesnt_have_13 ? "\n" : "\r\n" , 1 + !doesnt_have_13);

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (chrom_word_index, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

void piz_me23_uncompress_one_vb (VBlock *vb_)
{
    UNCOMPRESS_HEADER_AND_FIELDS (VBlockME23, true);  
    UNCOMPRESS_DATA_SECTION (SEC_NUMERIC_ID_DATA, id_numeric_data, char, false);    
    UNCOMPRESS_DATA_SECTION (SEC_HT_DATA, genotype_data, char, false);    
    UNCOMPRESS_DATA_SECTION (SEC_RANDOM_POS_DATA, random_pos_data, uint32_t, true);    

    piz_me23_reconstruct_vb (vb);
    UNCOMPRESS_DONE;
}

void piz_me23_read_one_vb (VBlock *vb_)
{ 
    PREPARE_TO_READ (VBlockME23, 5, SectionHeaderVbHeader);
    READ_FIELDS;
    READ_DATA_SECTION (SEC_NUMERIC_ID_DATA, false);
    READ_DATA_SECTION (SEC_HT_DATA, false);
    READ_DATA_SECTION (SEC_RANDOM_POS_DATA, true); // POS data that failed delta
    READ_DONE;
}
