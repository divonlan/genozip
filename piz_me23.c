// ------------------------------------------------------------------
//   piz_me23.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is the compute thread entry point. It receives all data of a variant block and processes it
// in memory to the uncompressed format. This thread then terminates the I/O thread writes the output.

#include "vblock.h"
#include "piz.h"
#include "endianness.h"
#include "zfile.h"
#include "buffer.h"
#include "move_to_front.h"
#include "regions.h"
#include "file.h"

void piz_me23_reconstruct_vb (VBlock *vb)
{
    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        uint32_t txt_line_i = vb->first_line + vb_line_i;

        DECLARE_SNIP;
        RECONSTRUCT_FROM_DICT (ME23_ID, true);

        uint32_t chrom_word_index = RECONSTRUCT_FROM_DICT (ME23_CHROM, true);
        RECONSTRUCT_FROM_DICT_POS (ME23_POS, vb->last_pos, true, NULL, true); // reconstruct from delta
        RECONSTRUCT_FROM_TABLESS_BUF (vb->mtf_ctx[ME23_GENOTYPE].local, vb->mtf_ctx[ME23_GENOTYPE].next_local, 2, false, "GENOTYPE");

        // remove the extra * added for ploidy=1 genotypes
        if (*LASTENT(char, vb->txt_data) == '*') vb->txt_data.len--; 

        // add the end-of-line
        LOAD_SNIP (ME23_HAS13);
        RECONSTRUCT (snip_len ? "\r\n" : "\n" , 1 + snip_len); // snip is "1" if there's \r\n and "" if not

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (chrom_word_index, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }
}
