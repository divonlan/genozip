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

    const char *snip;
    uint32_t snip_len, chrom_word_index;
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        uint32_t txt_line_i = vb->first_line + vb_line_i;
        bool doesnt_have_13=true; // note: for simplicity, we don't recover \r in case of --strip

        IFNOTSTRIP(".",1) {
            LOAD_SNIP (ME23_ID);        
            piz_reconstruct_id ((VBlockP)vb, &vb->id_numeric_data, &vb->next_numeric_id, snip, snip_len, &doesnt_have_13);
        }

        chrom_word_index = RECONSTRUCT_FROM_DICT (ME23_CHROM);
        RECONSTRUCT_FROM_DICT_POS (ME23_POS, true, NULL, true); // reconstruct from delta
        RECONSTRUCT_FROM_TABLESS_BUF (vb->genotype_data, vb->next_genotype, 2, false, "HT_DATA");

        // remove the extra * added for ploidy=1 genotypes
        if (*LASTENT(char, vb->txt_data) == '*') 
            vb->txt_data.len--; 

        // add the end-of-line
        buf_add (&vb->txt_data, doesnt_have_13 ? "\n" : "\r\n" , 1 + !doesnt_have_13);

        // after consuming sections' data, if this line is not to be outputed - shorten txt_data back to start of line
        if (flag_regions && !regions_is_site_included (chrom_word_index, vb->last_pos))
            vb->txt_data.len = txt_data_start; // remove excluded line
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

static void piz_me23_uncompress_all_sections (VBlockME23 *vb)
{
     ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->lines.len        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    unsigned section_i=1;

    // uncompress the fields     
    piz_uncompress_fields ((VBlockP)vb, section_index, &section_i);
    
    // uncompress the ID and Genotype data
    SectionHeader *id_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, id_header, &vb->id_numeric_data, "id_numeric_data", SEC_NUMERIC_ID_DATA);    

    SectionHeader *gt_header  = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, gt_header, &vb->genotype_data, "genotype_data", SEC_HT_DATA);    
}

void piz_me23_uncompress_one_vb (VBlock *vb_)
{
    START_TIMER;

    VBlockME23 *vb = (VBlockME23 *)vb_;

    piz_me23_uncompress_all_sections (vb);

    piz_me23_reconstruct_vb (vb);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway
    
    COPY_TIMER (vb->profile.compute);
}
