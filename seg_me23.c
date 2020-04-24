// ------------------------------------------------------------------
//   seg_me23.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "header.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"

#define DATA_LINE(i) ENT (ZipDataLineME23, vb->lines, i)

// called from seg_all_data_lines
void seg_me23_initialize (VBlock *vb_)
{
    VBlockME23 *vb = (VBlockME23 *)vb_;

    buf_alloc (vb, &vb->genotype_data, 2 * vb->lines.len, 1, "genotype_data", vb->vblock_i);
    buf_alloc (vb, &vb->rsid_data, 12 * vb->lines.len, 1, "rsid_data", vb->vblock_i);    
}             
             
const char *seg_me23_data_line (VBlock *vb_,   
                                const char *field_start_line)     // index in vb->txt_data where this line starts
{
    VBlockME23 *vb = (VBlockME23 *)vb_;

    const char *next_field, *field_start, *rsid_field_start;
    unsigned field_len=0, rsid_field_len=0;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    // RSID
    rsid_field_start = field_start_line;
    next_field = seg_get_next_item (vb, rsid_field_start, &len, false, true, false, &rsid_field_len, &separator, &has_13, "RSID");
    // wait before adding - we will tag a # to the ID if the row does not have a \r (it normally does)

    // CHROM
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "CHROM");
    uint32_t chrom_node_index = seg_one_field (vb, field_start, field_len, ME23_CHROM);

    random_access_update_chrom (vb_, chrom_node_index);

    // POS - store delta vs previous line
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, &has_13, "POS");
    vb->last_pos = seg_pos_field (vb_, vb->last_pos, ME23_POS, SEC_POS_B250, field_start, field_len, "POS");

    random_access_update_pos (vb_, vb->last_pos);

    // Genotype (a combination of two bases or "--")
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, true, false, false, &field_len, &separator, &has_13, "GENOTYPE");
    
    ASSERT (field_len == 1 || field_len == 2, "%s: Error in %s: expecting all genotype data to be 1 or 2 characters, but found one with %u: %.*s",
            global_cmd, txt_name, field_len, field_len, field_start);

    // for haplotypes (X, Y, MT) add a '*' (override the newline - no harm)
    if (field_len==1) *(char*)(field_start + 1) = '*';
    
    seg_add_to_data_buf (vb_, &vb->genotype_data, SEC_GT_DATA, field_start, 2, 0, field_len + 1 + has_13); 
    
    // Now, finalize RSID - if we DON'T have a \r (unexpected), then we add an extra bit.
    seg_id_field (vb_, &vb->rsid_data, ME23_ID, (char*)rsid_field_start, rsid_field_len, !has_13);
    //if (!has_13) ((char*)rsid_field_start)[rsid_field_len] = '#';
    
    //seg_add_to_data_buf (vb_, &vb->rsid_data, SEC_ID_DATA, rsid_field_start, rsid_field_len + !has_13, 
    //                     '\t', rsid_field_len+1 /* for the \t */); 
/*
    for (unsigned i=0; i < rsid_field_len; i++)
        if (IS_DIGIT (rsid_field_start[i])) { // first digit
            uint32_t id_num = BGEN32 (atoi (&rsid_field_start[i]));
            seg_add_to_data_buf (vb_, &vb->rsid_data, SEC_ID_DATA, (char*)&id_num, sizeof (id_num), 0, rsid_field_len+1);
            break;
        }
*/
    return next_field;
}
