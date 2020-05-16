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

    EXTENDED_FIELD_CTX (ME23_HT, dict_id_field (dict_id_make ("HT", 2)));

    buf_alloc (vb, &vb->genotype_data, 2 * vb->lines.len, 1, "genotype_data", vb->vblock_i);
    buf_alloc (vb, &vb->id_numeric_data, sizeof(uint32_t) * vb->lines.len, 1, "id_numeric_data", vb->vblock_i);    
    buf_alloc (vb, &vb->random_pos_data, vb->lines.len * sizeof (uint32_t), 1, "random_pos_data", vb->vblock_i);    
}             
             
const char *seg_me23_data_line (VBlock *vb_,   
                                const char *field_start_line)     // index in vb->txt_data where this line starts
{
    VBlockME23 *vb = (VBlockME23 *)vb_;

    const char *next_field, *field_start, *rsid_field_start;
    unsigned field_len=0, rsid_field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    // RSID
    rsid_field_start = field_start_line;
    next_field = seg_get_next_item (vb, rsid_field_start, &len, false, true, false, &rsid_field_len, &separator, NULL, "RSID");
    // wait before adding - we will tag a # to the ID if the row does not have a \r (it normally does)

    // CHROM
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, NULL, "CHROM");
    seg_chrom_field (vb_, field_start, field_len);

    // POS - store delta vs previous line
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, NULL, "POS");
    vb->last_pos = seg_pos_field (vb_, vb->last_pos, NULL, false, ME23_POS, field_start, field_len, "POS", true);

    random_access_update_pos (vb_, vb->last_pos);

    // Genotype (a combination of two bases or "--")
    field_start = next_field;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n
    next_field = seg_get_next_item (vb, field_start, &len, true, false, false, &field_len, &separator, &has_13, "GENOTYPE");
    
    ASSERT (field_len == 1 || field_len == 2, "%s: Error in %s: expecting all genotype data to be 1 or 2 characters, but found one with %u: %.*s",
            global_cmd, txt_name, field_len, field_len, field_start);

    // for haplotypes (X, Y, MT) add a '*' (override the newline - no harm)
    if (field_len==1) *(char*)(field_start + 1) = '*';
    
    seg_add_to_fixed_buf (vb_, &vb->genotype_data, field_start, 2); 
    vb->mtf_ctx[ME23_HT].txt_len += field_len + 1 + has_13;

    // Now, finalize RSID - if we DON'T have a \r (unexpected), then we add an extra bit.
    seg_id_field (vb_, &vb->id_numeric_data, (DictIdType)dict_id_fields[ME23_ID], SEC_ID_B250, SEC_NUMERIC_ID_DATA,
                  rsid_field_start, rsid_field_len, !has_13, true);
    
    return next_field;
}
