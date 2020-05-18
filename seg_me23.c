// ------------------------------------------------------------------
//   seg_me23.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"

#define DATA_LINE(i) ENT (ZipDataLineME23, vb->lines, i)

const char *seg_me23_data_line (VBlock *vb,   
                                const char *field_start_line)     // index in vb->txt_data where this line starts
{
    const char *next_field, *field_start, *rsid_field_start;
    unsigned field_len=0, rsid_field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    // RSID
    rsid_field_start = field_start_line;
    next_field = seg_get_next_item (vb, rsid_field_start, &len, false, true, false, &rsid_field_len, &separator, NULL, "RSID");
    seg_id_field (vb, (DictIdType)dict_id_fields[ME23_ID], rsid_field_start, rsid_field_len, true);

    // CHROM
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, NULL, "CHROM");
    seg_chrom_field (vb, field_start, field_len);

    // POS - store delta vs previous line (or full value in local)
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, false, true, false, &field_len, &separator, NULL, "POS");
    vb->last_pos = seg_pos_field (vb, vb->last_pos, NULL, false, ME23_POS, field_start, field_len, "POS", true);

    random_access_update_pos (vb, vb->last_pos);

    // Genotype (a combination of two bases or "--")
    field_start = next_field;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n
    next_field = seg_get_next_item (vb, field_start, &len, true, false, false, &field_len, &separator, &has_13, "GENOTYPE");
    
    ASSERT (field_len == 1 || field_len == 2, "%s: Error in %s: expecting all genotype data to be 1 or 2 characters, but found one with %u: %.*s",
            global_cmd, txt_name, field_len, field_len, field_start);

    // for haplotypes (X, Y, MT) add a '*' (override the newline - no harm)
    SAFE_ASSIGN (1, field_start + field_len, '*');
    seg_add_to_local_fixed (vb, &vb->mtf_ctx[ME23_GENOTYPE], field_start, 2); 
    SAFE_RESTORE (1);
    vb->mtf_ctx[ME23_GENOTYPE].txt_len += field_len + 1;

    // record whether there is a \r line ending
    seg_one_field (vb, has_13 ? "1" : "", has_13, ME23_HAS13, has_13);
    
    return next_field;
}
