// ------------------------------------------------------------------
//   me23.c
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
#include "piz.h"
#include "dict_id.h"

#define DATA_LINE(i) ENT (ZipDataLineME23, vb->lines, i)

void me23_seg_initialize (VBlock *vb)
{
    vb->contexts[ME23_CHROM].flags = CTX_FL_NO_STONS; // needs b250 node_index for random access
}

const char *me23_seg_txt_line (VBlock *vb, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    GET_NEXT_ITEM ("RSID");
    seg_id_field (vb, (DictIdType)dict_id_fields[ME23_ID], field_start, field_len, true);

    GET_NEXT_ITEM ("CHROM");
    seg_chrom_field (vb, field_start, field_len);

    GET_NEXT_ITEM ("POS");
    seg_pos_field (vb, ME23_POS, ME23_POS, false, field_start, field_len, true);
    random_access_update_pos (vb, ME23_POS);

    // Genotype (a combination of one or two bases or "--")
    GET_LAST_ITEM ("GENOTYPE");
    
    ASSERT (field_len == 1 || field_len == 2, "%s: Error in %s: expecting all genotype data to be 1 or 2 characters, but found one with %u: %.*s",
            global_cmd, txt_name, field_len, field_len, field_start);

    seg_add_to_local_fixed (vb, &vb->contexts[ME23_GENOTYPE], field_start, field_len); 
    vb->contexts[ME23_GENOTYPE].ltype = CTX_LT_SEQUENCE;
    
    char lookup[2] = { SNIP_LOOKUP, '0' + field_len };
    seg_by_did_i (vb, lookup, 2, ME23_GENOTYPE, field_len + 1);

    SEG_EOL (ME23_EOL, false);
    //seg_by_did_i (vb, *has_13 ? "\r\n" : "\n", (*has_13) + 1, ME23_EOL, *has_13);
    
    return next_field;
}

void me23_piz_reconstruct_vb (VBlock *vb)
{
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        vb->line_i = vb->first_line + vb_line_i;

        piz_reconstruct_from_ctx (vb, ME23_ID,    '\t');
        piz_reconstruct_from_ctx (vb, ME23_CHROM, '\t');
        piz_reconstruct_from_ctx (vb, ME23_POS,   '\t');
        piz_reconstruct_from_ctx (vb, ME23_GENOTYPE, 0);
        piz_reconstruct_from_ctx (vb, ME23_EOL, 0);

        // after consuming the line's data, if it is not to be outputted - trim txt_data back to start of line
        if (vb->dont_show_curr_line) vb->txt_data.len = txt_data_start; 
    }
}
