// ------------------------------------------------------------------
//   me23.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"
#include "piz.h"
#include "dict_id.h"

void me23_seg_initialize (VBlock *vb)
{
    vb->contexts[ME23_CHROM].inst = CTX_INST_NO_STONS; // needs b250 node_index for random access
    vb->contexts[ME23_GENOTYPE].ltype = LT_SEQUENCE;
}

void me23_seg_finalize (VBlockP vb)
{
    // top level snip
    Container top_level = { 
        .repeats   = vb->lines.len,
        .flags     = CONTAINER_TOPLEVEL,
        .num_items = 5,
        .items     = { { (DictId)dict_id_fields[ME23_ID],       DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[ME23_CHROM],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[ME23_POS],      DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[ME23_GENOTYPE], DID_I_NONE, ""   },
                       { (DictId)dict_id_fields[ME23_EOL],      DID_I_NONE, ""   } }
    };

    container_seg_by_ctx (vb, &vb->contexts[ME23_TOPLEVEL], &top_level, 0, 0, 0);
}

const char *me23_seg_txt_line (VBlock *vb, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    GET_NEXT_ITEM ("RSID");
    seg_id_field (vb, (DictId)dict_id_fields[ME23_ID], field_start, field_len, true);

    GET_NEXT_ITEM ("CHROM");
    seg_chrom_field (vb, field_start, field_len);

    GET_NEXT_ITEM ("POS");
    seg_pos_field (vb, ME23_POS, ME23_POS, false, field_start, field_len, 0, field_len+1);
    random_access_update_pos (vb, ME23_POS);

    // Genotype (a combination of one or two bases or "--")
    GET_LAST_ITEM ("GENOTYPE");
    
    ASSERT (field_len == 1 || field_len == 2, "%s: Error in %s: expecting all genotype data to be 1 or 2 characters, but found one with %u: %.*s",
            global_cmd, txt_name, field_len, field_len, field_start);

    seg_add_to_local_fixed (vb, &vb->contexts[ME23_GENOTYPE], field_start, field_len); 
        
    char lookup[2] = { SNIP_LOOKUP, '0' + field_len };
    seg_by_did_i (vb, lookup, 2, ME23_GENOTYPE, field_len + 1);

    SEG_EOL (ME23_EOL, false);
    
    return next_field;
}

