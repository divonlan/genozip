// ------------------------------------------------------------------
//   chain.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//
// handle UCSC Chain format: https://genome.ucsc.edu/goldenPath/help/chain.html

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "strings.h"
#include "piz.h"
#include "dict_id.h"
#include "stats.h"
#include "codec.h"
#include "version.h"
#include "chain.h"
#include "reconstruct.h"

//-----------------------
// TXTFILE stuff
//-----------------------

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t chain_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */)
{
    ARRAY (char, data, vb->txt_data)

    ASSERTE (*i >= 0 && *i < data_len, "*i=%d is out of range [0,%"PRIu64"]", *i, data_len);

    for (; *i >= (int32_t)first_i; (*i)--) {
        // in CHAIN - an "end of line" is two newlines (\n\n or \n\r\n)
        if (data[*i] == '\n' && 
            ((*i >= 1 && data[*i-1] == '\n') || (*i >= 2 && data[*i-1] == '\r' && data[*i-2] == '\n')))
            return data_len - *i - 1;
    }

    return -1; // cannot find end-of-line in the data starting first_i
}

//-----------------------
// Segmentation functions
//-----------------------

void chain_seg_initialize (VBlock *vb)
{
    vb->contexts[CHAIN_TOPLEVEL].no_stons = 
    vb->contexts[CHAIN_CHAIN]   .no_stons = 
    vb->contexts[CHAIN_TSTRAND] .no_stons = true; // keep in b250 so it can be eliminated as all_the_same

    vb->contexts[CHAIN_TSTART].flags.store = 
    vb->contexts[CHAIN_TEND].  flags.store = 
    vb->contexts[CHAIN_QSTART].flags.store = 
    vb->contexts[CHAIN_QEND].  flags.store = 
    vb->contexts[CHAIN_DT_DQ]. flags.store = STORE_INT;
}

void chain_seg_finalize (VBlockP vb)
{
    // top level snip
    SmallContainer top_level = { 
        .repeats     = vb->lines.len,
        .is_toplevel = true,
        .nitems_lo   = 16,
        .items       = { { .dict_id = (DictId)dict_id_fields[CHAIN_CHAIN],   .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_SCORE],   .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_TNAME],   .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_TQSIZE],  .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_TSTRAND], .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_TSTART],  .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_TEND],    .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_QNAME],   .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_TQSIZE],  .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_QSTRAND], .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_QSTART],  .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_QEND],    .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_ID]                          },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]                         },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_SET]                         }, 
                         { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]                         } },
    };

    container_seg_by_ctx (vb, &vb->contexts[CHAIN_TOPLEVEL], (ContainerP)&top_level, 0, 0, 0);
}

bool chain_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

static void chain_seg_size_field (VBlock *vb, const char *field_start, int field_len)
{
    int64_t size;
    ASSSEG (str_get_int (field_start, field_len, &size), field_start, "Expecting size to be an integer, but found %*s", field_len, field_start);

    Context *ctx = &vb->contexts[CHAIN_SIZE];

    if (vb->contexts[CHAIN_TEND].last_value.i - vb->contexts[CHAIN_TSTART].last_value.i == size) { // happens if the alignment set has only one alignment
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_SIZE }), 2, ctx, field_len + 1);
        ctx->numeric_only = false;
    }

    else
        seg_integer_or_not (vb, ctx, field_start, field_len, field_len + 1);
}

static void chain_seg_qend_field (VBlock *vb, const char *field_start, int field_len)
{
    Context *ctx = &vb->contexts[CHAIN_QEND];

    ASSSEG (str_get_int (field_start, field_len, &ctx->last_value.i), field_start, "Expecting qEnd to be an integer, but found %*s", field_len, field_start);

    if (vb->contexts[CHAIN_TEND].last_value.i - vb->contexts[CHAIN_TSTART].last_value.i ==
        ctx->last_value.i                     - vb->contexts[CHAIN_QSTART].last_value.i) {
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_QEND }), 2, ctx, field_len + 1);        
        ctx->numeric_only = false;
    }

    else
        seg_pos_field (vb, CHAIN_QEND, CHAIN_QSTART, false, field_start, field_len, 0, field_len + 1);
}

const char *chain_seg_txt_line (VBlock *vb, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    SEG_NEXT_ITEM_SP (CHAIN_CHAIN);
    SEG_NEXT_ITEM_SP (CHAIN_SCORE);
    SEG_NEXT_ITEM_SP (CHAIN_TNAME);    
    SEG_NEXT_ITEM_SP (CHAIN_TQSIZE); // only a handful a LN values in the file (one per chrom)
    SEG_NEXT_ITEM_SP (CHAIN_TSTRAND);

    GET_NEXT_ITEM_SP ("TSTART");
    seg_pos_field (vb, CHAIN_TSTART, CHAIN_TEND, false, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP ("TEND");
    seg_pos_field (vb, CHAIN_TEND, CHAIN_TSTART, false, field_start, field_len, 0,field_len + 1);

    SEG_NEXT_ITEM_SP (CHAIN_QNAME);    
    SEG_NEXT_ITEM_SP (CHAIN_TQSIZE); 
    SEG_NEXT_ITEM_SP (CHAIN_QSTRAND);

    GET_NEXT_ITEM_SP ("QSTART");
    seg_pos_field (vb, CHAIN_QSTART, CHAIN_QEND, false, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP ("QEND");
    chain_seg_qend_field (vb, field_start, field_len);

    // if ID is numeric, preferably store as delta, if not - normal snip
    GET_LAST_ITEM_SP ("ID");
    if (str_is_int (field_start, field_len))
        seg_pos_field (vb, CHAIN_ID, CHAIN_ID, false, field_start, field_len, 0, field_len + 1); // just a numeric delta    
    else
        seg_by_did_i (vb, field_start, field_len, CHAIN_ID, field_len + 1);

    SEG_EOL (CHAIN_EOL, true);

    // set containing 1 or more alignments of SIZE, DDST, SRC, last alignment in set is only SIZE
    uint32_t num_alignments=0;
    bool is_last_alignment=false;
    do {
        GET_MAYBE_LAST_ITEM_SP ("SIZE");
        chain_seg_size_field (vb, field_start, field_len);

        is_last_alignment = (separator == '\n');

        if (!is_last_alignment) {
            // note: we store both DIFFs in the same context as they are correlated (usually one of them is 0)
            SEG_NEXT_ITEM_SP (CHAIN_DT_DQ);
            SEG_LAST_ITEM_SP (CHAIN_DT_DQ);
        }
        // last line doesn't have DDST and DSRC - delete space seperator
        else {
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_DT_DQ, 0);
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_DT_DQ, 0);
        }
        SEG_EOL (CHAIN_EOL, false);

        num_alignments++;
    } 
    while (!is_last_alignment);

    // top level snip
    SmallContainer alignment_set = { 
        .repeats     = num_alignments,
        .nitems_lo   = 4,
        .keep_empty_item_sep = true, // avoid double-deleting the space - only chain_piz_special_BACKSPACE should delete it, not container_reconstruct_do
        .items       = { { .dict_id = (DictId)dict_id_fields[CHAIN_SIZE],  .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_DT_DQ], .seperator = {' '} },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_DT_DQ]                     },
                         { .dict_id = (DictId)dict_id_fields[CHAIN_EOL] }                     }
    };

    container_seg_by_ctx (vb, &vb->contexts[CHAIN_SET], (ContainerP)&alignment_set, 0, 0, 0);

    // Empty line after alignment set
    GET_LAST_ITEM ("EmptyLine");
    ASSSEG0 (!field_len, field_start, "Expecting an empty line after alignment set");
    SEG_EOL (CHAIN_EOL, false); 

    return next_field;
}

//--------------
// PIZ functions
//--------------

SPECIAL_RECONSTRUCTOR (chain_piz_special_BACKSPACE)
{
    ASSERTE0 (vb->txt_data.len, "vb->txt_data.len is 0");
    vb->txt_data.len--;
    return false;
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_QEND)
{
    new_value->i = vb->contexts[CHAIN_QSTART].last_value.i + 
                   vb->contexts[CHAIN_TEND].last_value.i - vb->contexts[CHAIN_TSTART].last_value.i;
    RECONSTRUCT_INT (new_value->i);
    return true; // new_value has been set
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_SIZE)
{
    new_value->i = vb->contexts[CHAIN_TEND].last_value.i - vb->contexts[CHAIN_TSTART].last_value.i;
    RECONSTRUCT_INT (new_value->i);
    return true; // new_value has been set
}

