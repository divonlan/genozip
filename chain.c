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
    vb->contexts[CHAIN_TOPLEVEL].no_stons  = 
    vb->contexts[CHAIN_CHAIN]   .no_stons  = 
    vb->contexts[CHAIN_TSTRAND] .no_stons  = true; // keep in b250 so it can be eliminated as all_the_same

    vb->contexts[CHAIN_TSTART].flags.store = 
    vb->contexts[CHAIN_TEND]  .flags.store = 
    vb->contexts[CHAIN_QSTART].flags.store = 
    vb->contexts[CHAIN_QEND]  .flags.store = 
    vb->contexts[CHAIN_SIZE]  .flags.store = 
    vb->contexts[CHAIN_DT_DQ] .flags.store = STORE_INT;

    vb->contexts[CHAIN_TNAME] .flags.store =
    vb->contexts[CHAIN_QNAME] .flags.store = STORE_INDEX;
}

void chain_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - if changing fields, update chain_piz_filter
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 16,
        .items        = { { .dict_id = (DictId)dict_id_fields[CHAIN_CHAIN],   .seperator = {' '} },
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
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_DT_DQ, 0); // dt
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_DT_DQ, 0); // dq
        }
        SEG_EOL (CHAIN_EOL, false);

        num_alignments++;
    } 
    while (!is_last_alignment);

    // alignment set - IMPORTNAT - if changing fields, update chain_piz_filter
    SmallContainer alignment_set = { 
        .repeats             = num_alignments,
        .nitems_lo           = 4,
        .keep_empty_item_sep = true, // avoid double-deleting the space - only chain_piz_special_BACKSPACE should delete it, not container_reconstruct_do
        .filter_items        = true,
        .items               = { { .dict_id = (DictId)dict_id_fields[CHAIN_SIZE],  .seperator = {' '} },
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_DT_DQ], .seperator = {' '} }, // dt
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_DT_DQ]                     }, // dq
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

    Context *dtdq_ctx = &vb->contexts[CHAIN_DT_DQ];
    dtdq_ctx->last_value.i = dtdq_ctx->local.param = 0; // no dt and dq

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

//-------------------------
// "chain" buffer functions
//-------------------------

typedef struct {
    const char *primary_chrom; // pointer into dict of CHAIN_QNAME
    PosType primary_first_pos, primary_last_pos;

    const char *secondary_chrom; // pointer into dict of CHAIN_TNAME
    PosType secondary_first_pos, secondary_last_pos;
} ChainAlignment;

static Buffer chain = EMPTY_BUFFER;
static Mutex chain_mutex = {};
static PosType next_primary_pos=0, next_secondary_pos=0;

static void chain_display_alignment (const ChainAlignment *aln)
{
    fprintf (info_stream, "Primary: %s %"PRId64"-%"PRId64" Secondary: : %s %"PRId64"-%"PRId64"\n",
             aln->primary_chrom,   aln->primary_first_pos,   aln->primary_last_pos,
             aln->secondary_chrom, aln->secondary_first_pos, aln->secondary_last_pos);
}

// initialize alignment set from alignment set header
static inline void chain_piz_filter_init_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);
    next_primary_pos   = vb->contexts[CHAIN_QSTART].last_value.i; 
    next_secondary_pos = vb->contexts[CHAIN_TSTART].last_value.i; 
    mutex_unlock (chain_mutex);
}

// store dt value in local.param, as last_value.i will be overridden by dq, as they share the same context
static inline void chain_piz_filter_save_dt (VBlock *vb)
{
    Context *dtdq_ctx = &vb->contexts[CHAIN_DT_DQ];
    dtdq_ctx->local.param = dtdq_ctx->last_value.i; 
}

// ingest one alignment
static inline void chain_piz_filter_ingest_alignmet (VBlock *vb)
{
    mutex_lock (chain_mutex);

    buf_alloc_more (vb, &chain, 1, 0, ChainAlignment, 2, "chain"); 

    int64_t size = vb->contexts[CHAIN_SIZE].last_value.i;

    NEXTENT (ChainAlignment, chain) = (ChainAlignment){ 
        .primary_chrom       = ctx_get_snip_by_word_index (&vb->contexts[CHAIN_QNAME], vb->contexts[CHAIN_QNAME].last_value.i, 0, 0),
        .primary_first_pos   = next_primary_pos,
        .primary_last_pos    = next_primary_pos + size - 1,

        .secondary_chrom     = ctx_get_snip_by_word_index (&vb->contexts[CHAIN_TNAME], vb->contexts[CHAIN_TNAME].last_value.i, 0, 0),
        .secondary_first_pos = next_secondary_pos,
        .secondary_last_pos  = next_secondary_pos + size - 1
    };

    Context *dtdq_ctx   = &vb->contexts[CHAIN_DT_DQ];
    next_primary_pos   += size + dtdq_ctx->last_value.i; // dq
    next_secondary_pos += size + dtdq_ctx->local.param;  // dt
    
    mutex_unlock (chain_mutex);

    //chain_display_alignment (LASTENT (ChainAlignment, chain));

    vb->dont_show_curr_line = true;
}

// verify that adding up all alignments and gaps, results in the end position specified in the header
static inline void chain_piz_filter_verify_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);

    ASSERTE (next_primary_pos == vb->contexts[CHAIN_QEND].last_value.i,
                "Expecting alignments to add up to qEND=%"PRId64", but they add up to %"PRId64":\n%*s",
                vb->contexts[CHAIN_QEND].last_value.i, next_primary_pos, 
                (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    ASSERTE (next_secondary_pos == vb->contexts[CHAIN_TEND].last_value.i,
                "Expecting alignments to add up to tEND=%"PRId64", but they add up to %"PRId64":\n%*s",
                vb->contexts[CHAIN_TEND].last_value.i, next_secondary_pos, 
                (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    mutex_unlock (chain_mutex);
}

CONTAINER_FILTER_FUNC (chain_piz_filter)
{
    if (!flag.reading_chain) goto done; // only relevant to ingesting the chain due to --chain

    // before alignment-set first EOL but before alignments - initialize next_primary_pos and next_secondary_pos
    if (dict_id.num == dict_id_fields[CHAIN_TOPLEVEL] && item == 13) 
        chain_piz_filter_init_alignment_set (vb);

    // save dt before reconstructing dq (dt was just processed)
    else if (dict_id.num == dict_id_fields[CHAIN_SET] && item == 2) 
        chain_piz_filter_save_dt (vb);

    // before EOF of each alignment, ingest alignment
    else if (dict_id.num == dict_id_fields[CHAIN_SET] && item == 3) 
        chain_piz_filter_ingest_alignmet (vb);

    // before alignment-set second EOL and after alignments - verify that numbers add up
    else if (dict_id.num == dict_id_fields[CHAIN_TOPLEVEL] && item == 15) 
        chain_piz_filter_verify_alignment_set (vb);

done:    
    return true;    
}

// zip: import chain file
void chain_load_chain (void)
{
    ASSERTE0 (flag.reading_chain, "reading_chain is NULL");
    SAVE_FLAGS;

    buf_alloc_more (evb, &chain, 0, 1000, ChainAlignment, 1, "chain"); // must be allocated by I/O thread
    
    mutex_initialize (chain_mutex);

    z_file = file_open (flag.reading_chain, READ, Z_FILE, DT_CHAIN);    
    z_file->basename = file_basename (flag.reading_chain, false, "(chain-file)", NULL, 0);

    // save and reset flags that are intended to operate on the compressed file rather than the reference file
    flag.test = flag.md5 = flag.show_memory = flag.show_stats= flag.no_header =
    flag.header_one = flag.header_only = flag.regions = flag.show_index = flag.show_dict = 
    flag.show_b250 = flag.show_ref_contigs = flag.list_chroms = 0;
    flag.grep = flag.show_time = flag.unbind = 0;
    flag.dict_id_show_one_b250 = flag.dump_one_b250_dict_id = flag.dump_one_local_dict_id = DICT_ID_NONE;
    flag.show_one_dict = NULL;

    TEMP_VALUE (command, PIZ);

    piz_one_file (0, false, false);

    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtfile_genozip_to_txt_header.
    
    mutex_destroy (chain_mutex);

    flag.reading_chain = NULL;
}

