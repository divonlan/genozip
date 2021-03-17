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
#include "liftover.h"

typedef struct {
    WordIndex query_chrom; // index in CHAIN_QNAME
    PosType query_first_pos, query_last_pos;

    WordIndex target_chrom; // index in CHAIN_TNAME - the nodes+dicts copied to the genozipped file (eg as VCF_oCHROM)
    PosType target_first_pos, target_last_pos;
} ChainAlignment;

static Buffer chain = EMPTY_BUFFER;   // immutable after loaded
static Mutex chain_mutex = {};        // protect chain wnile loading
static PosType next_query_pos=0, next_target_pos=0; // used for loading chain

bool chain_is_loaded = false; // global

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
    vb->contexts[CHAIN_QNAME]   .no_stons  =       
    vb->contexts[CHAIN_TNAME]   .no_stons  =       // (these 2 ^) need b250 node_index for reference
    vb->contexts[CHAIN_TOPLEVEL].no_stons  = 
    vb->contexts[CHAIN_CHAIN]   .no_stons  = 
    vb->contexts[CHAIN_TSTRAND] .no_stons  = true; // (these 3 ^) keep in b250 so it can be eliminated as all_the_same

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

    if (vb->last_int(CHAIN_TEND) - vb->last_int(CHAIN_TSTART) == size) { // happens if the alignment set has only one alignment
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

    if (vb->last_int(CHAIN_TEND) - vb->last_int(CHAIN_TSTART) ==
        ctx->last_value.i        - vb->last_int(CHAIN_QSTART)) {
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

void chain_piz_initialize (void)
{
    mutex_initialize (chain_mutex);
}

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
    new_value->i = vb->last_int(CHAIN_QSTART) + 
                   vb->last_int(CHAIN_TEND) - vb->last_int(CHAIN_TSTART);
    RECONSTRUCT_INT (new_value->i);
    return true; // new_value has been set
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_SIZE)
{
    new_value->i = vb->last_int(CHAIN_TEND) - vb->last_int(CHAIN_TSTART);
    RECONSTRUCT_INT (new_value->i);
    return true; // new_value has been set
}

// ---------------------------------------
// using the chain data in genozip --chain
// ---------------------------------------

static void chain_display_alignments (void) 
{
    ARRAY (ChainAlignment, aln, chain);

    for (uint32_t i=0; i < chain.len; i++) {
        const char *query_chrom  = ctx_get_words_snip (&z_file->contexts[CHAIN_QNAME], aln[i].query_chrom);
        const char *target_chrom = ctx_get_words_snip (&z_file->contexts[CHAIN_TNAME], aln[i].target_chrom);

        fprintf (info_stream, "Query: %s %"PRId64"-%"PRId64" Target: %s %"PRId64"-%"PRId64"\n",
                 query_chrom,  aln[i].query_first_pos,  aln[i].query_last_pos,
                 target_chrom, aln[i].target_first_pos, aln[i].target_last_pos);
    }
}

// initialize alignment set from alignment set header
static inline void chain_piz_filter_init_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);
    next_query_pos  = vb->last_int(CHAIN_QSTART); 
    next_target_pos = vb->last_int(CHAIN_TSTART); 
    mutex_unlock (chain_mutex);
}

// store dt value in local.param, as last_value.i will be overridden by dq, as they share the same context
static inline void chain_piz_filter_save_dt (VBlock *vb)
{
    vb->contexts[CHAIN_DT_DQ].local.param = vb->last_int(CHAIN_DT_DQ); 
}

// ingest one alignment
static inline void chain_piz_filter_ingest_alignmet (VBlock *vb)
{
    mutex_lock (chain_mutex);

    buf_alloc_more (vb, &chain, 1, 0, ChainAlignment, 2, "chain"); 

    int64_t size = vb->last_int(CHAIN_SIZE);

    NEXTENT (ChainAlignment, chain) = (ChainAlignment){ 
        .query_chrom      = vb->last_int(CHAIN_QNAME),
        .query_first_pos  = next_query_pos,
        .query_last_pos   = next_query_pos + size - 1,

        .target_chrom     = vb->last_int(CHAIN_TNAME),
        .target_first_pos = next_target_pos,
        .target_last_pos  = next_target_pos + size - 1
    };

    Context *dtdq_ctx   = &vb->contexts[CHAIN_DT_DQ];
    next_query_pos   += size + dtdq_ctx->last_value.i; // dq
    next_target_pos += size + dtdq_ctx->local.param;  // dt
    
    mutex_unlock (chain_mutex);

    vb->dont_show_curr_line = true;
}

// verify that adding up all alignments and gaps, results in the end position specified in the header
static inline void chain_piz_filter_verify_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);

    ASSINP (next_query_pos == vb->last_int(CHAIN_QEND),
            "Bad data in chain file %s: Expecting alignments to add up to qEND=%"PRId64", but they add up to %"PRId64":\n%*s",
            z_name, vb->last_int(CHAIN_QEND), next_query_pos, 
            (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    ASSINP (next_target_pos == vb->last_int(CHAIN_TEND),
            "Bad data in chain file %s: Expecting alignments to add up to tEND=%"PRId64", but they add up to %"PRId64":\n%*s",
            z_name, vb->last_int(CHAIN_TEND), next_target_pos, 
            (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    mutex_unlock (chain_mutex);
}

CONTAINER_FILTER_FUNC (chain_piz_filter)
{
    if (!flag.reading_chain) goto done; // only relevant to ingesting the chain due to --chain

    // before alignment-set first EOL and before alignments - initialize next_query_pos and next_target_pos
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

// sort the alignmants by (qname_index, qstart)
static int chain_sorter (const void *a_, const void *b_)
{   
    ChainAlignment *a = (ChainAlignment *)a_, *b = (ChainAlignment *)b_;

    if (a->query_chrom != b->query_chrom) 
        return a->query_chrom - b->query_chrom;
    else 
        return a->query_first_pos - b->query_first_pos;
}

// zip: load chain file as a result of genozip --chain
void chain_load (void)
{
    ASSERTE0 (flag.reading_chain, "reading_chain is NULL");
    SAVE_FLAGS;

    buf_alloc_more (evb, &chain, 0, 1000, ChainAlignment, 1, "chain"); // must be allocated by I/O thread
    
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

    // copy qname and tname dictionaries before z_file of the chain file is freed
    liftover_copy_data_from_chain_file();

    // sort the alignmants by (qname_index, qstart)
    qsort (chain.data, chain.len, sizeof (ChainAlignment), chain_sorter);
        
    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    if (flag.show_chain) {
        chain_display_alignments();
        exit_ok;
    }

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtfile_genozip_to_txt_header.
    
    flag.reading_chain = NULL;
    chain_is_loaded = true;
}

// get tname, tpos from qname, qpos (binary search on chain)
static LiftOverStatus chain_get_liftover_coords_do (WordIndex qname_index, PosType qpos, 
                                                    int32_t start, int32_t end,
                                                    WordIndex *tname_index, PosType *tpos) // out
{
    // case: no mapping
    if (end < start) {
        *tname_index = WORD_INDEX_NONE;
        *tpos = 0;
        return LO_NO_MAPPING;
    }

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = ENT (ChainAlignment, chain, mid);

    // case: success
    if (aln->query_chrom     == qname_index && 
        aln->query_first_pos <= qpos && 
        aln->query_last_pos >= qpos) {
            *tname_index = aln->target_chrom;
            *tpos = aln->target_first_pos + (qpos - aln->query_first_pos);
            return LO_OK;
        }

    // case: qname,qpos is less than aln - search lower half
    if (qname_index < aln->query_chrom ||
        (qname_index == aln->query_chrom && aln->query_first_pos > qpos))
        return chain_get_liftover_coords_do (qname_index, qpos, start, mid-1, tname_index, tpos);
    
    // case: qname,qpos is more than aln - search upper half
    else
        return chain_get_liftover_coords_do (qname_index, qpos, mid+1, end, tname_index, tpos);
}

LiftOverStatus chain_get_liftover_coords (WordIndex qname_index, PosType qpos, WordIndex *tname_index, PosType *tpos)
{
    return chain_get_liftover_coords_do (qname_index, qpos, 0, chain.len-1, tname_index, tpos);
}
