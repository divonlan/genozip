// ------------------------------------------------------------------
//   chain.c
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//
// handle UCSC Chain format: https://genome.ucsc.edu/goldenPath/help/chain.html
// also here: http://genomewiki.ucsc.edu/index.php/Chains_Nets#:~:text=5%20FAQs%3F-,Basic%20definitions,and%20mouse%20is%20the%20dst.
//
// Note on terminology: I find Chain file terminology very confusing an error-prone. I use a different terminology:
// Chain: Target Genozip: Src, Primary
// Chain: Query  Genozip: Dst, Laft (i.e. lifted-over, "laft" being an alternative past tense of "lift")

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
    WordIndex src_chrom; // index in CHAIN_NAMESRC
    PosType src_first_1pos, src_last_1pos;

    WordIndex dst_chrom; // index in CHAIN_NAMEDST - the nodes+dicts copied to the genozipped file (eg as VCF_oCHROM)
    PosType dst_first_1pos, dst_last_1pos;
} ChainAlignment;

static Buffer chain = EMPTY_BUFFER;   // immutable after loaded
static Mutex chain_mutex = {};        // protect chain wnile loading
static PosType next_dst_0pos=0, next_src_0pos=0; // used for loading chain

char *chain_filename = NULL; // global - chain filename

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
    vb->contexts[CHAIN_NAMEDST] .no_stons  =       
    vb->contexts[CHAIN_NAMESRC] .no_stons  =       // (these 2 ^) need b250 node_index for reference
    vb->contexts[CHAIN_TOPLEVEL].no_stons  = 
    vb->contexts[CHAIN_CHAIN]   .no_stons  = 
    vb->contexts[CHAIN_STRNDSRC].no_stons  = true; // (these 3 ^) keep in b250 so it can be eliminated as all_the_same

    vb->contexts[CHAIN_STARTSRC].flags.store = 
    vb->contexts[CHAIN_ENDSRC]  .flags.store = 
    vb->contexts[CHAIN_STARTDST].flags.store = 
    vb->contexts[CHAIN_ENDDST]  .flags.store = 
    vb->contexts[CHAIN_SIZE]    .flags.store = 
    vb->contexts[CHAIN_GAPS]    .flags.store = STORE_INT;

    vb->contexts[CHAIN_NAMESRC] .flags.store =
    vb->contexts[CHAIN_NAMEDST] .flags.store = STORE_INDEX;
}

void chain_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - if changing fields, update chain_piz_filter
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 16,
        .items        = { { .dict_id = (DictId)dict_id_fields[CHAIN_CHAIN],    .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_SCORE],    .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_NAMESRC],  .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_TQSIZE],   .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STRNDSRC], .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STARTSRC], .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_ENDSRC],   .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_NAMEDST],  .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_TQSIZE],   .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STRNDDST], .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STARTDST], .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_ENDDST],   .seperator = {' '} },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_ID]                           },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]                          },
                          { .dict_id = (DictId)dict_id_fields[CHAIN_SET]                          }, 
                          { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]                          } },
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

    if (vb->last_int(CHAIN_ENDSRC) - vb->last_int(CHAIN_STARTSRC) == size) { // happens if the alignment set has only one alignment
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_SIZE }), 2, ctx, field_len + 1);
        ctx->numeric_only = false;
    }

    else
        seg_integer_or_not (vb, ctx, field_start, field_len, field_len + 1);
}

static void chain_seg_dst_end_field (VBlock *vb, const char *field_start, int field_len)
{
    Context *ctx = &vb->contexts[CHAIN_ENDDST];

    ASSSEG (str_get_int (field_start, field_len, &ctx->last_value.i), field_start, "Expecting qEnd to be an integer, but found %*s", field_len, field_start);

    if (vb->last_int(CHAIN_ENDSRC) - vb->last_int(CHAIN_STARTSRC) ==
        ctx->last_value.i        - vb->last_int(CHAIN_STARTDST)) {
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_ENDDST }), 2, ctx, field_len + 1);        
        ctx->numeric_only = false;
    }

    else
        seg_pos_field (vb, CHAIN_ENDDST, CHAIN_STARTDST, false, field_start, field_len, 0, field_len + 1);
}

const char *chain_seg_txt_line (VBlock *vb, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    SEG_NEXT_ITEM_SP (CHAIN_CHAIN);
    SEG_NEXT_ITEM_SP (CHAIN_SCORE);
    SEG_NEXT_ITEM_SP (CHAIN_NAMESRC);    
    SEG_NEXT_ITEM_SP (CHAIN_TQSIZE); // only a handful a LN values in the file (one per chrom)
    SEG_NEXT_ITEM_SP (CHAIN_STRNDSRC);

    GET_NEXT_ITEM_SP ("tSTART");
    seg_pos_field (vb, CHAIN_STARTSRC, CHAIN_ENDSRC, false, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP ("tEND");
    seg_pos_field (vb, CHAIN_ENDSRC, CHAIN_STARTSRC, false, field_start, field_len, 0,field_len + 1);

    SEG_NEXT_ITEM_SP (CHAIN_NAMEDST);    
    SEG_NEXT_ITEM_SP (CHAIN_TQSIZE); 
    SEG_NEXT_ITEM_SP (CHAIN_STRNDDST);

    GET_NEXT_ITEM_SP ("qSTART");
    seg_pos_field (vb, CHAIN_STARTDST, CHAIN_ENDDST, false, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP ("qEND");
    chain_seg_dst_end_field (vb, field_start, field_len);

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
            SEG_NEXT_ITEM_SP (CHAIN_GAPS);
            SEG_LAST_ITEM_SP (CHAIN_GAPS);
        }
        // last line doesn't have DDST and DSRC - delete space seperator
        else {
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // src_gap
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // dst_gap
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
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_GAPS], .seperator = {' '} }, // src_gap
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_GAPS]                     }, // dst_gap
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

// called after reconstructing the txt header and before compute threads
void chain_piz_initialize (void)
{
    mutex_initialize (chain_mutex);
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_BACKSPACE)
{
    ASSERTE0 (vb->txt_data.len, "vb->txt_data.len is 0");
    vb->txt_data.len--;

    Context *gaps_ctx = &vb->contexts[CHAIN_GAPS];
    gaps_ctx->last_value.i = gaps_ctx->local.param = 0; // no src_gap and dst_gap

    return false;
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_ENDDST)
{
    new_value->i = vb->last_int(CHAIN_STARTDST) + 
                   vb->last_int(CHAIN_ENDSRC) - vb->last_int(CHAIN_STARTSRC);
    RECONSTRUCT_INT (new_value->i);
    return true; // new_value has been set
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_SIZE)
{
    new_value->i = vb->last_int(CHAIN_ENDSRC) - vb->last_int(CHAIN_STARTSRC);
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
        const char *dst_chrom = ctx_get_words_snip (&z_file->contexts[CHAIN_NAMEDST], aln[i].dst_chrom);
        const char *src_chrom = ctx_get_words_snip (&z_file->contexts[CHAIN_NAMESRC], aln[i].src_chrom);

        fprintf (info_stream, "Primary: %s %"PRId64"-%"PRId64" Laft: %s %"PRId64"-%"PRId64"\n",
                 src_chrom, aln[i].src_first_1pos, aln[i].src_last_1pos,
                 dst_chrom, aln[i].dst_first_1pos, aln[i].dst_last_1pos);
    }
}

// initialize alignment set from alignment set header
static inline void chain_piz_filter_init_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);
    next_dst_0pos = vb->last_int(CHAIN_STARTDST); 
    next_src_0pos = vb->last_int(CHAIN_STARTSRC); 
    mutex_unlock (chain_mutex);
}

// store src_gap value in local.param, as last_value.i will be overridden by dst_gap, as they share the same context
static inline void chain_piz_filter_save_src_gap (VBlock *vb)
{
    vb->contexts[CHAIN_GAPS].local.param = vb->last_int(CHAIN_GAPS); 
}

// ingest one alignment
static inline void chain_piz_filter_ingest_alignmet (VBlock *vb)
{
    mutex_lock (chain_mutex);

    buf_alloc_more (vb, &chain, 1, 0, ChainAlignment, 2, "chain"); 

    int64_t size = vb->last_int(CHAIN_SIZE);

    NEXTENT (ChainAlignment, chain) = (ChainAlignment){ 
        .src_chrom      = vb->last_int(CHAIN_NAMESRC),
        .src_first_1pos = 1 + next_src_0pos,  // +1 bc our alignments are 1-based vs the chain file that is 0-based
        .src_last_1pos  = 1 + next_src_0pos + size - 1,

        .dst_chrom      = vb->last_int(CHAIN_NAMEDST),
        .dst_first_1pos = 1 + next_dst_0pos,
        .dst_last_1pos  = 1 + next_dst_0pos + size - 1,
    };

    Context *gaps_ctx   = &vb->contexts[CHAIN_GAPS];
    next_src_0pos += size + gaps_ctx->local.param;  // src_gap
    next_dst_0pos += size + gaps_ctx->last_value.i; // dst_gap
    
    mutex_unlock (chain_mutex);

    vb->dont_show_curr_line = true;
}

// verify that adding up all alignments and gaps, results in the end position specified in the header
static inline void chain_piz_filter_verify_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);

    ASSINP (next_src_0pos == vb->last_int(CHAIN_ENDSRC),
            "Bad data in chain file %s: Expecting alignments to add up to ENDSRC=%"PRId64", but they add up to %"PRId64":\n%*s",
            z_name, vb->last_int(CHAIN_ENDSRC), next_src_0pos, 
            (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    ASSINP (next_dst_0pos == vb->last_int(CHAIN_ENDDST),
            "Bad data in chain file %s: Expecting alignments to add up to ENDDST=%"PRId64", but they add up to %"PRId64":\n%*s",
            z_name, vb->last_int(CHAIN_ENDDST), next_dst_0pos, 
            (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    mutex_unlock (chain_mutex);
}

CONTAINER_FILTER_FUNC (chain_piz_filter)
{
    // ingest alignments (and report Chain file data issues) only when consuming with --chain, not when merely pizzing
    if (!flag.reading_chain) goto done; 

    // before alignment-set first EOL and before alignments - initialize next_dst_0pos and next_src_0pos
    if (dict_id.num == dict_id_fields[CHAIN_TOPLEVEL] && item == 13) 
        chain_piz_filter_init_alignment_set (vb);

    // save src_gap before reconstructing dst_gap (src_gap was just processed)
    else if (dict_id.num == dict_id_fields[CHAIN_SET] && item == 2) 
        chain_piz_filter_save_src_gap (vb);

    // before EOF of each alignment, ingest alignment
    else if (dict_id.num == dict_id_fields[CHAIN_SET] && item == 3) 
        chain_piz_filter_ingest_alignmet (vb);

    // before alignment-set second EOL and after alignments - verify that numbers add up
    else if (dict_id.num == dict_id_fields[CHAIN_TOPLEVEL] && item == 15) 
        chain_piz_filter_verify_alignment_set (vb);

done:    
    return true;    
}

// sort the alignmants by (dst_contig_index, qstart)
static int chain_sorter (const void *a_, const void *b_)
{   
    ChainAlignment *a = (ChainAlignment *)a_, *b = (ChainAlignment *)b_;

    if (a->src_chrom != b->src_chrom) 
        return a->src_chrom - b->src_chrom;
    else 
        return a->src_first_1pos - b->src_first_1pos;
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

    TEMP_FLAG (quiet, true); // don't show progress indicator for the chain file - it is very fast 
    piz_one_file (0, false, false);
    RESTORE_FLAG (quiet);

    // copy src_contig and dst_contig dictionaries before z_file of the chain file is freed
    liftover_copy_data_from_chain_file();

    // sort the alignmants by (dst_contig_index, qstart)
    qsort (chain.data, chain.len, sizeof (ChainAlignment), chain_sorter);
        
    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    if (flag.show_chain) {
        chain_display_alignments();
        exit_ok;
    }

    chain_filename = file_make_unix_filename (z_name); // full-path unix-style filename, allocates memory

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtfile_genozip_to_txt_header.
    
    flag.reading_chain = NULL;
}

// get dst_contig, src_pos from src_contig, dst_pos (binary search on chain)
static LiftOverStatus chain_get_liftover_coords_do (WordIndex src_contig_index, PosType src_1pos, 
                                                    int32_t start, int32_t end,
                                                    WordIndex *dst_contig_index, PosType *dst_1pos) // out
{
    // case: no mapping
    if (end < start) {
        *dst_contig_index = WORD_INDEX_NONE;
        *dst_1pos = 0;
        return LO_NO_MAPPING;
    }

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = ENT (ChainAlignment, chain, mid);

    // case: success
    if (aln->src_chrom     == src_contig_index && 
        aln->src_first_1pos <= src_1pos && 
        aln->src_last_1pos >= src_1pos) {
            *dst_contig_index = aln->dst_chrom;
            *dst_1pos = aln->dst_first_1pos + (src_1pos - aln->src_first_1pos);
            return LO_OK;
        }

    // case: src_contig, dst_pos is less than aln - search lower half
    if (src_contig_index < aln->src_chrom ||
        (src_contig_index == aln->src_chrom && aln->src_first_1pos > src_1pos))
        return chain_get_liftover_coords_do (src_contig_index, src_1pos, start, mid-1, dst_contig_index, dst_1pos);
    
    // case: src_contig,dst_pos is more than aln - search upper half
    else
        return chain_get_liftover_coords_do (src_contig_index, src_1pos, mid+1, end, dst_contig_index, dst_1pos);
}

LiftOverStatus chain_get_liftover_coords (WordIndex src_contig_index, PosType src_1pos, 
                                          WordIndex *dst_contig_index, PosType *dst_1pos) // out
{
    return chain_get_liftover_coords_do (src_contig_index, src_1pos, 0, chain.len-1, dst_contig_index, dst_1pos);
}
