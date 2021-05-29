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
// Chain: Query  Genozip: Dst, Luft (i.e. lifted-over, "luft" being an alternative past tense of "lift")

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "reference.h"
#include "endianness.h"
#include "strings.h"
#include "piz.h"
#include "dict_id.h"
#include "stats.h"
#include "codec.h"
#include "version.h"
#include "chain.h"
#include "reconstruct.h"
#include "zfile.h"

static WordIndex num_ref_contigs=0; // ZIP of a chain file

typedef struct {
    WordIndex src_chrom; // index in CHAIN_NAMESRC
    PosType src_first_1pos, src_last_1pos;

    WordIndex dst_chrom; // index in CHAIN_NAMEDST - the nodes+dicts copied to the genozipped file (eg as VCF_oCHROM)
    PosType dst_first_1pos, dst_last_1pos;

    bool xstrand;
} ChainAlignment;

static Buffer chain = EMPTY_BUFFER;   // immutable after loaded
static Mutex chain_mutex = {};        // protect chain wnile loading
static PosType next_dst_0pos=0, next_src_0pos=0; // used for loading chain

// 1) When pizzing a chain as a result of genozip --chain: chain_piz_initialize: src_contig and dst_contig are generated 
//    from z_file->contexts of a chain file 
// 2) Contig lengths are added when pizzing each alignment set, in chain_piz_filter_add_contig_length
// 3) At the beginning of ZIP of a file with --chain: dst_contig copied to z_file->contexts of the file being compressed
// 4) vcf_lo_get_liftover_coords src_contig is consulted to get map the index of the primary chrom from the node_index
//    of the file being compressed to the word_index of src_contig in the chain data
static Buffer src_contig_dict = EMPTY_BUFFER;
static Buffer src_contigs     = EMPTY_BUFFER;
static Buffer dst_contig_dict = EMPTY_BUFFER;
static Buffer dst_contigs     = EMPTY_BUFFER;

char *chain_filename = NULL; // global - chain filename

//-----------------------
// TXTFILE stuff
//-----------------------

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t chain_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */)
{
    ARRAY (char, data, vb->txt_data)

    ASSERT (*i >= 0 && *i < data_len, "*i=%d is out of range [0,%"PRIu64"]", *i, data_len);

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

void chain_zip_initialize (void)
{
    ASSINP0 (flag.reference == REF_EXTERNAL, "Please specify the destination reference file with --reference");
    flag.reference = REF_MAKE_CHAIN; // since we're not attempting to use the data from the reference to compress the chain file itself

    if (z_file->num_txt_components_so_far == 1) { // run only for 1st component in a bound file

        ConstBufferP contig_dict, contigs;
        ref_contigs_get (&contig_dict, &contigs);

        num_ref_contigs = (WordIndex)contigs->len;
        ctx_build_zf_ctx_from_contigs (CHAIN_NAMEDST, contigs, contig_dict);
    }
}

void chain_seg_initialize (VBlock *vb)
{
    vb->contexts[CHAIN_NAMEDST] .no_stons  =       
    vb->contexts[CHAIN_NAMESRC] .no_stons  =       // (these 2 ^) need b250 node_index for reference
    vb->contexts[CHAIN_TOPLEVEL].no_stons  = 
    vb->contexts[CHAIN_CHAIN]   .no_stons  = 
    vb->contexts[CHAIN_SEP]     .no_stons  = 
    vb->contexts[CHAIN_STRNDSRC].no_stons  = true; // (these 4 ^) keep in b250 so it can be eliminated as all_the_same

    vb->contexts[CHAIN_STARTSRC].flags.store = 
    vb->contexts[CHAIN_ENDSRC]  .flags.store = 
    vb->contexts[CHAIN_SIZESRC] .flags.store = 
    vb->contexts[CHAIN_STARTDST].flags.store = 
    vb->contexts[CHAIN_ENDDST]  .flags.store = 
    vb->contexts[CHAIN_SIZEDST] .flags.store = 
    vb->contexts[CHAIN_SIZE]    .flags.store = 
    vb->contexts[CHAIN_GAPS]    .flags.store = STORE_INT;

    vb->contexts[CHAIN_NAMESRC] .flags.store =
    vb->contexts[CHAIN_NAMEDST] .flags.store = STORE_INDEX;
}

void chain_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - chain_piz_filter needs to changed and support backward compatability
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 16,
        .items        = { { .dict_id = (DictId)dict_id_fields[CHAIN_CHAIN],    .seperator = {' '} },  // 0
                          { .dict_id = (DictId)dict_id_fields[CHAIN_SCORE],    .seperator = {' '} },  // 1
                          { .dict_id = (DictId)dict_id_fields[CHAIN_NAMESRC],  .seperator = {' '} },  // 2
                          { .dict_id = (DictId)dict_id_fields[CHAIN_SIZESRC],  .seperator = {' '} },  // 3
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STRNDSRC], .seperator = {' '} },  // 4
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STARTSRC], .seperator = {' '} },  // 5
                          { .dict_id = (DictId)dict_id_fields[CHAIN_ENDSRC],   .seperator = {' '} },  // 6
                          { .dict_id = (DictId)dict_id_fields[CHAIN_NAMEDST],  .seperator = {' '} },  // 7
                          { .dict_id = (DictId)dict_id_fields[CHAIN_SIZEDST],  .seperator = {' '} },  // 8
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STRNDDST], .seperator = {' '} },  // 9
                          { .dict_id = (DictId)dict_id_fields[CHAIN_STARTDST], .seperator = {' '} },  // 10
                          { .dict_id = (DictId)dict_id_fields[CHAIN_ENDDST],   .seperator = {' '} },  // 11
                          { .dict_id = (DictId)dict_id_fields[CHAIN_ID]                           },  // 12
                          { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]                          },  // 13
                          { .dict_id = (DictId)dict_id_fields[CHAIN_SET]                          },  // 14
                          { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]                          } },// 15
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
        ctx->last_value.i          - vb->last_int(CHAIN_STARTDST)) {
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_ENDDST }), 2, ctx, field_len + 1);        
        ctx->numeric_only = false;
    }

    else
        seg_pos_field (vb, CHAIN_ENDDST, CHAIN_STARTDST, false, false, 0, field_start, field_len, 0, field_len + 1);
}

static void chain_seg_verify_size_dst (VBlock *vb, WordIndex ref_index, 
                                       const char *name_dst, unsigned name_dst_len,
                                       const char *size_dst_start, unsigned size_dst_len)
{
    PosType size_according_to_chain, size_according_to_ref;

    ASSSEG (str_get_int_range64 (size_dst_start, size_dst_len, 1, 1000000000000000, &size_according_to_chain),
            size_dst_start, "Invalid size field of \"%.*s\": %.*s", name_dst_len, name_dst, size_dst_len, size_dst_start);

    size_according_to_ref = ref_contigs_get_contig_length (ref_index, 0, 0, true);

    ASSERT (size_according_to_chain == size_according_to_ref, 
            "Size of \"%.*s\" in chain file is %"PRId64", but in the reference file it is %"PRId64,
             name_dst_len, name_dst, size_according_to_chain, size_according_to_ref);
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
    SEG_NEXT_ITEM_SP (CHAIN_SIZESRC);
    SEG_NEXT_ITEM_SP (CHAIN_STRNDSRC);

    GET_NEXT_ITEM_SP (CHAIN_STARTSRC);
    seg_pos_field (vb, CHAIN_STARTSRC, CHAIN_ENDSRC, false, false, 0, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_ENDSRC);
    seg_pos_field (vb, CHAIN_ENDSRC, CHAIN_STARTSRC, false, false, 0, field_start, field_len, 0,field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_NAMEDST);

    bool is_new;
    WordIndex name_dst_index = seg_by_did_i_ex (vb, CHAIN_NAMEDST_str, CHAIN_NAMEDST_len, CHAIN_NAMEDST, CHAIN_NAMEDST_len+1, &is_new);
    
    // case: chrom name doesn't match any ref contig. if it is an alternative name (eg 1->chr1), then we add
    // it to the CHAIN_NAMEDST dict (in addition to all the ref contigs that were copied in chain_zip_initialize),
    // and a mapping will be created by ref_alt_chroms_compress. 
    WordIndex ref_index = name_dst_index;
    if (name_dst_index >= num_ref_contigs) {
        static bool once = false;
        ref_index = ref_alt_chroms_zip_get_alt_index (CHAIN_NAMEDST_str, CHAIN_NAMEDST_len, WI_REF_CONTIG, WORD_INDEX_NONE);
        if (ref_index == WORD_INDEX_NONE && is_new) {
            if (!once)
                WARN ("\"%.*s\": Contig appears in chain file %s but is missing in reference file %s. When %s is used to liftover files, lines with this contig will not be lifted over.", 
                      CHAIN_NAMEDST_len, CHAIN_NAMEDST_str, txt_name, ref_filename, z_name);
            else
                WARN ("\"%.*s\": Contig appears in chain file but is missing in reference file", CHAIN_NAMEDST_len, CHAIN_NAMEDST_str);
            once = true;
        }
    }

    SEG_NEXT_ITEM_SP (CHAIN_SIZEDST);
    if (ref_index != WORD_INDEX_NONE)
        chain_seg_verify_size_dst (vb, ref_index, CHAIN_NAMEDST_str, CHAIN_NAMEDST_len, CHAIN_SIZEDST_str, CHAIN_SIZEDST_len);

    SEG_NEXT_ITEM_SP (CHAIN_STRNDDST);

    GET_NEXT_ITEM_SP (CHAIN_STARTDST);
    seg_pos_field (vb, CHAIN_STARTDST, CHAIN_ENDDST, false, false, 0, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_ENDDST);
    chain_seg_dst_end_field (vb, field_start, field_len);

    // if ID is numeric, preferably store as delta, if not - normal snip
    GET_LAST_ITEM_SP (CHAIN_ID);
    if (str_is_int (field_start, field_len))
        seg_pos_field (vb, CHAIN_ID, CHAIN_ID, false, false, 0, field_start, field_len, 0, field_len + 1); // just a numeric delta    
    else
        seg_by_did_i (vb, field_start, field_len, CHAIN_ID, field_len + 1);

    SEG_EOL (CHAIN_EOL, true);

    // set containing 1 or more alignments of SIZE, DDST, SRC, last alignment in set is only SIZE
    uint32_t num_alignments=0;
    bool is_last_alignment=false;
    char last_gap_sep = ' ';
    do {
        GET_MAYBE_LAST_ITEM_SP (CHAIN_SIZE);
        chain_seg_size_field (vb, field_start, field_len);

        is_last_alignment = (separator == '\n');

        if (!is_last_alignment) {
            last_gap_sep = separator;

            // note: we store both DIFFs in the same context as they are correlated (usually one of them is 0)
            seg_by_did_i (vb, &separator, 1, CHAIN_SEP, 0); // Chain format may have space or tab as a separator
            { SEG_NEXT_ITEM_SP (CHAIN_GAPS); }

            seg_by_did_i (vb, &separator, 1, CHAIN_SEP, 0); // space or tab
            { SEG_LAST_ITEM_SP (CHAIN_GAPS); }
        
        }
        // last line doesn't have DDST and DSRC - seg the "standard" separator and then delete it 
        // (this way we don't ruin the all_the_same of CHAIN_SEP)
        else {
            seg_by_did_i (vb, &last_gap_sep, 1, CHAIN_SEP, 0); 
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // src_gap

            seg_by_did_i (vb, &last_gap_sep, 1, CHAIN_SEP, 0); 
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // dst_gap
        }
        SEG_EOL (CHAIN_EOL, false);

        num_alignments++;
    } 
    while (!is_last_alignment);

    // alignment set - IMPORTNAT - if changing fields, update chain_piz_filter
    SmallContainer alignment_set = { 
        .repeats             = num_alignments,
        .nitems_lo           = 6,
        .keep_empty_item_sep = true, // avoid double-deleting the space - only chain_piz_special_BACKSPACE should delete it, not container_reconstruct_do
        .filter_items        = true,
        .items               = { { .dict_id = (DictId)dict_id_fields[CHAIN_SIZE] },
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_SEP]  }, 
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_GAPS] }, // src_gap
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_SEP]  }, 
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_GAPS] }, // dst_gap
                                 { .dict_id = (DictId)dict_id_fields[CHAIN_EOL]  } }
    };

    container_seg_by_ctx (vb, &vb->contexts[CHAIN_SET], (ContainerP)&alignment_set, 0, 0, 0);

    // Empty line after alignment set
    GET_LAST_ITEM (EmptyLine);
    ASSSEG0 (!field_len, field_start, "Expecting an empty line after alignment set");
    SEG_EOL (CHAIN_EOL, false); 

    return next_field;
}

//--------------
// PIZ functions
//--------------

static void chain_piz_initialize_copy_contigs (Buffer *contigs, Buffer *dict, DidIType name_did_i, bool get_LN_from_ref)
{
    // copy dictionary buffer
    buf_copy (evb, dict, &z_file->contexts[name_did_i].dict, char, 0, 0, "chain_contig_dict");

    ARRAY (CtxWord, contig_words, z_file->contexts[name_did_i].word_list);
    dict->param = contig_words_len; // dict param contains number of words

    // copy word_list into a contigs buffer - for convenience of using ctx_build_zf_ctx_from_contigs later
    buf_alloc (evb, contigs, 0, contig_words_len, RefContig, 0, "chain_contigs");
    contigs->len = contig_words_len;
    ARRAY (RefContig, ctg, *contigs);
    
    for (uint32_t i=0; i < contig_words_len; i++) {
        PosType LN=0;

        if (get_LN_from_ref) {
            const char *contig = ctx_get_snip_by_word_index (&z_file->contexts[name_did_i], i, 0, 0);
            LN = ref_contigs_get_contig_length (i, contig, contig_words[i].snip_len, false);
        }
        
        ctg[i] = (RefContig) { .char_index = contig_words[i].char_index,
                               .snip_len   = contig_words[i].snip_len,
                               .LN         = LN };
    }
}

// called after reconstructing the txt header and before compute threads
bool chain_piz_initialize (void)
{
    mutex_initialize (chain_mutex);

    if (flag.reading_chain) {
        chain_piz_initialize_copy_contigs (&src_contigs, &src_contig_dict, CHAIN_NAMESRC, false);
        chain_piz_initialize_copy_contigs (&dst_contigs, &dst_contig_dict, CHAIN_NAMEDST, true);
    }

    return true; // proceed with PIZ
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_BACKSPACE)
{
    ASSERT0 (vb->txt_data.len, "vb->txt_data.len is 0");
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

        iprintf ("Primary: %s %"PRId64"-%"PRId64" Luft: %s %"PRId64"-%"PRId64"\n",
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

    buf_alloc (NULL, &chain, 1, 0, ChainAlignment, 2, "chain");  // initial allocation is in chain_load

    int64_t size = vb->last_int(CHAIN_SIZE);

    NEXTENT (ChainAlignment, chain) = (ChainAlignment){ 
        .src_chrom      = vb->last_index(CHAIN_NAMESRC),
        .src_first_1pos = 1 + next_src_0pos,  // +1 bc our alignments are 1-based vs the chain file that is 0-based
        .src_last_1pos  = 1 + next_src_0pos + size - 1,

        // use alt chroms as they appear in the reference file, rather than their name in the chain file
        .dst_chrom      = ref_alt_get_final_index (vb->last_index(CHAIN_NAMEDST)), 
        .dst_first_1pos = 1 + next_dst_0pos,
        .dst_last_1pos  = 1 + next_dst_0pos + size - 1,

        .xstrand  = *last_txt(vb, CHAIN_STRNDSRC) != *last_txt(vb, CHAIN_STRNDDST)
    };

    Context *gaps_ctx   = &vb->contexts[CHAIN_GAPS];
    next_src_0pos += size + gaps_ctx->local.param;  // src_gap
    next_dst_0pos += size + gaps_ctx->last_value.i; // dst_gap
    
    mutex_unlock (chain_mutex);

    vb->drop_curr_line = "chain";
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

// set contig lengths according to tSize and qSize
static inline void chain_piz_filter_add_contig_length (VBlock *vb)
{
    mutex_lock (chain_mutex);

    RefContig *src_rc = ENT (RefContig, src_contigs, vb->last_index(CHAIN_NAMESRC));
    RefContig *dst_rc = ENT (RefContig, dst_contigs, vb->last_index(CHAIN_NAMEDST));

    src_rc->LN = src_rc->max_pos = vb->last_int(CHAIN_SIZESRC);
    dst_rc->LN = dst_rc->max_pos = vb->last_int(CHAIN_SIZEDST);

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
    else if (dict_id.num == dict_id_fields[CHAIN_SET] && item == 4) 
        chain_piz_filter_save_src_gap (vb);

    // before EOF of each alignment, ingest alignment
    else if (dict_id.num == dict_id_fields[CHAIN_SET] && item == 5) 
        chain_piz_filter_ingest_alignmet (vb);

    // before alignment-set second EOL and after alignments - verify that numbers add up, and also set contigs
    else if (dict_id.num == dict_id_fields[CHAIN_TOPLEVEL] && item == 15) {
        chain_piz_filter_verify_alignment_set (vb);
        chain_piz_filter_add_contig_length (vb);
    }

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

// main thread: ZIP with --chain, reading txt header: returns null-terminated string of contig, or NULL if contig_i is out of range
const char *chain_get_dst_contig (uint32_t contig_i, PosType *length)
{
    if (contig_i >= dst_contigs.len) return NULL;

    const RefContig *ctg = ENT (RefContig, dst_contigs, contig_i);

    *length = ctg->LN;
    return ENT (char, dst_contig_dict, ctg->char_index);
}

// returns the contig index in src_contigs, word WORD_INDEX_NONE if there is none
WordIndex chain_get_src_contig_index (const char *contig, unsigned contig_len)
{
    WordIndex ctg_i;
    ARRAY (RefContig, src_ctgs, src_contigs);
    
    for (ctg_i=0; ctg_i < src_ctgs_len; ctg_i++) {
        const char *src_contig = ENT (char, src_contig_dict, src_ctgs[ctg_i].char_index);

        if (src_ctgs[ctg_i].snip_len == contig_len && 
            !memcmp (contig, src_contig, contig_len)) break; // found
    }

    return (ctg_i == src_ctgs_len) ? WORD_INDEX_NONE : ctg_i;
}

// called by vcf_lo_zip_initialize when zipping with --chain to copy chain->dst_contigs to z_file->ochrom 
void chain_copy_dst_contigs_to_z_file (DidIType dst_contig_did_i)
{
    ctx_build_zf_ctx_from_contigs (dst_contig_did_i, &dst_contigs, &dst_contig_dict);
}

static void chain_contigs_show (const Buffer *contigs, const Buffer *dict, const char *type)
{
    ARRAY (const RefContig, cn, *contigs);

    iprintf ("\n%s contigs calculated from the chain file alignments:\n", type);
    for (uint32_t i=0; i < cn_len; i++) {

        const char *chrom_name = ENT (const char, *dict, cn[i].char_index);

        if (cn[i].snip_len)
            iprintf ("i=%u '%s' LN=%"PRId64" chrom_index=%d char_index=%"PRIu64" snip_len=%u\n",
                     i, chrom_name, cn[i].LN, cn[i].chrom_index, cn[i].char_index, cn[i].snip_len);
        else
            iprintf ("i=%u chrom_index=%d (unused - not present in txt data)\n", i, cn[i].chrom_index);
    }
}

// zip: load chain file as a result of genozip --chain
void chain_load (void)
{
    ASSERTNOTNULL (flag.reading_chain);
    SAVE_FLAGS_AUX;

    buf_alloc (evb, &chain, 0, 1000, ChainAlignment, 1, "chain"); // must be allocated by main thread
    
    z_file = file_open (flag.reading_chain, READ, Z_FILE, DT_CHAIN);    
    zfile_read_genozip_header (0);
    
    ASSINP (z_file->data_type == DT_CHAIN, "expected %s to be a genozip'ed chain file, but its a %s file. Tip: compress the chain with \"genozip --input chain\"", 
            z_name, dt_name (z_file->data_type));

    z_file->basename = file_basename (flag.reading_chain, false, "(chain-file)", NULL, 0);

    TEMP_VALUE (command, PIZ);

    if (flag.reference) flag.reference = REF_LIFTOVER;

    flag.no_writer = true;
    flag.genocat_no_reconstruct = false;
    flag.genocat_global_area_only = false;
    
    Dispatcher dispachter = piz_z_file_initialize (false);

    // load the reference, now that it is set (either explicitly from the command line, or implicitly from the chain GENOZIP_HEADER)
    SAVE_VALUE (z_file); // actually, read the reference first
    ref_load_external_reference (false);
    RESTORE_VALUE (z_file);

    flag.quiet = true; // don't show progress indicator for the chain file - it is very fast 
    flag.maybe_vb_modified_by_reconstructor = true; // we drop all the lines

    piz_one_txt_file (dispachter, false);

    // --show-chain-contigs
    if (flag.show_chain_contigs) {
        chain_contigs_show (&src_contigs, &src_contig_dict, "SRC");
        chain_contigs_show (&dst_contigs, &dst_contig_dict, "DST");
        if (exe_type == EXE_GENOCAT) exit_ok;  // in genocat this, not the data
    }

    // sort the alignmants by (dst_contig_index, qstart)
    qsort (chain.data, chain.len, sizeof (ChainAlignment), chain_sorter);
        
    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    flag.reference = REF_LIFTOVER; // will be needed by the files we are about to compress with --chain

    if (flag.show_chain) {
        chain_display_alignments();
        exit_ok;
    }

    chain_filename = file_make_unix_filename (z_name); // full-path unix-style filename, allocates memory

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtheader_piz_read_and_reconstruct.
    
    flag.reading_chain = NULL;
}

// reset to factory defaults
void chain_destroy (void)
{
    buf_destroy (&chain);
    buf_destroy (&src_contig_dict);
    buf_destroy (&src_contigs);
    buf_destroy (&dst_contig_dict);
    buf_destroy (&dst_contigs);
    mutex_destroy (chain_mutex);
    next_dst_0pos = next_src_0pos=0;
    chain_filename = NULL; 
}

// binary search for the first alignment of a src contig
static ChainAlignment *chain_get_first_aln (WordIndex src_contig_index, int32_t start, int32_t end) 
{
    // We know that this src_contig exists, so its unexpected if it doesn't
    ASSERT (end >= start, "Unable to find src_contig_index=%d in chain", src_contig_index);

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = ENT (ChainAlignment, chain, mid);

    // case: success - same src_contig, and previous aln has a different src_contig 
    if (aln->src_chrom == src_contig_index && (!mid || (aln-1)->src_chrom != src_contig_index))
        return aln;

    // case: src_contig less than aln's, OR is the same, but aln isn't the first with this src_contig - search lower half
    if (src_contig_index <= aln->src_chrom)
        return chain_get_first_aln (src_contig_index, start, mid-1);
    
    // case: src_contig is more than aln's - search upper half
    else
        return chain_get_first_aln (src_contig_index, mid+1, end);
}

// append dst_contigs buffer with all dst chroms indicies for which there exists alignment with src_chrom
void chain_append_all_dst_contig_index (const char *src_contig_name, unsigned src_contig_name_len, Buffer *dst_contigs)
{
    WordIndex src_contig_index = chain_get_src_contig_index (src_contig_name, src_contig_name_len);
    WordIndex prev_dst = WORD_INDEX_NONE;

    if (src_contig_index == WORD_INDEX_NONE) return; // this contig is not in the chain file

    for (ChainAlignment *aln = chain_get_first_aln (src_contig_index, 0, chain.len-1);
         aln < LASTENT (ChainAlignment, chain) && aln->src_chrom == src_contig_index;
         aln++)
         
         if (aln->dst_chrom != prev_dst) { // note: this prevents consecutive duplicates, but not non-consecutive duplicates
            buf_alloc (NULL, dst_contigs, 100, 1, WordIndex, 2, NULL);
            NEXTENT (WordIndex, *dst_contigs) = prev_dst = aln->dst_chrom;
         }
}

// get dst_contig, src_pos from src_contig, dst_pos (binary search on chain)
// returns true if successful
static bool chain_get_liftover_coords_do (WordIndex src_contig_index, PosType src_1pos, 
                                          int32_t start, int32_t end,
                                          WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i) // out
{
    // case: no mapping
    if (end < start) {
        if (dst_contig_index) *dst_contig_index = WORD_INDEX_NONE;
        if (dst_1pos) *dst_1pos = 0;
        if (xstrand) *xstrand = false;
        return false;
    }

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = ENT (ChainAlignment, chain, mid);

    // case: success
    if (aln->src_chrom == src_contig_index && 
        aln->src_first_1pos <= src_1pos && 
        aln->src_last_1pos >= src_1pos) {
            if (dst_contig_index) *dst_contig_index = aln->dst_chrom;
            if (dst_1pos) *dst_1pos = aln->dst_first_1pos + (src_1pos - aln->src_first_1pos);
            if (xstrand) *xstrand  = aln->xstrand;
            if (aln_i) *aln_i = mid;
            return true;
        }

    // case: src_contig, dst_pos is less than aln - search lower half
    if (src_contig_index < aln->src_chrom ||
        (src_contig_index == aln->src_chrom && src_1pos < aln->src_first_1pos))
        return chain_get_liftover_coords_do (src_contig_index, src_1pos, start, mid-1, dst_contig_index, dst_1pos, xstrand, aln_i);
    
    // case: src_contig,dst_pos is more than aln - search upper half
    else
        return chain_get_liftover_coords_do (src_contig_index, src_1pos, mid+1, end, dst_contig_index, dst_1pos, xstrand, aln_i);
}

bool chain_get_liftover_coords (WordIndex src_contig_index, PosType src_1pos, 
                                WordIndex *dst_contig_index, PosType *dst_1pos, bool *xstrand, uint32_t *aln_i) // out
{
    return chain_get_liftover_coords_do (src_contig_index, src_1pos, 0, chain.len-1, dst_contig_index, dst_1pos, xstrand, aln_i);
}
