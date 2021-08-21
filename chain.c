// ------------------------------------------------------------------
//   chain.c
//   Copyright (C) 2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt
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

typedef struct {
    WordIndex prim_chrom; // index in CHAIN_NAMEPRIM
    PosType prim_first_1pos, prim_last_1pos;

    WordIndex luft_chrom; // index in CHAIN_NAMELUFT - the nodes+dicts copied to the genozipped file (eg as VCF_oCHROM)
    PosType luft_first_1pos, luft_last_1pos;

    bool is_xstrand;
} ChainAlignment;

static Buffer chain = EMPTY_BUFFER;   // immutable after loaded
static Mutex chain_mutex = {};        // protect chain wnile loading
static PosType next_luft_0pos=0, next_prim_0pos=0; // used for loading chain

// 1) When pizzing a chain as a result of genozip --chain: chain_piz_initialize: prim_contig and luft_contig are generated 
//    from z_file->contexts of a chain file 
// 2) Contig lengths are added when pizzing each alignment set, in chain_piz_filter_add_contig_length
// 3) At the beginning of ZIP of a file with --chain: luft_contig copied to z_file->contexts of the file being compressed
// 4) vcf_lo_get_liftover_coords prim_contig is consulted to get map the index of the primary chrom from the node_index
//    of the file being compressed to the word_index of prim_contig in the chain data
static Buffer prim_contig_dict   = EMPTY_BUFFER;
static Buffer prim_contig_counts = EMPTY_BUFFER;
static Buffer prim_contigs       = EMPTY_BUFFER;
static Buffer luft_contig_dict   = EMPTY_BUFFER;
static Buffer luft_contigs       = EMPTY_BUFFER;

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
    ASSINP0 (flag.reference == REF_EXTERNAL || flag.force, "Please specify the destination reference file with --reference or use --force to override this");

    if (flag.reference != REF_EXTERNAL) return;

    flag.reference = REF_MAKE_CHAIN; // since we're not attempting to use the data from the reference to compress the chain file itself

    if (z_file->num_txt_components_so_far == 1) { // run only for 1st component in a bound file

        ConstBufferP contig_dict, contigs;

        ref_contigs_get (prim_ref, &contig_dict, &contigs);
        ctx_build_zf_ctx_from_contigs (CHAIN_NAMEPRIM, contigs, contig_dict);

        ref_contigs_get (gref, &contig_dict, &contigs);
        ctx_build_zf_ctx_from_contigs (CHAIN_NAMELUFT, contigs, contig_dict);
    }
}

void chain_seg_initialize (VBlock *vb)
{
    CTX(CHAIN_NAMELUFT)-> no_stons  =       
    CTX(CHAIN_NAMEPRIM)-> no_stons  =       // (these 2 ^) need b250 node_index for reference
    CTX(CHAIN_TOPLEVEL)-> no_stons  = 
    CTX(CHAIN_CHAIN)->    no_stons  = 
    CTX(CHAIN_SEP)->      no_stons  = 
    CTX(CHAIN_VERIFIED)-> no_stons  = 
    CTX(CHAIN_STRNDPRIM)->no_stons  = true; // (these 4 ^) keep in b250 so it can be eliminated as all_the_same

    CTX(CHAIN_STARTPRIM)->flags.store = 
    CTX(CHAIN_ENDPRIM)->  flags.store = 
    CTX(CHAIN_SIZEPRIM)-> flags.store = 
    CTX(CHAIN_STARTLUFT)->flags.store = 
    CTX(CHAIN_ENDLUFT)->  flags.store = 
    CTX(CHAIN_SIZELUFT)-> flags.store = 
    CTX(CHAIN_SIZE)->     flags.store = 
    CTX(CHAIN_VERIFIED)-> flags.store = 
    CTX(CHAIN_GAPS)->     flags.store = STORE_INT;

    CTX(CHAIN_NAMEPRIM)-> flags.store =
    CTX(CHAIN_NAMELUFT)-> flags.store = STORE_INDEX;

    CTX(CHAIN_NAMEPRIM)-> no_vb1_sort =
    CTX(CHAIN_NAMELUFT)-> no_vb1_sort = true;

    CTX(CHAIN_NAMEPRIM)-> counts_section = true;
}

void chain_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - chain_piz_filter needs to changed and support backward compatability
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 16,
        .items        = { { .dict_id = { _CHAIN_CHAIN },     .seperator = {' '} },  // 0
                          { .dict_id = { _CHAIN_SCORE },     .seperator = {' '} },  // 1
                          { .dict_id = { _CHAIN_NAMEPRIM },  .seperator = {' '} },  // 2
                          { .dict_id = { _CHAIN_SIZEPRIM },  .seperator = {' '} },  // 3
                          { .dict_id = { _CHAIN_STRNDPRIM }, .seperator = {' '} },  // 4
                          { .dict_id = { _CHAIN_STARTPRIM }, .seperator = {' '} },  // 5
                          { .dict_id = { _CHAIN_ENDPRIM },   .seperator = {' '} },  // 6
                          { .dict_id = { _CHAIN_NAMELUFT },  .seperator = {' '} },  // 7
                          { .dict_id = { _CHAIN_SIZELUFT },  .seperator = {' '} },  // 8
                          { .dict_id = { _CHAIN_STRNDLUFT }, .seperator = {' '} },  // 9
                          { .dict_id = { _CHAIN_STARTLUFT }, .seperator = {' '} },  // 10
                          { .dict_id = { _CHAIN_ENDLUFT },   .seperator = {' '} },  // 11
                          { .dict_id = { _CHAIN_ID }                            },  // 12
                          { .dict_id = { _CHAIN_EOL }                           },  // 13
                          { .dict_id = { _CHAIN_SET }                           },  // 14
                          { .dict_id = { _CHAIN_EOL }                           } },// 15
    };

    container_seg (vb, CTX(CHAIN_TOPLEVEL), (ContainerP)&top_level, 0, 0, 0);
}

bool chain_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return true; // contexts are expected to have small dictionaries
}

static void chain_seg_size_field (VBlock *vb, const char *field_start, int field_len)
{
    int64_t size;
    ASSSEG (str_get_int (field_start, field_len, &size), field_start, "Expecting size to be an integer, but found %*s", field_len, field_start);

    Context *ctx = CTX(CHAIN_SIZE);

    if (vb->last_int(CHAIN_ENDPRIM) - vb->last_int(CHAIN_STARTPRIM) == size) { // happens if the alignment set has only one alignment
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_SIZE }), 2, ctx, field_len + 1);
        ctx->numeric_only = false;
    }

    else
        seg_integer_or_not (vb, ctx, field_start, field_len, field_len + 1);
}

static void chain_seg_luft_end_field (VBlock *vb, const char *field_start, int field_len)
{
    Context *ctx = CTX(CHAIN_ENDLUFT);

    ASSSEG (str_get_int (field_start, field_len, &ctx->last_value.i), field_start, "Expecting qEnd to be an integer, but found %*s", field_len, field_start);

    if (vb->last_int(CHAIN_ENDPRIM) - vb->last_int(CHAIN_STARTPRIM) ==
        ctx->last_value.i          - vb->last_int(CHAIN_STARTLUFT)) {
        seg_by_ctx (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_ENDLUFT }), 2, ctx, field_len + 1);        
        ctx->numeric_only = false;
    }

    else
        seg_pos_field (vb, CHAIN_ENDLUFT, CHAIN_STARTLUFT, 0, 0, field_start, field_len, 0, field_len + 1);
}

static bool chain_seg_verify_contig (VBlock *vb, Reference ref, WordIndex name_index, bool is_new, bool *once,
                                     const char *name, unsigned name_len,
                                     const char *size_str, unsigned size_len,
                                     const char *last_name, unsigned last_name_len)
{
    // case: chrom name doesn't match any ref contig. if it is an alternative name (eg 1->chr1), then we add
    // it to the CHAIN_NAMELUFT dict (in addition to all the ref contigs that were copied in chain_zip_initialize),
    // and a mapping will be created by ref_alt_chroms_compress. 
    WordIndex ref_index = name_index;
    if (name_index >= ref_get_contigs (ref)->len) {
        ref_index = ref_alt_chroms_get_alt_index (ref, name, name_len, atoi (size_str), WORD_INDEX_NONE);
        if (ref_index == WORD_INDEX_NONE && is_new) {
            if (! *once)
                WARN ("\"%.*s\": Contig appears in chain file %s but is missing in reference file %s %s. When %s is used to liftover files, lines with this contig will not be lifted over.", 
                      name_len, name, txt_name, ref==gref ? "LUFT" : "PRIMARY", ref_get_filename (ref), z_name);
            else
                WARN ("\"%.*s\": Contig appears in chain file but is missing in %s reference file", name_len, name, ref==gref ? "LUFT" : "PRIMARY");
            *once = true;
        }
    }

    if (ref_index == WORD_INDEX_NONE) return false;

    PosType size_according_to_chain, size_according_to_ref;

    ASSSEG (str_get_int_range64 (size_str, size_len, 1, 1000000000000000, &size_according_to_chain),
            size_str, "Invalid size field of \"%.*s\": %.*s", name_len, name, size_len, size_str);

    size_according_to_ref = ref_contigs_get_contig_length (ref, ref_index, 0, 0, true);

    if (size_according_to_chain != size_according_to_ref) {
        // output error in the first alignment of this contig
        ASSERTW (str_issame (name, last_name),
                "Size of \"%.*s\" in chain file is %"PRId64", but in %s it is %"PRId64". Excluding it.",
                name_len, name, size_according_to_chain, ref_get_filename (ref), size_according_to_ref);
        return false;
    }

    return true; // verified
}

const char *chain_seg_txt_line (VBlock *vb, const char *field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    static bool once_src = false, once_dst = false;
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;
    bool is_new;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    SEG_NEXT_ITEM_SP (CHAIN_CHAIN);
    SEG_NEXT_ITEM_SP (CHAIN_SCORE);
    
    GET_NEXT_ITEM_SP (CHAIN_NAMEPRIM);
    WordIndex name_prim_index = seg_by_did_i_ex (vb, CHAIN_NAMEPRIM_str, CHAIN_NAMEPRIM_len, CHAIN_NAMEPRIM, CHAIN_NAMEPRIM_len+1, &is_new);
    
    SEG_NEXT_ITEM_SP (CHAIN_SIZEPRIM);
    bool verified = chain_seg_verify_contig (vb, prim_ref, name_prim_index, is_new, &once_src, 
                                             CHAIN_NAMEPRIM_str, CHAIN_NAMEPRIM_len, CHAIN_SIZEPRIM_str, CHAIN_SIZEPRIM_len,
                                             last_txt(vb, CHAIN_NAMEPRIM), vb->last_txt_len(CHAIN_NAMEPRIM));

    SEG_NEXT_ITEM_SP (CHAIN_STRNDPRIM);

    GET_NEXT_ITEM_SP (CHAIN_STARTPRIM);
    seg_pos_field (vb, CHAIN_STARTPRIM, CHAIN_ENDPRIM, 0, 0, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_ENDPRIM);
    seg_pos_field (vb, CHAIN_ENDPRIM, CHAIN_STARTPRIM, 0, 0, field_start, field_len, 0,field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_NAMELUFT);
    WordIndex name_luft_index = seg_by_did_i_ex (vb, CHAIN_NAMELUFT_str, CHAIN_NAMELUFT_len, CHAIN_NAMELUFT, CHAIN_NAMELUFT_len+1, &is_new);
    
    SEG_NEXT_ITEM_SP (CHAIN_SIZELUFT);
    verified &= chain_seg_verify_contig (vb, gref, name_luft_index, is_new, &once_dst, 
                                         CHAIN_NAMELUFT_str, CHAIN_NAMELUFT_len, CHAIN_SIZELUFT_str, CHAIN_SIZELUFT_len,
                                         last_txt(vb, CHAIN_NAMELUFT), vb->last_txt_len(CHAIN_NAMELUFT));

    seg_by_did_i (vb, verified ? "1" : "0", 1, CHAIN_VERIFIED, 0);

    seg_set_last_txt (vb, CTX(CHAIN_NAMEPRIM), CHAIN_NAMEPRIM_str, CHAIN_NAMEPRIM_len, STORE_NONE);
    seg_set_last_txt (vb, CTX(CHAIN_NAMELUFT), CHAIN_NAMELUFT_str, CHAIN_NAMELUFT_len, STORE_NONE);

    SEG_NEXT_ITEM_SP (CHAIN_STRNDLUFT);

    GET_NEXT_ITEM_SP (CHAIN_STARTLUFT);
    seg_pos_field (vb, CHAIN_STARTLUFT, CHAIN_ENDLUFT, 0, 0, field_start, field_len, 0, field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_ENDLUFT);
    chain_seg_luft_end_field (vb, field_start, field_len);

    // if ID is numeric, preferably store as delta, if not - normal snip
    GET_LAST_ITEM_SP (CHAIN_ID);
    if (str_is_int (field_start, field_len))
        seg_pos_field (vb, CHAIN_ID, CHAIN_ID, 0, 0, field_start, field_len, 0, field_len + 1); // just a numeric delta    
    else
        seg_by_did_i (vb, field_start, field_len, CHAIN_ID, field_len + 1);

    SEG_EOL (CHAIN_EOL, true);

    // set containing 1 or more alignments of SIZE, DLUFT, PRIM, last alignment in set is only SIZE
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
        // last line doesn't have DLUFT and DPRIM - seg the "standard" separator and then delete it 
        // (this way we don't ruin the all_the_same of CHAIN_SEP)
        else {
            seg_by_did_i (vb, &last_gap_sep, 1, CHAIN_SEP, 0); 
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // prim_gap

            seg_by_did_i (vb, &last_gap_sep, 1, CHAIN_SEP, 0); 
            seg_by_did_i (vb, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // luft_gap
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
        .items               = { { .dict_id = { _CHAIN_SIZE } },
                                 { .dict_id = { _CHAIN_SEP }  }, 
                                 { .dict_id = { _CHAIN_GAPS } }, // prim_gap
                                 { .dict_id = { _CHAIN_SEP }  }, 
                                 { .dict_id = { _CHAIN_GAPS } }, // luft_gap
                                 { .dict_id = { _CHAIN_EOL }  } }
    };

    container_seg (vb, CTX(CHAIN_SET), (ContainerP)&alignment_set, 0, 0, 0);

    // Empty line after alignment set
    GET_LAST_ITEM (EmptyLine);
    ASSSEG0 (!field_len, field_start, "Expecting an empty line after alignment set");
    SEG_EOL (CHAIN_EOL, false); 

    return next_field;
}

//--------------
// PIZ functions
//--------------

static void chain_piz_initialize_copy_contigs (Buffer *contigs, Buffer *dict, Buffer *counts, DidIType name_did_i, bool get_LN_from_ref)
{
    // copy dictionary buffer
    buf_copy (evb, dict, &ZCTX(name_did_i)->dict, char, 0, 0, "chain_contig_dict");

    // copy counts buffer if requested
    if (counts)
        buf_copy (evb, counts, &ZCTX(name_did_i)->counts, int64_t, 0, 0, "chain_contig_counts");

    ARRAY (CtxWord, contig_words, ZCTX(name_did_i)->word_list);
    dict->param = contig_words_len; // dict param contains number of words

    // copy word_list into a contigs buffer - for convenience of using ctx_build_zf_ctx_from_contigs later
    buf_alloc (evb, contigs, 0, contig_words_len, RefContig, 0, "chain_contigs");
    contigs->len = contig_words_len;
    ARRAY (RefContig, ctg, *contigs);
    
    for (uint32_t i=0; i < contig_words_len; i++) {
        PosType LN=0;

        if (get_LN_from_ref) {
            const char *contig = ctx_get_snip_by_word_index (&z_file->contexts[name_did_i], i, 0, 0);
            LN = ref_contigs_get_contig_length (gref, i, contig, contig_words[i].snip_len, false);
        }
        
        ctg[i] = (RefContig) { .char_index = contig_words[i].char_index,
                               .snip_len   = contig_words[i].snip_len,
                               .max_pos    = LN };
    }
}

// called after reconstructing the txt header and before compute threads
bool chain_piz_initialize (void)
{
    mutex_initialize (chain_mutex);

    if (flag.reading_chain) {
        chain_piz_initialize_copy_contigs (&prim_contigs, &prim_contig_dict, &prim_contig_counts, CHAIN_NAMEPRIM, false);
        chain_piz_initialize_copy_contigs (&luft_contigs, &luft_contig_dict, NULL, CHAIN_NAMELUFT, true);
    }

    return true; // proceed with PIZ
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_BACKSPACE)
{
    ASSERT0 (vb->txt_data.len, "vb->txt_data.len is 0");
    vb->txt_data.len--;

    Context *gaps_ctx = CTX(CHAIN_GAPS); // note: we can't rely on the constants CHAIN_* in PIZ
    gaps_ctx->last_value.i = gaps_ctx->local.param = 0; // no prim_gap and luft_gap

    return false;
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_ENDLUFT)
{
    new_value->i = CTX(CHAIN_STARTLUFT)->last_value.i +
                   CTX(CHAIN_ENDPRIM)  ->last_value.i -
                   CTX(CHAIN_STARTPRIM)->last_value.i;
                   
    RECONSTRUCT_INT (new_value->i);
    return true; // new_value has been set
}

SPECIAL_RECONSTRUCTOR (chain_piz_special_SIZE)
{
    new_value->i = CTX(CHAIN_ENDPRIM)->last_value.i - CTX(CHAIN_STARTPRIM)->last_value.i;
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
        const char *luft_chrom = ctx_get_words_snip (ZCTX(CHAIN_NAMELUFT), aln[i].luft_chrom);
        const char *prim_chrom = ctx_get_words_snip (ZCTX(CHAIN_NAMEPRIM), aln[i].prim_chrom);

        iprintf ("Primary: %s %"PRId64"-%"PRId64" Luft: %s %"PRId64"-%"PRId64" Xstrand=%c\n",
                 prim_chrom, aln[i].prim_first_1pos, aln[i].prim_last_1pos, 
                 luft_chrom, aln[i].luft_first_1pos, aln[i].luft_last_1pos,
                 aln[i].is_xstrand ? 'X' : '-');
    }
}

// initialize alignment set from alignment set header
static inline void chain_piz_filter_init_alignment_set (VBlock *vb)
{
    reconstruct_from_ctx (vb, CHAIN_VERIFIED, 0, false); // sets last_int

    mutex_lock (chain_mutex);
    next_luft_0pos = vb->last_int(CHAIN_STARTLUFT); 
    next_prim_0pos = vb->last_int(CHAIN_STARTPRIM); 
    mutex_unlock (chain_mutex);
}

// store prim_gap value in local.param, as last_value.i will be overridden by luft_gap, as they share the same context
static inline void chain_piz_filter_save_prim_gap (VBlock *vb)
{
    CTX(CHAIN_GAPS)->local.param = vb->last_int(CHAIN_GAPS); 
}

// ingest one alignment
static inline void chain_piz_filter_ingest_alignmet (VBlock *vb)
{
    if (!vb->last_int (CHAIN_VERIFIED)) goto done; // a line Seg couldn't verify - don't ingest

    mutex_lock (chain_mutex);

    buf_alloc (NULL, &chain, 1, 0, ChainAlignment, 2, "chain");  // initial allocation is in chain_load

    int64_t size = vb->last_int(CHAIN_SIZE);

    ASSINP (*last_txt(vb, CHAIN_STRNDPRIM) == '+', "Chain file contains alignments with tStrand=\"%c\" (strand of Primary reference), this is not support by Genozip.",
            *last_txt(vb, CHAIN_STRNDLUFT));

    bool is_xstrand = (*last_txt(vb, CHAIN_STRNDLUFT) == '-');
    PosType last_prim_0pos = next_prim_0pos + size - 1; // last POS of the dst alignment, 0-based
    PosType last_luft_0pos = next_luft_0pos + size - 1; // last POS of the dst alignment, 0-based
    PosType size_dst = vb->last_int (CHAIN_SIZELUFT);  // contig LN

    // note on negative strand: the source region (the source always has a positive strand), is aligned to a region on the destination
    // contig, but the positions given are counting starting from end of the contig (the last base of the contig is 0), and the sequence
    // is complemented. For our chain alignment representation, we convert the coordinates to positive-strand terms, and set XSTRAND=X
    // to record that this is a reverse strand

    NEXTENT (ChainAlignment, chain) = (ChainAlignment){ 
        .prim_chrom      = vb->last_index(CHAIN_NAMEPRIM),
        .prim_first_1pos = 1 + next_prim_0pos,  // +1 bc our alignments are 1-based vs the chain file that is 0-based
        .prim_last_1pos  = 1 + last_prim_0pos,

        // use alt chroms as they appear in the reference file, rather than their name in the chain file
        .luft_chrom      = ref_alt_get_final_index (vb->last_index(CHAIN_NAMELUFT)), 
        .luft_first_1pos = is_xstrand ? (size_dst - last_luft_0pos) : 1 + next_luft_0pos, // +1 to convert to 1-base. note that "size_dst" has a built-in +1.
        .luft_last_1pos  = is_xstrand ? (size_dst - next_luft_0pos) : 1 + last_luft_0pos,

        .is_xstrand     = is_xstrand
    };

    Context *gaps_ctx   = CTX(CHAIN_GAPS);
    next_prim_0pos += size + gaps_ctx->local.param;  // prim_gap
    next_luft_0pos += size + gaps_ctx->last_value.i; // luft_gap

    mutex_unlock (chain_mutex);

done:
    vb->drop_curr_line = "chain";
}

// verify that adding up all alignments and gaps, results in the end position specified in the header
static inline void chain_piz_filter_verify_alignment_set (VBlock *vb)
{
    mutex_lock (chain_mutex);

    ASSINP (next_prim_0pos == vb->last_int(CHAIN_ENDPRIM),
            "Bad data in chain file %s: Expecting alignments to add up to ENDPRIM=%"PRId64", but they add up to %"PRId64":\n%*s",
            z_name, vb->last_int(CHAIN_ENDPRIM), next_prim_0pos, 
            (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    ASSINP (next_luft_0pos == vb->last_int(CHAIN_ENDLUFT),
            "Bad data in chain file %s: Expecting alignments to add up to ENDLUFT=%"PRId64", but they add up to %"PRId64":\n%*s",
            z_name, vb->last_int(CHAIN_ENDLUFT), next_luft_0pos, 
            (int)(vb->txt_data.len - vb->line_start), ENT (char, vb->txt_data, vb->line_start));

    mutex_unlock (chain_mutex);
}

// set contig lengths according to tSize and qSize
static inline void chain_piz_filter_add_contig_length (VBlock *vb)
{
    mutex_lock (chain_mutex);

    RefContig *prim_rc = ENT (RefContig, prim_contigs, vb->last_index(CHAIN_NAMEPRIM));
    RefContig *luft_rc = ENT (RefContig, luft_contigs, vb->last_index(CHAIN_NAMELUFT));

    prim_rc->max_pos = vb->last_int(CHAIN_SIZEPRIM);
    luft_rc->max_pos = vb->last_int(CHAIN_SIZELUFT);

    mutex_unlock (chain_mutex);
}

static inline void chain_piz_filter_add_qName_chr (VBlock *vb)
{
    char c_3 = *(LASTENT (char, vb->txt_data)-3);
    char c_2 = *(LASTENT (char, vb->txt_data)-2);
    char c_1 = *(LASTENT (char, vb->txt_data)-1);
    
    // case: last reconstructed text look like " 2 " - single character chrom
    if (c_2 == ' ') { 
        vb->txt_data.len -= 2;
        RECONSTRUCT ("chr", 3);
        RECONSTRUCT1 (c_1);
        RECONSTRUCT1 (' ');
    }

    // case: last reconstructed text look like " 22 " - double digit chrom
    else if (c_3 == ' ' && IS_DIGIT (c_2) && IS_DIGIT (c_1)) { 
        vb->txt_data.len -= 3;
        RECONSTRUCT ("chr", 3);
        RECONSTRUCT1 (c_2);
        RECONSTRUCT1 (c_1);
        RECONSTRUCT1 (' ');
    }

    // case: last reconstructed text look like " MT " 
    else if (c_3 == ' ' && c_2 == 'M' && c_1 == 'T') { 
        vb->txt_data.len -= 3;
        RECONSTRUCT ("chrM ", 5);
    }
}

CONTAINER_FILTER_FUNC (chain_piz_filter)
{
    // ingest alignments (and report Chain file data issues) only when consuming with --chain, not when merely pizzing
    if (flag.reading_chain) { 

        // before alignment-set first EOL and before alignments - initialize next_luft_0pos and next_prim_0pos
        if (dict_id.num == _CHAIN_TOPLEVEL && item == 13) 
            chain_piz_filter_init_alignment_set (vb);

        // save prim_gap before reconstructing luft_gap (prim_gap was just processed)
        else if (dict_id.num == _CHAIN_SET && item == 4) 
            chain_piz_filter_save_prim_gap (vb);

        // before EOF of each alignment, ingest alignment (only if it passed verification during Seg)
        else if (dict_id.num == _CHAIN_SET && item == 5) 
            chain_piz_filter_ingest_alignmet (vb);

        // before alignment-set second EOL and after alignments - verify that numbers add up, and also set contigs
        else if (dict_id.num == _CHAIN_TOPLEVEL && item == 15 && vb->last_int (CHAIN_VERIFIED)) {
            chain_piz_filter_verify_alignment_set (vb);
            chain_piz_filter_add_contig_length (vb);
        }
    }

    // genocat of a chain file (not --chain)
    else {
        // rewrite NAMELUFT 1,...,Y -> chr1,...,chrY ; MT -> chrM
        if (flag.with_chr && dict_id.num == _CHAIN_TOPLEVEL && item == 8) 
            chain_piz_filter_add_qName_chr (vb);
    }
    
    return true;    
}

// sort the alignmants by (luft_contig_index, qstart)
static int chain_sorter (const void *a_, const void *b_)
{   
    ChainAlignment *a = (ChainAlignment *)a_, *b = (ChainAlignment *)b_;

    if (a->prim_chrom != b->prim_chrom) 
        return a->prim_chrom - b->prim_chrom;
    else 
        return a->prim_first_1pos - b->prim_first_1pos;
}

// main thread: ZIP with --chain, reading txt header: returns null-terminated string of contig, or NULL if contig_i is out of range
const char *chain_get_luft_contig (uint32_t contig_i, PosType *length)
{
    if (contig_i >= luft_contigs.len) return NULL;

    const RefContig *ctg = ENT (RefContig, luft_contigs, contig_i);

    *length = ctg->max_pos;
    return ENT (char, luft_contig_dict, ctg->char_index);
}

uint64_t chain_get_num_prim_contigs (void)
{
    return prim_contigs.len;
}

// returns the contig index in prim_contigs, word WORD_INDEX_NONE if there is none
WordIndex chain_get_prim_contig_index_by_name (const char *contig, unsigned contig_len, bool recursive)
{
    WordIndex ctg_i;
    ARRAY (RefContig, prim_ctgs, prim_contigs);
    
    for (ctg_i=0; ctg_i < prim_ctgs_len; ctg_i++) {
        const char *prim_contig = ENT (char, prim_contig_dict, prim_ctgs[ctg_i].char_index);

        if (prim_ctgs[ctg_i].snip_len == contig_len && 
            !memcmp (contig, prim_contig, contig_len)) break; // found
    }

    // case: not found (this can happen if VCF file doesn't have contigs in its header), check for alternative names.
    // similar logic to ref_alt_chroms_get_alt_index
    if (ctg_i == prim_ctgs_len && !recursive) {

        // 22 -> chr22 (1->22, X, Y, M, MT chromosomes)
        if ((contig_len == 1 && (IS_DIGIT (contig[0]) || contig[0]=='X' || contig[0]=='Y')) ||
            (contig_len == 2 && ((IS_DIGIT (contig[0]) && IS_DIGIT (contig[1]))))) {

            char chr_chrom[5] = "chr";
            chr_chrom[3] = contig[0];
            chr_chrom[4] = (contig_len == 2 ? contig[1] : 0);

            ctg_i = chain_get_prim_contig_index_by_name (chr_chrom, contig_len+3, true); 
        }

        // M, MT, chrM, chrMT 
        else if ((contig_len==4 && contig[0]=='c' && contig[1]=='h' && contig[2]=='r' && contig[3]=='M') || 
                 (contig_len==5 && contig[0]=='c' && contig[1]=='h' && contig[2]=='r' && contig[3]=='M' && contig[4]=='T') || 
                 (contig_len==1 && contig[0]=='M') ||
                 (contig_len==2 && contig[0]=='M' && contig[1]=='T')) {

                                          ctg_i = chain_get_prim_contig_index_by_name ("chrMT", 5, true); 
            if (ctg_i == WORD_INDEX_NONE) ctg_i = chain_get_prim_contig_index_by_name ("chrM",  4, true); 
            if (ctg_i == WORD_INDEX_NONE) ctg_i = chain_get_prim_contig_index_by_name ("MT",    2, true); 
            if (ctg_i == WORD_INDEX_NONE) ctg_i = chain_get_prim_contig_index_by_name ("M",     1, true); 
        }
        
        // Chr? or Chr?? -> ? or ??
        else if ((contig_len == 4 || contig_len == 5) && !memcmp (contig, "chr", 3))
            ctg_i = chain_get_prim_contig_index_by_name (&contig[3], contig_len-3, true); 
    }

    // note: PRIMNAME contains all contigs, copied from the primary reference. We check counts, to verify that the contig is actually used in the chain file.
    if (ctg_i != WORD_INDEX_NONE && ctg_i < prim_ctgs_len && *ENT (int64_t, prim_contig_counts, ctg_i) > 0)
        return ctg_i;
    else
        return WORD_INDEX_NONE;
}

// called by vcf_lo_zip_initialize when zipping with --chain to copy luft_contigs to z_file->ochrom and prim_contigs to z_file->CHROM
void chain_copy_contigs_to_z_file (DidIType luft_contig_did_i)
{
    ctx_build_zf_ctx_from_contigs (luft_contig_did_i, &luft_contigs, &luft_contig_dict);
    // note: primary contigs are copied in txtheader_zip_prepopulate_contig_ctxs
}

static void chain_contigs_show (const Buffer *contigs, const Buffer *dict, const char *type, const char *ref_filename)
{
    ARRAY (const RefContig, cn, *contigs);

    iprintf ("\n%s chain file contigs that also exist in the reference file:\n", type);
    for (uint32_t i=0; i < cn_len; i++) {

        const char *chrom_name = ENT (const char, *dict, cn[i].char_index);

        if (cn[i].snip_len && cn[i].max_pos > 0)
            iprintf ("%s %s length=%"PRId64"\n", type, chrom_name, cn[i].max_pos);
    }
}

// zip: load chain file as a result of genozip --chain
void chain_load (void)
{
    ASSERTNOTNULL (flag.reading_chain);
    SAVE_FLAGS_AUX;

    buf_alloc (evb, &chain, 0, 1000, ChainAlignment, 1, "chain"); // must be allocated by main thread
    
    z_file = file_open (flag.reading_chain, READ, Z_FILE, DT_CHAIN);    
    SectionHeaderGenozipHeader header;
    zfile_read_genozip_header (&header);
    
    ASSINP (z_file->data_type == DT_CHAIN, "expected %s to be a genozip'ed chain file, but its a %s file. Tip: compress the chain with \"genozip --input chain\"", 
            z_name, dt_name (z_file->data_type));

    ASSINP (ref_get_filename (gref) && ref_get_filename (prim_ref), 
            "%s is an invalid chain file: it was not compressed with source and destination references using genozip --reference", z_name);

    z_file->basename = file_basename (flag.reading_chain, false, "(chain-file)", NULL, 0);

    TEMP_VALUE (command, PIZ);

    if (flag.reference) flag.reference = REF_LIFTOVER;

    flag.no_writer = true;
    flag.genocat_no_reconstruct = false;
    flag.genocat_global_area_only = false;
    
    Dispatcher dispachter = piz_z_file_initialize (false);

    // load both references, now that it is set (either explicitly from the command line, or implicitly from the chain GENOZIP_HEADER)
    SAVE_VALUE (z_file); // actually, read the references first
    ref_load_external_reference (gref, NULL);
    ref_load_external_reference (prim_ref, NULL);
    RESTORE_VALUE (z_file);

    // test for matching MD5 between external references and reference in the chain file header (doing it here, because reference is read after chain file, if its explicitly specified)
    digest_verify_ref_is_equal (gref, header.ref_filename, header.ref_file_md5);
    digest_verify_ref_is_equal (prim_ref, header.dt_specific.chain.prim_filename, header.dt_specific.chain.prim_file_md5);

    flag.quiet = true; // don't show progress indicator for the chain file - it is very fast 
    flag.maybe_vb_modified_by_reconstructor = true; // we drop all the lines

    piz_one_txt_file (dispachter, false);

    // --show-chain-contigs
    if (flag.show_chain_contigs) {
        chain_contigs_show (&prim_contigs, &prim_contig_dict, "PRIMARY", ref_get_filename (prim_ref));
        chain_contigs_show (&luft_contigs, &luft_contig_dict, "LUFT",    ref_get_filename (gref));
        if (exe_type == EXE_GENOCAT) exit_ok();  // in genocat this, not the data
    }

    // sort the alignmants by (luft_contig_index, qstart)
    qsort (chain.data, chain.len, sizeof (ChainAlignment), chain_sorter);
        
    if (flag.show_chain) {
        chain_display_alignments();
        exit_ok();
    }

    chain_filename = file_make_unix_filename (z_name); // full-path unix-style filename, allocates memory

    file_close (&z_file, false, false);
    file_close (&txt_file, false, false); // close the txt_file object we created (even though we didn't open the physical file). it was created in file_open called from txtheader_piz_read_and_reconstruct.
    
    // recover globals
    RESTORE_VALUE (command);
    RESTORE_FLAGS;

    flag.reference     = REF_LIFTOVER; // will be needed by the files we are about to compress with --chain
    flag.reading_chain = NULL;
}

// reset to factory defaults
void chain_destroy (void)
{
    buf_destroy (&chain);
    buf_destroy (&prim_contig_dict);
    buf_destroy (&prim_contigs);
    buf_destroy (&luft_contig_dict);
    buf_destroy (&luft_contigs);
    mutex_destroy (chain_mutex);
    next_luft_0pos = next_prim_0pos=0;
    chain_filename = NULL; 
}

// binary search for the first alignment of a src contig
static ChainAlignment *chain_get_first_aln (WordIndex prim_contig_index, int32_t start, int32_t end) 
{
    if (end < start) return NULL; // prim_contig doesn't exist in chain file

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = ENT (ChainAlignment, chain, mid);

    // case: success - same prim_contig, and previous aln has a different prim_contig 
    if (aln->prim_chrom == prim_contig_index && (!mid || (aln-1)->prim_chrom != prim_contig_index))
        return aln;

    // case: prim_contig less than aln's, OR is the same, but aln isn't the first with this prim_contig - search lower half
    if (prim_contig_index <= aln->prim_chrom)
        return chain_get_first_aln (prim_contig_index, start, mid-1);
    
    // case: prim_contig is more than aln's - search upper half
    else
        return chain_get_first_aln (prim_contig_index, mid+1, end);
}

// append luft_contigs buffer with all dst chroms indicies for which there exists alignment with prim_chrom
void chain_append_all_luft_contig_index (const char *prim_contig_name, unsigned prim_contig_name_len, Buffer *luft_contigs)
{
    WordIndex prim_contig_index = chain_get_prim_contig_index_by_name (prim_contig_name, prim_contig_name_len, false);
    WordIndex prev_dst = WORD_INDEX_NONE;

    if (prim_contig_index == WORD_INDEX_NONE) return; // this contig is not in the chain file

    for (ChainAlignment *aln = chain_get_first_aln (prim_contig_index, 0, chain.len-1);
         aln && aln < LASTENT (ChainAlignment, chain) && aln->prim_chrom == prim_contig_index;
         aln++)
         
         if (aln->luft_chrom != prev_dst) { // note: this prevents consecutive duplicates, but not non-consecutive duplicates
            buf_alloc (NULL, luft_contigs, 100, 1, WordIndex, 2, NULL);
            NEXTENT (WordIndex, *luft_contigs) = prev_dst = aln->luft_chrom;
         }
}

// get luft_contig, prim_pos from prim_contig, luft_pos (binary search on chain)
// returns true if successful
static bool chain_get_liftover_coords_do (WordIndex prim_contig_index, PosType prim_1pos, 
                                          int32_t start, int32_t end,
                                          WordIndex *luft_contig_index, PosType *luft_1pos, bool *is_xstrand, uint32_t *aln_i) // out
{
    // case: no mapping
    if (end < start) {
        if (luft_contig_index) *luft_contig_index = WORD_INDEX_NONE;
        if (luft_1pos) *luft_1pos = 0;
        if (is_xstrand) *is_xstrand = false;
        return false;
    }

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = ENT (ChainAlignment, chain, mid);

    // case: success
    if (aln->prim_chrom == prim_contig_index && 
        aln->prim_first_1pos <= prim_1pos && 
        aln->prim_last_1pos >= prim_1pos) {
            if (luft_contig_index) *luft_contig_index = aln->luft_chrom;
            
            PosType pos_offset = prim_1pos - aln->prim_first_1pos; // offset of POS from the beginning of the src (src is always positive strand)
            if (luft_1pos) *luft_1pos = aln->is_xstrand ? (aln->luft_last_1pos  - pos_offset)  // dst is reverse strand - offset is from the end of the alignment, going back
                                                      : (aln->luft_first_1pos + pos_offset); // dst is positive strand - offset is from the start of the alignment
            
            if (is_xstrand) *is_xstrand = aln->is_xstrand;
            if (aln_i) *aln_i = mid;
            return true;
        }

    // case: prim_contig, luft_pos is less than aln - search lower half
    if (prim_contig_index < aln->prim_chrom ||
        (prim_contig_index == aln->prim_chrom && prim_1pos < aln->prim_first_1pos))
        return chain_get_liftover_coords_do (prim_contig_index, prim_1pos, start, mid-1, luft_contig_index, luft_1pos, is_xstrand, aln_i);
    
    // case: prim_contig,luft_pos is more than aln - search upper half
    else
        return chain_get_liftover_coords_do (prim_contig_index, prim_1pos, mid+1, end, luft_contig_index, luft_1pos, is_xstrand, aln_i);
}

bool chain_get_liftover_coords (WordIndex prim_contig_index, PosType prim_1pos, 
                                WordIndex *luft_contig_index, PosType *luft_1pos, bool *is_xstrand, uint32_t *aln_i) // out
{
    return chain_get_liftover_coords_do (prim_contig_index, prim_1pos, 0, chain.len-1, luft_contig_index, luft_1pos, is_xstrand, aln_i);
}

PosType chain_get_aln_gap_after (uint32_t aln_i)
{
    const ChainAlignment *aln = ENT (const ChainAlignment, chain, aln_i);

    // case: this is the last alignment of this contig
    if (aln_i == chain.len - 1 || aln->prim_chrom != (aln+1)->prim_chrom) { 
        const RefContig *prim_rc = ENT (const RefContig, prim_contigs, aln->prim_chrom);
        return prim_rc->max_pos - aln->prim_last_1pos;
    }
    else
        return (aln+1)->prim_first_1pos - aln->prim_last_1pos - 1;
}

PosType chain_get_aln_last_pos (uint32_t aln_i)
{
    return ENT (const ChainAlignment, chain, aln_i)->prim_last_1pos;
}
