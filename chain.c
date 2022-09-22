// ------------------------------------------------------------------
//   chain.c
//   Copyright (C) 2021-2022 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.
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
#include "segconf.h"
#include "version.h"
#include "website.h"
#include "random_access.h"
#include "contigs.h"
#include "chrom.h"

typedef struct {
    WordIndex prim_chrom; // index in CHAIN_NAMEPRIM
    PosType prim_first_1pos, prim_last_1pos;

    WordIndex luft_chrom; // index in CHAIN_NAMELUFT - the nodes+dicts copied to the genozipped file (eg as VCF_oCHROM)
    PosType luft_first_1pos, luft_last_1pos;

    bool is_xstrand;

    // used by chain_display_alignments
    int32_t aln_i;
    int32_t overlap_aln_i; // the luft range overlaps with with luft range
} ChainAlignment;

static Buffer chain = EMPTY_BUFFER;   // immutable after loaded
static Mutex chain_mutex = {};        // protect chain while PIZ: loading, SEG: accessing seg_alns_so_far

static Buffer seg_alns_so_far = EMPTY_BUFFER; // Seg: alignments encountered so far - used for de-dupping in --match

// 1) When pizzing a chain as a result of genozip --chain: chain_piz_initialize: prim_contig and luft_contig are generated 
//    from z_file->contexts of a chain file 
// 2) Contig lengths are added when pizzing each alignment set, in chain_piz_filter_add_contig_length
// 3) At the beginning of ZIP of a file with --chain: luft_contig copied to z_file->contexts of the file being compressed
// 4) vcf_lo_get_liftover_coords prim_contig is consulted to get map the index of the primary chrom from the node_index
//    of the file being compressed to the word_index of prim_contig in the chain data
static ContigPkg prim_ctgs = { .name = "chain_prim_ctgs" }, luft_ctgs = { .name = "chain_luft_ctgs" };

char *chain_filename = NULL; // global - chain filename

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockCHAIN {
    VBLOCK_COMMON_FIELDS
    PosType next_luft_0pos, next_prim_0pos; // used for loading chain
} VBlockCHAIN;

unsigned chain_vb_size (DataType dt) { return sizeof (VBlockCHAIN); }

void chain_vb_release_vb (VBlockCHAIN *vb)
{
    vb->next_luft_0pos = vb->next_prim_0pos = 0;
}

static void chain_display_alignments (void) 
{
    iprint0 ("##fileformat=GENOZIP-CHAIN\n");
    iprintf ("##primary_reference=%s\n", ref_get_filename (prim_ref));
    iprintf ("##luft_reference=%s\n", ref_get_filename (gref));
    iprint0 ("##documentation=" WEBSITE_CHAIN "\n");
    iprint0 ("#ALN_I\tPRIM_CONTIG\tPRIM_START\tPRIM_END\tLUFT_CONTIG\tLUFT_START\tLUFT_ENDS\tXSTRAND\tALN_OVERLAP\n");

    for_buf (ChainAlignment, aln, chain) {
        if (aln->aln_i == -1) continue; // duplicate

        rom luft_chrom = ctx_get_words_snip (ZCTX(CHAIN_NAMELUFT), aln->luft_chrom);
        rom prim_chrom = ctx_get_words_snip (ZCTX(CHAIN_NAMEPRIM), aln->prim_chrom);

        iprintf ("%u\t%s\t%"PRId64"\t%"PRId64"\t%s\t%"PRId64"\t%"PRId64"\t%c",
                 aln->aln_i, prim_chrom, aln->prim_first_1pos, aln->prim_last_1pos, 
                 luft_chrom, aln->luft_first_1pos, aln->luft_last_1pos,
                 aln->is_xstrand ? 'X' : '-');

        if (aln->overlap_aln_i)
            iprintf ("\t%u\n", aln->overlap_aln_i);
        else
            iprint0 ("\n");
    }
}

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

// main thread: writing data-type specific fields to genozip header
void chain_zip_genozip_header (SectionHeaderGenozipHeader *header)
{
    if (IS_REF_MAKE_CHAIN) {
        strncpy (header->chain.prim_filename, ref_get_filename (prim_ref), REF_FILENAME_LEN-1);
        header->chain.prim_file_md5 = ref_get_file_md5 (prim_ref);
    }
}

void chain_zip_initialize (void)
{
    ASSINP0 (IS_REF_EXTERNAL || flag.force, "Please specify the destination reference file with --reference or use --force to override this");

    if (IS_REF_EXTERNAL && z_file->num_txts_so_far == 1) {
        flag.reference = REF_MAKE_CHAIN; // since we're not attempting to use the data from the reference to compress the chain file itself
       
        // initialize NAMEPRIM and NAMELUFT from references, so chain and ref contigs have identical indices
        ctx_populate_zf_ctx_from_contigs (prim_ref, CHAIN_NAMEPRIM, ref_get_ctgs (prim_ref)); 
        ctx_populate_zf_ctx_from_contigs (gref,     CHAIN_NAMELUFT, ref_get_ctgs (gref));
    }
    else
        segconf.chain_mismatches_ref = true; // chain file can be compressed, but not used with --chain

    if (flag.match_chrom_to_reference) {
        mutex_initialize (chain_mutex);
        seg_alns_so_far.len = 0;     // free allocation by previous files
        buf_alloc (evb, &seg_alns_so_far, 0, 100, ChainAlignment, 1, "seg_alns_so_far"); // first allocation must be in main thread, so the buffer can be added to evb buf_list
    }
}

// Tests lengths of prim/luft contigs of first line in chain file vs the two references. Error if they don't match, but swapping them matches.
static void chain_verify_references_not_reversed (VBlockP vb)
{
    rom newline = memchr (vb->txt_data.data, '\n', vb->txt_data.len);
    if (!newline) return;

    str_split (vb->txt_data.data, newline - vb->txt_data.data, 13, ' ', item, true);
    if (!n_items) return;

    PosType prim_len_ref, prim_len_chain, luft_len_ref, luft_len_chain;
    if (!str_get_int (STRi(item, 3), &prim_len_chain)) return;
    if (!str_get_int (STRi(item, 8), &luft_len_chain)) return;

    prim_len_ref = ref_contigs_get_contig_length (prim_ref, WORD_INDEX_NONE, STRi(item,2), false);
    luft_len_ref = ref_contigs_get_contig_length (gref,     WORD_INDEX_NONE, STRi(item,7), false);
    
    // if lengths don't match (either because they differ, or because len_ref=-1 bc contig not found in ref file), check in reverse
    if (prim_len_ref != prim_len_chain || luft_len_ref != luft_len_chain) {
        prim_len_ref = ref_contigs_get_contig_length (gref,     WORD_INDEX_NONE, STRi(item,2), false);
        luft_len_ref = ref_contigs_get_contig_length (prim_ref, WORD_INDEX_NONE, STRi(item,7), false);
        
        if (prim_len_ref == prim_len_chain && luft_len_ref == luft_len_chain) 
            ABORTINP ("Reference files %s and %s are given in the command line in reverse order. Expecting Primary reference first and Luft reference second.\n", ref_get_filename (prim_ref), ref_get_filename (gref)); // not WARN bc segconf_calculate set flag.quiet=true
    }
}

void chain_seg_initialize (VBlockP vb)
{
    // check if the references are reversed
    if (segconf.running && flag.reference)
        chain_verify_references_not_reversed (vb);

    ctx_set_no_stons (vb, CHAIN_NAMELUFT, CHAIN_NAMEPRIM, CHAIN_TOPLEVEL, 
                      CHAIN_CHAIN, CHAIN_SEP, CHAIN_VERIFIED, CHAIN_STRNDPRIM,  // keep in b250 so it can be eliminated as all_the_same
                      CHAIN_STARTPRIM, CHAIN_ENDPRIM, CHAIN_STARTLUFT, CHAIN_ENDLUFT, CHAIN_ID, // required by seg_pos_field
                      DID_EOL);
                      
    ctx_set_store (VB, STORE_INT, CHAIN_STARTPRIM, CHAIN_ENDPRIM, CHAIN_SIZEPRIM,
                   CHAIN_STARTLUFT, CHAIN_ENDLUFT, CHAIN_SIZELUFT,
                   CHAIN_SIZE, CHAIN_VERIFIED, CHAIN_GAPS, CHAIN_ID, DID_EOL);

    ctx_set_store (VB, STORE_INDEX, CHAIN_NAMEPRIM, CHAIN_NAMELUFT, DID_EOL);

    CTX(CHAIN_NAMEPRIM)-> no_vb1_sort =
    CTX(CHAIN_NAMELUFT)-> no_vb1_sort = true;

    CTX(CHAIN_NAMEPRIM)-> counts_section = 
    CTX(CHAIN_NAMELUFT)-> counts_section = true;
}

void chain_seg_finalize (VBlockP vb)
{
    // top level snip - IMPORTNAT - chain_piz_filter needs to changed and support backward compatability
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 16,
        .items        = { { .dict_id = { _CHAIN_CHAIN },     .separator = {' '} },  // 0
                          { .dict_id = { _CHAIN_SCORE },     .separator = {' '} },  // 1
                          { .dict_id = { _CHAIN_NAMEPRIM },  .separator = {' '} },  // 2
                          { .dict_id = { _CHAIN_SIZEPRIM },  .separator = {' '} },  // 3
                          { .dict_id = { _CHAIN_STRNDPRIM }, .separator = {' '} },  // 4
                          { .dict_id = { _CHAIN_STARTPRIM }, .separator = {' '} },  // 5
                          { .dict_id = { _CHAIN_ENDPRIM },   .separator = {' '} },  // 6
                          { .dict_id = { _CHAIN_NAMELUFT },  .separator = {' '} },  // 7
                          { .dict_id = { _CHAIN_SIZELUFT },  .separator = {' '} },  // 8
                          { .dict_id = { _CHAIN_STRNDLUFT }, .separator = {' '} },  // 9
                          { .dict_id = { _CHAIN_STARTLUFT }, .separator = {' '} },  // 10
                          { .dict_id = { _CHAIN_ENDLUFT },   .separator = {' '} },  // 11
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

static void chain_seg_size_field (VBlockCHAIN *vb, rom field_start, int field_len)
{
    int64_t size;
    ASSSEG (str_get_int (field_start, field_len, &size), field_start, "Expecting size to be an integer, but found %*s", field_len, field_start);

    Context *ctx = CTX(CHAIN_SIZE);

    if (vb->last_int(CHAIN_ENDPRIM) - vb->last_int(CHAIN_STARTPRIM) == size)  // happens if the alignment set has only one alignment
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_SIZE }), 2, ctx, field_len + 1);

    else
        seg_integer_or_not (VB, ctx, field_start, field_len, field_len + 1);
}

static void chain_seg_luft_end_field (VBlockCHAIN *vb, rom field_start, int field_len)
{
    Context *ctx = CTX(CHAIN_ENDLUFT);

    ASSSEG (str_get_int (field_start, field_len, &ctx->last_value.i), field_start, "Expecting qEnd to be an integer, but found %*s", field_len, field_start);

    if (vb->last_int(CHAIN_ENDPRIM) - vb->last_int(CHAIN_STARTPRIM) ==
        ctx->last_value.i           - vb->last_int(CHAIN_STARTLUFT)) {
        seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_ENDLUFT }), 2, ctx, field_len + 1);        
    }

    else
        seg_pos_field (VB, CHAIN_ENDLUFT, CHAIN_STARTLUFT, 0, 0, field_start, field_len, 0, field_len + 1);
}

// returns false if length is inconsistent with reference, and true if length is the same or contig doesn't exist in reference
// returns *mismatch=true if contig doesn't exist with the exact name in the reference
static bool chain_seg_verify_contig (VBlockCHAIN *vb, Reference ref, WordIndex name_index, bool *once, STRp(name), STRp(last_name), PosType size_according_to_chain,
                                     bool *mismatch) // out
{
    WordIndex ref_index = name_index;
    if (flag.reference && ref_index >= ref_num_contigs(ref)) { 
        *mismatch = true;

        ref_index = contigs_get_matching (ref_get_ctgs(ref), STRa(name), size_according_to_chain, true, NULL);
        
        // case: this contig is not in the reference, even with an alternate name. 
        // note: if its in the reference with an alternate name, a warning was already displayed in chrom_seg_ex
        if (ref_index == WORD_INDEX_NONE && !flag.match_chrom_to_reference && !__atomic_test_and_set (once, __ATOMIC_RELAXED)) 
            iprintf ("\n\nWarning: Some contigs, for example \"%.*s\", appear the chain file %s, but don't appear in the %s reference %s. Use --match-chrom-to-reference to exclude these alignments from %s.\n", 
                        name_len, name, txt_name, ref==gref ? "LUFT" : "PRIMARY", ref_get_filename (ref), z_name);
    }

    // case: contig is in the reference - verify its size
    else if (name_index < ref_num_contigs (ref)) {
        
        PosType size_according_to_ref = ref_contigs_get_contig_length (ref, ref_index, 0, 0, true);

        if (size_according_to_chain != size_according_to_ref) {
            // output error in the first alignment of this contig
            ASSERTW (str_issame (name, last_name),
                    "Size of \"%.*s\" in chain file is %"PRId64", but in %s \"%s\" is %"PRId64". Excluding it.",
                    name_len, name, size_according_to_chain, ref_get_filename (ref), ref_contigs_get_name (ref, ref_index, NULL), size_according_to_ref);
            return false;
        }
    }

    return true; // verified
}

static WordIndex chain_seg_name_and_size (VBlockCHAIN *vb, Reference ref, bool is_luft, Did did_i_name, STRp(name), Did did_i_size, STRp(size), 
                                          bool *mismatch, bool *verified) // out
{
    PosType size_value=0;
    ASSSEG (str_get_int_range64 (STRa(size), 1, 1000000000000000, &size_value),
            size, "Invalid size field of \"%.*s\": %.*s", STRf(name), STRf(size));;

    seg_by_did (VB, STRa(size), did_i_size, size_len+1);
    
    bool is_alt;
    WordIndex node_index = chrom_seg_ex (VB, did_i_name, STRa (name), size_value, &is_alt, name_len+1, true, NULL);

    static bool once[2] = {};
    *verified &= is_alt || chain_seg_verify_contig (vb, ref, node_index, &once[is_luft], STRa (name), 
                                                    last_txt(VB, did_i_name), vb->last_txt_len(did_i_name), size_value, mismatch);

    seg_set_last_txt (VB, CTX(did_i_name), STRa(name)); // note: in some cases, this is already done by chrom_seg_ex. no harm.

    return node_index;
}

static bool chain_seg_is_duplicate (ChainAlignment aln)
{
    mutex_lock (chain_mutex);

    // check if this alignment was already encountered
    for (uint64_t i=0; i < seg_alns_so_far.len; i++) {
        ChainAlignment *old_aln =  B(ChainAlignment, seg_alns_so_far, i);
        if (!memcmp (old_aln, &aln, sizeof (ChainAlignment))) {
            mutex_unlock (chain_mutex);
            return true;  // a duplicate
        }
    }
    
    // first encounter with this alignment - add it to the buffer
    buf_alloc (NULL, &seg_alns_so_far, 1, 0, ChainAlignment, CTX_GROWTH, NULL);
    BNXT (ChainAlignment, seg_alns_so_far) = aln;

    mutex_unlock (chain_mutex);
    return false; // not duplicate
}

bool chain_zip_dts_flag (void)
{
    return segconf.chain_mismatches_ref;
}

rom chain_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockCHAIN *vb = (VBlockCHAIN *)vb_;
    rom next_field=field_start_line, field_start;
    unsigned field_len=0;
    char separator;
    bool mismatch=false;

    seg_create_rollback_point (VB, NULL, NUM_CHAIN_FIELDS, CHAIN_NAMELUFT, CHAIN_STRNDLUFT, CHAIN_STARTLUFT, CHAIN_ENDLUFT, CHAIN_SIZELUFT, 
                               CHAIN_NAMEPRIM, CHAIN_STRNDPRIM, CHAIN_STARTPRIM, CHAIN_ENDPRIM, CHAIN_SIZEPRIM, CHAIN_CHAIN, CHAIN_SCORE, 
                               CHAIN_ID, CHAIN_VERIFIED, CHAIN_SET, CHAIN_SIZE, CHAIN_GAPS, CHAIN_EOL, CHAIN_TOPLEVEL, CHAIN_SEP, CHAIN_DEBUG_LINES);
    typeof(vb->recon_size) save_recon_size = vb->recon_size;

    int32_t len = BAFTtxt - field_start_line;

    SEG_NEXT_ITEM_SP (CHAIN_CHAIN);
    SEG_NEXT_ITEM_SP (CHAIN_SCORE);
    
    GET_NEXT_ITEM_SP (CHAIN_NAMEPRIM);
    GET_NEXT_ITEM_SP (CHAIN_SIZEPRIM);
    bool verified = true; // start optimistically
    WordIndex prim_chrom = chain_seg_name_and_size (vb, prim_ref, false, CHAIN_NAMEPRIM, STRd(CHAIN_NAMEPRIM), CHAIN_SIZEPRIM, STRd(CHAIN_SIZEPRIM),
                                                    &mismatch, &verified);

    SEG_NEXT_ITEM_SP (CHAIN_STRNDPRIM);
    
    #define ASSERT_VALID_STRAND(f) ASSSEG (f##_len==1 && (*f##_str=='+' || *f##_str=='-'), f##_str, "Invalid strand character \"%.*s\", expecting '+' or '-'", f##_len, f##_str);
    ASSERT_VALID_STRAND (CHAIN_STRNDPRIM);

    GET_NEXT_ITEM_SP (CHAIN_STARTPRIM);
    seg_pos_field (vb_, CHAIN_STARTPRIM, CHAIN_ENDPRIM, 0, 0, STRd(CHAIN_STARTPRIM), 0, field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_ENDPRIM);
    seg_pos_field (vb_, CHAIN_ENDPRIM, CHAIN_STARTPRIM, 0, 0, STRd(CHAIN_ENDPRIM), 0,field_len + 1);

    random_access_update_first_last_pos (vb_, 0, prim_chrom, STRd(CHAIN_STARTPRIM), STRd (CHAIN_ENDPRIM));

    GET_NEXT_ITEM_SP (CHAIN_NAMELUFT);
    GET_NEXT_ITEM_SP (CHAIN_SIZELUFT);
    WordIndex luft_chrom = chain_seg_name_and_size (vb, gref, true, CHAIN_NAMELUFT, STRd(CHAIN_NAMELUFT), CHAIN_SIZELUFT, STRd(CHAIN_SIZELUFT), 
                                                    &mismatch, &verified);

    seg_by_did (VB, verified ? "1" : "0", 1, CHAIN_VERIFIED, 0);

    SEG_NEXT_ITEM_SP (CHAIN_STRNDLUFT);
    ASSERT_VALID_STRAND (CHAIN_STRNDLUFT);

    GET_NEXT_ITEM_SP (CHAIN_STARTLUFT);
    seg_pos_field (vb_, CHAIN_STARTLUFT, CHAIN_ENDLUFT, 0, 0, STRd(CHAIN_STARTLUFT), 0, field_len + 1);

    GET_NEXT_ITEM_SP (CHAIN_ENDLUFT);
    chain_seg_luft_end_field (vb, field_start, field_len);

    random_access_update_first_last_pos (vb_, 1, luft_chrom, STRd(CHAIN_STARTLUFT), STRd(CHAIN_ENDLUFT));

    // if ID is numeric, preferably store as delta, if not - normal snip
    GET_LAST_ITEM_SP (CHAIN_ID);
    if (str_is_int (STRd(CHAIN_ID)))
        seg_pos_field (vb_, CHAIN_ID, CHAIN_ID, 0, 0, STRd(CHAIN_ID), 0, field_len + 1); // just a numeric delta    
    else
        seg_by_did (VB, STRd(CHAIN_ID), CHAIN_ID, field_len + 1);

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
            seg_by_did (VB, &separator, 1, CHAIN_SEP, 0); // Chain format may have space or tab as a separator
            { SEG_NEXT_ITEM_SP (CHAIN_GAPS); }

            seg_by_did (VB, &separator, 1, CHAIN_SEP, 0); // space or tab
            { SEG_LAST_ITEM_SP (CHAIN_GAPS); }
        
        }
        // last line doesn't have DLUFT and DPRIM - seg the "standard" separator and then delete it 
        // (this way we don't ruin the all_the_same of CHAIN_SEP)
        else {
            seg_by_did (VB, &last_gap_sep, 1, CHAIN_SEP, 0); 
            seg_by_did (VB, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // prim_gap

            seg_by_did (VB, &last_gap_sep, 1, CHAIN_SEP, 0); 
            seg_by_did (VB, ((char[]){ SNIP_SPECIAL, CHAIN_SPECIAL_BACKSPACE }), 2, CHAIN_GAPS, 0); // luft_gap
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
                                 { .dict_id = { _CHAIN_SEP  } }, 
                                 { .dict_id = { _CHAIN_GAPS } }, // prim_gap
                                 { .dict_id = { _CHAIN_SEP  } }, 
                                 { .dict_id = { _CHAIN_GAPS } }, // luft_gap
                                 { .dict_id = { _CHAIN_EOL  } } }
    };

    container_seg (vb, CTX(CHAIN_SET), (ContainerP)&alignment_set, 0, 0, 0);

    // Empty line after alignment set
    GET_LAST_ITEM (EmptyLine);
    ASSSEG0 (!field_len, field_start, "Expecting an empty line after alignment set");
    SEG_EOL (CHAIN_EOL, false); 

    // case: alignment contig names mismatch the references - handling depends on --match-chrom-to-reference
    if (mismatch) {
        if (flag.match_chrom_to_reference) 
            goto rollback;
        else
            segconf.chain_mismatches_ref = true; // mark this chain file as unsuitable for use with --chain (we only ever set this, never reset, so no thread safety issues)
    }

    // in --match, we make sure to remove duplicate alignmets
    else if (!segconf.running && flag.match_chrom_to_reference &&
             chain_seg_is_duplicate ((ChainAlignment){ .prim_chrom      = prim_chrom,                          .luft_chrom      = luft_chrom,  // node_index==word_index bc ref contigs prepopulated in zctx(CHROM)
                                                       .prim_first_1pos = CTX(CHAIN_STARTPRIM)->last_value.i,  .luft_first_1pos = CTX(CHAIN_STARTLUFT)->last_value.i,
                                                       .prim_last_1pos  = CTX(CHAIN_ENDPRIM  )->last_value.i,  .luft_last_1pos  = CTX(CHAIN_ENDLUFT  )->last_value.i,
                                                       .is_xstrand      = *CHAIN_STRNDPRIM_str != *CHAIN_STRNDLUFT_str }))
        goto rollback;

    return next_field;

rollback:
    // remove the alignment completely
    seg_rollback (vb_); 
    vb->recon_size = save_recon_size - (next_field - field_start_line);
    vb->line_i--;

    return next_field;
}

//--------------
// PIZ functions
//--------------

// called after reconstructing the txt header and before compute threads
bool chain_piz_initialize (void)
{
    mutex_initialize (chain_mutex);

    if (flag.reading_chain) {
        contigs_build_contig_pkg_from_zctx (&prim_ctgs, ZCTX(CHAIN_NAMEPRIM), SORT_BY_NAME | SORT_BY_AC);
        contigs_build_contig_pkg_from_zctx (&luft_ctgs, ZCTX(CHAIN_NAMELUFT), SORT_BY_NONE);

        // add length from reference - same index as the contigs were copied from the reference to the CHAIN_NAMELUFT ctx during ZIP
        for_buf2 (Contig, ctg, i, luft_ctgs.contigs) 
            ctg->max_pos = ref_contigs_get_contig_length (gref, i, 0, 0, false);
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

#define CHAIN_SORT_BY(coord)                                        \
static int chain_sort_by_##coord (const void *a_, const void *b_)   \
{                                                                   \
    ChainAlignment *a = (ChainAlignment *)a_;                       \
    ChainAlignment *b = (ChainAlignment *)b_;                       \
                                                                    \
    if (a->coord##_chrom != b->coord##_chrom)                       \
        return a->coord##_chrom - b->coord##_chrom;                 \
                                                                    \
    /* note: don't use substraction - POS is 64b, return is 32b */  \
    if (a->coord##_first_1pos < b->coord##_first_1pos) return -1;   \
    if (b->coord##_first_1pos < a->coord##_first_1pos) return +1;   \
                                                                    \
    if (a->coord##_last_1pos  < b->coord##_last_1pos ) return -1;   \
    if (b->coord##_last_1pos  < a->coord##_last_1pos ) return +1;   \
                                                                    \
    return 0;                                                       \
}
CHAIN_SORT_BY(prim)
CHAIN_SORT_BY(luft)

static ASCENDING_SORTER (chain_sort_by_aln_i, ChainAlignment, aln_i)

static void chain_sort_mark_dups (void)
{
    ARRAY (ChainAlignment, aln, chain);

    qsort (STRb(chain), sizeof (ChainAlignment), chain_sort_by_prim);

    for (int32_t i=0; i < chain.len; i++) 
        aln[i].aln_i = i+1;

    qsort (STRb(chain), sizeof (ChainAlignment), chain_sort_by_luft);

    // mark duplicates alignments (can happen when chain file has both 1-> and chr-> and --match collapses them to the same contig)
    #define REMOVED 0x7fffffff // MAX_INT
    for (int32_t i=1; i < chain.len; i++) 
        if (aln[i].luft_chrom == aln[i-1].luft_chrom && aln[i].luft_first_1pos == aln[i-1].luft_first_1pos && aln[i].luft_last_1pos == aln[i-1].luft_last_1pos &&
            aln[i].prim_chrom == aln[i-1].prim_chrom && aln[i].prim_first_1pos == aln[i-1].prim_first_1pos && aln[i].prim_last_1pos == aln[i-1].prim_last_1pos)
            aln[i-1].aln_i = aln[i-1].luft_chrom = REMOVED;

    // sort again after de-dupping
    qsort (STRb(chain), sizeof (ChainAlignment), chain_sort_by_luft);

    // remove duplicates
    for (int32_t i=0; i < chain.len; i++) 
        if (aln[i].luft_chrom == REMOVED) {
            chain.len = i;
            break;
        }

    // mark overlapping alignment (only needed if we plan to display)
    if (flag.show_chain)
        for (int32_t i=1; i < chain.len; i++) 
            if (aln[i].luft_chrom == aln[i-1].luft_chrom && aln[i].luft_first_1pos <= aln[i-1].luft_last_1pos) {
                aln[i].overlap_aln_i   = aln[i-1].aln_i;
                aln[i-1].overlap_aln_i = aln[i].aln_i;
            }

    qsort (STRb(chain), sizeof (ChainAlignment), chain_sort_by_aln_i);
}

// initialize alignment set from alignment set header
static inline void chain_piz_filter_init_alignment_set (VBlockCHAIN *vb)
{
    reconstruct_from_ctx (vb, CHAIN_VERIFIED, 0, false); // sets last_int

    vb->next_luft_0pos = vb->last_int(CHAIN_STARTLUFT); 
    vb->next_prim_0pos = vb->last_int(CHAIN_STARTPRIM); 
}

// store prim_gap value in local.param, as last_value.i will be overridden by luft_gap, as they share the same context
static inline void chain_piz_filter_save_prim_gap (VBlockP vb)
{
    CTX(CHAIN_GAPS)->local.param = vb->last_int(CHAIN_GAPS); 
}

// ingest one alignment
static inline void chain_piz_filter_ingest_alignmet (VBlockCHAIN *vb)
{
    if (!vb->last_int (CHAIN_VERIFIED)) goto done; // a line Seg couldn't verify - don't ingest

    mutex_lock (chain_mutex);

    buf_alloc (NULL, &chain, 1, 0, ChainAlignment, 2, "chain");  // initial allocation is in chain_load

    int64_t size = vb->last_int(CHAIN_SIZE);

    ASSINP (*last_txt(VB, CHAIN_STRNDPRIM) == '+', "Chain file contains alignments with tStrand=\"%c\" (strand of Primary reference), this is not support by Genozip.",
            *last_txt(VB, CHAIN_STRNDLUFT));

    bool is_xstrand = (*last_txt(VB, CHAIN_STRNDLUFT) == '-');
    PosType last_prim_0pos = vb->next_prim_0pos + size - 1; // last POS of the dst alignment, 0-based
    PosType last_luft_0pos = vb->next_luft_0pos + size - 1; // last POS of the dst alignment, 0-based
    PosType size_dst = vb->last_int (CHAIN_SIZELUFT);  // contig LN

    // note on negative strand: the source region (the source always has a positive strand), is aligned to a region on the destination
    // contig, but the positions given are counting starting from end of the contig (the last base of the contig is 0), and the sequence
    // is complemented. For our chain alignment representation, we convert the coordinates to positive-strand terms, and set XSTRAND=X
    // to record that this is a reverse strand

    BNXT (ChainAlignment, chain) = (ChainAlignment){ 
        .prim_chrom      = vb->last_index(CHAIN_NAMEPRIM), // same index as chain contigs must be --match-chroms to be used for --chain
        .prim_first_1pos = 1 + vb->next_prim_0pos,         // +1 bc our alignments are 1-based vs the chain file that is 0-based
        .prim_last_1pos  = 1 + last_prim_0pos,
        .luft_chrom      = vb->last_index(CHAIN_NAMELUFT), // ditto
        .luft_first_1pos = is_xstrand ? (size_dst - last_luft_0pos)     : 1 + vb->next_luft_0pos, // +1 to convert to 1-base. note that "size_dst" has a built-in +1.
        .luft_last_1pos  = is_xstrand ? (size_dst - vb->next_luft_0pos) : 1 + last_luft_0pos,

        .is_xstrand     = is_xstrand
    };

    Context *gaps_ctx   = CTX(CHAIN_GAPS);
    vb->next_prim_0pos += size + gaps_ctx->local.param;  // prim_gap
    vb->next_luft_0pos += size + gaps_ctx->last_value.i; // luft_gap

    mutex_unlock (chain_mutex);

done:
    vb->drop_curr_line = "chain";
}

// verify that adding up all alignments and gaps, results in the end position specified in the header
static inline void chain_piz_filter_verify_alignment_set (VBlockCHAIN *vb)
{
    ASSINP (vb->next_prim_0pos == vb->last_int(CHAIN_ENDPRIM),
            "%.*s\nBad data ^^^ in chain file %s: Expecting alignments to add up to ENDPRIM=%s, but they add up to %s",
            (int)(vb->txt_data.len - vb->line_start), Bc (vb->txt_data, vb->line_start),
            z_name, str_int_commas (vb->last_int(CHAIN_ENDPRIM)).s, str_int_commas (vb->next_prim_0pos).s);

    ASSINP (vb->next_luft_0pos == vb->last_int(CHAIN_ENDLUFT),
            "%.*sBad data ^^^ in chain file %s: Expecting alignments to add up to ENDLUFT=%s, but they add up to %s",
            (int)(vb->txt_data.len - vb->line_start), Bc (vb->txt_data, vb->line_start),
            z_name, str_int_commas (vb->last_int(CHAIN_ENDLUFT)).s, str_int_commas (vb->next_luft_0pos).s);
}

// set contig lengths according to tSize and qSize
static inline void chain_piz_filter_add_contig_length (VBlockP vb)
{
    mutex_lock (chain_mutex);

    Contig *prim_ctg = B(Contig, prim_ctgs.contigs, vb->last_index(CHAIN_NAMEPRIM));
    Contig *luft_ctg = B(Contig, luft_ctgs.contigs, vb->last_index(CHAIN_NAMELUFT));

    prim_ctg->max_pos = vb->last_int(CHAIN_SIZEPRIM);
    luft_ctg->max_pos = vb->last_int(CHAIN_SIZELUFT);

    mutex_unlock (chain_mutex);
}

CONTAINER_FILTER_FUNC (chain_piz_filter)
{
    VBlockCHAIN *vb_ = (VBlockCHAIN *)vb;

    // ingest alignments (and report Chain file data issues) only when consuming with --chain, not when merely pizzing
    if (flag.reading_chain) { 

        // before alignment-set first EOL and before alignments - initialize next_luft_0pos and next_prim_0pos
        if (dict_id.num == _CHAIN_TOPLEVEL && item == 13) 
            chain_piz_filter_init_alignment_set (vb_);

        // save prim_gap before reconstructing luft_gap (prim_gap was just processed)
        else if (dict_id.num == _CHAIN_SET && item == 4) 
            chain_piz_filter_save_prim_gap (vb);

        // before EOF of each alignment, ingest alignment (only if it passed verification during Seg)
        else if (dict_id.num == _CHAIN_SET && item == 5) 
            chain_piz_filter_ingest_alignmet (vb_);

        // before alignment-set second EOL and after alignments - verify that numbers add up, and also set contigs
        else if (dict_id.num == _CHAIN_TOPLEVEL && item == 15 && vb->last_int (CHAIN_VERIFIED)) {
            chain_piz_filter_verify_alignment_set (vb_);
            chain_piz_filter_add_contig_length (vb);
        }
    }
    
    return true;    
}

// main thread: ZIP with --chain, reading txt header: returns null-terminated string of contig, or NULL if contig_i is out of range
rom chain_get_luft_contig (uint32_t contig_i, PosType *length)
{
    if (contig_i >= luft_ctgs.contigs.len) return NULL;

    const Contig *ctg = B(Contig, luft_ctgs.contigs, contig_i);

    *length = ctg->max_pos;
    return Bc (luft_ctgs.dict, ctg->char_index);
}

uint64_t chain_get_num_prim_contigs (void)
{
    return prim_ctgs.contigs.len;
}

static void chain_contigs_show (bool is_primary)
{
    ContextP zctx = ZCTX (is_primary ? CHAIN_NAMEPRIM : CHAIN_NAMELUFT);
    ASSERTISALLOCED (zctx->counts); // will fail with chain files compressed prior to 12.0.35

    ConstBufferP contigs = is_primary ? &prim_ctgs.contigs : &luft_ctgs.contigs;
    ARRAY (const Contig, cn, *contigs);

    iprintf ("\n%s chain file contigs:\n", is_primary ? "PRIMARY" : "LUFT");
    for (uint32_t i=0; i < contigs->len; i++) 
        if (*B64(zctx->counts, i)) // counts>0 means contig (copied from reference) appears in chain file
            iprintf ("%s index=%-2u %s length=%"PRId64" %s\n", is_primary ? "PRIMARY" : "LUFT", i, 
                     ctx_get_words_snip (zctx, i), cn[i].max_pos, display_acc_num (&cn[i].metadata.parsed.ac).s);
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
    
    bool chain_not_for_use = flag.show_chain || flag.show_chain_contigs;
    
    ASSINP (Z_DT(CHAIN), "expected %s to be a genozip'ed chain file, but its a %s file. Tip: compress the chain with \"genozip --input chain\"", 
            z_name, dt_name (z_file->data_type));

    ASSINP (chain_not_for_use || (ref_get_filename (gref) && ref_get_filename (prim_ref)),
            "%s is unsuitable for use with %s because it was not compressed with source and destination references using genozip --reference. See: " WEBSITE_CHAIN, 
            z_name, OT("chain", "C"));

    ASSINP (chain_not_for_use || !z_file->z_flags.dts_mismatch, 
            "%s is unsuitable for use with %s because its contigs mismatch those of the references. To fix, recreate it using the --match-chrom-to-reference option. See: " WEBSITE_CHAIN,
            z_name, OT("chain", "C"));

    z_file->basename = file_basename (flag.reading_chain, false, "(chain-file)", NULL, 0);

    TEMP_VALUE (command, PIZ);

    if (flag.reference) flag.reference = REF_LIFTOVER;

    flag.no_writer = flag.no_writer_thread = true;
    flag.genocat_no_reconstruct = false;
    flag.genocat_global_area_only = false;
    
    Dispatcher dispachter = piz_z_file_initialize();

    // load both references, now that it is set (either explicitly from the command line, or implicitly from the chain GENOZIP_HEADER)
    SAVE_VALUE (z_file); // actually, read the references first
    ref_load_external_reference (gref, NULL);
    ref_load_external_reference (prim_ref, NULL);
    RESTORE_VALUE (z_file);

    // test for matching MD5 between external references and reference in the chain file header (doing it here, because reference is read after chain file, if its explicitly specified)
    digest_verify_ref_is_equal (gref, header.ref_filename, header.ref_file_md5);
    digest_verify_ref_is_equal (prim_ref, header.chain.prim_filename, header.chain.prim_file_md5);

    flag.quiet = true; // don't show progress indicator for the chain file - it is very fast 
    flag.maybe_vb_modified_by_reconstructor = true; // we drop all the lines

    piz_one_txt_file (dispachter, false, false, COMP_NONE);

    // --show-chain-contigs
    if (flag.show_chain_contigs) {
        chain_contigs_show (true);
        chain_contigs_show (false);
        if (is_genocat) exit_ok();  // in genocat this, not the data
    }

    // sort the alignmants by (prim_ref_index, prim_start)
    chain_sort_mark_dups();
        
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
    buf_destroy (chain);
    contigs_destroy (&prim_ctgs);
    contigs_destroy (&luft_ctgs);
    mutex_destroy (chain_mutex);
    chain_filename = NULL; 
}

// binary search for the first alignment of a src contig
static ChainAlignment *chain_get_first_aln (WordIndex prim_ref_index, int32_t start, int32_t end) 
{
    if (end < start) return NULL; // prim_contig doesn't exist in chain file

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = B(ChainAlignment, chain, mid);

    // case: success - same prim_contig, and previous aln has a different prim_contig 
    if (aln->prim_chrom == prim_ref_index && (!mid || (aln-1)->prim_chrom != prim_ref_index))
        return aln;

    // case: prim_contig less than aln's, OR is the same, but aln isn't the first with this prim_contig - search lower half
    if (prim_ref_index <= aln->prim_chrom)
        return chain_get_first_aln (prim_ref_index, start, mid-1);
    
    // case: prim_contig is more than aln's - search upper half
    else
        return chain_get_first_aln (prim_ref_index, mid+1, end);
}

// append luft_contigs buffer with all dst chroms indicies for which there exists alignment with prim_chrom
void chain_append_all_luft_ref_index (rom prim_contig_name, unsigned prim_contig_name_len, PosType LN, BufferP luft_contigs)
{
    WordIndex prim_ref_index = contigs_get_matching (&prim_ctgs, STRa(prim_contig_name), LN, false, NULL);
    WordIndex prev_dst = WORD_INDEX_NONE;

    if (prim_ref_index == WORD_INDEX_NONE) return; // this contig is not in the chain file

    for (ChainAlignment *aln = chain_get_first_aln (prim_ref_index, 0, chain.len-1);
         aln && aln < BLST (ChainAlignment, chain) && aln->prim_chrom == prim_ref_index;
         aln++)
         
         if (aln->luft_chrom != prev_dst) { // note: this prevents consecutive duplicates, but not non-consecutive duplicates
            buf_alloc (NULL, luft_contigs, 100, 1, WordIndex, 2, NULL);
            BNXT (WordIndex, *luft_contigs) = prev_dst = aln->luft_chrom;
         }
}

// get luft_contig, prim_pos from prim_contig, luft_pos (binary search on chain)
// returns true if successful
static bool chain_get_liftover_coords_do (WordIndex prim_ref_index, PosType prim_1pos, 
                                          int32_t start, int32_t end,
                                          WordIndex *luft_ref_index, PosType *luft_1pos, bool *is_xstrand, uint32_t *aln_i) // out
{
    // case: no mapping
    if (end < start) {
        if (luft_ref_index) *luft_ref_index = WORD_INDEX_NONE;
        if (luft_1pos) *luft_1pos = 0;
        if (is_xstrand) *is_xstrand = false;
        return false;
    }

    uint32_t mid = (start + end) / 2;
    ChainAlignment *aln = B(ChainAlignment, chain, mid);

    // case: success
    if (aln->prim_chrom == prim_ref_index && 
        aln->prim_first_1pos <= prim_1pos && 
        aln->prim_last_1pos >= prim_1pos) {
            if (luft_ref_index) *luft_ref_index = aln->luft_chrom;
            
            PosType pos_offset = prim_1pos - aln->prim_first_1pos; // offset of POS from the beginning of the src (src is always positive strand)
            if (luft_1pos) *luft_1pos = aln->is_xstrand ? (aln->luft_last_1pos  - pos_offset)  // dst is reverse strand - offset is from the end of the alignment, going back
                                                      : (aln->luft_first_1pos + pos_offset); // dst is positive strand - offset is from the start of the alignment
            
            if (is_xstrand) *is_xstrand = aln->is_xstrand;
            if (aln_i) *aln_i = mid;
            return true;
        }

    // case: prim_contig, luft_pos is less than aln - search lower half
    if (prim_ref_index < aln->prim_chrom ||
        (prim_ref_index == aln->prim_chrom && prim_1pos < aln->prim_first_1pos))
        return chain_get_liftover_coords_do (prim_ref_index, prim_1pos, start, mid-1, luft_ref_index, luft_1pos, is_xstrand, aln_i);
    
    // case: prim_contig,luft_pos is more than aln - search upper half
    else
        return chain_get_liftover_coords_do (prim_ref_index, prim_1pos, mid+1, end, luft_ref_index, luft_1pos, is_xstrand, aln_i);
}

bool chain_get_liftover_coords (WordIndex prim_ref_index, PosType prim_1pos, 
                                WordIndex *luft_ref_index, PosType *luft_1pos, bool *is_xstrand, uint32_t *aln_i) // out
{
    return chain_get_liftover_coords_do (prim_ref_index, prim_1pos, 0, chain.len-1, luft_ref_index, luft_1pos, is_xstrand, aln_i);
}

PosType chain_get_aln_gap_after (uint32_t aln_i)
{
    const ChainAlignment *aln = B(const ChainAlignment, chain, aln_i);

    // case: this is the last alignment of this contig
    if (aln_i == chain.len - 1 || aln->prim_chrom != (aln+1)->prim_chrom) { 
        const Contig *prim_ctg = CONTIG (prim_ctgs, aln->prim_chrom);
        return prim_ctg->max_pos - aln->prim_last_1pos;
    }
    else
        return (aln+1)->prim_first_1pos - aln->prim_last_1pos - 1;
}

PosType chain_get_aln_prim_last_pos (uint32_t aln_i)
{
    return B(const ChainAlignment, chain, aln_i)->prim_last_1pos;
}
