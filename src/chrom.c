// ------------------------------------------------------------------
//   chrom.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <stdarg.h>
#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "segconf.h"
#include "buffer.h"
#include "context.h"
#include "file.h"
#include "random_access.h"
#include "chrom.h"
#include "zfile.h"
#include "endianness.h"
#include "contigs.h"
#include "reference.h"
#include "version.h"

static ContextP sorter_ctx = NULL; 
static Buffer chrom_sorter = {}; // ZIP/PIZ: index into sorter_ctx->nodes/word_list, sorted alphabetically by snip

// the data of SEC_CHROM2REF_MAP - this is part of the Genozip file format
typedef struct __attribute__ ((__packed__)) { WordIndex chrom_index, ref_index; } Chrom2Ref; 

//-------------------------
// chrom2ref mapping stuff
//-------------------------

// ZIP of a file with an external reference: 
// The CHROM dictionary includes: first - the txtheader contig names, then reference file contig names that are not in txtheader, 
// and finally names encountered in the file that are not in either the header or the reference (see ctx_populate_zf_ctx_from_contigs()). 
// Since the CHROM context is prepopulated from the txtheader and the reference, often not all the entries are used by the file data.
//
// Here, we create a mapping between those chrom entries that are used (count>0) and the index of the contig in the reference file
// against which this chrom data was compressed. We rely on contigs_get_matching() returning the same result here as it did in Seg.

void chrom_2ref_compress (Reference ref)
{
    if (flag.show_chrom2ref) 
        iprint0 ("\nAlternative chrom indices (output of --show-chrom2ref): chroms that are in the file and are mapped to a different index in the reference\n");
    
    ASSERTNOTINUSE (evb->scratch);
    buf_alloc (evb, &evb->scratch, 0, ZCTX(CHROM)->chrom2ref_map.len, Chrom2Ref, 1, "scratch");
    ContextP zctx = ZCTX(DTFZ(prim_chrom));

    for_buf2 (WordIndex, ref_index, chrom_node_index, ZCTX(CHROM)->chrom2ref_map) {
        if (flag.show_chrom2ref) {
            rom chrom_name = ctx_snip_from_zf_nodes (zctx, chrom_node_index, 0, 0);
            rom ref_name = *ref_index >= 0 ? ref_contigs_get_name (ref, *ref_index, NULL) : "(none)";

            if (*ref_index != WORD_INDEX_NONE) 
                iprintf ("In file: '%s' (%d)\tIn reference: '%s' (%d)\t%s\n", 
                         chrom_name, chrom_node_index, ref_name, *ref_index, chrom_node_index != *ref_index ? "INDEX_CHANGE" : "");
            else
                iprintf ("In file: '%s' (%d)\tNot in reference\n", chrom_name, chrom_node_index);
        }

        // adds the mapping if not identify and adds -1 if this chrom doesn't map to a ref contig.
        // note: we add only contigs that are used (count>0) except for aligner_available in which case we don't have counts (for REF_EXTERNAL, we have 
        // populated all contigs in zip_initialize, and for REF_EXT_STORE we add contigs with any bit set in is_set) 
        if (*ref_index != chrom_node_index && *ref_index != WORD_INDEX_NONE && (*B64(zctx->counts, chrom_node_index) || flag.aligner_available))
            BNXT (Chrom2Ref, evb->scratch) = (Chrom2Ref){ .chrom_index = BGEN32(chrom_node_index), 
                                                          .ref_index   = BGEN32(*ref_index)       };
    }

    if (evb->scratch.len) {
        evb->scratch.len *= sizeof (Chrom2Ref);
        zfile_compress_section_data_ex (evb, NULL, SEC_CHROM2REF_MAP, &evb->scratch, 0,0, CODEC_LZMA, SECTION_FLAGS_NONE, NULL); // compresses better with LZMA than BZLIB
    }

    buf_free (evb->scratch);
}

void chrom_2ref_load (Reference ref)
{
    Section sec = sections_last_sec (SEC_CHROM2REF_MAP, SOFT_FAIL);
    if (!sec) return;

    zfile_get_global_section (SectionHeader, sec, &evb->scratch, "scratch");

    if (flag.show_chrom2ref) 
        iprint0 ("\nAlternative chrom indices (output of --show-chrom2ref): chroms that are in the txt file and are mapped to a different index in the reference\n");

    evb->scratch.len /= sizeof (Chrom2Ref);
    Context *zctx = ZCTX(CHROM);

    // create mapping user index -> reference index
    buf_alloc (evb, &ZCTX(CHROM)->chrom2ref_map, 0, zctx->word_list.len, WordIndex, 1, "ZCTX(CHROM)->chrom2ref_map");
    ZCTX(CHROM)->chrom2ref_map.len = zctx->word_list.len;

    // initialize with unity mapping
    ARRAY (WordIndex, map, ZCTX(CHROM)->chrom2ref_map);
    for (uint32_t i=0; i < zctx->word_list.len32; i++)
        map[i] = i;

    // the indices of chroms that are NOT in the reference (they are only in the user file), will be mapped to ref chroms
    ConstContigPkgP ctgs = ref_get_ctgs (ref); 
    WordIndex num_ref_contigs = ctgs->contigs.len; // must be signed int

    for_buf2 (Chrom2Ref, ent, i, evb->scratch) {
        WordIndex chrom_index = BGEN32 (ent->chrom_index);
        WordIndex ref_index   = BGEN32 (ent->ref_index);

        ASSERT (chrom_index >= 0 && chrom_index < zctx->word_list.len, "chrom_index=%d out of range [0,%d]", chrom_index, (int32_t)zctx->word_list.len-1);
        ASSERT (!num_ref_contigs /* ref not loaded */ || (ref_index >= -1 && ref_index < num_ref_contigs), 
                "ref_index=%d out of range [-1,%u] (chrom_index=%u i=%u len=%u)", 
                ref_index, num_ref_contigs-1, chrom_index, i, evb->scratch.len32);

        map[chrom_index] = ref_index;

        if (flag.show_chrom2ref) {
            rom chrom_name = ctx_get_words_snip (zctx, chrom_index);
            rom ref_name   = ref_index >= 0 ? ref_contigs_get_name (ref, ref_index, NULL) : NULL;
            if (ref_name)
                iprintf ("In file: '%s' (%d)\tIn reference: '%s' (%d)\n", chrom_name, chrom_index, ref_name, ref_index);
            else
                iprintf ("In file: '%s' (%d)\tNot in reference\n", chrom_name, chrom_index);
        }
    }

    if (flag.show_chrom2ref && is_genocat) exit (EXIT_OK); // in genocat this, not the data

    buf_free (evb->scratch);
}

// ZIP: returns the ref index by the chrom index, works only after Segging of CHROM
WordIndex chrom_2ref_seg_get (Reference ref, ConstVBlockP vb, WordIndex chrom_index)
{ 
    ASSSEG (chrom_index >= WORD_INDEX_NONE, "invalid chrom_index=%d", chrom_index);

    int32_t ol_len = vb->ol_chrom2ref_map.len32;
    decl_const_ctx(CHROM);

    WordIndex ref_index = (chrom_index == WORD_INDEX_NONE)                  ? WORD_INDEX_NONE
                        : (chrom_index < ol_len)                            ? *B(WordIndex, vb->ol_chrom2ref_map, chrom_index)
                        : (chrom_index < ol_len + ctx->chrom2ref_map.len32) ? *B(WordIndex, ctx->chrom2ref_map, chrom_index - ol_len)
                        :                                                     WORD_INDEX_NONE; // not in reference

    ASSSEG (ref_index >= WORD_INDEX_NONE && ref_index < (WordIndex)ref_num_contigs (ref), 
            "ref_index=%d out of range: ref->ranges.len=%u, chrom_index=%d", ref_index, ref_num_contigs (ref), chrom_index);

    return ref_index;
}   

void chrom_calculate_ref2chrom (uint64_t num_ref_contigs)
{
    buf_alloc_255 (evb, &z_file->ref2chrom_map, 0, num_ref_contigs, WordIndex, 0, "ref2chrom_map");
    z_file->ref2chrom_map.len = num_ref_contigs;

    ARRAY (WordIndex, r2c, z_file->ref2chrom_map);
    ARRAY (WordIndex, c2r, ZCTX(CHROM)->chrom2ref_map);
    
    for (unsigned i=0; i < c2r_len; i++)
        if (c2r[i] != WORD_INDEX_NONE)
            r2c[c2r[i]] = i;
}

//-------------
// Seg stuff
//-------------

WordIndex chrom_seg_ex (VBlockP vb, Did did_i, 
                        STRp(chrom), 
                        PosType64 LN,     // Optional, if readily known
                        bool *is_alt_out, // need iff flag.match_chrom_to_reference.
                        int add_bytes,    // must be signed
                        bool recon_changes_if_match, // whether reconstruction changes in case of change in chrom name due to --match-chrom
                        bool *is_new_out) // optional out
{
    ASSERTNOTZERO (chrom_len);
    decl_ctx (did_i);
    bool is_primary = did_i == DTF(prim_chrom); // note: possibly neither primary or luft, eg SA_RNAME
    bool is_luft    = did_i == DTF(luft_chrom);
    bool has_chain  = chain_is_loaded || IS_REF_MAKE_CHAIN;

    WordIndex chrom_node_index = WORD_INDEX_NONE, ref_index = WORD_INDEX_NONE;
    int32_t chrom_name_growth=0;
    bool is_new, is_alt=false;
    rom save_chrom = chrom;
    uint32_t save_chrom_len = chrom_len;
    
    Reference ref = !(flag.reference & REF_ZIP_LOADED) ? NULL
                  : (has_chain && is_primary)          ? prim_ref // DT_VCF --chain or DT_CHAIN 
                  :                                      gref;
    
    // case match_chrom_to_reference: rather than segging the chrom as in the txt_data, we seg the matching name in the reference, if there is one.
    if (flag.match_chrom_to_reference && 
        (!is_luft || has_chain)) { // note: can't match the Luft reference of a DVCF (can only match Luft with --chain)

        // case: chrom is the same as previous line - use cached (note: might be different than ctx->last_snip since we match_chrom_to_reference)
        #define last_growth last_delta
        if (vb->line_i && chrom_len == ctx->last_txt.len && ctx->last_txt.index != INVALID_LAST_TXT_INDEX && 
            !memcmp (chrom, last_txtx(vb, ctx), chrom_len)) {
            if (is_alt_out) *is_alt_out = ctx->last_is_alt;
            if (is_new_out) *is_new_out = false;
            vb->recon_size += ctx->last_growth;
            
            return seg_duplicate_last (vb, ctx, add_bytes + ctx->last_growth);
        }

        // test for ref match
        ref_index = ref_contigs_get_matching (ref, LN, STRa(chrom), pSTRa(chrom), false, &is_alt, &chrom_name_growth); 
        if (ref_index != WORD_INDEX_NONE) { // chrom/chrom_len are modified in-place
            if (is_alt_out) *is_alt_out = is_alt;

            if (!recon_changes_if_match) chrom_name_growth = 0; // in BAM, chroms are stored as a ref_id (=chrom_index), so no change in reconstruction size

            vb->recon_size += chrom_name_growth;
            
            chrom_node_index = seg_by_ctx_ex (vb, STRa(chrom), ctx, add_bytes + chrom_name_growth, &is_new); // seg modified chrom
        }

        // update cache
        seg_set_last_txt (vb, ctx, STRa(save_chrom));
        ctx->last_is_alt = is_alt;
        ctx->last_growth = chrom_name_growth;  
        ctx->no_stons    = true; // needed for seg_duplicate_last

        // if a match was found, we're done
        if (ref_index != WORD_INDEX_NONE) goto finalize;
    }

    // case: either without --match-chrom-to-reference OR chrom not found in the reference
    chrom_node_index = seg_by_ctx_ex (vb, STRa(chrom), ctx, add_bytes, &is_new); // note: this is not the same as ref_index, bc ctx->nodes contains the header contigs first, followed by the reference contigs that are not already in the header
    
    STR0 (ref_contig);
    if (is_new && ref)
        ref_index = ref_contigs_get_matching (ref, LN, STRa(chrom), pSTRa(ref_contig), false, &is_alt, NULL);

    // warn if the file's contigs are called by a different name in the reference (eg 22/chr22)
    static bool once[2]={};
    if (ref_index != WORD_INDEX_NONE && is_alt && // a new chrom that matched to the reference with an alternative name
        (is_primary || is_luft) &&
        !segconf.running  &&              // segconf runs with flag.quiet so the user won't see the warning
        !flag.match_chrom_to_reference && // we didn't already attempt to match to the reference
        !__atomic_test_and_set (&once[is_primary], __ATOMIC_RELAXED)) {  // skip if we've shown the warning already
            
            if (VB_DT(CHAIN))
                WARN ("Warning: Contig name mismatch between %s and reference file %s. For example: %s file: \"%.*s\" Reference file: \"%.*s\". "
                      "You may use --match-chrom-to-reference to create %s with contigs matching those of the reference. This is required, if this file is for use with 'genozip --chain'. More info: %s\n",
                      txt_name, ref_get_filename (ref), dt_name (vb->data_type), STRf(chrom), STRf(ref_contig), z_name, WEBSITE_CHAIN);
            else
                WARN ("FYI: Contigs name mismatch between %s and reference file %s. For example: %s file: \"%.*s\" Reference file: \"%.*s\". "
                      "You may use --match-chrom-to-reference to create %s with contigs matching those of the reference. This makes no difference for the compression. More info: %s",
                      txt_name, ref_get_filename (ref), dt_name (vb->data_type), STRf(chrom), STRf(ref_contig), z_name, WEBSITE_MATCH_CHROM);
    } // we don't use WARN_ONCE bc we want the "once" to also include ref_contigs_get_matching

    if (is_alt_out) *is_alt_out = false;

finalize:
    if (is_new_out) *is_new_out = is_new;        

    if (is_primary || is_luft)
        random_access_update_chrom (vb, !is_primary, chrom_node_index, STRa(chrom)); 

    if (is_primary) {
        vb->chrom_node_index = chrom_node_index;

        if (chrom_node_index != WORD_INDEX_NONE) 
            STRset (vb->chrom_name, chrom);
    }

    if (is_new) { 
        buf_alloc_255 (vb, &ctx->chrom2ref_map, 0, chrom_node_index+1, WordIndex, CTX_GROWTH, "chrom2ref_map");
        *B(WordIndex, ctx->chrom2ref_map, chrom_node_index) = ref_index; // note: not a simple BNXT bc possibly nodes can be created outside of chrom_seg_ex, eg SNIP_SPECIAL
    }

    return chrom_node_index;
}

WordIndex chrom_seg_no_b250 (VBlockP vb, STRp(chrom_name), bool *is_new)
{
    WordIndex chrom_node_index = chrom_seg_ex (VB, CHROM, STRa(chrom_name), 0, NULL, 0, false, is_new); // also adds to random access etc
    ctx_decrement_count (VB, CTX(CHROM), chrom_node_index);
    CTX(CHROM)->b250.len--;

    return chrom_node_index;
}

bool chrom_seg_cb (VBlockP vb, ContextP ctx, STRp (chrom), uint32_t repeat)
{
    chrom_seg_ex (vb, ctx->did_i, STRa(chrom), 0, NULL, chrom_len, true, NULL);

    return true; // segged successfully
}

static SORTER (chrom_create_zip_sorter)
{
    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    CtxNode *word_a = B(CtxNode, sorter_ctx->nodes, index_a);
    CtxNode *word_b = B(CtxNode, sorter_ctx->nodes, index_b);

    return strcmp (Bc (sorter_ctx->dict, word_a->char_index),  
                   Bc (sorter_ctx->dict, word_b->char_index));
}

static SORTER (chrom_create_piz_sorter)
{
    uint32_t index_a = *(uint32_t *)a;
    uint32_t index_b = *(uint32_t *)b;

    CtxWord *word_a = B(CtxWord, sorter_ctx->word_list, index_a);
    CtxWord *word_b = B(CtxWord, sorter_ctx->word_list, index_b);
    
    return strcmp (Bc (sorter_ctx->dict, word_a->index),
                   Bc (sorter_ctx->dict, word_b->index));
}

// ZIP/PIZ MUST be run by the main thread only
void chrom_index_by_name (Did chrom_did_i)
{
    sorter_ctx = ZCTX(chrom_did_i);
    uint32_t num_words = (command==ZIP) ? sorter_ctx->nodes.len32 : sorter_ctx->word_list.len32;

    buf_free (chrom_sorter);

    // chrom_sorter - an array of uint32 of indexes into ZCTX(CHROM)->word_list - sorted by alphabetical order of the snip in ZCTX(CHROM)->dict
    buf_alloc (evb, &chrom_sorter, 0, num_words, uint32_t, 1, "chrom_sorter");
    
    if (IS_ZIP) {
        for_buf2 (CtxNode, node, i, sorter_ctx->nodes)
            if (node->snip_len) BNXT32(chrom_sorter) = i;
    }
    else
        for_buf2 (CtxWord, word, i, sorter_ctx->word_list)
            if (word->len) BNXT32(chrom_sorter) = i;

    qsort (STRb(chrom_sorter), sizeof(uint32_t), IS_ZIP ? chrom_create_zip_sorter : chrom_create_piz_sorter);
}

// binary search for this chrom in ZCTX(CHROM). we count on gcc tail recursion optimization to keep this fast.
static WordIndex chrom_zip_get_by_name_do (rom chrom_name, WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    STR (snip);
    WordIndex node_index = *B(WordIndex, chrom_sorter, mid_sorted_index);
    ctx_snip_from_zf_nodes (ZCTX(CHROM), node_index, pSTRa(snip));

    int cmp = strcmp (snip, chrom_name);
    if (cmp < 0) return chrom_zip_get_by_name_do (chrom_name, mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return chrom_zip_get_by_name_do (chrom_name, first_sorted_index, mid_sorted_index-1);

    return node_index;
}             

// binary search for this chrom in ZCTX(CHROM). we count on gcc tail recursion optimization to keep this fast.
static WordIndex chrom_piz_get_by_name_do (STRp (chrom_name), WordIndex first_sorted_index, WordIndex last_sorted_index)
{
    if (first_sorted_index > last_sorted_index) return WORD_INDEX_NONE; // not found

    WordIndex mid_sorted_index = (first_sorted_index + last_sorted_index) / 2;
    
    STR (snip);
    WordIndex word_index = *B(WordIndex, chrom_sorter, mid_sorted_index);
    ctx_get_snip_by_word_index (ZCTX(CHROM), word_index, snip);

    int cmp = strncmp (snip, chrom_name, chrom_name_len);
    if (!cmp && snip_len != chrom_name_len) // identical prefix but different length
        cmp = snip_len - chrom_name_len;

    if (cmp < 0) return chrom_piz_get_by_name_do (STRa(chrom_name), mid_sorted_index+1, last_sorted_index);
    if (cmp > 0) return chrom_piz_get_by_name_do (STRa(chrom_name), first_sorted_index, mid_sorted_index-1);

    return word_index;
}             

// note: search within the partial set of chroms that existed when chrom_index_by_name was called.
WordIndex chrom_get_by_name (STRp (chrom_name))
{
    if (!chrom_sorter.len) return WORD_INDEX_NONE;

    WordIndex wi;
    SAFE_NULT(chrom_name); 
    
    if (IS_ZIP) wi =  chrom_zip_get_by_name_do (chrom_name, 0, chrom_sorter.len-1); // not necessarily all of CHROM, just chrom_sorter.len
    else                wi =  chrom_piz_get_by_name_do (STRa(chrom_name), 0, chrom_sorter.len-1);
    
    SAFE_RESTORE;
    return wi;
}

void chrom_finalize (void)
{
    buf_destroy (chrom_sorter);
}
