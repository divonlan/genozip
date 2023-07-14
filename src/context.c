// ------------------------------------------------------------------
//   context.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include <stdarg.h>
#include "genozip.h"
#include "profiler.h"
#include "sections.h"
#include "vcf.h"
#include "vblock.h"
#include "context.h"
#include "zfile.h"
#include "endianness.h"
#include "file.h"
#include "hash.h"
#include "seg.h"
#include "dict_id.h"
#include "reference.h"
#include "mutex.h"
#include "progress.h"
#include "dispatcher.h"
#include "compressor.h"
#include "strings.h"
#include "codec.h"
#include "flags.h"
#include "zip.h"
#include "piz.h"
#include "reconstruct.h"
#include "contigs.h"
#include "segconf.h"
#include "chrom.h"
#include "buffer.h"
#include "version.h"
#include "qname.h"
#include "threads.h"
#include "aliases.h"

#define INITIAL_NUM_NODES 10000

// ZIP only: set a bit in counts of a node, so it is not removed by ctx_shorten_unused_dict_words
#define COUNT_PROTECTED_FROM_REMOVAL64 0x8000000000000000ULL 
#define COUNT_PROTECTED_FROM_REMOVAL32 0x80000000

#define EXCESSIVE_DICT_SIZE (512 MB) // warn if dict size goes beyond 512M

#define ZMUTEX(zctx) z_file->ctx_mutex[(zctx)->did_i] // ZIP only

// inserts dict_id->did_i to the map, if one of the two entries is available
static inline void set_d2d_map (DictIdtoDidMap d2d_map, DictId dict_id, Did did_i)
{
    // thread safety for z_file d2d_map: we don't bother with having a mutex, in the worst case scenario, two threads will test an entry
    // as empty and then both write to it, with one of them prevailing. that's fine (+ Likely its the same did_i anyway).

    if (d2d_map[dict_id.map_key[0]] == DID_NONE)    // d2d_map entry is free
        d2d_map[dict_id.map_key[0]] = did_i;

    else if (d2d_map[dict_id.map_key[0]] == did_i)  // already has requested value - nothing to do 
        {}
    
    else if (d2d_map[ALT_KEY(dict_id)] == DID_NONE) // fallback entry is free or we can override it
        d2d_map[ALT_KEY(dict_id)] = did_i;
}

// ZIP: add a snip to the dictionary the first time it is encountered in the txt file.
// the dictionary will be written to z_file.
typedef enum { DICT_VB, DICT_ZF, DICT_ZF_SINGLETON } DictType;
static inline CharIndex ctx_insert_to_dict (VBlockP vb_of_dict, ContextP ctx, DictType type, STRp(snip))
{
    BufferP dict = (type == DICT_ZF_SINGLETON) ? &ctx->stons : &ctx->dict; 

    static rom buf_name[3] = { "contexts->dict", "zctx->dict", "zctx->stons" };

    buf_alloc (vb_of_dict, dict, snip_len + 1, INITIAL_NUM_NODES * MIN_(10, snip_len), char, CTX_GROWTH, buf_name[type]);

    if (type == DICT_ZF) {
        ASSERT (snip_len <= Z_MAX_WORD_LEN, "ctx=%s has a snip with len=%u which is too big to fit in a dictionary. First 100 chars=\"%.*s\"",
                ctx->tag_name, snip_len, MIN_(100, snip_len), snip);

        ASSERT (dict->len  <= Z_MAX_DICT_LEN, "ctx=%s has a dictionary beyond maximum size of %"PRIu64". Example of a snip (first 100 chars)=\"%.*s\"",
                ctx->tag_name, (uint64_t)Z_MAX_DICT_LEN, MIN_(100, snip_len), snip);
    }
    
    CharIndex char_index = dict->len;
    char *dict_p = Bc (*dict, char_index);

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = 0; // dictionaries have a \0 separating snips, so that PIZ can generate word_list

    dict->len += (uint64_t)snip_len + 1;
    return char_index;
}

// ZIP only (PIZ doesn't have nodes) nodes index to node - possibly in ol_nodes, or in nodes
CtxNode *ctx_node_vb_do (ConstContextP vctx, WordIndex node_index, 
                         rom *snip_in_dict, uint32_t *snip_len,  // optional outs
                         FUNCLINE)
{
    ASSERT (vctx->dict_id.num, "this vctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index < vctx->nodes.len32 + vctx->ol_nodes.len32 && node_index >= 0, 
            "out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u. Caller: %s:%u",  
            vctx->tag_name, node_index, vctx->nodes.len32, vctx->ol_nodes.len32, func, code_line);

    bool is_ol = node_index < vctx->ol_nodes.len32; // is this entry from a previous vb (overlay buffer)

    CtxNode *node = is_ol ? B(CtxNode, vctx->ol_nodes, node_index)
                          : B(CtxNode, vctx->nodes, node_index - vctx->ol_nodes.len32);

    if (snip_in_dict) {
        ConstBufferP dict = is_ol ? &vctx->ol_dict : &vctx->dict;
        ASSERT0 (buf_is_alloc (dict), "dict not allocated");

        if (is_ol) {
            if (node->char_index + (uint64_t)node->snip_len >= dict->len) {
                ContextP zctx = ctx_get_zctx_from_vctx (vctx, false, true); // if is_ol, we know zctx exists

                ASSERT (node->char_index + (uint64_t)node->snip_len < dict->len, "Called from %s:%u: snip of %s out of range: node_index=%d is_ol=TRUE node->char_index=%"PRIu64" + node->snip_len=%u >= ol_dict->len=%"PRIu64". ol_dict->size=%"PRIu64". ol_nodes.len=%"PRIu64" zctx->nodes.len=%"PRIu64" zctx->dict.len=%"PRIu64" zctx->ston_nodes.len=%"PRIu64" zctx->stons.len=%"PRIu64,
                        func, code_line, vctx->tag_name, node_index, node->char_index, node->snip_len, dict->len, (uint64_t)dict->size, vctx->ol_nodes.len,
                        zctx->nodes.len, zctx->dict.len, zctx->ston_nodes.len, zctx->stons.len);
            }
        }
        else
            ASSERT (node->char_index + (uint64_t)node->snip_len < dict->len, "Called from %s:%u: snip of %s out of range: node_index=%d is_ol=FALSE node->char_index=%"PRIu64" + node->snip_len=%u >= dict->len=%"PRIu64". dict->size=%"PRIu64". ol_nodes.len=%"PRIu64,
                    func, code_line, vctx->tag_name, node_index, node->char_index, node->snip_len, dict->len, (uint64_t)dict->size, vctx->ol_nodes.len);

        *snip_in_dict = Bc (*dict, node->char_index);
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// ZIP only (PIZ doesn't have nodes) nodes index to node - possibly in ol_nodes, or in nodes
CtxNode *ctx_node_zf_do (ConstContextP zctx, int32_t node_index, 
                         rom *snip_in_dict, uint32_t *snip_len,  // optional outs
                         FUNCLINE)
{
    ASSERT (zctx->dict_id.num, "this zctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index > -2 - (int32_t)zctx->ston_nodes.len && node_index < (int32_t)zctx->nodes.len, 
            "out of range: dict=%s node_index=%d nodes.len=%u ston_nodes.len=%u. Caller: %s:%u",  
            zctx->tag_name, node_index, zctx->nodes.len32, zctx->ston_nodes.len32, func, code_line);

    bool is_singleton = node_index < 0; // is this entry from a previous vb (overlay buffer)

    CtxNode *node = is_singleton ? B(CtxNode, zctx->ston_nodes, -node_index - 2)
                                 : B(CtxNode, zctx->nodes, node_index);

    if (snip_in_dict) {
        ConstBufferP dict = is_singleton ? &zctx->stons : &zctx->dict;
        ASSERTISALLOCED(*dict);

        *snip_in_dict = Bc(*dict, node->char_index);
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// PIZ: search for a node matching this snip in a directory and return the node index. note that we do a linear
// search as PIZ doesn't have hash tables.
WordIndex ctx_search_for_word_index (ContextP ctx, STRp(snip))
{
    ASSERT (ctx->word_list.len, "word_list is empty for context %s", ctx->tag_name);

    for_buf2 (CtxWord, word, wi, ctx->word_list)
        if (word->len == snip_len && !memcmp (Bc(ctx->dict, word->index), snip, snip_len))
            return wi;

    return WORD_INDEX_NONE;
}

// PIZ
WordIndex ctx_decode_b250 (bytes *b, bool advance, B250Size b250_size, rom ctx_name)
{
    ASSERT (*b, "*b is NULL in ctx=%s", ctx_name);

    switch ((*b)[0]) {
        #define RETURN(res,n) if (advance) { *b += n; } return res
        
        case BASE250_MOST_FREQ0 : RETURN (0, 1);
        case BASE250_MOST_FREQ1 : RETURN (1, 1);
        case BASE250_MOST_FREQ2 : RETURN (2, 1);
        case BASE250_ONE_UP     : RETURN (WORD_INDEX_ONE_UP,  1);
        case BASE250_EMPTY_SF   : RETURN (WORD_INDEX_EMPTY,   1);
        case BASE250_MISSING_SF : RETURN (WORD_INDEX_MISSING, 1);
        default /* 0 ... 249 */ : {
            WordIndex value;
            switch (b250_size) {
                case B250_BYTES_1: 
                    value = (*b)[0]; 
                    RETURN (value, 1);
                case B250_BYTES_2:
                    value = ((uint32_t)(*b)[0] << 8) | (uint32_t)(*b)[1]; // careful not to use BGEN as string might not be aligned to word boundary
                    RETURN (value, 2);
                case B250_BYTES_3:
                    value = ((uint32_t)(*b)[0] << 16) | ((uint32_t)(*b)[1] << 8) | (uint32_t)(*b)[2]; 
                    RETURN (value, 3);
                case B250_BYTES_4:
                    value = ((uint32_t)(*b)[0] << 24) | ((uint32_t)(*b)[1] << 16) | ((uint32_t)(*b)[2] << 8) | (uint32_t)(*b)[3]; 
                    RETURN (value, 4);
                default:
                    ABORT ("Invalid b250_size=%u", b250_size);
            }
        }
        #undef RETURN
    }
}

// PIZ and reading Pair1 in ZIP (uses word_list): returns word index, and advances the iterator
WordIndex ctx_get_next_snip (VBlockP vb, ContextP ctx, bool is_pair, pSTRp (snip)/*optional out*/)
{
    ASSERT (ctx, "%s: ctx is NULL", VB_NAME);

    bool zip_pair = (is_pair && command==ZIP);

    ConstBufferP b250      = is_pair  ? &ctx->b250R1                 : &ctx->b250;
    B250Size b250_size     = is_pair  ? ctx->pair_b250_size          : ctx->b250_size;
    SnipIterator *iterator = is_pair  ? &ctx->pair_b250_iter         : &ctx->iterator;
    bool all_the_same      = is_pair  ? ctx->pair_flags.all_the_same : ctx->flags.all_the_same;
    Buffer *list           = zip_pair ? &ctx->ol_nodes               : &ctx->word_list;
    
    // if the entire b250 in a VB consisted of word_index=0, we don't output the b250 to the file, and just 
    // consider it to always emit 0
    if (!b250->len) {
        if (snip) {
            if (!zip_pair) {
                ASSERT (ctx->word_list.len32, "%s.word_list.len32=0", ctx->tag_name);
                
                CtxWord *dict_word = B1ST (CtxWord, ctx->word_list);

                ASSERT (dict_word->index + dict_word->len < ctx->dict.len, // < and not <= because seperator \0 is not included in word len
                        "expecting: %s.dict_word->index=%"PRIu64" + len=%"PRIu64" < %s->dict.len=%"PRIu64,
                        ctx->tag_name, (uint64_t)dict_word->index, (uint64_t)dict_word->len, ctx->tag_name, ctx->dict.len);

                *snip = Bc (ctx->dict, dict_word->index);
                *snip_len = dict_word->len;
            }
            else {
                ASSERT (ctx->ol_nodes.len32, "%s.ol_nodes.len32=0", ctx->tag_name);

                CtxNode *dict_word = B1ST (CtxNode, ctx->ol_nodes);

                ASSERT (dict_word->char_index + dict_word->snip_len < ctx->dict.len, // likewise, < and not <=
                        "expecting: %s.dict_word->char_index=%"PRIu64" + snip_len=%u < %s->dict.len=%"PRIu64,
                        ctx->tag_name, dict_word->char_index, dict_word->snip_len, ctx->tag_name, ctx->dict.len);

                *snip = Bc (ctx->dict, dict_word->char_index);
                *snip_len = dict_word->snip_len;
            }
        }
        return 0;
    }
    
    if (!iterator->next_b250)  // initialize
        *iterator = (SnipIterator){ .next_b250       = B1ST8 (*b250), 
                                    .prev_word_index = WORD_INDEX_NONE };

    WordIndex word_index = ctx_decode_b250 (&iterator->next_b250, !all_the_same, b250_size, ctx->tag_name);  // if this line has no non-GT subfields, it will not have a ctx 

    // we check after (no risk of segfault because of buffer overflow protector) - since b250 word is variable length
    if (IS_PIZ)
        ASSPIZ (iterator->next_b250 <= BAFT8 (*b250), 
                "while reconstructing (vb->lines.len=%u): iterator for %s(%u) %sreached end of b250. %s.len=%u, preprocessing=%s.%s", 
                vb->lines.len32, ctx->tag_name, ctx->did_i, is_pair ? "(PAIR) ": "", is_pair ? "pair" : "b250", b250->len32, TF(vb->preprocessing),
                b250->len32 ? "" : " Check Skip function (since b250.len=0).");
    else
        ASSERT (iterator->next_b250 <= BAFT8 (*b250), 
                "%s: while reconstructing (last line of vb is %d): iterator for %s(%u) %sreached end of b250. %s.len=%u", 
                LN_NAME, vb->lines.len32 - 1, ctx->tag_name, ctx->did_i, is_pair ? "(PAIR) ": "", is_pair ? "pair" : "b250", b250->len32);

    // case: a Container item is missing (eg a subfield in a Sample, a FORMAT or Samples items in a file)
    if (word_index == WORD_INDEX_MISSING) {
        if (snip) {
            *snip = NULL; // ignore this dict_id - don't even output a separator
            *snip_len = 0;
        }
    }

    // case: a subfield snip is empty, eg "AB::CD" (VCF GT Data) or "OA:Z:chr13,52863337,-,56S25M70S,0,;" (SAM OA optional field)
    else if (word_index == WORD_INDEX_EMPTY) { 
        if (snip) {
            *snip = ""; // pointer to static empty string
            *snip_len = 0;
        }
    }

    else {
        if (word_index == WORD_INDEX_ONE_UP) 
            word_index = iterator->prev_word_index + 1;

        ASSERT (word_index < list->len32, 
                "%s: word_index=%u but %s(%u).%s.len=%u (is_pair=%s, b250.len=%u, iterator(after)=%u)",
                LN_NAME, word_index, ctx->tag_name, ctx->did_i, zip_pair ? "ol_nodes" : "word_list", list->len32, 
                TF(is_pair), b250->len32, BNUM(*b250, iterator->next_b250));

        if (!zip_pair) {
            CtxWord *dict_word = B(CtxWord, ctx->word_list, word_index);
            if (snip) {
                *snip = Bc (ctx->dict, dict_word->index);
                *snip_len = dict_word->len;
            }
        }

        // case: iterating on pair-1 data while compressing pair-2 (FASTQ)
        else {
            CtxNode *dict_word = B(CtxNode, ctx->ol_nodes, word_index);
            if (snip) {
                *snip = Bc (ctx->ol_dict, dict_word->char_index);
                *snip_len = dict_word->snip_len;
            }
        }
    }
    
    iterator->prev_word_index = word_index;    

    return word_index;
}

// similar to ctx_get_next_snip, except iterator is left intact 
WordIndex ctx_peek_next_snip (VBlockP vb, ContextP ctx, pSTRp (snip)/*optional out*/)  
{
    SnipIterator save_iterator = ctx->iterator;
    WordIndex wi = ctx_get_next_snip (vb, ctx, false, STRa(snip));
    ctx->iterator = save_iterator;
    return wi;
}

// Process a snip durig merge - add it to zctx->dict or zctx->stons return its node index. Also used to pre-populate zctx.
static WordIndex ctx_commit_node (VBlockP vb, ContextP zctx, ContextP vctx, STRp(snip), int64_t count, CtxNode **node, bool *is_new)  // out
{
    // if this turns out to be a singelton - i.e. a new snip globally - where it goes depends on whether its a singleton in the VB 
    // and whether we are allowed to move singletons to local. if its a singleton:
    // 1. we keep it in ston_nodes and the index in the node is negative to indicate that
    // 2. we insert it to ston_nodes instead of dict - i.e. it doesn't get written the dict section
    // 3. we move it to the local section of this vb
    // 4. we set the word_index of its nodes to be the word_index of the SNIP_LOOKUP snip
    bool is_singleton_in_vb = (count == 1 && ctx_can_have_singletons (vctx)); // is singleton in this VB

    // attempt to get the node from the hash table 
    WordIndex node_index = hash_global_get_entry (zctx, STRa(snip), 
                                                  is_singleton_in_vb ? HASH_NEW_OK_SINGLETON_IN_VB : HASH_NEW_OK_NOT_SINGLETON, node);

    // case: existing non-singleton node. Possibly it was a singleton node before, that thanks to us hash_global_get_entry
    // converted to non-singleton - i.e. node_index is always >= 0.
    if (*node) { 
        *is_new = false;
        return node_index; // >= 0
    }
    
    // NEW SNIP globally - this snip was just added to the hash table - either as a regular or singleton node
    bool is_singleton_in_global = (node_index < 0);
    BufferP nodes = is_singleton_in_global ? &zctx->ston_nodes : &zctx->nodes;
    ASSERT (nodes->len <= MAX_WORDS_IN_CTX, 
            "too many words in ctx %s, max allowed number of words is is %u", zctx->tag_name, MAX_WORDS_IN_CTX);

    // set either singleton node or regular node with this snip
    buf_alloc (evb, nodes, 1, INITIAL_NUM_NODES, CtxNode, CTX_GROWTH, is_singleton_in_global ? "zctx->ston_nodes" : "zctx->nodes");
    *node = BLST (CtxNode, *nodes); // note: (*nodes).len was already incremented in hash_global_get_entry
    **node = (CtxNode){
        .snip_len   = snip_len,
        .char_index = ctx_insert_to_dict (evb, zctx, (is_singleton_in_global ? DICT_ZF_SINGLETON : DICT_ZF), STRa(snip))
    };

    // case: vb singleton turns out to be a global singleton - we add it to local and return the SNIP_LOOKUP node
    // instead of the singleton node (which is guaranteed to be non-singleton, and hence >= 0)
    // note: local is dedicated to singletons and contains nothing else, since inst.no_stons is not set
    if (node_index < 0) {
        seg_add_to_local_fixed_do (vb, vctx, STRa(snip), true, LOOKUP_NONE, true, 0); 

        if (HAS_DEBUG_SEG(vctx)) {
            char printable_snip[snip_len+20];
            iprintf ("commit_node: %s: SINGLETON replaced by SNIP_LOCAL: snip=%s snip_len=%u\n", 
                     vctx ? vctx->tag_name : "", str_print_snip (STRa(snip), printable_snip), snip_len);
        }

        static char lookup = SNIP_LOOKUP;
        return ctx_commit_node (vb, zctx, vctx, &lookup, 1, -1 /* not singleton */, node, is_new);
    }

    // case: not singleton - we return this (new) node
    else {
        (*node)->word_index = node_index;

        buf_alloc_zero (evb, &zctx->counts, 1, INITIAL_NUM_NODES, uint64_t, CTX_GROWTH, "zctx->counts");
        zctx->counts.len++; // actually assigned in ctx_merge_in_one_vctx 

        if (is_new) *is_new = true;

        return node_index; // >= 0
    }
}

// Seg: inserts snip into the hash, nodes and dictionary, if it was not already there, and returns its node index.
// Does NOT add the word index to b250.
WordIndex ctx_create_node_do (VBlockP vb, ContextP vctx, STRp(snip), bool *is_new/*optional out*/)
{
    ASSERTNOTNULL (vctx);
    ASSERT (vctx->dict_id.num, "vctx has no dict_id (did_i=%u)", (unsigned)(vctx - vb->contexts));
    ASSERT (snip || !snip_len, "vctx=%s: snip=NULL but snip_len=%u", vctx->tag_name, snip_len);
    
    WordIndex node_index;
    #define RETURN(ni) { node_index=(ni); goto done; }

    if (!snip_len) {
        if (is_new) *is_new = false;
        RETURN ((!snip || (vb->data_type == DT_VCF && dict_id_is_vcf_format_sf (vctx->dict_id) && *snip != ':')) 
                 ? WORD_INDEX_MISSING : WORD_INDEX_EMPTY);
    }

    // short-circuit if identical to previous snip
    if (vctx->last_snip && str_issame (snip, vctx->last_snip)) {
        (*B32(vctx->counts, vctx->last_snip_ni))++; 
        if (is_new) *is_new = false;
        RETURN (vctx->last_snip_ni);
    }
    
    WordIndex node_index_if_new = vctx->ol_nodes.len + vctx->nodes.len;
    
#ifdef DEBUG // time consuming and only needed during development
    unsigned actual_len = strnlen (snip, snip_len);
    ASSERT (actual_len == snip_len, "vb=%u ctx=%s: snip_len=%u but unexpectedly has an 0 at index %u: \"%.*s\"", 
            vb->vblock_i, vctx->tag_name, snip_len, actual_len, STRf(snip));
#endif

    ASSERT (node_index_if_new <= MAX_WORD_INDEX, 
            "ctx of %s is full (max allowed words=%u): ol_nodes.len=%u nodes.len=%u",
            vctx->tag_name, MAX_WORDS_IN_CTX, vctx->ol_nodes.len32, vctx->nodes.len32);

    // get the node from the hash table if it already exists, or add this snip to the hash table if not
    CtxNode *node;
    WordIndex existing_node_index = hash_get_entry_for_seg (vb, vctx, STRa(snip), node_index_if_new, &node);
    if (existing_node_index != WORD_INDEX_NONE) {
        (*B32(vctx->counts, existing_node_index))++; // note: counts.len = nodes.len + ol_nodes.len
        if (is_new) *is_new = false;
        
        // populate last_snip 
        bool is_ol = existing_node_index < vctx->ol_nodes.len; // is this entry from a previous vb (overlay buffer)
        ConstBufferP dict = is_ol ? &vctx->ol_dict : &vctx->dict;
        vctx->last_snip     = Bc (*dict, node->char_index);
        vctx->last_snip_len = snip_len;
        vctx->last_snip_ni  = existing_node_index;

        RETURN (existing_node_index); // snip found - we're done
    }
    
    // this snip isn't in the hash table - its a new snip
    ASSERT (vctx->nodes.len < MAX_WORD_INDEX, "too many words in dictionary %s (MAX_WORD_INDEX=%u)", vctx->tag_name, MAX_WORD_INDEX);

    buf_alloc (vb, &vctx->nodes,  1, INITIAL_NUM_NODES, CtxNode,  CTX_GROWTH, "contexts->nodes");
    buf_alloc (vb, &vctx->counts, 1, INITIAL_NUM_NODES, uint32_t, CTX_GROWTH, "contexts->counts");

    BNXT (CtxNode, vctx->nodes) = (CtxNode){
        .snip_len   = snip_len,
        .char_index = ctx_insert_to_dict (vb, vctx, DICT_VB, STRa(snip)),
        .node_index = node_index_if_new
    };

    // populate last_snip / last_snip_len if requested
    vctx->last_snip     = Bc (vctx->dict, BLST (CtxNode, vctx->nodes)->char_index);
    vctx->last_snip_len = snip_len;
    vctx->last_snip_ni  = node_index_if_new;

    BNXT32 (vctx->counts) = 1;

    ASSERT (vctx->counts.len == vctx->nodes.len + vctx->ol_nodes.len, "Expecting vctx->counts.len=%"PRId64" == vctx->nodes.len=%"PRId64" + vctx->ol_nodes.len=%"PRId64,
            vctx->counts.len, vctx->nodes.len, vctx->ol_nodes.len);

    if (is_new) *is_new = true;
    RETURN (node_index_if_new);

done:
    if (HAS_DEBUG_SEG(vctx)) { 
        char printable_snip[snip_len+20];
        if (snip) iprintf ("create_node: %s %s(%u): snip=%s snip_len=%u node_index=%d\n",  LN_NAME, vctx->tag_name, vctx->did_i, str_print_snip (STRa(snip), printable_snip), snip_len, node_index);
        else      iprintf ("create_node: %s %s(%u): snip=NULL snip_len=0 node_index=%d\n", LN_NAME, vctx->tag_name, vctx->did_i, node_index);
    }

    return node_index;
    #undef RETURN
}

// Seg only: create a node without adding to counts
WordIndex ctx_create_node_is_new (VBlockP vb, Did did_i, STRp (snip), bool *is_new)
{
    WordIndex node_index = ctx_create_node_do (vb, CTX(did_i), STRa(snip), is_new); 
    ctx_decrement_count (vb, CTX(did_i), node_index);

    return node_index;
}

uint32_t ctx_get_count (VBlockP vb, ContextP ctx, WordIndex node_index)
{
    ASSERT (node_index >= 0 && node_index < (WordIndex)ctx->counts.len, "node_index=%d out of range counts[%s].len=%"PRIu64, node_index, ctx->tag_name, ctx->counts.len);

    return *B32(ctx->counts, node_index);
}

// Seg only: if after ctx_create_node_do we don't add the snip to b250, we need to reduce its count
void ctx_decrement_count (VBlockP vb, ContextP ctx, WordIndex node_index)
{
    ASSERT (node_index < (WordIndex)ctx->counts.len, "node_index=%d out of range counts[%s].len=%"PRIu64, node_index, ctx->tag_name, ctx->counts.len);

    uint32_t *count_p = B32(ctx->counts, node_index);

    if (node_index < 0) return; // WORD_INDEX_EMPTY or WORD_INDEX_MISSING
     
    ASSERT (*count_p >= 1, "count[%s]=%u too low to be decremented", ctx->tag_name, *count_p);
    (*count_p)--;

    // edge case: if we removed the last b250 of this node_index, we now have an unused node and an unused word dict. 
    // this may in some cases cause the condition "ctx->nodes.len != ctx->b250.len" in zip_handle_unique_words_ctxs to
    // incorrectly fail, causing moving of an incorrect dict to local. to prevent this, we don't allow singletons in this case.
    if (! *count_p && node_index >= ctx->ol_nodes.len) ctx->no_stons = true;
}

// Seg only: if we add a b250 without evaluating (if node_index is known)
void ctx_increment_count (VBlockP vb, ContextP ctx, WordIndex node_index)
{
    ASSERT (node_index < ctx->counts.len, "%s: node_index=%d out of range counts[%s].len=%"PRIu64, LN_NAME, node_index, ctx->tag_name, ctx->counts.len);

    (*B32(ctx->counts, node_index))++;
}

// Seg only: mark this node as one that should NOT be removed by ctx_shorten_unused_dict_words if it is unused
void ctx_protect_from_removal (VBlockP vb, ContextP ctx, WordIndex node_index)
{
    ASSERT (node_index < ctx->counts.len, "node_index=%d out of range counts[%s].len=%"PRIu64, node_index, ctx->tag_name, ctx->counts.len);

    *B32(ctx->counts, node_index) |= COUNT_PROTECTED_FROM_REMOVAL32;
}

void ctx_append_b250 (VBlockP vb, ContextP vctx, WordIndex node_index)
{
    // case: context is currently all_the_same - (ie len>=1, but only one actual entry)...
    if (vctx->flags.all_the_same) {
        WordIndex first_node_index; 

        // ... and the new entry is more of the same
        if (node_index == (first_node_index = *B1ST (WordIndex,vctx->b250)))
            vctx->b250.len++; // increment len, but don't actually grow the array

        // ... and the new entry is different - add (len-1) copies of first_node, and the new one
        else {
            buf_alloc (vb, &vctx->b250, 1, AT_LEAST(vctx->did_i), WordIndex, CTX_GROWTH, CTX_TAG_B250); // add 1 more, meaning a total of len+1
        
            for (unsigned i=1; i < vctx->b250.len32; i++) *B(WordIndex, vctx->b250, i) = first_node_index;
            BNXT32 (vctx->b250) = node_index;

            vctx->flags.all_the_same = false; // no longer all_the_same 
        }
    }

    // b250 is not all_the_same: either it is empty, or it contains multiple different node_index values
    else {
        // case: this is the first entry - its all_the_same in a trivial way (still, we mark it so zip_generate_b250_section can potentially eliminate it)
        if (!vctx->b250.len) vctx->flags.all_the_same = true;

        buf_alloc (vb, &vctx->b250, 1, AT_LEAST(vctx->did_i), WordIndex, CTX_GROWTH, CTX_TAG_B250);
        BNXT32 (vctx->b250) = node_index;
    }
}

// ZIP compute thread: overlay and/or copy the current state of the global contexts to the vb, ahead of segging this vb.
void ctx_clone (VBlockP vb)
{
    Did z_num_contexts = __atomic_load_n (&z_file->num_contexts, __ATOMIC_ACQUIRE);

    START_TIMER; // including mutex wait time

    // note: because each dictionary has its own mutex, it is possible that we will see only a partial set
    // of dictionaries (eg some but not all of the fields) when we are arrive here while another thread is mid-way 
    // through merging and adding a bunch of dictionaries.
    // however z_num_contexts will always correctly state the number of dictionaries that are available.
    uint32_t num_cloned = 0;
    bool cloned[z_num_contexts];
    memset (cloned, 0, z_num_contexts);

    while (num_cloned < z_num_contexts) {

        bool achieved_something = false;
        for (Did did_i=0; did_i < z_num_contexts; did_i++) {
            if (cloned[did_i]) continue;

            ContextP vctx = CTX(did_i);
            ContextP zctx = ZCTX(did_i);
            ContextP dict_ctx = ZCTX(zctx->dict_did_i); // dict_did_i is equial to either did_i or, if its an alias, to the destination did_i

            // case: this context doesn't really exist (happens when incrementing num_contexts when adding RNAME in ctx_populate_zf_ctx_from_contigs)
            if (!ZMUTEX(dict_ctx).initialized) goto did_i_cloned;

            if (!mutex_trylock (ZMUTEX(dict_ctx)))
                continue;

            if (dict_ctx->dict.len) {  // there's something already in this dict
                // overlay the global dict and nodes - these will be treated as read-only by the VB
                buf_overlay (vb, &vctx->ol_dict,  &dict_ctx->dict,  "contexts->ol_dict");   
                buf_overlay (vb, &vctx->ol_nodes, &dict_ctx->nodes, "contexts->ol_nodes");   

                // overlay the hash table, that may still change by future vb's merging... this vb will only use
                // entries that are up to this merge_num
                buf_overlay (vb, &vctx->global_hash, &dict_ctx->global_hash, "contexts->global_hash");
                buf_overlay (vb, &vctx->global_ents, &dict_ctx->global_ents, "contexts->global_ents");

                vctx->merge_num = dict_ctx->merge_num;
                vctx->global_hash_prime = dict_ctx->global_hash_prime; // can never change
                vctx->num_new_entries_prev_merged_vb = dict_ctx->num_new_entries_prev_merged_vb;
            }

            vctx->did_i       = did_i;
            vctx->dict_id     = zctx->dict_id;
            vctx->st_did_i    = zctx->st_did_i;
            vctx->dict_did_i  = zctx->dict_did_i;
            vctx->luft_trans  = zctx->luft_trans;
            vctx->last_line_i = LAST_LINE_I_INIT;
            vctx->tag_i       = -1;

            // note: lcodec and bcodec are inherited in merge (see comment in codec_assign_best_codec)

            memcpy ((char*)vctx->tag_name, zctx->tag_name, sizeof (vctx->tag_name));
            
            // note: ZCTX(CHROM)->chrom2ref_map is protected by ZCTX[CHROM]->mutex
            if (chrom_2ref_seg_is_needed(did_i) && ZCTX(CHROM)->chrom2ref_map.len) 
                buf_overlay (vb, &vb->ol_chrom2ref_map, &ZCTX(CHROM)->chrom2ref_map, "ol_chrom2ref_map");

            mutex_unlock (ZMUTEX(dict_ctx));

            // stuff that doesn't require mutex
        
            if (vctx->ol_nodes.len) {
                buf_alloc_zero (vb, &vctx->counts, 0, vctx->ol_nodes.len, uint32_t, CTX_GROWTH, "contexts->counts");
                vctx->counts.len = vctx->ol_nodes.len;
            }

            set_d2d_map (vb->d2d_map, vctx->dict_id, did_i);

            ctx_init_iterator (vctx);

        did_i_cloned:
            cloned[did_i] = achieved_something = true;
            num_cloned++;
        }

        if (!achieved_something) sched_yield(); // all the contexts we still need to clone are locked - allow context switch.
    }

    vb->num_contexts = z_num_contexts;
       
    COPY_TIMER (ctx_clone);
}

static void ctx_initialize_ctx (ContextP ctx, Did did_i, DictId dict_id, DictIdtoDidMap d2d_map, STRp(tag_name))
{
    ctx->did_i            = did_i;
    ctx->dict_did_i       = did_i;
    ctx->st_did_i         = DID_NONE; // this is other_did_i in PIZ
    ctx->dict_id          = dict_id;
    ctx->last_line_i      = LAST_LINE_I_INIT;
    ctx->tag_i            = -1;
    ctx->dict.can_be_big  = true; // don't warn if dict buffers grow really big
    ctx->pair_assist_type = SEC_NONE;

    if (tag_name_len) {
        ASSINP (tag_name_len <= MAX_TAG_LEN-1, "Tag name \"%.*s\" is of length=%u beyond the maximum tag length supported by Genozip=%u",
                tag_name_len, tag_name, tag_name_len, MAX_TAG_LEN-1);

        memcpy ((char*)ctx->tag_name, tag_name, tag_name_len);
    }
    else
        strcpy ((char*)ctx->tag_name, dis_dict_id (dict_id).s);

    ctx_init_iterator (ctx);
    
    set_d2d_map (d2d_map, dict_id, did_i);

    // add a user-requested SEC_COUNTS section
    if (IS_ZIP) {
        if (flag.show_one_counts.num == dict_id_typeless (ctx->dict_id).num) 
            ctx->counts_section = ctx->no_stons = true;

        if (z_file && is_p_in_range (ctx, z_file->contexts, sizeof(ContextArray))) { // z_file context
            mutex_initialize (ZMUTEX(ctx));
            buf_set_shared (&ctx->dict);
            buf_set_shared (&ctx->nodes);
        }
    } 
}

WordIndex ctx_populate_zf_ctx (Did dst_did_i, STRp (contig_name), WordIndex ref_index)
{
    decl_zctx(dst_did_i);
    
    // zctx is already initialized, but this is the first call to populate it with a contig
    if (!buf_is_alloc (&zctx->global_hash)) { 
        zctx->no_stons = true;
        hash_alloc_global (zctx, 0);
    }

    CtxNode *zf_node;
    WordIndex zf_node_index = ctx_commit_node (NULL, zctx, NULL, STRa(contig_name), 0, &zf_node, NULL);
    zf_node->word_index = zf_node_index;

    if (dst_did_i == CHROM/*primary*/ && (flag.reference & REF_ZIP_LOADED)) {
        buf_alloc_255 (evb, &zctx->chrom2ref_map, 0, MAX_(INITIAL_NUM_NODES, zf_node_index+1), WordIndex, CTX_GROWTH, "ZCTX(CHROM)->chrom2ref_map");
        *B(WordIndex, zctx->chrom2ref_map, zf_node_index) = ref_index;
        zctx->chrom2ref_map.len32 = MAX_(zctx->chrom2ref_map.len32, zf_node_index+1);
    }

    if (zctx->dict.len > EXCESSIVE_DICT_SIZE)
        zctx->dict_len_excessive = true; // suppress warning - as this would not be an unexpectedly large dict size due to bad segging
        
    return zf_node_index;
}

// ZIP main thread: 
// 1. when starting to zip a new file, with pre-loaded external reference, we integrate the reference's FASTA CONTIG
//    dictionary as the chrom dictionary of the new file
// 2. in SAM, DENOVO: after creating contigs from SQ records, we copy them to the RNAME dictionary
// 3. When loading a chain file - copy tName word_list/dict read from the chain file to a context
void ctx_populate_zf_ctx_from_contigs (Reference ref, Did dst_did_i, ConstContigPkgP ctgs)
{
    if (!ctgs->contigs.len) return; // nothing to do

    decl_zctx(dst_did_i);
    
    // initialize if empty ctx. note: it might not be empty if contigs have been added already from the txt header
    if (!buf_is_alloc (&zctx->global_hash)) { 
        zctx->no_stons = true;
        hash_alloc_global (zctx, ctgs->contigs.len32);
    }

    if (flag.reference & REF_ZIP_LOADED)
        buf_alloc_255 (evb, &ZCTX(CHROM)->chrom2ref_map, ctgs->contigs.len, INITIAL_NUM_NODES, WordIndex, 1, "ZCTX(CHROM)->chrom2ref_map");

    bool is_primary = (dst_did_i == DTFZ (prim_chrom)); // note: primary is not CHROM in the case of DT_CHAIN

    for_buf2 (Contig, contig, i, ctgs->contigs) {
        CtxNode *zf_node;
        bool is_new;

        STR(contig_name);
        contig_name = contigs_get_name (ctgs, i, &contig_name_len);

        WordIndex zf_node_index = ctx_commit_node (NULL, zctx, NULL, STRa(contig_name), 0, &zf_node, &is_new);
        zf_node->word_index = zf_node_index;

        if (is_new && is_primary && (flag.reference & REF_ZIP_LOADED)) {
            if (contig->ref_index >= 0) buf_add_int (evb, ZCTX(CHROM)->chrom2ref_map, contig->ref_index); // header contigs also know the ref_index
            else                        buf_add_int (evb, ZCTX(CHROM)->chrom2ref_map, ref_contigs_get_matching (ref, 0, STRa(contig_name), NULL, NULL, false, NULL, NULL));
        }
    }
}

// Find the z_file context that corresponds to dict_id, or return NULL if there is none.
// ZIP note: It could be possibly a different did_i than in the vb - in case this dict_id is new to this vb, but another 
// vb already inserted it to z_file
ContextP ctx_get_existing_zctx (DictId dict_id)
{
    Did did_i = get_matching_did_i_from_map (z_file->contexts, z_file->d2d_map, dict_id);
    
    if (did_i != DID_NONE)
        return &z_file->contexts[did_i];
    
    else {
        Did z_num_contexts = IS_ZIP ? __atomic_load_n (&z_file->num_contexts, __ATOMIC_ACQUIRE)
                                    : z_file->num_contexts; // in PIZ, num_contexts is immutable so no need for atomic

        for (ContextP zctx=&z_file->contexts[0]; zctx < &z_file->contexts[z_num_contexts]; zctx++)
            if (dict_id.num == zctx->dict_id.num) 
                return zctx;
    }
    
    return NULL;
}

static ContextP ctx_add_new_zf_ctx_do (STRp (tag_name), DictId dict_id, TranslatorId luft_trans, 
                                       Did st_did_i, bool is_stats_parent)
{
    ASSERT (z_file->num_contexts+1 < MAX_DICTS, // load num_contexts - this time with mutex protection - it could have changed
            "z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    decl_zctx(z_file->num_contexts);
    zctx->did_i           = z_file->num_contexts; 
    zctx->dict_id         = dict_id;
    zctx->luft_trans      = luft_trans;
    zctx->st_did_i        = st_did_i;
    zctx->dict_did_i      = zctx->did_i; // this is a new context -> it is not a predefined context -> it is not an alias
    zctx->is_stats_parent = is_stats_parent;

    memcpy ((char*)zctx->tag_name, tag_name, tag_name_len);
    // note: lcodec is NOT copied here, see comment in codec_assign_best_codec

    mutex_initialize (ZMUTEX(zctx));

    // only when the new entry is finalized, do we increment num_contexts, atmoically, this is because
    // other threads might access it without a mutex when searching for a dict_id
    __atomic_fetch_add (&z_file->num_contexts, 1, __ATOMIC_ACQ_REL); 

    // only after updating num_contexts, we add it to the d2d_map. 
    set_d2d_map (z_file->d2d_map, zctx->dict_id, zctx->did_i);

    buf_set_shared (&zctx->dict);
    buf_set_shared (&zctx->nodes);

    return zctx;
}

// ZIP only: called by main thread when inspecting a txtheader for assigning liftover translators
ContextP ctx_add_new_zf_ctx_from_txtheader (STRp(tag_name), DictId dict_id, TranslatorId luft_translator)
{
    ASSINP (tag_name_len <= MAX_TAG_LEN-1, "Tag name \"%.*s\" is of length=%u beyond the maximum tag length supported by Genozip=%u",
            tag_name_len, tag_name, tag_name_len, MAX_TAG_LEN-1);

    mutex_lock (z_file->dicts_mutex); // note: mutex needed bc VBs of a previous component might still be merging

    ContextP zctx = ctx_get_existing_zctx (dict_id) 
                  ? NULL // not new - possibly a tag is duplicate in the header, or two tag names map to the same dict_id
                  : ctx_add_new_zf_ctx_do (STRa(tag_name), dict_id, luft_translator, DID_NONE, false);

    mutex_unlock (z_file->dicts_mutex);

    return zctx; 
}

// ZIP only: called by merging VBs to add a new dict to z_file - copying some stuff from vctx
// note: this is never called for predefined contexts
static ContextP ctx_add_new_zf_ctx (ConstContextP vctx)
{
    mutex_lock (z_file->dicts_mutex);

    // check if another thread raced and created this dict before us
    ContextP zctx = ctx_get_zctx_from_vctx (vctx, false, false);
    if (!zctx) 
        zctx = ctx_add_new_zf_ctx_do (vctx->tag_name, sizeof(zctx->tag_name), vctx->dict_id, 
                                      vctx->luft_trans, vctx->st_did_i, vctx->is_stats_parent);

    mutex_unlock (z_file->dicts_mutex);
    return zctx;
}

ContextP ctx_get_zctx_from_vctx (ConstContextP vctx, 
                                 bool create_if_missing, // if false, returns NULL if context doesn't exist
                                 bool follow_alias)      // if true, returns the ALIAS_DICT context if there is one
{
    if (vctx->did_i < DTFZ(num_fields)) 
        return follow_alias ? ZCTX(ZCTX(vctx->did_i)->dict_did_i) : ZCTX(vctx->did_i); // aliases only exist for predefined contexts
    
    // check dict_id->did_i map, and if not found search linearly
    Did did_i = ctx_get_existing_did_i_do (vctx->dict_id, z_file->contexts, z_file->d2d_map, NULL, z_file->num_contexts);
    
    return (did_i != DID_NONE) ? ZCTX(did_i) 
         : create_if_missing   ? ctx_add_new_zf_ctx (vctx) 
         :                       NULL;
}

// update zctx with codec as it is assigned - don't wait for merge, to increase the chance that subsequent
// VBs can get this codec and don't need to test for themselves.
void ctx_commit_codec_to_zf_ctx (VBlockP vb, ContextP vctx, bool is_lcodec, bool is_lcodec_inherited)
{
    // case: context might not exist yet in z_file, because no VB with it has merged yet - in which case we create it.
    ContextP zctx = ctx_get_zctx_from_vctx (vctx, true, false);

    mutex_lock (ZMUTEX(zctx));

    if (is_lcodec) {
        if (zctx->lcodec == vctx->lcodec) {
            if (zctx->lcodec_count < 255) zctx->lcodec_count++; // counts number of VBs in a row that set this codec
        }
        else if (is_lcodec_inherited) {
            zctx->lcodec_count = 0;
            zctx->lcodec = vctx->lcodec; 
        }
        else 
            zctx->lcodec_non_inherited = vctx->lcodec; 
    }
    else {
        if (zctx->bcodec == vctx->bcodec) {
            if (zctx->bcodec_count < 255) zctx->bcodec_count++; 
        }
        else {
            zctx->bcodec_count = 0;
            zctx->bcodec = vctx->bcodec; 
        }
    }           

    mutex_unlock (ZMUTEX(zctx));
}

void ctx_reset_codec_commits (void)
{
    for_zctx {
        if (zctx->lcodec && !zctx->lcodec_non_inherited) 
            zctx->lcodec_non_inherited = zctx->lcodec; // for stats
        
        zctx->lcodec = zctx->bcodec = CODEC_UNKNOWN; 
        zctx->lcodec_count = zctx->bcodec_count = 0;
    }
}

static inline void ctx_drop_all_the_same (VBlockP vb, ContextP zctx, ContextP vctx)
{
    rom reason=0;
    #define NO_DROP(s) { reason=(s); goto no_drop; }

    // note: pair_b250/local will causes zfile_compress_b250/local_data to set flags.paired, which will be different than vb_i=1 flags
    if (!vctx->flags.all_the_same) return;
    
    if (vctx->no_drop_b250) NO_DROP("no_drop_b250 is set");

    if (is_fastq_pair_2 (vb) && fastq_zip_use_pair_identical (vctx->dict_id)) {
        if (vctx->b250R1.len) 
            NO_DROP("pair-identical context: dropping R2.b250 would cause loading R1.b250 instead"); // even though not identical
        
        if (vctx->localR1.len && !vctx->local.len) 
            NO_DROP("pair-identical context: dropping b250 will cause recon of R1.local and not dict[0]"); // since PIZ will load R1.local as there is no R2.local
    }

    ASSERTISALLOCED (vctx->b250);
    WordIndex node_index = *B1ST (WordIndex, vctx->b250); // the only b250 in this context, as it is all_the_same
    
    // - if we have ol_node - we can only drop the existing word_index=0
    // - if we don't have ol_node - we can only drop new node_index=0 (possibly we have ctx_create_node'd several)
    // - a special WORD_INDEX_* - not droppable
    if (node_index != 0) NO_DROP ("node_index is not 0"); 

    bool is_new_word = node_index >= vctx->ol_nodes.len32; // its a new word if its not in the overlayed nodes (from previous VBs)
    
    if (!is_new_word && B(CtxNode, vctx->ol_nodes, node_index)->word_index > 0) NO_DROP ("existing word_index is not 0"); // old word - but not word_index=0

    if (is_new_word && 
        zctx->dict.len && strcmp (zctx->dict.data, vctx->dict.data)) // if it is_new_word - b250 can still be droppable if this new word (the only word in vctx->dict) is the same as word_index_i=0 in zctx->dict, created independently by two VBs in parallel 
        NO_DROP ("new word_index is not 0"); // new word - but not the first in the dictionary so not word_index=0

    STR(snip);
    ctx_get_vb_snip_ex (vctx, node_index, pSTRa(snip)); // get the snip of the only b250 we have

    if (*snip == SNIP_SELF_DELTA) NO_DROP ("snip is SNIP_SELF_DELTA");

    bool is_simple_lookup = str_is_1char (snip, SNIP_LOOKUP);

    // we can't drop b250 if we have local data (as reconstructor will reconstruct from local rather than dict), unless it is a simple SNIP_LOOKUP
    if (vctx->local.len && !is_simple_lookup) NO_DROP ("ctx has local, but snip is not SNIP_LOOKUP"); 

    struct FlagsCtx my_flags = vctx->flags;
    struct FlagsCtx vb_1_flags = (vb->vblock_i==1) ? (struct FlagsCtx){} : zctx->flags; // flags are set by vb_i=1 
    my_flags.all_the_same = vb_1_flags.all_the_same = false; // compare flags, except for all_the_same

    // we can't drop if vctx has flags that need to be passed to piz (unless its vb_i>1 and flags are same as vb_i=1)
    if (((SectionFlags)my_flags).flags != ((SectionFlags)vb_1_flags).flags) NO_DROP (vb->vblock_i==1 ? "this is vb_i=1, and vctx has flags" : "vctx has flags, different from those set by vb_i=1");

    // we survived all tests - we can now drop the b250 and maybe also dict

    if (flag.debug_generate) 
        iprintf ("%s: %s is \"all_the_same\" - dropped %sb250 b250.len=%"PRIu64"\n", VB_NAME, vctx->tag_name, (is_simple_lookup ? "dict AND " : ""), vctx->b250.len);
    
    buf_free (vctx->b250);
    if (is_simple_lookup) buf_free (vctx->dict); // if this is a SNIP_LOOKUP, we don't need the new dict entry added in this VB, as absent dict, lookup will happen 
    return;

no_drop:
    if (flag.debug_generate) 
        iprintf ("%s: %s is all_the_same but cannot drop b250 because %s\n", VB_NAME, vctx->tag_name, reason);
}

// in case we're dropping vctx after already merging it - substract txt_len added in ctx_merge_in_one_vctx
void ctx_substract_txt_len (VBlockP vb, ContextP vctx)
{
    ContextP zctx = ctx_get_zctx_from_vctx (vctx, false, false);
    if (!zctx) return;

    mutex_lock (ZMUTEX(zctx));
    zctx->txt_len -= vctx->txt_len;
    mutex_unlock (ZMUTEX(zctx));
}

static inline bool vctx_needs_merge (VBlockP vb, ContextP vctx)
{
    return vctx->dict_id.num && 
                (  vctx->nodes.len32                     // nodes need merging (some might pre-populated without b250)
                || vctx->ol_nodes.len32                  // possibly need to update counts 
                || vctx->b250.len32                      // some b250 might not have nodes (eg WORD_INDEX_MISSING)
                || vctx->txt_len                         // a context may have txt_len but the txt was segged to another context
                || (vctx->local.len && vb->vblock_i==1)  // zctx->flags needs setting
                || vctx->st_did_i != DID_NONE            // st_did_i needs setting
                || vctx->is_stats_parent);               // is_stats_parent needs setting
}

// run by the vb_i=1 thread
static inline uint8_t *ctx_get_vb_1_alias_count (VBlockP vb)
{
    uint8_t *alias_count = CALLOC (MAX_DICTS);

    // set alias count to the number of vctxs which have nodes that map to this zctx
    for_vctx_that (vctx_needs_merge (vb, vctx)) 
        alias_count[ctx_get_zctx_from_vctx (vctx, true, true/*follow aliases*/)->did_i]++;

    // unlock all vb_1 locks for which vb_1 has no nodes as subsequent VBs might have and need to merge
    for (Did did_i=0; did_i < MAX_DICTS; did_i++) 
        if (!alias_count[did_i])
            mutex_unlock (z_file->wait_for_vb_1_mutex[did_i]);

    return alias_count;
}

// increment counts, where increment may or may not have the protection bit. if it does, it sets the
// protection bit of the counter.
static inline void add_count (uint64_t *counter, uint32_t increment)
{
    if (increment & COUNT_PROTECTED_FROM_REMOVAL32) {
        *counter += increment & ~COUNT_PROTECTED_FROM_REMOVAL64; 
        *counter |= COUNT_PROTECTED_FROM_REMOVAL64; // this bit needs to be ORed, not added.
    }
    else
        *counter += increment;
}

// ZIP only: this is called towards the end of compressing one vb - merging its dictionaries into the z_file 
// each dictionary is protected by its own mutex, and there is one z_file mutex protecting num_dicts.
// we are careful never to hold two muteces at the same time to avoid deadlocks
static inline bool ctx_merge_in_one_vctx (VBlockP vb, ContextP vctx, uint8_t *vb_1_alias_count, ContextP *zctx_out)
{
    // get the zctx or create a new one. note: ctx_add_new_zf_ctx() must be called before mutex_lock() because it locks the z_file mutex (avoid a deadlock)
    ContextP zctx       = ctx_get_zctx_from_vctx (vctx, true, true);
    ContextP zctx_alias = ctx_get_zctx_from_vctx (vctx, true, false);

    if ((vb->vblock_i != 1 && !mutex_wait (z_file->wait_for_vb_1_mutex[zctx->did_i], false)) || // let vb_i=1 merge first, as it has the sorted dictionaries, other vbs can go in arbitrary order. 
        !mutex_trylock (ZMUTEX(zctx))) // also implies locking all its aliases including zctx_alias
        return false; // zctx is currently locked by another VB 

    //iprintf ( ("Merging dict_id=%.8s into z_file vb_i=%u vb_did_i=%u z_did_i=%u\n", dis_dict_id (vctx->dict_id).s, vb->vblock_i, did_i, z_did_i);

    ctx_drop_all_the_same (vb, zctx, vctx); // drop b250 and maybe also vctx->dict, if warranted
    
    zctx->merge_num++; // first merge is #1 (first clone which happens before the first merge, will get vb-)
    zctx->num_new_entries_prev_merged_vb = vctx->nodes.len; // number of new words in this dict from this VB
    
    zctx->counts_section |= vctx->counts_section;   // for use of ctx_compress_counts
    zctx_alias->counts.count += vctx->counts.count; // context-specific counter (not passed to PIZ)

    if (vb->vblock_i == 1 && (vctx->b250.len || vctx->local.len))  // vb=1 must have either a b250 or local section to carry the flags, otherwise the default flags are 0
        zctx_alias->flags = vctx->flags; // vb_1 flags will be the default flags for this context, used by piz in case there are no b250 or local sections due to all_the_same. see zip_generate_b250_section and piz_read_all_ctxs

    uint32_t ol_len = vctx->ol_nodes.len32;
    bool has_count = !DTP(zip_vb_has_count) || DTP(zip_vb_has_count)(vb); // normally we count, but callback can override

    if ((flag.show_stats_comp_i == COMP_NONE && (vb->data_type != DT_VCF || vb->comp_i != VCF_COMP_PRIM_ONLY)) // we don't include ##primary_only VBs as they are not in the primary reconstruction, but we do include ##luft_only
        || flag.show_stats_comp_i == vb->comp_i)
        zctx_alias->txt_len += vctx->txt_len; // for stats

    if (vctx->st_did_i != DID_NONE && zctx_alias->st_did_i == DID_NONE) {
        ContextP st_ctx = ctx_get_zctx_from_vctx (CTX(vctx->st_did_i), false, false);
        if (st_ctx) zctx_alias->st_did_i = st_ctx->did_i; // st_did_i is not necessarily the same for vb and zf
    }

    if (vctx->is_stats_parent)
        zctx_alias->is_stats_parent = true; // we set, but we never revert back

    // we assign VB a codec from zctx, if not already assigned by Seg. See comment in codec_assign_best_codec
    if (!vctx->lcodec) vctx->lcodec = zctx_alias->lcodec;
    if (!vctx->bcodec) vctx->bcodec = zctx_alias->bcodec;
    
    if (!buf_is_alloc (&vctx->dict)) goto finish; // no new snips introduced in this VB
 
    // case: this is first vb that is merging this dict: allocate hash tabl based on the statistics gathered by this VB. 
    if (!buf_is_alloc (&zctx->global_hash)) {
        uint32_t estimated_entries = hash_get_estimated_entries (vb, zctx, vctx);
        hash_alloc_global (zctx, estimated_entries);
    }

    // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
    // (vctx->nodes contains only new words, old words from previous vbs are in vctx->ol_nodes)
    for_buf2 (CtxNode, vb_node, i, vctx->nodes) {
        if (vb_node->node_index == WORD_INDEX_NONE) continue; // canceled in ctx_rollback

        rom snip = Bc (vctx->dict, vb_node->char_index);
        uint32_t count = *B32(vctx->counts, ol_len + i);
        bool is_new;

        // use evb and not vb because zf_context is z_file (which belongs to evb)
        CtxNode *zf_node;
        WordIndex zf_node_index = 
            ctx_commit_node (vb, zctx, vctx, snip, vb_node->snip_len, count, &zf_node, &is_new); // also allocs nodes, counts and chrom2ref_map

        ASSERT (zf_node_index >= 0 && zf_node_index < zctx->nodes.len32, 
                "zf_node_index=%d out of range - len=%i", zf_node_index, vctx->nodes.len32);

        if (has_count) 
            add_count (B64(zctx->counts, zf_node_index), count);

        // set word_index to be indexing the global dict - to be used by zip_generate_b250_section()
        if (is_new) {
            // note: chrom2ref_map is protected by ZMUTEX(zctx) (zctx is CHROM as we followed alias)
            if (chrom_2ref_seg_is_needed(zctx->did_i) && vctx->chrom2ref_map.len) { // SAM: it is ctx->chrom2ref_map.len=0 and nodes.len>0: in case of gencomp_len=2 will will have a SPECIAL snip node, but we don't use chrom2ref as all contigs must be in the SAM header
                buf_alloc_255 (evb, &ZCTX(CHROM)->chrom2ref_map, 0, MAX_(INITIAL_NUM_NODES, zf_node_index+1), WordIndex, CTX_GROWTH, "ZCTX(CHROM)->chrom2ref_map");
                *B(WordIndex, ZCTX(CHROM)->chrom2ref_map, zf_node_index) = *B(WordIndex, vctx->chrom2ref_map, i);
                ZCTX(CHROM)->chrom2ref_map.len32 = MAX_(ZCTX(CHROM)->chrom2ref_map.len32, zf_node_index + 1);
            }

            vb_node->word_index = zf_node->word_index = zf_node_index;
        } else 
            // a previous VB already already calculated the word index for this node. if it was done by vb_i=1,
            // then it is also re-sorted and the word_index is no longer the same as the node_index
            vb_node->word_index = zf_node->word_index;
    }

    // warn if dict size is excessive
    if (zctx->dict.len > EXCESSIVE_DICT_SIZE && !zctx->dict_len_excessive) {
        zctx->dict_len_excessive = true; // warn only once (per context)
        WARN ("WARNING: excessive zctx dictionary size - causing slow compression and decompression and reduced compression ratio. Please report this to " EMAIL_SUPPORT ".\n"
              "%s %s data_type=%s ctx=%s vb=%s vb_size=%"PRIu64" zctx->dict.len=%"PRIu64" version=%s. First 1000 bytes: ", 
              cond_str (VB_DT(BAM) || VB_DT(SAM), "sam_mapper=", segconf_sam_mapper_name()), 
              cond_str (VB_DT(BAM) || VB_DT(SAM) || VB_DT(FASTQ) || VB_DT(KRAKEN), "segconf_qf_name=", segconf_qf_name(QNAME1)), 
              z_dt_name(), zctx->tag_name, VB_NAME, segconf.vb_size, zctx->dict.len, GENOZIP_CODE_VERSION);
        str_print_dict (stderr, zctx->dict.data, 1000, false, false);
    }

    if (flag.deep) {
        if (VB_DT(BAM) || VB_DT(SAM)) zctx->dict_flags.deep_sam = true;
        else if (VB_DT(FASTQ))        zctx->dict_flags.deep_fastq = true;
    }

finish:
    // just update counts for ol_node (i.e. known to be existing) snips
    if (has_count) {
        ARRAY32 (uint32_t, vcounts, vctx->counts);
        ARRAY32 (uint64_t, zcounts, zctx->counts);

        ASSERT (zcounts_len >= ol_len, "%s ctx=%s: Expecting zcounts_len=%u >= ol_len=%u", 
                VB_NAME, zctx->tag_name, zcounts_len, ol_len);

        for (WordIndex ni=0; ni < ol_len; ni++) 
            add_count (&zcounts[ni], vcounts[ni]);
    }

    mutex_unlock (ZMUTEX(zctx));

    if (vb->vblock_i == 1) { 
        ASSERT (vb_1_alias_count[zctx->did_i], "Unexpectedly vb_1_alias_count[%s]=0", zctx->tag_name);

        if (--vb_1_alias_count[zctx->did_i] == 0) 
            mutex_unlock (z_file->wait_for_vb_1_mutex[zctx->did_i]); 
    }

    *zctx_out = zctx;
    return true; // merge was done
}

// ZIP only: merge new words added in this vb into the z_file.contexts, and compresses dictionaries.
void ctx_merge_in_vb_ctx (VBlockP vb)
{
    START_TIMER;

    // note: vb_1_alias_count is only accessed by the vb_i=1 thread
    uint8_t *vb_1_alias_count = (vb->vblock_i == 1) ? ctx_get_vb_1_alias_count (vb) : NULL;

    // merge all contexts 
    bool all_merged=false;
    bool custom_merge_pending = !!DTP(zip_custom_merge);
    
    ContextP v_did_i_to_zctx[vb->num_contexts];
    memset (v_did_i_to_zctx, 0, vb->num_contexts * sizeof(ContextP));

    while (!all_merged) {

        all_merged=true; // unless proven otherwise
        bool any_merged=false;

        for_vctx_that (!vctx->dict_merged && vctx_needs_merge (vb, vctx)) {
            vctx->dict_merged = ctx_merge_in_one_vctx (vb, vctx, vb_1_alias_count, &v_did_i_to_zctx[vctx - CTX(0)]); // false if zctx is locked by another VB
            all_merged &= vctx->dict_merged;
            any_merged |= vctx->dict_merged;
        }

        // data-type specific non-context merge. advantage of running logic here vs zip_after_compress is that it merges
        // contexts while the custom mutex is locked by another thread.
        if (custom_merge_pending) {
            if (mutex_trylock (z_file->custom_merge_mutex)) {
                DTP(zip_custom_merge)(vb);
                mutex_unlock (z_file->custom_merge_mutex);
                any_merged = true;
                custom_merge_pending = false;
            }

            else
                all_merged = false;
        }
        
        // case: couldn't merge any - perhaps the one large final vctx is busy - sleep a bit rather than busy-wait
        if (!any_merged && !all_merged) {
            START_TIMER;
            usleep (100); 
            COPY_TIMER (wait_for_merge);
        }
    }

    // note: z_file->num_contexts might be larger than vb->num_contexts at this point, for example:
    // vb_i=1 started, z_file is empty, created 20 contexts
    // vb_i=2 started, z_file is empty, created 10 contexts
    // vb_i=1 completes, merges 20 contexts to z_file, which has 20 contexts after
    // vb_i=2 completes, merges 10 contexts, of which 5 (for example) are shared with vb_i=1. Now z_file has 25 contexts after.

    if (vb->vblock_i == 1) {
        for (Did did_i=0; did_i < MAX_DICTS; did_i++)
            ASSERT (!vb_1_alias_count[did_i], "Expecting vb_1_alias_count[%s] to be 0, but it is %u", ZCTX(did_i)->tag_name, vb_1_alias_count[did_i]);
        FREE (vb_1_alias_count);
    }

    COPY_TIMER (ctx_merge_in_vb_ctx);
}

static Did ctx_did_i_search (const ContextIndex *ctx_index, Did num_contexts, DictId dict_id, Did first, Did last)
{
    if (last < first) return DID_NONE;

    Did mid = (first + last) / 2;
    DictId dict_id_mid = ctx_index[mid].dict_id;

    if (dict_id_mid.num == dict_id.num)    
        return ctx_index[mid].did_i;

    else if (dict_id_mid.num < dict_id.num) 
        return ctx_did_i_search (ctx_index, num_contexts, dict_id, mid+1, last);

    else
        return ctx_did_i_search (ctx_index, num_contexts, dict_id, first, mid-1);
}

// returns an existing did_i in this vb, or DID_NONE if there isn't one
Did ctx_get_unmapped_existing_did_i (const ContextArray contexts, const ContextIndex *ctx_index, Did num_contexts, DictId dict_id)
{
    int did_i; // signed

    // binary search if we have ctx_index (we will have it in PIZ compute threads)
    if (ctx_index) {
        if ((did_i = ctx_did_i_search (ctx_index, num_contexts, dict_id, 0, num_contexts-1)) != DID_NONE)
            return did_i;
    }

    else // linear search if not
        for (did_i=num_contexts-1; did_i >= 0 ; did_i--)  // Search backwards as unmapped ctxs are more likely to be towards the end.
            if (dict_id.num == contexts[did_i].dict_id.num) return contexts[did_i].did_i; // note: in case of an alias, contexts[did_i].did_i is the dst_did_i

    return DID_NONE; // not found
}

// gets did_id if the dictionary exists, and creates a new dictionary if its the first time dict_id is encountered
// threads: no issues - called by PIZ for vb and zf (but dictionaries are immutable) and by Seg (ZIP) on vctx only
ContextP ctx_get_unmapped_ctx (ContextArray contexts, DataType dt, DictIdtoDidMap d2d_map, 
                               Did *num_contexts, DictId dict_id, STRp(tag_name))
{
    // search to see if this dict_id has a context, despite not in the d2d_map (due to contention). 
    for (int/*signed*/ did_i=*num_contexts-1; did_i >= 0 ; did_i--)  // Search backwards as unmapped ctxs are more likely to be towards the end.
        if (dict_id.num == contexts[did_i].dict_id.num) 
            return &contexts[did_i];

    ContextP ctx = &contexts[*num_contexts]; 
 
    //iprintf ("New context: dict_id=%s in did_i=%u \n", dis_dict_id (dict_id).s, did_i);
    ASSERT (*num_contexts < MAX_DICTS, "cannot create a context for %.*s (dict_id=%s) because number of dictionaries would exceed MAX_DICTS=%u", 
            tag_name_len, tag_name, dis_dict_id (dict_id).s, MAX_DICTS);

    ctx_initialize_ctx (ctx, *num_contexts, dict_id, d2d_map, STRa (tag_name));

    // thread safety: the increment below MUST be AFTER the initialization of ctx, bc piz_get_line_subfields
    // might be reading this data at the same time as the piz dispatcher thread adding more dictionaries
    (*num_contexts)++;

    return ctx;
}

// called from seg_all_data_lines (ZIP) and ctx_piz_initialize_zctxs (PIZ) to initialize all
// primary field ctx's. these are not always used (e.g. when some are not read from disk due to genocat options)
// but we maintain their fixed positions anyway as the code relies on it
// Note: z_file Context Buffers are already initialized in file_initialize_z_file_data
void ctx_initialize_predefined_ctxs (ContextArray contexts, 
                                     DataType dt,
                                     DictIdtoDidMap d2d_map,
                                     Did *num_contexts)
{
    ASSERT0 (contexts == z_file->contexts, "expecting z_file contexts");
    ASSERTMAINTHREAD;

    *num_contexts = MAX_(dt_fields[dt].num_fields, *num_contexts);

    init_dict_id_to_did_map (d2d_map); // reset, in case data_type changed

    BufferP aliases = aliases_get(); // NULL if no aliases
    Did alias_dids[aliases->len]; // entry i corresponds to alias i
    memset (alias_dids, 0xff, sizeof (alias_dids)); // initialize to DID_NONE to detect duplicates

    for (Did did_i=0; did_i < dt_fields[dt].num_fields; did_i++) {
        DictId dict_id = dt_fields[dt].predefined[did_i].dict_id;
        ASSERT (dict_id.num, "No did_i->dict_id mapping is defined for predefined did_i=%u in dt=%s", did_i, dt_name (dt));

        // skip if its an alias 
        AliasType alias_type = ALIAS_NONE;
        for_buf2 (DictIdAlias, alias, alias_i, *aliases) 
            if (dict_id.num == alias->alias.num) {
                alias_dids[alias_i] = did_i;
                alias_type = alias->alias_type;
                goto breakout; // we can't "break" out of a for_buf2, using goto instead
            }
        breakout:

        // initialize, but in PIZ, not if ALIAS_CTX
        if (IS_ZIP || alias_type != ALIAS_CTX)
            ctx_initialize_ctx (&contexts[did_i], did_i, dict_id, d2d_map, 
                                dt_fields[dt].predefined[did_i].tag_name, dt_fields[dt].predefined[did_i].tag_name_len);
    }

    // initialize alias contexts (note: all aliases and their destinations are predefined)
    // note: all alias destinations that ever existed in previous versions of Genozip must be defined in #pragma GENDICT for this to work
    for_buf2 (DictIdAlias, alias, alias_i, *aliases) {
        ASSERT0 (alias_dids[alias_i] != DID_NONE, "Duplicate alias detected"); // a missing alias_did, means another entry got set twice
        
        ContextP alias_ctx = &z_file->contexts[alias_dids[alias_i]];
        ContextP dst_ctx   = ctx_get_existing_zctx (alias->dst); 
        
        ASSERT (dst_ctx, "destination %s of alias %s->%s is missing a #pragma GENDICT", 
                dis_dict_id (alias->dst).s, dis_dict_id (alias->alias).s, dis_dict_id (alias->dst).s);

        switch (alias->alias_type) {
            case ALIAS_CTX: // initialize as alias - instead of ctx_initialize_ctx
                if (IS_PIZ) { // note: in ZIP we seg directly to the destination, alias contexts are not used
                    alias_ctx->did_i        = dst_ctx->did_i; // dst did_i instead of our own did_i
                    alias_ctx->dict_id      = alias->alias;
                    alias_ctx->is_ctx_alias = true;
                    strcpy ((char*)alias_ctx->tag_name, dis_dict_id (alias_ctx->dict_id).s);

                    set_d2d_map (d2d_map, alias->alias, dst_ctx->did_i);
                }
                break;

            case ALIAS_DICT:
                alias_ctx->dict_did_i = dst_ctx->did_i;
                break;

            default:
                ABORT ("Invalid alias_type=%u", alias->alias_type);
        }
    }

    buf_destroy (*aliases);
}

// PIZ only: this is called by the main thread after it integrated all the dictionary fragment read from disk for one VB.
// Here we hand over the integrated dictionaries to the VB - in preparation for the Compute Thread to use them.
// We overlay the z_file's dictionaries and word lists to the vb. these data remain unchanged - neither
// the vb nor the dispatcher thread will ever change snips placed in these. the dispatcher thread may append
// the dictionary and word list as new fragments become available from subsequent VBs. If the memory is not 
// sufficient, the dispatcher thread will "abandon" this memory, leaving it to the VB to continue to use it
// while starting a larger dict/word_list on a fresh memory allocation.
void ctx_overlay_dictionaries_to_vb (VBlockP vb)
{
    for (Did did_i=0; did_i < z_file->num_contexts; did_i++) {
        ContextP zctx = ZCTX(did_i);
        ContextP vctx = CTX(did_i);

        if (!zctx->dict_id.num) continue;

        // we create a VB contexts even if there are no dicts (perhaps skipped due to flag) - needed for containers to work
        vctx->did_i        = zctx->did_i;      // note: in case of an alias, this will be did_i != zctx->did_i
        vctx->dict_id      = zctx->dict_id;
        vctx->is_loaded    = zctx->is_loaded;  // for now, we know that dictionary loaded. If not set here, it might still be set in piz_read_all_ctxs
        vctx->is_ctx_alias = zctx->is_ctx_alias;
        vctx->other_did_i  = DID_NONE;
        vctx->last_line_i  = LAST_LINE_I_INIT;
        vctx->pair_assist_type = SEC_NONE;

        memcpy ((char*)vctx->tag_name, zctx->tag_name, sizeof (vctx->tag_name));

        set_d2d_map (vb->d2d_map, vctx->dict_id, did_i);

        ctx_init_iterator (vctx);

        ContextP dict_ctx = ZCTX(zctx->dict_did_i); // this is either our zctx, or if we're a ALIAS_DICT - our destination's context

        if (buf_is_alloc (&dict_ctx->dict))
            buf_overlay (vb, &vctx->dict, &dict_ctx->dict, "ctx->dict");    
        
        if (buf_is_alloc (&dict_ctx->word_list))
            buf_overlay (vb, &vctx->word_list, &dict_ctx->word_list, "ctx->word_list");
    }

    vb->num_contexts = z_file->num_contexts;
}

// used by random_access_show_index
CtxNode *ctx_get_node_by_word_index (ConstContextP ctx, WordIndex word_index)
{
    for_buf (CtxNode, node, ctx->nodes)
        // if (node->word_index.n == word_index) return node;
        if (node->word_index == word_index) return node;

    ABORT ("ctx_get_node_by_word_index failed to find word_index=%d in did_i=%u", word_index, ctx->did_i);
}

// PIZ: get snip by normal word index (doesn't support WORD_INDEX_* - returns "")
rom ctx_get_snip_by_word_index_do (ConstContextP ctx, WordIndex word_index, STRp(*snip), FUNCLINE)
{
    ASSERT (buf_is_alloc (&ctx->word_list), "called from %s:%u: word_list is not allocated for ctx=%s", func, code_line, ctx->tag_name);

    ASSERT ((uint32_t)word_index < ctx->word_list.len32, "called from %s:%u: word_index=%d out of range: word_list.len=%u for ctx=%s",
            func, code_line, word_index, ctx->word_list.len32, ctx->tag_name);

    if (word_index < 0) {
        static rom empty="";
        if (snip) *snip = empty;
        if (snip_len) *snip_len = 0;
        return empty;
    }

    CtxWord *word = B(CtxWord, ctx->word_list, word_index);
    rom my_snip = B(const char, ctx->dict, word->index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = word->len;

    return my_snip; 
}

// PIZ: returns word index of the snip, or WORD_INDEX_NONE if it is not in the dictionary
WordIndex ctx_get_word_index_by_snip (VBlockP vb, ContextP ctx, STRp(snip))
{
    if (!snip_len) snip_len = strlen (snip);

    ARRAY (const CtxWord, words, ctx->word_list);
    ARRAY (const char, dict, ctx->dict);

    if (!ctx->piz_word_list_hash.data) 
        buf_alloc_255 (vb, &ctx->piz_word_list_hash, 0, 5*ctx->word_list.len, WordIndex, 0, "piz_word_list_hash");

    WordIndex *hash_ent = B(WordIndex, ctx->piz_word_list_hash, hash_do (5*ctx->word_list.len, STRa(snip)));
    if (*hash_ent != WORD_INDEX_NONE && str_issame_(STRa(snip), &dict[words[*hash_ent].index], words[*hash_ent].len))
        return *hash_ent;

    for (WordIndex wi=0; wi < words_len; wi++)
        if (str_issame_(STRa(snip), &dict[words[wi].index], words[wi].len)) {
            *hash_ent = wi;
            return wi;
        }

    return WORD_INDEX_NONE;
}

// ZIP
rom ctx_snip_from_zf_nodes (ConstContextP zctx, WordIndex node_index, pSTRp(snip))
{
    ASSERT (node_index >= 0 && node_index < zctx->nodes.len32, 
            "node_index=%d out of range, nodes.len=%u", node_index, zctx->nodes.len32);

    CtxNode *node = B(CtxNode, zctx->nodes, node_index);
    rom my_snip = Bc (zctx->dict, node->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = node->snip_len;

    return my_snip; 
}

// get the node_index=word_index of a snip from zctx. its an error if the snip doesn't exist.
WordIndex ctx_get_ol_node_index_by_snip (VBlockP vb, ContextP ctx, STRp(snip))
{
    CtxNode *node = NULL;
    return hash_get_entry_for_seg (vb, ctx, STRa(snip), WORD_INDEX_NONE/*read-only*/, &node);
}

static BufferP sorter_cmp_counts = NULL; // for use by sorter_cmp - used only in vblock_i=1, so no thread safety issues
static SORTER (sorter_cmp)  
{ 
    return DESCENDING_RAW (*B32(*sorter_cmp_counts, *(WordIndex *)a), 
                           *B32(*sorter_cmp_counts, *(WordIndex *)b));
}

void ctx_sort_dictionaries_vb_1 (VBlockP vb)
{
    START_TIMER;

    BufferP sorter = &vb->codec_bufs[0];
    BufferP unsorted_dict = &vb->codec_bufs[1];
    ASSERTNOTINUSE (*sorter);
    ASSERTNOTINUSE (*unsorted_dict);
    
    // thread safety note: no issues here, as this is run only by the compute thread of vblock_i=1
    for_ctx {
        if (ctx->no_vb1_sort || !ctx->nodes.len || chrom_2ref_seg_is_needed(ctx->did_i)) continue;
        
        // prepare sorter array containing indices into ctx->nodes. We are going to sort it rather than sort nodes directly
        // as the b250 data contains node indices into ctx->nodes.
        buf_alloc (vb, sorter, 0, ctx->nodes.len, WordIndex, CTX_GROWTH, "vb_1_sorter");
        for (WordIndex i=0; i < ctx->nodes.len32; i++)
            BNXT (WordIndex, *sorter) = i;

        // sort in ascending order of nodes->count
        sorter_cmp_counts = &ctx->counts; // communicate the ctx to sorter_cmp via a global var
        qsort (sorter->data, ctx->nodes.len, sizeof (WordIndex), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in nodes
        buf_move (vb, *unsorted_dict, "unsorted_dict", ctx->dict);
        buf_alloc_exact (vb, ctx->dict, unsorted_dict->len, char, "contexts->dict");

        // note: we sort dict and update nodes->char_index and word_index. nodes and counts are not sorted.
        char *next = B1STc (ctx->dict);
        for (WordIndex i=0; i < (WordIndex)ctx->nodes.len32; i++) {
            WordIndex node_index = *B(WordIndex, *sorter, i);
            CtxNode *node = B(CtxNode, ctx->nodes, node_index);

            if (node->node_index == WORD_INDEX_NONE) continue; // node was canceled in ctx_rollback

            rom snip = Bc(*unsorted_dict, node->char_index); 

            node->char_index = BNUM64 (ctx->dict, next);
            node->node_index = i;
            next = mempcpy (next, snip, node->snip_len + 1 /* +1 for SNIP_SEP */);

            if (HAS_DEBUG_SEG(ctx)) {
                char printable_snip[node->snip_len+20];
                iprintf ("ctx_sort_dictionaries_vb_1: %s: word_index=%u snip=%s snip_len=%u count=%u\n",
                          ctx->tag_name, i, str_print_snip (snip, node->snip_len, printable_snip), node->snip_len, 
                          *B32(ctx->counts, ctx->ol_nodes.len + node_index));
            }
        }

        buf_free (*sorter);
        buf_free (*unsorted_dict);
    }

    COPY_TIMER(ctx_sort_dictionaries_vb_1);
}

// ZIP only: run by main thread during zfile_output_processed_vb()
void ctx_update_stats (VBlockP vb)
{
    for_vctx {    
        ContextP zctx = ctx_get_zctx_from_vctx (vctx, false, false);
        if (!zctx) continue; // this can happen if FORMAT subfield appears, but no line has data for it

        zctx->b250.count      += vctx->b250.count;      // number of segs into b250
        zctx->local_num_words += vctx->local_num_words; // number of segs into local
        zctx->local.len       += vctx->local.len * lt_width(vctx); // uncompressed size of local

        // fields segged of this type in the file - if we have both, take the MAX. cases:
        // - all fields have b250, some with look up and local too, some without local (max is b250)
        // - all fields are lookup - b250 is 1 (all the same), all fields in local (take local)
        zctx->word_list.count += MAX_(vctx->local_num_words, vctx->b250.count);
    }
}

void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local /* true = local, false = b250 */)
{
    char dump_fn[MAX_TAG_LEN + 50];
    sprintf (dump_fn, "%s.%05u.%s", ctx->tag_name, vb->vblock_i, local ? "local" : "b250");
    
    bool success = local ? buf_dump_to_file (dump_fn, &ctx->local, lt_width(ctx), false, true, true, false)
                         : buf_dump_to_file (dump_fn, &ctx->b250, 1, false, true, true, false);

    ASSERTW (success, "Warning: ctx_dump_binary failed to output file %s: %s", dump_fn, strerror (errno));
}

// rewrite unused (i.e count=0) zctx dict words as "" (we don't delete them completely, so word_index's in b250 sections remain correct)
void ctx_shorten_unused_dict_words (Did did_i)
{
    ContextP zctx = ZCTX(did_i);
    ARRAY (CtxNode,  nodes,  zctx->nodes);
    ARRAY (uint64_t, counts, zctx->counts);
    ARRAY (char, dict, zctx->dict);

    // note that dict words are not in the order of the nodes array. 
    // pass 1: mark characters for deletion
    uint64_t deleted_so_far = 0;
    for (WordIndex ni=0; ni < nodes_len; ni++) {
        uint64_t new_char_index = nodes[ni].char_index - deleted_so_far;

        if (!counts[ni]) {
            memset (&dict[nodes[ni].char_index], SNIP_RESERVED, nodes[ni].snip_len); // A value guaranteed not to exist in dictionary data
            deleted_so_far += nodes[ni].snip_len; 
            nodes[ni].snip_len = 0;
        }

        nodes[ni].char_index = new_char_index;
    }
        
    // pass 2: delete characters
    char *next = dict;
    for (uint64_t i=0; i < dict_len; i++)
        if (dict[i] != SNIP_RESERVED)
            *next++ = dict[i];
    
    zctx->dict.len = next - dict;
}

// -----------------------------
// ZIP & PIZ: SEC_COUNTS sections
// -----------------------------

typedef struct {
    int64_t count;
    WordIndex word_index;
    rom snip;
} ShowCountsEnt;

static DESCENDING_SORTER (show_counts_cmp, ShowCountsEnt, count)

static void ctx_show_counts (ContextP zctx)
{
    static Buffer show_counts_buf = {};
    buf_free (show_counts_buf);
    buf_alloc (evb, &show_counts_buf, 0, zctx->counts.len, ShowCountsEnt, 0, "show_counts_buf");

    // QUAL counts store Longr value-to-bin mapping
    bool maybe_longr = ((Z_DT(BAM) || Z_DT(SAM)) && (zctx->dict_id.num == _SAM_DOMQRUNS || zctx->dict_id.num == _OPTION_OQ_DOMQRUNS || zctx->dict_id.num == _OPTION_U2_DOMQRUNS))
                    || (Z_DT(FASTQ) && zctx->dict_id.num == _FASTQ_DOMQRUNS);

    uint64_t total=0;
    for (uint32_t i=0; i < zctx->counts.len32; i++) {
        uint64_t count = *B64(zctx->counts, i) & ~COUNT_PROTECTED_FROM_REMOVAL64;
        if (!count && !maybe_longr) continue;

        total += count;
        BNXT (ShowCountsEnt, show_counts_buf) = (ShowCountsEnt){ 
            .count = count,
            .word_index = i,
            .snip  = maybe_longr     ? "" 
                   : (command==ZIP)  ? ctx_snip_from_zf_nodes (zctx, i, 0, 0) 
                   :                   ctx_get_words_snip (zctx, i)
        };
    }

    ARRAY (ShowCountsEnt, counts, show_counts_buf);
    
    if (!maybe_longr)
        qsort (counts, counts_len, sizeof (ShowCountsEnt), show_counts_cmp);

    iprintf ("Showing counts of %s (did_i=%u). Total items=%"PRIu64" Number of categories=%u\n", zctx->tag_name, zctx->did_i, total, (unsigned)counts_len);    

    if (total)
        for (uint32_t i=0; i < counts_len; i++) 
            iprintf ("\"%s\"(%d)\t%"PRIu64"\t%-4.2f%%\n", counts[i].snip, counts[i].word_index, counts[i].count, 
                        100 * (float)counts[i].count / (float)total);

    if (is_genocat) exit_ok;
}

rom ctx_get_snip_with_largest_count (Did did_i, int64_t *count)
{
    ContextP zctx = ZCTX(did_i);
    ARRAY (CtxNode,  nodes,  zctx->nodes);
    ARRAY (uint64_t, counts, zctx->counts);

    *count = -1;
    rom snip = "";

    for (WordIndex i=0; i < nodes_len; i++)
        if (counts[i] > *count) {
            *count = counts[i] & ~COUNT_PROTECTED_FROM_REMOVAL64;
            snip = Bc (zctx->dict, nodes[i].char_index);
        }

    return snip;
}

// ZIP
rom ctx_get_vb_snip_ex (ConstContextP vctx, WordIndex vb_node_index, pSTRp(snip))
{
    CtxNode *node;
    rom out_snip;

    if (vb_node_index < vctx->ol_nodes.len32) { // is this entry from a previous vb (overlay buffer)
        node = B(CtxNode, vctx->ol_nodes, vb_node_index);
        out_snip = Bc (vctx->ol_dict, node->char_index);
    }
    else {
        ASSERT (vb_node_index - vctx->ol_nodes.len32 < vctx->nodes.len32, 
                "Expecting node_index=%u < ol_nodes.len=%u + nodes.len=%u", 
                vb_node_index, vctx->ol_nodes.len32, vctx->nodes.len32);

        node = B(CtxNode, vctx->nodes, vb_node_index - vctx->ol_nodes.len32);
        out_snip = Bc (vctx->dict, node->char_index);
    }

    if (snip) *snip = out_snip;
    if (snip_len) *snip_len = node->snip_len;

    return out_snip;
}

void ctx_compress_counts (void)
{
    START_TIMER;

    for_zctx {
        if (flag.show_one_counts.num == dict_id_typeless (zctx->dict_id).num) 
            ctx_show_counts (zctx);

        if (zctx->counts_section && zctx->counts.len) {

            // remove protection bit
            for_buf (uint64_t, count, zctx->counts)
                *count &= ~COUNT_PROTECTED_FROM_REMOVAL64;
                
            BGEN_u64_buf (&zctx->counts, NULL);       

            zctx->counts.len *= sizeof (uint64_t);

            Codec codec = codec_assign_best_codec (evb, NULL, &zctx->counts, SEC_COUNTS);

            SectionHeaderCounts header = (SectionHeaderCounts){ 
                .magic                 = BGEN32 (GENOZIP_MAGIC),
                .section_type          = SEC_COUNTS,
                .data_uncompressed_len = BGEN32 (zctx->counts.len32),
                .codec                 = codec,
                .vblock_i              = 0,
                .nodes_param           = BGEN64 (zctx->nodes.param),
                .dict_id               = zctx->dict_id
            };

            comp_compress (evb, zctx, &evb->z_data, &header, zctx->counts.data, NO_CALLBACK, "SEC_COUNTS");
            zctx->counts.len /= sizeof (uint64_t);

            BGEN_u64_buf (&zctx->counts, NULL); // we need it for stats
        }
    }

    COPY_TIMER_EVB (ctx_compress_counts);
}

// PIZ: called by main thread after reading GENOZIP_HEADER - create predefined contexts as well as all
// contexts that exist in the file, and mark those with data in the file as "z_data_exists". 
// z_data_exists can be relied on in flags_update_piz_one_z_file and IS_SKIP functions.
void ctx_piz_initialize_zctxs (void)
{
    if (z_file->num_contexts) return; // already initialized (this happens, bc zfile_read_genozip_header is called multiple times)

    ctx_initialize_predefined_ctxs (z_file->contexts, z_file->data_type, z_file->d2d_map, &z_file->num_contexts);

    Section sec = NULL;
    while (sections_next_sec3 (&sec, SEC_B250, SEC_LOCAL, SEC_DICT)) {
        ContextP zctx = ctx_get_existing_zctx (sec->dict_id);

        // case: first encounter in z_file with this non-predefined dict_id - initialize a ctx for it
        if (!zctx) { 
            zctx = &z_file->contexts[z_file->num_contexts++];
            ctx_initialize_ctx (zctx, z_file->num_contexts-1, sec->dict_id, z_file->d2d_map, 0, 0);
        }

        if (!zctx->z_data_exists) // predefined not encountered before, or non-predefined just initialized
            zctx->z_data_exists = true; // store to memory only in rare cases of new encounter
    }
}

// PIZ: called by the main threads from piz_read_global_area
void ctx_read_all_counts (void)
{
    Section sec = NULL;
    bool counts_shown=false;
    while (sections_next_sec (&sec, SEC_COUNTS))  {

        ContextP zctx = ctx_get_existing_zctx (sec->dict_id);
        
        zfile_get_global_section (SectionHeaderCounts, sec, &zctx->counts, "counts");
        if (flag.only_headers || !zctx->counts.len) continue; // only show headers, or section skipped

        zctx->counts.len /= sizeof (uint64_t);
        BGEN_u64_buf (&zctx->counts, NULL);

        zctx->nodes.param = BGEN64 (header.nodes_param);
        
        if (flag.show_one_counts.num == dict_id_typeless (sec->dict_id).num) {
            ctx_show_counts (zctx);
            counts_shown=true;
        }
    }

    ASSINP (counts_shown || !flag.show_one_counts.num || !is_genocat, "There is no SEC_COUNTS section for %s", dis_dict_id (flag.show_one_counts).s);
}

TagNameEx ctx_tag_name_ex (ConstContextP ctx)
{
    TagNameEx s = {};

    if (!ctx) return (TagNameEx){ .s = "<none>" };

    unsigned start=0;
    if (z_file && (Z_DT(VCF) || Z_DT(BCF))) {
        if      (dict_id_is_vcf_format_sf (ctx->dict_id)) { memcpy (s.s, "FORMAT/", 7); start = 7; }
        else if (dict_id_is_vcf_info_sf   (ctx->dict_id)) { memcpy (s.s, "INFO/",   5); start = 5; }
    }

    memcpy (&s.s[start], ctx->tag_name, MAX_TAG_LEN);

    return s;
}

// rolls back a context to the rollback point registered in ctx_set_rollback
void ctx_rollback (VBlockP vb, ContextP ctx, bool override_id)
{
    ASSERT (override_id || ctx->rback_id == vb->rback_id, "Expected ctx->rback_id=%"PRId64" == vb->rback_id=%"PRId64, 
            ctx->rback_id, vb->rback_id);
             
    if (HAS_DEBUG_SEG(ctx)) 
        iprintf ("%s: ctx_rollback: %s: rolling back b250=%u nodes=%u\n", VB_NAME, ctx->tag_name, 
                 ctx->b250.len32 - ctx->rback_b250_len, ctx->nodes.len32 - ctx->rback_nodes_len);

    // if we evaluated this context since the rollback count - undo
    while (ctx->b250.len32 > ctx->rback_b250_len) { 
        WordIndex node_index = LASTb250(ctx);
        if (node_index >= 0)
            (*B32(ctx->counts, node_index))--;
        ctx->b250.len32--;
    }

    // update all_the_same - true, if len is down to 1, and false if it is down to 0
    if (ctx->b250.len32 <= 1) ctx->flags.all_the_same = (bool)ctx->b250.len32;

    // we can't actually remove nodes as they are refered to from the hash table, instead, we mark them as unused
    // this will prevent merging the word into the zctx->dict, and also prevent hash_get_entry_for_seg from returning this node if this word
    // is segged again in this VB - a new node will be created.
    for (WordIndex node_index_i = ctx->rback_nodes_len ; node_index_i < ctx->nodes.len32 ; node_index_i++) 
        B(CtxNode, ctx->nodes, node_index_i)->node_index = WORD_INDEX_NONE;
    
    ctx->local.len32     = ctx->rback_local_len;
    ctx->txt_len         = ctx->rback_txt_len;
    ctx->local_num_words = ctx->rback_local_num_words;
    ctx->b250.count      = ctx->rback_b250_count;
    ctx->last_value      = ctx->rback_last_value;
    ctx->last_delta      = ctx->rback_last_delta;
    ctx->last_txt        = ctx->rback_last_txt;
    ctx->last_line_i     = LAST_LINE_I_INIT; // undo "encountered in line"
    ctx->last_snip       = NULL; // cancel caching in ctx_create_node_do
    ctx->rback_id        = -1; // rolled back
}

// ZIP: get total length in VB's z_data of b250 and local of all contexts in a context group
uint64_t ctx_get_ctx_group_z_len (VBlockP vb, Did group_did_i)
{
    Did save = CTX(group_did_i)->st_did_i;
    CTX(group_did_i)->st_did_i = group_did_i; // so parent is also included in the loop

    uint64_t z_len = 0;
    for_ctx_that (ctx->st_did_i == group_did_i) 
        z_len += ctx->local_in_z_len + ctx->b250_in_z_len; 

    CTX(group_did_i)->st_did_i = save;

    return z_len;
}

// ZIP: one per file from vcf_after_vbs - set the st_did_i of the winning group, 
// and mark the directories of the losing group for deletion
void ctx_declare_winning_group (Did winning_group_did_i, Did losing_group_did_i, Did new_st_did_i) 
{
    ZCTX(winning_group_did_i)->st_did_i = winning_group_did_i; // so group head is also included in the loop
    ZCTX(losing_group_did_i) ->st_did_i = losing_group_did_i;

    for_zctx 
        if (zctx->st_did_i == winning_group_did_i) {
            zctx->st_did_i = new_st_did_i; 
            zctx->is_stats_parent = false;
        }
        else if (zctx->st_did_i == losing_group_did_i)             
            zctx->please_remove_dict = true;

    ZCTX(new_st_did_i)->is_stats_parent = true; // assuming it is predefined - no need for mutex - set once and never reset
}

// initialize "promiscuous" buffers - evb buffers that might be allocated by compute threads
// promiscuous buffers must be initialized by the main thread, and buffer.c does not verify their integrity.
void ctx_zip_init_promiscuous (ContextP ctx)
{
    ASSERTMAINTHREAD;

    #define INIT(buf) ({ buf_set_promiscuous (&ctx->buf, "contexts->" #buf);  })

    INIT(global_hash);
    INIT(global_ents);
    INIT(dict);
    INIT(nodes);
    INIT(counts);
    INIT(stons);
    INIT(ston_nodes);
}

// ZIP: called from compute thread or main thread: for predefined contexts only (i.e. ctx->did_i==zctx->did_i): 
// if we create any B250, LOCAL or DICT sections of this did_i, zctx->z_data_exists is set to true.
// Used to filter aliases written to the file.
void ctx_zip_z_data_exist (ContextP ctx)
{
    if (ctx->did_i < DTFZ(num_fields)) // case: ctx is a predefined context
        ZCTX(ctx->did_i)->z_data_exists = true;
}

// qsort comparison function
ASCENDING_SORTER (sort_by_dict_id, ContextIndex, dict_id.num)

#define MAX_MULTI_CTX_LEN 100
#define SET_MULTI_CTX(parameter_before,prop,value)      \
    va_list args;                                       \
    va_start (args, parameter_before);                  \
    unsigned d=0; for (; d < MAX_MULTI_CTX_LEN; d++) {  \
        Did did_i = (Did)va_arg (args, int);            \
        if (did_i == DID_EOL) break;                    \
        if (did_i != DID_NONE)                          \
            CTX(did_i)->prop = (value);                 \
    }                                                   \
    ASSERT (d != MAX_MULTI_CTX_LEN, "%s: Expecting last item in %s to be DID_EOL", VB_NAME, __FUNCTION__); \
    va_end (args)                                      

// ZIP, typically called from seg_initialize
void ctx_set_no_stons (VBlockP vb,                    ...) { SET_MULTI_CTX (vb, no_stons, true); }
void ctx_set_store_per_line (VBlockP vb,              ...) { SET_MULTI_CTX (vb, flags.store_per_line, true); }
void ctx_set_same_line (VBlockP vb,                   ...) { SET_MULTI_CTX (vb, flags.same_line, true); }
void ctx_set_store (VBlockP vb, int store_type, ...)       { SET_MULTI_CTX (store_type, flags.store, store_type); } // clang issues a warning if store_type is of type StoreType
void ctx_set_ltype (VBlockP vb, int ltype,      ...)       { SET_MULTI_CTX (ltype, ltype, ltype); }                 // clang issues a warning if ltype is of type LocalType
void ctx_consolidate_stats (VBlockP vb, int parent,   ...) { SET_MULTI_CTX (parent, st_did_i, parent); CTX(parent)->is_stats_parent = true;} // clang issues a warning if parent is of type Did

void ctx_consolidate_stats_(VBlockP vb, Did parent, unsigned num_deps, ContextP dep_ctxs[])
{
    for (unsigned d=0; d < num_deps; d++)
        if (dep_ctxs[d]->did_i != parent) 
            dep_ctxs[d]->st_did_i = parent;

    if (CTX(parent)->st_did_i == DID_NONE)
        CTX(parent)->is_stats_parent = true;
}

// consolidate a consecutive block of Dids
void ctx_consolidate_statsN (VBlockP vb, Did parent, Did first_dep, unsigned num_deps)
{
    for (ContextP ctx=CTX(first_dep); ctx < CTX(first_dep + num_deps); ctx++)
        if (ctx->did_i != parent) 
            ctx->st_did_i = parent;

    if (CTX(parent)->st_did_i == DID_NONE)
        CTX(parent)->is_stats_parent = true;
}

ContextP buf_to_ctx (ContextArray ca, ConstBufferP buf) 
{ 
    if (is_p_in_range (buf, ca, sizeof(ContextArray))) 
        return &ca[((rom)buf - (rom)(ca)) / ((sizeof(ContextArray)) / MAX_DICTS)]; // note "sizeof(ContextArray)) / MAX_DICTS" might be a bit bigger than sizeof(Context) due to alignment
    else
        return NULL;
}
