// ------------------------------------------------------------------
//   context.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include "genozip.h"
#include "profiler.h"
#include "sections.h"
#include "vcf.h"
#include "base250.h"
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

#define INITIAL_NUM_NODES 10000

// inserts dict_id->did_i to the map, if one of the two entries is available
static inline void set_dict_id_to_did_i_map (DidIType *map, DictId dict_id, DidIType did_i)
{
    // thread safety for z_file map: we don't bother with having a mutex, in the worst case scenario, two threads will test an entry
    // as empty and then both write to it, with one of them prevailing. that's fine (+ Likely its the same did_i anyway).

    if (map[dict_id.map_key] == DID_I_NONE)  // map is free
        map[dict_id.map_key] = did_i;

    else if (map[dict_id.map_key] == did_i)  // already has requested value - nothing to do 
        {}
    
    else if (map[ALT_KEY(dict_id)] == DID_I_NONE) // fallback entry is free or we can override it
        map[ALT_KEY(dict_id)] = did_i;
}

// show a dict_id if the requested string is a subset of it, excluding unprintable characters
bool ctx_is_show_dict_id (DictId dict_id)
{
    if (!flag.show_one_dict) return false;
    
    char dict_id_str[9] = "";
    unsigned s_len=0;
    dict_id = dict_id_typeless (dict_id);

    for (unsigned i=0; i < DICT_ID_LEN; i++)
        if (IS_NON_WS_PRINTABLE(dict_id.id[i]))
            dict_id_str[s_len++] = dict_id.id[i];

    unsigned len = strlen (flag.show_one_dict);
    if (len <= 8)
        return !memcmp (dict_id_str, flag.show_one_dict, len);
    else
        return !memcmp (dict_id_str, flag.show_one_dict, 4) && !memcmp(&dict_id_str[4], &flag.show_one_dict[len-4], 4);
}

// ZIP: add a snip to the dictionary the first time it is encountered in the VCF file.
// the dictionary will be written to GENOZIP and used to reconstruct the MTF during decompression
typedef enum { DICT_VB, DICT_ZF, DICT_ZF_SINGLETON } DictType;
static inline CharIndex ctx_insert_to_dict (VBlock *vb_of_dict, Context *ctx, DictType type, const char *snip, uint32_t snip_len)
{
    Buffer *dict = (type == DICT_ZF_SINGLETON) ? &ctx->ol_dict : &ctx->dict; // in z_file, ol_dict contains singletons

    static const char *buf_name[3] = { "contexts->dict", "zctx->dict", "zctx->ol_dict" };
    buf_alloc (vb_of_dict, dict, snip_len + 1, INITIAL_NUM_NODES * MIN_(10, snip_len), char,CTX_GROWTH, buf_name[type]);
    
    if (type == DICT_ZF) buf_set_overlayable (dict); // during merge
    
    CharIndex char_index = dict->len;
    char *dict_p = ENT (char, *dict, char_index);

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = 0; // dictionaries have a \0 separating snips, so that PIZ can generate word_list

    dict->len += snip_len + 1;
    return char_index;
}

// ZIP only (PIZ doesn't have nodes) nodes index to node - possibly in ol_nodes, or in nodes
CtxNode *ctx_node_vb_do (const Context *vctx, WordIndex node_index, 
                         const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                         const char *func, uint32_t code_line)
{
    ASSERT (vctx->dict_id.num, "this vctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index < vctx->nodes.len + vctx->ol_nodes.len, "out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u. Caller: %s:%u",  
            vctx->tag_name, node_index, (uint32_t)vctx->nodes.len, (uint32_t)vctx->ol_nodes.len, func, code_line);

    bool is_ol = node_index < vctx->ol_nodes.len; // is this entry from a previous vb (overlay buffer)

    CtxNode *node = is_ol ? ENT (CtxNode, vctx->ol_nodes, node_index)
                          : ENT (CtxNode, vctx->nodes, node_index - vctx->ol_nodes.len);

    if (snip_in_dict) {
        const Buffer *dict = is_ol ? &vctx->ol_dict : &vctx->dict;
        ASSERT0 (buf_is_alloc (dict), "dict not allocated");

        ASSERT (node->char_index + (uint64_t)node->snip_len < dict->len, "snip of %s out of range: node->char_index=%"PRIu64" + node->snip_len=%u >= %s->len=%"PRIu64,
                vctx->tag_name, node->char_index, node->snip_len, is_ol ? "ol_dict" : "dict", dict->len);

        *snip_in_dict = ENT (char, *dict, node->char_index);
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// ZIP only (PIZ doesn't have nodes) nodes index to node - possibly in ol_nodes, or in nodes
CtxNode *ctx_node_zf_do (const Context *zctx, int32_t node_index, 
                         const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                         const char *func, uint32_t code_line)
{
    ASSERT (zctx->dict_id.num, "this zctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index > -2 - (int32_t)zctx->ol_nodes.len && node_index < (int32_t)zctx->nodes.len, 
            "out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u. Caller: %s:%u",  
            zctx->tag_name, node_index, (uint32_t)zctx->nodes.len, (uint32_t)zctx->ol_nodes.len, func, code_line);

    bool is_singleton = node_index < 0; // is this entry from a previous vb (overlay buffer)

    CtxNode *node = is_singleton ? ENT (CtxNode, zctx->ol_nodes, -node_index - 2)
                                 : ENT (CtxNode, zctx->nodes, node_index);

    if (snip_in_dict) {
        const Buffer *dict = is_singleton ? &zctx->ol_dict : &zctx->dict;
        ASSERT0 (buf_is_alloc (dict), "dict not allocated");

        *snip_in_dict = &dict->data[node->char_index];
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// PIZ: search for a node matching this snip in a directory and return the node index. note that we do a linear
// search as PIZ doesn't have hash tables.
WordIndex ctx_search_for_word_index (Context *ctx, const char *snip, unsigned snip_len)
{
    CtxWord *words = (CtxWord *)ctx->word_list.data;
    ASSERT (words, "word_list is empty for context %s", ctx->tag_name);

    for (unsigned i=0; i < ctx->word_list.len; i++)
        if (words[i].snip_len == snip_len && !memcmp (&ctx->dict.data[words[i].char_index], snip, snip_len))
            return i;

    return WORD_INDEX_NONE;
}

// PIZ and reading Pair1 in ZIP (uses word_list): returns word index, and advances the iterator
WordIndex ctx_get_next_snip (VBlock *vb, Context *ctx, bool all_the_same, bool is_pair, pSTRp (snip)/*optional out*/)
{
    ASSERT (ctx, "ctx is NULL. vb_i=%u", vb->vblock_i);

    const Buffer *b250 = is_pair ? &ctx->pair : &ctx->b250;
    SnipIterator *iterator = is_pair ? &ctx->pair_b250_iter : &ctx->iterator;

    // if the entire b250 in a VB consisted of word_index=0, we don't output the b250 to the file, and just 
    // consider it to always emit 0
    if (!b250->len) {
        if (snip) {
            CtxWord *dict_word = FIRSTENT (CtxWord, ctx->word_list);
            *snip = &ctx->dict.data[dict_word->char_index];
            *snip_len = dict_word->snip_len;
        }
        return 0;
    }
    
    if (!iterator->next_b250)  // initialize
        *iterator = (SnipIterator){ .next_b250       = FIRSTENT (uint8_t, *b250), 
                                    .prev_word_index = WORD_INDEX_NONE };

    WordIndex word_index = base250_decode (&iterator->next_b250, !all_the_same, ctx->tag_name);  // if this line has no non-GT subfields, it will not have a ctx 

    // we check after (no risk of segfault because of buffer overflow protector) - since b250 word is variable length
    if (command == PIZ)
        ASSPIZ (iterator->next_b250 <= AFTERENT (uint8_t, *b250), 
                "while reconstructing (vb->lines.len=%"PRIu64"): iterator for %s (did_i=%u) %sreached end of b250. b250.len=%"PRIu64, 
                vb->lines.len, ctx->tag_name, ctx->did_i, is_pair ? "(PAIR) ": "", b250->len);
    else
        ASSERT (iterator->next_b250 <= AFTERENT (uint8_t, *b250), 
                "while reconstructing vb->line_i=%"PRIu64" (last line of vb is %"PRIu64") in vb_i=%u: iterator for %s (did_i=%u) %sreached end of b250. b250.len=%"PRIu64, 
                vb->line_i, vb->first_line + vb->lines.len - 1, vb->vblock_i, ctx->tag_name, ctx->did_i, is_pair ? "(PAIR) ": "", b250->len);

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

        ASSERT (word_index < ctx->word_list.len, 
                "while parsing vb=%u line=%"PRIu64": word_index=%u is out of bounds - %s %sdictionary (did=%u) has only %u entries. b250.len=%"PRId64" iterator(after)=%"PRId64,
                vb->vblock_i, vb->line_i, word_index, ctx->tag_name, is_pair ? "(PAIR) ": " ", 
                ctx->did_i, (uint32_t)ctx->word_list.len, b250->len, (uint64_t)((char*)iterator->next_b250 - (char*)b250->data));

        CtxWord *dict_word = ENT (CtxWord, ctx->word_list, word_index);

        if (snip) {
            *snip = &ctx->dict.data[dict_word->char_index];
            *snip_len = dict_word->snip_len;
        }
    }
    
    iterator->prev_word_index = word_index;    

    return word_index;
}

// similar to ctx_get_next_snip, except iterator is left intact 
WordIndex ctx_peek_next_snip (VBlock *vb, Context *ctx, bool all_the_same, pSTRp (snip)/*optional out*/)  
{
    SnipIterator save_iterator = ctx->iterator;
    WordIndex wi = ctx_get_next_snip (vb, ctx, all_the_same, false, STRa(snip));
    ctx->iterator = save_iterator;
    return wi;
}

// Process and snip - return its node index, and enter it into the directory if its not already there. Called
// 1. During segregate - as snips are encountered in the data. No base250 encoding yet
// 2. During ctx_merge_in_one_vctx() - to enter snips into z_file->contexts - also encoding in base250
static WordIndex ctx_commit_node (VBlock *vb, Context *zctx, Context *vctx, STRp(snip), int64_t count, CtxNode **node, bool *is_new)  // out
{
    // if this turns out to be a singelton - i.e. a new snip globally - where it goes depends on whether its a singleton in the VB 
    // and whether we are allowed to move singletons to local. if its a singleton:
    // 1. we keep it in ston_nodes and the index in the node is negative to indicate that
    // 2. we insert it to ston_nodes instead of dict - i.e. it doesn't get written the dict section
    // 3. we move it to the local section of this vb
    // 4. we set the word_index of its nodes to be the word_index of the SNIP_LOOKUP snip
    bool is_singleton_in_vb = (count == 1 && ctx_can_have_singletons (vctx)); // is singleton in this VB

    // attempt to get the node from the hash table 
    WordIndex node_index = hash_global_get_entry (zctx, snip, snip_len, 
                                                  is_singleton_in_vb ? HASH_NEW_OK_SINGLETON_IN_VB : HASH_NEW_OK_NOT_SINGLETON, node);

    // case: existing non-singleton node. Possibly it was a singleton node before, that thanks to us hash_global_get_entry
    // converted to non-singleton - i.e. node_index is always >= 0.
    if (*node) { 
        *is_new = false;
        return node_index; // >= 0
    }
    
    // NEW SNIP globally - this snip was just added to the hash table - either as a regular or singleton node
    bool is_singleton_in_global = (node_index < 0);
    Buffer *nodes = is_singleton_in_global ? &zctx->ston_nodes : &zctx->nodes;
    ASSERT (nodes->len <= MAX_WORDS_IN_CTX, 
            "too many words in ctx %s, max allowed number of words is is %u", zctx->tag_name, MAX_WORDS_IN_CTX);

    // set either singleton node or regular node with this snip
    buf_alloc (evb, nodes, 1, INITIAL_NUM_NODES, CtxNode, CTX_GROWTH, is_singleton_in_global ? "zctx->ston_nodes" : "zctx->nodes");
    *node = LASTENT (CtxNode, *nodes); // note: (*nodes).len was already incremented in hash_global_get_entry
    **node = (CtxNode){
        .snip_len   = snip_len,
        .char_index = ctx_insert_to_dict (evb, zctx, (is_singleton_in_global ? DICT_ZF_SINGLETON : DICT_ZF), snip, snip_len)
    };

    // case: vb singleton turns out to be a global singleton - we add it to local and return the SNIP_LOOKUP node
    // instead of the singleton node (which is guaranteed to be non-singleton, and hence >= 0)
    // note: local is dedicated to singletons and contains nothing else, since inst.no_stons is not set
    if (node_index < 0) {
        seg_add_to_local_text (vb, vctx, snip, snip_len, 0);

        if (HAS_DEBUG_SEG(vctx)) {
            char printable_snip[snip_len+20];
            iprintf ("commit_node: vb_i=%u %s: SINGLETON replaced by SNIP_LOCAL: snip=%s snip_len=%u\n", 
                     vb->vblock_i, vctx->tag_name, str_print_snip (snip, snip_len, printable_snip), snip_len);
        }

        static char lookup = SNIP_LOOKUP;
        return ctx_commit_node (vb, zctx, vctx, &lookup, 1, -1 /* not singleton */, node, is_new);
    }

    // case: not singleton - we return this (new) node
    else {
        buf_set_overlayable (&zctx->nodes);
        (*node)->word_index.n = node_index;

        buf_alloc_zero (evb, &zctx->counts, 1, INITIAL_NUM_NODES, int64_t, CTX_GROWTH, "zctx->counts");
        zctx->counts.len++; // actually assigned in ctx_merge_in_one_vctx 

        if (chrom_2ref_seg_is_needed (zctx->did_i)) { 
            buf_alloc_255 (evb, &z_file->chrom2ref_map, 1, INITIAL_NUM_NODES, WordIndex, CTX_GROWTH, "z_file->chrom2ref_map");
            z_file->chrom2ref_map.len++; // actually assigned in ctx_merge_in_one_vctx 
        }

        if (is_new) *is_new = true;

        return node_index; // >= 0
    }
}

// Seg: inserts snip into the hash, nodes and dictionary, if it was not already there, and returns its node index.
// Does NOT add the word index to b250.
WordIndex ctx_create_node_do (VBlockP vb, ContextP vctx, STRp(snip), bool *is_new /* out */)
{
    ASSERTNOTNULL (vctx);
    ASSERT (vctx->dict_id.num, "vctx has no dict_id (did_i=%u)", (unsigned)(vctx - vb->contexts));
    WordIndex node_index;
    #define RETURN(ni) { node_index=(ni); goto done; }

    if (!snip_len) {
        if (is_new) *is_new = false;
        RETURN ((!snip || (vb->data_type == DT_VCF && dict_id_is_vcf_format_sf (vctx->dict_id) && *snip != ':')) 
                 ? WORD_INDEX_MISSING : WORD_INDEX_EMPTY);
    }

    // short-circuit if identical to previous snip
    if (vctx->last_snip && str_issame (snip, vctx->last_snip)) {
        (*ENT (int64_t, vctx->counts, vctx->last_snip_ni))++; 
        if (is_new) *is_new = false;
        RETURN (vctx->last_snip_ni);
    }
    
    WordIndex node_index_if_new = vctx->ol_nodes.len + vctx->nodes.len;
    
#ifdef DEBUG // time consuming and only needed during development
    unsigned actual_len = strnlen (snip, snip_len);
    ASSERT (actual_len == snip_len, "vb=%u ctx=%s: snip_len=%u but unexpectedly has an 0 at index %u: \"%.*s\"", 
            vb->vblock_i, vctx->tag_name, snip_len, actual_len, snip_len, snip);
#endif

    ASSERT (node_index_if_new <= MAX_NODE_INDEX, 
            "ctx of %s is full (max allowed words=%u): ol_nodes.len=%u nodes.len=%u",
            vctx->tag_name, MAX_WORDS_IN_CTX, (uint32_t)vctx->ol_nodes.len, (uint32_t)vctx->nodes.len);

    // get the node from the hash table if it already exists, or add this snip to the hash table if not
    CtxNode *node;
    WordIndex existing_node_index = hash_get_entry_for_seg (vb, vctx, snip, snip_len, node_index_if_new, &node);
    if (existing_node_index != NODE_INDEX_NONE) {
        (*ENT (int64_t, vctx->counts, existing_node_index))++; // note: counts.len = nodes.len + ol_nodes.len
        if (is_new) *is_new = false;
        
        // populate last_snip 
        bool is_ol = existing_node_index < vctx->ol_nodes.len; // is this entry from a previous vb (overlay buffer)
        const Buffer *dict = is_ol ? &vctx->ol_dict : &vctx->dict;
        vctx->last_snip     = ENT (char, *dict, node->char_index);
        vctx->last_snip_len = snip_len;
        vctx->last_snip_ni  = existing_node_index;

        RETURN (existing_node_index); // snip found - we're done
    }
    
    // this snip isn't in the hash table - its a new snip
    ASSERT (vctx->nodes.len < MAX_NODE_INDEX, "too many words in dictionary %s (MAX_NODE_INDEX=%u)", vctx->tag_name, MAX_NODE_INDEX);

    buf_alloc (vb, &vctx->nodes,  1, INITIAL_NUM_NODES, CtxNode, CTX_GROWTH, "contexts->nodes");
    buf_alloc (vb, &vctx->counts, 1, INITIAL_NUM_NODES, int64_t, CTX_GROWTH, "contexts->counts");

    NEXTENT (CtxNode, vctx->nodes) = (CtxNode){
        .snip_len     = snip_len,
        .char_index   = ctx_insert_to_dict (vb, vctx, DICT_VB, snip, snip_len),
        .word_index.n = node_index_if_new
    };

    // populate last_snip / last_snip_len if requested
    vctx->last_snip     = ENT (char, vctx->dict, LASTENT (CtxNode, vctx->nodes)->char_index);
    vctx->last_snip_len = snip_len;
    vctx->last_snip_ni  = node_index_if_new;

    NEXTENT (int64_t, vctx->counts) = 1;

    ASSERT (vctx->counts.len == vctx->nodes.len + vctx->ol_nodes.len, "Expecting vctx->counts.len=%"PRId64" == vctx->nodes.len=%"PRId64" + vctx->ol_nodes.len=%"PRId64,
            vctx->counts.len, vctx->nodes.len, vctx->ol_nodes.len);

    if (is_new) *is_new = true;
    RETURN (node_index_if_new);

done:
    if (HAS_DEBUG_SEG(vctx)) { 
        char printable_snip[snip_len+20];
        if (snip) iprintf ("create_node: vb_i=%u %s(%u): snip=%s snip_len=%u node_index=%d\n", vb->vblock_i, vctx->tag_name, vctx->did_i, str_print_snip (snip, snip_len, printable_snip), snip_len, node_index);
        else      iprintf ("create_node: vb_i=%u %s(%u): snip=NULL snip_len=0 node_index=%d\n", vb->vblock_i, vctx->tag_name, vctx->did_i, node_index);
    }

    return node_index;
    #undef RETURN
}

// Seg only: create a node without adding to counts
WordIndex ctx_create_node (VBlockP vb, DidIType did_i, STRp (snip))
{
    WordIndex node_index = ctx_create_node_do (vb, CTX(did_i), STRa(snip), NULL); 
    ctx_decrement_count (vb, CTX(did_i), node_index);

    return node_index;
}

// Seg only: if after ctx_create_node_do we don't add the snip to b250, we need to reduce its count
int64_t ctx_decrement_count (VBlock *vb, Context *ctx, WordIndex node_index)
{
    ASSERT (node_index < (WordIndex)ctx->counts.len, "node_index=%d out of range counts[%s].len=%"PRIu64, node_index, ctx->tag_name, ctx->counts.len);

    int64_t *count_p = ENT (int64_t, ctx->counts, node_index);

    if (node_index < 0) return 0; // WORD_INDEX_EMPTY or WORD_INDEX_MISSING
     
    ASSERT (*count_p >= 1, "count[%s]=%"PRId64" too low to be decremented", ctx->tag_name, *count_p);
    (*count_p)--;

    // edge case: if we removed the last b250 of this node_index, we now have an unused node and an unused word dict. 
    // this may in some cases cause the condition "ctx->nodes.len != ctx->b250.len" in zip_handle_unique_words_ctxs to
    // incorrectly fail, causing moving of an incorrect dict to local. to prevent this, we don't allow singletons in this case.
    if (! *count_p && node_index >= ctx->ol_nodes.len) ctx->no_stons = true;

    return *count_p;
}

// Seg only: if we add a b250 without evaluating (if node_index is known)
void ctx_increment_count (VBlock *vb, Context *ctx, WordIndex node_index)
{
    ASSERT (node_index < ctx->counts.len, "node_index=%d out of range counts[%s].len=%"PRIu64, node_index, ctx->tag_name, ctx->counts.len);

    (*ENT (int64_t, ctx->counts, node_index))++;
}

void ctx_append_b250 (VBlockP vb, ContextP vctx, WordIndex node_index)
{
    // case: context is currently all_the_same - (ie len>=1, but only one actual entry)...
    if (vctx->flags.all_the_same) {
        WordIndex first_node_index; 

        // ... and the new entry is more of the same
        if (node_index == (first_node_index = *FIRSTENT (WordIndex,vctx->b250)))
            vctx->b250.len++; // increment len, but don't actually grow the array

        // ... and the new entry is different - add (len-1) copies of first_node, and the new one
        else {
            buf_alloc (vb, &vctx->b250, 1, AT_LEAST(vctx->did_i), WordIndex, CTX_GROWTH, "contexts->b250"); // add 1 more, meaning a total of len+1
        
            for (unsigned i=1; i < vctx->b250.len; i++) *ENT (WordIndex, vctx->b250, i) = first_node_index;
            NEXTENT (uint32_t, vctx->b250) = node_index;

            vctx->flags.all_the_same = false; // no longer all_the_same 
        }
    }

    // b250 is not all_the_same: either it is empty, or it contains multiple different node_index values
    else {
        // case: this is the first entry - its all_the_same in a trivial way (still, we mark it so zip_generate_b250_section can potentially eliminate it)
        if (!vctx->b250.len) vctx->flags.all_the_same = true;

        buf_alloc (vb, &vctx->b250, 1, AT_LEAST(vctx->did_i), WordIndex, CTX_GROWTH, "contexts->b250");
        NEXTENT (uint32_t, vctx->b250) = node_index;
    }
}

// ZIP only: overlay and/or copy the current state of the global contexts to the vb, ahead of segging this vb.
void ctx_clone (VBlock *vb)
{
    unsigned z_num_contexts = __atomic_load_n (&z_file->num_contexts, __ATOMIC_RELAXED);

    START_TIMER; // including mutex wait time

    // note: because each dictionary has its own mutex, it is possible that we will see only a partial set
    // of dictionaries (eg some but not all of the fields) when we are arrive here while another thread is mid-way 
    // through merging and adding a bunch of dictionaries.
    // however z_num_contexts will always correctly state the number of dictionaries that are available.

    for (DidIType did_i=0; did_i < z_num_contexts; did_i++) {
        Context *vctx = CTX(did_i);
        Context *zctx = ZCTX(did_i);

        // case: this context doesn't really exist (happens when incrementing num_contexts when adding RNAME and RNEXT in ctx_populate_zf_ctx_from_contigs)
        if (!zctx->mutex.initialized) continue;

        mutex_lock (zctx->mutex);

        if (buf_is_alloc (&zctx->dict)) {  // something already for this dict_id

            // overlay the global dict and nodes - these will not change by this (or any other) VB
            //iprintf ( ("ctx_clone: overlaying old dict %.8s, to vb_i=%u vb_did_i=z_did_i=%u\n", dis_dict_id (zctx->dict_id).s, vb->vblock_i, did_i);
            buf_overlay (vb, &vctx->ol_dict,  &zctx->dict,  "contexts->ol_dict");   
            buf_overlay (vb, &vctx->ol_nodes, &zctx->nodes, "contexts->ol_nodes");   

            // overlay the hash table, that may still change by future vb's merging... this vb will only use
            // entries that are up to this merge_num
            buf_overlay (vb, &vctx->global_hash, &zctx->global_hash, "contexts->global_hash");
            vctx->merge_num = zctx->merge_num;
            vctx->global_hash_prime = zctx->global_hash_prime; // can never change
            vctx->num_new_entries_prev_merged_vb = zctx->num_new_entries_prev_merged_vb;

            buf_alloc_zero (vb, &vctx->counts, 0, vctx->ol_nodes.len, int64_t, CTX_GROWTH, "contexts->counts");
            vctx->counts.len = vctx->ol_nodes.len;

        }

        vctx->did_i       = did_i;
        vctx->dict_id     = zctx->dict_id;
        vctx->st_did_i    = zctx->st_did_i;
        vctx->luft_trans  = zctx->luft_trans;
        vctx->last_line_i = LAST_LINE_I_INIT;
        vctx->tag_i       = -1;

        // note: lcodec and bcodec are inherited in merge (see comment in zip_assign_best_codec)

        memcpy ((char*)vctx->tag_name, zctx->tag_name, sizeof (vctx->tag_name));

        set_dict_id_to_did_i_map (vb->dict_id_to_did_i_map, vctx->dict_id, did_i);
        
        ctx_init_iterator (vctx);

        // note: chrom2ref_map is protected by ZCTX[CHROM]->mutex
        if (chrom_2ref_seg_is_needed(did_i) && z_file->chrom2ref_map.len) {
            buf_set_overlayable (&z_file->chrom2ref_map);
            buf_overlay (vb, &vb->ol_chrom2ref_map, &z_file->chrom2ref_map, "ol_chrom2ref_map");
        }

        mutex_unlock (zctx->mutex);
    }

    vb->num_contexts = z_num_contexts;
    
    COPY_TIMER (ctx_clone);
}

static void ctx_initialize_ctx (Context *ctx, DidIType did_i, DictId dict_id, DidIType *dict_id_to_did_i_map,
                                const char *tag_name, unsigned tag_name_len)
{
    ctx->did_i       = did_i;
    ctx->st_did_i    = DID_I_NONE;
    ctx->dict_id     = dict_id;
    ctx->last_line_i = LAST_LINE_I_INIT;
    ctx->tag_i       = -1;
    ctx->dict.can_be_big = true; // don't warn if dict buffers grow really big

    if (tag_name_len) {
        ASSINP (tag_name_len <= MAX_TAG_LEN-1, "Tag name \"%.*s\" is of length=%u beyond the maximum tag length supported by Genozip=%u",
                tag_name_len, tag_name, tag_name_len, MAX_TAG_LEN-1);

        memcpy ((char*)ctx->tag_name, tag_name, tag_name_len);
    }
    else
        strcpy ((char*)ctx->tag_name, dis_dict_id (dict_id).s);

    ctx_init_iterator (ctx);
    
    set_dict_id_to_did_i_map (dict_id_to_did_i_map, dict_id, did_i);

    bool is_zf_ctx = z_file && (ctx - z_file->contexts) >= 0 && (ctx - z_file->contexts) <= (sizeof(z_file->contexts)/sizeof(z_file->contexts[0]));

    // add a user-requested SEC_COUNT section
    if (command == ZIP) {
        if (flag.show_one_counts.num == dict_id_typeless (ctx->dict_id).num) 
            ctx->counts_section = ctx->no_stons = true;

        if (is_zf_ctx) mutex_initialize (ctx->mutex);

        // a new non-field context that is not defined in the header - see if we have a default translator for it
        if (z_dual_coords && !dict_id_is_field (dict_id)) 
            ctx->luft_trans = vcf_lo_luft_trans_id (dict_id, '.'); // TO DO: make data-type agnostic (bug 359)
    } 
}

WordIndex ctx_populate_zf_ctx (DidIType dst_did_i, STRp (contig_name), WordIndex ref_index)
{
    Context *zctx = ZCTX(dst_did_i);
    
    if (!buf_is_alloc (&zctx->global_hash)) { // first call
        zctx->no_stons = true;
        zctx->st_did_i = DID_I_NONE;    
        hash_alloc_global (zctx, 0);
    }

    if (flag.reference & REF_ZIP_LOADED)
        buf_alloc_255 (evb, &z_file->chrom2ref_map, 1, INITIAL_NUM_NODES, WordIndex, 1, "z_file->chrom2ref_map");

    bool is_primary = (dst_did_i == DTFZ (prim_chrom));
    
    CtxNode *zf_node;
    WordIndex zf_node_index = ctx_commit_node (NULL, zctx, NULL, STRa(contig_name), 0, &zf_node, NULL);
    zf_node->word_index = base250_encode (zf_node_index);

    if (is_primary && ref_index != WORD_INDEX_NONE && (flag.reference & REF_ZIP_LOADED))
        *ENT (WordIndex, z_file->chrom2ref_map, zf_node_index) = ref_index;

    return zf_node_index;
}

// ZIP main thread: 
// 1. when starting to zip a new file, with pre-loaded external reference, we integrate the reference's FASTA CONTIG
//    dictionary as the chrom dictionary of the new file
// 2. in SAM, DENOVO: after creating contigs from SQ records, we copy them to the RNAME dictionary
// 3. When loading a chain file - copy tName word_list/dict read from the chain file to a context
void ctx_populate_zf_ctx_from_contigs (Reference ref, DidIType dst_did_i, ConstContigPkgP ctgs)
{
    if (!ctgs->contigs.len) return; // nothing to do

    Context *zctx = ZCTX(dst_did_i);
    
    if (!buf_is_alloc (&zctx->global_hash)) { // first call
        zctx->no_stons = true;
        zctx->st_did_i = DID_I_NONE;
        hash_alloc_global (zctx, 0);
    }

    ARRAY (Contig, contigs, ctgs->contigs);

    if (flag.reference & REF_ZIP_LOADED)
        buf_alloc_255 (evb, &z_file->chrom2ref_map, contigs_len, INITIAL_NUM_NODES, WordIndex, 1, "z_file->chrom2ref_map");

    bool is_primary = (dst_did_i == DTFZ (prim_chrom));

    for (uint64_t i=0; i < contigs_len; i++) {
        CtxNode *zf_node;
        bool is_new;

        STR(contig_name);
        contig_name = contigs_get_name (ctgs, i, &contig_name_len);

        WordIndex zf_node_index = ctx_commit_node (NULL, zctx, NULL, STRa(contig_name), 0, &zf_node, &is_new);
        zf_node->word_index = base250_encode (zf_node_index);

        if (is_new && is_primary && (flag.reference & REF_ZIP_LOADED))
            *ENT (WordIndex, z_file->chrom2ref_map, zf_node_index) = (contigs[i].ref_index != WORD_INDEX_NONE) 
                ? contigs[i].ref_index // header contigs also know the ref_index
                : ref_contigs_get_matching (ref, 0, STRa(contig_name), NULL, NULL, false, NULL, NULL);
    }
}

// find the z_file context that corresponds to dict_id. It could be possibly a different did_i
// than in the vb - in case this dict_id is new to this vb, but another vb already inserted
// it to z_file
static Context *ctx_get_zctx (DictId dict_id, bool check_predefined)
{
    DidIType z_num_contexts = __atomic_load_n (&z_file->num_contexts, __ATOMIC_RELAXED);

    for (ContextP zctx=ZCTX(check_predefined ? 0 : DTFZ(num_fields)); zctx < ZCTX(z_num_contexts); zctx++)
        if (dict_id.num == zctx->dict_id.num) 
            return zctx;

    return NULL;
}

ContextP ctx_get_zctx_from_vctx (ConstContextP vctx)  // returns NULL if context doesn't exist
{
    if (vctx->did_i < DTFZ(num_fields)) 
        return ZCTX(vctx->did_i); 
    
    // check dict_id->did_i map, and if not found search linearly
    DidIType did_i = ctx_get_existing_did_i_do (vctx->dict_id, z_file->contexts, z_file->dict_id_to_did_i_map, NULL, z_file->num_contexts);
    return did_i != DID_I_NONE ? ZCTX(did_i) : NULL;
}

// ZIP only: called by merging VBs to add a new dict to z_file - copying some stuff from vctx
static Context *ctx_add_new_zf_ctx (VBlock *vb, const Context *vctx)
{
    // adding a new dictionary is proctected by a mutex. note that z_file->num_contexts is accessed by other threads
    // without mutex proction when searching for a dictionary - that's why we update it at the end, after the new
    // zctx is set up with the new dict_id (ready for another thread to search it)
    mutex_lock (z_file->dicts_mutex);

    // check if another thread raced and created this dict before us
    Context *zctx = ctx_get_zctx_from_vctx (vctx);
    if (zctx) goto finish;

    ASSERT (z_file->num_contexts+1 < MAX_DICTS, // load num_contexts - this time with mutex protection - it could have changed
            "z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    zctx = ZCTX(z_file->num_contexts);

    mutex_initialize (zctx->mutex);

    zctx->did_i           = z_file->num_contexts; 
    zctx->st_did_i        = vctx->st_did_i;
    zctx->is_stats_parent = vctx->is_stats_parent;
    zctx->dict_id         = vctx->dict_id;
    zctx->luft_trans      = vctx->luft_trans;
    memcpy ((char*)zctx->tag_name, vctx->tag_name, sizeof(zctx->tag_name));
    // note: lcodec is NOT copied here, see comment in zip_assign_best_codec

    // only when the new entry is finalized, do we increment num_contexts, atmoically , this is because
    // other threads might access it without a mutex when searching for a dict_id
    __atomic_store_n (&z_file->num_contexts, z_file->num_contexts+1, __ATOMIC_RELAXED); 

    // only after updating num_contexts, we add it to the map. 
    set_dict_id_to_did_i_map (z_file->dict_id_to_did_i_map, zctx->dict_id, zctx->did_i);

finish:
    mutex_unlock (z_file->dicts_mutex);
    return zctx;
}

// ZIP only: called by main thread when inspecting a txtheader for assigning liftover translators
ContextP ctx_add_new_zf_ctx_from_txtheader (const char *tag_name, unsigned tag_name_len, DictId dict_id, TranslatorId luft_translator)
{
    ASSINP (tag_name_len <= MAX_TAG_LEN-1, "Tag name \"%.*s\" is of length=%u beyond the maximum tag length supported by Genozip=%u",
            tag_name_len, tag_name, tag_name_len, MAX_TAG_LEN-1);

    // adding a new dictionary is proctected by a mutex. note that z_file->num_contexts is accessed by other threads
    // without mutex proction when searching for a dictionary - that's why we update it at the end, after the new
    // zctx is set up with the new dict_id (ready for another thread to search it)
    mutex_lock (z_file->dicts_mutex); // note: mutex needed bc VBs of a previous component might still be merging

    // possibly a tag is duplicate in the header, OR two tag names map to the same dict_id
    Context *zctx = ctx_get_zctx (dict_id, true);
    if (zctx) {
        zctx = NULL; // not new
        goto finish;
    }

    ASSERT (z_file->num_contexts+1 < MAX_DICTS, // load num_contexts - this time with mutex protection - it could have changed
            "z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    zctx = ZCTX(z_file->num_contexts);

    mutex_initialize (zctx->mutex);

    zctx->did_i      = z_file->num_contexts; 
    zctx->st_did_i   = DID_I_NONE;
    zctx->dict_id    = dict_id;
    zctx->luft_trans = luft_translator;
    memcpy ((char*)zctx->tag_name, tag_name, tag_name_len);

    // only when the new entry is finalized, do we increment num_contexts, atmoically , this is because
    // other threads might access it without a mutex when searching for a dict_id
    __atomic_store_n (&z_file->num_contexts, z_file->num_contexts+1, __ATOMIC_RELAXED); 

finish:
    mutex_unlock (z_file->dicts_mutex);

    return zctx; // NULL if not new
}

// update zctx with codec as it is assigned - don't wait for merge, to increase the chance that subsequent
// VBs can get this codec and don't need to test for themselves.
void ctx_commit_codec_to_zf_ctx (VBlock *vb, Context *vctx, bool is_lcodec)
{
    Context *zctx  = ctx_get_zctx_from_vctx (vctx);
    if (!zctx) zctx = ctx_add_new_zf_ctx (vb, vctx); // this context doesn't exist yet in z_file, because no VB with it has merged yet. create it now.

    // note: no need to mutex_lock zctx - these are only used for cloning to vctx, so the then-current one wins
    if (is_lcodec) zctx->lcodec = vctx->lcodec;
    else           zctx->bcodec = vctx->bcodec;
}

static inline void ctx_drop_all_the_same (VBlock *vb, Context *zctx, Context *vctx)
{
    const char *reason=0;
    #define NO_DROP(s) { reason=(s); goto no_drop; }

    // note: pair_b250/local will causes zfile_compress_b250/local_data to set flags.paired, which will be different than vb_i=1 flags
    if (!vctx->flags.all_the_same) return;
    
    if (vctx->no_drop_b250) NO_DROP("no_drop_b250 is set");
    if (vctx->pair_b250)    NO_DROP("pair_b250 is set");
    if (vctx->pair_local)   NO_DROP("pair_local is set");

    ASSERTISALLOCED (vctx->b250);
    WordIndex node_index = *FIRSTENT (WordIndex, vctx->b250); // the only b250 in this context, as it is all_the_same
    if (node_index < 0) NO_DROP ("node_index=WORD_INDEX_*"); // a special WORD_INDEX_* - not droppable

    bool is_new_word = node_index >= vctx->ol_nodes.len; // its a new word if its not in the overlayed nodes (from previous VBs)
    
    if (!is_new_word && ENT (CtxNode, vctx->ol_nodes, node_index)->word_index.n > 0) NO_DROP ("existing word_index is not 0"); // old word - but not word_index=0

    if (is_new_word && 
        zctx->dict.len && strcmp (zctx->dict.data, vctx->dict.data)) // if it is_new_word - b250 can still be droppable if this new word (the only word in vctx->dict) is the same as word_index_i=0 in zctx->dict, created independently by two VBs in parallel 
        NO_DROP ("new word_index is not 0"); // new word - but not the first in the dictionary so not word_index=0

    STR(snip);
    ctx_get_vb_snip_ex (vctx, node_index, pSTRa(snip)); // get the the snip of the only b250 we have

    if (*snip == SNIP_SELF_DELTA) NO_DROP ("snip is SNIP_SELF_DELTA");

    bool is_simple_lookup = (*snip == SNIP_LOOKUP) && (snip_len == 1);

    // we can't drop b250 if we have local data (as reconstructor will reconstruct from local rather than dict), unless it is a simple SNIP_LOOKUP
    if (vctx->local.len && !is_simple_lookup) NO_DROP ("ctx has local, but snip is not SNIP_LOOKUP"); 

    struct FlagsCtx my_flags = vctx->flags;
    struct FlagsCtx vb_1_flags = (vb->vblock_i==1) ? (struct FlagsCtx){} : zctx->flags; // flags are set by vb_i=1 
    my_flags.all_the_same = vb_1_flags.all_the_same = false; // compare flags, except for all_the_same

    // we can't drop if vctx has flags that need to be passed to piz (unless its vb_i>1 and flags are same as vb_i=1)
    if (((SectionFlags)my_flags).flags != ((SectionFlags)vb_1_flags).flags) NO_DROP (vb->vblock_i==1 ? "this is vb_i=1, and vctx has flags" : "vctx has flags, different from those set by vb_i=1");

    // we survived all tests - we can now drop the b250 and maybe also dict

    if (flag.debug_generate) 
        iprintf ("%s vb_i=%u is \"all_the_same\" - dropped %sb250 b250.len=%"PRIu64"\n", vctx->tag_name, vb->vblock_i, (is_simple_lookup ? "dict AND " : ""), vctx->b250.len);
    
    buf_free (&vctx->b250);
    if (is_simple_lookup) buf_free (&vctx->dict); // if this is a SNIP_LOOKUP, we don't need the new dict entry added in this VB, as absent dict, lookup will happen 
    return;

no_drop:
    if (flag.debug_generate) 
        iprintf ("%s vb_i=%u is all_the_same but cannot drop b250 because %s\n", vctx->tag_name, vb->vblock_i, reason);
}

// ZIP only: this is called towards the end of compressing one vb - merging its dictionaries into the z_file 
// each dictionary is protected by its own mutex, and there is one z_file mutex protecting num_dicts.
// we are careful never to hold two muteces at the same time to avoid deadlocks
static bool ctx_merge_in_one_vctx (VBlock *vb, ContextP vctx)
{
    // get the ctx or create a new one. note: ctx_add_new_zf_ctx() must be called before mutex_lock() because it locks the z_file mutex (avoid a deadlock)
    Context *zctx  = ctx_get_zctx_from_vctx (vctx);
    zctx = ctx_add_new_zf_ctx (vb, vctx);

    if (!mutex_trylock (zctx->mutex)) 
        return false; // zctx is currently locked by another VB 

    //iprintf ( ("Merging dict_id=%.8s into z_file vb_i=%u vb_did_i=%u z_did_i=%u\n", dis_dict_id (vctx->dict_id).s, vb->vblock_i, did_i, z_did_i);

    ctx_drop_all_the_same (vb, zctx, vctx); // drop b250 and maybe also dict, if warranted
    
    zctx->merge_num++; // first merge is #1 (first clone which happens before the first merge, will get vb-)
    zctx->num_new_entries_prev_merged_vb = vctx->nodes.len; // number of new words in this dict from this VB
    zctx->num_singletons += vctx->num_singletons; // add singletons created by seg_by_ctx_ex (SNIP_LOOKUP b250 and snip in local)
    zctx->counts_section |= vctx->counts_section; // for use of ctx_compress_counts
    
    if (vb->vblock_i == 1)
        zctx->flags = vctx->flags; // vb_1 flags will be the default flags for this context, used by piz in case there are no b250 or local sections due to all_the_same. see zip_generate_b250_section and piz_read_all_ctxs

    uint64_t ol_len = vctx->ol_nodes.len;
    bool has_count = !vb->is_rejects_vb; // don't count rejects VB - these are duplicate lines counted in the normal VBs.

    if (flag.rejects_coord != DC_PRIMARY) // we don't include ##primary_only VBs as they are not in the primary reconstruction, but we do include ##luft_only
        zctx->txt_len += vctx->txt_len; // for stats

    if (vctx->st_did_i != DID_I_NONE && zctx->st_did_i == DID_I_NONE) {
        Context *st_ctx = ctx_get_zctx_from_vctx (CTX(vctx->st_did_i));
        if (st_ctx) zctx->st_did_i = st_ctx->did_i; // st_did_i is not necessarily the same for vb and zf
    }

    if (vctx->is_stats_parent)
        zctx->is_stats_parent = true; // we set, but we never revert back

    // we assign VB a codec from zctx, if not already assigned by Seg. See comment in zip_assign_best_codec
    if (!vctx->lcodec) vctx->lcodec = zctx->lcodec;
    if (!vctx->bcodec) vctx->bcodec = zctx->bcodec;
    
    if (!buf_is_alloc (&vctx->dict)) goto finish; // no new snips introduced in this VB
 
    if (!buf_is_alloc (&zctx->dict)) {
        // allocate hash table, based on the statistics gather by this first vb that is merging this dict and 
        // populate the hash table without needing to reevalate the snips (we know none are in the hash table, but all are in nodes and dict)
        if (zctx->global_hash.size <= 1) { // only initial allocation in zip_dict_data_initialize
            uint32_t estimated_entries = hash_get_estimated_entries (vb, zctx, vctx);
            hash_alloc_global (zctx, estimated_entries);
        }
    }

    // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
    // (vctx->nodes contains only new words, old words from previous vbs are in vctx->ol_nodes)
    for (uint64_t i=0; i < vctx->nodes.len; i++) {
        CtxNode *vb_node = ENT (CtxNode, vctx->nodes, i), *zf_node;

        if (vb_node->word_index.n == WORD_INDEX_NONE) continue; // canceled in ctx_rollback

        const char *snip = ENT (char, vctx->dict, vb_node->char_index);
        int64_t count = *ENT (int64_t, vctx->counts, ol_len + i);
        bool is_new;

        // use evb and not vb because zf_context is z_file (which belongs to evb)
        WordIndex zf_node_index = 
            ctx_commit_node (vb, zctx, vctx, snip, vb_node->snip_len, count, &zf_node, &is_new);

        ASSERT (zf_node_index >= 0 && zf_node_index < zctx->nodes.len, 
                "zf_node_index=%d out of range - len=%i", zf_node_index, (uint32_t)vctx->nodes.len);

        if (has_count) 
            *ENT (int64_t, zctx->counts, zf_node_index) += *ENT (int64_t, vctx->counts, ol_len + i); // also allocs nodes, counts and chrom2ref_map

        // set word_index to be indexing the global dict - to be used by vcf_zip_generate_genotype_one_section() and zip_generate_b250_section()
        if (is_new) {
            vb_node->word_index = zf_node->word_index = base250_encode (zf_node_index);

            // note: chrom2ref_map is protected by ZCTX[CHROM]->mutex
            if (chrom_2ref_seg_is_needed(vctx->did_i)) 
                *ENT (WordIndex, z_file->chrom2ref_map, zf_node_index) = *ENT (WordIndex, vb->chrom2ref_map, i);

        } else 
            // a previous VB already already calculated the word index for this node. if it was done by vb_i=1,
            // then it is also re-sorted and the word_index is no longer the same as the node_index
            vb_node->word_index = zf_node->word_index;
    }

finish:
    // just update counts for ol_node (i.e. known to be existing) snips
    if (has_count) 
        for (uint64_t i=0; i < ol_len; i++) 
            *ENT (int64_t, zctx->counts, i) += *ENT (int64_t, vctx->counts, i);

    mutex_unlock (zctx->mutex);
    return true; // merge was done
}

// ZIP only: merge new words added in this vb into the z_file.contexts, and compresses dictionaries.
void ctx_merge_in_vb_ctx (VBlock *vb)
{
    START_TIMER;
    
    // merge all contexts 
    bool all_merged=false;
    while (!all_merged) {

        all_merged=true; // unless proven otherwise
        bool any_merged=false;
        for (ContextP vctx=CTX(0); vctx < CTX(vb->num_contexts); vctx++) 
            if (!vctx->dict_merged) {
                vctx->dict_merged = ctx_merge_in_one_vctx (vb, vctx); // false if zctx is locked by another VB
                all_merged &= vctx->dict_merged;
                any_merged |= vctx->dict_merged;
            }

        // case: couldn't merge any - perhaps the one large final vctx is busy - sleep a bit rather than busy-wait
        if (!any_merged && !all_merged) usleep (100); 
    }

    // note: z_file->num_contexts might be larger than vb->num_contexts at this point, for example:
    // vb_i=1 started, z_file is empty, created 20 contexts
    // vb_i=2 started, z_file is empty, created 10 contexts
    // vb_i=1 completes, merges 20 contexts to z_file, which has 20 contexts after
    // vb_i=2 completes, merges 10 contexts, of which 5 (for example) are shared with vb_i=1. Now z_file has 25 contexts after.

    COPY_TIMER_VB (vb, ctx_merge_in_vb_ctx);
}

// ZIP / PIZ : called from main thread
void ctx_add_compressor_time_to_zf_ctx (VBlockP vb)
{
    for (DidIType vb_did_i=0; vb_did_i < vb->num_contexts; vb_did_i++) {
        ContextP vctx = CTX(vb_did_i); 
        if (vctx->compressor_time) {
            ContextP zctx = vctx->st_did_i != DID_I_NONE ? ctx_get_zctx_from_vctx (CTX(vctx->st_did_i)) // we accumulate at stats parent context if there is one
                                                         : ctx_get_zctx_from_vctx (vctx);
            zctx->compressor_time += vctx->compressor_time;
        }
    }
}

static DidIType ctx_did_i_search (const ContextIndex *ctx_index, DidIType num_contexts, DictId dict_id, DidIType first, DidIType last)
{
    if (last < first) return DID_I_NONE;

    DidIType mid = (first + last) / 2;
    DictId dict_id_mid = ctx_index[mid].dict_id;

    if (dict_id_mid.num == dict_id.num)    
        return ctx_index[mid].did_i;

    else if (dict_id_mid.num < dict_id.num) 
        return ctx_did_i_search (ctx_index, num_contexts, dict_id, mid+1, last);

    else
        return ctx_did_i_search (ctx_index, num_contexts, dict_id, first, mid-1);
}

// returns an existing did_i in this vb, or DID_I_NONE if there isn't one
DidIType ctx_get_unmapped_existing_did_i (const Context *contexts, const ContextIndex *ctx_index, DidIType num_contexts, DictId dict_id)
{
    int did_i; // signed

    // binary search if we have ctx_index (we will have it in PIZ compute threads)
    if (ctx_index) {
        if ((did_i = ctx_did_i_search (ctx_index, num_contexts, dict_id, 0, num_contexts-1)) != DID_I_NONE)
            return did_i;
    }

    else // linear search if not
        for (did_i=num_contexts-1; did_i >= 0 ; did_i--)  // Search backwards as unmapped ctxs are more likely to be towards the end.
            if (dict_id.num == contexts[did_i].dict_id.num) return did_i;

    // PIZ only: check if its an alias that's not mapped in ctx_initialize_predefined_ctxs (due to contention)
    if (command != ZIP && dict_id_aliases) {
        for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
            if (dict_id.num == dict_id_aliases[alias_i].alias.num) { // yes! its an alias

                if (ctx_index)
                    return ctx_did_i_search (ctx_index, num_contexts, dict_id_aliases[alias_i].dst, 0, num_contexts-1);

                else 
                    for (did_i=0; did_i < num_contexts; did_i++) 
                        if (dict_id_aliases[alias_i].dst.num == contexts[did_i].dict_id.num) return did_i;
            }
    }

    return DID_I_NONE; // not found
}

// gets did_id if the dictionary exists, and creates a new dictionary if its the first time dict_id is encountered
// threads: no issues - called by PIZ for vb and zf (but dictionaries are immutable) and by Seg (ZIP) on vctx only
ContextP ctx_get_unmapped_ctx (Context *contexts/* array */, DataType dt, DidIType *dict_id_to_did_i_map, 
                               DidIType *num_contexts, DictId dict_id, STRp(tag_name))
{
    // search to see if this dict_id has a context, despite not in the map (due to contention). 
    for (int/*signed*/ did_i=*num_contexts-1; did_i >= 0 ; did_i--)  // Search backwards as unmapped ctxs are more likely to be towards the end.
        if (dict_id.num == contexts[did_i].dict_id.num) 
            return &contexts[did_i];

    Context *ctx = &contexts[*num_contexts]; 
 
    //iprintf ("New context: dict_id=%s in did_i=%u \n", dis_dict_id (dict_id).s, did_i);
    ASSERT (*num_contexts < MAX_DICTS, "cannot create a context for %.*s because number of dictionaries would exceed MAX_DICTS=%u", 
            tag_name_len, tag_name, MAX_DICTS);

    ctx_initialize_ctx (ctx, *num_contexts, dict_id, dict_id_to_did_i_map, STRa (tag_name));

    // thread safety: the increment below MUST be AFTER the initialization of ctx, bc piz_get_line_subfields
    // might be reading this data at the same time as the piz dispatcher thread adding more dictionaries
    (*num_contexts)++;

    return ctx;
}

// called from seg_all_data_lines (ZIP) and ctx_read_all_dictionaries (PIZ) to initialize all
// primary field ctx's. these are not always used (e.g. when some are not read from disk due to genocat options)
// but we maintain their fixed positions anyway as the code relies on it
void ctx_initialize_predefined_ctxs (Context *contexts /* an array */, 
                                     DataType dt,
                                     DidIType *dict_id_to_did_i_map,
                                     DidIType *num_contexts)
{
    *num_contexts = MAX_(dt_fields[dt].num_fields, *num_contexts);

    for (int did_i=0; did_i < dt_fields[dt].num_fields; did_i++) {
        DictId dict_id = dt_fields[dt].predefined[did_i].dict_id;
        ASSERT (dict_id.num, "No did_i->dict_id mapping is defined for predefined did_i=%u in dt=%s", did_i, dt_name (dt));

        // check if its an alias (PIZ only)
        Context *dst_ctx = NULL;
        if (command != ZIP && dict_id_aliases) 
            for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
                if (dict_id.num == dict_id_aliases[alias_i].alias.num) 
                    // note: all alias destinations that ever existed previous versions of Genozip must be defined in #pragma GENDICT for this to work
                    dst_ctx = ctx_get_zctx (dict_id_aliases[alias_i].dst, true); 

        if (!dst_ctx) // normal field, not an alias
            ctx_initialize_ctx (&contexts[did_i], did_i, dict_id, dict_id_to_did_i_map, 
                                dt_fields[dt].predefined[did_i].tag_name, dt_fields[dt].predefined[did_i].tag_name_len);

        else { // an alias
            contexts[did_i].did_i = dst_ctx->did_i;
            contexts[did_i].dict_id = DICT_ID_NONE; // this is how reconstruct_from_ctx_do identifies it is an alias
            
            set_dict_id_to_did_i_map (dict_id_to_did_i_map, dict_id, dst_ctx->did_i);
        }
    }
}

// PIZ only: this is called by the main thread after it integrated all the dictionary fragment read from disk for one VB.
// Here we hand over the integrated dictionaries to the VB - in preparation for the Compute Thread to use them.
// We overlay the z_file's dictionaries and word lists to the vb. these data remain unchanged - neither
// the vb nor the dispatcher thread will ever change snips placed in these. the dispatcher thread may append
// the dictionary and word list as new fragments become available from subsequent VBs. If the memory is not 
// sufficient, the dispatcher thread will "abandon" this memory, leaving it to the VB to continue to use it
// while starting a larger dict/word_list on a fresh memory allocation.
void ctx_overlay_dictionaries_to_vb (VBlock *vb)
{
    for (DidIType did_i=0; did_i < MAX_DICTS; did_i++) {
        Context *zctx = ZCTX(did_i);
        Context *vctx = CTX(did_i);

        if (!zctx->dict_id.num) continue;

        // we create a VB contexts even if there are now dicts (perhaps skipped due to flag) - needed for containers to work
        vctx->did_i       = did_i;
        vctx->dict_id     = zctx->dict_id;
        vctx->last_line_i = LAST_LINE_I_INIT;
        memcpy ((char*)vctx->tag_name, zctx->tag_name, sizeof (vctx->tag_name));

        set_dict_id_to_did_i_map (vb->dict_id_to_did_i_map, vctx->dict_id, did_i);

        ctx_init_iterator (vctx);

        if (buf_is_alloc (&zctx->dict))
            buf_overlay (vb, &vctx->dict, &zctx->dict, "ctx->dict");    
        
        if (buf_is_alloc (&zctx->word_list))
            buf_overlay (vb, &vctx->word_list, &zctx->word_list, "ctx->word_list");
    }
    vb->num_contexts = z_file->num_contexts;
}

// used by random_access_show_index
CtxNode *ctx_get_node_by_word_index (ConstContextP ctx, WordIndex word_index)
{
    ARRAY (CtxNode, nodes, ctx->nodes);

    for (uint64_t i=0; i < ctx->nodes.len; i++)
        if (nodes[i].word_index.n == word_index) return &nodes[i];

    ABORT_R ("ctx_get_node_by_word_index failed to find word_index=%d in did_i=%u", word_index, ctx->did_i);
}

// PIZ: get snip by normal word index (doesn't support WORD_INDEX_* - returns "")
const char *ctx_get_snip_by_word_index (ConstContextP ctx, WordIndex word_index, STRp(*snip))
{
    ASSERT (buf_is_alloc (&ctx->word_list), "word_list is not allocated for ctx=%s", ctx->tag_name);

    ASSERT ((uint32_t)word_index < ctx->word_list.len, "word_index=%d out of range: word_list.len=%u for ctx=%s",
            word_index, (uint32_t)ctx->word_list.len, ctx->tag_name);

    if (word_index < 0) {
        static const char *empty="";
        if (snip) *snip = empty;
        if (snip_len) *snip_len = 0;
        return empty;
    }

    CtxWord *word = ENT (CtxWord, ctx->word_list, word_index);
    const char *my_snip = ENT (const char, ctx->dict, word->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = word->snip_len;

    return my_snip; 
}

// PIZ: returns word index of the snip, or WORD_INDEX_NONE if it is not in the dictionary
WordIndex ctx_get_word_index_by_snip (ConstContextP ctx, const char *snip)
{
    unsigned snip_len = strlen (snip);
    ARRAY (const CtxWord, words, ctx->word_list);
    ARRAY (const char, dict, ctx->dict);

    for (WordIndex i=0; i < words_len; i++)
        if (words[i].snip_len == snip_len && !memcmp (&dict[words[i].char_index], snip, snip_len))
            return i;

    return WORD_INDEX_NONE;
}

// ZIP
const char *ctx_get_snip_by_zf_node_index (const Buffer *nodes, const Buffer *dict, WordIndex node_index, 
                                           const char **snip, uint32_t *snip_len)
{
    ASSERT (node_index >= 0 && node_index < nodes->len, "node_index=%d out of range, nodes.len=%"PRIu64, node_index, nodes->len);

    CtxNode *node = ENT (CtxNode, *nodes, node_index);
    const char *my_snip = ENT (const char, *dict, node->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = node->snip_len;

    return my_snip; 
}

static Buffer *sorter_cmp_counts = NULL; // for use by sorter_cmp - used only in vblock_i=1, so no thread safety issues
static int sorter_cmp(const void *a, const void *b)  
{ 
    return (int)(*ENT (int64_t, *sorter_cmp_counts, *(WordIndex *)b) -
                 *ENT (int64_t, *sorter_cmp_counts, *(WordIndex *)a));
}

void ctx_sort_dictionaries_vb_1(VBlock *vb)
{
    // thread safety note: no issues here, as this is run only by the compute thread of vblock_i=1
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {

        Context *ctx = CTX(did_i);

        if (ctx->no_vb1_sort || !ctx->nodes.len) continue;
        
        // prepare sorter array containing indices into ctx->nodes. We are going to sort it rather than sort nodes directly
        // as the b250 data contains node indices into ctx->nodes.
        static Buffer sorter = EMPTY_BUFFER;
        buf_alloc (vb, &sorter, 0, ctx->nodes.len, WordIndex, CTX_GROWTH, "sorter");
        for (WordIndex i=0; i < ctx->nodes.len; i++)
            NEXTENT (WordIndex, sorter) = i;

        // sort in ascending order of nodes->count
        sorter_cmp_counts = &ctx->counts; // communicate the ctx to sorter_cmp via a global var
        qsort (sorter.data, ctx->nodes.len, sizeof (WordIndex), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in nodes
        static Buffer unsorted_dict   = EMPTY_BUFFER;
        buf_move (vb, &unsorted_dict, vb, &ctx->dict);
        buf_alloc (vb, &ctx->dict, 0, unsorted_dict.len, char, CTX_GROWTH, "contexts->dict");
        ctx->dict.len = unsorted_dict.len;

        // note: we sort dict and update nodes->char_index and word_index. nodes and counts are not sorted.
        char *next = FIRSTENT (char, ctx->dict);
        for (WordIndex i=0; i < (WordIndex)ctx->nodes.len; i++) {
            WordIndex node_index = *ENT (WordIndex, sorter, i);
            CtxNode *node = ENT (CtxNode, ctx->nodes, node_index);

            if (node->word_index.n == WORD_INDEX_NONE) continue; // node was canceled in ctx_rollback

            const char *snip = ENT (char, unsorted_dict, node->char_index); 
            memcpy (next, snip, node->snip_len + 1 /* +1 for SNIP_SEP */);
            node->char_index = ENTNUM (ctx->dict, next);
            node->word_index.n = i;
            next += node->snip_len + 1;

            if (HAS_DEBUG_SEG(ctx)) {
                char printable_snip[node->snip_len+20];
                iprintf ("ctx_sort_dictionaries_vb_1: %s: word_index=%u snip=%s snip_len=%u count=%"PRId64"\n",
                          ctx->tag_name, i, str_print_snip (snip, node->snip_len, printable_snip), node->snip_len, 
                          *ENT (int64_t, ctx->counts, ctx->ol_nodes.len + node_index));
            }
        }

        buf_destroy (&sorter); // destroy and not free as it is first allocated by vb=0 and then again vb=1
        buf_destroy (&unsorted_dict);
    }
}

// ZIP only: run by main thread during zfile_output_processed_vb()
void ctx_update_stats (VBlock *vb)
{
    // zctx doesn't store b250, but we just use b250.len as a counter for displaying in genozip_show_sections
    for (ContextP vctx=CTX(0); vctx < CTX(vb->num_contexts); vctx++) {
        Context *zctx = ctx_get_zctx_from_vctx (vctx);
        if (!zctx) continue; // this can happen if FORMAT subfield appears, but no line has data for it

        zctx->b250.num_ctx_words += vctx->b250.num_ctx_words; // thread safety: no issues, this only updated only by the main thread
        
        if (vctx->ltype == LT_TEXT) // textual - count words added in seg_add_to_local_text
            zctx->local.num_ctx_words += zctx->local.num_ctx_words;
        else if (vctx->ltype != LT_BITMAP && vctx->ltype != LT_SEQUENCE && vctx->ltype != LT_CODEC) // numeric types - count numbers
            zctx->local.num_ctx_words += zctx->local.len;
    }
}

void ctx_free_context (Context *ctx, DidIType did_i)
{
    ctx_foreach_buffer (ctx, false, buf_free_do);

    memset ((char*)ctx->tag_name, 0, sizeof(ctx->tag_name));
    ctx->did_i = 0; 
    ctx->st_did_i = 0;
    ctx->ltype = 0;
    ctx->flags = (struct FlagsCtx){};
    ctx->pair_flags = (struct FlagsCtx){};
    ctx->dict_id.num = 0;
    ctx->pair_b250_iter = (SnipIterator){};
    ctx->lcodec = ctx->bcodec = ctx->lsubcodec_piz = 0;

    ctx->no_stons = ctx->pair_local = ctx->pair_b250 = ctx->stop_pairing = ctx->no_callback = ctx->line_is_luft_trans =
    ctx->local_param = ctx->no_vb1_sort = ctx->local_always = ctx->counts_section = ctx->no_drop_b250 =
    ctx->dynamic_size_local = ctx->is_stats_parent = ctx->local_compressed = ctx->dict_merged = 0;
    ctx->local_dep = 0;
    ctx->local_hash_prime = 0;
    ctx->num_new_entries_prev_merged_vb = 0;
    ctx->nodes_len_at_1_3 = ctx->nodes_len_at_2_3 = 0;
    ctx->luft_trans = 0;

    ctx->global_hash_prime = 0;
    ctx->merge_num = 0;
    ctx->txt_len = ctx->num_singletons = ctx->num_failed_singletons = 0;
    ctx->rback_b250_len = ctx->rback_local_len = ctx->rback_txt_len = ctx->rback_nodes_len = 0;    
    ctx->rback_last_txt_index = ctx->rback_last_txt_len = ctx->rback_num_singletons = 0;
    ctx->rback_last_value.i = 0;
    ctx->rback_last_delta = 0;

    mutex_destroy (ctx->mutex);

    ctx->iterator = (SnipIterator){};
    ctx->next_local = 0;

    ctx->last_line_i = ctx->last_sample_i = ctx->ctx_specific = 0; 
    ctx->last_value.i = 0;
    ctx->last_delta = 0;
    ctx->last_txt_index = ctx->last_txt_len = 0;
    ctx->semaphore = 0;
    ctx->compressor_time = 0;
    ctx->tag_i = 0;

    ctx->last_snip = 0; 
    ctx->last_snip_len = 0;
    ctx->last_snip_ni = 0;

    memset (ctx->tag_name, 0, sizeof(ctx->tag_name));
}

// Called by file_close ahead of freeing File memory containing contexts
void ctx_destroy_context (Context *ctx, DidIType did_i)
{
    ctx_foreach_buffer (ctx, false, buf_destroy_do);
    mutex_destroy (ctx->mutex);

    // test that ctx_free_context indeed frees everything
    #ifdef DEBUG
        ctx_free_context (ctx, did_i);

        #define REL_LOC(field) ((char*)(&ctx->field) - (char*)ctx)
        for (char *c=(char*)ctx; c  < (char*)(ctx+1); c++)
            if (*c) {
                fprintf (stderr, "relative location for debugging: %"PRIu64"\n", REL_LOC(num_singletons)); // help find the offending field 
                ABORT ("ctx_free_context didn't fully clear the context did_i=%u, byte %u != 0", did_i, (unsigned)(c - (char*)ctx)); 
            }
    #endif

    memset (ctx, 0, sizeof (Context));
}

void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local /* true = local, false = b250 */)
{
    char dump_fn[MAX_TAG_LEN + 50];
    sprintf (dump_fn, "%s.%05u.%s", ctx->tag_name, vb->vblock_i, local ? "local" : "b250");
    
    bool success = local ? buf_dump_to_file (dump_fn, &ctx->local, lt_desc[ctx->ltype].width, false, true, true, false)
                         : buf_dump_to_file (dump_fn, &ctx->b250, 1, false, true, true, false);

    ASSERTW (success, "Warning: ctx_dump_binary failed to output file %s: %s", dump_fn, strerror (errno));
}

// -------------------------------------
// ZIP: Compress and output dictionaries
// -------------------------------------

static Context *frag_ctx;
static const CtxNode *frag_next_node;
static Codec frag_codec = CODEC_UNKNOWN;
static unsigned frag_size = 0;

// compress the dictionary fragment - either an entire dict, or divide it to fragments if large to allow multi-threaded
// compression and decompression
static void ctx_prepare_for_dict_compress (VBlockP vb)
{
    while (frag_ctx < ZCTX(z_file->num_contexts)) {

        if (!frag_next_node) {
            if (!frag_ctx->nodes.len) {
                frag_ctx++;
                continue; // unused context
            }
            frag_next_node = FIRSTENT (const CtxNode, frag_ctx->nodes);

            ASSERT (frag_next_node->char_index + frag_next_node->snip_len <= frag_ctx->dict.len, 
                    "Corrupt nodes in ctx=%.8s did_i=%u", frag_ctx->tag_name, (int)(frag_ctx - z_file->contexts));

            ASSERT (frag_ctx->dict_id.num, "dict_id=0 for context \"%s\" did_i=%u", frag_ctx->tag_name, (unsigned)(frag_ctx - ZCTX(0)));

            // frag_size is set by default 1MB, but must be at least double the longest snip, or it will cause
            // mis-calculation of size_upper_bound in ctx_dict_read_one_vb
            frag_size = 1<<20; // 1MB
            for (const CtxNode *node=FIRSTENT (CtxNode, frag_ctx->nodes); node <= LASTENT (CtxNode, frag_ctx->nodes); node++)
                if (node->snip_len * 2 > frag_size)
                    frag_size = 2 * roundup2pow (node->snip_len); // must be power of 2 for ctx_dict_read_one_vb
        }

        vb->fragment_ctx   = frag_ctx;
        vb->fragment_start = ENT (char, frag_ctx->dict, frag_next_node->char_index);
        vb->fragment_codec = frag_codec;

        while (frag_next_node < AFTERENT (CtxNode, frag_ctx->nodes) && 
               vb->fragment_len + frag_next_node->snip_len + 1 < frag_size) {

            vb->fragment_len += frag_next_node->snip_len + 1;
            vb->fragment_num_words++;
            frag_next_node++;
        }

        if (frag_next_node == AFTERENT (CtxNode, frag_ctx->nodes)) {
            frag_ctx++;
            frag_next_node = NULL;
            frag_codec = CODEC_UNKNOWN; // we assign a fresh codec for each dictionary
            frag_size = 0;
        }

        if (vb->fragment_len) {

            // if its the first fragment - assign a codec
            if (vb->fragment_codec == CODEC_UNKNOWN)
                vb->fragment_codec = frag_codec =
                    codec_assign_best_codec (vb, vb->fragment_ctx, NULL, SEC_DICT);
            vb->ready_to_dispatch = true;
            break;
        }
    }
}

static void ctx_compress_one_dict_fragment (VBlockP vb)
{
    START_TIMER;

    SectionHeaderDictionary header = (SectionHeaderDictionary){ 
        .h.magic                 = BGEN32 (GENOZIP_MAGIC),
        .h.section_type          = SEC_DICT,
        .h.data_uncompressed_len = BGEN32 (vb->fragment_len),
        .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary)),
        .h.codec                 = vb->fragment_codec,
        .h.vblock_i              = BGEN32 (vb->vblock_i),
        .num_snips               = BGEN32 (vb->fragment_num_words),
        .dict_id                 = vb->fragment_ctx->dict_id
    };

    if (flag.show_dict) 
        iprintf ("%s (vb_i=%u, did=%u, num_snips=%u)\n", 
                 vb->fragment_ctx->tag_name, vb->vblock_i, vb->fragment_ctx->did_i, vb->fragment_num_words);
    
    if (ctx_is_show_dict_id (vb->fragment_ctx->dict_id))
        str_print_dict (vb->fragment_start, vb->fragment_len, flag.show_dict, false);

    if (flag.list_chroms && vb->fragment_ctx->did_i == CHROM)
        str_print_dict (vb->fragment_start, vb->fragment_len, false, VB_DT(DT_SAM) || VB_DT(DT_BAM));

    if (flag.show_time) codec_show_time (vb, "DICT", vb->fragment_ctx->tag_name, vb->fragment_codec);

    comp_compress (vb, &vb->z_data, (SectionHeader*)&header, vb->fragment_start, NULL);

    COPY_TIMER (ctx_compress_one_dict_fragment)    

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

// called by main thread in zip_write_global_area
void ctx_compress_dictionaries (void)
{
    frag_ctx = ZCTX(0);
    frag_next_node = NULL;

    dispatcher_fan_out_task ("compress_dicts", NULL, PROGRESS_MESSAGE, "Writing dictionaries...", false, true, true, false, 0, 20000,
                             ctx_prepare_for_dict_compress, 
                             ctx_compress_one_dict_fragment, 
                             zfile_output_processed_vb);
}

// -------------------------------------
// PIZ: Read and decompress dictionaries
// -------------------------------------
static Section dict_sl = NULL; 
static Context *dict_ctx;

static void ctx_dict_read_one_vb (VBlockP vb)
{
    if (!sections_next_sec (&dict_sl, SEC_DICT))
        return; // we're done - no more SEC_DICT sections

    // create context if if section is skipped, for containters to work (skipping a section should be mirror in 
    // a container filter)
    bool new_ctx = (!dict_ctx || dict_sl->dict_id.num != dict_ctx->dict_id.num);
    if (new_ctx)
        dict_ctx = ctx_get_ctx_do (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts, dict_sl->dict_id, 0, 0);

    if (piz_is_skip_sectionz (SEC_DICT, dict_sl->dict_id)) goto done;
    
    int32_t offset = zfile_read_section (z_file, vb, dict_sl->vblock_i, &vb->z_data, "z_data", SEC_DICT, dict_sl);    
    SectionHeaderDictionary *header = 
        (offset != SECTION_SKIPPED) ? (SectionHeaderDictionary *)vb->z_data.data : NULL;

    vb->fragment_len = header ? BGEN32 (header->h.data_uncompressed_len) : 0;

    ASSERT (!header || header->dict_id.num == dict_sl->dict_id.num, "Expecting dictionary fragment with DictId=%s but found one with DictId=%s",
            dis_dict_id (dict_sl->dict_id).s, dis_dict_id (header->dict_id).s);

    // new context
    // note: in v9+ same-dict fragments are consecutive in the file, and all but the last are frag_size or a bit less, allowing pre-allocation
    if (header && new_ctx && z_file->genozip_version >= 9) {
        unsigned num_fragments=0; 
        for (Section sl=dict_sl; sl->dict_id.num == dict_ctx->dict_id.num; sl++) num_fragments++;

        // get size: for multi-fragment dictionaries, first fragment will be at or less than a power of 2, but more than the previous power of two.
        // this allows us to calculate the frag_size with which this dictionary was compressed and hence an upper bound on the size
        uint64_t size_upper_bound = (num_fragments == 1) ? vb->fragment_len : ((uint64_t)roundup2pow (vb->fragment_len) * (uint64_t)num_fragments);
        
        buf_alloc (evb, &dict_ctx->dict, 0, size_upper_bound, char, 0, "contexts->dict");
        dict_ctx->dict.param = ((uint64_t)num_fragments) << 32 | (uint64_t)vb->fragment_len; // for error reporting

        buf_set_overlayable (&dict_ctx->dict);
    }

    // when pizzing a v8 file, we run in single-thread since we need to do the following dictionary enlargement with which fragment
    if (z_file->genozip_version == 8) {
        buf_alloc (evb, &dict_ctx->dict, vb->fragment_len, 0, char, 0, "contexts->dict");
        buf_set_overlayable (&dict_ctx->dict);
    }

    if (header) {
        vb->fragment_ctx         = dict_ctx;
        vb->fragment_start       = ENT (char, dict_ctx->dict, dict_ctx->dict.len);
        dict_ctx->word_list.len += header ? BGEN32 (header->num_snips) : 0;
        dict_ctx->dict.len      += (uint64_t)vb->fragment_len;

        ASSERT (dict_ctx->dict.len <= dict_ctx->dict.size, "Dict %s len=%"PRIu64" exceeds allocated size=%"PRIu64" (num_fragments=%u fragment_1_len=%u)", 
                dict_ctx->tag_name, dict_ctx->dict.len, dict_ctx->dict.size, (unsigned)(dict_ctx->dict.param >> 32), (unsigned)(dict_ctx->dict.param & 0xffffffffULL));
    }

done: 
    // note: in cases we just "goto" here, no data is read, and a thread is needlessly created to decompress it
    // this is because the vb_i of the section needs to match the vb_i of the thread
    vb->ready_to_dispatch = true;
}

// entry point of compute thread of dictionary decompression
static void ctx_dict_uncompress_one_vb (VBlockP vb)
{
    if (!vb->fragment_ctx || flag.only_headers) goto done; // nothing to do in this thread
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)vb->z_data.data;

    ASSERT (vb->fragment_start + BGEN32 (header->h.data_uncompressed_len) <= AFTERENT (char, vb->fragment_ctx->dict), 
            "Buffer overflow when uncompressing dict=%s", vb->fragment_ctx->tag_name);

    // a hack for uncompressing to a location within the buffer - while multiple threads are uncompressing into 
    // non-overlappying regions in the same buffer in parallel
    Buffer copy = vb->fragment_ctx->dict;
    copy.data   = vb->fragment_start;
    zfile_uncompress_section (vb, header, &copy, NULL, 0, SEC_DICT); // NULL name prevents buf_alloc

done:
    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

static void ctx_dict_build_word_lists (void)
{    
    START_TIMER;

    Context *last = ZCTX((flag.show_headers && exe_type == EXE_GENOCAT) ? CHROM : z_file->num_contexts);

    for (Context *ctx=z_file->contexts; ctx <= last; ctx++) {

        if (!ctx->word_list.len || ctx->word_list.data) continue; // skip if 1. no words, or 2. already built

        buf_alloc (evb, &ctx->word_list, 0, ctx->word_list.len, CtxWord, 0, "contexts->word_list");
        buf_set_overlayable (&ctx->word_list);

        const char *word_start = ctx->dict.data;
        for (uint64_t snip_i=0; snip_i < ctx->word_list.len; snip_i++) {

            const char *c=word_start; while (*c) c++;

            *ENT (CtxWord, ctx->word_list, snip_i) = (CtxWord) {
                .snip_len   = c - word_start,
                .char_index = word_start - ctx->dict.data
            };

            word_start = c+1; // skip over the \0 separator
        }
    }

    COPY_TIMER_VB (evb, ctx_dict_build_word_lists);
}

// PIZ main thread
void ctx_read_all_dictionaries (void)
{
    START_TIMER;

    ctx_initialize_predefined_ctxs (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts);

    dict_sl = NULL;
    dict_ctx = NULL;

    dispatcher_fan_out_task (flag.reading_reference ? "read_dicts_ref" 
                            :flag.reading_chain     ? "read_dicts_chain" 
                            :flag.reading_kraken    ? "read_dicts_kraken" 
                            :                         "read_dicts",
                             NULL, PROGRESS_NONE, "Reading dictionaries...", 
                             flag.test, true, true, 
                             z_file->genozip_version == 8, // For v8 files, we read all fragments in the main thread as was the case in v8. This is because they are very small, and also we can't easily calculate the totel size of each dictionary.
                             0, 
                             10, // must be short - many dictionaries are just a few bytes
                             ctx_dict_read_one_vb, 
                             ctx_dict_uncompress_one_vb, 
                             NULL);

    // build word lists in z_file->contexts with dictionary data 
    ctx_dict_build_word_lists();

    // output the dictionaries if we're asked to
    if (flag.show_dict || flag.show_one_dict || flag.list_chroms) {
        for (uint32_t did_i=0; did_i < z_file->num_contexts; did_i++) {
            Context *ctx = ZCTX(did_i);
            if (!ctx->dict.len) continue;

            if (flag.list_chroms && ((!flag.luft && ctx->did_i == CHROM) || (flag.luft && ctx->did_i == VCF_oCHROM)))
                str_print_dict (STRb(ctx->dict), true, Z_DT(DT_SAM));
            
            if (flag.show_dict) 
                iprintf ("%s (did_i=%u, num_snips=%u, dict_size=%u bytes)\n", 
                         ctx->tag_name, did_i, (uint32_t)ctx->word_list.len, (uint32_t)ctx->dict.len);

            if (ctx_is_show_dict_id (ctx->dict_id))
                str_print_dict (STRb(ctx->dict), true, false);
        }
        iprint0 ("\n");

        if (exe_type == EXE_GENOCAT) exit_ok(); // if this is genocat - we're done
    }

    COPY_TIMER_VB (evb, ctx_read_all_dictionaries);
}

// -----------------------------
// ZIP & PIZ: SEC_COUNT sections
// -----------------------------

typedef struct {
    int64_t count;
    const char *snip;
} ShowCountsEnt;

static int show_counts_cmp (const void *a_, const void *b_)
{
    int64_t count_a = ((ShowCountsEnt *)a_)->count,
            count_b = ((ShowCountsEnt *)b_)->count;
    
    return count_b > count_a ? 1
         : count_a > count_b ? -1
         :                     0;  // don't use minus as numbers are 64b and return value is int
}

static void ctx_show_counts (Context *zctx)
{
    static Buffer show_counts_buf = EMPTY_BUFFER;
    buf_free (&show_counts_buf);
    buf_alloc (evb, &show_counts_buf, 0, zctx->counts.len, ShowCountsEnt, 0, "show_counts_buf");

    int64_t total=0;
    for (uint32_t i=0; i < zctx->counts.len; i++) {
        int64_t count = *ENT (int64_t, zctx->counts, i);
        if (!count) continue;

        total += count;
        NEXTENT (ShowCountsEnt, show_counts_buf) = (ShowCountsEnt){ 
            .count = count,
            .snip  = (command==ZIP) ? ctx_get_zf_nodes_snip (zctx, i) : ctx_get_words_snip (zctx, i)
        };
    }

    ARRAY (ShowCountsEnt, counts, show_counts_buf);
    qsort (counts, counts_len, sizeof (ShowCountsEnt), show_counts_cmp);

    iprintf ("Showing counts of %s (did_i=%u). Total items=%"PRId64" Number of categories=%u\n", zctx->tag_name, zctx->did_i, total, (unsigned)counts_len);    

    if (total)
        for (uint32_t i=0; i < counts_len; i++) 
            iprintf ("%s\t%"PRId64"\t%-4.2f%%\n", counts[i].snip, counts[i].count, 
                        100 * (float)counts[i].count / (float)total);

    if (exe_type == EXE_GENOCAT) exit_ok();
}

const char *ctx_get_snip_with_largest_count (DidIType did_i, int64_t *count)
{
    Context *ctx = ZCTX(did_i);
    ARRAY (CtxNode, nodes, ctx->nodes);
    ARRAY (int64_t, counts, ctx->counts);

    *count = -1;
    const char *snip = "";

    for (WordIndex i=0; i < nodes_len; i++)
        if (counts[i] > *count) {
            *count = counts[i];
            snip = ENT (char, ctx->dict, nodes[i].char_index);
        }

    return snip;
}

const char *ctx_get_vb_snip_ex (ConstContextP vctx, WordIndex vb_node_index, const char **snip, uint32_t *snip_len)
{
    CtxNode *node;
    const char *out_snip;

    if (vb_node_index < vctx->ol_nodes.len) { // is this entry from a previous vb (overlay buffer)
        node = ENT (CtxNode, vctx->ol_nodes, vb_node_index);
        out_snip = ENT (char, vctx->ol_dict, node->char_index);
    }
    else {
        node = ENT (CtxNode, vctx->nodes, vb_node_index - vctx->ol_nodes.len);
        out_snip = ENT (char, vctx->dict, node->char_index);
    }

    if (snip) *snip = out_snip;
    if (snip_len) *snip_len = node->snip_len;

    return out_snip;
}

void ctx_compress_counts (void)
{
    for (DidIType did_i=0; did_i < z_file->num_contexts; did_i++) {
        Context *ctx = ZCTX(did_i);

        if (flag.show_one_counts.num == dict_id_typeless (ctx->dict_id).num) 
            ctx_show_counts (ctx);

        if (ctx->counts_section && ctx->counts.len) {

            BGEN_u64_buf (&ctx->counts, NULL);       

            ctx->counts.len *= sizeof (int64_t);

            Codec codec = codec_assign_best_codec (evb, NULL, &ctx->counts, SEC_COUNTS);

            SectionHeaderCounts header = (SectionHeaderCounts){ 
                .h.magic                 = BGEN32 (GENOZIP_MAGIC),
                .h.section_type          = SEC_COUNTS,
                .h.data_uncompressed_len = BGEN32 (ctx->counts.len),
                .h.compressed_offset     = BGEN32 (sizeof(SectionHeaderCounts)),
                .h.codec                 = codec,
                .h.vblock_i              = 0,
                .nodes_param             = BGEN64 (ctx->nodes.param),
                .dict_id                 = ctx->dict_id
            };

            comp_compress (evb, &evb->z_data, (SectionHeader*)&header, ctx->counts.data, NULL);
            ctx->counts.len /= sizeof (int64_t);

            BGEN_u64_buf (&ctx->counts, NULL); // we need it for stats
        }
    }
}

// PIZ: called by the main threads from piz_read_global_area
void ctx_read_all_counts (void)
{
    Section sl = NULL;
    bool counts_shown=false;
    while (sections_next_sec (&sl, SEC_COUNTS))  {

        if (piz_is_skip_sectionz (SEC_COUNTS, sl->dict_id)) continue;

        Context *ctx = ctx_get_zctx (sl->dict_id, true);
        
        zfile_get_global_section (SectionHeaderCounts, SEC_COUNTS, sl, &ctx->counts, "counts");
        if (flag.only_headers) continue; // only show headers

        ctx->counts.len /= sizeof (int64_t);
        BGEN_u64_buf (&ctx->counts, NULL);

        ctx->nodes.param = BGEN64 (header.nodes_param);
        
        if (flag.show_one_counts.num == dict_id_typeless (sl->dict_id).num) {
            ctx_show_counts (ctx);
            counts_shown=true;
        }
    }

    ASSINP (counts_shown || !flag.show_one_counts.num || exe_type != EXE_GENOCAT, "There is no SEC_COUNTS section for %s", dis_dict_id (flag.show_one_counts).s);
}

TagNameEx ctx_tag_name_ex (ConstContextP ctx)
{
    TagNameEx s = {};

    if (!ctx) return (TagNameEx){ .s = "<none>" };

    unsigned start=0;
    if (z_file && (Z_DT(DT_VCF) || Z_DT(DT_BCF))) {
        if      (dict_id_is_vcf_format_sf (ctx->dict_id)) { memcpy (s.s, "FORMAT/", 7); start = 7; }
        else if (dict_id_is_vcf_info_sf   (ctx->dict_id)) { memcpy (s.s, "INFO/",   5); start = 5; }
    }

    memcpy (&s.s[start], ctx->tag_name, MAX_TAG_LEN);

    return s;
}

// rolls back a context to the rollback point registered in ctx_create_rollback_point
void ctx_rollback (VBlockP vb, Context *ctx)
{
    if (HAS_DEBUG_SEG(ctx)) 
        iprintf ("ctx_rollback: vb_i=%u %s: rolling back b250=%u nodes=%u\n", vb->vblock_i, ctx->tag_name, 
                 (int)(ctx->b250.len - ctx->rback_b250_len), (int)(ctx->nodes.len - ctx->rback_nodes_len));

    // if we evaluated this context since the rollback count - undo
    while (ctx->b250.len > ctx->rback_b250_len) { 
        (*ENT (int64_t, ctx->counts, LASTb250(ctx)))--;
        ctx->b250.len--;
    }

    // update all_the_same - true, if len is down to 1, and false if it is down to 0
    if (ctx->b250.len <= 1) ctx->flags.all_the_same = (bool)ctx->b250.len;

    // we can't actually remove nodes as they are refered to from the hash table, instead, we mark them as unused
    // this will prevent merging the word into the zctx->dict, and also prevent hash_get_entry_for_seg from returning this node if this word
    // is segged again in this VB - a new node will be created.
    for (WordIndex node_index_i = ctx->rback_nodes_len ; node_index_i < ctx->nodes.len ; node_index_i++) 
        ENT (CtxNode, ctx->nodes, node_index_i
        )->word_index.n = WORD_INDEX_NONE;
    
    ctx->local.len      = ctx->rback_local_len;
    ctx->txt_len        = ctx->rback_txt_len;
    ctx->num_singletons = ctx->rback_num_singletons;
    ctx->last_value     = ctx->rback_last_value;
    ctx->last_delta     = ctx->rback_last_delta;
    ctx->last_txt_index = ctx->rback_last_txt_index;
    ctx->last_txt_len   = ctx->rback_last_txt_len;
    ctx->last_line_i    = LAST_LINE_I_INIT; // undo "encountered in line"
    ctx->last_snip      = NULL; // cancel caching in ctx_create_node_do
}

void ctx_foreach_buffer(ContextP ctx, bool set_name, void (*func)(BufferP buf, const char *func, unsigned line)) 
{
    { Buffer *buf = &(ctx)->dict;        if (set_name) buf->name = "contexts->dict"        ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->b250;        if (set_name) buf->name = "contexts->b250"        ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->local;       if (set_name) buf->name = "contexts->local"       ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->pair;        if (set_name) buf->name = "contexts->pair"        ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->ol_dict;     if (set_name) buf->name = "contexts->ol_dict"     ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->ol_nodes;    if (set_name) buf->name = "contexts->ol_nodes"    ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->nodes;       if (set_name) buf->name = "contexts->nodes"       ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->counts;      if (set_name) buf->name = "contexts->counts"      ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->local_hash;  if (set_name) buf->name = "contexts->local_hash"  ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->global_hash; if (set_name) buf->name = "contexts->global_hash" ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->word_list;   if (set_name) buf->name = "contexts->word_list"   ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->con_cache;   if (set_name) buf->name = "contexts->con_cache"   ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->con_index;   if (set_name) buf->name = "contexts->con_index"   ; func (buf, __FUNCTION__, __LINE__); }  
    { Buffer *buf = &(ctx)->con_len;     if (set_name) buf->name = "contexts->con_len"     ; func (buf, __FUNCTION__, __LINE__); }  
}
