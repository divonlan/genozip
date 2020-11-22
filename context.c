// ------------------------------------------------------------------
//   context.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

/*
zip:
    1) during segregate - build ctx_context + dictionary for each dict_id

    2) during generate - convert snips in vcf to indexes into ctx->nodes

    3) merge back into the main (z_file) dictionaries - we use thread synchronization to make
    sure this happens in the sequencial order of variant blocks. this merging will also
    causes update of the word and char indices in ctx->nodes
*/

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

#define INITIAL_NUM_NODES 10000

MUTEX (wait_for_vb_1_mutex);

static inline void ctx_lock_do (VBlock *vb, pthread_mutex_t *mutex, const char *func, uint32_t code_line, const char *name, uint32_t param)
{
    //printf ("thread %u vb_i=%u LOCKING %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
    mutex_lock (*mutex);
    //printf ("thread %u vb_i=%u LOCKED %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
}
#define ctx_lock(vb, mutex, name, param) ctx_lock_do (vb, mutex, __FUNCTION__, __LINE__, name, param)

static inline void ctx_unlock_do (VBlock *vb, pthread_mutex_t *mutex, const char *func, uint32_t code_line, const char *name, uint32_t param)
{
    mutex_unlock (*mutex);
    //printf ("thread %u vb_i=%u UNLOCKED %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
}
#define ctx_unlock(vb, mutex, name, param) ctx_unlock_do (vb, mutex, __FUNCTION__, __LINE__, name, param)

void ctx_vb_1_lock (VBlockP vb)
{
    ASSERT0 (vb->vblock_i == 1, "Error: Only vb_i=1 can call ctx_vb_1_lock");

    ctx_lock (vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", 1);
}

void ctx_vb_1_unlock (VBlockP vb)
{
    ASSERT0 (vb->vblock_i == 1, "Error: Only vb_i=1 can call ctx_vb_1_unlock");

    ctx_unlock (vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", 1);
}

// ZIP: add a snip to the dictionary the first time it is encountered in the VCF file.
// the dictionary will be written to GENOZIP and used to reconstruct the MTF during decompression
typedef enum { DICT_VB, DICT_ZF, DICT_ZF_SINGLETON } DictType;
static inline CharIndex ctx_insert_to_dict (VBlock *vb_of_dict, Context *ctx, DictType type, const char *snip, uint32_t snip_len)
{
    Buffer *dict = (type == DICT_ZF_SINGLETON) ? &ctx->ol_dict : &ctx->dict;

    static const char *buf_name[3] = { "contexts->dict", "zf_ctx->dict", "zf_ctx->ol_dict" };
    buf_alloc (vb_of_dict, dict, MAX ((dict->len + snip_len + 1), INITIAL_NUM_NODES * MIN (10, snip_len)), 
               CTX_GROWTH, buf_name[type]);
    
    if (type == DICT_ZF) buf_set_overlayable (dict); // during merge
    
    CharIndex char_index = dict->len;
    char *dict_p = ENT (char, *dict, char_index);

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = SNIP_SEP; // dictionary have a SNIP_SEP separating snips, so that PIZ can generate word_list

    dict->len += snip_len + 1;
    return char_index;
}

// ZIP only (PIZ doesn't have nodes) nodes index to node - possibly in ol_nodes, or in nodes
CtxNode *ctx_node_vb_do (const Context *vb_ctx, WordIndex node_index, 
                         const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                         const char *func, uint32_t code_line)
{
    ASSERT (vb_ctx->dict_id.num, "Error in ctx_node_do: this vb_ctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index < vb_ctx->nodes.len + vb_ctx->ol_nodes.len, "Error in ctx_node_do: out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u. Caller: %s:%u",  
            vb_ctx->name, node_index, (uint32_t)vb_ctx->nodes.len, (uint32_t)vb_ctx->ol_nodes.len, func, code_line);

    bool is_ol = node_index < vb_ctx->ol_nodes.len; // is this entry from a previous vb (overlay buffer)

    CtxNode *node = is_ol ? ENT (CtxNode, vb_ctx->ol_nodes, node_index)
                          : ENT (CtxNode, vb_ctx->nodes, node_index - vb_ctx->ol_nodes.len);

    if (snip_in_dict) {
        const Buffer *dict = is_ol ? &vb_ctx->ol_dict : &vb_ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error in ctx_node_do: dict not allocated");

        ASSERT (node->char_index + (uint64_t)node->snip_len < dict->len, "Error in ctx_node_vb_do: snip of %s out of range: node->char_index=%"PRIu64" + node->snip_len=%u >= %s->len=%"PRIu64,
                vb_ctx->name, node->char_index, node->snip_len, is_ol ? "ol_dict" : "dict", dict->len);

        *snip_in_dict = ENT (char, *dict, node->char_index);
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// ZIP only (PIZ doesn't have nodes) nodes index to node - possibly in ol_nodes, or in nodes
CtxNode *ctx_node_zf_do (const Context *zf_ctx, int32_t node_index, 
                         const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                         const char *func, uint32_t code_line)
{
    ASSERT (zf_ctx->dict_id.num, "Error in ctx_node_do: this zf_ctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index > -2 - (int32_t)zf_ctx->ol_nodes.len && node_index < (int32_t)zf_ctx->nodes.len , "Error in ctx_node_do: out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u. Caller: %s:%u",  
            zf_ctx->name, node_index, (uint32_t)zf_ctx->nodes.len, (uint32_t)zf_ctx->ol_nodes.len, func, code_line);

    bool is_singleton = node_index < 0; // is this entry from a previous vb (overlay buffer)

    CtxNode *node = is_singleton ? ENT (CtxNode, zf_ctx->ol_nodes, -node_index - 2)
                                 : ENT (CtxNode, zf_ctx->nodes, node_index);

    if (snip_in_dict) {
        const Buffer *dict = is_singleton ? &zf_ctx->ol_dict : &zf_ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error in ctx_node_do: dict not allocated");

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

    for (unsigned i=0; i < ctx->word_list.len; i++)
        if (words[i].snip_len == snip_len && !memcmp (&ctx->dict.data[words[i].char_index], snip, snip_len))
            return i;

    return WORD_INDEX_NONE;
}

// PIZ only (uses word_list): returns word index, and advances the iterator
WordIndex ctx_get_next_snip (VBlock *vb, Context *ctx, uint8_t ctx_flags,
                             SnipIterator *override_iterator,   // if NULL, defaults to ctx->iterator
                             const char **snip, uint32_t *snip_len) // optional out 
{
    WordIndex word_index;
    ASSERT (ctx || override_iterator, "Error in ctx_get_next_snip: ctx is NULL. vb_i=%u", vb->vblock_i);

    // if the entire b250 in a VB consisted of word_index=0, we don't output the b250 to the file, and just 
    // consider it to always emit 0
    if (!buf_is_allocated (&ctx->b250)) {
        if (snip) {
            CtxWord *dict_word = FIRSTENT (CtxWord, ctx->word_list);
            *snip = &ctx->dict.data[dict_word->char_index];
            *snip_len = dict_word->snip_len;
        }
        return 0;
    }
    
    SnipIterator *iterator = override_iterator ? override_iterator : &ctx->iterator;
    
    if (!override_iterator && !iterator->next_b250) 
        iterator->next_b250 = FIRSTENT (uint8_t, ctx->b250); // initialize (GT data initializes to the beginning of each sample rather than the beginning of the data)

    // an imperfect test for overflow, but this should never happen anyway 
    ASSERT (override_iterator || iterator->next_b250 <= LASTENT (uint8_t, ctx->b250), "Error while reconstrucing line %u vb_i=%u: iterator for %s reached end of data",
            vb->line_i, vb->vblock_i, ctx->name);
            
    word_index = base250_decode (&iterator->next_b250, !(ctx_flags & CTX_FL_ALL_THE_SAME));  // if this line has no non-GT subfields, it will not have a ctx 

    // case: a Container item is missing (eg a subfield in a Sample, a FORMAT or Samples items in a file)
    if (word_index == WORD_INDEX_MISSING_SF) {
        if (snip) {
            *snip = NULL; // ignore this dict_id - don't even output a separator
            *snip_len = 0;
        }
    }

    // case: a subfield snip is empty, eg "AB::CD" (VCF GT Data) or "OA:Z:chr13,52863337,-,56S25M70S,0,;" (SAM OA optional field)
    else if (word_index == WORD_INDEX_EMPTY_SF) { 
        if (snip) {
            *snip = ""; // pointer to static empty string
            *snip_len = 0;
        }
    }

    else {
        if (word_index == WORD_INDEX_ONE_UP) 
            word_index = iterator->prev_word_index + 1;

        ASSERT (word_index < ctx->word_list.len, "Error while parsing line %u: word_index=%u is out of bounds - %s dictionary has only %u entries",
                vb->line_i, word_index, ctx->name, (uint32_t)ctx->word_list.len);

        CtxWord *dict_word = ENT (CtxWord, ctx->word_list, word_index);

        if (snip) {
            *snip = &ctx->dict.data[dict_word->char_index];
            *snip_len = dict_word->snip_len;
        }
    }
    
    iterator->prev_word_index = word_index;    

    return word_index;
}

// get next snip without advancing the iterator
const char *ctx_peek_next_snip (VBlock *vb, Context *ctx)
{
    SnipIterator save = ctx->iterator;

    const char *snip;
    uint32_t snip_len;
    ctx_get_next_snip (vb, ctx, ctx->flags, NULL, &snip, &snip_len);
    
    ctx->iterator = save; // restore

    return snip;
}

// Process and snip - return its node index, and enter it into the directory if its not already there. Called
// 1. During segregate - as snips are encountered in the data. No base250 encoding yet
// 2. During ctx_merge_in_vb_ctx_one_dict_id() - to enter snips into z_file->contexts - also encoding in base250
static WordIndex ctx_evaluate_snip_merge (VBlock *merging_vb, Context *zf_ctx, Context *vb_ctx, 
                                          const char *snip, uint32_t snip_len, uint32_t count,
                                          CtxNode **node, bool *is_new)  // out
{
    // if this turns out to be a singelton - i.e. a new snip globally - where it goes depends on whether its a singleton in the VB 
    // and whether we are allowed to move singletons to local. if its a singleton:
    // 1. we keep it in ol_nodes and the index in the node is negative to indicate that
    // 2. we insert it to ol_dict instead of dict - i.e. it doesn't get written the dict section
    // 3. we move it to the local section of this vb
    // 4. we set the word_index of its nodes to be the word_index of the SNIP_LOOKUP snip
    bool is_singleton_in_vb = (count == 1 && (vb_ctx->ltype == LT_TEXT) && !(vb_ctx->inst & CTX_INST_NO_STONS)); // is singleton in this VB

    // attempt to get the node from the hash table
    WordIndex node_index = hash_global_get_entry (zf_ctx, snip, snip_len, is_singleton_in_vb ? HASH_NEW_OK_SINGLETON_IN_VB : HASH_NEW_OK_NOT_SINGLETON, node);
    
    // case: existing non-singleton node. Possibly it was a singleton node before, that thanks to us hash_global_get_entry
    // converted to non-singleton - i.e. node_index is always >= 0.
    if (*node) { 
        *is_new = false;
        return node_index; // >= 0
    }
    
    // NEW SNIP globally - this snip was just added to the hash table - either as a regular or singleton node
    bool is_singleton_in_global = (node_index < 0);
    Buffer *nodes = is_singleton_in_global ? &zf_ctx->ol_nodes : &zf_ctx->nodes;
    ASSERT (nodes->len <= MAX_WORDS_IN_CTX, 
            "Error: too many words in ctx %s, max allowed number of words is is %u", zf_ctx->name, MAX_WORDS_IN_CTX);

    buf_alloc (evb, nodes, sizeof (CtxNode) * MAX(INITIAL_NUM_NODES, nodes->len), CTX_GROWTH, 
               is_singleton_in_global ? "zf_ctx->ol_nodes" : "zf_ctx->nodes");

    // set either singleton node or regular node with this snip
    *node = LASTENT (CtxNode, *nodes);
    memset (*node, 0, sizeof(CtxNode)); // safety
    (*node)->snip_len   = snip_len;
    (*node)->char_index = ctx_insert_to_dict (evb, zf_ctx, (is_singleton_in_global ? DICT_ZF_SINGLETON : DICT_ZF), snip, snip_len);

    // case: vb singleton turns out to be a global singleton - we add it to local and return the SNIP_LOOKUP node
    // instead of the singleton node (which is guaranteed to be non-singleton, and hence >= 0)
    // note: local is dedicated to singletons and contains nothing else, since CTX_INST_NO_STONS is not set
    if (node_index < 0) {
        seg_add_to_local_text (merging_vb, vb_ctx, snip, snip_len, 0);
        
        static char lookup = SNIP_LOOKUP;
        return ctx_evaluate_snip_merge (merging_vb, zf_ctx, vb_ctx, &lookup, 1, 2 /* not singleton */, node, is_new);
    }

    // case: not singleton - we return this (new) node
    else {
        buf_set_overlayable (&zf_ctx->nodes);
        (*node)->word_index.n = node_index;
        *is_new = true;
        return node_index; // >= 0
    }
}

WordIndex ctx_evaluate_snip_seg (VBlock *segging_vb, Context *vb_ctx, 
                                 const char *snip, uint32_t snip_len,
                                 bool *is_new /* out */)
{
    ASSERT0 (vb_ctx, "Error in ctx_evaluate_snip_seg: vb_ctx is NULL");

    if (!snip_len) {
        if (is_new) *is_new = false;
        return (!snip || (segging_vb->data_type == DT_VCF && *snip != ':')) ? WORD_INDEX_MISSING_SF : WORD_INDEX_EMPTY_SF;
    }

    WordIndex node_index_if_new = vb_ctx->ol_nodes.len + vb_ctx->nodes.len;
    
#ifdef DEBUG // time consuming and only needed during development
    ASSERT (strnlen (snip, snip_len) == snip_len, "Error in ctx_evaluate_snip_seg in vb=%u ctx=%s: snip_len=%u but unexpectedly has an 0 in its midst", 
            segging_vb->vblock_i, vb_ctx->name, snip_len);
#endif

    ASSERT (node_index_if_new <= MAX_NODE_INDEX, 
            "Error: ctx of %s is full (max allowed words=%u): ol_nodes.len=%u nodes.len=%u",
            vb_ctx->name, MAX_WORDS_IN_CTX, (uint32_t)vb_ctx->ol_nodes.len, (uint32_t)vb_ctx->nodes.len)

    // get the node from the hash table if it already exists, or add this snip to the hash table if not
    CtxNode *node;
    WordIndex existing_node_index = hash_get_entry_for_seg (segging_vb, vb_ctx, snip, snip_len, node_index_if_new, &node);
    if (existing_node_index != NODE_INDEX_NONE) {
        node->count++;
        if (is_new) *is_new = false;
        return existing_node_index; // snip found - we're done
    }
    
    // this snip isn't in the hash table - its a new snip
    ASSERT (vb_ctx->nodes.len < MAX_NODE_INDEX, "Error: too many words in dictionary %s", vb_ctx->name);

    buf_alloc (segging_vb, &vb_ctx->nodes, sizeof (CtxNode) * MAX(INITIAL_NUM_NODES, 1+vb_ctx->nodes.len), CTX_GROWTH, 
               "contexts->nodes");

    vb_ctx->nodes.len++; // new hash entry or extend linked list

    node = ctx_node_vb (vb_ctx, node_index_if_new, NULL, NULL);
    memset (node, 0, sizeof(CtxNode)); // safety
    node->snip_len     = snip_len;
    node->char_index   = ctx_insert_to_dict (segging_vb, vb_ctx, DICT_VB, snip, snip_len);
    node->word_index.n = node_index_if_new;
    node->count++;

    if (is_new) *is_new = true;
    return node_index_if_new;
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
        Context *vb_ctx = &vb->contexts[did_i];
        Context *zf_ctx = &z_file->contexts[did_i];

        // case: this context doesn't really exist (happens when incrementing num_contexts when adding RNAME and RNEXT in ctx_copy_ref_contigs_to_zf)
        if (!zf_ctx->mutex_initialized) continue;

        ctx_lock (vb, &zf_ctx->mutex, "zf_ctx", did_i);

        if (buf_is_allocated (&zf_ctx->dict)) {  // something already for this dict_id

            // overlay the global dict and nodes - these will not change by this (or any other) VB
            //fprintf (stderr,  ("ctx_clone: overlaying old dict %.8s, to vb_i=%u vb_did_i=z_did_i=%u\n", dict_id_printable (zf_ctx->dict_id).id, vb->vblock_i, did_i);
            buf_overlay (vb, &vb_ctx->ol_dict, &zf_ctx->dict, "ctx->ol_dict", did_i);   
            buf_overlay (vb, &vb_ctx->ol_nodes, &zf_ctx->nodes, "ctx->ol_nodes", did_i);   

            // overlay the hash table, that may still change by future vb's merging... this vb will only use
            // entries that are up to this merge_num
            buf_overlay (vb, &vb_ctx->global_hash, &zf_ctx->global_hash, "contexts->global_hash", did_i);
            vb_ctx->merge_num = zf_ctx->merge_num;
            vb_ctx->global_hash_prime = zf_ctx->global_hash_prime; // can never change
            vb_ctx->num_new_entries_prev_merged_vb = zf_ctx->num_new_entries_prev_merged_vb;
        }

        vb_ctx->did_i    = did_i;
        vb_ctx->dict_id  = zf_ctx->dict_id;
        // note: lcodec and bcodec are inherited in merge (see comment in zip_assign_best_codec)

        memcpy ((char*)vb_ctx->name, zf_ctx->name, sizeof (vb_ctx->name));

        vb->dict_id_to_did_i_map[vb_ctx->dict_id.map_key] = did_i;
        
        ctx_init_iterator (vb_ctx);

        ctx_unlock (vb, &zf_ctx->mutex, "zf_ctx", did_i);
    }

    vb->num_contexts = z_num_contexts;

    COPY_TIMER (ctx_clone);
}

static void ctx_initialize_ctx (Context *ctx, DidIType did_i, DictId dict_id, DidIType *dict_id_to_did_i_map)
{
    ctx->did_i    = did_i;
    ctx->st_did_i = DID_I_NONE;
    ctx->dict_id  = dict_id;
    
    memcpy ((char*)ctx->name, dict_id_printable (dict_id).id, DICT_ID_LEN);
    ((char*)ctx->name)[DICT_ID_LEN] = 0;

    ctx_init_iterator (ctx);
    
    if (dict_id_to_did_i_map[dict_id.map_key] == DID_I_NONE)
        dict_id_to_did_i_map[dict_id.map_key] = did_i;
}

// ZIP I/O thread: 
// 1. when starting to zip a new file, with pre-loaded external reference, we integrate the reference FASTA CONTIG
//    dictionary as the chrom dictionary of the new file
// 2. in SAM, DENOVO: after creating loaded_contigs from SQ records, we copy them to the RNAME dictionary
void ctx_copy_ref_contigs_to_zf (DidIType dst_did_i, ConstBufferP contigs_buf, ConstBufferP contigs_dict_buf)
{
    // note: in REF_INTERNAL it is possible that there are no contigs - unaligned SAM
    if (flag.reference == REF_INTERNAL && (!contigs_buf || !contigs_buf->len)) return;

    ASSERT0 (buf_is_allocated (contigs_buf) && buf_is_allocated (contigs_dict_buf),
             "Error in ctx_copy_ref_contigs_to_zf: expecting contigs and contigs_dict to be allocated");
    
    Context *zf_ctx = &z_file->contexts[dst_did_i];
    zf_ctx->inst = CTX_INST_NO_STONS;
    mutex_initialize (zf_ctx->mutex);

    // copy reference dict
    ARRAY (RefContig, contigs, *contigs_buf);

    buf_copy (evb, &zf_ctx->dict, contigs_dict_buf, 0,0,0, "z_file->contexts->dict");
    buf_set_overlayable (&zf_ctx->dict);

    // build nodes from reference word_list
    buf_alloc (evb, &zf_ctx->nodes, sizeof (CtxNode) * contigs_buf->len, 1, "z_file->contexts->nodes");
    buf_set_overlayable (&zf_ctx->nodes);
    zf_ctx->nodes.len = contigs_buf->len;

    for (unsigned i=0 ; i < zf_ctx->nodes.len; i++) {
        CtxNode *node = ENT (CtxNode, zf_ctx->nodes, i);
        node->char_index = contigs[i].char_index;
        node->snip_len   = contigs[i].snip_len;
        node->word_index = base250_encode (i);
        node->count      = 0;
    }
    
    // allocate and populate hash from zf_ctx->nodes
    hash_alloc_global (zf_ctx, zf_ctx->nodes.len);

    z_file->num_contexts = MAX (z_file->num_contexts, dst_did_i+1);
}

void ctx_initialize_for_zip (void)
{
    mutex_initialize (wait_for_vb_1_mutex);
}

// find the z_file context that corresponds to dict_id. It could be possibly a different did_i
// than in the vb - in case this dict_id is new to this vb, but another vb already inserted
// it to z_file
static Context *ctx_get_zf_ctx (DictId dict_id)
{
    DidIType z_num_contexts = __atomic_load_n (&z_file->num_contexts, __ATOMIC_RELAXED);

    for (DidIType did_i=0; did_i < z_num_contexts; did_i++)
        if (dict_id.num == z_file->contexts[did_i].dict_id.num) 
            return &z_file->contexts[did_i];

    return NULL;
}

// ZIP only: called by merging VBs to add a new dict to z_file - copying some stuff from vb_ctx
static Context *ctx_add_new_zf_ctx (VBlock *merging_vb, const Context *vb_ctx)
{
    // adding a new dictionary is proctected by a mutex. note that z_file->num_contexts is accessed by other threads
    // without mutex proction when searching for a dictionary - that's why we update it at the end, after the new
    // zf_ctx is set up with the new dict_id (ready for another thread to search it)
    ctx_lock (merging_vb, &z_file->dicts_mutex, "dicts_mutex", 0);

    // check if another thread raced and created this dict before us
    Context *zf_ctx = ctx_get_zf_ctx (vb_ctx->dict_id);
    if (zf_ctx) goto finish;

    ASSERT (z_file->num_contexts+1 < MAX_DICTS, // load num_contexts - this time with mutex protection - it could have changed
            "Error: z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    zf_ctx = &z_file->contexts[z_file->num_contexts];

    mutex_initialize (zf_ctx->mutex);

    zf_ctx->did_i    = z_file->num_contexts; 
    zf_ctx->st_did_i = vb_ctx->st_did_i;
    zf_ctx->dict_id  = vb_ctx->dict_id;
    memcpy ((char*)zf_ctx->name, vb_ctx->name, sizeof(zf_ctx->name));
    // note: lcodec is NOT copied here, see comment in zip_assign_best_codec

    // only when the new entry is finalized, do we increment num_contexts, atmoically , this is because
    // other threads might access it without a mutex when searching for a dict_id
    __atomic_store_n (&z_file->num_contexts, z_file->num_contexts+1, __ATOMIC_RELAXED); // stamp our merge_num as the ones that set the node_i

finish:
    ctx_unlock (merging_vb, &z_file->dicts_mutex, "dicts_mutex", 0);
    return zf_ctx;
}

void ctx_commit_codec_to_zf_ctx (VBlock *vb, Context *vb_ctx, bool is_lcodec)
{
    Context *zf_ctx  = ctx_get_zf_ctx (vb_ctx->dict_id);
    ASSERT (zf_ctx, "Error in ctx_commit_codec_to_zf_ctx: zf_ctx is missing for %s in vb=%u", vb_ctx->name, vb->vblock_i); // zf_ctx is expected to exist as this is called after merge

    { START_TIMER; 
      ctx_lock (vb, &zf_ctx->mutex, "zf_ctx", zf_ctx->did_i);
      COPY_TIMER_VB (vb, lock_mutex_zf_ctx);  
    }

    if (is_lcodec) zf_ctx->lcodec = vb_ctx->lcodec;
    else           zf_ctx->bcodec = vb_ctx->bcodec;

    ctx_unlock (vb, &zf_ctx->mutex, "zf_ctx->mutex", zf_ctx->did_i);
}

// ZIP only: this is called towards the end of compressing one vb - merging its dictionaries into the z_file 
// each dictionary is protected by its own mutex, and there is one z_file mutex protecting num_dicts.
// we are careful never to hold two muteces at the same time to avoid deadlocks
static void ctx_merge_in_vb_ctx_one_dict_id (VBlock *merging_vb, unsigned did_i)
{
    Context *vb_ctx = &merging_vb->contexts[did_i];

    // get the ctx or create a new one. note: ctx_add_new_zf_ctx() must be called before ctx_lock() because it locks the z_file mutex (avoid a deadlock)
    Context *zf_ctx  = ctx_get_zf_ctx (vb_ctx->dict_id);
    if (!zf_ctx) zf_ctx = ctx_add_new_zf_ctx (merging_vb, vb_ctx); 

    { START_TIMER; 
      ctx_lock (merging_vb, &zf_ctx->mutex, "zf_ctx", zf_ctx->did_i);
      COPY_TIMER_VB (merging_vb, lock_mutex_zf_ctx);  
    }

    START_TIMER; // note: careful not to count time spent waiting for the mutex
    //fprintf (stderr,  ("Merging dict_id=%.8s into z_file vb_i=%u vb_did_i=%u z_did_i=%u\n", dict_id_printable (vb_ctx->dict_id).id, merging_vb->vblock_i, did_i, z_did_i);

    zf_ctx->merge_num++; // first merge is #1 (first clone which happens before the first merge, will get vb-)
    zf_ctx->txt_len += vb_ctx->txt_len; // for stats
    zf_ctx->num_new_entries_prev_merged_vb = vb_ctx->nodes.len; // number of new words in this dict from this VB
    zf_ctx->num_singletons += vb_ctx->num_singletons; // add singletons created by seg (i.e. SNIP_LOOKUP_* in b250, and snip in local)

    if (vb_ctx->st_did_i != DID_I_NONE) 
        zf_ctx->st_did_i = vb_ctx->st_did_i; // we assign stats consolidation, but never revert back to DID_I_NONE

    // we assign VB a codec from zf_ctx, if not already assigned by Seg. See comment in zip_assign_best_codec
    if (!vb_ctx->lcodec) vb_ctx->lcodec = zf_ctx->lcodec;
    if (!vb_ctx->bcodec) vb_ctx->bcodec = zf_ctx->bcodec;
    
    if (!buf_is_allocated (&vb_ctx->dict)) goto finish; // nothing yet for this dict_id
 
    if (!buf_is_allocated (&zf_ctx->dict)) {
        // allocate hash table, based on the statitics gather by this first vb that is merging this dict and 
        // populate the hash table without needing to reevalate the snips (we know none are in the hash table, but all are in nodes and dict)
        if (zf_ctx->global_hash.size <= 1) { // only initial allocation in zip_dict_data_initialize
            uint32_t estimated_entries = hash_get_estimated_entries (merging_vb, zf_ctx, vb_ctx);
            hash_alloc_global (zf_ctx, estimated_entries);
        }
    }

    // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
    // (vb_ctx->nodes contains only new words, old words from previous vbs are in vb_ctx->ol_nodes)
    for (unsigned i=0; i < vb_ctx->nodes.len; i++) {
        CtxNode *vb_node = &((CtxNode *)vb_ctx->nodes.data)[i];
        
        CtxNode *zf_node;
        bool is_new;
        // use evb and not vb because zf_context is z_file (which belongs to evb)
        WordIndex zf_node_index = 
            ctx_evaluate_snip_merge (merging_vb, zf_ctx, vb_ctx, &vb_ctx->dict.data[vb_node->char_index], 
                                     vb_node->snip_len, vb_node->count, &zf_node, &is_new);

        ASSERT (zf_node_index >= 0 && zf_node_index < zf_ctx->nodes.len, "Error: zf_node_index=%d out of range - len=%i", zf_node_index, (uint32_t)vb_ctx->nodes.len);

        // set word_index to be indexing the global dict - to be used by vcf_zip_generate_genotype_one_section() and zip_generate_b250_section()
        if (is_new)
            vb_node->word_index = zf_node->word_index = base250_encode (zf_node_index);
        else 
            // a previous VB already already calculated the word index for this node. if it was done by vb_i=1,
            // then it is also re-sorted and the word_index is no longer the same as the node_index
            vb_node->word_index = zf_node->word_index;
    }

finish:
    COPY_TIMER_VB (merging_vb, ctx_merge_in_vb_ctx_one_dict_id)
    ctx_unlock (merging_vb, &zf_ctx->mutex, "zf_ctx->mutex", zf_ctx->did_i);
}

// ZIP only: merge new words added in this vb into the z_file.contexts, and compresses dictionaries.
void ctx_merge_in_vb_ctx (VBlock *merging_vb)
{
    START_TIMER;

    // vb_i=1 goes first, as it has the sorted dictionaries, other vbs can go in
    // arbitrary order. at the end of this function, vb_i releases the mutex it locked along time ago,
    // while the other vbs wait for vb_1 by attempting to lock the mutex
    if (merging_vb->vblock_i != 1) {
        ctx_lock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", merging_vb->vblock_i);
        ctx_unlock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", merging_vb->vblock_i);
    }
    
    ctx_verify_field_ctxs (merging_vb); // this was useful in the past to catch nasty thread issues

    // merge all contexts
    for (DidIType did_i=0; did_i < merging_vb->num_contexts; did_i++) 
        ctx_merge_in_vb_ctx_one_dict_id (merging_vb, did_i);

    // note: z_file->num_contexts might be larger than merging_vb->num_contexts at this point, for example:
    // vb_i=1 started, z_file is empty, created 20 contexts
    // vb_i=2 started, z_file is empty, created 10 contexts
    // vb_i=1 completes, merges 20 contexts to z_file, which has 20 contexts after
    // vb_i=2 completes, merges 10 contexts, of which 5 (for example) are shared with vb_i=1. Now z_file has 25 contexts after.

    COPY_TIMER_VB (merging_vb, ctx_merge_in_vb_ctx);
}

// PIZ: add aliases to dict_id_to_did_i_map
void ctx_map_aliases (VBlockP vb)
{
    if (!dict_id_aliases) return;

    for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
        for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) 
            if (dict_id_aliases[alias_i].dst.num == vb->contexts[did_i].dict_id.num && 
                vb->dict_id_to_did_i_map[dict_id_aliases[alias_i].alias.map_key] == DID_I_NONE)

                vb->dict_id_to_did_i_map[dict_id_aliases[alias_i].alias.map_key] = did_i;    
}

// returns an existing did_i in this vb, or DID_I_NONE if there isn't one
DidIType ctx_get_existing_did_i_if_not_found_by_inline (VBlockP vb, DictId dict_id)
{
    // a different dict_id is in the map, perhaps a hash-clash...
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) 
        if (dict_id.num == vb->contexts[did_i].dict_id.num) return did_i;

    // PIZ only: check if its an alias that's not mapped in ctx_map_aliases (due to contention)
    if (command != ZIP && dict_id_aliases) {
        for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
            if (dict_id.num == dict_id_aliases[alias_i].alias.num) { // yes! its an alias
                for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) 
                    if (dict_id_aliases[alias_i].dst.num == vb->contexts[did_i].dict_id.num) return did_i;
            }
    }

    return DID_I_NONE; // not found
}

ContextP ctx_get_existing_ctx_do (VBlockP vb, DictId dict_id) 
{
    DidIType did_i = ctx_get_existing_did_i(vb, dict_id); 
    return (did_i == DID_I_NONE) ? NULL : &vb->contexts[did_i]; 
}

// gets did_id if the dictionary exists, and creates a new dictionary if its the first time dict_id is encountered
// threads: no issues - called by PIZ for vb and zf (but dictionaries are immutable) and by Segregate (ZIP) on vb_ctx only
Context *ctx_get_ctx_if_not_found_by_inline (
    Context *contexts /* an array */, 
    DataType dt, 
    DidIType *dict_id_to_did_i_map, 
    DidIType did_i,
    DidIType *num_contexts, 
    DictId dict_id)
{
    // case: its not in mapper - mapper is occupied by another - perhaps it exists
    // and missing the opportunity to enter mapper - search for it
    if (did_i != DID_I_NONE)    
        for (did_i=0; did_i < *num_contexts; did_i++) 
            if (dict_id.num == contexts[did_i].dict_id.num) goto done;

    did_i = *num_contexts; // note: *num_contexts cannot be updated until ctx is initialized, see comment below
    
    Context *ctx = &contexts[did_i]; 

    //fprintf (stderr, "New context: dict_id=%.8s in did_i=%u \n", dict_id_print (dict_id), did_i);
    ASSERT (*num_contexts+1 < MAX_DICTS, 
            "Error: number of dictionaries is greater than MAX_DICTS=%u", MAX_DICTS);

    ctx_initialize_ctx (ctx, did_i, dict_id, dict_id_to_did_i_map);

    // thread safety: the increment below MUST be AFTER the initialization of ctx, bc piz_get_line_subfields
    // might be reading this data at the same time as the piz dispatcher thread adding more dictionaries
    (*num_contexts) = did_i + 1; 

done:
    ctx = &contexts[did_i];
    return ctx;
}

// called from seg_all_data_lines (ZIP) and ctx_read_all_dictionaries (PIZ) to initialize all
// primary field ctx's. these are not always used (e.g. when some are not read from disk due to genocat options)
// but we maintain their fixed positions anyway as the code relies on it
void ctx_initialize_primary_field_ctxs (Context *contexts /* an array */, 
                                        DataType dt,
                                        DidIType *dict_id_to_did_i_map,
                                        DidIType *num_contexts)
{
    *num_contexts = MAX (dt_fields[dt].num_fields, *num_contexts);

    for (int f=0; f < dt_fields[dt].num_fields; f++) {
        const char *fname  = dt_fields[dt].names[f];
        DictId dict_id = dict_id_field (dict_id_make (fname, strlen(fname)));
        Context *dst_ctx  = NULL;

        // check if its an alias (PIZ only)
        if (command != ZIP && dict_id_aliases) 
            for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
                if (dict_id.num == dict_id_aliases[alias_i].alias.num) 
                    dst_ctx = ctx_get_zf_ctx (dict_id_aliases[alias_i].dst);

        if (!dst_ctx) // normal field, not an alias
            ctx_initialize_ctx (&contexts[f], f, dict_id, dict_id_to_did_i_map);

        else { // an alias
            contexts[f].did_i = dst_ctx->did_i;
            dict_id_to_did_i_map[dict_id.map_key] = dst_ctx->did_i;        
        }
    }
}

// PIZ only: this is called by the I/O thread after it integrated all the dictionary fragment read from disk for one VB.
// Here we hand over the integrated dictionaries to the VB - in preparation for the Compute Thread to use them.
// We overlay the z_file's dictionaries and word lists to the vb. these data remain unchanged - neither
// the vb nor the dispatcher thread will ever change snips placed in these. the dispatcher thread may append
// the dictionary and word list as new fragments become available from subsequent VBs. If the memory is not 
// sufficient, the dispatcher thread will "abandon" this memory, leaving it to the VB to continue to use it
// while starting a larger dict/word_list on a fresh memory allocation.
void ctx_overlay_dictionaries_to_vb (VBlock *vb)
{
    for (DidIType did_i=0; did_i < MAX_DICTS; did_i++) {
        Context *zf_ctx = &z_file->contexts[did_i];
        Context *vb_ctx = &vb->contexts[did_i];

        if (!zf_ctx->dict_id.num) continue;

        if (buf_is_allocated (&zf_ctx->dict) && buf_is_allocated (&zf_ctx->word_list)) { 
            
            vb_ctx->did_i    = did_i;
            vb_ctx->dict_id  = zf_ctx->dict_id;
            memcpy ((char*)vb_ctx->name, zf_ctx->name, sizeof (vb_ctx->name));

            if (vb->dict_id_to_did_i_map[vb_ctx->dict_id.map_key] == DID_I_NONE)
                vb->dict_id_to_did_i_map[vb_ctx->dict_id.map_key] = did_i;

            ctx_init_iterator (vb_ctx);

            buf_overlay (vb, &vb_ctx->dict, &zf_ctx->dict, "ctx->dict", did_i);    
            buf_overlay (vb, &vb_ctx->word_list, &zf_ctx->word_list, "ctx->word_list", did_i);
        }
    }
    vb->num_contexts = z_file->num_contexts;
}

// used by random_access_show_index
CtxNode *ctx_get_node_by_word_index (Context *ctx, WordIndex word_index)
{
    ARRAY (CtxNode, nodes, ctx->nodes);

    for (uint64_t i=0; i < ctx->nodes.len; i++)
        if (nodes[i].word_index.n == word_index) return &nodes[i];

    ABORT ("ctx_get_node_by_word_index failed to find word_index=%d in did_i=%u", word_index, ctx->did_i);
    return NULL; // never reaches here
}

// PIZ: get snip by normal word index (doesn't support WORD_INDEX_*)
const char *ctx_get_snip_by_word_index (const Buffer *word_list, const Buffer *dict, WordIndex word_index, 
                                        const char **snip, uint32_t *snip_len)
{
    CtxWord *word = ENT (CtxWord, *word_list, word_index);
    const char *my_snip = ENT (const char, *dict, word->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = word->snip_len;

    return my_snip; 
}

// ZIP
const char *ctx_get_snip_by_node_index (const Buffer *nodes, const Buffer *dict, WordIndex node_index, 
                                        const char **snip, uint32_t *snip_len)
{
    CtxNode *node = ENT (CtxNode, *nodes, node_index);
    const char *my_snip = ENT (const char, *dict, node->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = node->snip_len;

    return my_snip; 
}

static Buffer *sorter_cmp_mtf = NULL; // for use by sorter_cmp - used only in vblock_i=1, so no thread safety issues
static int sorter_cmp(const void *a, const void *b)  
{ 
    return ENT (CtxNode, *sorter_cmp_mtf, *(uint32_t *)b)->count -
           ENT (CtxNode, *sorter_cmp_mtf, *(uint32_t *)a)->count;
}

void ctx_sort_dictionaries_vb_1(VBlock *vb)
{
    // thread safety note: no issues here, as this is run only by the compute thread of vblock_i=1
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {

        Context *ctx = &vb->contexts[did_i];

        if (ctx->inst & CTX_INST_NO_VB1_SORT) continue;
        
        // prepare sorter array containing indices into ctx->nodes. We are going to sort it rather than sort nodes directly
        // as the b250 data contains node indices into ctx->nodes.
        static Buffer sorter = EMPTY_BUFFER;
        buf_alloc (vb, &sorter, ctx->nodes.len * sizeof (int32_t), CTX_GROWTH, "sorter");
        for (WordIndex i=0; i < ctx->nodes.len; i++)
            NEXTENT (WordIndex, sorter) = i;

        // sort in ascending order of nodes->count
        sorter_cmp_mtf = &ctx->nodes; // communicate the ctx to sorter_cmp via a global var
        qsort (sorter.data, ctx->nodes.len, sizeof (WordIndex), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in nodes
        static Buffer old_dict = EMPTY_BUFFER;
        buf_move (vb, &old_dict, vb, &ctx->dict);

        buf_alloc (vb, &ctx->dict, old_dict.len, CTX_GROWTH, "contexts->dict");
        ctx->dict.len = old_dict.len;

        char *next = ctx->dict.data;
        for (WordIndex i=0; i < (WordIndex)ctx->nodes.len; i++) {
            WordIndex node_index = *ENT (WordIndex, sorter, i);
            CtxNode *node = ENT (CtxNode, ctx->nodes, node_index);
            memcpy (next, &old_dict.data[node->char_index], node->snip_len + 1 /* +1 for SNIP_SEP */);
            node->char_index   = next - ctx->dict.data;
            node->word_index.n = i;

            next += node->snip_len + 1;
        }

        buf_free (&sorter);
        buf_destroy (&old_dict);
    }
}

// for safety, verify that field ctxs are what they say they are. we had bugs in the past where they got mixed up due to
// delicate thread logic.
void ctx_verify_field_ctxs_do (VBlock *vb, const char *func, uint32_t code_line)
{
    for (int f=0; f < DTF(num_fields); f++) {

            Context *ctx = &vb->contexts[f];

            ASSERT (dict_id_fields[f] == ctx->dict_id.num,
                    "ctx_verify_field_ctxs called from %s:%u: dict_id mismatch with section type: f=%s ctx->dict_id=%s vb_i=%u",
                    func, code_line, (char*)DTF(names)[f], ctx->name, vb->vblock_i);
    }
}

// ZIP only: run by I/O thread during zfile_output_processed_vb()
void ctx_update_stats (VBlock *vb)
{
    // zf_ctx doesn't store node_i, but we just use node_i.len as a counter for displaying in genozip_show_sections
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        Context *vb_ctx = &vb->contexts[did_i];
    
        Context *zf_ctx = ctx_get_zf_ctx (vb_ctx->dict_id);
        if (!zf_ctx) continue; // this can happen if FORMAT subfield appears, but no line has data for it

        zf_ctx->node_i.len += vb_ctx->node_i.len; // thread safety: no issues, this only updated only by the I/O thread
    }
}

void ctx_free_context (Context *ctx)
{
    buf_free (&ctx->ol_dict);
    buf_free (&ctx->ol_nodes);
    buf_free (&ctx->dict);
    buf_free (&ctx->nodes);
    buf_free (&ctx->word_list);
    buf_free (&ctx->local_hash);
    buf_free (&ctx->global_hash);
    buf_free (&ctx->node_i);
    buf_free (&ctx->b250);
    buf_free (&ctx->local);
    buf_free (&ctx->con_cache);
    buf_free (&ctx->con_index);
    buf_free (&ctx->con_len);
    buf_free (&ctx->pair);
    
    ctx->node_i.len = 0; // VCF stores FORMAT length in here for stats, even if node_i is not allocated (and therefore buf_free will not cleanup)
    ctx->local.len = 0; // For callback ctxs, length is stored, but data is not copied to local and is kept in vb->txt_data
    ctx->dict_id.num = 0;
    ctx->iterator.next_b250 = ctx->pair_b250_iter.next_b250 = NULL;
    ctx->iterator.prev_word_index = ctx->pair_b250_iter.prev_word_index = 0;
    ctx->local_hash_prime = 0;
    ctx->global_hash_prime = 0;
    ctx->merge_num = 0;
    ctx->nodes_len_at_1_3 = ctx->nodes_len_at_2_3 = 0;
    ctx->txt_len = ctx->next_local = ctx->num_singletons = ctx->num_failed_singletons = 0;
    ctx->last_delta = 0;
    ctx->last_value.i = 0;
    ctx->last_line_i = 0;
    ctx->did_i = ctx->flags = ctx->inst = 0;
    ctx->ltype = 0;
    ctx->lcodec = ctx->bcodec = 0;
    memset ((char*)ctx->name, 0, sizeof(ctx->name));
    mutex_destroy (ctx->mutex);
}

// Called by file_close ahead of freeing File memory containing contexts
void ctx_destroy_context (Context *ctx)
{
    buf_destroy (&ctx->ol_dict);
    buf_destroy (&ctx->ol_nodes);
    buf_destroy (&ctx->dict);
    buf_destroy (&ctx->nodes);
    buf_destroy (&ctx->word_list);
    buf_destroy (&ctx->local_hash);
    buf_destroy (&ctx->global_hash);
    buf_destroy (&ctx->node_i);
    buf_destroy (&ctx->b250);
    buf_destroy (&ctx->con_cache);
    buf_destroy (&ctx->con_index);
    buf_destroy (&ctx->con_len);

    mutex_destroy (ctx->mutex);
}

void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local /* true = local, false = b250 */)
{
    char dump_fn[50];
    sprintf (dump_fn, "%s.%05u.%s", ctx->name, vb->vblock_i, local ? "local" : "b250");
    
    bool success = local ? file_put_buffer (dump_fn, &ctx->local, lt_desc[ctx->ltype].width)
                         : file_put_buffer (dump_fn, &ctx->b250, 1);

    ASSERTW (success, "Warning: ctx_dump_binary failed to output file %s: %s", dump_fn, strerror (errno));
}

// -------------------------------------
// ZIP: Compress and output dictionaries
// -------------------------------------

static Context *frag_ctx;
static const CtxNode *frag_next_node;
static Codec frag_codec = CODEC_UNKNOWN;

// compress the dictionary fragment - either an entire dict, or divide it to fragments if large to allow multi-threaded
// compression and decompression
static void ctx_prepare_for_dict_compress (VBlockP vb)
{
    // max fragment size - 1MB - a relatively small size to enable utilization of more cores, as only a handful of dictionaries
    // are expected to be big enough to have multiple fragments
    #define FRAGMENT_SIZE (1<<20)

    while (frag_ctx < &z_file->contexts[z_file->num_contexts]) {

        if (!frag_next_node) {
            if (!frag_ctx->nodes.len) {
                frag_ctx++;
                continue; // unused context
            }
            frag_next_node = FIRSTENT (const CtxNode, frag_ctx->nodes);
        }

        vb->fragment_ctx = frag_ctx;
        vb->fragment_start = ENT (char, frag_ctx->dict, frag_next_node->char_index);
        vb->fragment_codec = frag_codec;

        while (frag_next_node < AFTERENT (CtxNode, frag_ctx->nodes) && 
               vb->fragment_len + frag_next_node->snip_len < FRAGMENT_SIZE) {

            // we allow snips to be so large that it will cause the fragment to be FRAGMENT_SIZE/2 or less, which will cause
            // mis-calculation of size_upper_bound in ctx_dict_read_one_vb (if this ever becomes a problem, we can set FRAGMENT_SIZE
            // dynamically based on the largest snip in the dictionary)
            ASSERT (frag_next_node->snip_len < FRAGMENT_SIZE/2,
                    "Error: found a word in dict=%s that is larger than %u, the maximum supported by genozip", frag_ctx->name, FRAGMENT_SIZE/2);

            vb->fragment_len += frag_next_node->snip_len + 1;
            vb->fragment_num_words++;
            frag_next_node++;
        }

        if (frag_next_node == AFTERENT (CtxNode, frag_ctx->nodes)) {
            frag_ctx++;
            frag_next_node = NULL;
            frag_codec = CODEC_UNKNOWN;
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

    if (flag.show_dict) {
        fprintf (stderr, "%s (vb_i=%u, did=%u, num_snips=%u):\t", 
                 vb->fragment_ctx->name, vb->vblock_i, vb->fragment_ctx->did_i, vb->fragment_num_words);
        str_print_null_seperated_data (vb->fragment_start, vb->fragment_len, true, false);
    }
    
    if (dict_id_printable (vb->fragment_ctx->dict_id).num == flag.dict_id_show_one_dict.num)
        str_print_null_seperated_data (vb->fragment_start, vb->fragment_len, false, false);

    if (flag.list_chroms && vb->fragment_ctx->did_i == CHROM)
        str_print_null_seperated_data (vb->fragment_start, vb->fragment_len, false, vb->data_type == DT_SAM);

    if (flag.show_time) codec_show_time (vb, "DICT", vb->fragment_ctx->name, vb->fragment_codec);

    vb->z_data.name = "z_data"; // comp_compress requires that it is set in advance
    comp_compress (vb, &vb->z_data, false, (SectionHeader*)&header, vb->fragment_start, NULL);

    COPY_TIMER (ctx_compress_one_dict_fragment)    

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined.
}

// called by I/O thread in zip_write_global_area
void ctx_compress_dictionaries (void)
{
    frag_ctx = &z_file->contexts[0];
    frag_next_node = NULL;

    dispatcher_fan_out_task (NULL, PROGRESS_MESSAGE, "Writing dictionaries...", false, false,
                             ctx_prepare_for_dict_compress, 
                             ctx_compress_one_dict_fragment, 
                             zfile_output_processed_vb);
}

// -------------------------------------
// PIZ: Read and decompress dictionaries
// -------------------------------------
static const SectionListEntry *dict_sl = NULL; 
static ReadChromeType read_chrom_policy;
static Context *dict_ctx;

static void ctx_dict_read_one_vb (VBlockP vb)
{
    buf_alloc (vb, &vb->z_section_headers, 1 * sizeof(int32_t), 0, "z_section_headers"); // room for 1 section header

    if (!sections_get_next_section_of_type (&dict_sl, SEC_DICT, false, false))
        return; // we're done - no more SEC_DICT sections

    // cases where we can skip reading these dictionaries because we don't be using them
    bool is_chrom = (dict_sl->dict_id.num == dict_id_fields[CHROM]);
    if (read_chrom_policy == DICTREAD_CHROM_ONLY  && !is_chrom) return; // we're done - since CHROM is always the first dictionary
    if (read_chrom_policy == DICTREAD_EXCEPT_CHROM && is_chrom) goto done;

    if (piz_is_skip_sectionz (SEC_DICT, dict_sl->dict_id)) goto done;
    
    zfile_read_section (z_file, vb, dict_sl->vblock_i, &vb->z_data, "z_data", SEC_DICT, dict_sl);    
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)vb->z_data.data;

    vb->fragment_len = BGEN32 (header->h.data_uncompressed_len);

    // new context
    if (!dict_ctx || header->dict_id.num != dict_ctx->dict_id.num) {
        dict_ctx = ctx_get_ctx_do (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts, header->dict_id);

        // in v9+ same-dict fragments are consecutive in the file, and all but the last are FRAGMENT_SIZE or a bit less, allowing pre-allocation
        if (z_file->genozip_version >= 9) {
            unsigned num_fragments=0; 
            for (const SectionListEntry *sl=dict_sl; sl->dict_id.num == dict_ctx->dict_id.num; sl++) num_fragments++;

            // get size: for multi-fragment dictionaries, first fragment will be at or slightly less than FRAGMENT_SIZE, which is a power of 2.
            // this allows us to calculate the FRAGMENT_SIZE with which this file was compressed and hence an upper bound on the size
            uint32_t size_upper_bound = (num_fragments == 1) ? vb->fragment_len : roundup2pow (vb->fragment_len) * num_fragments;
            
            buf_alloc (evb, &dict_ctx->dict, size_upper_bound, 0, "context->dict");
            buf_set_overlayable (&dict_ctx->dict);
        }
    }

    // when pizzing a v8 file, we run in single-thread since we need to do the following dictionary enlargement with which fragment
    if (z_file->genozip_version == 8) {
        buf_alloc_more (evb, &dict_ctx->dict, vb->fragment_len, 0, char, 0, "context->dict");
        buf_set_overlayable (&dict_ctx->dict);
    }

    vb->fragment_ctx         = dict_ctx;
    vb->fragment_start       = ENT (char, dict_ctx->dict, dict_ctx->dict.len);
    dict_ctx->word_list.len += BGEN32 (header->num_snips);
    dict_ctx->dict.len      += vb->fragment_len;

done: 
    // note: in cases we just "goto" here, no data is read, and a thread is needlessly created to decompress it
    // this is because the vb_i of the section needs to match the vb_i of the thread
    vb->ready_to_dispatch = true;
}

// entry point of compute thread of dictionary decompression
static void ctx_dict_uncompress_one_vb (VBlockP vb)
{
    if (!vb->fragment_ctx || (flag.show_headers && exe_type == EXE_GENOCAT)) goto done; // nothing to do in this thread

    SectionHeaderDictionary *header = (SectionHeaderDictionary *)vb->z_data.data;

    // a hack for uncompressing to a location withing the buffer - while multiple threads are uncompressing into 
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

    for (Context *ctx=z_file->contexts; ctx < &z_file->contexts[z_file->num_contexts]; ctx++) {

        if (!ctx->word_list.len || ctx->word_list.data) continue; // skip if 1. no words, or 2. already built

        buf_alloc (evb, &ctx->word_list, ctx->word_list.len * sizeof (CtxWord), 0, "word_list");
        buf_set_overlayable (&ctx->word_list);

        const char *word_start = ctx->dict.data;
        for (uint32_t snip_i=0; snip_i < ctx->word_list.len; snip_i++) {

            const char *c=word_start; while (*c) c++;

            *ENT (CtxWord, ctx->word_list, snip_i) = (CtxWord) {
                .snip_len   = c - word_start,
                .char_index = word_start - ctx->dict.data
            };

            word_start = c+1; // skip over the \0 seperator
        }
    }

    COPY_TIMER_VB (evb, ctx_dict_build_word_lists);
}

void ctx_read_all_dictionaries (ReadChromeType read_chrom)
{
    START_TIMER;

    ctx_initialize_primary_field_ctxs (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_contexts);

    dict_sl = NULL;
    read_chrom_policy = read_chrom;
    dict_ctx = NULL;

    dispatcher_fan_out_task (NULL, PROGRESS_NONE, "Reading dictionaries...", 
                             flag.test, 
                             z_file->genozip_version == 8, // For v8 files, we read all fragments in the I/O thread as was the case in v8. This is because they are very small, and also we can't easily calculate the totel size of each dictionary.
                             ctx_dict_read_one_vb, 
                             ctx_dict_uncompress_one_vb, 
                             NULL);

    // build word lists in z_file->contexts with dictionary data 
    if (!(flag.show_headers && exe_type == EXE_GENOCAT))
        ctx_dict_build_word_lists();

    // output the dictionaries if we're asked to
    if (flag.show_dict || flag.dict_id_show_one_dict.num || flag.list_chroms) {
        for (uint32_t did_i=0; did_i < z_file->num_contexts; did_i++) {
            Context *ctx = &z_file->contexts[did_i];

#define MAX_PRINTABLE_DICT_LEN 100000

            if (dict_id_printable (ctx->dict_id).num == flag.dict_id_show_one_dict.num) 
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), false, false);
            
            if (flag.list_chroms && ctx->did_i == CHROM)
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), false, z_file->data_type == DT_SAM);
            
            if (flag.show_dict) {
                fprintf (stderr, "%s (did_i=%u, num_snips=%u):\t", ctx->name, did_i, (uint32_t)ctx->word_list.len);
                str_print_null_seperated_data (ctx->dict.data, (uint32_t)MIN(ctx->dict.len,MAX_PRINTABLE_DICT_LEN), true, false);
            }
        }
        fprintf (stderr, "\n");

        if (exe_type == EXE_GENOCAT) exit_ok; // if this is genocat - we're done
    }

    COPY_TIMER_VB (evb, ctx_read_all_dictionaries);
}
