// ------------------------------------------------------------------
//   move-to-front.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

/*
zip:
    1) during segregate - build mtf_context + dictionary for each dict_id

    2) during generate - convert snips in vcf to indexes into ctx->mtf

    3) merge back into the main (z_file) dictionaries - we use thread synchronization to make
    sure this happens in the sequencial order of variant blocks. this merging will also
    causes update of the word and char indices in ctx->mtf

    4) compress the incremental part of the dictionaries added by this VB

unzip:
    1) Dispatcher thread integrates the dictionaries fragments added by this VB

    2) Create MTF array mapping word indices to char indices (one array in z-file)

    3) Re-create genotype data by looking up words in the dictionaries
*/

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

#define INITIAL_NUM_NODES 10000

MUTEX (wait_for_vb_1_mutex);
MUTEX (compress_dictionary_data_mutex);

static inline void mtf_lock_do (VBlock *vb, pthread_mutex_t *mutex, const char *func, uint32_t code_line, const char *name, uint32_t param)
{
    //printf ("thread %u vb_i=%u LOCKING %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
    mutex_lock (*mutex);
    //printf ("thread %u vb_i=%u LOCKED %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
}
#define mtf_lock(vb, mutex, name, param) mtf_lock_do (vb, mutex, __FUNCTION__, __LINE__, name, param)

static inline void mtf_unlock_do (VBlock *vb, pthread_mutex_t *mutex, const char *func, uint32_t code_line, const char *name, uint32_t param)
{
    mutex_unlock (*mutex);
    //printf ("thread %u vb_i=%u UNLOCKED %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
}
#define mtf_unlock(vb, mutex, name, param) mtf_unlock_do (vb, mutex, __FUNCTION__, __LINE__, name, param)

void mtf_vb_1_lock (VBlockP vb)
{
    ASSERT0 (vb->vblock_i == 1, "Error: Only vb_i=1 can call mtf_vb_1_lock");

    mtf_lock (vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", 1);
}

// ZIP: add a snip to the dictionary the first time it is encountered in the VCF file.
// the dictionary will be written to GENOZIP and used to reconstruct the MTF during decompression
typedef enum { DICT_VB, DICT_ZF, DICT_ZF_SINGLTON } DictType;
static inline CharIndex mtf_insert_to_dict (VBlock *vb_of_dict, Context *ctx, DictType type, const char *snip, uint32_t snip_len)
{
    Buffer *dict = (type == DICT_ZF_SINGLTON) ? &ctx->ol_dict : &ctx->dict;

    static const char *buf_name[3] = { "contexts->dict", "zf_ctx->dict", "zf_ctx->ol_dict" };
    buf_alloc (vb_of_dict, dict, MAX ((dict->len + snip_len + 1), INITIAL_NUM_NODES * MIN (10, snip_len)), 
               CTX_GROWTH, buf_name[type] , ctx->did_i);
    
    if (type == DICT_ZF) buf_set_overlayable (dict); // during merge
    
    CharIndex char_index = dict->len;
    char *dict_p = ENT (char, *dict, char_index);

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = SNIP_SEP; // dictionary have a SNIP_SEP separating snips, so that PIZ can generate word_list

    dict->len += snip_len + 1;
    return char_index;
}

// ZIP only (PIZ doesn't have mtf) mtf index to node - possibly in ol_mtf, or in mtf
MtfNode *mtf_node_vb_do (const Context *vb_ctx, WordIndex node_index, 
                         const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                         const char *func, uint32_t code_line)
{
    ASSERT (vb_ctx->dict_id.num, "Error in mtf_node_do: this vb_ctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index < vb_ctx->mtf.len + vb_ctx->ol_mtf.len, "Error in mtf_node_do: out of range: dict=%s node_index=%d mtf.len=%u ol_mtf.len=%u. Caller: %s:%u",  
            vb_ctx->name, node_index, (uint32_t)vb_ctx->mtf.len, (uint32_t)vb_ctx->ol_mtf.len, func, code_line);

    bool is_ol = node_index < vb_ctx->ol_mtf.len; // is this entry from a previous vb (overlay buffer)

    MtfNode *node = is_ol ? ENT (MtfNode, vb_ctx->ol_mtf, node_index)
                          : ENT (MtfNode, vb_ctx->mtf, node_index - vb_ctx->ol_mtf.len);

    if (snip_in_dict) {
        const Buffer *dict = is_ol ? &vb_ctx->ol_dict : &vb_ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error in mtf_node_do: dict not allocated");

        ASSERT (node->char_index + (uint64_t)node->snip_len < dict->len, "Error in mtf_node_vb_do: snip of %s out of range: node->char_index=%"PRIu64" + node->snip_len=%u >= %s->len=%"PRIu64,
                vb_ctx->name, node->char_index, node->snip_len, is_ol ? "ol_dict" : "dict", dict->len);

        *snip_in_dict = ENT (char, *dict, node->char_index);
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// ZIP only (PIZ doesn't have mtf) mtf index to node - possibly in ol_mtf, or in mtf
MtfNode *mtf_node_zf_do (const Context *zf_ctx, int32_t node_index, 
                         const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                         const char *func, uint32_t code_line)
{
    ASSERT (zf_ctx->dict_id.num, "Error in mtf_node_do: this zf_ctx is not initialized (dict_id.num=0) - called from %s:%u", func, code_line);
    
    ASSERT (node_index > -2 - (int32_t)zf_ctx->ol_mtf.len && node_index < (int32_t)zf_ctx->mtf.len , "Error in mtf_node_do: out of range: dict=%s node_index=%d mtf.len=%u ol_mtf.len=%u. Caller: %s:%u",  
            zf_ctx->name, node_index, (uint32_t)zf_ctx->mtf.len, (uint32_t)zf_ctx->ol_mtf.len, func, code_line);

    bool is_singleton = node_index < 0; // is this entry from a previous vb (overlay buffer)

    MtfNode *node = is_singleton ? ENT (MtfNode, zf_ctx->ol_mtf, -node_index - 2)
                                 : ENT (MtfNode, zf_ctx->mtf, node_index);

    if (snip_in_dict) {
        const Buffer *dict = is_singleton ? &zf_ctx->ol_dict : &zf_ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error in mtf_node_do: dict not allocated");

        *snip_in_dict = &dict->data[node->char_index];
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// PIZ: search for a node matching this snip in a directory and return the node index. note that we do a linear
// search as PIZ doesn't have hash tables.
WordIndex mtf_search_for_word_index (Context *ctx, const char *snip, unsigned snip_len)
{
    MtfWord *words = (MtfWord *)ctx->word_list.data;

    for (unsigned i=0; i < ctx->word_list.len; i++)
        if (words[i].snip_len == snip_len && !memcmp (&ctx->dict.data[words[i].char_index], snip, snip_len))
            return i;

    return WORD_INDEX_NONE;
}

// PIZ only (uses word_list): returns word index, and advances the iterator
WordIndex mtf_get_next_snip (VBlock *vb, Context *ctx, 
                             SnipIterator *override_iterator,   // if NULL, defaults to ctx->iterator
                             const char **snip, uint32_t *snip_len) // optional out 
{
    ASSERT (ctx || override_iterator, "Error in mtf_get_next_snip: ctx is NULL. vb_i=%u", vb->vblock_i);

    SnipIterator *iterator = override_iterator ? override_iterator : &ctx->iterator;
    
    if (!override_iterator && !iterator->next_b250) { 
        ASSERT (buf_is_allocated (&ctx->b250), "Error in mtf_get_next_snip: b250 is unallocated. dict_id=%s", ctx->name);
        
        iterator->next_b250 = FIRSTENT (uint8_t, ctx->b250); // initialize (GT data initializes to the beginning of each sample rather than the beginning of the data)
    }

    // an imperfect test for overflow, but this should never happen anyway 
    ASSERT (override_iterator || iterator->next_b250 <= LASTENT (uint8_t, ctx->b250), "Error while reconstrucing line %u vb_i=%u: iterator for %s reached end of data",
            vb->line_i, vb->vblock_i, ctx->name);
            
    WordIndex word_index = base250_decode (&iterator->next_b250);  // if this line has no non-GT subfields, it will not have a ctx 

    // case: a subfield snip is missing - either the genotype data has less subfields than declared in FORMAT, or not provided at all for some (or all) samples.
    if (word_index == WORD_INDEX_MISSING_SF) {
        ASSERT (!ctx || (vb->data_type == DT_VCF && dict_id_is_vcf_format_sf (ctx->dict_id)), 
                "Error while reconstructing line %u vb_i=%u: BASE250_MISSING_SF unexpectedly found in b250 data of %s",
                vb->line_i, vb->vblock_i, ctx->dict_id.num ? ctx->name : "(no ctx)"); // there will be no context if this GT subfield was always missing - never appeared on any sample

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

        MtfWord *dict_word = ENT (MtfWord, ctx->word_list, word_index);

        if (snip) {
            *snip = &ctx->dict.data[dict_word->char_index];
            *snip_len = dict_word->snip_len;
        }
    }
    
    iterator->prev_word_index = word_index;    

    return word_index;
}

// get next snip without advancing the iterator
const char *mtf_peek_next_snip (VBlock *vb, Context *ctx)
{
    SnipIterator save = ctx->iterator;

    const char *snip;
    uint32_t snip_len;
    mtf_get_next_snip (vb, ctx, NULL, &snip, &snip_len);
    
    ctx->iterator = save; // restore

    return snip;
}

// Process and snip - return its node index, and enter it into the directory if its not already there. Called
// 1. During segregate - as snips are encountered in the data. No base250 encoding yet
// 2. During mtf_merge_in_vb_ctx_one_dict_id() - to enter snips into z_file->contexts - also encoding in base250
static WordIndex mtf_evaluate_snip_merge (VBlock *merging_vb, Context *zf_ctx, Context *vb_ctx, 
                                          const char *snip, uint32_t snip_len, uint32_t count,
                                          MtfNode **node, bool *is_new)  // out
{
    // if this turns out to be a singelton - i.e. a new snip globally - where it goes depends on whether its a singleton in the VB 
    // and whether we are allowed to move singletons to local. if its a singleton:
    // 1. we keep it in ol_mtf and the index in the node is negative to indicate that
    // 2. we insert it to ol_dict instead of dict - i.e. it doesn't get written the dict section
    // 3. we move it to the local section of this vb
    // 4. we set the word_index of its mtf to be the word_index of the SNIP_LOOKUP_TEXT snip
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
    Buffer *mtf = is_singleton_in_global ? &zf_ctx->ol_mtf : &zf_ctx->mtf;
    ASSERT (mtf->len <= MAX_WORDS_IN_CTX, 
            "Error: too many words in ctx %s, max allowed number of words is is %u", zf_ctx->name, MAX_WORDS_IN_CTX);

    buf_alloc (evb, mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, mtf->len), CTX_GROWTH, 
               is_singleton_in_global ? "zf_ctx->ol_mtf" : "zf_ctx->mtf", zf_ctx->did_i);

    // set either singleton node or regular node with this snip
    *node = LASTENT (MtfNode, *mtf);
    memset (*node, 0, sizeof(MtfNode)); // safety
    (*node)->snip_len   = snip_len;
    (*node)->char_index = mtf_insert_to_dict (evb, zf_ctx, (is_singleton_in_global ? DICT_ZF_SINGLTON : DICT_ZF), snip, snip_len);

    // case: vb singleton turns out to be a global singleton - we add it to local and return the SNIP_LOOKUP_TEXT node
    // instead of the singleton node (which is guaranteed to be non-singleton, and hence >= 0)
    if (node_index < 0) {
        seg_add_to_local_text (merging_vb, vb_ctx, snip, snip_len, 0);
        
        char lookup = SNIP_LOOKUP;
        return mtf_evaluate_snip_merge (merging_vb, zf_ctx, vb_ctx, &lookup, 1, 2 /* not singleton */, node, is_new);
    }

    // case: not singleton - we return this (new) node
    else {
        buf_set_overlayable (&zf_ctx->mtf);
        (*node)->word_index.n = node_index;
        *is_new = true;
        return node_index; // >= 0
    }
}

WordIndex mtf_evaluate_snip_seg (VBlock *segging_vb, Context *vb_ctx, 
                                 const char *snip, uint32_t snip_len,
                                 bool *is_new /* out */)
{
    ASSERT0 (vb_ctx, "Error in mtf_evaluate_snip_seg: vb_ctx is NULL");

    if (!snip_len) 
        return (!snip || (segging_vb->data_type == DT_VCF && *snip != ':')) ? WORD_INDEX_MISSING_SF : WORD_INDEX_EMPTY_SF;

    WordIndex node_index_if_new = vb_ctx->ol_mtf.len + vb_ctx->mtf.len;
    
#ifdef DEBUG // time consuming and only needed during development
    ASSERT (strnlen (snip, snip_len) == snip_len, "Error in mtf_evaluate_snip_seg: snip_len=%u but unexpectedly has an 0 in its midst", snip_len);
#endif

    ASSERT (node_index_if_new <= MAX_NODE_INDEX, 
            "Error: ctx of %s is full (max allowed words=%u): ol_mtf.len=%u mtf.len=%u",
            vb_ctx->name, MAX_WORDS_IN_CTX, (uint32_t)vb_ctx->ol_mtf.len, (uint32_t)vb_ctx->mtf.len)

    // get the node from the hash table if it already exists, or add this snip to the hash table if not
    MtfNode *node;
    WordIndex existing_node_index = hash_get_entry_for_seg (segging_vb, vb_ctx, snip, snip_len, node_index_if_new, &node);
    if (existing_node_index != NODE_INDEX_NONE) {
        node->count++;
        if (is_new) *is_new = false;
        return existing_node_index; // snip found - we're done
    }
    
    // this snip isn't in the hash table - its a new snip
    ASSERT (vb_ctx->mtf.len < MAX_NODE_INDEX, "Error: too many words in dictionary %s", vb_ctx->name);

    buf_alloc (segging_vb, &vb_ctx->mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, 1+vb_ctx->mtf.len), CTX_GROWTH, 
               "contexts->mtf", vb_ctx->did_i);

    vb_ctx->mtf.len++; // new hash entry or extend linked list

    node = mtf_node_vb (vb_ctx, node_index_if_new, NULL, NULL);
    memset (node, 0, sizeof(MtfNode)); // safety
    node->snip_len     = snip_len;
    node->char_index   = mtf_insert_to_dict (segging_vb, vb_ctx, DICT_VB, snip, snip_len);
    node->word_index.n = node_index_if_new;
    node->count++;

    if (is_new) *is_new = true;
    return node_index_if_new;
}

// ZIP only: overlay and/or copy the current state of the global context to the vb, ahead of compressing this vb.
void mtf_clone_ctx (VBlock *vb)
{
    unsigned z_num_dict_ids = __atomic_load_n (&z_file->num_dict_ids, __ATOMIC_RELAXED);

    START_TIMER; // including mutex wait time

    // note: because each dictionary has its own mutex, it is possible that we will see only a partial set
    // of dictionaries (eg some but not all of the fields) when we are arrive here while another thread is mid-way 
    // through merging and adding a bunch of dictionaries.
    // however z_num_dict_ids will always correctly state the number of dictionaries that are available.

    for (DidIType did_i=0; did_i < z_num_dict_ids; did_i++) {
        Context *vb_ctx = &vb->contexts[did_i];
        Context *zf_ctx = &z_file->contexts[did_i];

        ASSERT (zf_ctx->mutex_initialized, "Error: expected zf_ctx->mutex_initialized for did_i=%u", did_i);
        mtf_lock (vb, &zf_ctx->mutex, "zf_ctx", did_i);

        if (buf_is_allocated (&zf_ctx->dict)) {  // something already for this dict_id

            // overlay the global dict and mtf - these will not change by this (or any other) VB
            //fprintf (stderr,  ("mtf_clone_ctx: overlaying old dict %.8s, to vb_i=%u vb_did_i=z_did_i=%u\n", dict_id_printable (zf_ctx->dict_id).id, vb->vblock_i, did_i);
            buf_overlay (vb, &vb_ctx->ol_dict, &zf_ctx->dict, "ctx->ol_dict", did_i);   
            buf_overlay (vb, &vb_ctx->ol_mtf, &zf_ctx->mtf, "ctx->ol_mtf", did_i);   

            // overlay the hash table, that may still change by future vb's merging... this vb will only use
            // entries that are up to this merge_num
            buf_overlay (vb, &vb_ctx->global_hash, &zf_ctx->global_hash, "contexts->global_hash", did_i);
            vb_ctx->merge_num = zf_ctx->merge_num;
            vb_ctx->global_hash_prime = zf_ctx->global_hash_prime; // can never change
            vb_ctx->num_new_entries_prev_merged_vb = zf_ctx->num_new_entries_prev_merged_vb;
        }

        vb_ctx->did_i    = did_i;
        vb_ctx->dict_id  = zf_ctx->dict_id;
        vb_ctx->flags    = zf_ctx->flags;
        vb_ctx->ltype    = zf_ctx->ltype;
        vb_ctx->inst     = zf_ctx->inst;
        vb_ctx->lcomp = zf_ctx->lcomp;

        memcpy ((char*)vb_ctx->name, zf_ctx->name, sizeof (vb_ctx->name));

        vb->dict_id_to_did_i_map[vb_ctx->dict_id.map_key] = did_i;
        
        mtf_init_iterator (vb_ctx);

        mtf_unlock (vb, &zf_ctx->mutex, "zf_ctx", did_i);
    }

    vb->num_dict_ids = z_num_dict_ids;

    COPY_TIMER (vb->profile.mtf_clone_ctx)
}

static void mtf_initialize_ctx (Context *ctx, DataType dt, DidIType did_i, DictId dict_id, DidIType *dict_id_to_did_i_map)
{
    ctx->did_i   = did_i;
    ctx->dict_id = dict_id;
    
    if (command == ZIP) ctx->lcomp = COMP_BZ2; // default, may be changed in seg

    memcpy ((char*)ctx->name, dict_id_printable (dict_id).id, DICT_ID_LEN);
    ((char*)ctx->name)[DICT_ID_LEN] = 0;

    mtf_init_iterator (ctx);
    
    if (dict_id_to_did_i_map[dict_id.map_key] == DID_I_NONE)
        dict_id_to_did_i_map[dict_id.map_key] = did_i;
}

static Context *mtf_add_new_zf_ctx (VBlock *merging_vb, const Context *vb_ctx); // forward declaration

// ZIP I/O thread: when starting to zip a new file, with pre-loaded external reference, we integrate the reference FASTA CONTIG
// dictionary as the chrom dictionary of the new file
static void mtf_copy_reference_contig_to_chrom_ctx (void)
{
    ConstBufferP ref_contigs, ref_config_dict;
    ref_contigs_get (&ref_config_dict, &ref_contigs);

    ASSERT0 (buf_is_allocated (ref_contigs) && buf_is_allocated (ref_config_dict), 
             "Error in mtf_copy_reference_contig_to_chrom_ctx: expecting ref_contigs and ref_config_dict to be allocated");

    // Create chrom context, this is the first context so it will be did_i=0, hence the requirement that chrom is always the first field
    ASSERT0 (CHROM == 0, "Error: CHROM must be 0");
    
    Context copy_from_ctx;
    memset (&copy_from_ctx, 0, sizeof (Context));
    copy_from_ctx.dict_id = (DictId)dict_id_fields[CHROM];
    copy_from_ctx.inst    = CTX_INST_NO_STONS; // needs b250 node_index for random access;
    strcpy ((char*)copy_from_ctx.name, DTFZ(names)[CHROM]);
    
    mtf_add_new_zf_ctx (evb, &copy_from_ctx);

    // copy reference dict
    Context *zf_ctx = &z_file->contexts[CHROM];
    ARRAY (RefContig, contigs, *ref_contigs);

    buf_copy (evb, &zf_ctx->dict, ref_config_dict, 0,0,0, "z_file->contexts->dict", zf_ctx->did_i);
    buf_set_overlayable (&zf_ctx->dict);

    // build mtf from reference word_list
    buf_alloc (evb, &zf_ctx->mtf, sizeof (MtfNode) * ref_contigs->len, 1, "z_file->contexts->mtf", zf_ctx->did_i);
    buf_set_overlayable (&zf_ctx->mtf);
    zf_ctx->mtf.len = ref_contigs->len;

    for (unsigned i=0 ; i < zf_ctx->mtf.len; i++) {
        MtfNode *node = ENT (MtfNode, zf_ctx->mtf, i);
        node->char_index = contigs[i].char_index;
        node->snip_len   = contigs[i].snip_len;
        node->word_index = base250_encode (i);
        node->count      = 0;
    }
    
    // allocate and populate hash from zf_ctx->mtf
    hash_alloc_global (zf_ctx, zf_ctx->mtf.len);

    zfile_compress_dictionary_data (evb, zf_ctx, zf_ctx->mtf.len, zf_ctx->dict.data, zf_ctx->dict.len);
}

void mtf_initialize_for_zip (void)
{
    if (z_file->dicts_mutex_initialized) return;

    mutex_initialize (z_file->dicts_mutex);
    mutex_initialize (wait_for_vb_1_mutex);
    mutex_initialize (compress_dictionary_data_mutex);

    if (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE)
        mtf_copy_reference_contig_to_chrom_ctx();
}

// find the z_file context that corresponds to dict_id. It could be possibly a different did_i
// than in the vb - in case this dict_id is new to this vb, but another vb already inserted
// it to z_file
static Context *mtf_get_zf_ctx (DictId dict_id)
{
    DidIType z_num_dict_ids = __atomic_load_n (&z_file->num_dict_ids, __ATOMIC_RELAXED);

    for (DidIType did_i=0; did_i < z_num_dict_ids; did_i++)
        if (dict_id.num == z_file->contexts[did_i].dict_id.num) 
            return &z_file->contexts[did_i];

    return NULL;
}

// ZIP only: called by merging VBs to add a new dict to z_file - copying some stuff from vb_ctx
static Context *mtf_add_new_zf_ctx (VBlock *merging_vb, const Context *vb_ctx)
{
    // adding a new dictionary is proctected by a mutex. note that z_file->num_dict_ids is accessed by other threads
    // without mutex proction when searching for a dictionary - that's why we update it at the end, after the new
    // zf_ctx is set up with the new dict_id (ready for another thread to search it)
    mtf_lock (merging_vb, &z_file->dicts_mutex, "dicts_mutex", 0);

    // check if another thread raced and created this dict before us
    Context *zf_ctx = mtf_get_zf_ctx (vb_ctx->dict_id);
    if (zf_ctx) goto finish;

    ASSERT (z_file->num_dict_ids+1 < MAX_DICTS, // load num_dict_ids - this time with mutex protection - it could have changed
            "Error: z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    zf_ctx = &z_file->contexts[z_file->num_dict_ids];

    mutex_initialize (zf_ctx->mutex);

    zf_ctx->did_i      = z_file->num_dict_ids; 
    zf_ctx->dict_id    = vb_ctx->dict_id;
    zf_ctx->flags      = vb_ctx->flags;
    zf_ctx->inst       = vb_ctx->inst;
    zf_ctx->lcomp = vb_ctx->lcomp;
    zf_ctx->ltype      = vb_ctx->ltype;
    memcpy ((char*)zf_ctx->name, vb_ctx->name, sizeof(zf_ctx->name));

    // only when the new entry is finalized, do we increment num_dict_ids, atmoically , this is because
    // other threads might access it without a mutex when searching for a dict_id
    __atomic_store_n (&z_file->num_dict_ids, z_file->num_dict_ids+1, __ATOMIC_RELAXED); // stamp our merge_num as the ones that set the mtf_i

finish:
    mtf_unlock (merging_vb, &z_file->dicts_mutex, "dicts_mutex", 0);
    return zf_ctx;
}

// ZIP only: this is called towards the end of compressing one vb - merging its dictionaries into the z_file 
// each dictionary is protected by its own mutex, and there is one z_file mutex protecting num_dicts.
// we are careful never to hold two muteces at the same time to avoid deadlocks
static void mtf_merge_in_vb_ctx_one_dict_id (VBlock *merging_vb, unsigned did_i)
{
    Context *vb_ctx = &merging_vb->contexts[did_i];

    // get the ctx or create a new one. note: mtf_add_new_zf_ctx() must be called before mtf_lock() because it locks the z_file mutex (avoid a deadlock)
    Context *zf_ctx  = mtf_get_zf_ctx (vb_ctx->dict_id);
    if (!zf_ctx) zf_ctx = mtf_add_new_zf_ctx (merging_vb, vb_ctx); 

    { START_TIMER; 
      mtf_lock (merging_vb, &zf_ctx->mutex, "zf_ctx", zf_ctx->did_i);
      COPY_TIMER (merging_vb->profile.lock_mutex_zf_ctx);  
    }

    START_TIMER; // note: careful not to count time spent waiting for the mutex
    //fprintf (stderr,  ("Merging dict_id=%.8s into z_file vb_i=%u vb_did_i=%u z_did_i=%u\n", dict_id_printable (vb_ctx->dict_id).id, merging_vb->vblock_i, did_i, z_did_i);

    zf_ctx->merge_num++; // first merge is #1 (first clone which happens before the first merge, will get vb-)
    zf_ctx->txt_len += vb_ctx->txt_len; // for stats
    zf_ctx->num_new_entries_prev_merged_vb = vb_ctx->mtf.len; // number of new words in this dict from this VB
    zf_ctx->num_singletons += vb_ctx->num_singletons; // add singletons created by seg (i.e. SNIP_LOOKUP_* in b250, and snip in local)

    if (!buf_is_allocated (&vb_ctx->dict)) goto finish; // nothing yet for this dict_id

    uint64_t start_dict_len = zf_ctx->dict.len;
    uint64_t start_mtf_len  = zf_ctx->mtf.len;
 
    if (!buf_is_allocated (&zf_ctx->dict)) {
        // allocate hash table, based on the statitics gather by this first vb that is merging this dict and 
        // populate the hash table without needing to reevalate the snips (we know none are in the hash table, but all are in mtf and dict)
        if (zf_ctx->global_hash.size <= 1) { // only initial allocation in zip_dict_data_initialize
            uint32_t estimated_entries = hash_get_estimated_entries (merging_vb, zf_ctx, vb_ctx);
            hash_alloc_global (zf_ctx, estimated_entries);
        }
    }

    // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
    // (vb_ctx->mtf contains only new words, old words from previous vbs are in vb_ctx->ol_mtf)
    for (unsigned i=0; i < vb_ctx->mtf.len; i++) {
        MtfNode *vb_node = &((MtfNode *)vb_ctx->mtf.data)[i];
        
        MtfNode *zf_node;
        bool is_new;
        // use evb and not vb because zf_context is z_file (which belongs to evb)
        WordIndex zf_node_index = 
            mtf_evaluate_snip_merge (merging_vb, zf_ctx, vb_ctx, &vb_ctx->dict.data[vb_node->char_index], 
                                     vb_node->snip_len, vb_node->count, &zf_node, &is_new);

        ASSERT (zf_node_index >= 0 && zf_node_index < zf_ctx->mtf.len, "Error: zf_node_index=%d out of range - len=%i", zf_node_index, (uint32_t)vb_ctx->mtf.len);

        // set word_index to be indexing the global dict - to be used by vcf_zip_generate_genotype_one_section() and zip_generate_b250_section()
        if (is_new)
            vb_node->word_index = zf_node->word_index = base250_encode (zf_node_index);
        else 
            // a previous VB already already calculated the word index for this node. if it was done by vb_i=1,
            // then it is also re-sorted and the word_index is no longer the same as the node_index
            vb_node->word_index = zf_node->word_index;
    }

    // we now compress the dictionaries directly from z_file. note: we must continue to hold
    // the mutex during compression, lest another thread re-alloc the dictionary.
    const char *start_dict = &zf_ctx->dict.data[start_dict_len]; // we take the pointer AFTER the evaluate, since dict can be reallocted
    uint32_t added_chars   = (uint32_t)(zf_ctx->dict.len - start_dict_len);
    uint32_t added_words   = (uint32_t)(zf_ctx->mtf.len  - start_mtf_len);

    // compress incremental part of dictionary added by this vb. note: dispatcher calls this function in the correct order of VBs.
    if (added_chars) {
        // special optimization for the GL dictionary (it is ineffective with --optimize-GL that already optimizes GL)
        if (zf_ctx->dict_id.num == dict_id_FORMAT_GL && !flag_optimize_GL) 
            start_dict = gl_optimize_dictionary (merging_vb, &zf_ctx->dict, ENT (MtfNode, zf_ctx->mtf, start_mtf_len), start_dict_len, added_words);

        // we need to protect z_file->dict_data while we're writing to it. this ensures a single writer
        // to this data. we also need this mutex embedded in the zf_ctx->mutex, so that fragments of
        // a dictionary are written in the order they are created. 
        {   START_TIMER; 
            mtf_lock (merging_vb, &compress_dictionary_data_mutex, "compress_dictionary_data_mutex", merging_vb->vblock_i);
            COPY_TIMER (merging_vb->profile.lock_mutex_compress_dict);
        }  
        zfile_compress_dictionary_data (merging_vb, zf_ctx, added_words, start_dict, added_chars);
        mtf_unlock (merging_vb, &compress_dictionary_data_mutex, "compress_dictionary_data_mutex", merging_vb->vblock_i);
    }

finish:
    COPY_TIMER (merging_vb->profile.mtf_merge_in_vb_ctx_one_dict_id)
    mtf_unlock (merging_vb, &zf_ctx->mutex, "zf_ctx->mutex", zf_ctx->did_i);
}

// ZIP only: merge new words added in this vb into the z_file.contexts, and compresses dictionaries.
void mtf_merge_in_vb_ctx (VBlock *merging_vb)
{
    // vb_i=1 goes first, as it has the sorted dictionaries, other vbs can go in
    // arbitrary order. at the end of this function, vb_i releases the mutex it locked along time ago,
    // while the other vbs wait for vb_1 by attempting to lock the mutex
    if (merging_vb->vblock_i != 1) {
        mtf_lock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", merging_vb->vblock_i);
        mtf_unlock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", merging_vb->vblock_i);
    }
    
    mtf_verify_field_ctxs (merging_vb); // this was useful in the past to catch nasty thread issues

    // merge all contexts
    for (DidIType did_i=0; did_i < merging_vb->num_dict_ids; did_i++) 
        mtf_merge_in_vb_ctx_one_dict_id (merging_vb, did_i);

    // note: z_file->num_dict_ids might be larger than merging_vb->num_dict_ids at this point, for example:
    // vb_i=1 started, z_file is empty, created 20 contexts
    // vb_i=2 started, z_file is empty, created 10 contexts
    // vb_i=1 completes, merges 20 contexts to z_file, which has 20 contexts after
    // vb_i=2 completes, merges 10 contexts, of which 5 (for example) are shared with vb_i=1. Now z_file has 25 contexts after.

    if (merging_vb->vblock_i == 1)  
        mtf_unlock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", 1);
}

// PIZ only (no thread issues - dictionaries are immutable) - gets did_id if the dictionary exists, 
// or returns DID_I_NONE, if not. Called by subfield mapper functions only. NOTE: if ctx is local only (not dict)
// then it will NOT be found, as it is not in z_file
DidIType mtf_get_existing_did_i_from_z_file (DictId dict_id)
{
    // attempt to get did_i from dict_id mapper
    DidIType did_i = z_file->dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && z_file->contexts[did_i].dict_id.num == dict_id.num) return did_i;

    for (did_i=0; did_i < z_file->num_dict_ids; did_i++) 
        if (dict_id.num == z_file->contexts[did_i].dict_id.num) return did_i;

    return DID_I_NONE; // not found
}

// returns an existing did_i in this vb, or DID_I_NONE if there isn't one
DidIType mtf_get_existing_did_i (VBlockP vb, DictId dict_id)
{
    DidIType did_i = vb->dict_id_to_did_i_map[dict_id.map_key];
    if (did_i == DID_I_NONE || vb->contexts[did_i].dict_id.num == dict_id.num) return did_i;

    // a different dict_id is in the map, perhaps a hash-clash...
    for (did_i=0; did_i < vb->num_dict_ids; did_i++) 
        if (dict_id.num == vb->contexts[did_i].dict_id.num) return did_i;

    return DID_I_NONE; // not found
}

ContextP mtf_get_existing_ctx (VBlockP vb, DictId dict_id) 
{
    DidIType did_i = mtf_get_existing_did_i(vb, dict_id); 
    return (did_i == DID_I_NONE) ? NULL : &vb->contexts[did_i]; 
}


// gets did_id if the dictionary exists, and creates a new dictionary if its the first time dict_id is encountered
// threads: no issues - called by PIZ for vb and zf (but dictionaries are immutable) and by Segregate (ZIP) on vb_ctx only
Context *mtf_get_ctx_if_not_found_by_inline (
    Context *contexts /* an array */, 
    DataType dt, 
    DidIType *dict_id_to_did_i_map, 
    DidIType did_i,
    DidIType *num_dict_ids, 
    DictId dict_id)
{
    // case: its not in mapper - mapper is occupied by another - perhaps it exists
    // and missing the opportunity to enter mapper - search for it
    if (did_i != DID_I_NONE)    
        for (did_i=0; did_i < *num_dict_ids; did_i++) 
            if (dict_id.num == contexts[did_i].dict_id.num) goto done;

    // PIZ only: check if its an alias
    if (command != ZIP && dict_id_aliases) {
        for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
            if (dict_id.num == dict_id_aliases[alias_i].alias.num) { // yes! its an alias

                // get the context of the dst dictionary (+ create it if it doesn't exist yet)
                Context *dst_ctx = mtf_get_ctx_do (contexts, dt, dict_id_to_did_i_map, num_dict_ids, dict_id_aliases[alias_i].dst);

                ASSERT (dst_ctx, "Error in mtf_get_ctx_if_not_found_by_inline: cannot find ctx of %.*s which is the destination of the alias %.*s",
                        DICT_ID_LEN, dict_id_printable (dict_id_aliases[alias_i].dst).id, DICT_ID_LEN, dict_id_printable (dict_id_aliases[alias_i].alias).id);
                        
                // if have room in the map - add the alias there pointing to the dst did_i
                if (did_i == DID_I_NONE) 
                    dict_id_to_did_i_map[dict_id_aliases[alias_i].alias.map_key] = dst_ctx->did_i;

                return dst_ctx;
            }
    }

    did_i = *num_dict_ids; // note: *num_dict_ids cannot be updated until ctx is initialized, see comment below
    
    Context *ctx = &contexts[did_i]; 

    //fprintf (stderr, "New context: dict_id=%.8s in did_i=%u \n", dict_id_printable (dict_id).id, did_i);
    ASSERT (*num_dict_ids+1 < MAX_DICTS, 
            "Error: number of dictionaries is greater than MAX_DICTS=%u", MAX_DICTS);

    mtf_initialize_ctx (ctx, dt, did_i, dict_id, dict_id_to_did_i_map);

    // thread safety: the increment below MUST be AFTER the initialization of ctx, bc piz_get_line_subfields
    // might be reading this data at the same time as the piz dispatcher thread adding more dictionaries
    (*num_dict_ids) = did_i + 1; 

done:
    ctx = &contexts[did_i];
    return ctx;
}

// called from seg_all_data_lines (ZIP) and zfile_read_all_dictionaries (PIZ) to initialize all
// primary field ctx's. these are not always used (e.g. when some are not read from disk due to genocat options)
// but we maintain their fixed positions anyway as the code relies on it
void mtf_initialize_primary_field_ctxs (Context *contexts /* an array */, 
                                        DataType dt,
                                        DidIType *dict_id_to_did_i_map,
                                        DidIType *num_dict_ids)
{
    *num_dict_ids = MAX (dt_fields[dt].num_fields, *num_dict_ids);

    for (int f=0; f < dt_fields[dt].num_fields; f++) {
        const char *fname  = dt_fields[dt].names[f];
        DictId dict_id = dict_id_field (dict_id_make (fname, strlen(fname)));
        Context *dst_ctx  = NULL;

        // check if its an alias (PIZ only)
        if (command != ZIP && dict_id_aliases) 
            for (uint32_t alias_i=0; alias_i < dict_id_num_aliases; alias_i++)
                if (dict_id.num == dict_id_aliases[alias_i].alias.num) 
                    dst_ctx = mtf_get_zf_ctx (dict_id_aliases[alias_i].dst);

        if (!dst_ctx) // normal field, not an alias
            mtf_initialize_ctx (&contexts[f], dt, f, dict_id, dict_id_to_did_i_map);

        else { // an alias
            contexts[f].did_i = dst_ctx->did_i;
            dict_id_to_did_i_map[dict_id.map_key] = dst_ctx->did_i;        
        }
    }
}

// PIZ only: this is called by the I/O thread after reading a dictionary section 
void mtf_integrate_dictionary_fragment (VBlock *vb, char *section_data)
{    
    START_TIMER;

    // thread safety note: this function is called only from the piz dispatcher thread,
    // so no thread safety issues with this static buffer.
    static Buffer fragment = EMPTY_BUFFER;

    // thread-safety note: while the dispatcher thread is integrating new dictionary fragments,
    // compute threads might be using these dictionaries. This is ok, bc the dispatcher thread makes
    // sure we integrate dictionaries from vbs by order - so that running compute threads never
    // need to access the new parts of dictionaries. We also pre-allocate the dictionaries in
    // txtfile_genozip_to_txt_header() so that they don't need to be realloced. dict.len may be accessed
    // by compute threads, but its change is assumed to be atomic, so that no weird things will happen
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)section_data;

    ASSERT (header->h.section_type == SEC_DICT,
            "Error in mtf_integrate_dictionary_fragment: header->h.section_type=%s is not SEC_DICT", st_name(header->h.section_type));

    uint32_t num_snips = BGEN32 (header->num_snips);

    zfile_uncompress_section (vb, section_data, &fragment, "fragment", header->h.section_type);

    // special treatment if this is GL - de-optimize
    if (header->dict_id.num == dict_id_FORMAT_GL)
        gl_deoptimize (fragment.data, fragment.len);

    // in piz, the same did_i is used for z_file and vb contexts, meaning that in vbs there could be
    // a non-contiguous array of contexts (some are missing if not used by this vb)

    Context *zf_ctx = mtf_get_ctx_do (z_file->contexts, z_file->data_type, z_file->dict_id_to_did_i_map, &z_file->num_dict_ids, header->dict_id);
    zf_ctx->flags = header->h.flags;

    // append fragment to dict. If there is no room - old memory is abandoned (so that VBs that are overlaying
    // it continue to work uninterrupted) and a new memory is allocated, where the old dict is joined by the new fragment
    uint64_t dict_old_len = zf_ctx->dict.len;
    buf_alloc (evb, &zf_ctx->dict, zf_ctx->dict.len + fragment.len, CTX_GROWTH, "z_file->contexts->dict", header->h.section_type);
    buf_set_overlayable (&zf_ctx->dict);

    memcpy (AFTERENT (char, zf_ctx->dict), fragment.data, fragment.len);
    zf_ctx->dict.len += fragment.len;

    // extend word list memory - and calculate the new words. If there is no room - old memory is abandoned 
    // (so that VBs that are overlaying it continue to work uninterrupted) and a new memory is allocated
    buf_alloc (evb, &zf_ctx->word_list, (zf_ctx->word_list.len + (uint64_t)num_snips) * sizeof (MtfWord), CTX_GROWTH, 
               "z_file->contexts->word_list", zf_ctx->did_i);
    buf_set_overlayable (&zf_ctx->word_list);

    char *start = fragment.data;
    for (uint32_t snip_i=0; snip_i < num_snips; snip_i++) {

        MtfWord *word = &NEXTENT (MtfWord, zf_ctx->word_list);

        char *c=start; while (*c != SNIP_SEP) c++;

        word->snip_len   = c - start;
        word->char_index = dict_old_len + (start - fragment.data);

        start = c+1; // skip over the SNIP_SEP
    }

    buf_free (&fragment);

    COPY_TIMER(vb->profile.mtf_integrate_dictionary_fragment);
}

// PIZ only: this is called by the I/O thread after it integrated all the dictionary fragment read from disk for one VB.
// Here we hand over the integrated dictionaries to the VB - in preparation for the Compute Thread to use them.
// We overlay the z_file's dictionaries and word lists to the vb. these data remain unchanged - neither
// the vb nor the dispatcher thread will ever change snips placed in these. the dispatcher thread may append
// the dictionary and word list as new fragments become available from subsequent VBs. If the memory is not 
// sufficient, the dispatcher thread will "abandon" this memory, leaving it to the VB to continue to use it
// while starting a larger dict/word_list on a fresh memory allocation.
void mtf_overlay_dictionaries_to_vb (VBlock *vb)
{
    for (DidIType did_i=0; did_i < MAX_DICTS; did_i++) {
        Context *zf_ctx = &z_file->contexts[did_i];
        Context *vb_ctx = &vb->contexts[did_i];

        if (!zf_ctx->dict_id.num) continue;

        if (buf_is_allocated (&zf_ctx->dict) && buf_is_allocated (&zf_ctx->word_list)) { 
            
            vb_ctx->did_i   = did_i;
            vb_ctx->dict_id = zf_ctx->dict_id;
            vb_ctx->flags   = zf_ctx->flags;
            memcpy ((char*)vb_ctx->name, zf_ctx->name, sizeof (vb_ctx->name));

            if (vb->dict_id_to_did_i_map[vb_ctx->dict_id.map_key] == DID_I_NONE)
                vb->dict_id_to_did_i_map[vb_ctx->dict_id.map_key] = did_i;

            mtf_init_iterator (vb_ctx);

            buf_overlay (vb, &vb_ctx->dict, &zf_ctx->dict, "ctx->dict", did_i);    
            buf_overlay (vb, &vb_ctx->word_list, &zf_ctx->word_list, "ctx->word_list", did_i);
        }
    }
    vb->num_dict_ids = z_file->num_dict_ids;
}

// used by random_access_show_index
MtfNode *mtf_get_node_by_word_index (Context *ctx, WordIndex word_index)
{
    ARRAY (MtfNode, mtf, ctx->mtf);

    for (uint64_t i=0; i < ctx->mtf.len; i++)
        if (mtf[i].word_index.n == word_index) return &mtf[i];

    ABORT ("mtf_get_node_by_word_index failed to find word_index=%d in did_i=%u", word_index, ctx->did_i);
    return NULL; // never reaches here
}

const char *mtf_get_snip_by_word_index (const Buffer *word_list, const Buffer *dict, WordIndex word_index, 
                                        const char **snip, uint32_t *snip_len)
{
    MtfWord *word = ENT (MtfWord, *word_list, word_index);
    const char *my_snip = ENT (const char, *dict, word->char_index);
    
    if (snip) *snip = my_snip;
    if (snip_len) *snip_len = word->snip_len;

    return my_snip; 
}

static Buffer *sorter_cmp_mtf = NULL; // for use by sorter_cmp - used only in vblock_i=1, so no thread safety issues
static int sorter_cmp(const void *a, const void *b)  
{ 
    return ENT (MtfNode, *sorter_cmp_mtf, *(uint32_t *)b)->count -
           ENT (MtfNode, *sorter_cmp_mtf, *(uint32_t *)a)->count;
}

void mtf_sort_dictionaries_vb_1(VBlock *vb)
{
    // thread safety note: no issues here, as this is run only by the compute thread of vblock_i=1
    for (DidIType did_i=0; did_i < vb->num_dict_ids; did_i++) {

        Context *ctx = &vb->contexts[did_i];

        // prepare sorter array containing indices into ctx->mtf. We are going to sort it rather than sort mtf directly
        // as the b250 data contains node indices into ctx->mtf.
        static Buffer sorter = EMPTY_BUFFER;
        buf_alloc (vb, &sorter, ctx->mtf.len * sizeof (int32_t), CTX_GROWTH, "sorter", ctx->did_i);
        for (WordIndex i=0; i < ctx->mtf.len; i++)
            NEXTENT (WordIndex, sorter) = i;

        // sort in ascending order of mtf->count
        sorter_cmp_mtf = &ctx->mtf; // communicate the ctx to sorter_cmp via a global var
        qsort (sorter.data, ctx->mtf.len, sizeof (WordIndex), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in mtf
        static Buffer old_dict = EMPTY_BUFFER;
        buf_move (vb, &old_dict, vb, &ctx->dict);

        buf_alloc (vb, &ctx->dict, old_dict.len, CTX_GROWTH, "contexts->dict", did_i);
        ctx->dict.len = old_dict.len;

        char *next = ctx->dict.data;
        for (WordIndex i=0; i < (WordIndex)ctx->mtf.len; i++) {
            WordIndex node_index = *ENT (WordIndex, sorter, i);
            MtfNode *node = ENT (MtfNode, ctx->mtf, node_index);
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
void mtf_verify_field_ctxs_do (VBlock *vb, const char *func, uint32_t code_line)
{
    for (int f=0; f < DTF(num_fields); f++) {

            Context *ctx = &vb->contexts[f];

            ASSERT (dict_id_fields[f] == ctx->dict_id.num,
                    "mtf_verify_field_ctxs called from %s:%u: dict_id mismatch with section type: f=%s ctx->dict_id=%s vb_i=%u",
                    func, code_line, (char*)DTF(names)[f], ctx->name, vb->vblock_i);
    }
}

// ZIP only: run by I/O thread during zip_output_processed_vb()
void mtf_update_stats (VBlock *vb)
{
    // zf_ctx doesn't store mtf_i, but we just use mtf_i.len as a counter for displaying in genozip_show_sections
    for (DidIType did_i=0; did_i < vb->num_dict_ids; did_i++) {
        Context *vb_ctx = &vb->contexts[did_i];
    
        Context *zf_ctx = mtf_get_zf_ctx (vb_ctx->dict_id);
        if (!zf_ctx) continue; // this can happen if FORMAT subfield appears, but no line has data for it

        zf_ctx->mtf_i.len += vb_ctx->mtf_i.len; // thread safety: no issues, this only updated only by the I/O thread
    }
}

void mtf_free_context (Context *ctx)
{
    buf_free (&ctx->ol_dict);
    buf_free (&ctx->ol_mtf);
    buf_free (&ctx->dict);
    buf_free (&ctx->mtf);
    buf_free (&ctx->word_list);
    buf_free (&ctx->local_hash);
    buf_free (&ctx->global_hash);
    buf_free (&ctx->mtf_i);
    buf_free (&ctx->b250);
    buf_free (&ctx->local);
    ctx->mtf_i.len = 0; // VCF stores FORMAT length in here for stats, even if mtf_i is not allocated (and therefore buf_free will not cleanup)
    ctx->local.len = 0; // For callback ctxs, length is stored, but data is not copied to local and is kept in vb->txt_data
    ctx->dict_id.num = 0;
    ctx->iterator.next_b250 = NULL;
    ctx->iterator.prev_word_index =0;
    ctx->local_hash_prime = 0;
    ctx->global_hash_prime = 0;
    ctx->merge_num = 0;
    ctx->mtf_len_at_1_3 = ctx->mtf_len_at_2_3 = 0;
    ctx->txt_len = ctx->next_local = ctx->num_singletons = ctx->num_failed_singletons = 0;
    ctx->last_delta = 0;
    ctx->last_value.i = 0;
    ctx->last_line_i = 0;
    memset ((char*)ctx->name, 0, sizeof(ctx->name));
    mutex_destroy (ctx->mutex);
}

// Called by file_close ahead of freeing File memory containing contexts
void mtf_destroy_context (Context *ctx)
{
    buf_destroy (&ctx->ol_dict);
    buf_destroy (&ctx->ol_mtf);
    buf_destroy (&ctx->dict);
    buf_destroy (&ctx->mtf);
    buf_destroy (&ctx->word_list);
    buf_destroy (&ctx->local_hash);
    buf_destroy (&ctx->global_hash);
    buf_destroy (&ctx->mtf_i);
    buf_destroy (&ctx->b250);
    mutex_destroy (ctx->mutex);
}

static FILE *dump_file; // used by --dump-one-b250 and --dump-one-local
MUTEX (dump_mutex);

void mtf_initialize_binary_dump (const char *field, DictId *dict_id, const char *filename_ext)
{
    *dict_id = dict_id_make (field, strlen (field)); 
    
    char dump_fn[strlen(field) + strlen (filename_ext) + 2];
    sprintf (dump_fn, "%s.%s", field, filename_ext);
    
    dump_file = fopen (dump_fn, "wb"); // it will be closed implicitly when the process terminates

    fprintf (stderr, "Writing %s data to %s\n", field, dump_fn);

    mutex_initialize (dump_mutex);
}

void mtf_dump_local (ContextP ctx, bool local /* true = local, false = b250 */)
{
    mutex_lock (dump_mutex);

    if (local)
        fwrite (ctx->local.data, 1, ctx->local.len * lt_sizeof_one[ctx->ltype], dump_file);
    else
        fwrite (ctx->b250.data, 1, ctx->b250.len, dump_file);

    mutex_unlock (dump_mutex);
}
