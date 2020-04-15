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
#include "gloptimize_vcf.h"
#include "base250.h"
#include "vblock.h"
#include "move_to_front.h"
#include "zfile.h"
#include "endianness.h"
#include "file.h"
#include "hash.h"

#define INITIAL_NUM_NODES 10000

#define CTX_GROWTH 1.75

static pthread_mutex_t wait_for_vb_1_mutex;
static pthread_mutex_t compress_dictionary_data_mutex;

static inline void mtf_lock_do (VBlock *vb, pthread_mutex_t *mutex, const char *func, uint32_t code_line, const char *name, uint32_t param)
{
    //printf ("thread %u vb_i=%u LOCKING %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
    pthread_mutex_lock (mutex);
    //printf ("thread %u vb_i=%u LOCKED %s:%u from %s:%u\n", (unsigned)pthread_self(), vb->vblock_i, name, param, func, code_line);
}
#define mtf_lock(vb, mutex, name, param) mtf_lock_do (vb, mutex, __FUNCTION__, __LINE__, name, param)

static inline void mtf_unlock_do (VBlock *vb, pthread_mutex_t *mutex, const char *func, uint32_t code_line, const char *name, uint32_t param)
{
    pthread_mutex_unlock (mutex);
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
static inline uint32_t mtf_insert_to_dict (VBlock *vb_of_dict, MtfContext *ctx, bool is_zf_ctx, const char *snip, uint32_t snip_len)
{
    buf_alloc (vb_of_dict, &ctx->dict, MAX ((ctx->dict.len + snip_len + 1), INITIAL_NUM_NODES * MIN (10, snip_len)), 
               CTX_GROWTH, is_zf_ctx ? "z_file->mtf_ctx->dict" : "mtf_ctx->dict", ctx->did_i);
    
    if (is_zf_ctx) buf_set_overlayable (&ctx->dict); // during merge
    
    unsigned char_index = ctx->dict.len;
    char *dict_p = &ctx->dict.data[char_index];

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = '\t'; // dictionary snips have a \t separator within dictionary string

    ctx->dict.len += snip_len + 1;
    return char_index;
}

// ZIP only (PIZ doesn't have mtf) mtf index to node - possibly in ol_mtf, or in mtf
MtfNode *mtf_node_do (const MtfContext *ctx, uint32_t mtf_i, 
                      const char **snip_in_dict, uint32_t *snip_len,  // optional outs
                      const char *func, uint32_t code_line)
{
    ASSERT0 (ctx->dict_id.num, "Error: this ctx is not initialized");
    
    ASSERT (mtf_i < ctx->mtf.len + ctx->ol_mtf.len, "Error in mtf_node: out of range: dict=%.*s %s mtf_i=%d mtf.len=%u ol_mtf.len=%u. Caller: %s:%u",  
            DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, st_name (ctx->dict_section_type),
            mtf_i, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len, func, code_line);

    bool is_ol = mtf_i < ctx->ol_mtf.len; // is this entry from a previous vb (overlay buffer)

    MtfNode *node = is_ol ? ENT (MtfNode, ctx->ol_mtf, mtf_i)
                          : ENT (MtfNode, ctx->mtf, mtf_i - ctx->ol_mtf.len);

    if (snip_in_dict) {
        const Buffer *dict = is_ol ? &ctx->ol_dict : &ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error: dict not allocated");

        *snip_in_dict = &dict->data[node->char_index];
    }

    if (snip_len) *snip_len = node->snip_len;
    
    return node;
}

// PIZ: search for a node matching this snip in a directory and return the node index. note that we do a linear
// search as PIZ doesn't have hash tables.
int32_t mtf_search_for_node_index (MtfContext *ctx, const char *snip, unsigned snip_len)
{
    MtfWord *words = (MtfWord *)ctx->word_list.data;

    for (unsigned i=0; i < ctx->word_list.len; i++)
        if (words[i].snip_len == snip_len && !memcmp (&ctx->dict.data[words[i].char_index], snip, snip_len))
            return i;

    return NIL;
}

// PIZ only (uses word_list): called when pizzing a genotype section - returns snip and index, and advances the index
uint32_t mtf_get_next_snip (VBlock *vb, MtfContext *ctx, 
                            SnipIterator *override_iterator,   // if NULL, taken from ctx
                            const char **snip, uint32_t *snip_len, // optional out
                            uint32_t vcf_line) 
{
    SnipIterator *iterator = override_iterator ? override_iterator : &ctx->iterator;
    
    if (!override_iterator && !iterator->next_b250) // INFO and Field1-9 data (GT data uses override_next_b250)
        iterator->next_b250 = (uint8_t *)ctx->b250.data; // initialize (GT data initializes to the beginning of each sample rather than the beginning of the data)

    uint32_t word_index = z_file->genozip_version >= 2 ? base250_decode (&iterator->next_b250)  // if this line has no non-GT subfields, it will not have a ctx 
                                                       : v1_base250_decode (&iterator->next_b250);

    // case: a subfield snip is missing - either the genotype data has less subfields than declared in FORMAT, or not provided at all for some (or all) samples.
    if (word_index == WORD_INDEX_MISSING_SF) {
        ASSERT (!ctx || ctx->b250_section_type == SEC_VCF_GT_DATA, "Error while reconstrucing line %u vb_i=%u: BASE250_MISSING_SF unexpectedly found in b250 data of %.*s (%s)",
                vcf_line, vb->vblock_i, DICT_ID_LEN, dict_id_printable(ctx->dict_id).id, st_name(ctx->b250_section_type)); // there will be no context if this GT subfield was always missing - never appeared on any sample

        if (snip) {
            *snip = NULL; // ignore this dict_id - don't even output a separator
            *snip_len = 0;
        }
    }

    // case: a subfield snip is empty, eg AB::CD
    else if (word_index == WORD_INDEX_EMPTY_SF) { 
        ASSERT (ctx->b250_section_type == SEC_VCF_GT_DATA, "Error while reconstrucing line %u: BASE250_EMPTY_SF unexpectedly found in b250 data of %.*s",
                vcf_line, DICT_ID_LEN, dict_id_printable(ctx->dict_id).id);
        if (snip) {
            *snip = ""; // pointer to static empty string
            *snip_len = 0;
        }
    }

    else {
        if (word_index == WORD_INDEX_ONE_UP) 
            word_index = ctx->iterator.prev_word_index + 1;

        ASSERT (word_index < ctx->word_list.len, "Error while parsing line %u: word_index=%u is out of bounds - %s%s \"%.*s\" dictionary has only %u entries",
                vcf_line, word_index, 
                ctx->dict_section_type == SEC_VCF_INFO_SF_DICT ? "INFO" : "",
                ctx->dict_section_type == SEC_VCF_FRMT_SF_DICT ? "FORMAT" : "",
                DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, (uint32_t)ctx->word_list.len);

        //MtfWord *dict_word = &((MtfWord*)ctx->word_list.data)[word_index];
        MtfWord *dict_word = mtf_get_word (ctx, word_index);

        if (snip) {
            *snip = &ctx->dict.data[dict_word->char_index];
            *snip_len = dict_word->snip_len;
        }
    }
    
    iterator->prev_word_index = word_index;    

    return word_index;
}

// Process and snip - return its node index, and enter it into the directory if its not already there. Called
// 1. During segregate - as snips are encountered in the data. No base250 encoding yet
// 2. During mtf_merge_in_vb_ctx_one_dict_id() - to enter snips into z_file->mtf_ctx - also encoding in base250
static uint32_t mtf_evaluate_snip_merge (VBlock *merging_vb, MtfContext *zf_ctx, MtfContext *vb_ctx, 
                                         const char *snip, uint32_t snip_len,
                                         MtfNode **node, bool *is_new)  // out
{
    // attempt to get the node from the hash table
    uint32_t new_mtf_i_if_no_old_one = zf_ctx->mtf.len;

    int32_t existing_mtf_i = hash_get_entry_for_merge (zf_ctx, snip, snip_len, new_mtf_i_if_no_old_one, node);
    if (existing_mtf_i != NIL) {
        *is_new = false;
        return existing_mtf_i; // snip was found in hash table - we're done
    }
    
    // this snip was just added to the hash table - its a new snip

    zf_ctx->mtf.len++; // we have a new dictionary item - this snip

    ASSERT (zf_ctx->mtf.len < MAX_WORDS_IN_CTX, 
            "Error: too many words in directory %.*s, max allowed number of words is is %u", DICT_ID_LEN, 
            dict_id_printable (zf_ctx->dict_id).id, MAX_WORDS_IN_CTX);

    buf_alloc (evb, &zf_ctx->mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, zf_ctx->mtf.len), CTX_GROWTH, 
               "z_file->mtf_ctx->mtf", zf_ctx->did_i);
    buf_set_overlayable (&zf_ctx->mtf);

    *node = mtf_node (zf_ctx, new_mtf_i_if_no_old_one, NULL, NULL);
    memset (*node, 0, sizeof(MtfNode)); // safety
    (*node)->snip_len     = snip_len;
    (*node)->char_index   = mtf_insert_to_dict (evb, zf_ctx, true, snip, snip_len);
    (*node)->word_index.n = new_mtf_i_if_no_old_one;

    *is_new = true;

    return new_mtf_i_if_no_old_one;
}

uint32_t mtf_evaluate_snip_seg (VBlock *segging_vb, MtfContext *vb_ctx, 
                                const char *snip, uint32_t snip_len,
                                MtfNode **node /* out */, bool *is_new /* out */)
{
    ASSERT0 (vb_ctx, "Error in mtf_evaluate_snip_seg: vb_ctx is NULL");

    segging_vb->z_section_entries[vb_ctx->b250_section_type]++; 

    if (!snip_len) 
        return (!snip || (z_file->data_type == DATA_TYPE_VCF && *snip != ':')) ? WORD_INDEX_MISSING_SF : WORD_INDEX_EMPTY_SF;

    uint32_t new_mtf_i_if_no_old_one = vb_ctx->ol_mtf.len + vb_ctx->mtf.len;
    
    ASSERT (new_mtf_i_if_no_old_one <= MAX_WORDS_IN_CTX, 
            "Error: ctx of %.*s is full (max allowed words=%u): ol_mtf.len=%u mtf.len=%u",
            DICT_ID_LEN, dict_id_printable(vb_ctx->dict_id).id, MAX_WORDS_IN_CTX, (uint32_t)vb_ctx->ol_mtf.len, (uint32_t)vb_ctx->mtf.len)

    // get the node from the hash table if it already exists, or add this snip to the hash table if not
    int32_t existing_mtf_i = hash_get_entry_for_seg (segging_vb, vb_ctx,snip, snip_len, new_mtf_i_if_no_old_one, node);
    if (existing_mtf_i != NIL) {

        if (segging_vb->vblock_i == 1) 
            ENT (SorterEnt, vb_ctx->sorter, existing_mtf_i)->count++;

        if (is_new) *is_new = false;
        return existing_mtf_i; // snip found - we're done
    }
    
    // this snip isn't in the hash table - its a new snip
    ASSERT (vb_ctx->mtf.len < 0x7fffffff, "Error: too many words in directory %.*s", DICT_ID_LEN, dict_id_printable (vb_ctx->dict_id).id);

    segging_vb->z_section_entries[vb_ctx->dict_section_type]++; 

    buf_alloc (segging_vb, &vb_ctx->mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, 1+vb_ctx->mtf.len), CTX_GROWTH, 
               "mtf_ctx->mtf", vb_ctx->did_i);

    vb_ctx->mtf.len++; // new hash entry or extend linked list

    *node = mtf_node (vb_ctx, new_mtf_i_if_no_old_one, NULL, NULL);
    memset (*node, 0, sizeof(MtfNode)); // safety
    (*node)->snip_len     = snip_len;
    (*node)->char_index   = mtf_insert_to_dict (segging_vb, vb_ctx, false, snip, snip_len);
    (*node)->word_index.n = new_mtf_i_if_no_old_one;

    // if this is the first variant block - allocate/grow sorter to contain exactly the same number of entries as mtf
    if (segging_vb->vblock_i == 1) {
        unsigned prev_size = vb_ctx->sorter.size;
        buf_alloc (segging_vb, &vb_ctx->sorter, sizeof (SorterEnt) * (vb_ctx->mtf.size / sizeof(MtfNode)), 1, "mtf_ctx->sorter", 0);
        if (vb_ctx->sorter.size > prev_size) memset (&vb_ctx->sorter.data[prev_size], 0, vb_ctx->sorter.size - prev_size);

        SorterEnt *sorter_ent = ENT (SorterEnt, vb_ctx->sorter, new_mtf_i_if_no_old_one);
        sorter_ent->mtf_i = new_mtf_i_if_no_old_one;
        sorter_ent->count = 1;
    }

    if (is_new) *is_new = true;
    return new_mtf_i_if_no_old_one;
}

static void mtf_init_mapper (VBlock *vb, VcfFields field_i, Buffer *mapper_buf, const char *name)
{
    if (!buf_is_allocated (&vb->mtf_ctx[field_i].ol_mtf)) return;
        
    mapper_buf->len = vb->mtf_ctx[field_i].ol_mtf.len;
    
    buf_alloc (vb, mapper_buf,  mapper_buf->len * sizeof (SubfieldMapper), 2, name, 0);
    
    for (unsigned i=0; i < mapper_buf->len; i++) 
        ((SubfieldMapper *)mapper_buf->data)[i].num_subfields = (uint8_t)NIL;
}

// ZIP only: overlay and/or copy the current state of the global context to the vb, ahead of compressing this vb.
void mtf_clone_ctx (VBlock *vb)
{
    unsigned z_num_dict_ids = __atomic_load_n (&z_file->num_dict_ids, __ATOMIC_RELAXED);

    START_TIMER; // including mutex wait time

    for (unsigned did_i=0; did_i < z_num_dict_ids; did_i++) {
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];
        MtfContext *zf_ctx = &z_file->mtf_ctx[did_i];

        ASSERT (zf_ctx->mutex_initialized, "Error: expected zf_ctx->mutex_initialized for did_i=%u", did_i);
        mtf_lock (vb, &zf_ctx->mutex, "zf_ctx", did_i);

        if (buf_is_allocated (&zf_ctx->dict)) {  // something already for this dict_id

            // overlay the global dict and mtf - these will not change by this (or any other) VB
            //fprintf (stderr,  ("mtf_clone_ctx: overlaying old dict %.8s, to vb_i=%u vb_did_i=z_did_i=%u\n", dict_id_printable (zf_ctx->dict_id).id, vb->vblock_i, did_i);
            buf_overlay (vb, &vb_ctx->ol_dict, &zf_ctx->dict, "ctx->ol_dict", did_i);   
            buf_overlay (vb, &vb_ctx->ol_mtf, &zf_ctx->mtf, "ctx->ol_mtf", did_i);   

            // overlay the hash table, that may still change by future vb's merging... this vb will only use
            // entries that are up to this merge_num
            buf_overlay (vb, &vb_ctx->global_hash, &zf_ctx->global_hash, "mtf_ctx->global_hash", did_i);
            vb_ctx->merge_num = zf_ctx->merge_num;
            vb_ctx->global_hash_prime = zf_ctx->global_hash_prime; // can never change
            vb_ctx->num_new_entries_prev_merged_vb = zf_ctx->num_new_entries_prev_merged_vb;
        }

        vb_ctx->did_i             = did_i;
        vb_ctx->dict_id           = zf_ctx->dict_id;
        vb_ctx->dict_section_type = zf_ctx->dict_section_type;
        vb_ctx->b250_section_type = zf_ctx->b250_section_type;
        mtf_init_iterator (vb_ctx);

        mtf_unlock (vb, &zf_ctx->mutex, "zf_ctx", did_i);
    }

    vb->num_dict_ids = z_num_dict_ids;

    // For VCF: initialize mappers for VCF_FORMAT and VCF_INFO
    if (vb->data_type == DATA_TYPE_VCF) {    
        mtf_init_mapper (vb, VCF_FORMAT, &((VBlockVCF *)vb)->format_mapper_buf, "format_mapper_buf");    
        mtf_init_mapper (vb, VCF_INFO,   &((VBlockVCF *)vb)->iname_mapper_buf, "iname_mapper_buf");    
    }

    COPY_TIMER (vb->profile.mtf_clone_ctx)
}

void mtf_initialize_mutex (void)
{
    if (z_file->dicts_mutex_initialized) return;

    unsigned ret = pthread_mutex_init (&z_file->dicts_mutex, NULL);
    ASSERT0 (!ret, "pthread_mutex_init failed for z_file->dicts_mutex");

    ret = pthread_mutex_init (&wait_for_vb_1_mutex, NULL);
    ASSERT0 (!ret, "pthread_mutex_init failed for wait_for_vb_1_mutex");

    ret = pthread_mutex_init (&compress_dictionary_data_mutex, NULL);
    ASSERT0 (!ret, "pthread_mutex_init failed for compress_dictionary_data_mutex");
}

// find the z_file context that corresponds to dict_id. It could be possibly a different did_i
// than in the vb - in case this dict_id is new to this vb, but another vb already inserted
// it to z_file
static MtfContext *mtf_get_zf_ctx (DictIdType dict_id)
{
    unsigned z_num_dict_ids = __atomic_load_n (&z_file->num_dict_ids, __ATOMIC_RELAXED);

    for (unsigned did_i=0; did_i < z_num_dict_ids; did_i++)
        if (dict_id.num == z_file->mtf_ctx[did_i].dict_id.num) 
            return &z_file->mtf_ctx[did_i];

    return NULL;
}

// ZIP only: called by merging VBs to add a new dict to z_file - copying some stuff from vb_ctx
static MtfContext *mtf_add_new_zf_ctx (VBlock *merging_vb, const MtfContext *vb_ctx)
{
    // adding a new dictionary is proctected by a mutex. note that z_file->num_dict_ids is accessed by other threads
    // without mutex proction when searching for a dictionary - that's why we update it at the end, after the new
    // zf_ctx is set up with the new dict_id (ready for another thread to search it)
    mtf_lock (merging_vb, &z_file->dicts_mutex, "dicts_mutex", 0);

    // check if another thread raced and created this dict before us
    MtfContext *zf_ctx = mtf_get_zf_ctx (vb_ctx->dict_id);
    if (zf_ctx) goto finish;

    ASSERT (z_file->num_dict_ids+1 < MAX_DICTS, // load num_dict_ids - this time with mutex protection - it could have changed
            "Error: z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    zf_ctx = &z_file->mtf_ctx[z_file->num_dict_ids];

    ASSERT (!pthread_mutex_init (&zf_ctx->mutex, NULL), 
            "pthread_mutex_init failed for zf_ctx->mutex did_i=%u", zf_ctx->did_i);

    zf_ctx->mutex_initialized = true;
    zf_ctx->did_i             = z_file->num_dict_ids; 
    zf_ctx->dict_id           = vb_ctx->dict_id;
    zf_ctx->b250_section_type = vb_ctx->b250_section_type;
    zf_ctx->dict_section_type = vb_ctx->dict_section_type;

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
    MtfContext *vb_ctx = &merging_vb->mtf_ctx[did_i];

    // get the ctx or create a new one. note: mtf_add_new_zf_ctx() must be called before mtf_lock() because it locks the z_file mutex (avoid a deadlock)
    MtfContext *zf_ctx = mtf_get_zf_ctx (vb_ctx->dict_id);
    if (!zf_ctx) zf_ctx = mtf_add_new_zf_ctx (merging_vb, vb_ctx); 

    { START_TIMER; 
      mtf_lock (merging_vb, &zf_ctx->mutex, "zf_ctx", zf_ctx->did_i);
      COPY_TIMER (merging_vb->profile.lock_mutex_zf_ctx);  
    }

    START_TIMER; // note: careful not to count time spent waiting for the mutex
    //fprintf (stderr,  ("Merging dict_id=%.8s into z_file vb_i=%u vb_did_i=%u z_did_i=%u\n", dict_id_printable (vb_ctx->dict_id).id, merging_vb->vblock_i, did_i, z_did_i);

    zf_ctx->merge_num++; // first merge is #1 (first clone which happens before the first merge, will get vb-)

    if (!buf_is_allocated (&vb_ctx->dict)) goto finish; // nothing yet for this dict_id

    uint32_t start_dict_len = zf_ctx->dict.len;
    uint32_t start_mtf_len  = zf_ctx->mtf.len;
    
    if (!buf_is_allocated (&zf_ctx->dict)) {

        // first data for this dict (usually, but not always, vb_i=1) - move to zf_ctx and leave overlay
        zf_ctx->num_new_entries_prev_merged_vb = vb_ctx->mtf.len; // number of new words in this dict from this VB

        buf_move (evb, &zf_ctx->dict, merging_vb, &vb_ctx->dict);
        buf_set_overlayable (&zf_ctx->dict);
        buf_overlay (merging_vb, &vb_ctx->ol_dict, &zf_ctx->dict, "ctx->ol_dict", did_i);

        buf_move (evb, &zf_ctx->mtf,  merging_vb, &vb_ctx->mtf);
        buf_set_overlayable (&zf_ctx->mtf);
        buf_overlay (merging_vb, &vb_ctx->ol_mtf, &zf_ctx->mtf, "ctx->ol_mtf", did_i);

        // allocate hash table, based on the statitics gather by this first vb that is merging this dict and 
        // populate the hash table without needing to reevalate the snips (we know none are in the hash table, but all are in mtf and dict)
        if (!buf_is_allocated (&zf_ctx->global_hash)) 
            hash_alloc_global (merging_vb, zf_ctx, vb_ctx);

        // encode in base250 - to be used by zip_vcf_generate_genotype_one_section() and zip_generate_b250_section()
        for (unsigned i=0; i < zf_ctx->mtf.len; i++) {
            MtfNode *zf_node = &((MtfNode *)zf_ctx->mtf.data)[i];
            zf_node->word_index = base250_encode (zf_node->word_index.n); // note that vb overlays this. also, vb_1 has been sorted so word_index != node_index

            ASSERT (zf_node->word_index.n < zf_ctx->mtf.len, // sanity check
                    "Error: word_index=%u out of bound - mtf.len=%u, in dictionary %.*s", 
                    (uint32_t)zf_node->word_index.n, (uint32_t)zf_ctx->mtf.len, DICT_ID_LEN, dict_id_printable (zf_ctx->dict_id).id);
        }
    }
    else {
        // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
        // (vb_ctx->mtf contains only new words, old words from previous vbs are in vb_ctx->ol_mtf)
        for (unsigned i=0; i < vb_ctx->mtf.len; i++) {
            MtfNode *vb_node = &((MtfNode *)vb_ctx->mtf.data)[i];
            
            MtfNode *zf_node;
            bool is_new;
            // use evb and not vb because zf_context is z_file (which belongs to evb)
            int32_t zf_node_index = mtf_evaluate_snip_merge (merging_vb, zf_ctx, vb_ctx, 
                                                             &vb_ctx->dict.data[vb_node->char_index], vb_node->snip_len,
                                                             &zf_node, &is_new);

            ASSERT (zf_node_index >= 0 && zf_node_index < zf_ctx->mtf.len, "Error: zf_node_index=%d out of range - len=%i", zf_node_index, (uint32_t)vb_ctx->mtf.len);

            // set word_index to be indexing the global dict - to be used by zip_vcf_generate_genotype_one_section() and zip_generate_b250_section()
            if (is_new)
                vb_node->word_index = zf_node->word_index = base250_encode (zf_node_index);
            else 
                // a previous VB already already calculated the word index for this node. if it was done by vb_i=1,
                // then it is also re-sorted and the word_index is no longer the same as the node_index
                vb_node->word_index = zf_node->word_index;
        }
    }

    // we now compress the dictionaries directly from z_file. note: we must continue to hold
    // the mutex during compression, lest another thread re-alloc the dictionary.
    const char *start_dict = &zf_ctx->dict.data[start_dict_len]; // we take the pointer AFTER the evaluate, since dict can be reallocted
    unsigned added_chars   = zf_ctx->dict.len - start_dict_len;
    unsigned added_words   = zf_ctx->mtf.len  - start_mtf_len;

    // compress incremental part of dictionary added by this vb. note: dispatcher calls this function in the correct order of VBs.
    if (added_chars) {
        // special optimization for the GL dictionary (it is ineffective with --optimize that already optimizes GL)
        if (zf_ctx->dict_id.num == dict_id_FORMAT_GL && !flag_optimize) 
            start_dict = gl_optimize_dictionary ((VBlockVCFP)merging_vb, &zf_ctx->dict, &((MtfNode *)zf_ctx->mtf.data)[start_mtf_len], start_dict_len, added_words);

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

// ZIP only: merge new words added in this vb into the z_file.mtf_ctx, and compresses dictionaries.
void mtf_merge_in_vb_ctx (VBlock *merging_vb)
{
    // vb_i=1 goes first, as it has the sorted dictionaries, other vbs can go in
    // arbitrary order. at the end of this function, vb_i releases the mutex it locked along time ago,
    // while the other vbs wait for vb_1 by attempting to lock the mutex
    if (merging_vb->vblock_i != 1) {
        mtf_lock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", merging_vb->vblock_i);
        mtf_unlock (merging_vb, &wait_for_vb_1_mutex, "wait_for_vb_1_mutex", merging_vb->vblock_i);
    }
    
    // first, all field dictionaries    
    for (unsigned did_i=0; did_i < merging_vb->num_dict_ids; did_i++) {
        if (!buf_is_allocated (&merging_vb->mtf_ctx[did_i].dict)) continue;

        SectionType dict_sec_type = merging_vb->mtf_ctx[did_i].dict_section_type;

        ASSERT (section_type_is_dictionary(dict_sec_type), "Error: dict_sec_type=%s is not a dictionary section", st_name(dict_sec_type));

        if (dict_sec_type != SEC_VCF_INFO_SF_DICT && dict_sec_type != SEC_VCF_FRMT_SF_DICT) 
            mtf_merge_in_vb_ctx_one_dict_id (merging_vb, did_i);
    }

    // second, all the info subfield dictionaries
    for (unsigned did_i=0; did_i < merging_vb->num_dict_ids; did_i++)         
        if (buf_is_allocated (&merging_vb->mtf_ctx[did_i].dict) && 
            merging_vb->mtf_ctx[did_i].dict_section_type == SEC_VCF_INFO_SF_DICT) 
            mtf_merge_in_vb_ctx_one_dict_id (merging_vb, did_i);

    // third, all the genotype subfield dictionaries
    for (unsigned did_i=0; did_i < merging_vb->num_dict_ids; did_i++)         
        if (buf_is_allocated (&merging_vb->mtf_ctx[did_i].dict) && 
            merging_vb->mtf_ctx[did_i].dict_section_type == SEC_VCF_FRMT_SF_DICT) 
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
// or returns NIL, if not
uint8_t mtf_get_existing_did_i_by_dict_id (VBlock *vb, DictIdType dict_id)
{
    //MtfContext *mtf_ctx = z_file->mtf_ctx;

    for (uint8_t did_i=0; did_i < z_file->num_dict_ids; did_i++) 
        if (dict_id.num == vb->mtf_ctx[did_i].dict_id.num) return did_i;

    return DID_I_NONE; // not found
}

// gets did_id if the dictionary exists, and creates a new dictionary if its the first time dict_id is encountered
// threads: no issues - called by PIZ for vb and zf (but dictionaries are immutable) and by Segregate (ZIP) on vb_ctx only
MtfContext *mtf_get_ctx_by_dict_id (MtfContext *mtf_ctx /* an array */, 
                                    unsigned *num_dict_ids, 
                                    uint8_t *num_subfields, // variable to increment if a new context is added
                                    DictIdType dict_id,
                                    SectionType dict_section_type)
{
    int did_i = dict_id_get_field (dict_id); // will return the field if its a main field dict (even if the dictinary doesn't exit yet), or -1 if not

    // case: not a main field dict - find the did_i if we have it already
    if (did_i < 0) { // did_i is minus the next field
        did_i=-did_i ; for (; did_i < *num_dict_ids; did_i++) 
            if (dict_id.num == mtf_ctx[did_i].dict_id.num) break;
    }

    MtfContext *ctx = &mtf_ctx[did_i]; // existing or new

    // case: dict_id encountered for this first time - initialize a mtf_ctx
    if (!ctx->dict_id.num) {

        //fprintf (stderr,  ("Inserting new vb dict_id=%.8s in did_i=num_dict_ids=%u \n", dict_id_printable (dict_id).id, did_i);
        ASSERT (*num_dict_ids+1 < MAX_DICTS, 
                "Error: number of dictionary types is greater than MAX_DICTS=%u", MAX_DICTS);

        ctx->did_i             = did_i;
        ctx->dict_id           = dict_id;
        ctx->dict_section_type = dict_section_type;
        ctx->b250_section_type = dict_section_type + 1; // the b250 is 1 after the dictionary for all dictionary sections
        mtf_init_iterator (ctx);

        // thread safety: the increment below MUST be AFTER memcpy, bc piz_get_line_subfields
        // might be reading this data at the same time as the piz dispatcher thread adding more dictionaries
        (*num_dict_ids) = MAX (datatype_last_field[z_file->data_type], did_i) + 1; 

        if (num_subfields) (*num_subfields)++;
    }

    ASSERT (ctx->dict_section_type == dict_section_type, "Error: mismatch in dict_id=%.*s dict_section_type: requested %s but in the ctx says: %s",
            DICT_ID_LEN, dict_id_printable (dict_id).id, st_name(dict_section_type), st_name(ctx->dict_section_type));

    return ctx;
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
    // header_genozip_to_txt() so that they don't need to be realloced. dict.len may be accessed
    // by compute threads, but its change is assumed to be atomic, so that no weird things will happen
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)section_data;

    ASSERT (section_type_is_dictionary(header->h.section_type),
            "Error: header->h.section_type=%s is not a dictionary section", st_name(header->h.section_type));

    uint32_t num_snips = BGEN32 (header->num_snips);

    zfile_uncompress_section (vb, section_data, &fragment, "fragment", header->h.section_type);

    // special treatment if this is GL - de-optimize
    if (header->dict_id.num == dict_id_FORMAT_GL)
        gl_deoptimize_dictionary (fragment.data, fragment.len);

    // in piz, the same did_i is used for z_file and vb contexts, meaning that in vbs there could be
    // a non-contiguous array of contexts (some are missing if not used by this vb)

    MtfContext *zf_ctx = mtf_get_ctx_by_dict_id (z_file->mtf_ctx, &z_file->num_dict_ids, NULL, header->dict_id, header->h.section_type);
    
    // append fragment to dict. If there is no room - old memory is abandoned (so that VBs that are overlaying
    // it continue to work uninterrupted) and a new memory is allocated, where the old dict is joined by the new fragment
    unsigned dict_old_len = zf_ctx->dict.len;
    buf_alloc (evb, &zf_ctx->dict, zf_ctx->dict.len + fragment.len, CTX_GROWTH, "z_file->mtf_ctx->dict", header->h.section_type);
    buf_set_overlayable (&zf_ctx->dict);

    memcpy (&zf_ctx->dict.data[zf_ctx->dict.len], fragment.data, fragment.len);
    zf_ctx->dict.len += fragment.len;

    // extend word list memory - and calculate the new words. If there is no room - old memory is abandoned 
    // (so that VBs that are overlaying it continue to work uninterrupted) and a new memory is allocated
    buf_alloc (evb, &zf_ctx->word_list, (zf_ctx->word_list.len + num_snips) * sizeof (MtfWord), CTX_GROWTH, "z_file->mtf_ctx->word_list", header->h.section_type);
    buf_set_overlayable (&zf_ctx->word_list);

    bool is_ref_alt = !strncmp ((char*)dict_id_printable (header->dict_id).id, vcf_field_names[VCF_REFALT], MIN (strlen(vcf_field_names[VCF_REFALT]+1), DICT_ID_LEN)); // compare inc. \0 terminator

    char *start = fragment.data;
    for (unsigned snip_i=0; snip_i < num_snips; snip_i++) {

        MtfWord *word = mtf_get_word (zf_ctx, zf_ctx->word_list.len++);

        char *c=start; while (*c != '\t') c++;

        // special case of REFALT - there is always one \t in the middle of the snip, eg "A\tC"
        if (is_ref_alt) {
            c++;
            while (*c != '\t') c++;
        }

        word->snip_len   = c - start;
        word->char_index = dict_old_len + (start - fragment.data);

        start = c+1; // skip over the \t
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
    for (unsigned did_i=0; did_i < MAX_DICTS; did_i++) {
        MtfContext *zf_ctx = &z_file->mtf_ctx[did_i];
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];

        if (!zf_ctx->dict_id.num) continue;

        if (buf_is_allocated (&zf_ctx->dict) && buf_is_allocated (&zf_ctx->word_list)) { 
            
            vb_ctx->did_i             = did_i;
            vb_ctx->dict_id           = zf_ctx->dict_id;
            vb_ctx->b250_section_type = zf_ctx->b250_section_type;
            vb_ctx->dict_section_type = zf_ctx->dict_section_type;
            mtf_init_iterator (vb_ctx);

            buf_overlay (vb, &vb_ctx->dict, &zf_ctx->dict, "ctx->dict", did_i);    
            buf_overlay (vb, &vb_ctx->word_list, &zf_ctx->word_list, "ctx->word_list", did_i);

            // count dictionaries of genotype data subfields
            // note: this is needed only for V1 files...
            if (vb->data_type == DATA_TYPE_VCF && dict_id_is_vcf_format_sf (vb_ctx->dict_id)) {

                ASSERT (++((VBlockVCFP)vb)->num_format_subfields <= MAX_SUBFIELDS, 
                        "Error: number of subfields in %s exceeds MAX_SUBFIELDS=%u, while reading vb_i=%u", 
                        z_name, MAX_SUBFIELDS, vb->vblock_i);
            }
        }
    }
    vb->num_dict_ids = z_file->num_dict_ids;
}

static int sorter_cmp(const void *a_, const void *b_)  
{ 
    SorterEnt *a = (SorterEnt *)a_;
    SorterEnt *b = (SorterEnt *)b_;
    
    return (int)b->count - (int)a->count;
}

void mtf_sort_dictionaries_vb_1(VBlock *vb)
{
    // thread safety note: no issues here, as this is run only by the compute thread of vblock_i=1
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {

        MtfContext *ctx = &vb->mtf_ctx[did_i];

        // sort in ascending order of mtf->count
        qsort (ctx->sorter.data, ctx->mtf.len, sizeof(SorterEnt), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in mtf
        static Buffer old_dict = EMPTY_BUFFER;
        buf_move (vb, &old_dict, vb, &ctx->dict);

        buf_alloc (vb, &ctx->dict, old_dict.size, 1, "mtf_ctx->dict", did_i);
        ctx->dict.len = old_dict.len;

        char *next = ctx->dict.data;
        for (unsigned i=0; i < ctx->mtf.len; i++) {
            int32_t mtf_i = ((SorterEnt *)ctx->sorter.data)[i].mtf_i;
            MtfNode *node = &((MtfNode *)ctx->mtf.data)[mtf_i];
            memcpy (next, &old_dict.data[node->char_index], node->snip_len + 1 /* +1 for \t */);
            node->char_index   = next - ctx->dict.data;
            node->word_index.n = i;

            next += node->snip_len + 1;
        }

        buf_destroy (&old_dict);
    }
}

// zero all sorters - this is called in case of a re-do of the first VB due to ploidy overflow
void mtf_zero_all_sorters (VBlock *vb)
{
    ASSERT (vb->vblock_i == 1, "Error in mtf_zero_all_sorters: expected vb_i==1, but vb_i==%u", vb->vblock_i);
    
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {
        MtfContext *ctx = &vb->mtf_ctx[did_i];
        buf_zero (&ctx->sorter);
    }
}

// ZIP only: run by I/O thread during zip_output_processed_vb()
void mtf_update_stats (VBlock *vb)
{
    // zf_ctx doesn't store mtf_i, but we just use mtf_i.len as a counter for displaying in genozip_show_sections
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];
    
        MtfContext *zf_ctx = mtf_get_zf_ctx (vb_ctx->dict_id);
        if (!zf_ctx) continue; // this can happen if FORMAT subfield appears, but no line has data for it

        zf_ctx->mtf_i.len += vb_ctx->mtf_i.len; // thread safety: no issues, this only updated only by the I/O thread
    }
}

void mtf_free_context (MtfContext *ctx)
{
    buf_free (&ctx->ol_dict);
    buf_free (&ctx->ol_mtf);
    buf_free (&ctx->dict);
    buf_free (&ctx->mtf);
    buf_free (&ctx->word_list);
    buf_free (&ctx->local_hash);
    buf_free (&ctx->global_hash);
    buf_free (&ctx->sorter);
    buf_free (&ctx->mtf_i);
    buf_free (&ctx->b250);
    ctx->dict_id.num = 0;
    ctx->dict_section_type = ctx->b250_section_type = 0;
    ctx->iterator.next_b250 = NULL;
    ctx->iterator.prev_word_index =0;
    ctx->local_hash_prime = 0;
    ctx->global_hash_prime = 0;
    ctx->merge_num = 0;
    ctx->num_lines_at_half = 0;
    ctx->mtf_len_at_half = 0;

    if (ctx->mutex_initialized) {
        pthread_mutex_destroy (&ctx->mutex);
        ctx->mutex_initialized = false;
    }
}

// Called by file_close ahead of freeing File memory containing mtf_ctx
void mtf_destroy_context (MtfContext *ctx)
{
    buf_destroy (&ctx->ol_dict);
    buf_destroy (&ctx->ol_mtf);
    buf_destroy (&ctx->dict);
    buf_destroy (&ctx->mtf);
    buf_destroy (&ctx->word_list);
    buf_destroy (&ctx->local_hash);
    buf_destroy (&ctx->global_hash);
    buf_destroy (&ctx->sorter);
    buf_destroy (&ctx->mtf_i);
    buf_destroy (&ctx->b250);

    if (ctx->mutex_initialized) {
        pthread_mutex_destroy (&ctx->mutex);
        ctx->mutex_initialized = false;
    }
}
