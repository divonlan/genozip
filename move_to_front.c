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
#include "gloptimize.h"
#include "base250.h"
#include "vb.h"
#include "move_to_front.h"
#include "zfile.h"
#include "endianness.h"
#include "file.h"

#define INITIAL_NUM_NODES 10000

// tested hash table sizes up to 5M. turns out smaller tables (up to a point) are faster, despite having longer
// average linked lists. probably bc the CPU can store the entire hash and mtf arrays in L1 or L2
// memory cache during segmentation
#define MTF_HASH_TABLE_LEN 19997         // a prime number, see: https://primes.utm.edu/lists/small/millions/
static inline int32_t mtf_hash (const MtfContext *ctx, const char *snip, unsigned snip_len)
{
    uint64_t result=0;
    for (unsigned i=0; i < snip_len; i++) 
        result = ((result << 23) | (result >> 41)) ^ (uint64_t)((uint8_t)snip[i]);

    return (uint32_t)(result % MTF_HASH_TABLE_LEN);
}

// add a snip to the dictionary the first time it is encountered in the VCF file.
// the dictionary will be written to GENOZIP and used to reconstruct the MTF during decompression
static inline uint32_t mtf_insert_to_dict (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len, bool overlayable)
{
    buf_alloc (vb, &ctx->dict, MAX ((ctx->dict.len + snip_len + 1), INITIAL_NUM_NODES*snip_len), 2, "mtf_ctx->dict", 0);
    if (overlayable) buf_set_overlayable (&ctx->dict);

    unsigned char_index = ctx->dict.len;
    char *dict_p = &ctx->dict.data[char_index];

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = '\t'; // dictionary snips have a \t separator within dictionary string

    ctx->dict.len += snip_len + 1;

    return char_index;
}

// mtf index to node - possibly in ol_mtf, or in mtf
MtfNode *mtf_node (const MtfContext *ctx, uint32_t mtf_i, const char **snip_in_dict /* optional out */)
{
    ASSERT (mtf_i >= 0 && mtf_i < ctx->mtf.len + ctx->ol_mtf.len, "Error: mtf_i=%u out of range: mtf.len=%u ol_mtf.len=%u", mtf_i, ctx->mtf.len, ctx->ol_mtf.len);
    bool is_ol = mtf_i < ctx->ol_mtf.len; // is this entry from a previous vb (overlay buffer)

    MtfNode *node = is_ol ? &((MtfNode *)ctx->ol_mtf.data)[mtf_i] 
                          : &((MtfNode *)ctx->mtf.data)[mtf_i - ctx->ol_mtf.len];

    if (snip_in_dict) {
        const Buffer *dict = is_ol ? &ctx->ol_dict : &ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error: dict not allocated");

        *snip_in_dict = &dict->data[node->char_index];
    }
    
    return node;
}

static inline int32_t mtf_get_mtf_i (VariantBlock *vb, MtfContext *ctx, const char *snip, unsigned snip_len, 
                                     MtfNode **node,        // out - node if node is found, NULL if not
                                     HashEnt **new_hashent) // optional out - the hash entry for a new node
{
    HashEnt head, *hashent = &head;
    hashent->next = mtf_hash (ctx, snip, snip_len); // entry in hash table determined by hash function on snip
    int32_t hashent_i;

    while (hashent->next != NIL) {

        ASSERT (hashent->next >= 0 && hashent->next < ctx->hash.len, "Error: hashent->next=%d out of range, hash.len=%u", hashent->next, ctx->hash.len);
        hashent_i = hashent->next;
        hashent = &((HashEnt*)ctx->hash.data)[hashent_i];

        // case: snip is not in hash table and also no other snip occupies the slot
        if (hashent->mtf_i == NIL) { // unoccupied space in core hash table
            if (new_hashent) *new_hashent = hashent;
            return NIL;
        }

        const char *snip_in_dict;
        *node = mtf_node (ctx, hashent->mtf_i, &snip_in_dict);

        // case: snip is in the hash table 
        if (snip_len == (*node)->snip_len && !memcmp (snip, snip_in_dict, snip_len)) {
            if (new_hashent) *new_hashent = NULL; // not a new node
            return hashent->mtf_i;
        }
    }

    // case: not found in hash table, and we are required to provide a new hash entry on the linked list
    if (new_hashent) {
        if (ctx->hash.size < sizeof (HashEnt) * (1 + ctx->hash.len))
            buf_alloc (vb, &ctx->hash, sizeof (HashEnt) * (1 + ctx->hash.len), 2, "mtf_ctx->hash" , 0);

        hashent = &((HashEnt*)ctx->hash.data)[hashent_i]; // might have changed after realloc
        hashent->next = ctx->hash.len++;

        *new_hashent = &((HashEnt*)ctx->hash.data)[hashent->next];
    }

    *node = NULL;
    return NIL;
}

typedef struct {
    uint32_t char_index;
    uint32_t snip_len;
} MtfWord;

// called when pizzing a genotype section
void mtf_get_snip_by_word_index (VariantBlock *vb, MtfContext *ctx, const uint8_t *word_index_base250,
                                 char **snip, uint32_t *snip_len) // out
{
    if (word_index_base250[0] == BASE250_MISSING_SF) 
        *snip = NULL; // ignore this dict_id - don't even output a separator

    else if (word_index_base250[0] == BASE250_EMPTY_SF) {
        *snip = ""; // pointer to static empty string
        *snip_len = 0;
    }
    else {
        uint32_t word_index = base250_decode (word_index_base250);

        ASSERT (word_index < ctx->word_list.len, "Error: word_index=%u is out of bounds - %.*s dictionary has only %u entries",
                word_index, DICT_ID_LEN, ctx->dict_id.id, ctx->word_list.len);

        MtfWord *dict_word = &((MtfWord*)ctx->word_list.data)[word_index];

        *snip = &ctx->dict.data[dict_word->char_index];
        *snip_len = dict_word->snip_len;
    }
}

// process a snip as it is segregated during zip. if its the first time its seen, it is added to the dictionary
// also used for adding snips to z_file->mtf_ctx during mtf_merge_in_vb_ctx_one_dict_id()
int32_t mtf_evaluate_snip (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len, bool overlayable,
                           MtfNode **node /* out */) 
{
    if (!snip_len) 
        return (!snip || *snip != ':') ? SEG_MISSING_SF : SEG_EMPTY_SF;
    
    // allocated fixed-size hash table if not already allocated
    if (!buf_is_allocated (&ctx->hash)) {
        buf_alloc (vb, &ctx->hash, sizeof(HashEnt) * MTF_HASH_TABLE_LEN, 1, "mtf_ctx->hash", 0);

        // we set all entries to {NIL, NIL} == {-1, -1} == {0xffffffff, 0xffffffff} (note: HashEnt is packed)
        memset (ctx->hash.data, 0xff, sizeof(HashEnt) * MTF_HASH_TABLE_LEN);
        ctx->hash.len = MTF_HASH_TABLE_LEN;
    }

    // attempt to get the node from the hash table
    HashEnt *new_hashent = NULL;
    int32_t mtf_i = mtf_get_mtf_i (vb, ctx, snip, snip_len, node, &new_hashent);
    if (mtf_i != NIL) return mtf_i; // snip found - we're done

    // this snip isn't in the hash table - its a new snip
    buf_alloc (vb, &ctx->mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, 1+ctx->mtf.len), 2, "mtf_ctx->mtf" , 0);
    if (overlayable) buf_set_overlayable (&ctx->mtf);

    new_hashent->mtf_i = ctx->ol_mtf.len + ctx->mtf.len++; // new hash entry or extend linked list
    new_hashent->next  = NIL;

    *node = mtf_node (ctx, new_hashent->mtf_i, NULL);
    (*node)->snip_len   = snip_len;
    (*node)->word_index = base250_encode (new_hashent->mtf_i);    
    (*node)->char_index = mtf_insert_to_dict (vb, ctx, snip, snip_len, overlayable);
    
    // if this is the first variant block - allocate/grow sorter to contain exactly the same number of entries as mtf
    if (vb->variant_block_i == 1) {
        unsigned prev_size = ctx->sorter.size;
        buf_alloc (vb, &ctx->sorter, sizeof (SorterEnt) * (ctx->mtf.size / sizeof(MtfNode)), 1, "mtf_ctx->sorter", 0);
        if (ctx->sorter.size > prev_size) memset (&ctx->sorter.data[prev_size], 0, ctx->sorter.size - prev_size);

        ((SorterEnt *)ctx->sorter.data)[new_hashent->mtf_i].mtf_i = new_hashent->mtf_i;
    }

    return new_hashent->mtf_i;
}

// overlay and/or copy the current state of the global context to the vb, ahead of compressing this vb.
void mtf_clone_ctx (VariantBlock *vb)
{
    pthread_mutex_lock (&vb->z_file->mutex);

    START_TIMER;

    for (unsigned did_i=0; did_i < vb->z_file->num_dict_ids; did_i++) {
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];
        MtfContext *zf_ctx = &vb->z_file->mtf_ctx[did_i];

        if (buf_is_allocated (&zf_ctx->dict)) {  // something already for this dict_id
            // overlay the global dict and mtf - these will not change by this (or any other) VB
            buf_overlay (&vb_ctx->ol_dict, &zf_ctx->dict, 0,0,0,0);   
            buf_overlay (&vb_ctx->ol_mtf, &zf_ctx->mtf, 0,0,0,0);   

            // copy the hash table (core + used entries only of extension) - we're copying and not overlaying as we might be changing it
            buf_copy (vb, &vb_ctx->hash, &zf_ctx->hash, sizeof(HashEnt), 0, zf_ctx->hash.len, 0, 0);
        }

        vb_ctx->dict_id = zf_ctx->dict_id;
    }

    vb->num_dict_ids = vb->z_file->num_dict_ids;

    COPY_TIMER (vb->profile.mtf_clone_ctx)

    pthread_mutex_unlock (&vb->z_file->mutex);
}

// find the z_file context that corresponds to dict_id. It could be possibly a different did_i
// than in the vb - in case this dict_id is new to this vb, but another vb already inserted
// it to z_file
static unsigned mtf_get_z_file_did_i (VariantBlock *vb, DictIdType dict_id)
{
    for (unsigned did_i=0; did_i < vb->z_file->num_dict_ids; did_i++)
        if (dict_id.num == vb->z_file->mtf_ctx[did_i].dict_id.num)
            return did_i;

    // z_file doesn't yet have this dict_id - add it now
    ASSERT (vb->z_file->num_dict_ids+1 < MAX_DICTS, 
            "Error: z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    vb->z_file->mtf_ctx[vb->z_file->num_dict_ids].dict_id = dict_id;
    vb->z_file->num_dict_ids++;

    return vb->z_file->num_dict_ids-1;
}

void mtf_initialize_mutex (File *z_file, unsigned next_variant_i_to_merge)
{
    unsigned ret = pthread_mutex_init (&z_file->mutex, NULL);
    z_file->mutex_initialized = true;
    ASSERT0 (!ret, "pthread_mutex_init failed");

    z_file->next_variant_i_to_merge = next_variant_i_to_merge;
}

// we need to add "our" new words to the global dictionaries in the correct order of VBs
static void mtf_wait_for_my_turn(VariantBlock *vb)
{
    for (unsigned i=0; ; i++) {
        pthread_mutex_lock (&vb->z_file->mutex);

        if (vb->z_file->next_variant_i_to_merge == vb->variant_block_i) 
            return; // our turn now

        else {
            pthread_mutex_unlock (&vb->z_file->mutex);
            usleep (100000); // wait 100 millisec and try again
        }

        ASSERT0 (i < 600, "Error: timeout while waiting for mutex")
    }
}

// this is called towards the end of compressing one vb - merging its dictionaries into the z_file 
static bool mtf_merge_in_vb_ctx_one_dict_id (VariantBlock *vb, unsigned did)
{
    MtfContext *vb_ctx = &vb->mtf_ctx[did];

    unsigned z_did_i = mtf_get_z_file_did_i (vb, vb_ctx->dict_id);
    MtfContext *zf_ctx = &vb->z_file->mtf_ctx[z_did_i];

    if (!buf_is_allocated (&vb_ctx->dict)) return false; // nothing yet for this dict_id

    uint32_t start_dict_len = zf_ctx->dict.len;
    uint32_t start_mtf_len  = zf_ctx->mtf.len;

    if (!buf_is_allocated (&zf_ctx->dict)) {
        
        // first data - move to zf_ctx and leave overlay
        buf_move (vb, &zf_ctx->dict, &vb_ctx->dict);
        buf_set_overlayable (&zf_ctx->dict);
        buf_overlay (&vb_ctx->ol_dict, &zf_ctx->dict, 0,0,0,0);

        buf_move (vb, &zf_ctx->mtf,  &vb_ctx->mtf);
        buf_set_overlayable (&zf_ctx->mtf);
        buf_overlay (&vb_ctx->ol_mtf, &zf_ctx->mtf, 0,0,0,0);

        buf_move (vb, &zf_ctx->hash, &vb_ctx->hash); // vb_ctx no longer needs the hash table - zf_ctx can take it

        zf_ctx->dict_id = vb_ctx->dict_id;
    }
    else {
        // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
        for (unsigned i=0; i < vb_ctx->mtf.len; i++) {
            MtfNode *vb_node = &((MtfNode *)vb_ctx->mtf.data)[i];
            
            MtfNode *zf_node;
            uint32_t zf_node_index = mtf_evaluate_snip (vb, zf_ctx, &vb_ctx->dict.data[vb_node->char_index], vb_node->snip_len, true, &zf_node);
            ASSERT (zf_node_index < zf_ctx->mtf.len, "Error: zf_node_index=%u out of range - len=%i", zf_node_index, vb_ctx->mtf.len);

            // set word_index to be indexing the global dict - to be used by zip_generate_genotype_one_section()
            vb_node->word_index = zf_node->word_index;
        }        
    }

    // we now compress the dictionaries directly from z_file. note: we must continue to hold
    // the mutex during compression, lest another thread re-alloc the dictionary.

    const char *start_dict = &zf_ctx->dict.data[start_dict_len]; // we take the pointer AFTER the evaluate, since dict can be reallocted
    unsigned added_chars = zf_ctx->dict.len - start_dict_len;
    unsigned added_words = zf_ctx->mtf.len  - start_mtf_len;

    // compress incremental part of dictionary added by this vb. note: dispatcher calls this function in the correct order of VBs.
    if (added_chars) {

        // special optimization for the GL dictionary
        if (!memcmp (&zf_ctx->dict_id.id, "GL\0\0\0\0\0\0", DICT_ID_LEN)) 
            start_dict = gl_optimize_dictionary (vb, &zf_ctx->dict, &((MtfNode *)zf_ctx->mtf.data)[start_mtf_len], start_dict_len, added_words);
     
        zfile_compress_dictionary_data (vb, zf_ctx->dict_id, added_words, start_dict, added_chars);

        vb->add_bytes[SEC_DICTIONARY] += added_chars;
    }

    return added_chars > 0;
}

// merge new words added in this vb into the z_file.mtf_ctx, and compresses dictionaries
// while holding exclusive access to the z_file dictionaries. returns num_dictionary_sections
unsigned mtf_merge_in_vb_ctx (VariantBlock *vb)
{
    mtf_wait_for_my_turn(vb); // we grab the mutex in the sequencial order of VBs

    START_TIMER; // note: careful not to count time spent waiting for the mutex

    unsigned num_dictionary_sections = 0;
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) 
        num_dictionary_sections += mtf_merge_in_vb_ctx_one_dict_id (vb, did_i);

    vb->z_file->next_variant_i_to_merge++;

    COPY_TIMER (vb->profile.mtf_merge_in_vb_ctx)

    pthread_mutex_unlock (&vb->z_file->mutex);

    return num_dictionary_sections;
}

unsigned mtf_get_did_i_by_dict_id (MtfContext *mtf_ctx, unsigned *num_dict_ids, 
                                   unsigned *num_subfields, // used only if dict_id is a subfield
                                   DictIdType dict_id)
{
    // check if we have this dict_id already
    unsigned did_i=0 ; for (; did_i < *num_dict_ids; did_i++) 
        if (dict_id.num == mtf_ctx[did_i].dict_id.num) break;

    // case: dict_id encountered for this first time - initialize a mtf_ctx
    if (did_i == *num_dict_ids) {

        ASSERT (*num_dict_ids+1 < MAX_DICTS, 
                "Error: number of dictionary types is greater than MAX_DICTS=%u", MAX_DICTS);


        mtf_ctx[*num_dict_ids].dict_id = dict_id;

        // thread safety: the increment below MUST be AFTER memcpy, bc piz_get_line_subfields
        // might be reading this data at the same time as the piz dispatcher thread adding more
        // dictionaries
        (*num_dict_ids)++; 

        if (num_subfields) (*num_subfields)++;
    }

    return did_i;
}

// this is called by the piz dispatcher thread after reading a dictionary section 
void mtf_integrate_dictionary_fragment (VariantBlock *vb, char *section_data)
{    
    START_TIMER;

    // thread safety note: this function is called only from the piz dispatcher thread,
    // so no thread safety issues with this static buffer.
    static Buffer fragment = EMPTY_BUFFER;

    // thread-safety note: while the dispatcher thread is integrating new dictionary fragments,
    // compute threads might be using these dictionaries. This is ok, bc the dispatcher thread makes
    // sure we integrate dictionaries from vbs by order - so that running compute threads never
    // need to access the new parts of dictionaries. We also pre-allocate the dictionaries in
    // vcf_header_genozip_to_vcf() so that they don't need to be realloced. dict.len may be accessed
    // by compute threads, but its change is assumed to be atomic, so that no weird things will happen
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)section_data;
    uint32_t num_snips = BGEN32 (header->num_snips);

    zfile_uncompress_section (vb, section_data, &fragment, SEC_DICTIONARY);

    // special treatment if this is GL - de-optimize
    if (!memcmp (header->dict_id.id, "GL\0\0\0\0\0\0", DICT_ID_LEN))
        gl_deoptimize_dictionary (fragment.data, fragment.len);

    if (flag_show_dict) 
        fprintf (stderr, "%*.*s (vb_i=%u):\t%.*s\n", -DICT_ID_LEN, DICT_ID_LEN, header->dict_id.id, vb->variant_block_i, fragment.len, fragment.data);

    // in piz, the same did_i is used for z_file and vb contexts, meaning that in vbs there could be
    // a non-contiguous array of contexts (some are missing if not used by this vb)
    unsigned did_i = mtf_get_did_i_by_dict_id (vb->z_file->mtf_ctx, &vb->z_file->num_dict_ids, NULL, header->dict_id);
    
    Buffer *dict = &vb->z_file->mtf_ctx[did_i].dict;

    // append fragment to dict. If there is no room - old memory is abandoned (so that VBs that are overlaying
    // it continue to work uninterrupted) and a new memory is allocated, where the old dict is joined by the new fragment
    unsigned dict_old_len = dict->len;
    buf_alloc (vb, dict, dict->len + fragment.len, 2, "mtf_ctx->dict", did_i);
    buf_set_overlayable (dict);

    memcpy (&dict->data[dict->len], fragment.data, fragment.len);
    dict->len += fragment.len;

    Buffer *word_list_buf = &vb->z_file->mtf_ctx[did_i].word_list;

    // extend word list memory - and calculate the new words. If there is no room - old memory is abandoned 
    // (so that VBs that are overlaying it continue to work uninterrupted) and a new memory is allocated
    buf_alloc (vb, word_list_buf, (word_list_buf->len + num_snips) * sizeof (MtfWord), 2, "mtf_ctx->word_list", did_i);
    buf_set_overlayable (word_list_buf);

    char *start = fragment.data;
    for (unsigned snip_i=0; snip_i < num_snips; snip_i++) {
        MtfWord *word = &((MtfWord *)word_list_buf->data)[word_list_buf->len++];

        char *c=start; while (*c != '\t') c++;
        
        word->snip_len   = c - start;
        word->char_index = dict_old_len + (start - fragment.data);

        start = c+1; // skip over the \t
    }

    buf_free (&fragment);

    COPY_TIMER(vb->profile.mtf_integrate_dictionary_fragment);
}

// here we overlay the z_file's dictionaries and word lists to the vb. these data remain unchanged - neither
// the vb nor the dispatcher thread will ever change snips placed in these. the dispatcher thread may append
// the dictionary and word list as new fragments become available from subsequent VBs. If the memory is not 
// sufficient, the dispatcher thread will "abandon" this memory, leaving it to the VB to continue to use it
// while starting a larger dict/word_list on a fresh memory allocation.
void mtf_overlay_dictionaries_to_vb (VariantBlock *vb)
{
    for (unsigned did_i=0; did_i < MAX_DICTS; did_i++) {
        MtfContext *zf_ctx = &vb->z_file->mtf_ctx[did_i];
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];

        if (!zf_ctx->dict_id.id[0]) break;

        if (buf_is_allocated (&zf_ctx->dict) && buf_is_allocated (&zf_ctx->word_list)) { 
            vb_ctx->dict_id = zf_ctx->dict_id;
            buf_overlay (&vb_ctx->dict, &zf_ctx->dict, 0,0,0,0);    
            buf_overlay (&vb_ctx->word_list, &zf_ctx->word_list, 0,0,0,0);

            // count dictionaries of genotype data subfields
            if (dict_id_is_gt_subfield (vb_ctx->dict_id)) {
                vb->num_subfields++;
                ASSERT (vb->num_subfields <= MAX_SUBFIELDS, 
                        "Error: number of subfields in %s exceeds MAX_SUBFIELDS=%u, while reading vb_i=%u", 
                        file_printname (vb->z_file), MAX_SUBFIELDS, vb->variant_block_i);
            }
        }
    }
}

static int sorter_cmp(const void *a_, const void *b_)  
{ 
    SorterEnt *a = (SorterEnt *)a_;
    SorterEnt *b = (SorterEnt *)b_;
    
    return (int)b->count - (int)a->count;
}

void mtf_sort_dictionaries_vb_1(VariantBlock *vb)
{
    // thread safety note: no issues here, as this is run only by the compute thread of variant_block_i=1
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {

        MtfContext *ctx = &vb->mtf_ctx[did_i];

        // sort in ascending order of mtf->count
        qsort (ctx->sorter.data, ctx->mtf.len, sizeof(SorterEnt), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in mtf
        static Buffer old_dict = EMPTY_BUFFER;
        buf_move (vb, &old_dict, &ctx->dict);

        buf_alloc (vb, &ctx->dict, old_dict.size, 1, "mtf_ctx->dict", did_i);
        ctx->dict.len = old_dict.len;

        char *next = ctx->dict.data;
        for (unsigned i=0; i < ctx->mtf.len; i++) {
            int32_t mtf_i = ((SorterEnt *)ctx->sorter.data)[i].mtf_i;
            MtfNode *node = &((MtfNode *)ctx->mtf.data)[mtf_i];
            memcpy (next, &old_dict.data[node->char_index], node->snip_len + 1 /* +1 for \t */);
            
            node->word_index = base250_encode (i);
            node->char_index = next - ctx->dict.data;
            
            next += node->snip_len + 1;
        }

        buf_free (&old_dict);
    }
}

void mtf_free_context (MtfContext *ctx)
{
    buf_free (&ctx->ol_dict);
    buf_free (&ctx->ol_mtf);
    buf_free (&ctx->dict);
    buf_free (&ctx->mtf);
    buf_free (&ctx->word_list);
    buf_free (&ctx->hash);
    buf_free (&ctx->sorter);
    ctx->dict_id.num = 0;
}

