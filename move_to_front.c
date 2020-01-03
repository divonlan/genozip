// ------------------------------------------------------------------
//   move-to-front.c
//   Copyright (C) 2019 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

/*
zip:
    1) during segregate - build mtf_context + dictionary for each subfield

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

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "genozip.h"

#define INITIAL_NUM_NODES 10000

// tested hash table sizes up to 5M. turns out smaller tables (up to a point) are faster, despite having longer
// average linked lists. probably bc the CPU can store the entire hash and mtf arrays in L1 or L2
// memory cache during segmentation
#define MTF_HASH_TABLE_LEN 19997         // a prime number, see: https://primes.utm.edu/lists/small/millions/
static inline int *mtf_hash (const MtfContext *ctx, const char *snip, unsigned snip_len)
{
    uint64_t result=0;
    for (unsigned i=0; i < snip_len; i++) 
        result = ((result << 23) | (result >> 41)) ^ (uint64_t)((uint8_t)snip[i]);

    uint32_t hash = (uint32_t)(result % MTF_HASH_TABLE_LEN);
    
    return &((int *)ctx->hash.data)[hash]; // pointer to the hash table entry containing mtf_i or NIL
}

static inline unsigned mtf_snip_len_t (const char *snip)
{
    const char *s; for (s=snip ; *s != '\t' ; s++);
    return s - snip;
}

// add a snip to the dictionary the first time it is encountered in the VCF file.
// the dictionary will be written to GENOZIP and used to reconstruct the MTF during decompression
static inline uint32_t mtf_insert_to_dict (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len)
{
    if (!buf_is_allocated (&ctx->dict) || ctx->dict.len + snip_len + 1 > ctx->dict.size)  // full or not allocated 
        buf_alloc (vb, &ctx->dict, MAX ((ctx->dict.len + snip_len + 1), INITIAL_NUM_NODES*snip_len), 2, "mtf_ctx->dict", 0);
        
    unsigned char_index = ctx->dict.len;
    char *dict_p = &ctx->dict.data[char_index];

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = '\t'; // dictionary snips have a \t separator within dictionary string

    ctx->dict.len += snip_len + 1;

    return char_index;
}

static inline int mtf_get_mtf_i (MtfContext *ctx, const char *snip, unsigned snip_len, 
                                 int **next) // optional out
{
    // check if its in the hash table
    int *hash_ent = mtf_hash (ctx, snip, snip_len);
    int mtf_i = *hash_ent;
    MtfNode *node = NULL;

    while (mtf_i != NIL) {

        node = &((MtfNode *)ctx->mtf.data)[mtf_i];

        // case: snip is in the hash table 
        if (snip_len == node->snip_len && !memcmp (snip, &ctx->dict.data[node->char_index], snip_len)) 
            return mtf_i;

        mtf_i = node->next; // move through the linked list for this hash value
    }

    // not found in hash table
    if (next) *next = (node ? &node->next : hash_ent); // the place to add the index of a new node
    return NIL;
}

// gets the index of a snip in the dictionary 
// returning the index in base-250. called when generating genotype sections
Base250 mtf_get_index_by_snip (VariantBlock *vb, MtfContext *ctx, char **src, unsigned snip_len)
{
    // get the node from the hash table
    int mtf_i = mtf_get_mtf_i (ctx, *src, snip_len, NULL);
    ASSERT (mtf_i != NIL, "Error: snip %.*s not found", snip_len, *src);

    *src += snip_len;

    return (&((MtfNode *)ctx->mtf.data)[mtf_i])->word_index;
}

// called when pizzing a genotype section
void mtf_get_snip_by_word_index (VariantBlock *vb, MtfContext *ctx, const uint8_t *word_index_base250,
                                 char **snip, uint32_t *snip_len) // out
{
    if (word_index_base250[0] == BASE250_MISSING_SF) 
        *snip = NULL; // ignore this subfield - don't even output a separator

    else if (word_index_base250[0] == BASE250_EMPTY_SF) {
        *snip = ""; // pointer to static empty string
        *snip_len = 0;
    }
    else {
        uint32_t word_index = base250_decode (word_index_base250);

        ASSERT (word_index < ctx->word_list.len, "Error: word_index=%u is out of bounds - %.*s dictionary has only %u entries",
                word_index, SUBFIELD_ID_LEN, ctx->subfield.id, ctx->word_list.len);

        MtfWord *dict_word = &((MtfWord*)ctx->word_list.data)[word_index];

        *snip = &ctx->dict.data[dict_word->char_index];
        *snip_len = dict_word->snip_len;
    }
}

// process a snip as it is segregated during zip. if its the first time its seen,
// it is added to the dictionary, tree and hash table.
uint32_t mtf_evaluate_snip (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len) 
{
    if (!snip_len) 
        return (!snip || *snip != ':') ? SEG_MISSING_SF : SEG_EMPTY_SF;
    
    // allocate memory - start with 1% of entries, we will increase if its not enough
    if (!buf_is_allocated (&ctx->mtf) || ctx->mtf.size < (ctx->mtf.len+1) * sizeof (MtfNode)) // full or not allocated 
        // start with INITIAL_NUM_NODES nodes, and double if we run out
        buf_alloc (vb, &ctx->mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, 1+ctx->mtf.len), 
                   2, "mtf_ctx->mtf" , 0);

    if (!buf_is_allocated (&ctx->hash)) {
        buf_alloc (vb, &ctx->hash, sizeof(int) * MTF_HASH_TABLE_LEN, 1, "mtf_ctx->hash", 0);
        int *hash = (int *)ctx->hash.data;
        for (unsigned i=0; i < MTF_HASH_TABLE_LEN; i++) hash[i] = NIL;
        ctx->hash.len = MTF_HASH_TABLE_LEN;
    }

    // attempt to get the node from the hash table
    int *next;
    int mtf_i = mtf_get_mtf_i (ctx, snip, snip_len, &next);
    if (mtf_i != NIL) return mtf_i;

    // this snip isn't in the hash table - its a new snip
    *next = ctx->mtf.len++; // new hash entry or extend linked list

    MtfNode *node = N(*next); // new snip
    node->next       = NIL;
    node->snip_len   = snip_len;
    node->word_index = base250_encode (*next);
    node->char_index = mtf_insert_to_dict (vb, ctx, snip, snip_len);
    node->count      = 0;

    return *next;
}

static void mtf_copy_ctx (VariantBlock *vb, MtfContext *dst, MtfContext *src)
{
    if (buf_is_allocated (&src->dict)) {  // something already for this subfield
        buf_copy (vb, &dst->dict, &src->dict, 1, 0, 0, "mtf_ctx->dict", vb->variant_block_i);
        buf_copy (vb, &dst->hash, &src->hash, sizeof(int), 0, 0, "mtf_ctx->hash", vb->variant_block_i);
        buf_copy (vb, &dst->mtf,  &src->mtf,  sizeof(MtfNode), 0, 0, "mtf_ctx->mtf", vb->variant_block_i);
    }

    dst->subfield = src->subfield;
}

// copy the current state of the global context to the vb, ahead of compressing this vb.
void mtf_clone_ctx (VariantBlock *vb)
{
    pthread_mutex_lock (&vb->z_file->mutex);

    START_TIMER;

    for (unsigned sf=0; sf < vb->z_file->num_subfields; sf++) {
        MtfContext *dst = &vb->mtf_ctx[sf];
        MtfContext *src = &vb->z_file->mtf_ctx[sf];

        mtf_copy_ctx (vb, dst, src);
        dst->cloned_mtf_len = dst->mtf.len;
    }

    vb->num_subfields = vb->z_file->num_subfields;

    COPY_TIMER (vb->profile.mtf_clone_ctx)

    pthread_mutex_unlock (&vb->z_file->mutex);
}

// find the z_file context that corresponds to subfield. It could be possibly a different sf_i
// than in the vb - in case this subfield is new to this vb, but another vb already inserted
// it to z_file
static unsigned mtf_get_z_file_sf_i (VariantBlock *vb, SubfieldIdType subfield)
{
    for (unsigned sf_i=0; sf_i < vb->z_file->num_subfields; sf_i++)
        if (!memcmp (subfield.id, &vb->z_file->mtf_ctx[sf_i].subfield.id, SUBFIELD_ID_LEN))
            return sf_i;

    // z_file doesn't yet have this subfield - add it now
    ASSERT (vb->z_file->num_subfields+1 < MAX_SUBFIELDS, 
            "Error: z_file has more subfield types than MAX_SUBFIELDS=%u", MAX_SUBFIELDS);

    vb->z_file->mtf_ctx[vb->z_file->num_subfields].subfield = subfield;
    vb->z_file->num_subfields++;

    return vb->z_file->num_subfields-1;
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
bool mtf_merge_in_vb_ctx_one_subfield(VariantBlock *vb, unsigned sf)
{
    MtfContext *vb_ctx = &vb->mtf_ctx[sf];

    unsigned z_sf_i = mtf_get_z_file_sf_i (vb, vb_ctx->subfield);
    MtfContext *zf_ctx = &vb->z_file->mtf_ctx[z_sf_i];
    
    if (!buf_is_allocated (&vb_ctx->dict)) return false; // nothing yet for this subfield

    uint32_t start_dict_len = zf_ctx->dict.len;
    uint32_t start_mtf_len  = zf_ctx->mtf.len;

    if (!buf_is_allocated (&zf_ctx->dict)) {
        // first data - just copy
        mtf_copy_ctx (vb, zf_ctx, vb_ctx);
    }
    else {
        // merge in words that are potentially new (but may have been already added
        // by other VBs since we cloned for this VB)
        for (unsigned i=vb_ctx->cloned_mtf_len; i < vb_ctx->mtf.len; i++) {
            MtfNode *node = &((MtfNode *)vb_ctx->mtf.data)[i];
            
            uint32_t zf_node_index = mtf_evaluate_snip (vb, zf_ctx, &vb_ctx->dict.data[node->char_index], node->snip_len);
            ASSERT (zf_node_index < zf_ctx->mtf.len, 
                    "Error: zf_node_index=%u out of range - len=%i", zf_node_index, vb_ctx->mtf.len);

            MtfNode *zf_node = &((MtfNode *)zf_ctx->mtf.data)[zf_node_index];

            // update indexes to be indexing the global dict
            node->char_index = zf_node->char_index;
            node->word_index = zf_node->word_index;
        }        
    }

    // we now compress the dictionaries directly from z_file. note: we must continue to hold
    // the mutex during compression, lest another thread re-alloc the dictionary.

    char *start_dict = &zf_ctx->dict.data[start_dict_len]; // we take the pointer AFTER the evaluate, since dict can be reallocted
    unsigned added_chars = zf_ctx->dict.len - start_dict_len;
    unsigned added_words = zf_ctx->mtf.len  - start_mtf_len;

    // compress incremental part of dictionary added by this vb. note: dispatcher
    // calls this function in the correct order of VBs.
    if (added_chars) {

        // special optimization for the GL dictionary
        if (!memcmp (&zf_ctx->subfield.id, "GL\0\0\0\0\0\0", SUBFIELD_ID_LEN)) {
            gl_optimize_dictionary (vb, &zf_ctx->dict, &((MtfNode *)zf_ctx->mtf.data)[start_mtf_len],
                                    start_dict_len, added_words);

            zfile_compress_dictionary_data (vb, zf_ctx->subfield, added_words, vb->optimized_gl_dict.data, added_chars);
        }
        else
            zfile_compress_dictionary_data (vb, zf_ctx->subfield, added_words, start_dict, added_chars);

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
    for (unsigned sf=0; sf < vb->num_subfields; sf++) 
        num_dictionary_sections += mtf_merge_in_vb_ctx_one_subfield (vb, sf);

    vb->z_file->next_variant_i_to_merge++;

    COPY_TIMER (vb->profile.mtf_merge_in_vb_ctx)

    pthread_mutex_unlock (&vb->z_file->mutex);

    return num_dictionary_sections;
}

unsigned mtf_get_sf_i_by_subfield (MtfContext *mtf_ctx, unsigned *num_subfields, SubfieldIdType subfield)
{
    // check if we have this subfield already
    unsigned sf_i=0 ; for (; sf_i < *num_subfields; sf_i++) {
        if (!memcmp (subfield.id, mtf_ctx[sf_i].subfield.id, SUBFIELD_ID_LEN)) break;
    }

    // case: subfield encountered for this first time - initialize a mtf_ctx
    if (sf_i == *num_subfields) {

        ASSERT (*num_subfields+1 < MAX_SUBFIELDS, 
                "Error: number of subfield types than MAX_SUBFIELDS=%u", MAX_SUBFIELDS);

        memcpy (mtf_ctx[*num_subfields].subfield.id, subfield.id, SUBFIELD_ID_LEN); 

        // thread safety: the increment below MUST be AFTER memcpy, bc piz_get_line_subfields
        // might be reading this data at the same time as the piz dispatcher thread adding more
        // dictionaries
        (*num_subfields)++; 
    }

    return sf_i;
}

// this is called by the piz dispatcher thread after reading a dictionary section 
void mtf_integrate_dictionary_fragment (VariantBlock *vb, const char *section_data)
{
    START_TIMER;

    // thread safety note: this function is called only from the piz dispatcher thread,
    // so no thread safety issues with this static buffer.
    static Buffer fragment = EMPTY_BUFFER;

    // thread-safety note: while the dispatcher thread is integrating new dictionary fragments,
    // compute threads might be using these dictionaries. This is ok, bc the dispatcher thread makes
    // sure we integrate dictionaries from vbs by order - so that running compute threads never
    // need to access the new parts of dictionaries. We also pre-allocate the dictionaries in
    // vcf_header_vcz_to_vcf() so that they don't need to be realloced. dict.len may be accessed
    // by compute threads, but its change is assumed to be atomic, so that no weird things will happen
    SectionHeaderDictionary *header = (SectionHeaderDictionary *)section_data;
    uint32_t num_snips = ENDN32 (header->num_snips);
    SubfieldIdType subfield; memcpy (subfield.id, header->subfield_id, SUBFIELD_ID_LEN);

    zfile_uncompress_section (vb, section_data, &fragment, SEC_DICTIONARY);

    // special treatment if this is GL - de-optimize
    if (!memcmp (header->subfield_id, "GL\0\0\0\0\0\0", SUBFIELD_ID_LEN))
        gl_deoptimize_dictionary (fragment.data, fragment.len);

    // in piz, the same sf_i is used for z_file and vb contexts, meaning that in vbs there could be
    // a non-contiguous array of contexts (some are missing if not used by this vb)
    unsigned sf_i = mtf_get_sf_i_by_subfield (vb->z_file->mtf_ctx, &vb->z_file->num_subfields, subfield);
    
    Buffer *dict = &vb->z_file->mtf_ctx[sf_i].dict;

    // append fragment to dict. If there is no room - old memory is abandoned (so that VBs that are overlaying
    // it continue to work uninterrupted) and a new memory is allocated, where the old dict is joined by the new fragment
    unsigned dict_old_len = dict->len;
    buf_append (vb, dict, &fragment, 1, "ctx->dict", sf_i);

    Buffer *word_list_buf = &vb->z_file->mtf_ctx[sf_i].word_list;

    // extend word list memory - and calculate the new words. If there is no room - old memory is abandoned 
    // (so that VBs that are overlaying it continue to work uninterrupted) and a new memory is allocated
    buf_extend (vb, word_list_buf, num_snips, sizeof (MtfWord), "word_list_buf", sf_i);
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
    for (unsigned sf=0; sf < MAX_SUBFIELDS; sf++) {
        MtfContext *zf_ctx = &vb->z_file->mtf_ctx[sf];
        MtfContext *vb_ctx = &vb->mtf_ctx[sf];

        if (!zf_ctx->subfield.id[0]) break;

        if (buf_is_allocated (&zf_ctx->dict) && buf_is_allocated (&zf_ctx->word_list)) { 
            vb_ctx->subfield = zf_ctx->subfield;
            buf_overlay (&vb_ctx->dict, &zf_ctx->dict, NULL, NULL, "mtf_ctx->dict", sf);    
            buf_overlay (&vb_ctx->word_list, &zf_ctx->word_list, NULL, NULL, "mtf_ctx->word_list", sf);
        }
    }
}

static int sorter_cmp(const void *a_, const void *b_)  
{ 
    MtfNode *a = *(MtfNode**)a_;
    MtfNode *b = *(MtfNode**)b_;
    
    return (int)b->count - (int)a->count;
}

void mtf_sort_dictionaries_vb_1(VariantBlock *vb)
{
    // thread safety note: no issues here, as this is run only by the compute thread of variant_block_i=1
    for (unsigned sf=0; sf < vb->num_subfields; sf++) {

        MtfContext *ctx = &vb->mtf_ctx[sf];

        static Buffer sorter_buf = EMPTY_BUFFER; // indices of ctx->mtf
        buf_alloc (vb, &sorter_buf, ctx->mtf.len * sizeof(MtfNode *), 2, "sorter_buf", 0);
        MtfNode **sorter = (MtfNode **)sorter_buf.data;

        // initialize - each entry contains and index into ctx->mtf
        for (unsigned i=0; i < ctx->mtf.len; i++) 
            sorter[i] = &((MtfNode *)ctx->mtf.data)[i];

        // sort in ascending order of mtf->count
        qsort (sorter, ctx->mtf.len, sizeof(MtfNode *), sorter_cmp);

        // rebuild dictionary is the sorted order, and update char and word indices in mtf
        static Buffer new_dict;
        buf_alloc (vb, &new_dict, ctx->dict.size, 1, "new_dict", sf);

        char *next = new_dict.data;
        for (unsigned i=0; i < ctx->mtf.len; i++) {
            memcpy (next, &ctx->dict.data[sorter[i]->char_index], sorter[i]->snip_len + 1 /* \t */);

            sorter[i]->word_index = base250_encode (i);
            sorter[i]->char_index = next - new_dict.data;

            next += sorter[i]->snip_len + 1;
        }

        // copy update dictionary to its place
        memcpy (ctx->dict.data, new_dict.data, ctx->dict.size);

        buf_free (&sorter_buf);
        buf_free (&new_dict);
    }
}

void mtf_free_context (MtfContext *ctx)
{
    buf_free (&ctx->dict);
    buf_free (&ctx->mtf);
    buf_free (&ctx->word_list);
    buf_free (&ctx->hash);
    memset (&ctx->subfield.id, 0, SUBFIELD_ID_LEN);
    ctx->cloned_mtf_len = 0;
}

