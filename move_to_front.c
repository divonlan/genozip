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
#include "random_access.h"

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
static inline uint32_t mtf_insert_to_dict (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len)
{
    buf_alloc (vb, &ctx->dict, MAX ((ctx->dict.len + snip_len + 1), INITIAL_NUM_NODES * MIN (10, snip_len)), 2, "mtf_ctx->dict", ctx->did_i);
    if (ctx->encoding != B250_ENC_NONE) buf_set_overlayable (&ctx->dict); // during merge

    unsigned char_index = ctx->dict.len;
    char *dict_p = &ctx->dict.data[char_index];

    memcpy (dict_p, snip, snip_len);
    dict_p[snip_len] = '\t'; // dictionary snips have a \t separator within dictionary string

    ctx->dict.len += snip_len + 1;

    return char_index;
}

// ZIP only (PIZ doesn't have mtf) mtf index to node - possibly in ol_mtf, or in mtf
MtfNode *mtf_node (const MtfContext *ctx, uint32_t mtf_i, 
                   const char **snip_in_dict, uint32_t *snip_len) // optional outs
{
    ASSERT0 (ctx->dict_id.num, "Error: this ctx is not initialized");
    ASSERT (mtf_i < ctx->mtf.len + ctx->ol_mtf.len, "Error: mtf_i=%u out of range: mtf.len=%u ol_mtf.len=%u", mtf_i, ctx->mtf.len, ctx->ol_mtf.len);

    bool is_ol = mtf_i < ctx->ol_mtf.len; // is this entry from a previous vb (overlay buffer)

    MtfNode *node = is_ol ? &((MtfNode *)ctx->ol_mtf.data)[mtf_i] 
                          : &((MtfNode *)ctx->mtf.data)[mtf_i - ctx->ol_mtf.len];

    if (snip_in_dict) {
        const Buffer *dict = is_ol ? &ctx->ol_dict : &ctx->dict;
        ASSERT0 (buf_is_allocated (dict), "Error: dict not allocated");

        *snip_in_dict = &dict->data[node->char_index];
    }

    if (snip_len) *snip_len = node->snip_len;
    
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
        uint32_t snip_len_in_dict;
        *node = mtf_node (ctx, hashent->mtf_i, &snip_in_dict, &snip_len_in_dict);

        // case: snip is in the hash table 
        if (snip_len == snip_len_in_dict && !memcmp (snip, snip_in_dict, snip_len)) {
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

// PIZ only (uses word_list): called when pizzing a genotype section - returns snip and index, and advances the index
uint32_t mtf_get_next_snip (VariantBlock *vb, MtfContext *ctx, 
                            SnipIterator *override_iterator,   // if NULL, taken from ctx
                            const char **snip, uint32_t *snip_len, // optional out
                            uint32_t vcf_line) 
{
    SnipIterator *iterator = override_iterator ? override_iterator : &ctx->iterator;
    
    if (!override_iterator && !iterator->next_b250) // INFO and Field1-9 data (GT data uses override_next_b250)
        iterator->next_b250 = (uint8_t *)ctx->b250.data; // initialize (GT data initializes to the beginning of each sample rather than the beginning of the data)

    uint32_t word_index = base250_decode (&iterator->next_b250, ctx ? ctx->encoding : B250_ENC_NONE); // if this line has no non-GT subfields, it will not have a ctx 

    // case: a subfield snip is missing - either the genotype data has less subfields than declared in FORMAT, or not provided at all for some (or all) samples.
    if (word_index == WORD_INDEX_MISSING_SF) {
        ASSERT (!ctx || ctx->b250_section_type == SEC_GENOTYPE_DATA, "Error while reconstrucing line %u vb_i=%u: BASE250_MISSING_SF unexpectedly found in b250 data of %.*s (%s)",
                vcf_line, vb->variant_block_i, DICT_ID_LEN, dict_id_printable(ctx->dict_id).id, st_name(ctx->b250_section_type)); // there will be no context if this GT subfield was always missing - never appeared on any sample

        if (snip) {
            *snip = NULL; // ignore this dict_id - don't even output a separator
            *snip_len = 0;
        }
    }

    // case: a subfield snip is empty, eg AB::CD
    else if (word_index == WORD_INDEX_EMPTY_SF) { 
        ASSERT (ctx->b250_section_type == SEC_GENOTYPE_DATA, "Error while reconstrucing line %u: BASE250_EMPTY_SF unexpectedly found in b250 data of %.*s",
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
                ctx->dict_section_type == SEC_INFO_SUBFIELD_DICT ? "INFO" : "",
                ctx->dict_section_type == SEC_FRMT_SUBFIELD_DICT ? "FORMAT" : "",
                DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, ctx->word_list.len);

        MtfWord *dict_word = &((MtfWord*)ctx->word_list.data)[word_index];

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
int32_t mtf_evaluate_snip (VariantBlock *vb, MtfContext *ctx, const char *snip, uint32_t snip_len,
                           MtfNode **node /* out */, bool *is_new /* optional out */) 
{
    vb->z_section_entries[ctx->b250_section_type]++; 

    if (!snip_len) 
        return (!snip || *snip != ':') ? WORD_INDEX_MISSING_SF : WORD_INDEX_EMPTY_SF;

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
    if (mtf_i != NIL) {
        if (is_new) *is_new = false;

        if (vb->variant_block_i == 1) {
            SorterEnt *sorter_ent = &((SorterEnt *)ctx->sorter.data)[mtf_i];
            sorter_ent->count++;
        }

        return mtf_i; // snip found - we're done
    }
    
    // this snip isn't in the hash table - its a new snip

    ASSERT (ctx->mtf.len < 0x7fffffff, "Error: too many words in directory %.*s", DICT_ID_LEN, dict_id_printable (ctx->dict_id).id);

    vb->z_section_entries[ctx->dict_section_type]++; 

    buf_alloc (vb, &ctx->mtf, sizeof (MtfNode) * MAX(INITIAL_NUM_NODES, 1+ctx->mtf.len), 2, "mtf_ctx->mtf", ctx->did_i);
    if (ctx->encoding != B250_ENC_NONE) buf_set_overlayable (&ctx->mtf); // when called from merge

    new_hashent->mtf_i = ctx->ol_mtf.len + ctx->mtf.len++; // new hash entry or extend linked list
    new_hashent->next  = NIL;

    *node = mtf_node (ctx, new_hashent->mtf_i, NULL, NULL);
    memset (*node, 0, sizeof(MtfNode)); // safety
    (*node)->snip_len   = snip_len;
    (*node)->char_index = mtf_insert_to_dict (vb, ctx, snip, snip_len);
    (*node)->word_index.n = new_hashent->mtf_i;

    // if this is the first variant block - allocate/grow sorter to contain exactly the same number of entries as mtf
    if (vb->variant_block_i == 1) {
        unsigned prev_size = ctx->sorter.size;
        buf_alloc (vb, &ctx->sorter, sizeof (SorterEnt) * (ctx->mtf.size / sizeof(MtfNode)), 1, "mtf_ctx->sorter", 0);
        if (ctx->sorter.size > prev_size) memset (&ctx->sorter.data[prev_size], 0, ctx->sorter.size - prev_size);

        SorterEnt *sorter_ent = &((SorterEnt *)ctx->sorter.data)[new_hashent->mtf_i];
        sorter_ent->mtf_i = new_hashent->mtf_i;
        sorter_ent->count = 1;
    }

    if (is_new) *is_new = true;

    return new_hashent->mtf_i;
}

static void mtf_init_mapper (VariantBlock *vb, VcfFields field_i, Buffer *mapper_buf, const char *name)
{
    if (!buf_is_allocated (&vb->mtf_ctx[field_i].ol_mtf)) return;
        
    mapper_buf->len = vb->mtf_ctx[field_i].ol_mtf.len;
    
    buf_alloc (vb, mapper_buf,  mapper_buf->len * sizeof (SubfieldMapperZip), 2, name, 0);
    
    for (unsigned i=0; i < mapper_buf->len; i++) 
        ((SubfieldMapperZip *)mapper_buf->data)[i].num_subfields = NIL;
}

// ZIP only: overlay and/or copy the current state of the global context to the vb, ahead of compressing this vb.
void mtf_clone_ctx (VariantBlock *vb)
{
    pthread_mutex_lock (&vb->z_file->mutex);

    START_TIMER; // careful not to include the mutex waiting in the time measured

    for (unsigned did_i=0; did_i < vb->z_file->num_dict_ids; did_i++) {
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];
        MtfContext *zf_ctx = &vb->z_file->mtf_ctx[did_i];

        if (buf_is_allocated (&zf_ctx->dict)) {  // something already for this dict_id

            // overlay the global dict and mtf - these will not change by this (or any other) VB
            //printf ("mtf_clone_ctx: overlaying old dict %.8s, to vb_i=%u vb_did_i=z_did_i=%u\n", dict_id_printable (zf_ctx->dict_id).id, vb->variant_block_i, did_i);
            buf_overlay (&vb_ctx->ol_dict, &zf_ctx->dict, 0,0,0,0);   
            buf_overlay (&vb_ctx->ol_mtf, &zf_ctx->mtf, 0,0,0,0);   

            // copy the hash table (core + used entries only of extension) - we're copying and not overlaying as we might be changing it
            buf_copy (vb, &vb_ctx->hash, &zf_ctx->hash, sizeof(HashEnt), 0, zf_ctx->hash.len, 0, 0);
        }

        vb_ctx->did_i             = did_i;
        vb_ctx->dict_id           = zf_ctx->dict_id;
        vb_ctx->dict_section_type = zf_ctx->dict_section_type;
        vb_ctx->b250_section_type = zf_ctx->b250_section_type;
        vb_ctx->encoding          = zf_ctx->encoding; // minimum encoding for this dict, VB might increase it further in mtf_decide_encoding
        mtf_init_iterator (vb_ctx);
    }

    vb->num_dict_ids = vb->z_file->num_dict_ids;

    COPY_TIMER (vb->profile.mtf_clone_ctx)

    pthread_mutex_unlock (&vb->z_file->mutex);

    // initialize mappers for FORMAT and INFO
    mtf_init_mapper (vb, FORMAT, &vb->format_mapper_buf, "format_mapper_buf");    
    mtf_init_mapper (vb, INFO, &vb->iname_mapper_buf, "iname_mapper_buf");    
}

// find the z_file context that corresponds to dict_id. It could be possibly a different did_i
// than in the vb - in case this dict_id is new to this vb, but another vb already inserted
// it to z_file
static unsigned mtf_get_z_file_did_i (VariantBlock *vb, DictIdType dict_id)
{
    for (unsigned did_i=0; did_i < vb->z_file->num_dict_ids; did_i++)
        if (dict_id.num == vb->z_file->mtf_ctx[did_i].dict_id.num) {
            //printf ("Inserting new z_file dict_id=%.8s in did_i=%u\n", dict_id_printable (dict_id).id, did_i);
            return did_i;
        }

    // z_file doesn't yet have this dict_id - add it now
    ASSERT (vb->z_file->num_dict_ids+1 < MAX_DICTS, 
            "Error: z_file has more dict_id types than MAX_DICTS=%u", MAX_DICTS);

    vb->z_file->mtf_ctx[vb->z_file->num_dict_ids].did_i   = vb->z_file->num_dict_ids;
    vb->z_file->mtf_ctx[vb->z_file->num_dict_ids].dict_id = dict_id;
    vb->z_file->num_dict_ids++;

    return vb->z_file->num_dict_ids-1;
}

// ZIP only: decide for each ctx whether we will use 8 or 16 bit encoding. Note: different VBs might use different
// encodings for b250 data based on the same dictionary - because of the logic
static void mtf_decide_encodings (MtfContext *vb_ctx, MtfContext *zf_ctx)
{
    if (vb_ctx->mtf.len + vb_ctx->ol_mtf.len == 0) return; // nothing to do

    // calculate enc - the minimum encoding level required for this VB
    Base250Encoding enc = B250_ENC_NONE;

    // note: this is a heuristic - after the dictionary is merged into z_file, mtf.len might end up being
    // more than 250 even in 8 bit, and as a result some base250 might have more than 1 numeral
    if (vb_ctx->b250_section_type == SEC_GENOTYPE_DATA || flag_encode_8)
        // all genotype dictionaries are 8bit - for now (bc we don't have a place in the section headers to 
        // indicate each subfield encoding). On the PIZ size, this is similarly set in piz_uncompress_all_sections()
        enc = B250_ENC_8; 

    else if (flag_encode_16)
        enc = B250_ENC_16; 

    else if (flag_encode_24)
        enc = B250_ENC_24; 

    else if (vb_ctx->mtf.len + vb_ctx->ol_mtf.len <= 250) // this condition is checked only if not flag_encode_16/24
        enc = B250_ENC_8; 

    else if (vb_ctx->mtf.len + vb_ctx->ol_mtf.len <= 62500)
        enc = B250_ENC_16; 

    else
        enc = B250_ENC_24; 

    // calculate the actual encoding to be used - possibly higher than enc, if previous VBs have already used higher
    // (and hence Base250 of previously entered nodes doesn't contain lower bit values). 
    // if enc is higher than the previously highest, then now the minimum encryption is increased
    vb_ctx->encoding /* to be used for this VB */ = zf_ctx->encoding /* minimum for all future VBs */ = MAX (enc, zf_ctx->encoding);
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

// ZIP only: this is called towards the end of compressing one vb - merging its dictionaries into the z_file 
static void mtf_merge_in_vb_ctx_one_dict_id (VariantBlock *vb, unsigned did_i)
{
    MtfContext *vb_ctx = &vb->mtf_ctx[did_i];

    unsigned z_did_i = mtf_get_z_file_did_i (vb, vb_ctx->dict_id);
    MtfContext *zf_ctx = &vb->z_file->mtf_ctx[z_did_i];

    //printf ("Merging dict_id=%.8s into z_file vb_i=%u vb_did_i=%u z_did_i=%u\n", dict_id_printable (vb_ctx->dict_id).id, vb->variant_block_i, did_i, z_did_i);

    if (!buf_is_allocated (&vb_ctx->dict)) return; // nothing yet for this dict_id

    uint32_t start_dict_len = zf_ctx->dict.len;
    uint32_t start_mtf_len  = zf_ctx->mtf.len;

    if (!buf_is_allocated (&zf_ctx->dict)) {

        // first data for this dict (usually, but not always, vb_i=1) - move to zf_ctx and leave overlay

        zf_ctx->b250_section_type = vb_ctx->b250_section_type;
        zf_ctx->dict_section_type = vb_ctx->dict_section_type;
        zf_ctx->dict_id           = vb_ctx->dict_id;
        zf_ctx->encoding          = B250_ENC_NONE; // initialize - mtf_decide_encodings might update this

        mtf_decide_encodings (vb_ctx, zf_ctx); // set vb encoding level, and minimum zf encoding level for all future VBs

        buf_move (vb, &zf_ctx->dict, &vb_ctx->dict);
        buf_set_overlayable (&zf_ctx->dict);
        buf_overlay (&vb_ctx->ol_dict, &zf_ctx->dict, 0,0,0,0);

        buf_move (vb, &zf_ctx->mtf,  &vb_ctx->mtf);
        buf_set_overlayable (&zf_ctx->mtf);
        buf_overlay (&vb_ctx->ol_mtf, &zf_ctx->mtf, 0,0,0,0);

        buf_move (vb, &zf_ctx->hash, &vb_ctx->hash); // vb_ctx no longer needs the hash table - zf_ctx can take it

        // encode in base250 - to be used by zip_generate_genotype_one_section() and zip_generate_b250_section()
        for (unsigned i=0; i < zf_ctx->mtf.len; i++) {
            MtfNode *zf_node = &((MtfNode *)zf_ctx->mtf.data)[i];
            zf_node->word_index = base250_encode (zf_node->word_index.n, zf_ctx->encoding); // note that vb overlays this. also, vb_1 has been sorted so word_index != node_index

            ASSERT (zf_node->word_index.n < zf_ctx->mtf.len, // sanity check
                    "Error: word_index=%u out of bound - mtf.len=%u, in dictionary %.*s", 
                    zf_node->word_index.n, zf_ctx->mtf.len, DICT_ID_LEN, dict_id_printable (zf_ctx->dict_id).id);
        }
    }
    else {

        // set vb encoding level, and possibly increase the minimum zf encoding level for all future VBs
        mtf_decide_encodings (vb_ctx, zf_ctx); 

        // merge in words that are potentially new (but may have been already added by other VBs since we cloned for this VB)
        for (unsigned i=0; i < vb_ctx->mtf.len; i++) {
            MtfNode *vb_node = &((MtfNode *)vb_ctx->mtf.data)[i];
            
            MtfNode *zf_node;
            bool is_new;
            uint32_t zf_node_index = mtf_evaluate_snip (vb, zf_ctx, &vb_ctx->dict.data[vb_node->char_index], vb_node->snip_len, &zf_node, &is_new);
            ASSERT (zf_node_index < zf_ctx->mtf.len, "Error: zf_node_index=%u out of range - len=%i", zf_node_index, vb_ctx->mtf.len);

            // set word_index to be indexing the global dict - to be used by zip_generate_genotype_one_section() and zip_generate_b250_section()
            // note that encoding is private to the vb - different vbs might encoding their b250 of a certain dictionary with
            // different encoding
            if (is_new)
                vb_node->word_index = zf_node->word_index = base250_encode (zf_node_index, zf_ctx->encoding);
            else 
                // a previous VB already already calculated the word index for this node. if it was done by vb_i=1,
                // then it is also resorted and the word_index is no longer the same as the node_index
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
       if (dict_id_is (zf_ctx->dict_id, "GL")) 
            start_dict = gl_optimize_dictionary (vb, &zf_ctx->dict, &((MtfNode *)zf_ctx->mtf.data)[start_mtf_len], start_dict_len, added_words);
     
        zfile_compress_dictionary_data (vb, zf_ctx, added_words, start_dict, added_chars);
    }
}

// ZIP only: merge new words added in this vb into the z_file.mtf_ctx, and compresses dictionaries
// while holding exclusive access to the z_file dictionaries. returns num_dictionary_sections
void mtf_merge_in_vb_ctx (VariantBlock *vb)
{
    mtf_wait_for_my_turn(vb); // we grab the mutex in the sequencial order of VBs

    START_TIMER; // note: careful not to count time spent waiting for the mutex

    // first, all field dictionaries    
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {
        if (!buf_is_allocated (&vb->mtf_ctx[did_i].dict)) continue;

        SectionType dict_sec_type = vb->mtf_ctx[did_i].dict_section_type;

        ASSERT (section_type_is_dictionary(dict_sec_type), "Error: dict_sec_type=%s is not a dictionary section", st_name(dict_sec_type));

        if (dict_sec_type != SEC_INFO_SUBFIELD_DICT && dict_sec_type != SEC_FRMT_SUBFIELD_DICT) 
            mtf_merge_in_vb_ctx_one_dict_id (vb, did_i);
    }

    // second, all the info subfield dictionaries
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++)         
        if (buf_is_allocated (&vb->mtf_ctx[did_i].dict) && 
            vb->mtf_ctx[did_i].dict_section_type == SEC_INFO_SUBFIELD_DICT) 
            mtf_merge_in_vb_ctx_one_dict_id (vb, did_i);

    // third, all the genotype subfield dictionaries
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++)         
        if (buf_is_allocated (&vb->mtf_ctx[did_i].dict) && 
            vb->mtf_ctx[did_i].dict_section_type == SEC_FRMT_SUBFIELD_DICT) 
            mtf_merge_in_vb_ctx_one_dict_id (vb, did_i);

    vb->z_file->next_variant_i_to_merge++;

    // note: vb->z_file->num_dict_ids might be larger than vb->num_dict_ids at this point, for example:
    // vb_i=1 started, z_file is empty, created 20 contexts
    // vb_i=2 started, z_file is empty, created 10 contexts
    // vb_i=1 completes, merges 20 contexts to z_file, which has 20 contexts after
    // vb_i=2 completes, merges 10 contexts, of which 5 (for example) are shared with vb_i=1. Now z_file has 25 contexts after.

    // now, we merge vb->ra_buf into z_file->ra_buf
    random_access_merge_in_vb (vb);

    pthread_mutex_unlock (&vb->z_file->mutex);

    COPY_TIMER (vb->profile.mtf_merge_in_vb_ctx)
}

// gets did_id if the dictionary exists, or returns NIL, if not
int mtf_get_existing_did_i_by_dict_id (VariantBlock *vb, DictIdType dict_id)
{
    //MtfContext *mtf_ctx = vb->z_file->mtf_ctx;

    for (unsigned did_i=0; did_i < vb->z_file->num_dict_ids; did_i++) 
        if (dict_id.num == vb->mtf_ctx[did_i].dict_id.num) return did_i;

    return NIL; // not found
}

// gets did_id if the dictionary exists, and creates a new dictionary if its the first time dict_id is encountered
MtfContext *mtf_get_ctx_by_dict_id (MtfContext *mtf_ctx /* an array */, 
                                    unsigned *num_dict_ids, 
                                    unsigned *num_subfields, // variable to increment if a new context is added
                                    DictIdType dict_id,
                                    SectionType dict_section_type)
{
    // check if we have this dict_id already
    unsigned did_i=0 ; for (; did_i < *num_dict_ids; did_i++) 
        if (dict_id.num == mtf_ctx[did_i].dict_id.num) break;

    MtfContext *ctx = &mtf_ctx[did_i];

    // case: dict_id encountered for this first time - initialize a mtf_ctx
    if (did_i == *num_dict_ids) {

        //printf ("Inserting new vb dict_id=%.8s in did_i=num_dict_ids=%u \n", dict_id_printable (dict_id).id, did_i);
        ASSERT (*num_dict_ids+1 < MAX_DICTS, 
                "Error: number of dictionary types is greater than MAX_DICTS=%u", MAX_DICTS);

        ctx->did_i             = did_i;
        ctx->dict_id           = dict_id;
        ctx->dict_section_type = dict_section_type;
        ctx->b250_section_type = dict_section_type + 1; // the b250 is 1 after the dictionary for all dictionary sections
        ctx->encoding          = ENCRYPTION_TYPE_NONE;         // encoding will be decided at the end of segregate
        mtf_init_iterator (ctx);

        // thread safety: the increment below MUST be AFTER memcpy, bc piz_get_line_subfields
        // might be reading this data at the same time as the piz dispatcher thread adding more dictionaries
        (*num_dict_ids)++; 

        if (num_subfields) (*num_subfields)++;
    }

    ASSERT (ctx->dict_section_type == dict_section_type, "Error: mismatch in dict_id=%.*s dict_section_type: requested %s but in the ctx says: %s",
            DICT_ID_LEN, dict_id_printable (dict_id).id, st_name(dict_section_type), st_name(ctx->dict_section_type));

    return ctx;
}

// PIZ only: this is called by the I/O thread after reading a dictionary section 
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

    ASSERT (section_type_is_dictionary(header->h.section_type),
            "Error: header->h.section_type=%s is not a dictionary section", st_name(header->h.section_type));

    uint32_t num_snips = BGEN32 (header->num_snips);

    zfile_uncompress_section (vb, section_data, &fragment, "fragment", header->h.section_type);

    // special treatment if this is GL - de-optimize
    if (dict_id_is (header->dict_id, "GL"))
        gl_deoptimize_dictionary (fragment.data, fragment.len);

    // in piz, the same did_i is used for z_file and vb contexts, meaning that in vbs there could be
    // a non-contiguous array of contexts (some are missing if not used by this vb)

    MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->z_file->mtf_ctx, &vb->z_file->num_dict_ids, NULL, header->dict_id, header->h.section_type);
    
    // append fragment to dict. If there is no room - old memory is abandoned (so that VBs that are overlaying
    // it continue to work uninterrupted) and a new memory is allocated, where the old dict is joined by the new fragment
    unsigned dict_old_len = ctx->dict.len;
    buf_alloc (vb, &ctx->dict, ctx->dict.len + fragment.len, 2, "mtf_ctx->dict", header->h.section_type);
    buf_set_overlayable (&ctx->dict);

    memcpy (&ctx->dict.data[ctx->dict.len], fragment.data, fragment.len);
    ctx->dict.len += fragment.len;

    // extend word list memory - and calculate the new words. If there is no room - old memory is abandoned 
    // (so that VBs that are overlaying it continue to work uninterrupted) and a new memory is allocated
    buf_alloc (vb, &ctx->word_list, (ctx->word_list.len + num_snips) * sizeof (MtfWord), 2, "mtf_ctx->word_list", header->h.section_type);
    buf_set_overlayable (&ctx->word_list);

    bool is_ref_alt = !strncmp ((char*)dict_id_printable (header->dict_id).id, vcf_field_names[REFALT], MIN (strlen(vcf_field_names[REFALT]+1), DICT_ID_LEN)); // compare inc. \0 terminator

    char *start = fragment.data;
    for (unsigned snip_i=0; snip_i < num_snips; snip_i++) {
        MtfWord *word = &((MtfWord *)ctx->word_list.data)[ctx->word_list.len++];

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
void mtf_overlay_dictionaries_to_vb (VariantBlock *vb)
{
    for (unsigned did_i=0; did_i < MAX_DICTS; did_i++) {
        MtfContext *zf_ctx = &vb->z_file->mtf_ctx[did_i];
        MtfContext *vb_ctx = &vb->mtf_ctx[did_i];

        if (!zf_ctx->dict_id.id[0]) break;

        if (buf_is_allocated (&zf_ctx->dict) && buf_is_allocated (&zf_ctx->word_list)) { 
            
            vb_ctx->did_i             = did_i;
            vb_ctx->dict_id           = zf_ctx->dict_id;
            vb_ctx->b250_section_type = zf_ctx->b250_section_type;
            vb_ctx->dict_section_type = zf_ctx->dict_section_type;
            mtf_init_iterator (vb_ctx);

            buf_overlay (&vb_ctx->dict, &zf_ctx->dict, 0,0,0,0);    
            buf_overlay (&vb_ctx->word_list, &zf_ctx->word_list, 0,0,0,0);

            // count dictionaries of genotype data subfields
            if (dict_id_is_gtdata_subfield (vb_ctx->dict_id)) {
                vb->num_format_subfields++;
                ASSERT (vb->num_format_subfields <= MAX_SUBFIELDS, 
                        "Error: number of subfields in %s exceeds MAX_SUBFIELDS=%u, while reading vb_i=%u", 
                        file_printname (vb->z_file), MAX_SUBFIELDS, vb->variant_block_i);
            }
        }
    }
    vb->num_dict_ids = vb->z_file->num_dict_ids;
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
            node->char_index   = next - ctx->dict.data;
            node->word_index.n = i;

            next += node->snip_len + 1;
        }

        buf_free (&old_dict);
    }
}

// zero all sorters - this is called in case of a re-do of the first VB due to ploidy overflow
void mtf_zero_all_sorters (VariantBlock *vb)
{
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {
        MtfContext *ctx = &vb->mtf_ctx[did_i];
        buf_zero (&ctx->sorter);
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
    buf_free (&ctx->mtf_i);
    buf_free (&ctx->b250);
    ctx->dict_id.num = 0;
    ctx->dict_section_type = ctx->b250_section_type = 0;
    ctx->iterator.next_b250 = NULL;
    ctx->iterator.prev_word_index =0;
}

