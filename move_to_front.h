// ------------------------------------------------------------------
//   move_to_front.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef MOVE_TO_FRONT_INCLUDED
#define MOVE_TO_FRONT_INCLUDED

#include <pthread.h>
#include "genozip.h"
#include "buffer.h"
#include "base250.h"
#include "dict_id.h"
#include "section_types.h"

#define MAX_WORDS_IN_CTX 0x7fffffff // limit on mtf.len, word_list.len - partly because hash uses signed int32_t

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define WORD_INDEX_MAX_INDEX  0xfffffffcUL // the number just smaller than all the special values below
#define WORD_INDEX_ONE_UP     0xfffffffdUL // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define WORD_INDEX_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

// Tell PIZ to replace this character by something else (can appear in any part of a snip in a dictionary, or even multiple times in a snip)
#define SNIP_SEP           '\0'   // Seperator between snips - both for dict and local 
#define PIZ_SNIP_SEP       (is_v5_or_above ? SNIP_SEP : '\t') // for use by PIZ
#define SNIP_LOOKUP_UINT32 '\1'   // Lookup from local containing big endian uint32   
#define SNIP_LOOKUP_TEXT   '\2'   // Lookup from local containing snips separated by SNIP_SEP
#define SNIP_VERBTIM       '\3'   // Appears as first character in the SNIP, tell PIZ to copy the remainder of the snip as is without special handing

#ifndef DID_I_NONE // also defined in vblock.h
#define DID_I_NONE   255
#endif
#define NIL ((int32_t)-1)
typedef struct MtfNode {
    uint32_t char_index;      // character index into dictionary array
    uint32_t snip_len;        // not including SNIP_SEP terminator present in dictionary array
    Base250  word_index;      // word index into dictionary 
    int32_t  count;           // number of times this snip has been evaluated in this VB
} MtfNode;

typedef struct {
    uint32_t char_index;
    uint32_t snip_len;
} MtfWord;

typedef struct { // initialize with mtf_init_iterator()
    const uint8_t *next_b250;  // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    int32_t prev_word_index;   // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initalized to -1)
} SnipIterator;

typedef struct MtfContext {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    unsigned did_i;            // the index of this ctx within the array vb->mtf_ctx
    DictIdType dict_id;        // which dict_id is this MTF dealing with
    const char name[DICT_ID_LEN+1]; // null-terminated printable dict_id
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // The buffer of b250 data containing indeces (in b250) to word_list. 
    Buffer local;              // VB only (not z_file): Data private to this VB that is not in the dictionary
    
    // ----------------------------
    // ZIP only fields
    // ----------------------------
    Buffer ol_dict;            // tab-delimited list of all unique snips - overlayed all previous VB dictionaries
    Buffer ol_mtf;             // MTF nodes - overlayed all previous VB dictionaries. char/word indeces are into ol_dict.
    Buffer mtf;                // array of MtfNode - in this VB that don't exist in ol_mtf. char/word indeces are into dict.
    Buffer mtf_i;              // contains 32bit indeces into the ctx->mtf - this is an intermediate step before generating b250 or genotype_data 
    
    // hash stuff 
    Buffer local_hash;         // hash table for entries added by this VB that are not yet in the global (until merge_number)
                               // obtained by hash function hash(snip) and the rest of linked to them by linked list
    uint32_t local_hash_prime; // prime number - size of the core (without extensions) has table 
    int32_t num_new_entries_prev_merged_vb; // zf_ctx: updated in every merge - how many new words happened in this VB
                               // vb_ctx: copied from zf_ctx during clone, and used to initialize the size of local_hash
                               //         0 means no VB merged yet with this. if a previous vb had 0 new words, it will still be 1.
    Buffer global_hash;        // global hash table that is populated during merge in zf_ctx and is overlayed to vb_ctx during clone.
    uint32_t global_hash_prime; // prime number - size of the core (without extensions) has table 

    uint32_t merge_num;        // in vb_ctx: the merge_num when global_hash was cloned. only entries with merge_num <= this number 
                               // are valid. other entries may be added by later merges and should be ignored.
                               // in zf_ctx: incremented with every merge into this ctx.
    // the next 2 are used in merge to set the size of the global hash table, when the first vb to create a ctx does so
    uint32_t mtf_len_at_1_3, mtf_len_at_2_3;  // value of mtf->len after an estimated 1/3 + 2/3 of the lines have been segmented

    // used by zf_ctx only in ZIP only    
    pthread_mutex_t mutex;     // MtfContext in z_file (only) is protected by a mutex 
    bool mutex_initialized;
    
    uint64_t txt_len;          // How many characters in the txt file are accounted for by snips in this ctx (for stats)
    
    // ----------------------------
    // PIZ only fields
    // ----------------------------
    Buffer word_list;          // PIZ only. word list. an array of MtfWord - listing the snips in dictionary
    SnipIterator iterator;     // PIZ only: used to iterate on the context, reading one b250 word_index at a time
    uint32_t next_local;       // PIZ only: iterator on MtfContext.local

} MtfContext;

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75

static inline void mtf_init_iterator (MtfContext *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; }

extern uint32_t mtf_evaluate_snip_seg (VBlockP segging_vb, MtfContextP vb_ctx, const char *snip, uint32_t snip_len, bool *is_new);
extern uint32_t mtf_get_next_snip (VBlockP vb, MtfContext *ctx, SnipIterator *override_iterator, const char **snip, uint32_t *snip_len, uint32_t vcf_line);
extern int32_t mtf_search_for_word_index (MtfContext *ctx, const char *snip, unsigned snip_len);
extern void mtf_clone_ctx (VBlockP vb);
extern MtfNode *mtf_node_do (const MtfContext *ctx, uint32_t mtf_i, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define mtf_node(ctx, mtf_i, snip_in_dict, snip_len) mtf_node_do(ctx, mtf_i, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void mtf_merge_in_vb_ctx (VBlockP vb);

extern MtfContext *mtf_get_ctx_if_not_found_by_inline (MtfContext *mtf_ctx, uint8_t *dict_id_to_did_i_map, uint8_t map_did_i, unsigned *num_dict_ids, uint8_t *num_subfields, DictIdType dict_id);

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define mtf_get_ctx(vb,dict_id) mtf_get_ctx_do (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, (dict_id))
#define mtf_get_ctx_sf(vb,num_subfields,dict_id) mtf_get_ctx_do (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, (num_subfields), (dict_id))
static inline MtfContext *mtf_get_ctx_do (MtfContext *mtf_ctx, uint8_t *dict_id_to_did_i_map, unsigned *num_dict_ids, uint8_t *num_subfields, DictIdType dict_id)
{
    uint8_t did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && mtf_ctx[did_i].dict_id.num == dict_id.num) 
        return &mtf_ctx[did_i];
    else    
        return mtf_get_ctx_if_not_found_by_inline (mtf_ctx, dict_id_to_did_i_map, did_i, num_dict_ids, num_subfields, dict_id);
}

extern uint8_t mtf_get_existing_did_i_by_dict_id (DictIdType dict_id);

extern void mtf_integrate_dictionary_fragment (VBlockP vb, char *data);
extern void mtf_overlay_dictionaries_to_vb (VBlockP vb);
extern void mtf_sort_dictionaries_vb_1(VBlockP vb);
extern void mtf_verify_field_ctxs_do (VBlockP vb, const char *func, uint32_t code_line);
#define mtf_verify_field_ctxs(vb) mtf_verify_field_ctxs_do(vb, __FUNCTION__, __LINE__);

extern void mtf_initialize_for_zip (void);
extern void mtf_update_stats (VBlockP vb);
extern void mtf_free_context (MtfContext *ctx);
extern void mtf_destroy_context (MtfContext *ctx);

extern void mtf_vb_1_lock (VBlockP vb);
extern MtfNode *mtf_get_node_by_word_index (MtfContext *ctx, uint32_t word_index);
extern void mtf_initialize_primary_field_ctxs (MtfContext *mtf_ctx /* an array */, DataType dt, uint8_t *dict_id_to_did_i_map, unsigned *num_dict_ids);

#endif