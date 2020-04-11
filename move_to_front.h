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
#include "base250.h"
#include "sections.h"

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define WORD_INDEX_MAX_INDEX  0xfffffffcUL // the number just smaller than all the special values below
#define WORD_INDEX_ONE_UP     0xfffffffdUL // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define WORD_INDEX_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

#define NIL ((int32_t)-1)
typedef struct mtfnode_ {
    uint32_t char_index;      // character index into dictionary array
    uint32_t snip_len;        // not including \t terminator present in dictionary array
    Base250 word_index;       // word index into dictionary 
} MtfNode;

typedef struct {
    uint32_t char_index;
    uint32_t snip_len;
} MtfWord;

// used by vblock_i=1 to collect frequency statistics
typedef struct {
    int32_t mtf_i;             // index into MtfContext.mtf
    int32_t count;             // number of times this snip has been encoutered so far
} SorterEnt;

typedef struct {
    const uint8_t *next_b250;  // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    int32_t prev_word_index;   // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initalized to -1)
} SnipIterator;

typedef struct mtfcontext_ {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    unsigned did_i;            // the index of this ctx within the array vb->mtf_ctx
    DictIdType dict_id;        // which dict_id is this MTF dealing with
    SectionType b250_section_type; // section type where the the corresponding b250 data goes
    SectionType dict_section_type; // section type for dictionary statistics (this is a "fake" section type, not one that we write to disk)
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // The buffer of b250 data containing indeces (in b250) to word_list. Used by vardata fields 
                               // and INFO subfields. Not used by FORMAT subfields - those are stored in genotype_data.
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
    uint32_t mtf_len_at_half;  // value of mtf->len after an estimated half of the lines have been segmented
    uint32_t num_lines_at_half;// the number of lines segmented when mtf_len_at_half was set

    Buffer sorter;             // used by the vb_i==1 only of ZIP to sort the dictionary - entries of SorterEnt

    // used by zf_ctx only in ZIP only    
    pthread_mutex_t mutex;     // MtfContext in z_file (only) is protected by a mutex 
    bool mutex_initialized;
    
    // ----------------------------
    // UNZIP only fields
    // ----------------------------
    Buffer word_list;          // PIZ only. word list. an array of MtfWord - referring to the snips in dictionary
    SnipIterator iterator;     // PIZ only: used to iterate on the context, reading one b250 word_index at a time

} MtfContext;

static inline void mtf_init_iterator (MtfContext *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; }

extern uint32_t mtf_evaluate_snip_seg (VBlockP segging_vb, MtfContextP vb_ctx, const char *snip, uint32_t snip_len, MtfNode **node, bool *is_new);
extern uint32_t mtf_get_next_snip (VBlockP vb, MtfContext *ctx, SnipIterator *override_iterator, const char **snip, uint32_t *snip_len, uint32_t vcf_line);
extern int32_t mtf_search_for_node_index (MtfContext *ctx, const char *snip, unsigned snip_len);
extern void mtf_clone_ctx (VBlockP vb);
extern MtfNode *mtf_node_do (const MtfContext *ctx, uint32_t mtf_i, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define mtf_node(ctx, mtf_i, snip_in_dict, snip_len) mtf_node_do(ctx, mtf_i, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void mtf_merge_in_vb_ctx (VBlockP vb);
extern MtfContext *mtf_get_ctx_by_dict_id (MtfContext *mtf_ctx, unsigned *num_dict_ids, uint8_t *num_subfields, DictIdType dict_id, SectionType dict_section_type);
extern uint8_t mtf_get_existing_did_i_by_dict_id (VBlockP vb, DictIdType dict_id);
extern void mtf_integrate_dictionary_fragment (VBlockP vb, char *data);
extern void mtf_overlay_dictionaries_to_vb (VBlockP vb);
extern void mtf_sort_dictionaries_vb_1(VBlockP vb);
extern void mtf_zero_all_sorters (VBlockP vb);

extern void mtf_initialize_mutex (void);
extern void mtf_update_stats (VBlockP vb);
extern void mtf_free_context (MtfContext *ctx);
extern void mtf_destroy_context (MtfContext *ctx);

extern void mtf_vb_1_lock (VBlockP vb);

#define mtf_get_word(ctx, word_index) (&((MtfWord*)(ctx)->word_list.data)[(word_index)])

#endif