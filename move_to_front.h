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
#include "section_types.h"
#include "data_types.h"

#define MAX_WORDS_IN_CTX 0x7ffffff0 // limit on mtf.len, word_list.len - partly because hash uses signed int32_t + 2 for singlton using index-2

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define WORD_INDEX_MAX_INDEX  0xfffffffcUL // the number just smaller than all the special values below
#define WORD_INDEX_ONE_UP     0xfffffffdUL // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define WORD_INDEX_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

// Tell PIZ to replace this character by something else (can appear in any part of a snip in a dictionary, or even multiple times in a snip)
#define SNIP_SEP                 '\0'   // Seperator between snips - both for dict and local 
#define PIZ_SNIP_SEP             (is_v5_or_above ? SNIP_SEP : '\t') // PIZ_SNIP_SEP must be used instead of SNIP_SEP, if there's any chance the code might be executed for PIZ of a VCF of an older version.
#define SNIP_LOOKUP              '\1'   // Lookup from local 
#define SNIP_OTHER_LOOKUP        '\2'   // Lookup from local of other dict_id (possibly with length for sequence storage)
#define SNIP_STRUCTURED          '\3'   // Appears as first character in the SNIP, followed by a specification of a structured field
#define SNIP_SELF_DELTA          '\4'   // The value is a uint32_t which is a result of the last value + the positive or negative textual int32_t value following this character
#define SNIP_OTHER_DELTA         '\5'   // The value is a uint32_t which is a result of the last value of another field + the delta value. following this char, {DictIdType dict_id, int32_t delta, bool update_other} in base64)
#define SNIP_SPECIAL             '\6'   // Special algorithm followed by ID of the algorithm 
#define SNIP_REDIRECTION         '\7'   // Get the data from another dict_id (can be in b250, local...)

// structured snip: it starts with SNIP_STRUCTURED, following by a base64 of a big endian Structured
#pragma pack(1)
#define STRUCTURED_DROP_LAST_SEP_OF_LAST_ELEMENT 0x01
#define STRUCTURED_MAX_REPEATS 4294967294UL // one less than maxuint32 to make it easier to loop with st.repeats without overflow 
#define STRUCTURED_MAX_PREFIXES_LEN 1000 // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="

typedef struct StructuredItem {
        DictIdType dict_id;  
        uint8_t did_i;    // Used only in PIZ, must remain DID_I_NONE in ZIP
        char seperator;
} StructuredItem;

typedef struct Structured {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // 1 to MAX_STRUCTURED_ITEMS
    uint8_t flags;
    char repsep[2];       // repeat seperator - two bytes that appear at the end of each repeat (ignored if 0)
    StructuredItem items[MAX_SUBFIELDS];
} Structured;

#pragma pack()
#define sizeof_structured(st) (sizeof(st) - sizeof((st).items) + (st).num_items * sizeof((st).items[0]))

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

// these values and flags are part of the file format (in SectionHeaderCtx.flags) - so values cannot be changed easily
// CTX_LT_* values are consistent with BAM optional 'B' types (and extend them)
#define CTX_LT_TEXT        0
#define CTX_LT_INT8        1    
#define CTX_LT_UINT8       2
#define CTX_LT_INT16       3
#define CTX_LT_UINT16      4
#define CTX_LT_INT32       5
#define CTX_LT_UINT32      6 
#define CTX_LT_INT64       7    // ffu
#define CTX_LT_UINT64      8    // ffu
#define CTX_LT_FLOAT32     9    // ffu
#define CTX_LT_FLOAT64     10   // ffu
#define CTX_LT_SEQUENCE    11   // length of data extracted is determined by vb->seq_len
#define NUM_CTX_LT         12
extern const char ctx_lt_to_sam_map[NUM_CTX_LT];
extern const int ctx_lt_sizeof_one[NUM_CTX_LT];
extern const bool ctx_lt_is_signed[NUM_CTX_LT];
extern const int64_t ctx_lt_min[NUM_CTX_LT], ctx_lt_max[NUM_CTX_LT];

#define CTX_FL_NO_STONS    0x01 // don't attempt to move singletons to local (singletons are never moved anyway if ltype!=CTX_LT_TEXT)
#define CTX_FL_LOCAL_LZMA  0x02 // compress local with lzma
#define CTX_FL_STORE_VALUE 0x04 // the values of this ctx are uint32_t, and are a basis for a delta calculation (by this field or another one)
#define CTX_FL_STRUCTURED  0x08 // snips usually contain Structured
#define CTX_FL_POS         0x03 // A POS field that stores a delta vs. a different field
#define CTX_FL_POS_BASE    0x07 // A POS field that is the base for delta calculations (with itself and/or other fields)
#define CTX_FL_ID          0x03 // An ID field that is split between a numeric component in local and a textual component in b250

typedef struct MtfContext {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    const char name[DICT_ID_LEN+1]; // null-terminated printable dict_id
    uint8_t did_i;             // the index of this ctx within the array vb->contexts
    uint8_t ltype;        // CTX_LT_*
    uint8_t flags;             // CTX_*
    DictIdType dict_id;        // which dict_id is this MTF dealing with
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // The buffer of b250 data containing indeces (in b250) to word_list. 
    Buffer local;              // VB: Data private to this VB that is not in the dictionary

    // ----------------------------
    // ZIP only fields
    // ----------------------------
    Buffer ol_dict;            // VB: tab-delimited list of all unique snips - overlayed all previous VB dictionaries
                               // zfile: singletons are stored here
    Buffer ol_mtf;             // MTF nodes - overlayed all previous VB dictionaries. char/word indeces are into ol_dict.
                               // zfile: nodes of singletons
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
    
    // stats
    uint64_t txt_len;          // How many characters in the txt file are accounted for by snips in this ctx (for stats)
    uint32_t num_singletons;   // True singletons that appeared exactly once in the entire file
    uint32_t num_failed_singletons;// Words that we wrote into local in one VB only to discover later that they're not a singleton, and wrote into the global dict too

    // ----------------------------
    // ZIP in z_file only
    // ----------------------------
    pthread_mutex_t mutex;     // MtfContext in z_file (only) is protected by a mutex 
    bool mutex_initialized;
    
    // ----------------------------
    // PIZ only fields
    // ----------------------------
    Buffer word_list;          // PIZ only. word list. an array of MtfWord - listing the snips in dictionary
    SnipIterator iterator;     // PIZ only: used to iterate on the context, reading one b250 word_index at a time
    uint32_t next_local;       // PIZ only: iterator on MtfContext.local

    uint32_t last_line_i;      // PIZ only: the last line_i this ctx was encountered
    int64_t last_value;        // PIZ only: last value from which to conduct a delta. 
    int64_t last_delta;        // PIZ only: last delta value calculated
} MtfContext;

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75

static inline void mtf_init_iterator (MtfContext *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern uint32_t mtf_evaluate_snip_seg (VBlockP segging_vb, MtfContextP vb_ctx, const char *snip, uint32_t snip_len, bool *is_new);
extern uint32_t mtf_get_next_snip (VBlockP vb, MtfContext *ctx, SnipIterator *override_iterator, const char **snip, uint32_t *snip_len);
extern int32_t mtf_search_for_word_index (MtfContext *ctx, const char *snip, unsigned snip_len);
extern void mtf_clone_ctx (VBlockP vb);
extern MtfNode *mtf_node_vb_do (const MtfContext *ctx, uint32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define mtf_node_vb(ctx, node_index, snip_in_dict, snip_len) mtf_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern MtfNode *mtf_node_zf_do (const MtfContext *ctx, int32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define mtf_node_zf(ctx, node_index, snip_in_dict, snip_len) mtf_node_zf_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void mtf_merge_in_vb_ctx (VBlockP vb);

extern MtfContext *mtf_get_ctx_if_not_found_by_inline (MtfContext *contexts, DataType dt, uint8_t *dict_id_to_did_i_map, uint8_t map_did_i, unsigned *num_dict_ids, DictIdType dict_id);

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define mtf_get_ctx(vb,dict_id) mtf_get_ctx_do (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_dict_ids, (dict_id))

static inline MtfContext *mtf_get_ctx_do (MtfContext *contexts, DataType dt, uint8_t *dict_id_to_did_i_map, unsigned *num_dict_ids, DictIdType dict_id)
{
    uint8_t did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return &contexts[did_i];
    else    
        return mtf_get_ctx_if_not_found_by_inline (contexts, dt, dict_id_to_did_i_map, did_i, num_dict_ids, dict_id);
}

extern uint8_t mtf_get_existing_did_i (VBlockP vb, DictIdType dict_id);
extern uint8_t mtf_get_existing_did_i_from_z_file (DictIdType dict_id);

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
extern void mtf_initialize_primary_field_ctxs (MtfContext *contexts /* an array */, DataType dt, uint8_t *dict_id_to_did_i_map, unsigned *num_dict_ids);

#endif