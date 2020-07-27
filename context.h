// ------------------------------------------------------------------
//   context.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef MOVE_TO_FRONT_INCLUDED
#define MOVE_TO_FRONT_INCLUDED

#include <pthread.h>
#include "genozip.h"
#include "buffer.h"
#include "base250.h"
#include "sections.h"
#include "data_types.h"

#define MAX_WORDS_IN_CTX 0x7ffffff0 // limit on mtf.len, word_list.len - partly because hash uses signed int32_t + 2 for singlton using index-2

// fake mtf index values that go into genotype_data after segregation if subfields are missing
#define WORD_INDEX_MAX_INDEX  0xfffffffcUL // the number just smaller than all the special values below
#define WORD_INDEX_ONE_UP     0xfffffffdUL // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   0xfffffffeUL // subfield is missing, terminating : present
#define WORD_INDEX_MISSING_SF 0xffffffffUL // subfield is missing at end of cell, no :

// Tell PIZ to replace this character by something else (can appear in any part of a snip in a dictionary, or even multiple times in a snip)
// We use characters that cannot appear in a snip - i.e. other than ASCII 32-127, \t (\x9) \n (\xA) \r (\xD)
#define SNIP_SEP                 '\x0'   // Seperator between snips - both for dict and local 
#define SNIP_LOOKUP              '\x1'   // Lookup from local (optionally followed by a snip - interpreted differently by local type, see piz_reconstruct_one_snip)
#define SNIP_OTHER_LOOKUP        '\x2'   // Lookup from local of other dict_id (possibly with length for sequence storage)
#define SNIP_PAIR_LOOKUP         '\x3'   // Lookup from paired file (when using --pair)  
#define SNIP_STRUCTURED          '\x4'   // Appears as first character in the SNIP, followed by a specification of a structured field
#define SNIP_SELF_DELTA          '\x5'   // The value is a uint32_t which is a result of the last value + the positive or negative textual int32_t value following this character
#define SNIP_OTHER_DELTA         '\x6'   // The value is a uint32_t which is a result of the last value of another field + the delta value. following this char, {DictId dict_id, int32_t delta, bool update_other} in base64)
#define SNIP_PAIR_DELTA          '\x7'   // The value is a uint32_t which is a result of the equivalent value in the paired file + the delta value (when using --pair)
#define SNIP_SPECIAL             '\x8'   // Special algorithm followed by ID of the algorithm 
#define SNIP_REDIRECTION         '\xB'   // Get the data from another dict_id (can be in b250, local...)
#define SNIP_DONT_STORE          '\xC'   // Reconcstruct the following value, but don't store it in last_value (overriding CTX_FL_STORE_INT)

// structured snip: it starts with SNIP_STRUCTURED, following by a base64 of a big endian Structured
#pragma pack(1)
#define STRUCTURED_DROP_LAST_SEP_OF_LAST_ELEMENT 0x01
#define STRUCTURED_MAX_REPEATS 4294967294UL // one less than maxuint32 to make it easier to loop with st.repeats without overflow 
#define STRUCTURED_MAX_PREFIXES_LEN 1000 // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="

typedef struct StructuredItem {
        DictId dict_id;  
        DidIType did_i;    // Used only in PIZ, must remain DID_I_NONE in ZIP
        char seperator[2];
        uint8_t ffu;
} StructuredItem;

typedef struct Structured {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // 1 to MAX_STRUCTURED_ITEMS
    uint8_t flags;
    char repsep[2];       // repeat seperator - two bytes that appear at the end of each repeat (ignored if 0)
    StructuredItem items[MAX_SUBFIELDS];
} Structured;

// identical to Structured, but with only one item
typedef struct MiniStructured {
    uint32_t repeats;     // number of "repeats" (array elements)
    uint8_t num_items;    // must be 1
    uint8_t flags;
    char repsep[2];       // repeat seperator - two bytes that appear at the end of each repeat (ignored if 0)
    StructuredItem items[1];
} MiniStructured;

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

// flags written to the genozip file (header.h.flags)
#define CTX_FL_STORE_INT    0x01 // the values of this ctx.local are uint32_t, and should be stored, eg because they are a basis for a delta calculation (by this field or another one)
#define CTX_FL_STORE_FLOAT  0x02 // the values of this ctx.local are float
#define CTX_FL_PAIRED       0x04 // appears in header.flags of b250 or local sections, indicating to PIZ that the the same section from the same vb of the previous (paired) file should be read from disk too
#define CTX_FL_STRUCTURED   0x08 // snips usually contain Structured
#define CTX_FL_COPY_PARAM   0x10 // appears in header.flags of b250 or local sections, that ctx.b250/local.param should be copied from SectionHeaderCtx.param

#define ctx_is_store(ctx, store_flag) (((ctx)->flags & 0x3) == (store_flag))

// ZIP-only instructions NOT written to the genozip file
#define CTX_INST_NO_STONS    0x01 // don't attempt to move singletons to local (singletons are never moved anyway if ltype!=LT_TEXT)
#define CTX_INST_PAIR_LOCAL  0x02 // this is the 2nd file of a pair - compare vs the first file, and set CTX_FL_PAIRED in the header of SEC_LOCAL
#define CTX_INST_PAIR_B250   0x04 // this is the 2nd file of a pair - compare vs the first file, and set CTX_FL_PAIRED in the header of SEC_B250
#define CTX_INST_NO_CALLBACK 0x08 // don't use LOCAL_GET_LINE_CALLBACK for compressing, despite it being defined
#define CTX_INST_LOCAL_PARAM 0x10 // copy local.param to SectionHeaderCtx

typedef union {
    int64_t i;
    double d;
} LastValueType;

typedef struct Context {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    const char name[DICT_ID_LEN+1]; // null-terminated printable dict_id
    DidIType did_i;            // the index of this ctx within the array vb->contexts
    LocalType ltype;           // LT_* - type of local data - included in the section header
    int8_t flags;              // CTX_FL_* - flags to be included in section header (8 bits)
    DictId dict_id;            // which dict_id is this MTF dealing with
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // The buffer of b250 data containing indeces (in b250) to word_list. 
    Buffer local;              // VB: Data private to this VB that is not in the dictionary
    Buffer pair;               // Used if this file is a PAIR_2 - contains a copy of either b250 or local of the PAIR_1 (if CTX_INST_PAIR_B250 or CTX_INST_PAIR_LOCAL is set)
    SnipIterator pair_b250_iter; // Iterator on pair, if it contains b250 data

    // ----------------------------
    // ZIP only fields
    // ----------------------------
    Buffer ol_dict;            // VB: tab-delimited list of all unique snips - overlayed all previous VB dictionaries
                               // zfile: singletons are stored here
    Buffer ol_mtf;             // MTF nodes - overlayed all previous VB dictionaries. char/word indeces are into ol_dict.
                               // zfile: nodes of singletons
    Buffer mtf;                // array of MtfNode - in this VB that don't exist in ol_mtf. char/word indeces are into dict.
    Buffer mtf_i;              // contains 32bit indeces into the ctx->mtf - this is an intermediate step before generating b250 or genotype_data 
    
    // settings
    CompressionAlg local_comp; // algorithm used to compress local
    uint8_t inst;              // instructions for seg/zip - ORed CTX_INST_ values

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
    pthread_mutex_t mutex;     // Context in z_file (only) is protected by a mutex 
    bool mutex_initialized;
    
    // ----------------------------
    // PIZ only fields
    // ----------------------------
    Buffer word_list;          // PIZ only. word list. an array of MtfWord - listing the snips in dictionary
    SnipIterator iterator;     // PIZ only: used to iterate on the context, reading one b250 word_index at a time
    uint32_t next_local;       // PIZ only: iterator on Context.local

    uint32_t last_line_i;      // PIZ only: the last line_i this ctx was encountered
    LastValueType last_value;  // PIZ only: last value from which to conduct a delta. 
    int64_t last_delta;        // PIZ only: last delta value calculated
} Context;

#define NEXTLOCAL(type, ctx) (*ENT (type, (ctx)->local, (ctx)->next_local++))
static inline bool PAIRBIT(Context *ctx)      { BitArray *b = buf_get_bitarray (&ctx->pair);  bool ret = bit_array_get (b, ctx->next_local);                    return ret; } // we do it like this and not in a #define to avoid anti-aliasing warning when casting part of a Buffer structure into a BitArray structure
static inline bool NEXTLOCALBIT(Context *ctx) { BitArray *b = buf_get_bitarray (&ctx->local); bool ret = bit_array_get (b, ctx->next_local); ctx->next_local++; return ret; }

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75

static inline void mtf_init_iterator (Context *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern uint32_t mtf_evaluate_snip_seg (VBlockP segging_vb, ContextP vb_ctx, const char *snip, uint32_t snip_len, bool *is_new);
extern uint32_t mtf_get_next_snip (VBlockP vb, Context *ctx, SnipIterator *override_iterator, const char **snip, uint32_t *snip_len);
extern const char *mtf_peek_next_snip (VBlockP vb, Context *ctx);
extern int32_t mtf_search_for_word_index (Context *ctx, const char *snip, unsigned snip_len);
extern void mtf_clone_ctx (VBlockP vb);
extern MtfNode *mtf_node_vb_do (const Context *ctx, uint32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define mtf_node_vb(ctx, node_index, snip_in_dict, snip_len) mtf_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern MtfNode *mtf_node_zf_do (const Context *ctx, int32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define mtf_node_zf(ctx, node_index, snip_in_dict, snip_len) mtf_node_zf_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void mtf_merge_in_vb_ctx (VBlockP vb);

extern Context *mtf_get_ctx_if_not_found_by_inline (Context *contexts, DataType dt, uint8_t *dict_id_to_did_i_map, uint8_t map_did_i, DidIType *num_dict_ids, DictId dict_id);

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define mtf_get_ctx(vb,dict_id) mtf_get_ctx_do (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_dict_ids, (DictId)(dict_id))

static inline Context *mtf_get_ctx_do (Context *contexts, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_dict_ids, DictId dict_id)
{
    DidIType did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return &contexts[did_i];
    else    
        return mtf_get_ctx_if_not_found_by_inline (contexts, dt, dict_id_to_did_i_map, did_i, num_dict_ids, dict_id);
}

extern DidIType mtf_get_existing_did_i (VBlockP vb, DictId dict_id);
extern DidIType mtf_get_existing_did_i_from_z_file (DictId dict_id);
extern ContextP mtf_get_existing_ctx (VBlockP vb, DictId dict_id); // returns NULL if context doesn't exist

extern void mtf_integrate_dictionary_fragment (VBlockP vb, char *data);
extern void mtf_overlay_dictionaries_to_vb (VBlockP vb);
extern void mtf_sort_dictionaries_vb_1(VBlockP vb);
extern void mtf_verify_field_ctxs_do (VBlockP vb, const char *func, uint32_t code_line);
#define mtf_verify_field_ctxs(vb) mtf_verify_field_ctxs_do(vb, __FUNCTION__, __LINE__);

extern void mtf_initialize_for_zip (void);
extern void mtf_update_stats (VBlockP vb);
extern void mtf_free_context (Context *ctx);
extern void mtf_destroy_context (Context *ctx);

extern void mtf_vb_1_lock (VBlockP vb);
extern MtfNode *mtf_get_node_by_word_index (Context *ctx, uint32_t word_index);
extern const char *mtf_get_snip_by_word_index (const Buffer *word_list, const Buffer *dict, int32_t word_index, 
                                               const char **snip, uint32_t *snip_len);

extern void mtf_initialize_primary_field_ctxs (Context *contexts /* an array */, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_dict_ids);

#endif