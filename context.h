// ------------------------------------------------------------------
//   context.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CONTEXT_INCLUDED
#define CONTEXT_INCLUDED

#include <pthread.h>
#include "genozip.h"
#include "buffer.h"
#include "base250.h"
#include "sections.h"
#include "bit_array.h"
#include "mutex.h"

#define MAX_WORDS_IN_CTX 0x7ffffff0 // limit on nodes.len, word_list.len - partly because hash uses signed int32_t + 2 for singlton using index-2

#define NODE_INDEX_NONE       -1
#define MAX_NODE_INDEX (MAX_WORDS_IN_CTX-1)

#define MAX_WORD_INDEX (MAX_WORDS_IN_CTX-1)
#define WORD_INDEX_NONE       -1
#define WORD_INDEX_ONE_UP     -2 // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   -3 // empty string
#define WORD_INDEX_MISSING_SF -4 // container item missing, remove preceding separator

// Tell PIZ to replace this character by something else (can appear in any part of a snip in a dictionary, or even multiple times in a snip)
// We use characters that cannot appear in a snip - i.e. other than ASCII 32-127, \t (\x9) \n (\xA) \r (\xD)
#define SNIP_SEP                 '\x0'   // Seperator between snips - both for dict and local 
#define SNIP_LOOKUP              '\x1'   // Lookup from local (optionally followed by a snip - interpreted differently by local type, see reconstruct_one_snip)
#define SNIP_OTHER_LOOKUP        '\x2'   // Lookup from local of other dict_id (possibly with length for sequence storage)
#define SNIP_PAIR_LOOKUP         '\x3'   // Lookup from paired file (when using --pair)  
#define SNIP_CONTAINER           '\x4'   // Appears as first character in the SNIP, followed by a specification of a container field
#define SNIP_SELF_DELTA          '\x5'   // The value is a uint32_t which is a result of the last value + the positive or negative textual int32_t value following this character
#define SNIP_OTHER_DELTA         '\x6'   // The value is a uint32_t which is a result of the last value of another field + the delta value. following this char, {DictId dict_id, int32_t delta, bool update_other} in base64)
#define SNIP_PAIR_DELTA          '\x7'   // The value is a uint32_t which is a result of the equivalent value in the paired file + the delta value (when using --pair)
#define SNIP_SPECIAL             '\x8'   // Special algorithm followed by ID of the algorithm 
#define SNIP_REDIRECTION         '\xB'   // Get the data from another dict_id (can be in b250, local...)
#define SNIP_DONT_STORE          '\xC'   // Reconcstruct the following value, but don't store it in last_value (overriding flags.store)

#define DECLARE_SNIP const char *snip=NULL; uint32_t snip_len=0

typedef struct CtxNode {
    CharIndex char_index; // character index into dictionary array
    uint32_t snip_len;    // not including SNIP_SEP terminator present in dictionary array
    int32_t  count;       // number of times this snip has been evaluated in this VB
    Base250  word_index;  // word index into dictionary 
} CtxNode;

typedef struct {
    CharIndex char_index;
    uint32_t snip_len;
} CtxWord;

typedef struct { // initialize with ctx_init_iterator()
    const uint8_t *next_b250;  // Pointer into b250 of the next b250 to be read (must be initialized to NULL)
    WordIndex prev_word_index; // When decoding, if word_index==BASE250_ONE_UP, then make it prev_word_index+1 (must be initalized to -1)
} SnipIterator;

// SIGNED NUMBERS ARE NOT UNTEST YET! NOT USE YET BY ANY SEG
// for signed numbers, we store them in our "interlaced" format rather than standard ISO format 
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// compressing better in an array that contains both positive and negative
#define SAFE_NEGATE(type,n) ((u##type)(-((int64_t)n))) // careful negation to avoid overflow eg -(-128)==0 in int8_t
#define INTERLACE(type,n) ((((type)n) < 0) ? ((SAFE_NEGATE(type,n) << 1) - 1) : (((u##type)n) << 1))
#define DEINTERLACE(signedtype,unum) (((unum) & 1) ? -(signedtype)(((unum)>>1)+1) : (signedtype)((unum)>>1))

typedef union LastValueType { // 64 bit
    int64_t i;
    double f;
} LastValueType;

typedef struct Context {
    // ----------------------------
    // common fields for ZIP & PIZ
    // ----------------------------
    const char name[DICT_ID_LEN+1]; // nul-terminated printable dict_id
    DidIType did_i;            // the index of this ctx within the array vb->contexts
    DidIType st_did_i;         // in --stats, consolidate this context into st_did_i
    LocalType ltype;           // LT_* - type of local data - included in the section header
    struct FlagsCtx flags;     // flags to be included in section header
    struct FlagsCtx pair_flags;// Used if this file is a PAIR_2 - contains ctx->flags of the PAIR_1
    DictId dict_id;            // which dict_id is this MTF dealing with
    Buffer dict;               // tab-delimited list of all unique snips - in this VB that don't exist in ol_dict
    Buffer b250;               // The buffer of b250 data containing indices (in b250) to word_list. 
    Buffer local;              // VB: Data private to this VB that is not in the dictionary
    Buffer pair;               // Used if this file is a PAIR_2 - contains a copy of either b250 or local of the PAIR_1 (if inst.pair_b250 or inst.pair_local is set)
    LastValueType last_value;  // last value of this context (it can be a basis for a delta, used for BAM translation, and other uses)
    int64_t last_delta;        // last delta value calculated
    SnipIterator pair_b250_iter; // Iterator on pair, if it contains b250 data

    // ----------------------------
    // ZIP only fields
    // ----------------------------
    Buffer ol_dict;            // VB: tab-delimited list of all unique snips - overlayed all previous VB dictionaries
                               // zfile: singletons are stored here
    Buffer ol_nodes;           // MTF nodes - overlayed all previous VB dictionaries. char/word indices are into ol_dict.
                               // zfile: nodes of singletons
    Buffer nodes;              // array of CtxNode - in this VB that don't exist in ol_nodes. char/word indices are into dict.
    Buffer node_i;             // contains 32bit indices into the ctx->nodes - this is an intermediate step before generating b250 or genotype_data 
    
    // settings
    Codec lcodec, bcodec;      // codec used to compress local and b250
    Codec lsubcodec_piz;       // piz to decompress with this codec, AFTER decompressing with lcodec

    // ZIP-only instructions NOT written to the genozip file
    bool no_stons;             // don't attempt to move singletons to local (singletons are never moved anyway if ltype!=LT_TEXT)
    bool pair_local;           // this is the 2nd file of a pair - compare vs the first file, and set flags.paired in the header of SEC_LOCAL
    bool pair_b250;            // this is the 2nd file of a pair - compare vs the first file, and set flags.paired in the header of SEC_B250
    bool stop_pairing;         // this is the 2nd file of a pair - don't use SNIP_PAIR_LOOKUP/DELTA anymore until the end of this VB
    bool no_callback;          // don't use LOCAL_GET_LINE_CALLBACK for compressing, despite it being defined
    bool local_param;          // copy local.param to SectionHeaderCtx
    bool no_vb1_sort;          // don't sort the dictionary in ctx_sort_dictionaries_vb_1
    bool local_always;         // always create a local section in zfile, even if it is empty 

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
    uint32_t nodes_len_at_1_3, nodes_len_at_2_3;  // value of nodes->len after an estimated 1/3 + 2/3 of the lines have been segmented
    
    // stats
    uint64_t txt_len;          // How many characters in the txt file are accounted for by snips in this ctx (for stats)
    uint32_t num_singletons;   // True singletons that appeared exactly once in the entire file
    uint32_t num_failed_singletons;// Words that we wrote into local in one VB only to discover later that they're not a singleton, and wrote into the global dict too

    // ----------------------------
    // ZIP in z_file only
    // ----------------------------
    Mutex mutex;               // Context in z_file (only) is protected by a mutex 
    
    // ----------------------------
    // PIZ only fields
    // ----------------------------
    Buffer word_list;          // PIZ only. word list. an array of CtxWord - listing the snips in dictionary
    SnipIterator iterator;     // PIZ only: used to iterate on the context, reading one b250 word_index at a time
    uint32_t next_local;       // PIZ only: iterator on Context.local

    uint32_t last_line_i;      // PIZ only: the last line_i this ctx was encountered

    // PIZ-only instructions
    bool semaphore;            // valid within the context of reconstructing a single line. MUST be reset ahead of completing the line.

    // Container cache 
    Buffer con_cache;       // An array of Container which includes the did_i. Each struct is truncated to used items, followed by prefixes. 
    Buffer con_index;       // Array of uint32_t - index into con_cache. Each item corresponds to word_index (PIZ) or node_index (ZIP)
    Buffer con_len;         // Array of uint16_t - length of item in cache
    
} Context;

#define CTX &vb->contexts
#define Z_CTX &z_file->contexts

#define NEXTLOCAL(type, ctx) (*ENT (type, (ctx)->local, (ctx)->next_local++))
static inline bool PAIRBIT(Context *ctx)      { BitArrayP b = buf_get_bitarray (&ctx->pair);  bool ret = bit_array_get (b, ctx->next_local);                    return ret; } // we do it like this and not in a #define to avoid anti-aliasing warning when casting part of a Buffer structure into a BitArray structure
static inline bool NEXTLOCALBIT(Context *ctx) { BitArrayP b = buf_get_bitarray (&ctx->local); bool ret = bit_array_get (b, ctx->next_local); ctx->next_local++; return ret; }

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75

static inline void ctx_init_iterator (Context *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern WordIndex ctx_evaluate_snip_seg (VBlockP segging_vb, ContextP vb_ctx, const char *snip, uint32_t snip_len, bool *is_new);
extern WordIndex ctx_get_next_snip (VBlockP vb, Context *ctx, bool all_the_same, SnipIterator *override_iterator, const char **snip, uint32_t *snip_len);
extern const char *ctx_peek_next_snip (VBlockP vb, Context *ctx);
extern WordIndex ctx_search_for_word_index (Context *ctx, const char *snip, unsigned snip_len);
extern void ctx_clone (VBlockP vb);
extern CtxNode *ctx_node_vb_do (const Context *ctx, WordIndex node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define ctx_node_vb(ctx, node_index, snip_in_dict, snip_len) ctx_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern CtxNode *ctx_node_zf_do (const Context *ctx, int32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define ctx_node_zf(ctx, node_index, snip_in_dict, snip_len) ctx_node_zf_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void ctx_merge_in_vb_ctx (VBlockP vb);
extern void ctx_commit_codec_to_zf_ctx (VBlockP vb, ContextP vb_ctx, bool is_lcodec);

extern Context *ctx_get_ctx_if_not_found_by_inline (Context *contexts, DataType dt, uint8_t *dict_id_to_did_i_map, uint8_t map_did_i, DidIType *num_contexts, DictId dict_id);

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define ctx_get_ctx(vb,dict_id) ctx_get_ctx_do (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts, (DictId)(dict_id))
static inline Context *ctx_get_ctx_do (Context *contexts, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_contexts, DictId dict_id)
{
    DidIType did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return &contexts[did_i];
    else    
        return ctx_get_ctx_if_not_found_by_inline (contexts, dt, dict_id_to_did_i_map, did_i, num_contexts, dict_id);
}

extern DidIType ctx_get_existing_did_i_if_not_found_by_inline (VBlockP vb, DictId dict_id);
#define ctx_get_existing_did_i(vb,dict_id) ctx_get_existing_did_i_do (vb, (DictId)dict_id, vb->contexts, vb->dict_id_to_did_i_map)
static inline DidIType ctx_get_existing_did_i_do (VBlockP vb, DictId dict_id, Context *contexts, DidIType *dict_id_to_did_i_map)
{
    DidIType did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i == DID_I_NONE || contexts[did_i].dict_id.num == dict_id.num) return did_i;
    return ctx_get_existing_did_i_if_not_found_by_inline (vb, dict_id);
}    
extern ContextP ctx_get_existing_ctx_do (VBlockP vb, DictId dict_id); // returns NULL if context doesn't exist
#define ctx_get_existing_ctx(vb,dict_id) ctx_get_existing_ctx_do ((VBlockP)vb, (DictId)dict_id)

extern void ctx_overlay_dictionaries_to_vb (VBlockP vb);
extern void ctx_sort_dictionaries_vb_1(VBlockP vb);
extern void ctx_verify_field_ctxs_do (VBlockP vb, const char *func, uint32_t code_line);
#define ctx_verify_field_ctxs(vb) ctx_verify_field_ctxs_do(vb, __FUNCTION__, __LINE__);

extern void ctx_update_stats (VBlockP vb);
extern void ctx_free_context (Context *ctx);
extern void ctx_destroy_context (Context *ctx);
extern void ctx_map_aliases (VBlockP vb);
extern CtxNode *ctx_get_node_by_word_index (Context *ctx, WordIndex word_index);
extern const char *ctx_get_snip_by_word_index (const Buffer *word_list, const Buffer *dict, WordIndex word_index, 
                                               const char **snip, uint32_t *snip_len);
extern const char *ctx_get_snip_by_node_index (const Buffer *nodes, const Buffer *dict, WordIndex node_index, 
                                               const char **snip, uint32_t *snip_len);

extern void ctx_initialize_primary_field_ctxs (Context *contexts /* an array */, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_contexts);

typedef enum {DICTREAD_ALL, DICTREAD_CHROM_ONLY, DICTREAD_EXCEPT_CHROM} ReadChromeType;
extern void ctx_read_all_dictionaries (ReadChromeType read_chrom);
extern void ctx_compress_dictionaries (void);
extern void ctx_copy_ref_contigs_to_zf (DidIType dst_did_i, ConstBufferP contigs, ConstBufferP contigs_dict);

extern void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local);

#endif