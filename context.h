// ------------------------------------------------------------------
//   context.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CONTEXT_INCLUDED
#define CONTEXT_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "base250.h"
#include "sections.h"
#include "bit_array.h"
#include "mutex.h"
#include "context_struct.h"
#include "vblock.h"

#define MAX_WORDS_IN_CTX 0x7ffffff0 // limit on nodes.len, word_list.len - partly because hash uses signed int32_t + 2 for singlton using index-2

#define NODE_INDEX_NONE       -1
#define MAX_NODE_INDEX (MAX_WORDS_IN_CTX-1)

#define MAX_WORD_INDEX (MAX_WORDS_IN_CTX-1)
#define WORD_INDEX_NONE       -1
#define WORD_INDEX_ONE_UP     -2 // the value is the one more than the previous value
#define WORD_INDEX_EMPTY_SF   -3 // empty string
#define WORD_INDEX_MISSING_SF -4 // container item missing, remove preceding separator
#define MIN_WORD_INDEX        -4

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
#define SNIP_OTHER_COPY          '\xE'   // Copy the last_txt of another dict_id 

#define DECLARE_SNIP const char *snip=NULL; uint32_t snip_len=0
#define SNIP(len) uint32_t snip_len=(len); char snip[len]

typedef struct CtxNode {
    CharIndex char_index; // character index into dictionary array
    uint32_t snip_len;    // not including SNIP_SEP terminator present in dictionary array
    Base250  word_index;  // word index into dictionary 
} CtxNode;

typedef struct {
    CharIndex char_index;
    uint32_t snip_len;
} CtxWord;

// SIGNED NUMBERS ARE NOT UNTEST YET! NOT USE YET BY ANY SEG
// for signed numbers, we store them in our "interlaced" format rather than standard ISO format 
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// compressing better in an array that contains both positive and negative
#define SAFE_NEGATE(type,n) ((u##type)(-((int64_t)n))) // careful negation to avoid overflow eg -(-128)==0 in int8_t
#define INTERLACE(type,n) ((((type)n) < 0) ? ((SAFE_NEGATE(type,n) << 1) - 1) : (((u##type)n) << 1))
#define DEINTERLACE(signedtype,unum) (((unum) & 1) ? -(signedtype)(((unum)>>1)+1) : (signedtype)((unum)>>1))

#define CTX &vb->contexts
#define Z_CTX &z_file->contexts

#define NEXTLOCAL(type, ctx) (*ENT (type, (ctx)->local, (ctx)->next_local++))
static inline bool PAIRBIT(Context *ctx)      { BitArrayP b = buf_get_bitarray (&ctx->pair);  bool ret = bit_array_get (b, ctx->next_local);                    return ret; } // we do it like this and not in a #define to avoid anti-aliasing warning when casting part of a Buffer structure into a BitArray structure
static inline bool NEXTLOCALBIT(Context *ctx) { BitArrayP b = buf_get_bitarray (&ctx->local); bool ret = bit_array_get (b, ctx->next_local); ctx->next_local++; return ret; }

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75  

#define ctx_node_vb(ctx, node_index, snip_in_dict, snip_len) ctx_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
#define node_word_index(vb,did_i,index) ((index)!=WORD_INDEX_NONE ? ctx_node_vb (&(vb)->contexts[did_i], (index), 0,0)->word_index.n : WORD_INDEX_NONE)

#define last_int(did_i)     contexts[did_i].last_value.i
#define last_index(did_i)   contexts[did_i].last_value.i
#define last_float(did_i)   contexts[did_i].last_value.f
#define last_delta(did_i)   contexts[did_i].last_delta
#define last_txtx(vb, ctx)  ENT (char, (vb)->txt_data, (ctx)->last_txt_index)
#define last_txt(vb, did_i) last_txtx (vb, &(vb)->contexts[did_i])
#define last_txt_len(did_i) contexts[did_i].last_txt_len
#define is_last_txt_valid(ctx) ((ctx)->last_txt_index != INVALID_LAST_TXT_INDEX)

static inline void ctx_init_iterator (Context *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern WordIndex ctx_evaluate_snip_seg (VBlockP segging_vb, ContextP vb_ctx, const char *snip, uint32_t snip_len, bool *is_new);
extern int64_t ctx_decrement_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_increment_count (VBlockP vb, ContextP ctx, WordIndex node_index);

extern WordIndex ctx_get_next_snip (VBlockP vb, Context *ctx, bool all_the_same, bool is_pair, const char **snip, uint32_t *snip_len);

extern WordIndex ctx_search_for_word_index (Context *ctx, const char *snip, unsigned snip_len);
extern void ctx_clone (VBlockP vb);
extern CtxNode *ctx_node_vb_do (const Context *ctx, WordIndex node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
extern CtxNode *ctx_node_zf_do (const Context *ctx, int32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define ctx_node_zf(ctx, node_index, snip_in_dict, snip_len) ctx_node_zf_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void ctx_merge_in_vb_ctx (VBlockP vb);
extern void ctx_commit_codec_to_zf_ctx (VBlockP vb, ContextP vb_ctx, bool is_lcodec);

extern Context *ctx_get_ctx_if_not_found_by_inline (Context *contexts, DataType dt, DidIType *dict_id_to_did_i_map, DidIType map_did_i, DidIType *num_contexts, DictId dict_id);

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define ctx_get_ctx(vb,dict_id) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->dict_id_to_did_i_map, &((VBlockP)(vb))->num_contexts, (DictId)(dict_id))
static inline Context *ctx_get_ctx_do (Context *contexts, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_contexts, DictId dict_id)
{
    DidIType did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return &contexts[did_i];
    else    
        return ctx_get_ctx_if_not_found_by_inline (contexts, dt, dict_id_to_did_i_map, did_i, num_contexts, dict_id);
}

extern DidIType ctx_get_existing_did_i_if_not_found_by_inline (VBlockP vb, DictId dict_id);

static inline DidIType ctx_get_existing_did_i_do (VBlockP vb, DictId dict_id, Context *contexts, DidIType *dict_id_to_did_i_map)
{
    DidIType did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i == DID_I_NONE || contexts[did_i].dict_id.num == dict_id.num) return did_i;
    return ctx_get_existing_did_i_if_not_found_by_inline (vb, dict_id);
}    
#define ctx_get_existing_did_i(vb,dict_id) ctx_get_existing_did_i_do (vb, dict_id, vb->contexts, vb->dict_id_to_did_i_map)

static inline ContextP ctx_get_existing_ctx_do (VBlockP vb, DictId dict_id)  // returns NULL if context doesn't exist
{
    DidIType did_i = ctx_get_existing_did_i(vb, dict_id); 
    return (did_i == DID_I_NONE) ? NULL : &vb->contexts[did_i]; 
}
#define ctx_get_existing_ctx(vb,dict_id) ctx_get_existing_ctx_do ((VBlockP)(vb), (DictId)(dict_id))

extern struct FlagsCtx ctx_get_zf_ctx_flags (DictId dict_id);

extern void ctx_add_new_zf_ctx_from_txtheader (DictId dict_id, TranslatorId luft_translator);

extern void ctx_overlay_dictionaries_to_vb (VBlockP vb);
extern void ctx_sort_dictionaries_vb_1(VBlockP vb);
extern void ctx_verify_field_ctxs_do (VBlockP vb, const char *func, uint32_t code_line);
#define ctx_verify_field_ctxs(vb) ctx_verify_field_ctxs_do(vb, __FUNCTION__, __LINE__);

extern void ctx_update_stats (VBlockP vb);
extern void ctx_free_context (Context *ctx);
extern void ctx_destroy_context (Context *ctx);
extern void ctx_map_aliases (VBlockP vb);
extern CtxNode *ctx_get_node_by_word_index (Context *ctx, WordIndex word_index);
extern const char *ctx_get_snip_by_word_index (const Context *ctx, WordIndex word_index, 
                                               const char **snip, uint32_t *snip_len);
#define ctx_get_words_snip(ctx, word_index) ctx_get_snip_by_word_index ((ctx), (word_index), 0, 0)

extern const char *ctx_get_snip_by_zf_node_index (const Buffer *nodes, const Buffer *dict, WordIndex node_index, 
                                                  const char **snip, uint32_t *snip_len);
#define ctx_get_zf_nodes_snip(ctx, node_index) ctx_get_snip_by_zf_node_index (&(ctx)->nodes, &(ctx)->dict, (node_index), 0, 0)

extern WordIndex ctx_get_word_index_by_snip (const Context *ctx, const char *snip);

extern void ctx_initialize_primary_field_ctxs (Context *contexts /* an array */, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_contexts);

extern void ctx_read_all_dictionaries (void);
extern void ctx_compress_dictionaries (void);
extern void ctx_read_all_counts (void);
extern void ctx_compress_counts (void);
extern const char *ctx_get_snip_with_largest_count (DidIType did_i, int64_t *count);
extern void ctx_build_zf_ctx_from_contigs (DidIType dst_did_i, ConstBufferP contigs, ConstBufferP contigs_dict);

extern void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local);

// returns true if dict_id was *previously* segged on this line, and we stored a valid last_value (int or float)
#define ctx_has_value_in_line_(vb, ctx) ((ctx)->last_line_i == (vb)->line_i)
static inline bool ctx_has_value_in_line_do (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    Context *ctx = ctx_get_existing_ctx (vb, dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_has_value_in_line_(vb, ctx);
}
#define ctx_has_value_in_line(vb, dict_id, p_ctx) ctx_has_value_in_line_do ((VBlockP)(vb), (DictId)(dict_id), (p_ctx))

static inline void ctx_set_last_value_do (VBlockP vb, ContextP ctx, LastValueType last_value)
{
    ctx->last_value    = last_value;
    ctx->last_line_i   = vb->line_i;
}
#define ctx_set_last_value(vb, ctx, last_value) ctx_set_last_value_do ((VBlockP)(vb), (ctx), (last_value))

// returns true if dict_id was *previously* segged on this line (last_value may be valid or not)
#define ctx_encountered_in_line_(vb, ctx) (((ctx)->last_line_i == (vb)->line_i) || ((ctx)->last_line_i == -(int32_t)(vb)->line_i - 1))
static inline bool ctx_encountered_in_line_do (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    Context *ctx = ctx_get_existing_ctx (vb, dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_encountered_in_line_(vb, ctx);
}
#define ctx_encountered_in_line(vb, dict_id, p_ctx) ctx_encountered_in_line_do ((VBlockP)(vb), (DictId)(dict_id), (p_ctx))

#define ctx_set_encountered_in_line(ctx)  /* set encountered if not already ctx_set_last_value */  \
    if ((ctx)->last_line_i != vb->line_i) \
        (ctx)->last_line_i   = -(int32_t)vb->line_i - 1; 

#endif