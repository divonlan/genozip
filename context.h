// ------------------------------------------------------------------
//   context.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef CONTEXT_INCLUDED
#define CONTEXT_INCLUDED

#include "genozip.h"
#include "base250.h"
#include "sections.h"
#include "bit_array.h"
#include "mutex.h"
#include "context_struct.h"
#include "vblock.h"
#include "segconf.h"

#define MAX_WORDS_IN_CTX 0x7ffffff0 // limit on nodes.len, word_list.len - partly because hash uses signed int32_t + 2 for singlton using index-2

#define NODE_INDEX_NONE    -1
#define MAX_NODE_INDEX (MAX_WORDS_IN_CTX-1)

#define MAX_WORD_INDEX (MAX_WORDS_IN_CTX-1)
#define WORD_INDEX_NONE    -1
#define WORD_INDEX_ONE_UP  -2 // the value is the one more than the previous value
#define WORD_INDEX_EMPTY   -3 // empty string
#define WORD_INDEX_MISSING -4 // container item missing, remove preceding separator
#define MIN_WORD_INDEX     -4

// Tell PIZ to replace this character by something else (can appear in any part of a snip in a dictionary, or even multiple times in a snip)
// We use characters that cannot appear in a snip - i.e. other than ApCII 32-127, \t (\x9) \n (\xA) \r (\xD)
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
#define SNIP_COPY                '\xE'   // Copy the last_txt of dict_id (same or other)
#define SNIP_DUAL                '\xF'   // A snip containing two snips separated by a SNIP_DUAL - for Primary and Luft reconstruction respectively
#define SNIP_LOOKBACK            '\x10'  // Copy an earlier snip in the same context. Snip is dict_id from which to take the lookback offset, and an optional delta to be applied to the retrieved numeric value. note: line number of the previous snip is variable, but its offset back is fixed (introduced 12.0.41)
#define SNIP_COPY_BUDDY          '\x11'  // Copy a snip on an earlier "buddy" line in the same or another context (note: offset back to the previous snip is variable, but its line number is fixed) (introduced 12.0.41)
#define NUM_SNIP_CODES           18

#define SNIP_CODES { "SNIP_SEP", "SNIP_LOOKUP", "SNIP_OTHER_LOOKUP", "SNIP_PAIR_LOOKUP",\
                     "SNIP_CONTAINER", "SNIP_SELF_DELTA", "SNIP_OTHER_DELTA", "SNIP_PAIR_DELTA", \
                     "SNIP_SPECIAL", "<TAB>", "<NL>", "SNIP_REDIRECTION",\
                     "SNIP_DONT_STORE", "<LF>", "SNIP_COPY", "SNIP_DUAL", "SNIP_LOOKBACK" }
#define SNIP(len) uint32_t snip_len=(len); char snip[len]

typedef struct CtxNode {
    CharIndex char_index; // character index into dictionary array
    uint32_t snip_len;    // not including \0 terminator present in dictionary array
    Base250  word_index;  // word index into dictionary 
} CtxNode;

typedef struct {
    CharIndex char_index;
    uint32_t snip_len;
} CtxWord;
#define WORD_IN_TXT_DATA(snip) ((CtxWord){ .char_index = ENTNUM (vb->txt_data, snip), .snip_len = snip##_len }) // get coordinates in txt_data

// SIGNED NUMBERS ARE NOT UNTEST YET! NOT USE YET BY ANY SEG
// for signed numbers, we store them in our "interlaced" format rather than standard ISO format 
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// compressing better in an array that contains both positive and negative
#define SAFE_NEGATE(type,n) ((u##type)(-((int64_t)n))) // careful negation to avoid overflow eg -(-128)==0 in int8_t
#define INTERLACE(type,n) ((((type)n) < 0) ? ((SAFE_NEGATE(type,n) << 1) - 1) : (((u##type)n) << 1))
#define DEINTERLACE(signedtype,unum) (((unum) & 1) ? -(signedtype)(((unum)>>1)+1) : (signedtype)((unum)>>1))

#define NEXTLOCAL(type, ctx) (*ENT (type, (ctx)->local, (ctx)->next_local++))
#define PEEKNEXTLOCAL(type, ctx, offset) (*ENT (type, (ctx)->local, (ctx)->next_local + offset))
static inline bool PAIRBIT(Context *ctx)      { BitArrayP b = buf_get_bitarray (&ctx->pair);  bool ret = bit_array_get (b, ctx->next_local);                    return ret; } // we do it like this and not in a #define to avoid anti-aliasing warning when casting part of a Buffer structure into a BitArray structure
static inline bool NEXTLOCALBIT(Context *ctx) { BitArrayP b = buf_get_bitarray (&ctx->local); bool ret = bit_array_get (b, ctx->next_local); ctx->next_local++; return ret; }

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75  

#define ctx_node_vb(ctx, node_index, snip_in_dict, snip_len) ctx_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
#define node_word_index(vb,did_i,index) ((index)!=WORD_INDEX_NONE ? ctx_node_vb (&(vb)->contexts[did_i], (index), 0,0)->word_index.n : WORD_INDEX_NONE)

#define CTX(did_i)          (&vb->contexts[did_i])
#define ZCTX(did_i)         (&z_file->contexts[did_i])
#define last_int(did_i)     contexts[did_i].last_value.i
#define last_index(did_i)   contexts[did_i].last_value.i
#define last_float(did_i)   contexts[did_i].last_value.f
#define last_delta(did_i)   contexts[did_i].last_delta
#define last_txtx(vb, ctx)  ENT (char, (vb)->txt_data, (ctx)->last_txt_index)
#define last_txt(vb, did_i) last_txtx (vb, &(vb)->contexts[did_i])
#define last_txt_len(did_i) contexts[did_i].last_txt_len
static inline bool is_last_txt_valid(ContextP ctx) { return ctx->last_txt_index != INVALID_LAST_TXT_INDEX; }
static inline bool is_same_last_txt(VBlockP vb, ContextP ctx, STRp(str)) { return str_len == ctx->last_txt_len && !memcmp (str, last_txtx(vb, ctx), str_len); }

static inline void ctx_init_iterator (Context *ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern WordIndex ctx_evaluate_snip_seg (VBlockP segging_vb, ContextP vctx, STRp (snip), bool *is_new);
extern WordIndex ctx_create_node (VBlockP vb, DidIType did_i, STRp (snip));

#define LASTb250(ctx) ((ctx)->flags.all_the_same ? *FIRSTENT(WordIndex, (ctx)->b250) : *LASTENT(WordIndex, (ctx)->b250))
extern void ctx_append_b250 (VBlockP vb, ContextP vctx, WordIndex node_index);

extern int64_t ctx_decrement_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_increment_count (VBlockP vb, ContextP ctx, WordIndex node_index);

extern WordIndex ctx_get_next_snip (VBlockP vb, ContextP ctx, bool all_the_same, bool is_pair, pSTRp (snip));
extern WordIndex ctx_peek_next_snip (VBlockP vb, ContextP ctx, bool all_the_same, pSTRp (snip));  

extern WordIndex ctx_search_for_word_index (Context *ctx, const char *snip, unsigned snip_len);
extern void ctx_clone (VBlockP vb);
extern CtxNode *ctx_node_vb_do (const Context *ctx, WordIndex node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
extern CtxNode *ctx_node_zf_do (const Context *ctx, int32_t node_index, const char **snip_in_dict, uint32_t *snip_len, const char *func, uint32_t code_line);
#define ctx_node_zf(ctx, node_index, snip_in_dict, snip_len) ctx_node_zf_do(ctx, node_index, snip_in_dict, snip_len, __FUNCTION__, __LINE__)
extern void ctx_merge_in_vb_ctx (VBlockP vb);
extern void ctx_add_compressor_time_to_zf_ctx (VBlockP vb);
extern void ctx_commit_codec_to_zf_ctx (VBlockP vb, ContextP vctx, bool is_lcodec);

extern Context *ctx_get_ctx_if_not_found_by_inline (Context *contexts, DataType dt, DidIType *dict_id_to_did_i_map, DidIType map_did_i, DidIType *num_contexts, DictId dict_id, const char *tag_name, unsigned tag_name_len);

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define ctx_get_ctx(vb,dict_id) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->dict_id_to_did_i_map, &((VBlockP)(vb))->num_contexts, (DictId)(dict_id), 0, 0)
#define ctx_get_ctx_tag(vb,dict_id,tag_name,tag_name_len) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->dict_id_to_did_i_map, &((VBlockP)(vb))->num_contexts, (DictId)(dict_id), (tag_name), (tag_name_len))
static inline Context *ctx_get_ctx_do (Context *contexts, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_contexts, DictId dict_id, const char *tag_name, unsigned tag_name_len)
{
    DidIType did_i = dict_id_to_did_i_map[dict_id.map_key];
    if (did_i != DID_I_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return &contexts[did_i];
    else    
        return ctx_get_ctx_if_not_found_by_inline (contexts, dt, dict_id_to_did_i_map, did_i, num_contexts, dict_id, tag_name, tag_name_len);
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
#define ECTX(dict_id) ctx_get_existing_ctx_do ((VBlockP)(vb), (DictId)(dict_id))

extern ContextP ctx_add_new_zf_ctx_from_txtheader (const char *tag_name, unsigned tag_name_len, DictId dict_id, TranslatorId luft_translator);

extern void ctx_overlay_dictionaries_to_vb (VBlockP vb);
extern void ctx_sort_dictionaries_vb_1(VBlockP vb);

extern void ctx_update_stats (VBlockP vb);
extern void ctx_free_context (Context *ctx, DidIType did_i);
extern void ctx_destroy_context (Context *ctx, DidIType did_i);
extern bool ctx_is_show_dict_id (DictId dict_id);

extern CtxNode *ctx_get_node_by_word_index (ConstContextP ctx, WordIndex word_index);
extern const char *ctx_get_snip_by_word_index (ConstContextP ctx, WordIndex word_index, 
                                               const char **snip, uint32_t *snip_len);

extern const char *ctx_get_vb_snip_ex (ConstContextP vctx, WordIndex vb_node_index, const char **snip, uint32_t *snip_len); 
static inline const char *ctx_get_vb_snip (ConstContextP vctx, WordIndex vb_node_index) { return ctx_get_vb_snip_ex (vctx, vb_node_index, 0, 0); }
 
static inline const char *ctx_get_words_snip(ConstContextP ctx, WordIndex word_index) 
    { return ctx_get_snip_by_word_index (ctx, word_index, 0, 0); }

extern const char *ctx_get_snip_by_zf_node_index (ConstBufferP nodes, ConstBufferP dict, WordIndex node_index, 
                                                  const char **snip, uint32_t *snip_len);
static inline const char *ctx_get_zf_nodes_snip(ConstContextP ctx, WordIndex node_index) 
    { return ctx_get_snip_by_zf_node_index (&ctx->nodes, &ctx->dict, node_index, 0, 0); }

extern WordIndex ctx_get_word_index_by_snip (ConstContextP ctx, const char *snip);

extern void ctx_initialize_predefined_ctxs (Context *contexts /* an array */, DataType dt, DidIType *dict_id_to_did_i_map, DidIType *num_contexts);

extern void ctx_read_all_dictionaries (void);
extern void ctx_compress_dictionaries (void);
extern void ctx_read_all_counts (void);
extern void ctx_compress_counts (void);
extern const char *ctx_get_snip_with_largest_count (DidIType did_i, int64_t *count);
extern void ctx_populate_zf_ctx_from_contigs (Reference ref, DidIType dst_did_i, ConstContigPkgP ctgs);
extern WordIndex ctx_populate_zf_ctx (DidIType dst_did_i, STRp (contig_name), WordIndex ref_index);

extern void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local);

typedef struct { char s[MAX_TAG_LEN+8]; } TagNameEx;
TagNameEx ctx_tag_name_ex (ConstContextP ctx);

// called before seg, to store the point to which we might roll back
static inline void ctx_create_rollback_point (ContextP ctx)
{
    ctx->rback_b250_len       = ctx->b250.len;
    ctx->rback_local_len      = ctx->local.len;
    ctx->rback_txt_len        = ctx->txt_len;
    ctx->rback_num_singletons = ctx->num_singletons;
    ctx->rback_last_value     = ctx->last_value;
    ctx->rback_last_delta     = ctx->last_delta;
    ctx->rback_last_txt_index = ctx->last_txt_index;
    ctx->rback_last_txt_len   = ctx->last_txt_len;
}
extern void ctx_rollback (VBlockP vb, ContextP ctx);

// returns true if dict_id was *previously* segged on this line, and we stored a valid last_value (int or float)
#define ctx_has_value_in_line_(vb, ctx) ((ctx)->last_line_i == (vb)->line_i)
static inline bool ctx_has_value_in_line_do (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    Context *ctx = ECTX (dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_has_value_in_line_(vb, ctx);
}
#define ctx_has_value_in_line(vb, dict_id, p_ctx) ctx_has_value_in_line_do ((VBlockP)(vb), (DictId)(dict_id), (p_ctx))

static inline void ctx_set_last_value (VBlockP vb, ContextP ctx, LastValueType last_value)
{
    ctx->last_value    = last_value;
    ctx->last_line_i   = vb->line_i;
}

// returns true if dict_id was *previously* segged on this line (last_value may be valid or not)
#define ctx_encountered_in_line_(vb, ctx) (((ctx)->last_line_i == (vb)->line_i) || ((ctx)->last_line_i == -(int32_t)(vb)->line_i - 1))
static inline bool ctx_encountered_in_line_do (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    Context *ctx = ECTX (dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_encountered_in_line_(vb, ctx);
}
#define ctx_encountered_in_line(vb, dict_id, p_ctx) ctx_encountered_in_line_do ((VBlockP)(vb), (DictId)(dict_id), (p_ctx))

#define ctx_set_encountered_in_line(ctx)  /* set encountered if not already ctx_set_last_value */  \
    if ((ctx)->last_line_i != vb->line_i) \
        (ctx)->last_line_i   = -(int32_t)vb->line_i - 1; 

extern void ctx_foreach_buffer(ContextP ctx, bool set_name, void (*func)(BufferP buf, const char *func, unsigned line));

#endif