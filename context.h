// ------------------------------------------------------------------
//   context.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "sections.h"
#include "bits.h"
#include "mutex.h"
#include "context_struct.h"
#include "vblock.h"
#include "segconf.h"
#include "strings.h"

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
// We can use non-printable characters that - i.e. other than ASCII 32-126, \t (\x9) \n (\xA) \r (\xD)
#define SNIP_SEP                  '\x0'   // Seperator between snips - both for dict and local 
#define SNIP_LOOKUP               '\x1'   // Lookup from local (optionally followed by a snip - interpreted differently by local type, see reconstruct_one_snip)
#define SNIP_OTHER_LOOKUP         '\x2'   // Lookup from local of other dict_id (possibly with length for sequence storage)
#define v13_SNIP_MATE_LOOKUP      '\x3'   // up to v13: (replaced by FASTQ_SPECIAL_mate_lookup) Lookup from paired file (when using --pair)  
#define SNIP_CONTAINER            '\x4'   // Appears as first character in the SNIP, followed by a specification of a container field
#define SNIP_SELF_DELTA           '\x5'   // The value is a uint32_t which is a result of the last value + the positive or negative textual int32_t value following this character
#define SNIP_OTHER_DELTA          '\x6'   // The value is a uint32_t which is a result of the last value of another field + the delta value. following this char, {DictId dict_id, int32_t delta, bool update_other} in base64)
#define v13_SNIP_FASTQ_PAIR2_GPOS '\x7'   // up to v13: (replaced by FASTQ_SPECIAL_PAIR2_GPOS) The value is a uint32_t which is a result of the equivalent value in the paired file + the delta value (when using --pair)
#define SNIP_SPECIAL              '\x8'   // Special algorithm followed by ID of the algorithm 
#define SNIP_REDIRECTION          '\xB'   // Get the data from another dict_id (can be in b250, local...)
#define SNIP_DONT_STORE           '\xC'   // Reconstruct the following value, but don't store it in last_value (overriding flags.store)
#define SNIP_COPY                 '\xE'   // Copy the last_txt of dict_id (same or other)
#define SNIP_DUAL                 '\xF'   // A snip containing two snips separated by a SNIP_DUAL - for Primary and Luft reconstruction respectively
#define SNIP_LOOKBACK             '\x10'  // Copy an earlier snip in the same context. Snip is dict_id from which to take the lookback offset, and an optional delta to be applied to the retrieved numeric value. note: line number of the previous snip is variable, but its offset back is fixed (introduced 12.0.41)
#define SNIP_COPY_BUDDY           '\x11'  // Copy a snip on an earlier "buddy" line in the same or another context (note: offset back to the previous snip is variable, but its line number is fixed) (introduced 12.0.41)
#define SNIP_XOR_DIFF             '\x12'  // XOR a string vs. previous string (introduced 13.0.5)    
#define SNIP_RESERVED             '\x13'  // A value guaranteed not to exist in dictionary data. Used internally by ctx_shorten_unused_dict_words. (13.0.7)
#define NUM_SNIP_CODES            20

#define SNIP_CODES { "SNIP_SEP", "SNIP_LOOKUP", "SNIP_OTHER_LOOKUP", "SNIP_MATE_LOOKUP",\
                     "SNIP_CONTAINER", "SNIP_SELF_DELTA", "SNIP_OTHER_DELTA", "v13_SNIP_FASTQ_PAIR2_GPOS", \
                     "SNIP_SPECIAL", "<TAB>", "<NL>", "SNIP_REDIRECTION",\
                     "SNIP_DONT_STORE", "<LF>", "SNIP_COPY", "SNIP_DUAL", "SNIP_LOOKBACK",\
                     "SNIP_COPY_BUDDY", "SNIP_XOR_DIFF", "SNIP_RESERVED" }

// Format on data in Context.b250: Each entry is either a single-byte special-code value 0xFA-0xFF, OR a 1, 2 or 4 big-endian integer.
// The number of bytes is determined by Context.b250_size transmitted via SectionHeaderCtx.b250_size, its selection is done separately for each VB.
// In FASTQ, for b250-pairing, b250 of the pair is stored in Context.pair and its size in Context.pair_b250_size
#define BASE250_EMPTY_SF   0xFA // 250 empty string
#define BASE250_MISSING_SF 0xFB // 251 container item missing, remove preceding separator
#define BASE250_ONE_UP     0xFC // 252 value is one higher than previous value. 
#define BASE250_MOST_FREQ0 0xFD // 253 this translates to 0,1,2 representing the most frequent values (according to vb_i=1 sorting).
#define BASE250_MOST_FREQ1 0xFE // 254
#define BASE250_MOST_FREQ2 0xFF // 255

#define B250_MAX_WI_1BYTE  (0xFA         - 1) // maximal value in Context.b250 of non-special-code entries, for b250_size to be set to 1 byte
#define B250_MAX_WI_2BYTES ((0xFA << 8)  - 1) // same, for 2 bytes
#define B250_MAX_WI_3BYTES ((0xFA << 16) - 1) // same, for 3 bytes

// ZIP only
typedef struct CtxNode {
    CharIndex char_index;     // character index into dictionary array
    uint32_t snip_len;        // not including \0 terminator present in dictionary array
    union {
        WordIndex word_index; // zctx, vctx->ol_nodes and vctx->nodes after ctx_merge_in_one_vctx 
        WordIndex node_index; // vctx->nodes before ctx_merge_in_one_vctx
    };
} CtxNode;

typedef ZWord CtxWord;

// for signed numbers, we store them in our "interlaced" format rather than standard ISO format 
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// compressing better in an array that contains both positive and negative
// Note: the full range of the signed type is supported, eg -128 to 127 for int8_t
#define SAFE_NEGATE(signedtype,n) ((u##signedtype)(-((int64_t)n))) // careful negation to avoid overflow eg -(-128)==0 in int8_t
#define INTERLACE(signedtype,num) ({ signedtype n=(signedtype)(num); (n < 0) ? ((SAFE_NEGATE(signedtype,n) << 1) - 1) : (((u##signedtype)n) << 1); })
#define DEINTERLACE(signedtype,unum) (((unum) & 1) ? -(signedtype)(((unum)>>1)+1) : (signedtype)((unum)>>1))

static inline bool NEXTLOCALBIT(ContextP ctx)      { bool ret = bits_get ((BitsP)&ctx->local, ctx->next_local); ctx->next_local++; return ret; }
static inline uint8_t NEXTLOCAL2BITS(ContextP ctx) { uint8_t ret = bits_get ((BitsP)&ctx->local, ctx->next_local) | (bits_get ((BitsP)&ctx->local, ctx->next_local+1) << 1); ctx->next_local += 2; return ret; }

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75  

#define ctx_node_vb(ctx, node_index, snip_in_dict, snip_len) ctx_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCLINE)
#define node_word_index(vb,did_i,index) ((index)!=WORD_INDEX_NONE ? ctx_node_vb (&(vb)->contexts[did_i], (index), 0,0)->word_index : WORD_INDEX_NONE)

#define CTX(did_i)   ({ Did my_did_i = (did_i); /* could be an expression */\
                        ASSERT (my_did_i < MAX_DICTS, "CTX(): did_i=%u out of range", my_did_i); /* will be optimized-out for constant did_i */ \
                        (&vb->contexts[my_did_i]); })

#define LOADED_CTX(did_i) ({ ContextP ctx=CTX(did_i); ASSISLOADED(ctx); ctx; })

#define ZCTX(did_i)  ({ Did my_did_i = (did_i);\
                        ASSERT (my_did_i < MAX_DICTS, "ZCTX(): did_i=%u out of range", my_did_i);  \
                        &z_file->contexts[my_did_i]; })

#define last_int(did_i)     contexts[did_i].last_value.i
#define last_index(did_i)   contexts[did_i].last_value.i
#define last_float(did_i)   contexts[did_i].last_value.f
#define last_delta(did_i)   contexts[did_i].last_delta
#define last_txtx(vb, ctx)  Btxt ((ctx)->last_txt.index)
static inline char *last_txt (VBlockP vb, Did did_i) { return last_txtx (vb, CTX(did_i)); }
#define last_txt_len(did_i) contexts[did_i].last_txt.len

#define history64(did_i, line_i) (*B(int64_t, CTX(did_i)->history, (line_i)))

static inline bool is_last_txt_valid(ContextP ctx) { return ctx->last_txt.index != INVALID_LAST_TXT_INDEX; }
static inline bool is_same_last_txt(VBlockP vb, ContextP ctx, STRp(str)) { return str_len == ctx->last_txt.len && !memcmp (str, last_txtx(vb, ctx), str_len); }

static inline void ctx_init_iterator (ContextP ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern WordIndex ctx_create_node_do (VBlockP segging_vb, ContextP vctx, STRp (snip), bool *is_new);
extern WordIndex ctx_create_node (VBlockP vb, Did did_i, STRp (snip));

#define LASTb250(ctx) ((ctx)->flags.all_the_same ? *B1ST(WordIndex, (ctx)->b250) : *BLST(WordIndex, (ctx)->b250))
extern void ctx_append_b250 (VBlockP vb, ContextP vctx, WordIndex node_index);

extern uint32_t ctx_get_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_decrement_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_increment_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_protect_from_removal (VBlockP vb, ContextP ctx, WordIndex node_index);

extern WordIndex ctx_decode_b250 (bytes *b, bool advance, B250Size b250_size, rom ctx_name);
extern WordIndex ctx_get_next_snip (VBlockP vb, ContextP ctx, bool is_pair, pSTRp (snip));
extern WordIndex ctx_peek_next_snip (VBlockP vb, ContextP ctx, pSTRp (snip));  

extern WordIndex ctx_search_for_word_index (ContextP ctx, STRp(snip));
extern void ctx_clone (VBlockP vb);
extern CtxNode *ctx_node_vb_do (ConstContextP ctx, WordIndex node_index, rom *snip_in_dict, uint32_t *snip_len, FUNCLINE);
extern CtxNode *ctx_node_zf_do (ConstContextP ctx, int32_t node_index, rom *snip_in_dict, uint32_t *snip_len, FUNCLINE);
#define ctx_node_zf(ctx, node_index, snip_in_dict, snip_len) ctx_node_zf_do(ctx, node_index, snip_in_dict, snip_len, __FUNCLINE)
extern void ctx_merge_in_vb_ctx (VBlockP vb);
extern void ctx_substract_txt_len (VBlockP vb, ContextP vctx);
extern void ctx_add_compressor_time_to_zf_ctx (VBlockP vb);
extern void ctx_commit_codec_to_zf_ctx (VBlockP vb, ContextP vctx, bool is_lcodec);
extern void ctx_convert_generated_b250_to_mate_lookup (VBlockP vb, ContextP vctx);

extern ContextP ctx_get_unmapped_ctx (ContextP contexts, DataType dt, Did *dict_id_to_did_i_map, Did *num_contexts, DictId dict_id, STRp(tag_name));

// returns did_i of dict_id if it is found in the map, or DID_NONE if not
static inline Did get_matching_did_i_from_map (const Context *contexts, const Did *map, DictId dict_id)
{
    Did did_i = map[dict_id.map_key];
    if (did_i != DID_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return did_i;

    did_i = map[ALT_KEY(dict_id)];
    if (did_i != DID_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return did_i;

    return DID_NONE;
}

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define ctx_get_ctx(vb,dict_id) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->dict_id_to_did_i_map, &((VBlockP)(vb))->num_contexts, (dict_id), 0, 0)
#define ctx_get_ctx_tag(vb,dict_id,tag_name,tag_name_len) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->dict_id_to_did_i_map, &((VBlockP)(vb))->num_contexts, (dict_id), (tag_name), (tag_name_len))
static inline ContextP ctx_get_ctx_do (Context *contexts, DataType dt, Did *dict_id_to_did_i_map, Did *num_contexts, DictId dict_id, STRp(tag_name))
{
    Did did_i = get_matching_did_i_from_map (contexts, dict_id_to_did_i_map, dict_id);
    if (did_i != DID_NONE) 
        return &contexts[did_i];
    else    
        return ctx_get_unmapped_ctx (contexts, dt, dict_id_to_did_i_map, num_contexts, dict_id, STRa(tag_name));
}

extern Did ctx_get_unmapped_existing_did_i (const Context *contexts, const ContextIndex *ctx_index, Did num_contexts, DictId dict_id);

static inline Did ctx_get_existing_did_i_do (DictId dict_id, const Context *contexts, Did *dict_id_to_did_i_map, const ContextIndex *ctx_index, Did num_contexts)
{
    Did did_i = get_matching_did_i_from_map (contexts, dict_id_to_did_i_map, dict_id);
    if (did_i != DID_NONE)
        return did_i;
    else
        return ctx_get_unmapped_existing_did_i (contexts, ctx_index, num_contexts, dict_id);
}    
#define ctx_get_existing_did_i(vb,dict_id) ctx_get_existing_did_i_do (dict_id, vb->contexts, vb->dict_id_to_did_i_map, (vb->has_ctx_index ? vb->ctx_index : NULL), vb->num_contexts)

static inline ContextP ctx_get_existing_ctx_do (VBlockP vb, DictId dict_id)  // returns NULL if context doesn't exist
{
    Did did_i = ctx_get_existing_did_i(vb, dict_id); 
    return (did_i == DID_NONE) ? NULL : CTX(did_i); 
}
#define ECTX(dict_id) ctx_get_existing_ctx_do ((VBlockP)(vb), (dict_id))

extern ContextP ctx_get_zctx_from_vctx (ConstContextP vctx);

extern ContextP ctx_add_new_zf_ctx_from_txtheader (STRp(tag_name), DictId dict_id, TranslatorId luft_translator);

extern void ctx_overlay_dictionaries_to_vb (VBlockP vb);
extern void ctx_sort_dictionaries_vb_1(VBlockP vb);

extern void ctx_update_stats (VBlockP vb);
extern void ctx_free_context (ContextP ctx, Did did_i);
extern void ctx_destroy_context (ContextP ctx, Did did_i);

extern CtxNode *ctx_get_node_by_word_index (ConstContextP ctx, WordIndex word_index);
extern rom ctx_get_snip_by_word_index_do (ConstContextP ctx, WordIndex word_index, pSTRp(snip), FUNCLINE);
#define ctx_get_snip_by_word_index(ctx,word_index,snip) ctx_get_snip_by_word_index_do ((ctx), (word_index), &snip, &snip##_len, __FUNCLINE)
#define ctx_get_snip_by_word_index0(ctx,word_index) ctx_get_snip_by_word_index_do ((ctx), (word_index), 0,0, __FUNCLINE)
                                                
extern rom ctx_get_vb_snip_ex (ConstContextP vctx, WordIndex vb_node_index, pSTRp(snip)); 
static inline rom ctx_get_vb_snip (ConstContextP vctx, WordIndex vb_node_index) { return ctx_get_vb_snip_ex (vctx, vb_node_index, 0, 0); }
 
static inline rom ctx_get_words_snip(ConstContextP ctx, WordIndex word_index)  // PIZ
    { return ctx_get_snip_by_word_index0 (ctx, word_index); }

extern WordIndex ctx_get_word_index_by_snip (VBlockP vb, ContextP ctx, STRp(snip));

extern rom ctx_snip_from_zf_nodes (ConstContextP zctx, WordIndex node_index, pSTRp(snip));

extern WordIndex ctx_get_ol_node_index_by_snip (VBlockP vb, ContextP ctx, STRp(snip)); 

extern void ctx_initialize_predefined_ctxs (ContextP contexts /* an array */, DataType dt, Did *dict_id_to_did_i_map, Did *num_contexts);

extern void ctx_shorten_unused_dict_words (Did did_i);
extern void ctx_read_all_counts (void);
extern void ctx_compress_counts (void);
extern rom ctx_get_snip_with_largest_count (Did did_i, int64_t *count);
extern void ctx_populate_zf_ctx_from_contigs (Reference ref, Did dst_did_i, ConstContigPkgP ctgs);
extern WordIndex ctx_populate_zf_ctx (Did dst_did_i, STRp (contig_name), WordIndex ref_index);

extern void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local);

typedef struct { char s[MAX_TAG_LEN+8]; } TagNameEx;
TagNameEx ctx_tag_name_ex (ConstContextP ctx);

// Seg: called before seg, to store the point to which we might roll back
static inline bool ctx_set_rollback (VBlockP vb, ContextP ctx, bool override_id)
{
    if (ctx->rback_id == vb->rback_id && !override_id) return false; // ctx already in rollback point 
    
    ctx->rback_b250_len        = ctx->b250.len32;
    ctx->rback_local_len       = ctx->local.len32;
    ctx->rback_nodes_len       = ctx->nodes.len32;
    ctx->rback_txt_len         = ctx->txt_len;
    ctx->rback_local_num_words = ctx->local_num_words;
    ctx->rback_b250_count      = ctx->b250.count;
    ctx->rback_last_value      = ctx->last_value;
    ctx->rback_last_delta      = ctx->last_delta;
    ctx->rback_last_txt        = ctx->last_txt;
    ctx->rback_id              = vb->rback_id;
    return true;
}

// Needed only if we intend to call ctx_set_rollback again without advancing vb->rback_id
static inline void ctx_unset_rollback (ContextP ctx) { ctx->rback_id = 1; }

extern void ctx_rollback (VBlockP vb, ContextP ctx, bool override_id);

static inline bool ctx_can_have_singletons (ContextP ctx) 
    { return (ctx->ltype == LT_TEXT) && !ctx->no_stons && 
              (ctx->flags.store != STORE_INDEX) && !ctx->flags.all_the_same; } 

// returns true if dict_id was *previously* segged on this line, and we stored a valid last_value (int or float)
#define ctx_has_value_in_line_(vb, ctx) ((ctx)->last_line_i == (vb)->line_i)
#define ctx_has_value_in_prev_line_(vb, ctx) ((ctx)->last_line_i+1 == (vb)->line_i)

static inline bool ctx_has_value_in_line_do (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    ContextP ctx = ECTX (dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_has_value_in_line_(vb, ctx);
}
#define ctx_has_value_in_line(vb, dict_id, p_ctx) ctx_has_value_in_line_do ((VBlockP)(vb), (dict_id), (p_ctx))

static inline void ctx_set_last_value (VBlockP vb, ContextP ctx, ValueType last_value)
{
    ctx->last_value    = last_value;
    ctx->last_line_i   = vb->line_i;
    ctx->last_sample_i = vb->sample_i; // used for VCF/FORMAT. otherwise meaningless but harmless.
}

// returns true if value is set
static inline bool ctx_set_last_value_from_str (VBlockP vb, ContextP ctx, STRp(str))
{
    if ((ctx->flags.store == STORE_INT   && str_get_int   (STRa(str), &ctx->last_value.i)) ||
        (ctx->flags.store == STORE_FLOAT && str_get_float (STRa(str), &ctx->last_value.f, 0, 0))) {

        ctx->last_line_i   = vb->line_i;
        ctx->last_sample_i = vb->sample_i; // used for VCF/FORMAT. otherwise meaningless but harmless.
        return true;
    }
    else
        return false;
}

// set encountered if not already ctx_set_last_value (encounted = seen, but without setting last_value)
static inline void ctx_set_encountered (VBlockP vb, ContextP ctx)
{
    if (ctx->last_line_i != vb->line_i || ctx->last_sample_i != vb->sample_i) // not already ctx_set_last_value in this line/sample
        ctx->last_line_i = -(int64_t)vb->line_i - 1; 

    ctx->last_sample_i = vb->sample_i; 
}

// after calling this, all these are false: ctx_encountered*, ctx_has_value*
static inline void ctx_unset_encountered (VBlockP vb, ContextP ctx)
{
    ctx->last_line_i = LAST_LINE_I_INIT;
}

// returns true if dict_id was *previously* segged on this line (last_value may be valid or not)
static inline bool ctx_encountered_in_line (VBlockP vb, Did did_i) 
{ 
    ContextP ctx = CTX(did_i); 
    return ((ctx->last_line_i == vb->line_i) || (ctx->last_line_i == -(int64_t)vb->line_i - 1)); 
}

static inline bool ctx_encountered_in_prev_line (VBlockP vb, Did did_i) 
{ 
    ContextP ctx = CTX(did_i); 
    return ((ctx->last_line_i == vb->line_i-1) || (ctx->last_line_i == -(int64_t)(vb->line_i-1) - 1)); 
}

static inline bool ctx_encountered_in_line_by_dict_id (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    ContextP ctx = ECTX (dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_encountered_in_line(vb, ctx->did_i);
}

// note: these macros are good for contexts in VCF/FORMAT or not. Use the *_in_line version in case
// we are sure its not VCF/FORMAT, they are a bit more efficient.
static inline bool ctx_encountered (VBlockP vb, Did did_i) 
{ 
    ContextP ctx = CTX(did_i);
    return ctx_encountered_in_line(vb, ctx->did_i) && ctx->last_sample_i == vb->sample_i; 
}

static inline bool ctx_encountered_by_dict_id (VBlockP vb, DictId dict_id, ContextP *p_ctx) 
    { return ctx_encountered_in_line_by_dict_id (vb, dict_id, p_ctx) && (*p_ctx)->last_sample_i == vb->sample_i; }
    
static inline bool ctx_has_value (VBlockP vb, Did did_i) 
{   
    ContextP ctx = CTX(did_i);
    return ctx_has_value_in_line_(vb, ctx) && ctx->last_sample_i == vb->sample_i; 
}

static inline bool ctx_has_value_by_dict_id (VBlockP vb, DictId dict_id, ContextP *p_ctx) 
    { return ctx_has_value_in_line_do (vb, dict_id, p_ctx) && (*p_ctx)->last_sample_i == vb->sample_i; }

extern uint64_t ctx_get_ctx_group_z_len (VBlockP vb, Did group_did_i);

typedef enum { KR_KEEP, KR_REMOVE } CtxKeepRemove;
extern void ctx_declare_winning_group (Did winning_group_did_i, Did losing_group_did_i, Did new_st_did_i);

extern void ctx_set_store (VBlockP vb, StoreType store_type, unsigned num_ctxs, ...);
extern void ctx_set_no_stons (VBlockP vb, unsigned num_ctxs, ...);
extern void ctx_set_delta_peek (VBlockP vb, unsigned num_ctxs, ...);
extern void ctx_set_store_per_line (VBlockP vb, unsigned num_ctxs, ...);
extern void ctx_set_ltype (VBlockP vb, LocalType ltype, unsigned num_ctxs, ...);
extern void ctx_consolidate_stats (VBlockP vb, Did parent, unsigned num_deps, ...);
extern void ctx_consolidate_statsN(VBlockP vb, Did parent, Did first_dep, unsigned num_deps);
extern void ctx_consolidate_stats_(VBlockP vb, Did parent, unsigned num_deps, ContextP *dep_ctxs);

extern rom dyn_type_name (DynType dyn_type);

extern void ctx_foreach_buffer(ContextP ctx, bool set_name, void (*func)(BufferP buf, FUNCLINE));

extern SORTER (sort_by_dict_id);

#define for_zctx for (ContextP zctx=ZCTX(0); zctx < ZCTX(z_file->num_contexts); zctx++) 
#define for_vctx for (ContextP vctx=CTX(0);  vctx < CTX(vb->num_contexts);      vctx++) 
#define for_ctx  for (ContextP ctx=CTX(0);   ctx  < CTX(vb->num_contexts);      ctx++) 
