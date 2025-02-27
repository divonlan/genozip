// ------------------------------------------------------------------
//   context.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "sections.h"
#include "bits.h"
#include "mutex.h"
#include "context_struct.h"
#include "vblock.h"
#include "segconf.h"
#include "strings.h"
#include "dict_io.h"
#include "file.h"
#include "sorter.h"

#define MAX_WORDS_IN_CTX (1<<29) // limit on nodes.len, word_list.len - limited by VARL_MAX_4B

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
#define SNIP_REDIRECTION          '\xB'   // Get the data from another context (can be in b250, local...)
#define SNIP_DONT_STORE           '\xC'   // Reconstruct the following value, but don't store it in last_value (overriding flags.store)
#define SNIP_COPY                 '\xE'   // Copy the last_txt of (optional) dict_id (self if no dict_id)
#define dvcf_SNIP_DUAL            '\xF'   // up to 15.0.41: A snip containing two snips separated by a SNIP_DUAL - for Primary and Luft reconstruction respectively
#define SNIP_LOOKBACK             '\x10'  // Copy an earlier snip in the same context. Snip is dict_id from which to take the lookback offset, and an optional delta to be applied to the retrieved numeric value. note: line number of the previous snip is variable, but its offset back is fixed (introduced 12.0.41)
#define v13_SNIP_COPY_BUDDY       '\x11'  // up to v13: Copy a snip on an earlier "buddy" line in the same or another context (note: offset back to the previous snip is variable, but its line number is fixed) (introduced 12.0.41)
#define SNIP_DIFF                 '\x12'  // XOR a string vs. previous string (introduced 13.0.5)    
#define SNIP_RESERVED             '\x13'  // A value guaranteed not to exist in dictionary data. Used internally by ctx_shorten_unused_dict_words. (13.0.7)
#define SNIP_NUMERIC              '\x14'  // Lookup for local, and format output (introduced v14)
#define NUM_SNIP_CODES            21

#define SNIP_CODES { "SNIP_SEP", "SNIP_LOOKUP", "SNIP_OTHER_LOOKUP", "SNIP_MATE_LOOKUP",\
                     "SNIP_CONTAINER", "SNIP_SELF_DELTA", "SNIP_OTHER_DELTA", "v13_SNIP_FASTQ_PAIR2_GPOS", \
                     "SNIP_SPECIAL", "<TAB>", "<NL>", "SNIP_REDIRECTION",\
                     "SNIP_DONT_STORE", "<LF>", "SNIP_COPY", "SNIP_DUAL", "SNIP_LOOKBACK",\
                     "SNIP_COPY_BUDDY", "SNIP_DIFF", "SNIP_RESERVED", "SNIP_NUMERIC" }

#define MAX_SNIP_DICTS 3 // maximum dicts for PLUS snip. Can be increased if needed, but never decreased.

#define decl_ctx(did_i) ContextP ctx = CTX(did_i)
#define decl_const_ctx(did_i) ConstContextP ctx = CTX(did_i)
#define decl_zctx(did_i) ContextP zctx = ZCTX(did_i)

// PIZ

// type of elements of ctx->word_list 
typedef struct { // 8 bytes 
    uint64_t char_index : 39; // up to CTX_MAX_DICT_LEN (was 40 up to 15.0.37)
    uint64_t snip_len   : 24; // up to CTX_MAX_SNIP_LEN
} CtxWord, *CtxWordP; 

// ZIP 
#pragma pack(4) // note: MUST be at least 4 or 8 for prevs_next assignment in hash_global_add_node() to be thread-safe
#define INITIAL_NUM_NODES 1024

#define CTX_MAX_DICT_LEN (512 GB - 1ULL) // (39 bit) maximum length of a context.dict (was: 1TB (40bit) from v14 to 15.0.37) 
#define CTX_MAX_SNIP_LEN (16 MB - 1ULL)  // (24 bit) maximum length of any snip in a context.dict (excluding its \0 separator) (v14)

// type of elements of zctx->nodes, vctx->ol_nodes and vctx->nodes entries before conversion to WordIndex. 
typedef struct {              // 12 bytes
    uint64_t char_index : 39; // up to CTX_MAX_DICT_LEN
    uint64_t snip_len   : 24; // up to CTX_MAX_SNIP_LEN
    uint64_t canceled   : 1;  // set if node was canceled (only happens in vctx->nodes)
    uint32_t next;            // linked list - index in buffer of next node in linked list or NO_NEXT (note: in vctx->nodes - this is NOT node_index, since we don't add ol_nodes.len)
} CtxNode, *CtxNodeP;

typedef struct {              // 8 bytes
    uint32_t digest;          // crc32 digest of the singleton
    uint32_t next;            // linked list - index of next entry on list in zctx->ston_ents or NO_NEXT
} SingletonEnt, *SingletonEntP;

#pragma pack()

// Interlaced integers are used for storing integers that might be positive or negative, while keeping
// as many higher bits zero as possible.
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// Note: the full range of the signed type is supported, eg -128 to 127 for int8_t
#define SAFE_NEGATE(signedtype,n) ((u##signedtype)(-((int64_t)n))) // careful negation to avoid overflow eg -(-128)==0 in int8_t
#define INTERLACE(signedtype,num) ({ signedtype n=(signedtype)(num); (n < 0) ? ((SAFE_NEGATE(signedtype,n) << 1) - 1) : (((u##signedtype)n) << 1); })
#define DEINTERLACE(signedtype,unum) (((unum) & 1) ? -(signedtype)(((unum)>>1)+1) : (signedtype)((unum)>>1))

// factor in which we grow buffers in CTX upon realloc
#define CTX_GROWTH 1.75

#define ctx_node_vb(ctx, node_index, snip_in_dict, snip_len) ctx_node_vb_do(ctx, node_index, snip_in_dict, snip_len, __FUNCLINE)

#define node_index_to_word_index(vb, vctx, vb_node_index) /* use with vctx after conversion */  \
    (((vb_node_index) >= (int32_t)(vctx)->ol_nodes.len32)                                       \
        ? *B(WordIndex, (vctx)->nodes, (vb_node_index) - (vctx)->ol_nodes.len32)                \
        : (vb_node_index))

#define CTX(did_i)   ({ Did my_did_i = (did_i); /* evaluate did_i only once */\
                        ASSERT (my_did_i < MAX_DICTS, "CTX(): did_i=%u out of range", my_did_i); /* optimized out for constant did_i */ \
                        (&vb->contexts[my_did_i]); })

#define LOADED_CTX(did_i) ({ ContextP ctx=CTX(did_i); ASSISLOADED(ctx); ctx; })

#define ZCTX(did_i)  ({ Did my_did_i = (did_i);\
                        ASSERT (my_did_i < MAX_DICTS, "ZCTX(): did_i=%u out of range", my_did_i); /* optimized out for constant did_i */ \
                        &z_file->contexts[my_did_i]; })

#define last_int(did_i)     contexts[did_i].last_value.i
#define last_index(did_i)   contexts[did_i].last_value.i
#define last_float(did_i)   contexts[did_i].last_value.f
#define last_delta(did_i)   contexts[did_i].last_delta
#define last_txtx(vb, ctx)  Btxt ((ctx)->last_txt.index)
static inline char *last_txt (VBlockP vb, Did did_i) { return last_txtx (vb, CTX(did_i)); }
#define last_txt_len(did_i) contexts[did_i].last_txt.len

#define set_last_txtC(ctx, value, value_len) (ctx)->last_txt = (TxtWord){ .index = BNUMtxt(value), .len = (value_len) }
#define set_last_txt_(did_i, value, value_len) set_last_txtC(CTX(did_i), (value), (value_len))
#define set_last_txt(did_i, value) set_last_txt_((did_i), value, value##_len)

#define history64(did_i, line_i) (*B(int64_t, CTX(did_i)->history, (line_i)))

static inline bool is_last_txt_valid(ContextP ctx) { return ctx->last_txt.index != INVALID_LAST_TXT_INDEX; }
static inline bool is_same_last_txt(VBlockP vb, ContextP ctx, STRp(str)) { return str_issame_(STRa(str), STRlst_(ctx)); }

static inline void ctx_init_iterator (ContextP ctx) { ctx->iterator.next_b250 = NULL ; ctx->iterator.prev_word_index = -1; ctx->next_local = 0; }

extern WordIndex ctx_create_node_do (VBlockP segging_vb, ContextP vctx, STRp (snip), bool *is_new);
extern WordIndex ctx_create_node_is_new (VBlockP vb, Did did_i, STRp (snip), bool *is_new);
static inline WordIndex ctx_create_node (VBlockP vb, Did did_i, STRp (snip)) { return ctx_create_node_is_new (vb, did_i, STRa(snip), NULL); }

extern uint32_t ctx_get_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_decrement_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_increment_count (VBlockP vb, ContextP ctx, WordIndex node_index);
extern void ctx_protect_from_removal (VBlockP vb, ContextP ctx, WordIndex node_index);

extern WordIndex ctx_get_next_snip (VBlockP vb, ContextP ctx, bool is_pair, pSTRp (snip));
extern uint32_t ctx_get_next_snip_from_local (VBlockP vb, ContextP ctx, pSTRp (snip));
extern WordIndex ctx_peek_next_snip (VBlockP vb, ContextP ctx, pSTRp (snip));  

extern WordIndex ctx_search_for_word_index (ContextP ctx, STRp(snip));
extern void ctx_clone (VBlockP vb);
extern CtxNode ctx_node_vb_do (ConstContextP ctx, WordIndex node_index, rom *snip_in_dict, uint32_t *snip_len, FUNCLINE);
extern void ctx_merge_in_vb_ctx (VBlockP vb);
extern void ctx_update_zctx_txt_len (VBlockP vb, ContextP vctx, int64_t increment);
extern void ctx_commit_codec_to_zf_ctx (VBlockP vb, ContextP vctx, bool is_lcodec, bool is_lcodec_inherited);
extern void ctx_reset_codec_commits (void);
extern void ctx_segconf_set_hard_coded_lcodec (Did did_i, Codec codec);
extern void ctx_get_z_codecs (ContextP zctx, Codec *lcodec, Codec *bcodec, uint8_t *lcodec_count, uint8_t *bcodec_count, bool *lcodec_hard_coded);

extern ContextP ctx_get_unmapped_ctx (ContextArray contexts, DataType dt, DictIdtoDidMap d2d_map, Did *num_contexts, DictId dict_id, STRp(tag_name));

// returns did_i of dict_id if it is found in the map, or DID_NONE if not
static inline Did get_matching_did_i_from_map (const ContextArray contexts, const Did *map, DictId dict_id)
{
    Did did_i = map[dict_id.map_key[0]];
    if (did_i != DID_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return did_i;

    did_i = map[ALT_KEY(dict_id)];
    if (did_i != DID_NONE && contexts[did_i].dict_id.num == dict_id.num) 
        return did_i;

    return DID_NONE;
}

// inline function for quick operation typically called several billion times in a typical file and > 99.9% can be served by the inline
#define ctx_get_ctx(vb,dict_id) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->d2d_map, &((VBlockP)(vb))->num_contexts, (DictId)(dict_id), 0, 0)
#define ctx_get_ctx_tag(vb,dict_id,tag_name,tag_name_len) ctx_get_ctx_do (((VBlockP)(vb))->contexts, ((VBlockP)(vb))->data_type, ((VBlockP)(vb))->d2d_map, &((VBlockP)(vb))->num_contexts, (DictId)(dict_id), (tag_name), (tag_name_len))
static inline ContextP ctx_get_ctx_do (ContextArray contexts, DataType dt, DictIdtoDidMap d2d_map, Did *num_contexts, DictId dict_id, STRp(tag_name))
{
    Did did_i = get_matching_did_i_from_map (contexts, d2d_map, dict_id);
    if (did_i != DID_NONE) 
        return &contexts[did_i];
    else    
        return ctx_get_unmapped_ctx (contexts, dt, d2d_map, num_contexts, dict_id, STRa(tag_name));
}

extern Did ctx_get_unmapped_existing_did_i (const ContextArray contexts, ConstBufferP ctx_index, Did num_contexts, DictId dict_id);

static inline Did ctx_get_existing_did_i_do (DictId dict_id, const ContextArray contexts, DictIdtoDidMap d2d_map, ConstBufferP ctx_index, Did num_contexts)
{
    Did did_i = get_matching_did_i_from_map (contexts, d2d_map, dict_id);
    if (did_i != DID_NONE)
        return did_i;
    else
        return ctx_get_unmapped_existing_did_i (contexts, ctx_index, num_contexts, dict_id);
}    
#define ctx_get_existing_did_i(vb,dict_id) ctx_get_existing_did_i_do ((dict_id), (vb)->contexts, (vb)->d2d_map, ((vb)->ctx_index.len ? &(vb)->ctx_index : NULL), (vb)->num_contexts)
#define zctx_get_existing_did_i(dict_id) ctx_get_existing_did_i_do ((dict_id), z_file->contexts, z_file->d2d_map, NULL, z_file->num_contexts)

static inline ContextP ctx_get_existing_ctx_do (VBlockP vb, DictId dict_id, FUNCLINE)  // returns NULL if context doesn't exist
{
    Did did_i = ctx_get_existing_did_i (vb, dict_id); 
    ASSERT (IN_RANGE((int16_t)did_i, -1, vb->num_contexts), "%s:%u: ECTX: did_i=%d ∉ [-1,%d). dict_id=%s", 
            func, code_line, did_i, vb->num_contexts, dis_dict_id (dict_id).s); 

    return (did_i == DID_NONE) ? NULL : &vb->contexts[did_i]; 
}
#define ECTX(dict_id) ctx_get_existing_ctx_do ((VBlockP)(vb), (DictId)(dict_id), __FUNCTION__, __LINE__)

extern ContextP ctx_get_existing_zctx (DictId dict_id);

extern ContextP ctx_get_zctx_from_vctx (ConstContextP vctx, bool create_if_missing, bool follow_alias);

extern ContextP ctx_add_new_zf_ctx_at_init (STRp(tag_name), DictId dict_id);

extern void ctx_overlay_dictionaries_to_vb (VBlockP vb);

extern void ctx_update_stats (VBlockP vb);

// PIZ - get snip
extern rom ctx_get_snip_by_word_index_do (ConstContextP ctx, WordIndex word_index, pSTRp(snip), FUNCLINE);
#define ctx_get_snip_by_word_index(ctx,word_index,snip) ctx_get_snip_by_word_index_do ((ctx), (word_index), &snip, &snip##_len, __FUNCLINE)
#define ctx_get_snip_by_word_index0(ctx,word_index) ctx_get_snip_by_word_index_do ((ctx), (word_index), 0,0, __FUNCLINE)

static inline rom ctx_get_words_snip(ConstContextP ctx, WordIndex word_index)  // PIZ
    { return ctx_get_snip_by_word_index0 (ctx, word_index); }

// ZIP - get snip
extern rom ctx_get_vb_snip_ex (ConstContextP vctx, WordIndex vb_node_index, pSTRp(snip)); 
extern StrTextMegaLong ctx_get_vb_snip (ConstContextP vctx, WordIndex vb_node_index);
extern rom ctx_get_z_snip_ex (ConstContextP zctx, WordIndex z_node_index, pSTRp(snip));
extern StrTextMegaLong ctx_get_z_snip (ConstContextP zctx, WordIndex z_node_index); 

extern WordIndex ctx_get_word_index_by_snip (VBlockP vb, ContextP ctx, STRp(snip)); // PIZ

extern WordIndex ctx_get_ol_node_index_by_snip (VBlockP vb, ContextP ctx, STRp(snip)); 

extern CharIndex ctx_get_char_index_of_snip (ContextP zctx, STRp(snip), bool soft_fail); // ZIP

extern void ctx_initialize_predefined_ctxs (ContextArray contexts, DataType dt, DictIdtoDidMap d2d_map, Did *num_contexts);
extern void ctx_zip_init_promiscuous (ContextP zctx);
extern void ctx_zip_z_data_exist (ContextP ctx);

extern void ctx_shorten_unused_dict_words (Did did_i);
extern void ctx_piz_initialize_zctxs (void);
extern void ctx_read_all_counts (void);
extern void ctx_compress_counts (void);
extern void ctx_read_all_subdicts (void);
extern void ctx_compress_subdicts (void);
extern rom ctx_get_snip_with_largest_count (Did did_i, int64_t *count);
extern void ctx_populate_zf_ctx_from_contigs (Did dst_did_i, ConstContigPkgP ctgs);
extern WordIndex ctx_populate_zf_ctx (Did dst_did_i, STRp (snip), WordIndex ref_index);

extern void ctx_dump_binary (VBlockP vb, ContextP ctx, bool local);

StrText ctx_tag_name_ex (ConstContextP ctx);

// Needed only if we intend to call ctx_set_rollback again without advancing vb->rback_id
static inline void ctx_unset_rollback (ContextP ctx) { ctx->rback_id = 1; }

extern void ctx_rollback (VBlockP vb, ContextP ctx, bool override_id);

static inline bool ctx_can_have_singletons (ContextP ctx) 
    { return (ctx->ltype == LT_SINGLETON) && !ctx->no_stons && 
              (ctx->flags.store != STORE_INDEX) && !ctx->flags.all_the_same; } 

// returns true if dict_id was *previously* segged on this line, and we stored a valid last_value (int or float)
#define ctx_has_value_in_line_(vb, ctx) ((ctx)->last_line_i == (vb)->line_i)
#define ctx_has_value_in_prev_line_(vb, ctx) ((ctx)->last_line_i != NO_LINE && (ctx)->last_line_i+1 == (vb)->line_i)

static inline bool ctx_has_value_in_line_do (VBlockP vb, DictId dict_id, ContextP *p_ctx /* optional out */) 
{ 
    ContextP ctx = ECTX (dict_id);
    if (p_ctx) *p_ctx = ctx;
    return ctx && ctx_has_value_in_line_(vb, ctx);
}
#define ctx_has_value_in_line(vb, dict_id, p_ctx) ctx_has_value_in_line_do ((VBlockP)(vb), (DictId)(dict_id), (p_ctx))

static inline void ctx_set_last_value (VBlockP vb, ContextP ctx, ValueType last_value)
{
    ctx->last_value    = last_value;
    ctx->last_line_i   = vb->line_i;
    ctx->last_sample_i = vb->sample_i; // used for VCF/FORMAT. otherwise meaningless but harmless.
}

#define ENCOUNTERED(line_i) (-(int32_t)(line_i) - 2) 

// set encountered if not already ctx_set_last_value (encounted = seen, but without setting last_value)
static inline void ctx_set_encountered (VBlockP vb, ContextP ctx)
{
    if (ctx->last_line_i != vb->line_i || ctx->last_sample_i != vb->sample_i) // not already ctx_set_last_value in this line/sample
        ctx->last_line_i = ENCOUNTERED (vb->line_i); 

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
    decl_ctx (did_i); 
    return ((ctx->last_line_i == vb->line_i) || (ctx->last_line_i == ENCOUNTERED (vb->line_i))); 
}

static inline bool ctx_encountered_in_prev_line (VBlockP vb, Did did_i) 
{ 
    decl_ctx (did_i); 
    return vb->line_i && ((ctx->last_line_i == vb->line_i-1) || (ctx->last_line_i == ENCOUNTERED (vb->line_i-1))); 
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
    decl_ctx (did_i);
    return ctx_encountered_in_line(vb, ctx->did_i) && ctx->last_sample_i == vb->sample_i; 
}

static inline bool ctx_encountered_by_dict_id (VBlockP vb, DictId dict_id, ContextP *p_ctx) 
    { return ctx_encountered_in_line_by_dict_id (vb, dict_id, p_ctx) && (*p_ctx)->last_sample_i == vb->sample_i; }
    
static inline bool ctx_has_value (VBlockP vb, Did did_i) 
{   
    decl_ctx (did_i);
    return ctx_has_value_in_line_(vb, ctx) && ctx->last_sample_i == vb->sample_i; 
}

static inline bool ctx_has_value_by_dict_id (VBlockP vb, DictId dict_id, ContextP *p_ctx) 
    { return ctx_has_value_in_line_do (vb, dict_id, p_ctx) && (*p_ctx)->last_sample_i == vb->sample_i; }

extern uint64_t ctx_get_ctx_group_z_len (VBlockP vb, Did group_did_i);

typedef enum { KR_KEEP, KR_REMOVE } CtxKeepRemove;
extern void ctx_declare_winning_group (Did winning_group_did_i, Did losing_group_did_i, Did new_st_did_i);

#define T(cond, did_i) ((cond) ? (did_i) : DID_NONE)
extern void ctx_set_store (VBlockP vb, int store_type, ...);
extern void ctx_set_no_stons (VBlockP vb, ...);
extern void ctx_set_same_line (VBlockP vb, ...);
extern void ctx_set_store_per_line (VBlockP vb, ...);
extern void ctx_set_dyn_int (VBlockP vb, ...);
extern void ctx_set_ltype (VBlockP vb, int ltype, ...);
extern void ctx_consolidate_stats (VBlockP vb, int parent, ...);
extern void ctx_consolidate_statsN(VBlockP vb, Did parent, Did first_dep, unsigned num_deps);
extern void ctx_consolidate_statsA(VBlockP vb, Did parent, ContextP ctxs[], unsigned num_deps);
extern void ctx_consolidate_stats_(VBlockP vb, ContextP parent_ctx, ContainerP con);
extern ContextP buf_to_ctx (ContextArray ca, ConstBufferP buf);
extern void ctx_show_zctx_big_consumers (FILE *out);

extern rom dyn_type_name (DynType dyn_type);

extern SORTER (sort_by_dict_id);

#define for_zctx for (ContextP zctx=ZCTX(0), fc_after=ZCTX(z_file->num_contexts); zctx < fc_after; zctx++) 
#define for_vctx for (ContextP vctx=CTX(0),  fc_after=CTX(vb->num_contexts);      vctx < fc_after; vctx++) 
#define for_ctx  for (ContextP ctx=CTX(0),   fc_after=CTX(vb->num_contexts);      ctx  < fc_after; ctx++) 
#define for_zctx_that for (ContextP zctx=ZCTX(0), fc_after=ZCTX(z_file->num_contexts); zctx < fc_after; zctx++) if /* for each context for which the condition is true */
#define for_vctx_that for (ContextP vctx=CTX(0),  fc_after=CTX(vb->num_contexts);      vctx < fc_after; vctx++) if 
#define for_ctx_that  for (ContextP ctx=CTX(0),   fc_after=CTX(vb->num_contexts);      ctx  < fc_after; ctx++)  if 
