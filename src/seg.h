// ------------------------------------------------------------------
//   seg.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include <stdint.h>
#include "genozip.h"
#include "sections.h"
#include "context.h"
#include "container.h"

typedef enum { ERR_SEG_NO_ERROR=0, ERR_SEG_OUT_OF_RANGE, ERR_SEG_NOT_INTEGER } SegError;

extern void seg_all_data_lines (VBlockP vb); 

typedef enum { GN_FORBIDEN, GN_SEP, GN_IGNORE } GetNextAllow;
extern rom seg_get_next_item (void *vb, rom str, int *str_len, GetNextAllow newline, GetNextAllow tab, GetNextAllow space,
                              unsigned *len, char *separator, bool *has_13, // out
                              rom item_name);
extern rom seg_get_next_line (void *vb_, rom str, int *str_len, unsigned *len, bool must_have_newline, bool *has_13 /* out */, rom item_name);

extern WordIndex seg_by_ctx_ex (VBlockP vb, STRp(snip), ContextP ctx, uint32_t add_bytes, bool *is_new);
static inline WordIndex seg_by_ctx (VBlockP vb, STRp(snip), ContextP ctx, unsigned add_bytes)                      { return seg_by_ctx_ex (vb, STRa(snip), ctx, add_bytes, NULL); }
static inline WordIndex seg_by_dict_id (VBlockP vb, STRp(snip), DictId dict_id, unsigned add_bytes)                { return seg_by_ctx_ex (vb, STRa(snip), ctx_get_ctx (vb, dict_id), add_bytes, NULL); }
static inline WordIndex seg_by_did_i_ex (VBlockP vb, STRp(snip), Did did_i, unsigned add_bytes, bool *is_new) { return seg_by_ctx_ex (vb, STRa(snip), CTX(did_i), add_bytes, is_new); }
static inline WordIndex seg_by_did (VBlockP vb, STRp(snip), Did did_i, unsigned add_bytes)                  { return seg_by_ctx_ex (vb, STRa(snip), CTX(did_i), add_bytes, NULL); }

extern WordIndex seg_known_node_index (VBlockP vb, ContextP ctx, WordIndex node_index, unsigned add_bytes);
extern WordIndex seg_duplicate_last (VBlockP vb, ContextP ctx, unsigned add_bytes);

extern void seg_integer (VBlockP vb, ContextP ctx, int64_t n, bool with_lookup, unsigned add_bytes);

extern WordIndex seg_integer_as_text_do (VBlockP vb, ContextP ctx, int64_t n, unsigned add_bytes); // segs integer as normal textual snip
#define seg_integer_as_text(vb,did_i,n,add_sizeof_n) seg_integer_as_text_do((VBlockP)(vb), &vb->contexts[did_i], (n), (add_sizeof_n) ? sizeof(n) : 0)

extern WordIndex seg_self_delta (VBlockP vb, ContextP ctx, int64_t value, char format, uint32_t add_bytes);

extern void seg_simple_lookup (VBlockP vb, ContextP ctx, unsigned add_bytes);
extern void seg_lookup_with_length (VBlockP vb, ContextP ctx, int32_t length, unsigned add_bytes);

extern bool seg_integer_or_not (VBlockP vb, ContextP ctx, STRp(this_value), unsigned add_bytes); // segs integer in local if possible
extern bool seg_integer_or_not_cb (VBlockP vb, ContextP ctx, STRp(int_str), uint32_t repeat);
extern bool seg_float_or_not (VBlockP vb, ContextP ctx, STRp(this_value), unsigned add_bytes);
extern void seg_numeric_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned field_width, unsigned add_bytes);

#define MAX_POS_DELTA   16384  // the max delta (in either direction) that we will put in a dictionary - above this it goes to local. This number can be changed at any time without affecting backward compatability - it is used only by ZIP, not PIZ
#define SPF_BAD_SNIPS_TOO   1  // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
#define SPF_ZERO_IS_BAD     2  // whether 0 is considered a bad POS (if true and POS is 0, to be handled according to seg_bad_snips_too)
#define SPF_UNLIMITED_DELTA 4  // Always use delta, even if larger than MAX_POS_DELTA
#define SPF_NO_DELTA        8  // All integer data goes into local
extern PosType seg_pos_field (VBlockP vb, Did snip_did_i, Did base_did_i, unsigned opt, 
                              char missing, STRp(pos_str), PosType this_pos, unsigned add_bytes);
extern bool seg_pos_field_cb (VBlockP vb, ContextP ctx, STRp(pos_str), uint32_t repeat);

extern void seg_id_field_init (ContextP ctx);
extern void seg_id_field_do (VBlockP vb, ContextP ctx, STRp(id_snip));
#define seg_id_field(vb, ctx, id, id_len, account_for_separator) \
    ({ seg_id_field_do(VB, (ctx), (id), (id_len)); (ctx)->txt_len += !!(account_for_separator); })
extern void seg_id_field2 (VBlockP vb, ContextP ctx, STRp(id), unsigned add_bytes);

extern bool seg_id_field_cb (VBlockP vb, ContextP ctx, STRp(id), uint32_t repeat);
extern bool seg_id_field2_cb (VBlockP vb, ContextP ctx, STRp(id), uint32_t repeat);

typedef enum { LOOKUP_NONE, LOOKUP_SIMPLE, LOOKUP_WITH_LENGTH } Lookup;
extern void seg_add_to_local_fixed_do (VBlockP vb, ContextP ctx, STRp(data), bool add_nul, Lookup lookup_type, bool is_singleton, unsigned add_bytes);

static inline void seg_add_to_local_text (VBlockP vb, ContextP ctx, STRp(snip), Lookup lookup_type, unsigned add_bytes) 
{ 
#ifdef DEBUG
    ASSERT (segconf.running || ctx->no_stons, "%s: Expecting ctx->no_stons=true in ctx=%s", LN_NAME, ctx->tag_name);
#endif
    seg_add_to_local_fixed_do (vb, ctx, STRa(snip), lookup_type != LOOKUP_WITH_LENGTH, lookup_type, false, add_bytes); 
}

static inline void seg_add_to_local_fixed (VBlockP vb, ContextP ctx, STRp(data), Lookup lookup_type, unsigned add_bytes)
    { seg_add_to_local_fixed_do (vb, ctx, STRa(data), false, lookup_type, false, add_bytes); }

extern void seg_integer_fixed (VBlockP vb, Context *ctx, void *number, bool with_lookup, unsigned add_bytes);

// requires setting ltype=LT_DYN_INT* in seg_initialize, but not need to set ltype as it will be set in zip_resize_local
static inline void seg_add_to_local_resizable (VBlockP vb, ContextP ctx, int64_t value, unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || ctx->ltype==LT_DYN_INT || ctx->ltype==LT_DYN_INT_h || ctx->ltype==LT_DYN_INT_H, 
            "%s: Expecting ctx->ltype=%s to be LT_DYN_INT* in ctx=%s ", LN_NAME, lt_name(ctx->ltype), ctx->tag_name);
#endif

    // TO DO: find a way to better estimate the size, see b250_per_line
    buf_alloc (vb, &ctx->local, 1, vb->lines.len, int64_t, CTX_GROWTH, "contexts->local");
    BNXT (int64_t, ctx->local) = value;
    if (add_bytes) ctx->txt_len += add_bytes;
    ctx->local_num_words++;
}

extern WordIndex seg_delta_vs_other_do (VBlockP vb, ContextP ctx, ContextP other_ctx, STRp(value), int64_t value_n, int64_t max_delta, unsigned add_bytes);
static inline WordIndex seg_delta_vs_other (VBlockP vb, ContextP ctx, ContextP other_ctx, STRp(value))
    { return seg_delta_vs_other_do (vb, ctx, other_ctx, STRa(value), ctx->last_value.i, -1, value_len); }

extern void seg_diff (VBlockP vb, ContextP ctx, ContextP base_ctx, STRp(value), bool entire_snip_if_same, unsigned add_bytes);

extern WordIndex seg_array (VBlockP vb, ContextP container_ctx, Did stats_conslidation_did_i, rom value, int32_t value_len, char sep, char subarray_sep, bool use_integer_delta, StoreType store_in_local, DictId arr_dict_id, int add_bytes);

typedef bool (*SegCallback) (VBlockP vb, ContextP ctx, STRp(value), uint32_t repeat);
extern int32_t seg_array_of_struct (VBlockP vb, ContextP ctx, MediumContainer con, STRp(snip), const SegCallback *callbacks, unsigned add_bytes);
extern bool seg_do_nothing_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep);

extern bool seg_by_container (VBlockP vb, ContextP ctx, ContainerP con, STRp(value), STRp(container_snip), bool normal_seg_if_fail, unsigned add_bytes);

extern void seg_prepare_snip_other_do (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int64_t int_param, char char_param,
                                       char *snip, unsigned *snip_len /* in / out */);
#define seg_prepare_snip_other(snip_code, other_dict_id, has_parameter, parameter, snip) \
    snip##_len = sizeof (snip);\
    seg_prepare_snip_other_do ((snip_code), (other_dict_id), (has_parameter), (parameter), 0, (snip), &snip##_len)

#define seg_prepare_snip_other_char(snip_code, other_dict_id, char_param, snip) \
    ({ snip##_len = sizeof (snip);\
       seg_prepare_snip_other_do ((snip_code), (other_dict_id), true, 0, (char_param), (snip), &snip##_len); })

#define seg_prepare_snip_other_chari(snip_code, other_dict_id, char_param, snip, i) \
    ({ snip##_lens[i] = sizeof (snip##s);\
       seg_prepare_snip_other_do ((snip_code), (other_dict_id), true, 0, (char_param), (snip##s)[i], &snip##_lens[i]); })

#define seg_prepare_snip_special_other(special_code, other_dict_id, snip, char_param) \
    ({ snip[0]=SNIP_SPECIAL; snip##_len=sizeof(snip)-1; \
       seg_prepare_snip_other_do ((special_code), (other_dict_id), char_param, 0, char_param, &snip[1], &snip##_len); \
       snip##_len++; })

extern void seg_prepare_multi_dict_id_special_snip (uint8_t special_code, unsigned num_dict_ids, DictId *dict_ids, char *out_snip, unsigned *out_snip_len);

extern void seg_prepare_minus_snip_do (DictId dict_id_a, DictId dict_id_b, uint8_t special_code, char *snip, unsigned *snip_len);
#define seg_prepare_minus_snip(dt, dict_id_a, dict_id_b, snip) \
    ({ snip##_len = sizeof (snip);\
       seg_prepare_minus_snip_do ((dict_id_a), (dict_id_b), dt##_SPECIAL_MINUS, (snip), &snip##_len); })
#define seg_prepare_minus_snip_i(dt, dict_id_a, dict_id_b, snip, i) \
    ({ snip##_lens[i] = sizeof (snip##s[i]);\
       seg_prepare_minus_snip_do ((dict_id_a), (dict_id_b), dt##_SPECIAL_MINUS, snip##s[i], &snip##_lens[i]); })

static void inline seg_set_last_txt (VBlockP vb, ContextP ctx, STRp(value))
{
    bool is_value_in_txt_data = value >= B1STtxt &&
                                value <= BLST  (char, vb->txt_data);

    ctx->last_txt = (TxtWord){ .index = is_value_in_txt_data ? BNUMtxt (value) : INVALID_LAST_TXT_INDEX,
                               .len   = value_len };

    ctx_set_encountered (vb, ctx);
}

bool seg_set_last_txt_store_value (VBlockP vb, ContextP ctx, STRp(value), StoreType store_type);

extern void seg_create_rollback_point (VBlockP vb, ContextP *ctxs, unsigned num_ctxs, ...); // list of did_i
extern void seg_add_ctx_to_rollback_point (VBlockP vb, ContextP ctx);
extern void seg_rollback (VBlockP vb);

// Multiplexers
#define BASE64_DICT_ID_LEN 14
#define MULTIPLEXER(n_channels)                     \
struct __attribute__ ((__packed__)) {               \
    /* all 32b/64b fields are word-aligned */       \
    ContextP ctx;                                   \
    bool no_stons;                                  \
    uint8_t special_code;                           \
    uint16_t num_channels;                          \
    DictId dict_ids[n_channels];                    \
    ContextP channel_ctx[n_channels];               \
    uint32_t snip_len;                              \
    char snip[BASE64_DICT_ID_LEN * (n_channels)];   \
}                               
#define MUX ((MultiplexerP)mux)
#define MUX_CAPACITY(mux) (sizeof((mux).dict_ids)/sizeof(DictId)) // max number of channels this mux can contain
#define MUX_CHANNEL_CTX(mux) ((ContextP *)((mux)->dict_ids + (mux)->num_channels))
#define MUX_SNIP_LEN(mux)    (*(uint32_t*)(MUX_CHANNEL_CTX(mux) + (mux)->num_channels))
#define MUX_SNIP(mux)        ((char*)(&MUX_SNIP_LEN(mux) + 1))

typedef MULTIPLEXER(1000) *MultiplexerP;
typedef const MULTIPLEXER(1000) *ConstMultiplexerP;

typedef MULTIPLEXER(2) Multiplexer2, *Multiplexer2P;
typedef MULTIPLEXER(3) Multiplexer3, *Multiplexer3P;
typedef MULTIPLEXER(4) Multiplexer4, *Multiplexer4P;
typedef MULTIPLEXER(5) Multiplexer5, *Multiplexer5P;
typedef MULTIPLEXER(6) Multiplexer6, *Multiplexer6P;
typedef MULTIPLEXER(7) Multiplexer7, *Multiplexer7P;

extern void seg_mux_init (VBlockP vb, ContextP ctx, unsigned num_channels, uint8_t special_code, bool no_stons, MultiplexerP mux, rom channel_letters);
extern ContextP seg_mux_get_channel_ctx (VBlockP vb, Did did_i, MultiplexerP mux, uint32_t channel_i);

// --------------------
// handling binary data
// --------------------

// loading a Little Endian uint32_t from an unaligned buffer
#define GET_UINT8(p)   ((uint8_t)(((uint8_t*)(p))[0]))
#define GET_UINT16(p)  ((uint16_t)(((uint8_t*)(p))[0] | (((uint8_t*)(p))[1] << 8)))
#define GET_UINT32(p)  ((uint32_t)(((uint8_t*)(p))[0] | (((uint8_t*)(p))[1] << 8) | (((uint8_t*)(p))[2] << 16) | (((uint8_t*)(p))[3] << 24)))
#define GET_UINT64(p)  ((uint64_t)(((uint8_t*)(p))[0] | ((uint64_t)((uint8_t*)(p))[1] << 8) | ((uint64_t)((uint8_t*)(p))[2] << 16) | ((uint64_t)((uint8_t*)(p))[3] << 24) | ((uint64_t)((uint8_t*)(p))[4] << 32) | ((uint64_t)((uint8_t*)(p))[5] << 40) | ((uint64_t)((uint8_t*)(p))[6] << 48) | ((uint64_t)((uint8_t*)(p))[7] << 56)))
#define GET_FLOAT32(p) ({ union { uint32_t i; float f; } n= {.i = GET_UINT32(p)}; n.f; })

// getting integers from the BAM data
#define NEXT_UINT8    ({ uint8_t  value = GET_UINT8   (next_field); next_field += sizeof (uint8_t ); value; })
#define NEXT_UINT16   ({ uint16_t value = GET_UINT16  (next_field); next_field += sizeof (uint16_t); value; })
#define NEXT_UINT32   ({ uint32_t value = GET_UINT32  (next_field); next_field += sizeof (uint32_t); value; })
#define NEXT_UINT64   ({ uint64_t value = GET_UINT64  (next_field); next_field += sizeof (uint64_t); value; })
#define NEXT_FLOAT32  ({ float    value = GET_FLOAT32 (next_field); next_field += sizeof (float);    value; })

#define NEXTP_UINT8   ({ uint8_t  value = GET_UINT8   (*next_field_p); *next_field_p += sizeof (uint8_t);  value; })
#define NEXTP_UINT16  ({ uint16_t value = GET_UINT16  (*next_field_p); *next_field_p += sizeof (uint16_t); value; })
#define NEXTP_UINT32  ({ uint32_t value = GET_UINT32  (*next_field_p); *next_field_p += sizeof (uint32_t); value; })
#define NEXTP_UINT64  ({ uint64_t value = GET_UINT64  (*next_field_p); *next_field_p += sizeof (uint64_t); value; })
#define NEXTP_FLOAT32 ({ uint32_t value = GET_FLOAT32 (*next_field_p); *next_field_p += sizeof (float);    value; })

// ------------------
// Seg utilities
// ------------------

// TAB separator between fields

#define FIELD(f) \
    rom f##_str __attribute__((unused)) = field_start;  \
    unsigned    f##_len __attribute__((unused)) = field_len

#define GET_NEXT_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_FORBIDEN, GN_SEP, GN_IGNORE, &field_len, &separator, NULL, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM(f) \
    GET_NEXT_ITEM (f); \
    seg_by_did (VB, field_start, field_len, f, field_len+1)

#define GET_LAST_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_FORBIDEN, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_LAST_ITEM(f) \
    GET_LAST_ITEM (f);\
    seg_by_did (VB, field_start, field_len, f, field_len+1)

#define GET_MAYBE_LAST_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_MAYBE_LAST_ITEM(f)  \
    GET_MAYBE_LAST_ITEM (f); \
    seg_by_did (VB, field_start, field_len, f, field_len+1)

// SPACE separator between fields

#define GET_NEXT_ITEM_SP(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_FORBIDEN, GN_SEP, GN_SEP, &field_len, &separator, NULL, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM_SP(f) \
    GET_NEXT_ITEM_SP (f); \
    seg_by_did (VB, field_start, field_len, f, field_len+1); 

#define GET_LAST_ITEM_SP(f)  \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_FORBIDEN, GN_FORBIDEN, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_LAST_ITEM_SP(f)  \
    GET_LAST_ITEM_SP (f); \
    seg_by_did (VB, field_start, field_len, f, field_len+1)

#define GET_MAYBE_LAST_ITEM_SP(f)  \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_SEP, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_MAYBE_LAST_ITEM_SP(f)  \
    GET_MAYBE_LAST_ITEM_SP (f); \
    seg_by_did (VB, field_start, field_len, f, field_len+1)

// NEWLINE separator

#define GET_NEXT_ITEM_NL(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_IGNORE, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM_NL(f)  \
    GET_NEXT_ITEM_NL (f); \
    seg_by_did (VB, field_start, field_len, f, field_len+1);

#define SEG_EOL(f,account_for_ascii10) ({ seg_by_did (VB, *(has_13) ? "\r\n" : "\n", 1 + *(has_13), (f), (account_for_ascii10) + *(has_13)); })

#define ASSSEG(condition, p_into_txt, format, ...)                                                                  \
    ASSINP ((condition), "%s %s/%s: "format "\nTechnical details: pos_in_vb: %"PRIi64" pos_in_file: %"PRIi64        \
                         "\n%d characters before to %d characters after (in quotes): \"%.*s\""                      \
                         "\nTo get vblock: genozip --biopsy %u %s",                                                 \
                                     str_time().s, txt_name, LN_NAME, __VA_ARGS__,                                  \
            /* pos_in_vb:         */ (PosType)(p_into_txt ? ((rom)p_into_txt - vb->txt_data.data) : -1),            \
            /* pos_in_file:       */ (PosType)(p_into_txt ? (vb->vb_position_txt_file + ((rom)p_into_txt - vb->txt_data.data)) : -1),\
            /* +- 10 char snip    */                                                                                \
            /* chars before:      */ p_into_txt ? MIN_(10, (unsigned)((rom)p_into_txt - vb->txt_data.data)) : -1,   \
            /* chars after:       */ p_into_txt ? MIN_(10, (unsigned)(vb->txt_data.data + vb->txt_data.len - (rom)p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN_((rom)p_into_txt+11, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX_((rom)p_into_txt-10, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && ((rom)p_into_txt >= vb->txt_data.data) && ((rom)p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX_((rom)p_into_txt-30, vb->txt_data.data) : "(inaccessible)"),            \
            /* biopsy parameters: */ vb->vblock_i, txt_name)

#define ASSSEG0(condition, p_into_txt, err_str) ASSSEG (condition, p_into_txt, err_str "%s", "")

#define ABOSEG(p_into_txt, format, ...) ASSSEG(false, p_into_txt, format, __VA_ARGS__)

#define ABOSEG0(p_into_txt, err_str) ABOSEG(false, p_into_txt, format, err_str "%s", "")
