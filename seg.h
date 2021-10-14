// ------------------------------------------------------------------
//   seg.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SEGREGATE_INCLUDED
#define SEGREGATE_INCLUDED

#include <stdint.h>
#include "genozip.h"
#include "sections.h"
#include "context.h"
#include "container.h"

typedef enum { ERR_SEG_NO_ERROR=0, ERR_SEG_OUT_OF_RANGE, ERR_SEG_NOT_INTEGER } SegError;

extern void seg_all_data_lines (VBlockP vb); 

typedef enum { GN_FORBIDEN, GN_SEP, GN_IGNORE } GetNextAllow;
extern const char *seg_get_next_item (void *vb, const char *str, int *str_len, 
                                      GetNextAllow newline, GetNextAllow tab, GetNextAllow space,
                                      unsigned *len, char *separator, bool *has_13, // out
                                      const char *item_name);
extern const char *seg_get_next_line (void *vb_, const char *str, int *str_len, unsigned *len, bool must_have_newline, bool *has_13 /* out */, const char *item_name);

extern WordIndex seg_by_ctx_ex (VBlockP vb, const char *snip, unsigned snip_len, ContextP ctx, uint32_t add_bytes, bool *is_new);
static inline WordIndex seg_by_ctx (VBlockP vb, STRp(snip), ContextP ctx, unsigned add_bytes)                      { return seg_by_ctx_ex (vb, STRa(snip), ctx, add_bytes, NULL); }
static inline WordIndex seg_by_dict_id (VBlockP vb, STRp(snip), DictId dict_id, unsigned add_bytes)                { return seg_by_ctx_ex (vb, STRa(snip), ctx_get_ctx (vb, dict_id), add_bytes, NULL); }
static inline WordIndex seg_by_did_i_ex (VBlockP vb, STRp(snip), DidIType did_i, unsigned add_bytes, bool *is_new) { return seg_by_ctx_ex (vb, STRa(snip), CTX(did_i), add_bytes, is_new); }
static inline WordIndex seg_by_did_i (VBlockP vb, STRp(snip), DidIType did_i, unsigned add_bytes)                  { return seg_by_ctx_ex (vb, STRa(snip), CTX(did_i), add_bytes, NULL); }

extern WordIndex seg_known_node_index (VBlockP vb, ContextP ctx, WordIndex node_index, unsigned add_bytes);
#define seg_duplicate_last(vb,ctx,add_bytes) seg_known_node_index ((VBlockP)(vb), (ctx), *LASTENT (uint32_t, (ctx)->b250), (add_bytes))

extern WordIndex seg_integer_do (VBlockP vb, DidIType did_i, int64_t n, unsigned add_bytes); // segs integer as normal textual snip
#define seg_integer(vb,did_i,n,add_sizeof_n) seg_integer_do((VBlockP)(vb), (did_i), (n), (add_sizeof_n) ? sizeof(n) : 0)

extern WordIndex seg_self_delta (VBlockP vb, ContextP ctx, int64_t value, uint32_t value_str_len);

extern void seg_simple_lookup (VBlockP vb, ContextP ctx, unsigned add_bytes);
extern bool seg_integer_or_not (VBlockP vb, ContextP ctx, const char *this_value, unsigned this_value_len, unsigned add_bytes); // segs integer in local if possible
extern bool seg_integer_or_not_cb (VBlockP vb, ContextP ctx, STRp(int_str), uint32_t repeat);
extern bool seg_float_or_not (VBlockP vb, ContextP ctx, const char *this_value, unsigned this_value_len, unsigned add_bytes);

#define MAX_POS_DELTA 32000    // the max delta (in either direction) that we will put in a dictionary - above this it goes to random_pos. This number can be changed at any time without affecting backward compatability - it is used only by ZIP, not PIZ
#define SPF_BAD_SNIPS_TOO   1  // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
#define SPF_ZERO_IS_BAD     2  // whether 0 is considered a bad POS (if true and POS is 0, to be handled according to seg_bad_snips_too)
#define SPF_UNLIMITED_DELTA 4  // Always use delta, even if larger than MAX_POS_DELTA
#define SPF_NO_DELTA        8  // All integer data goes into local
extern PosType seg_pos_field (VBlockP vb, DidIType snip_did_i, DidIType base_did_i, unsigned opt, 
                              char missing, const char *pos_str, unsigned pos_len, PosType this_pos, unsigned add_bytes);
extern bool seg_pos_field_cb (VBlockP vb, ContextP ctx, const char *pos_str, unsigned pos_len, uint32_t repeat);

extern void seg_id_field_do (VBlockP vb, ContextP ctx, STRp(id_snip));
#define seg_id_field(vb, ctx, id_snip, id_snip_len, account_for_separator) \
    do { seg_id_field_do(VB, (ctx), (id_snip), (id_snip_len)); (ctx)->txt_len += !!(account_for_separator); } while(0)
extern bool seg_id_field_cb (VBlockP vb, ContextP ctx, STRp(id_snip), uint32_t repeat);

#define MAX_ARRAY_ITEMS 36
extern Container seg_initialize_container_array_do (DictId dict_id, bool type_1_items, bool comma_sep);
#define seg_initialize_container_array(dict_id, type_1_items, comma_sep) seg_initialize_container_array_do ((DictId)dict_id, type_1_items, comma_sep)

extern void seg_add_to_local_text   (VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len, unsigned add_bytes);
extern void seg_add_to_local_fixed  (VBlockP vb, ContextP ctx, const void *data, unsigned data_len);
extern void seg_add_to_local_uint8  (VBlockP vb, ContextP ctx, uint8_t  value, unsigned add_bytes);
extern void seg_add_to_local_uint (VBlockP vb, ContextP ctx, uint32_t value, unsigned add_bytes);

extern WordIndex seg_delta_vs_other (VBlockP vb, Context *ctx, Context *other_ctx, const char *value, unsigned value_len, int64_t max_delta);

extern WordIndex seg_array (VBlockP vb, ContextP container_ctx, DidIType stats_conslidation_did_i, const char *value, int32_t value_len, char sep, char subarray_sep, bool use_integer_delta, bool store_int_in_local, bool items_are_id);

typedef bool (*SegCallback) (VBlockP vb, ContextP ctx, STRp(value), uint32_t repeat);
extern int32_t seg_array_of_struct (VBlockP vb, ContextP ctx, MediumContainer con, const char *snip, unsigned snip_len, const SegCallback *callbacks);

typedef void (*SegOptimize)(const char **snip, unsigned *snip_len, char *space_for_new_str);

extern void seg_prepare_snip_other_do (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int32_t parameter, 
                                       char *snip, unsigned *snip_len /* in / out */);
#define seg_prepare_snip_other(snip_code, other_dict_id, has_parameter, parameter, snip) \
    snip##_len = sizeof (snip);\
    seg_prepare_snip_other_do ((snip_code), (DictId)(other_dict_id), (has_parameter), (parameter), (snip), &snip##_len)

#define seg_prepare_snip_special_other(special_code, snip, other_dict_id) do { \
    snip[0]=SNIP_SPECIAL; snip##_len=sizeof(snip)-1; \
    seg_prepare_snip_other_do ((special_code), (DictId)(other_dict_id), 0, 0, &snip[1], &snip##_len); \
    snip##_len++;\
} while(0)

bool seg_set_last_txt (VBlockP vb, ContextP ctx, const char *value, unsigned value_len, StoreType store_type);

extern void seg_create_rollback_point (VBlockP vb, unsigned num_ctxs, ...); // list of did_i
extern void seg_rollback (VBlockP vb);

// ------------------
// Seg utilities
// ------------------

// TAB separator between fields

#define FIELD(f) \
    const char *f##_str __attribute__((unused)) = field_start;  \
    unsigned    f##_len __attribute__((unused)) = field_len

#define GET_NEXT_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_FORBIDEN, GN_SEP, GN_IGNORE, &field_len, &separator, NULL, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM(f) \
    GET_NEXT_ITEM (f); \
    seg_by_did_i (VB, field_start, field_len, f, field_len+1)

#define GET_LAST_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_FORBIDEN, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_LAST_ITEM(f) do \
    GET_LAST_ITEM (f);\
    seg_by_did_i (VB, field_start, field_len, f, field_len+1)

#define GET_MAYBE_LAST_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_MAYBE_LAST_ITEM(f)  \
    GET_MAYBE_LAST_ITEM (f); \
    seg_by_did_i (VB, field_start, field_len, f, field_len+1)

// SPACE separator between fields

#define GET_NEXT_ITEM_SP(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_FORBIDEN, GN_SEP, GN_SEP, &field_len, &separator, NULL, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM_SP(f) ({ \
    GET_NEXT_ITEM_SP (f); \
    seg_by_did_i (VB, field_start, field_len, f, field_len+1); })

#define GET_LAST_ITEM_SP(f)  \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_FORBIDEN, GN_FORBIDEN, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_LAST_ITEM_SP(f)  \
    GET_LAST_ITEM_SP (f); \
    seg_by_did_i (VB, field_start, field_len, f, field_len+1)

#define GET_MAYBE_LAST_ITEM_SP(f)  \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_SEP, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_MAYBE_LAST_ITEM_SP(f)  \
    GET_MAYBE_LAST_ITEM_SP (f); \
    seg_by_did_i (VB, field_start, field_len, f, field_len+1)

// NEWLINE separator

#define GET_NEXT_ITEM_NL(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_IGNORE, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM_NL(f)  \
    GET_NEXT_ITEM_NL (f); \
    seg_by_did_i (VB, field_start, field_len, f, field_len+1);

#define SEG_EOL(f,account_for_ascii10) do { seg_by_did_i (VB, *(has_13) ? "\r\n" : "\n", 1 + *(has_13), (f), (account_for_ascii10) + *(has_13)); } while (0)

#define ASSSEG(condition, p_into_txt, format, ...) \
    ASSINP (condition, "Error in file %s: "format "\n\nvb_line_i:%"PRIu64" vb_i:%u pos_in_vb: %"PRIi64" pos_in_file: %"PRIi64\
                       "\nvb pos in file (0-based):%"PRIu64" - %"PRIu64" (length %"PRIu64")" \
                       "\n%d characters before to %d characters after (in quotes): \"%.*s\""\
                       "\n%d characters before to %d characters after (in quotes): \"%.*s\""\
                       "\nTo get vblock: %s %s | head -c %"PRIu64" | tail -c %u > vb.%u%s"\
                       "\nDumped bad vblock from memory: %s", \
                                     txt_name, __VA_ARGS__, vb->line_i, vb->vblock_i, \
            /* pos_in_vb:         */ (PosType)(p_into_txt ? (p_into_txt - vb->txt_data.data) : -1), \
            /* pos_in_file:       */ (PosType)(p_into_txt ? (vb->vb_position_txt_file + (p_into_txt - vb->txt_data.data)) : -1),\
            /* vb start pos file: */ vb->vb_position_txt_file, \
            /* vb end pos file:   */ vb->vb_position_txt_file + vb->txt_data.len-1, \
            /* vb length:         */ vb->txt_data.len,\
            /* +- 30 char snip    */\
            /* chars before:      */ p_into_txt ? MIN_(30, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN_(30, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN_(p_into_txt+31, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX_(p_into_txt-30, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX_(p_into_txt-30, vb->txt_data.data) : "(inaccessible)"),\
            /* +- 2 char snip     */\
            /* chars before:      */ p_into_txt ? MIN_(2, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN_(2, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN_(p_into_txt+3, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX_(p_into_txt-2, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX_(p_into_txt-3, vb->txt_data.data) : "(inaccessible)"),\
            /* head/tail params:  */ codec_args[txt_file->codec].viewer, txt_name, vb->vb_position_txt_file + vb->txt_data.len, (uint32_t)vb->txt_data.len,\
            /* output filename:   */ vb->vblock_i, file_plain_ext_by_dt (vb->data_type),\
            /* dump filename:     */ txtfile_dump_vb (VB, txt_name))

#define ASSSEG0(condition, p_into_txt, err_str) ASSSEG (condition, p_into_txt, err_str "%s", "")

#define ABOSEG(p_into_txt, format, ...) ASSSEG(false, p_into_txt, format, __VA_ARGS__)

#define ABOSEG0(p_into_txt, err_str) ABOSEG(false, p_into_txt, format, err_str "%s", "")

#endif