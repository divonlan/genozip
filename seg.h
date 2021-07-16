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
extern const char *seg_get_next_line (void *vb_, const char *str, int *str_len, unsigned *len, bool *has_13 /* out */, const char *item_name);

extern WordIndex seg_by_ctx_do (VBlockP vb, const char *snip, unsigned snip_len, ContextP ctx, uint32_t add_bytes, bool *is_new);
#define seg_by_ctx(vb,snip,snip_len,ctx,add_bytes)               seg_by_ctx_do ((VBlockP)(vb), (snip), (snip_len), (ctx), (add_bytes), NULL)
#define seg_by_dict_id(vb,snip,snip_len,dict_id,add_bytes)       seg_by_ctx_do ((VBlockP)(vb), (snip), (snip_len), ctx_get_ctx ((vb), (dict_id)), (add_bytes), NULL)
#define seg_by_did_i_ex(vb,snip,snip_len,did_i,add_bytes,is_new) seg_by_ctx_do ((VBlockP)(vb), (snip), (snip_len), CTX(did_i), (add_bytes), (is_new))
#define seg_by_did_i(vb,snip,snip_len,did_i,add_bytes)           seg_by_ctx_do ((VBlockP)(vb), (snip), (snip_len), CTX(did_i), (add_bytes), NULL)

extern WordIndex seg_chrom_field (VBlockP vb, const char *chrom_str, unsigned chrom_str_len);

extern WordIndex seg_integer_do (VBlockP vb, DidIType did_i, int64_t n, unsigned add_bytes); // segs integer as normal textual snip
#define seg_integer(vb,did_i,n,add_sizeof_n) seg_integer_do((VBlockP)(vb), (did_i), (n), (add_sizeof_n) ? sizeof(n) : 0)

extern void seg_simple_lookup (VBlockP vb, ContextP ctx, unsigned add_bytes);
extern bool seg_integer_or_not_do (VBlockP vb, ContextP ctx, const char *this_value, unsigned this_value_len, unsigned add_bytes); // segs integer in local if possible
#define seg_integer_or_not(vb, ctx, this_value, this_value_len, add_bytes) seg_integer_or_not_do ((VBlockP)(vb), (ctx), (this_value), (this_value_len), (add_bytes))

extern bool seg_float_or_not_do (VBlockP vb, ContextP ctx, const char *this_value, unsigned this_value_len, unsigned add_bytes);
#define seg_float_or_not(vb, ctx, this_value, this_value_len, add_bytes) seg_float_or_not_do ((VBlockP)(vb), (ctx), (this_value), (this_value_len), (add_bytes))

#define MAX_POS_DELTA 32000   // the max delta (in either direction) that we will put in a dictionary - above this it goes to random_pos. This number can be changed at any time without affecting backward compatability - it is used only by ZIP, not PIZ
#define SPF_BAD_SNIPS_TOO  1  // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
#define SPF_ZERO_IS_BAD    2  // whether 0 is considered a bad POS (if true and POS is 0, to be handled according to seg_bad_snips_too)
#define SPF_UNLIMIED_DELTA 4  // Always use delta, even if larger than MAX_POS_DELTA
extern PosType seg_pos_field (VBlockP vb, DidIType snip_did_i, DidIType base_did_i, unsigned opt, 
                              char missing, const char *pos_str, unsigned pos_len, PosType this_pos, unsigned add_bytes);

extern void seg_id_field_do (VBlockP vb, DictId dict_id, const char *id_snip, unsigned id_snip_len, bool account_for_separator);
#define seg_id_field(vb, dict_id,id_snip, id_snip_len, account_for_separator) \
    seg_id_field_do((VBlockP)vb, (DictId)dict_id, (id_snip), (id_snip_len), (account_for_separator))

extern Container seg_initialize_container_array_do (DictId dict_id, bool type_1_items, bool comma_sep);
#define seg_initialize_container_array(dict_id, type_1_items, comma_sep) seg_initialize_container_array_do ((DictId)dict_id, type_1_items, comma_sep)

extern void seg_add_to_local_text   (VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len, unsigned add_bytes);
extern void seg_add_to_local_fixed  (VBlockP vb, ContextP ctx, const void *data, unsigned data_len);
extern void seg_add_to_local_uint8  (VBlockP vb, ContextP ctx, uint8_t  value, unsigned add_bytes);
extern void seg_add_to_local_uint16 (VBlockP vb, ContextP ctx, uint16_t value, unsigned add_bytes);
extern void seg_add_to_local_uint32 (VBlockP vb, ContextP ctx, uint32_t value, unsigned add_bytes);
extern void seg_add_to_local_uint64 (VBlockP vb, ContextP ctx, uint64_t value, unsigned add_bytes);

extern WordIndex seg_delta_vs_other_do (VBlockP vb, Context *ctx, Context *other_ctx, const char *value, unsigned value_len, int64_t max_delta);
#define seg_delta_vs_other(vb, ctx, other_ctx, value, value_len, max_delta) seg_delta_vs_other_do ((VBlockP)(vb), (ctx), (other_ctx), (value), (value_len), (max_delta))

extern WordIndex seg_array (VBlockP vb, ContextP container_ctx, DidIType stats_conslidation_did_i, const char *value, int32_t value_len, char sep, char subarray_sep, bool use_integer_delta, bool store_int_in_local);

typedef struct { bool slash, pipe, colon, dot, whitespace; /* seperators */ } SegCompoundArg; 
extern void seg_compound_field (VBlockP vb, ContextP field_ctx, const char *field, unsigned field_len, 
                                SegCompoundArg arg, unsigned nonoptimized_len, unsigned add_for_eol);

typedef void (*SegOptimize)(const char **snip, unsigned *snip_len, char *space_for_new_str);

extern void seg_prepare_snip_other_do (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int32_t parameter, 
                                       char *snip, unsigned *snip_len /* in / out */);
#define seg_prepare_snip_other(snip_code, other_dict_id, has_parameter, parameter, snip, snip_len) \
    seg_prepare_snip_other_do ((snip_code), (DictId)(other_dict_id), (has_parameter), (parameter), (snip), (snip_len))

bool seg_set_last_txt_do (VBlockP vb, ContextP ctx, const char *value, unsigned value_len, StoreType store_type);
#define seg_set_last_txt(vb, ctx, value, value_len, store_type) seg_set_last_txt_do ((VBlockP)(vb), (ctx), (value), (value_len), (store_type))

extern void seg_create_rollback_point (ContextP ctx);
extern void seg_rollback (VBlockP vb, ContextP ctx);

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
    seg_by_did_i (vb, field_start, field_len, f, field_len+1)

#define GET_LAST_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_FORBIDEN, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_LAST_ITEM(f) do \
    GET_LAST_ITEM (f);\
    seg_by_did_i (vb, field_start, field_len, f, field_len+1)

#define GET_MAYBE_LAST_ITEM(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_IGNORE, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_MAYBE_LAST_ITEM(f)  \
    GET_MAYBE_LAST_ITEM (f); \
    seg_by_did_i (vb, field_start, field_len, f, field_len+1)

// SPACE separator between fields

#define GET_NEXT_ITEM_SP(f) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_FORBIDEN, GN_SEP, GN_SEP, &field_len, &separator, NULL, #f); \
    FIELD (f)

#define SEG_NEXT_ITEM_SP(f)  \
    GET_NEXT_ITEM_SP (f); \
    seg_by_did_i (vb, field_start, field_len, f, field_len+1);

#define GET_LAST_ITEM_SP(f)  \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_FORBIDEN, GN_FORBIDEN, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_LAST_ITEM_SP(f)  \
    GET_LAST_ITEM_SP (f); \
    seg_by_did_i (vb, field_start, field_len, f, field_len+1)

#define GET_MAYBE_LAST_ITEM_SP(f)  \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, GN_SEP, GN_SEP, GN_SEP, &field_len, &separator, has_13, #f); \
    FIELD (f)

#define SEG_MAYBE_LAST_ITEM_SP(f)  \
    GET_MAYBE_LAST_ITEM_SP (f); \
    seg_by_did_i (vb, field_start, field_len, f, field_len+1)

#define SEG_EOL(f,account_for_ascii10) seg_by_did_i (vb, *(has_13) ? "\r\n" : "\n", 1 + *(has_13), (f), (account_for_ascii10) + *(has_13)); 

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
            /* chars before:      */ p_into_txt ? MIN (30, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN (30, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN (p_into_txt+31, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX (p_into_txt-30, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX (p_into_txt-30, vb->txt_data.data) : "(inaccessible)"),\
            /* +- 2 char snip     */\
            /* chars before:      */ p_into_txt ? MIN (2, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN (2, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN (p_into_txt+3, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX (p_into_txt-2, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX (p_into_txt-3, vb->txt_data.data) : "(inaccessible)"),\
            /* head/tail params:  */ codec_args[txt_file->codec].viewer, txt_name, vb->vb_position_txt_file + vb->txt_data.len, (uint32_t)vb->txt_data.len,\
            /* output filename:   */ vb->vblock_i, file_plain_ext_by_dt (vb->data_type),\
            /* dump filename:     */ txtfile_dump_vb ((VBlockP)vb, txt_name))

#define ASSSEG0(condition, p_into_txt, err_str) ASSSEG (condition, p_into_txt, err_str "%s", "")

#define ABOSEG(p_into_txt, format, ...) ASSSEG(false, p_into_txt, format, __VA_ARGS__)

#define ABOSEG0(p_into_txt, err_str) ABOSEG(false, p_into_txt, format, err_str "%s", "")

#endif