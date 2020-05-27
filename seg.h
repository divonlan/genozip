// ------------------------------------------------------------------
//   seg.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SEGREGATE_INCLUDED
#define SEGREGATE_INCLUDED

#include <stdint.h>
#include "genozip.h"
#include "section_types.h"
#include "dict_id.h"
#include "move_to_front.h"

extern void seg_all_data_lines (VBlockP vb); 

extern void seg_init_mapper (VBlockP vb, int field_i, BufferP mapper_buf, const char *name);

extern const char *seg_get_next_item (void *vb, const char *str, int *str_len, 
                                      bool allow_newline, bool allow_tab, bool allow_colon, 
                                      unsigned *len, char *separator, bool *has_13, // out
                                      const char *item_name);
extern const char *seg_get_next_line (void *vb_, const char *str, int *str_len, unsigned *len, bool *has_13 /* out */, const char *item_name);

extern uint32_t seg_by_ctx (VBlockP vb, const char *snip, unsigned snip_len, MtfContextP ctx, uint32_t add_bytes, bool *is_new);
#define seg_by_dict_id(vb,str,len,dict_id,add_bytes)       seg_by_ctx ((VBlockP)vb, str, len, mtf_get_ctx (vb, dict_id), add_bytes, NULL)
#define seg_by_did_i_ex(vb,str,len,did_i,add_bytes,is_new) seg_by_ctx ((VBlockP)vb, str, len, &vb->mtf_ctx[did_i], add_bytes, is_new);
#define seg_by_did_i(vb,str,len,did_i,add_bytes)           seg_by_ctx ((VBlockP)vb, str, len, &vb->mtf_ctx[did_i], add_bytes, NULL);

extern uint32_t seg_chrom_field (VBlockP vb, const char *chrom_str, unsigned chrom_str_len);

extern int64_t seg_scan_pos_snip (VBlockP vb, const char *snip, unsigned snip_len, bool allow_nonsense);

#define MAX_POS_DELTA 32000 // the max delta (in either direction) that we will put in a dictionary - above this it goes to random_pos. This number can be changed at any time without affecting backward compatability - it is used only by ZIP, not PIZ
extern void seg_pos_field (VBlockP vb, 
                           uint8_t snip_did_i,    // mandatory: the ctx the snip belongs to
                           uint8_t base_did_i,    // mandatory: base for delta
                           bool allow_non_number,      // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
                           const char *pos_str, unsigned pos_len, 
                           bool account_for_separator);

extern void seg_id_field (VBlockP vb, DictIdType dict_id, const char *id_snip, unsigned id_snip_len, bool account_for_separator);

typedef bool (*SegSpecialInfoSubfields)(VBlockP vb, MtfContextP ctx, const char **this_value, unsigned *this_value_len, char *optimized_snip);

extern void seg_structured_by_ctx (VBlockP vb, MtfContextP ctx, StructuredP st, unsigned add_bytes);
#define seg_structured_by_dict_id(vb,dict_id,st,add_bytes) seg_structured_by_ctx ((VBlockP)vb, mtf_get_ctx (vb, dict_id), st, add_bytes)

extern void seg_info_field (VBlockP vb, uint32_t *dl_info_mtf_i, BufferP iname_mapper_buf, uint8_t *num_info_subfields,
                            SegSpecialInfoSubfields seg_special_subfields,
                            const char *info_str, unsigned info_len, 
                            bool this_field_has_13, // this is the last field in the line, and it ends with a Windows-style \r\n - we account for it in txt_len
                            bool this_line_has_13); // this line ends with \r\n (this field may or may not be the last field) - we store this information as an info subfield for PIZ to recover

extern void seg_add_to_local_text   (VBlockP vb, MtfContextP ctx, const char *snip, unsigned snip_len, unsigned add_bytes);
extern void seg_add_to_local_fixed  (VBlockP vb, MtfContextP ctx, const void *data, unsigned data_len);
extern void seg_add_to_local_uint8  (VBlockP vb, MtfContextP ctx, uint8_t  value, unsigned add_bytes);
extern void seg_add_to_local_uint16 (VBlockP vb, MtfContextP ctx, uint16_t value, unsigned add_bytes);
extern void seg_add_to_local_uint32 (VBlockP vb, MtfContextP ctx, uint32_t value, unsigned add_bytes);
extern void seg_add_to_local_uint64 (VBlockP vb, MtfContextP ctx, uint64_t value, unsigned add_bytes);

extern void seg_initialize_compound_structured (VBlockP vb, char *name_template, Structured *st);
extern void seg_compound_field (VBlockP vb, MtfContextP field_ctx, const char *field, unsigned field_len, 
                                SubfieldMapperP mapper, Structured st, bool ws_is_sep, unsigned add_for_eol);

typedef void (*SegOptimize)(const char **snip, unsigned *snip_len, char *space_for_new_str);
extern void seg_array_field (VBlockP vb, DictIdType dict_id, const char *value, unsigned value_len, SegOptimize optimize);

extern void seg_prepare_snip_other (uint8_t snip_code, DictIdType other_dict_id, uint32_t lookup_len, 
                                    char *snip, unsigned *snip_len);

// ------------------
// Seg utilities
// ------------------

#define GET_NEXT_ITEM(item_name,has_13) \
    field_start = next_field; \
    next_field = seg_get_next_item (vb, field_start, &len, !!has_13, !has_13, false, &field_len, &separator, has_13, item_name);

// create extendent field contexts in the correct order of the fields
#define EXTENDED_FIELD_CTX(extended_field, dict_id_num) { \
    MtfContext *ctx = mtf_get_ctx (vb, (DictIdType)dict_id_num); \
    ASSERT (ctx->did_i == extended_field, "Error: expecting ctx->did_i=%u to be %u", ctx->did_i, extended_field); \
    dict_id_fields[ctx->did_i] = ctx->dict_id.num; \
}

#define SAFE_ASSIGN(reg,addr,char_val) /* we are careful to evaluate addr, char_val only once, less they contain eg ++ */ \
    char *__addr##reg = (char*)(addr); \
    char __save##reg = *__addr##reg; \
    *__addr##reg = (char_val);

#define SAFE_RESTORE(reg) *__addr##reg = __save##reg; 

#define SEG_EOL(f) seg_by_did_i (vb, has_13 ? "\r\n" : "\n", 1 + has_13, (f), 1 + has_13); \
                   *has_special_eol = *has_special_eol || has_13;

#define ASSSEG(condition, p_into_txt, format, ...) \
    ASSERT (condition, format "\nFile: %s vb_line_i:%u vb_i:%u pos_in_vb: %"PRIi64" pos_in_file: %"PRIi64\
                              "\nvb pos in file (0-based):%"PRIu64" - %"PRIu64" (length %"PRIu64")" \
                              "\n%d characters before to %d characters after (in quotes): \"%.*s\""\
                              "\n%d characters before to %d characters after (in quotes): \"%.*s\""\
                              "\nTo get vblock: %s %s | head -c %"PRIu64" | tail -c %"PRIu64 " > vb%s", \
            __VA_ARGS__, txt_name, vb->line_i, vb->vblock_i, \
            /* pos_in_vb:         */ (int64_t)(p_into_txt ? (p_into_txt - vb->txt_data.data) : -1), \
            /* pos_in_file:       */ (int64_t)(p_into_txt ? (vb->vb_position_txt_file + (p_into_txt - vb->txt_data.data)) : -1),\
            /* vb start pos file: */ vb->vb_position_txt_file, \
            /* vb end pos file:   */ vb->vb_position_txt_file + vb->txt_data.len-1, \
            /* vb length:         */ vb->txt_data.len,\
            /* +- 30 char snip    */\
            /* chars before:      */ p_into_txt ? MIN (30, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN (30, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN (p_into_txt+31, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX (p_into_txt-30, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX (p_into_txt-30, vb->txt_data.data) : "(inaccessible)"),\
            /* +- 2 char snip    */\
            /* chars before:      */ p_into_txt ? MIN (2, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN (2, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN (p_into_txt+3, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX (p_into_txt-2, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX (p_into_txt-3, vb->txt_data.data) : "(inaccessible)"),\
            /* head, tail params: */ file_viewer (txt_file), txt_name, vb->vb_position_txt_file + vb->txt_data.len, vb->txt_data.len,\
            /* plain file ext:    */ file_plain_ext_by_dt (vb->data_type))

#define ASSSEG0(condition, p_into_txt, err_str) ASSSEG (condition, p_into_txt, err_str "%s", "")

#define ABOSEG(p_into_txt, format, ...) ASSSEG(false, p_into_txt, format, __VA_ARGS__)

#define ABOSEG0(p_into_txt, err_str) ABOSEG(false, p_into_txt, format, err_str "%s", "")

#endif