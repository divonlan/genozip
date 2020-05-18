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

extern void seg_all_data_lines (VBlockP vb); 

extern void seg_init_mapper (VBlockP vb, int field_i, BufferP mapper_buf, const char *name);

extern const char *seg_get_next_item (void *vb, const char *str, int *str_len, bool allow_newline, bool allow_tab, bool allow_colon, 
                                      unsigned *len, char *separator, bool *has_13, // out
                                      const char *item_name);
extern const char *seg_get_next_line (void *vb_, const char *str, int *str_len, unsigned *len, bool *has_13 /* out */, const char *item_name);

extern void seg_store (VBlockP vb, 
                       bool *dst_is_spillover, uint32_t *dst_start, uint32_t *dst_len, // out
                       BufferP src_buf, uint32_t size, // Either src_buf OR size must be given
                       const char *limit_txt_data, // we cannot store in txt starting here. if NULL always allocates in txt_data_spillover
                       bool align32); // does start address need to be 32bit aligned to prevent aliasing issues

extern uint32_t seg_one_subfield (VBlockP vb, const char *str, unsigned len, DictIdType dict_id, uint32_t add_bytes);

extern uint32_t seg_one_snip (VBlockP vb, const char *str, unsigned len, int did_i, uint32_t add_bytes,
                              bool *is_new); // optional out

#define seg_one_field(vb,str,len,f, add_bytes) \
    seg_one_snip ((VBlockP)(vb), (str), (len), (f), add_bytes, NULL)

extern uint32_t seg_chrom_field (VBlockP vb, const char *chrom_str, unsigned chrom_str_len);

extern uint32_t seg_add_to_random_pos_data (VBlockP vb, const char *snip, unsigned snip_len, unsigned add_bytes, const char *field_name);

#define MAX_POS_DELTA 32000 // the max delta (in either direction) that we will put in a dictionary - above this it goes to random_pos. This number can be changed at any time without affecting backward compatability - it is used only by ZIP, not PIZ
extern int32_t seg_pos_field (VBlockP vb, int32_t last_pos, int32_t *last_pos_delta, bool allow_non_number, 
                              int pos_field, const char *pos_str, unsigned pos_len, 
                              const char *field_name, bool account_for_separator);

extern void seg_id_field (VBlockP vb, DictIdType dict_id, const char *id_snip, unsigned id_snip_len, bool account_for_separator);

typedef bool (*SegSpecialInfoSubfields)(VBlockP vb, MtfContextP ctx, const char **this_value, unsigned *this_value_len, char *optimized_snip);

extern void seg_info_field (VBlockP vb, uint32_t *dl_info_mtf_i, BufferP iname_mapper_buf, uint8_t *num_info_subfields,
                            SegSpecialInfoSubfields seg_special_subfields,
                            const char *info_str, unsigned info_len, 
                            bool this_field_has_13, // this is the last field in the line, and it ends with a Windows-style \r\n - we account for it in txt_len
                            bool this_line_has_13); // this line ends with \r\n (this field may or may not be the last field) - we store this information as an info subfield for PIZ to recover

extern void seg_add_to_local_text  (VBlockP vb, MtfContextP ctx, const char *snip, unsigned snip_len, unsigned add_bytes);
extern void seg_add_to_local_fixed (VBlockP vb, MtfContextP ctx, const void *data, unsigned data_len);

extern void seg_compound_field (VBlockP vb, MtfContextP field_ctx, const char *field, unsigned field_len, 
                                SubfieldMapperP mapper, DictIdType sf_dict_id, bool ws_is_sep, bool account_for_13);
                               
// ---------
// VCF Stuff
// ---------
extern const char *seg_vcf_data_line (VBlockP vb_, const char *field_start_line);
extern void seg_vcf_initialize  (VBlockP vb_);
extern void seg_vcf_complete_missing_lines (VBlockVCFP vb);

// ---------
// SAM Stuff
// ---------
extern const char *seg_sam_data_line (VBlockP vb_, const char *field_start_line);
extern uint32_t seg_sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len);
extern uint32_t seg_sam_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len, bool *is_numeric);

// THIS FOLLOWING CONSTANTS CANNOT CHANGE - THEY ARE PART OF THE FILE FORMAT
#define MAX_SAM_MD_LEN 1000 // maximum length of MD that is shortened.

// ---------------------------
// FASTA and FASTQ Stuff
// ---------------------------
extern void seg_fasta_initialize();
extern const char *seg_fastq_data_line();
extern const char *seg_fasta_data_line();

// ------------------
// GFF3 Stuff
// ------------------
extern const char *seg_gff3_data_line (VBlockP vb_, const char *field_start_line);
extern void seg_gff3_initialize (VBlockP vb_);

#define AOS_NUM_ENTRIES '\1' // first char of dictionary word in case of a valid AoS- followed by the number of entries. this value is part of the file format so cannot be (easily) changed.
extern void seg_gff3_array_of_struct_ctxs (VBlockGFF3P vb, DictIdType dict_id, unsigned num_items, 
                                           MtfContextP *ctx_array, MtfContextP *enst_ctx); // out

// ------------------
// ME23 Stuff
// ------------------
extern const char *seg_me23_data_line (VBlockP vb_, const char *field_start_line);

// ------------------
// Seg utilities
// ------------------

// create extendent field contexts in the correct order of the fields
#define EXTENDED_FIELD_CTX(extended_field, dict_id_num) { \
    MtfContext *ctx = mtf_get_ctx_by_dict_id (vb, (DictIdType)dict_id_num); \
    ASSERT (ctx->did_i == extended_field, "Error: expecting ctx->did_i=%u to be %u", ctx->did_i, extended_field); \
    dict_id_fields[ctx->did_i] = ctx->dict_id.num; \
}

#define SAFE_ASSIGN(reg,addr,char_val) /* we are careful to evaluate addr, char_val only once, less they contain eg ++ */ \
    char *__addr##reg = (char*)(addr); \
    char __save##reg = *__addr##reg; \
    *__addr##reg = (char_val);

#define SAFE_RESTORE(reg) *__addr##reg = __save##reg; 

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