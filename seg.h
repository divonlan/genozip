// ------------------------------------------------------------------
//   seg.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SEGREGATE_INCLUDED
#define SEGREGATE_INCLUDED

#include <stdint.h>
#include "genozip.h"
#include "sections.h"

typedef const char *SegDataLineFuncType (VBlockP vb, const char *field_start_line);
typedef void SegInitializer (VBlockP vb);

extern void seg_all_data_lines (VBlockP vb, SegDataLineFuncType seg_data_line, SegInitializer seg_initialize, unsigned sizeof_line); 

extern DictIdType seg_vcf_get_format_subfield (const char **data, uint32_t *len);

extern const char *seg_get_next_item (void *vb, const char *str, int *str_len, bool allow_newline, bool allow_tab, bool allow_colon, 
                                      unsigned *len, char *separator, bool *has_13, // out
                                      const char *item_name);
extern const char *seg_get_next_line (void *vb_, const char *str, int *str_len, unsigned *len, bool *has_13 /* out */, const char *item_name);

extern void seg_store (VBlockP vb, 
                       bool *dst_is_spillover, uint32_t *dst_start, uint32_t *dst_len, // out
                       BufferP src_buf, uint32_t size, // Either src_buf OR size must be given
                       const char *limit_txt_data, // we cannot store in txt starting here. if NULL always allocates in txt_data_spillover
                       bool align32); // does start address need to be 32bit aligned to prevent aliasing issues

extern uint32_t seg_one_subfield (VBlockP vb, const char *str, unsigned len,
                                  DictIdType dict_id, SectionType sec_b250, int accounts_for_chars);

extern uint32_t seg_one_snip (VBlockP vb, const char *str, unsigned len, int did_i, SectionType sec_b250,
                              bool *is_new); // optional out

#define seg_one_field(vb,str,len,f) seg_one_snip ((VBlockP)(vb), (str), (len), (f), FIELD_TO_B250_SECTION((vb)->data_type, f), NULL)

extern int32_t seg_pos_field (VBlockP vb, int32_t last_pos, int32_t *last_pos_delta, bool allow_non_number, 
                              int pos_field, SectionType sec_pos_b250,
                              const char *pos_str, unsigned pos_len, const char *field_name);

extern void seg_id_field (VBlockP vb, BufferP id_buf, int id_field, char *id_snip, unsigned id_snip_len, bool extra_bit);

extern void seg_add_to_data_buf (VBlockP vb, BufferP buf, SectionType sec, 
                                 const char *snip, unsigned snip_len, char add_separator, unsigned add_bytes);

extern void seg_compound_field (VBlockP vb, MtfContextP field_ctx, const char *field, unsigned field_len, 
                                SubfieldMapperP mapper, DictIdType sf_dict_id, bool ws_is_sep, bool account_for_13,
                                SectionType field_b250_sec, SectionType sf_b250_sec);
                               
// ---------
// VCF Stuff
// ---------
extern SegDataLineFuncType seg_vcf_data_line;
extern SegInitializer seg_vcf_initialize;
extern void seg_vcf_complete_missing_lines (VBlockVCFP vb);

// ---------
// SAM Stuff
// ---------
extern SegDataLineFuncType seg_sam_data_line;
extern SegInitializer seg_sam_initialize;
extern uint32_t seg_sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len);
extern uint32_t seg_sam_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len, bool *is_numeric);

// THIS FOLLOWING CONSTANTS CANNOT CHANGE - THEY ARE PART OF THE FILE FORMAT
#define MAX_SAM_MD_LEN 1000 // maximum length of MD that is shortened.

// ---------------------------
// FASTA and FASTQ Stuff
// ---------------------------
extern SegDataLineFuncType seg_fastq_data_line;
extern SegDataLineFuncType seg_fasta_data_line;
extern SegInitializer seg_fasta_initialize;

// ------------------
// ME23 Stuff
// ------------------
extern SegDataLineFuncType seg_me23_data_line;
extern SegInitializer seg_me23_initialize;

// ------------------
// Seg utilities
// ------------------

//#define ASSSEG(condition, p_into_txt, format, ...) ASSERT(condition, format "ext %u", __VA_ARGS__, 6)

#define ASSSEG(condition, p_into_txt, format, ...) \
    ASSERT (condition, format "\nFile: %s vb_line_i:%u vb_i:%u pos_in_vb: %"PRIi64" pos_in_file: %"PRIi64\
                              "\nvb pos in file (0-based):%"PRIu64" - %"PRIu64" (length %"PRIu64")" \
                              "\n%d characters before to %d characters after (in quotes): \"%.*s\""\
                              "\nTo get vblock: %s %s | head -c %"PRIu64" | tail -c %"PRIu64 " > vb", \
            __VA_ARGS__, txt_name, vb->line_i, vb->vblock_i, \
            /* pos_in_vb:         */ p_into_txt ? (p_into_txt - vb->txt_data.data) : -1LL, \
            /* pos_in_file:       */ p_into_txt ? (vb->vb_position_txt_file + (p_into_txt - vb->txt_data.data)) : -1LL,\
            /* vb start pos file: */ vb->vb_position_txt_file, \
            /* vb end pos file:   */ vb->vb_position_txt_file + vb->txt_data.len-1, \
            /* vb length:         */ vb->txt_data.len,\
            /* chars before:      */ p_into_txt ? MIN (30, (unsigned)(p_into_txt - vb->txt_data.data)) : -1, \
            /* chars after:       */ p_into_txt ? MIN (30, (unsigned)(vb->txt_data.data + vb->txt_data.len - p_into_txt)) : -1,\
            /* snip len:          */ p_into_txt ? (unsigned)(MIN (p_into_txt+31, vb->txt_data.data + vb->txt_data.len) /* end pos */ - MAX (p_into_txt-30, vb->txt_data.data) /* start_pos */) : -1,\
            /* condition for snip */ (vb->txt_data.data && p_into_txt && (p_into_txt >= vb->txt_data.data) && (p_into_txt <= /* = too */ vb->txt_data.data + vb->txt_data.len) ? \
            /* snip start:        */    MAX (p_into_txt-30, vb->txt_data.data) : "(inaccessible)"),\
            /* head, tail params: */ file_viewer (txt_file), txt_name, vb->vb_position_txt_file + vb->txt_data.len, vb->txt_data.len)

#define ASSSEG0(condition, p_into_txt, err_str) ASSSEG (condition, p_into_txt, err_str "%s", "")

#define ABOSEG(p_into_txt, format, ...) ASSSEG(false, p_into_txt, format, __VA_ARGS__)

#define ABOSEG0(p_into_txt, err_str) ABOSEG(false, p_into_txt, format, err_str "%s", "")

#endif