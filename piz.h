// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"

extern bool piz_dispatcher (const char *z_basename, unsigned max_threads, bool is_first_vcf_component, bool is_last_file);

extern int32_t piz_decode_pos (VBlockP vb, uint32_t txt_line_i,
                               int32_t last_pos, const char *delta_snip, unsigned delta_snip_len, 
                               int32_t *last_delta, char *pos_str, unsigned *pos_len);
extern void piz_map_iname_subfields (void);

// ----------------------
// VCF stuff
// ----------------------
extern void piz_vcf_read_one_vb  (VBlockP vb);
extern void piz_vcf_uncompress_one_vb (VBlockP vb);
extern void v2v3_piz_vcf_map_iname_subfields (BufferP vb);

// ----------------------
// SAM stuff
// ----------------------
extern void piz_sam_read_one_vb  (VBlockP vb);
extern void piz_sam_uncompress_one_vb (VBlockP vb);

// ----------------------
// FASTQ + FASTA stuff
// ----------------------
extern void piz_fast_read_one_vb (VBlockP vb);
extern void piz_fast_uncompress_one_vb (VBlockP vb);
extern bool piz_fast_test_grep (VBlockFASTP vb);

// ----------------------
// GFF3 stuff
// ----------------------
extern void piz_gff3_read_one_vb (VBlockP vb);
extern void piz_gff3_uncompress_one_vb
 (VBlockP vb);

// ----------------------
// 23andMe stuff
// ----------------------
extern void piz_me23_read_one_vb (VBlockP vb);
extern void piz_me23_uncompress_one_vb (VBlockP vb);

// ----------------------------------------------
// utilities for use by piz_*_read_one_vb
// ----------------------------------------------

#define PREPARE_TO_READ(vbblock_type,max_sections,sec_type_vb_header)  \
    START_TIMER; \
    vbblock_type *vb = (vbblock_type *)vb_; \
    SectionListEntry *sl = sections_vb_first (vb->vblock_i); \
    int vb_header_offset = zfile_read_section (vb_, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data", \
                                               sizeof(sec_type_vb_header), SEC_VB_HEADER, sl++); \
    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);\
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ \
    buf_alloc (vb, &vb->z_section_headers, (max_sections) * sizeof(char*), 0, "z_section_headers", 1); /* room for section headers */ \
    *FIRSTENT (unsigned, vb->z_section_headers) = vb_header_offset; /* vb header is at index 0 */ \
    vb->z_section_headers.len =1;

#define READ_SB_SECTION(sec,header_type,sb_i) \
    { *ENT (unsigned, vb->z_section_headers, vb->z_section_headers.len++) = vb->z_data.len; \
      zfile_read_section ((VBlockP)vb, vb->vblock_i, sb_i, &vb->z_data, "z_data", sizeof(header_type), sec, sl++); }

#define READ_SECTION(sec,header_type,is_optional) {  \
    if (!(is_optional) || sl->section_type == (sec)) \
        READ_SB_SECTION(sec, header_type, NO_SB_I);  \
}
#define READ_DATA_SECTION(sec,is_optional) READ_SECTION((sec), SectionHeader, (is_optional))

#define READ_FIELDS { for (int f=0; f < DTFZ(num_fields); f++) \
                          READ_SECTION (FIELD_TO_B250_SECTION(z_file->data_type, f), SectionHeaderBase250, false); }

#define READ_SUBFIELDS(num,sec_b250) { \
    num = sections_count_sec_type (vb->vblock_i, (sec_b250)); \
    for (uint8_t sf_i=0; sf_i < (num); sf_i++) \
        READ_SECTION (sec_b250, SectionHeaderBase250, false);  \
}

#define READ_DONE \
    COPY_TIMER (vb->profile.piz_read_one_vb); \
    vb->ready_to_dispatch = true; /* all good */ 

// --------------------------------------------------
// utilities for use by piz_*_uncompress_all_sections
// --------------------------------------------------

extern void piz_uncompress_compound_field (VBlockP vb, SectionType field_b250_sec, SectionType sf_b250_sec, SubfieldMapperP mapper, unsigned *section_i);

#define UNCOMPRESS_HEADER_AND_FIELDS(vb_block_type,also_uncompress_fields) \
    START_TIMER;\
    vb_block_type *vb = (vb_block_type *)vb_; \
    ARRAY (const unsigned, section_index, vb->z_section_headers); \
    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);\
    vb->first_line       = BGEN32 (header->first_line);      \
    vb->lines.len        = BGEN32 (header->num_lines);       \
    vb->vb_data_size     = BGEN32 (header->vb_data_size);    \
    vb->longest_line_len = BGEN32 (header->longest_line_len);\
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i); /* in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher because the dispatcher is re-initialized for every component */ \
    unsigned section_i=1;\
    if (also_uncompress_fields) UNCOMPRESS_FIELDS;

#define UNCOMPRESS_FIELDS { \
    for (int f=0 ; f < DTF(num_fields) ; f++) { \
        SectionType b250_sec = FIELD_TO_B250_SECTION(vb->data_type, f); \
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]); \
        zfile_uncompress_section ((VBlockP)vb, header, &vb->mtf_ctx[f].b250, "mtf_ctx.b250", b250_sec);\
    }\
}

#define UNCOMPRESS_SUBFIELDS(num_subfields,sec_b250) {\
    for (uint8_t sf_i=0; sf_i < num_subfields ; sf_i++) {\
        SectionHeaderBase250 *header = (SectionHeaderBase250 *)(vb->z_data.data + section_index[section_i++]); \
        if (zfile_is_skip_section (vb, (sec_b250), header->dict_id)) continue; \
        MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, header->dict_id, (sec_b250)-1); \
        zfile_uncompress_section ((VBlockP)vb, header, &ctx->b250, "mtf_ctx.b250", (sec_b250)); \
    } \
}

#define UNCOMPRESS_DATA_SECTION(sec, vb_buf_name, type, is_optional) { \
    SectionHeader *data_header = (SectionHeader *)(vb->z_data.data + section_index[section_i]); \
    if (!(is_optional) || (section_i < vb->z_section_headers.len && data_header->section_type == (sec))) { \
        zfile_uncompress_section ((VBlockP)vb, data_header, &vb->vb_buf_name, #vb_buf_name, (sec)); \
        vb->vb_buf_name.len /= sizeof(type); \
        section_i++;\
    }\
}

#define UNCOMPRESS_DONE \
    vb->is_processed = true; /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ \
    COPY_TIMER (vb->profile.compute);


// ----------------------------------------------
// utilities for use by piz_*_reconstruct_vb
// ----------------------------------------------

extern void piz_reconstruct_compound_field (VBlockP vb, SubfieldMapperP mapper, const char *separator, unsigned separator_len, 
                                            const char *template, unsigned template_len, uint32_t txt_line_i);

extern void piz_reconstruct_seq_qual (VBlockP vb, uint32_t seq_len, ConstBufferP data, uint32_t *next, 
                                      SectionType sec, uint32_t txt_line_i, bool grepped_out);

extern void piz_reconstruct_id (VBlockP vb, BufferP id_buf, uint32_t *next_id, 
                                const char *id_snip, unsigned id_snip_len, bool *extra_bit, bool add_tab);

typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, uint8_t did_i, DictIdType dict_id, uint32_t txt_line_i);

extern void piz_reconstruct_info (VBlockP vb, uint32_t iname_word_index, 
                                  const char *iname_snip, unsigned iname_snip_len, 
                                  PizReconstructSpecialInfoSubfields reconstruct_special_info_subfields,
                                  uint32_t txt_line_i, bool *has_13);

typedef struct PizSubfieldMapper {
    uint8_t num_subfields;        
    uint8_t did_i[MAX_SUBFIELDS]; // NIL for a subfield that has no values (no ctx) but still has a dict_id
    DictIdType dict_id[MAX_SUBFIELDS];
} PizSubfieldMapper;

#define DECLARE_SNIP const char *snip=NULL; uint32_t snip_len=0

// gets snip, snip_len from b250 data
#define LOAD_SNIP(did_i) mtf_get_next_snip ((VBlockP)vb, &vb->mtf_ctx[(did_i)], NULL, &snip, &snip_len, txt_line_i); 

#define RECONSTRUCT(s,len) buf_add (&vb->txt_data, (s), (len))
#define RECONSTRUCT1(s) buf_add (&vb->txt_data, (s), 1)
#define RECONSTRUCT_TABBED(s,len) { RECONSTRUCT((s), (len)); RECONSTRUCT ("\t", 1); }

#define RECONSTRUCT_FROM_DICT(did_i,add_tab)  /* we don't put in a {} so that caller can use index=RECONSTRUCT_FROM_DICT() */ \
    LOAD_SNIP(did_i) \
    RECONSTRUCT (snip, snip_len); \
    if (add_tab) RECONSTRUCT ("\t", 1); 

#define RECONSTRUCT_FROM_DICT_POS(did_i,last_pos,update_last_pos,last_delta,add_tab) { \
    if ((did_i) != DID_I_NONE) LOAD_SNIP(did_i);\
    char pos_str[30];\
    uint32_t new_pos = piz_decode_pos ((VBlockP)vb, txt_line_i, last_pos, snip, snip_len, last_delta, pos_str, &snip_len); \
    if (update_last_pos) last_pos = new_pos;\
    RECONSTRUCT (pos_str, snip_len);\
    if (add_tab) RECONSTRUCT ("\t", 1); }

#define LOAD_SNIP_FROM_BUF(buf,next,field_name) { \
    uint32_t start = next; \
    ARRAY (char, data, buf);\
    while (next < buf.len && data[next] != '\n') next++;\
    ASSERT (next < buf.len, \
            "Error reconstructing txt_line=%u: unexpected end of %s data (len=%u)", txt_line_i, field_name, (uint32_t)buf.len); \
    snip = &data[start];\
    snip_len = next - start; \
    next++; /* skip the tab */ }

// reconstructs from the buffer up to a tab    
#define RECONSTRUCT_FROM_BUF(buf,next,field_name,reconst_sep_str,reconst_sep_str_len) { \
    DECLARE_SNIP;\
    LOAD_SNIP_FROM_BUF(buf,next,field_name) \
    RECONSTRUCT (snip, snip_len); \
    if (reconst_sep_str_len) RECONSTRUCT (reconst_sep_str, reconst_sep_str_len);  }

// reconstructs a fix number of characters from a tab-less buffer
#define RECONSTRUCT_FROM_TABLESS_BUF(buf,next,fixed_len,add_tab,field_name) { \
    ARRAY (char, data, buf);\
    ASSERT ((next) + (fixed_len) <= buf.len, \
            "Error reconstructing txt_line=%u: unexpected end of " field_name " data (buf.len=%u next=%u fixed_len=%u)", \
            txt_line_i, (uint32_t)buf.len, (next), (fixed_len)); \
    RECONSTRUCT (&data[(next)], (fixed_len)); \
    if (add_tab) RECONSTRUCT ("\t", 1);  \
    (next) += (fixed_len); }

#define RECONSTRUCT_ID(did_i,id_data,next_id,extra_bit,add_tab) { \
    DECLARE_SNIP; \
    LOAD_SNIP (did_i);         \
    piz_reconstruct_id ((VBlockP)vb, id_data, next_id, snip, snip_len, extra_bit,add_tab);\
}

#endif

