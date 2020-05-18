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
extern void piz_map_iname_subfields (VBlockP vb);

// ----------------------
// VCF stuff
// ----------------------
extern bool piz_vcf_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern void piz_vcf_uncompress_vb(); // no parameter - implicit casting of VBlockP to VBlockVCFP
extern void v2v3_piz_vcf_map_iname_subfields (BufferP vb);

// ----------------------
// SAM stuff
// ----------------------
extern void piz_sam_reconstruct_vb ();

// ----------------------
// FASTQ + FASTA stuff
// ----------------------
extern bool piz_fast_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern bool piz_fast_test_grep (VBlockFASTP vb);
extern void piz_fasta_reconstruct_vb(); // no parameter - implicit casting of VBlockP to VBlockFASTP
extern void piz_fastq_reconstruct_vb();

// ----------------------
// GFF3 stuff
// ----------------------
extern bool piz_gff3_read_one_vb (VBlockP vb, SectionListEntryP sl);
extern void piz_gff3_reconstruct_vb(); // no parameter - implicit casting of VBlockP to VBlockGFF3P

// ----------------------
// 23andMe stuff
// ----------------------
extern void piz_me23_reconstruct_vb (VBlockP vb);

// ----------------------------------------------
// utilities for use by piz_*_read_one_vb
// ----------------------------------------------

#define PREPARE_TO_READ(vbblock_type,max_sections,sec_type_vb_header)  \
    START_TIMER; \
    vbblock_type *vb = (vbblock_type *)vb_; \
    SectionListEntryP sl = sections_vb_first (vb->vblock_i); \
    int vb_header_offset = zfile_read_section (vb_, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data", \
                                               sizeof(sec_type_vb_header), SEC_VB_HEADER, sl++); \
    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);\
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ \
    buf_alloc (vb, &vb->z_section_headers, (max_sections) * sizeof(char*), 0, "z_section_headers", 1); /* room for section headers */ \
    *FIRSTENT (unsigned, vb->z_section_headers) = vb_header_offset; /* vb header is at index 0 */ \
    vb->z_section_headers.len =1;

#define READ_SB_SECTION(sec,header_type,sb_i) \
    { NEXTENT (unsigned, vb->z_section_headers) = (uint32_t)vb->z_data.len; \
      zfile_read_section ((VBlockP)vb, vb->vblock_i, sb_i, &vb->z_data, "z_data", sizeof(header_type), sec, sl++); }

#define READ_SECTION(sec,header_type,is_optional) {  \
    if (!(is_optional) || sl->section_type == (sec)) \
        READ_SB_SECTION(sec, header_type, NO_SB_I);  \
}
#define READ_DATA_SECTION(sec,is_optional) READ_SECTION((sec), SectionHeader, (is_optional))

#define READ_DONE \
    COPY_TIMER (vb->profile.piz_read_one_vb); \
    vb->ready_to_dispatch = true; /* all good */ 

// --------------------------------------------------
// utilities for use by piz_*_uncompress_all_sections
// --------------------------------------------------

extern uint32_t piz_uncompress_all_ctxs (VBlockP vb);

extern void piz_map_compound_field (VBlockP vb, bool (*predicate)(DictIdType), SubfieldMapperP mapper);

// ----------------------------------------------
// utilities for use by piz_*_reconstruct_vb
// ----------------------------------------------

extern uint32_t piz_reconstruct_from_ctx (VBlockP vb, uint8_t did_i, const char *sep, unsigned sep_len, uint32_t txt_line_i);

extern void piz_reconstruct_compound_field (VBlockP vb, SubfieldMapperP mapper, const char *separator, unsigned separator_len, 
                                            const char *template, unsigned template_len, uint32_t txt_line_i);

extern void piz_reconstruct_seq_qual (VBlockP vb, MtfContext *ctx, uint32_t seq_len, uint32_t txt_line_i, bool grepped_out);

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

#define RECONSTRUCT_FROM_DICT(did_i,add_tab) piz_reconstruct_from_ctx ((VBlockP)vb, (did_i), (add_tab) ? "\t" : NULL, add_tab, txt_line_i)

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
    while (next < buf.len && data[next] != LOCAL_BUF_TEXT_SEP) next++;\
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

#endif

