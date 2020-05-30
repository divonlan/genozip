// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"

extern bool piz_dispatcher (const char *z_basename, unsigned max_threads, bool is_first_vcf_component, bool is_last_file);

extern int32_t piz_decode_pos (VBlockP vb, int32_t last_pos, const char *delta_snip, unsigned delta_snip_len, 
                               int32_t *last_delta, char *pos_str, unsigned *pos_len);

#define piz_is_skip_section(vb,st,dict_id) (vb->data_type != DT_NONE     && DTP(is_skip_secetion)  && DTP (is_skip_secetion)((VBlockP)(vb), (st), (dict_id)))
#define piz_is_skip_sectionz(st,dict_id)   (z_file->data_type != DT_NONE && DTPZ(is_skip_secetion) && DTPZ(is_skip_secetion)(NULL, (st), (dict_id)))

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

extern uint32_t piz_reconstruct_from_ctx_do (VBlockP vb, uint8_t did_i, char sep);
#define piz_reconstruct_from_ctx(vb,did_i,sep) piz_reconstruct_from_ctx_do ((VBlockP)(vb),(did_i),(sep))

extern void piz_reconstruct_one_snip (VBlockP vb, MtfContextP ctx, const char *snip, unsigned snip_len);

typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, uint8_t did_i, DictIdType dict_id);

extern void piz_reconstruct_structured_do (VBlockP vb, ConstStructuredP st, const char *prefixes, uint32_t prefixes_len);

typedef struct PizSubfieldMapper {
    uint8_t num_subfields;        
    DictIdType dict_id[MAX_SUBFIELDS];
} PizSubfieldMapper;

#define DECLARE_SNIP const char *snip=NULL; uint32_t snip_len=0

// gets snip, snip_len from b250 data
#define LOAD_SNIP(did_i) mtf_get_next_snip ((VBlockP)vb, &vb->contexts[(did_i)], NULL, &snip, &snip_len); 

#define RECONSTRUCT(s,len) buf_add (&vb->txt_data, (s), (len))
#define RECONSTRUCT1(c) NEXTENT (char, vb->txt_data) = c
#define RECONSTRUCT_SEP(s,len,sep) { RECONSTRUCT((s), (len)); RECONSTRUCT1 (sep); }
#define RECONSTRUCT_TABBED(s,len) RECONSTRUCT_SEP (s, len, '\t')

#define RECONSTRUCT_INT(n) unsigned n_len = str_int ((n), AFTERENT (char, vb->txt_data)); /* not in a block because some need access to n_len */ \
                           vb->txt_data.len += n_len; 

#define RECONSTRUCT_FROM_DICT(did_i,add_tab) /* not a block so caller get the return value of mtf_get_next_snip */ \
    LOAD_SNIP (did_i);\
    RECONSTRUCT (snip, snip_len);\
    if (add_tab) RECONSTRUCT1 ('\t');

#define RECONSTRUCT_FROM_DICT_POS(did_i,last_pos,update_last_pos,last_delta,add_tab) { \
    if ((did_i) != DID_I_NONE) LOAD_SNIP(did_i);\
    char pos_str[30];\
    uint32_t new_pos = piz_decode_pos ((VBlockP)vb, last_pos, snip, snip_len, last_delta, pos_str, &snip_len); \
    if (update_last_pos) last_pos = new_pos;\
    RECONSTRUCT (pos_str, snip_len);\
    if (add_tab) RECONSTRUCT1 ('\t'); }

#endif

