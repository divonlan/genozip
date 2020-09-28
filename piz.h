// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"

extern bool piz_dispatcher (bool is_first_component, bool is_last_file);

#define piz_is_skip_section(vb,st,dict_id) (vb->data_type != DT_NONE     && DTP(is_skip_secetion)  && DTP (is_skip_secetion)((VBlockP)(vb), (st), (dict_id)))
#define piz_is_skip_sectionz(st,dict_id)   (z_file->data_type != DT_NONE && DTPZ(is_skip_secetion) && DTPZ(is_skip_secetion)(NULL, (st), (dict_id)))

// --------------------------------------------------
// utilities for use by piz_*_uncompress_all_sections
// --------------------------------------------------

extern uint32_t piz_uncompress_all_ctxs (VBlockP vb, uint32_t pair_vb_i);

// ----------------------------------------------
// utilities for use by piz_*_reconstruct_vb
// ----------------------------------------------

extern int32_t piz_reconstruct_from_ctx_do (VBlockP vb, DidIType did_i, char sep);
#define piz_reconstruct_from_ctx(vb,did_i,sep) piz_reconstruct_from_ctx_do ((VBlockP)(vb),(did_i),(sep))

extern void piz_reconstruct_one_snip (VBlockP vb, ContextP ctx, WordIndex word_index, const char *snip, unsigned snip_len);

typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, DidIType did_i, DictId dict_id);

extern void piz_reconstruct_structured_do (VBlockP vb, DictId dict_id, ConstStructuredP st, const char *prefixes, uint32_t prefixes_len);

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

#endif

