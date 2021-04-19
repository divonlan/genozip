// ------------------------------------------------------------------
//   reconstruct.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef RECONSTRUCT_INCLUDED
#define RECONSTRUCT_INCLUDED

#include "genozip.h"

extern int32_t reconstruct_from_ctx_do (VBlockP vb, DidIType did_i, char sep, bool reconstruct, const char *func);
#define reconstruct_from_ctx(vb,did_i,sep,reconstruct) reconstruct_from_ctx_do ((VBlockP)(vb),(did_i),(sep),(reconstruct), __FUNCTION__)

extern void reconstruct_one_snip (VBlockP vb, ContextP ctx, WordIndex word_index, const char *snip, unsigned snip_len, bool reconstruct);

extern int64_t reconstruct_from_local_int (VBlockP vb, ContextP ctx, char seperator /* 0 if none */, bool reconstruct);

typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, DidIType did_i, DictId dict_id);

// gets snip, snip_len from b250 data
#define LOAD_SNIP(did_i) ctx_get_next_snip ((VBlockP)vb, &vb->contexts[(did_i)], vb->contexts[(did_i)].flags.all_the_same, false, &snip, &snip_len); 

#define RECONSTRUCT(s,len) buf_add (&vb->txt_data, (char*)(s), (len))
#define RECONSTRUCT1(c) NEXTENT (char, vb->txt_data) = c
#define RECONSTRUCT_SEP(s,len,sep) do { RECONSTRUCT((s), (len)); RECONSTRUCT1 (sep); } while(0)
#define RECONSTRUCT_TABBED(s,len) RECONSTRUCT_SEP (s, len, '\t')

#define RECONSTRUCT_INT(n) unsigned n_len = str_int ((n), AFTERENT (char, vb->txt_data)); /* not in a block because some need access to n_len */ \
                           vb->txt_data.len += n_len; 

#define RECONSTRUCT_FROM_DICT(did_i,add_tab) /* not a block so caller get the return value of ctx_get_next_snip */ \
    LOAD_SNIP (did_i);\
    RECONSTRUCT (snip, snip_len);\
    if (add_tab) RECONSTRUCT1 ('\t');

// binary reconstructions
#define RECONSTRUCT_BIN8(n)  RECONSTRUCT (&n, 1)
#define RECONSTRUCT_BIN16(n) do { uint16_t lten = (uint16_t)(n); lten = LTEN16(lten); RECONSTRUCT (&lten, sizeof (uint16_t)); } while (0)
#define RECONSTRUCT_BIN32(n) do { uint32_t lten = (uint32_t)(n); lten = LTEN32(lten); RECONSTRUCT (&lten, sizeof (uint32_t)); } while (0)

#endif