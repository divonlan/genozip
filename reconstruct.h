// ------------------------------------------------------------------
//   reconstruct.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// goes into ctx->history if not STORE_INT
typedef struct {
    CharIndex char_index;
    uint32_t snip_len;
    enum { LookupTxtData, LookupDict, LookupLocal, LookupPerLine } lookup;
} HistoryWord;

extern void reconstruct_initialize (void);

extern int32_t reconstruct_from_ctx_do (VBlockP vb, DidIType did_i, char sep, bool reconstruct, const char *func);
#define reconstruct_from_ctx(vb,did_i,sep,reconstruct) reconstruct_from_ctx_do ((VBlockP)(vb),(did_i),(sep),(reconstruct), __FUNCTION__)

extern void reconstruct_one_snip (VBlockP vb, ContextP ctx, WordIndex word_index, STRp(snip), bool reconstruct);
extern void reconstruct_from_local_sequence (VBlockP vb, ContextP ctx, STRp(snip));

extern ContextP reconstruct_get_other_ctx_from_snip (VBlockP vb, ContextP ctx, pSTRp (snip));

extern ContextP recon_multi_dict_id_get_ctx_first_time (VBlockP vb, ContextP ctx, STRp(snip), unsigned ctx_i);
#define MCTX(ctx_i,snip,snip_len) ((ctx->con_cache.len && *ENT(ContextP, ctx->con_cache, ctx_i)) \
                                        ? *ENT(ContextP, ctx->con_cache, ctx_i)                  \
                                        : recon_multi_dict_id_get_ctx_first_time ((VBlockP)vb, ctx, (snip), (snip_len), (ctx_i)))

// use SCTX if we are certain that ctx can only be one other_dict_id in its snips 
// snip is expected to be : 1-char-code + base64-dict_id + other stuff. snip is modified to be after the dict_id
#define SCTX(snip) ({ ContextP sctx;                                \
                      if (ctx->other_did_i != DID_I_NONE)  {        \
                          snip       += base64_sizeof (DictId) + 1; \
                          snip##_len -= base64_sizeof (DictId) + 1; \
                          sctx = CTX(ctx->other_did_i);             \
                      }                                             \
                      else                                          \
                          sctx = reconstruct_get_other_ctx_from_snip (VB, ctx, &snip, &snip##_len); \
                      sctx;                                         \
                   })

// a version of SCTX with where just a base64-dict_id
#define SCTX0(snip) ({ const char *snip0 = (snip)-1;                    \
                       uint32_t snip0_len = base64_sizeof (DictId) + 1; \
                       SCTX(snip0); })

extern ValueType reconstruct_peek (VBlockP vb, ContextP ctx, pSTRp(txt));
extern ValueType reconstruct_peek_do (VBlockP vb, DictId dict_id, pSTRp(txt));
#define reconstruct_peek_(vb, dict_id, txt, txt_len) reconstruct_peek_do ((VBlockP)(vb), (DictId)(dict_id), (txt), (txt_len))

extern void reconstruct_set_buddy (VBlockP vb);
extern bool reconstruct_from_buddy (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct, ValueType *new_value);
extern void reconstruct_from_buddy_get_textual_snip (VBlockP vb, ContextP ctx, pSTRp(snip));

typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, DidIType did_i, DictId dict_id);

// gets snip, snip_len from b250 data
#define LOAD_SNIP(did_i) ctx_get_next_snip (VB, CTX(did_i), CTX(did_i)->flags.all_the_same, false, &snip, &snip_len); 
#define PEEK_SNIP(did_i) ctx_peek_next_snip (VB, CTX(did_i), CTX(did_i)->flags.all_the_same, &snip, &snip_len); 

#define LOAD_SNIP_FROM_LOCAL(ctx) ( {           \
    uint32_t start = ctx->next_local;           \
    ARRAY (char, data, ctx->local);             \
    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != 0) ctx->next_local++; \
    snip = &data[start];                        \
    snip_len = ctx->next_local - start;         \
    ctx->next_local++; /* skip the separator */ \
    start;                                      \
} ) 

#define RECONSTRUCT(s,len) buf_add (&vb->txt_data, (char*)(s), (len))
#define RECONSTRUCT1(c) NEXTENT (char, vb->txt_data) = c
#define RECONSTRUCT_SEP(s,len,sep) do { RECONSTRUCT((s), (len)); RECONSTRUCT1 (sep); } while(0)
#define RECONSTRUCT_TABBED(s,len) RECONSTRUCT_SEP (s, len, '\t')
#define RECONSTRUCT_BUF(buf) RECONSTRUCT((buf).data,(buf).len)

#define RECONSTRUCT_INT(n) ({ unsigned n_len = str_int ((n), AFTERENT (char, vb->txt_data)); \
                              vb->txt_data.len += n_len; n_len; })

#define RECONSTRUCT_FROM_DICT(did_i,add_tab) /* not a block so caller get the return value of ctx_get_next_snip */ \
    LOAD_SNIP (did_i);\
    RECONSTRUCT (snip, snip_len);\
    if (add_tab) RECONSTRUCT1 ('\t');

// binary reconstructions
#define RECONSTRUCT_BIN8(n)  RECONSTRUCT (&n, 1)
#define RECONSTRUCT_BIN16(n) do { uint16_t lten = (uint16_t)(n); lten = LTEN16(lten); RECONSTRUCT (&lten, sizeof (uint16_t)); } while (0)
#define RECONSTRUCT_BIN32(n) do { uint32_t lten = (uint32_t)(n); lten = LTEN32(lten); RECONSTRUCT (&lten, sizeof (uint32_t)); } while (0)
