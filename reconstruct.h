// ------------------------------------------------------------------
//   reconstruct.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// goes into ctx->history if not STORE_INT
typedef enum __attribute__ ((__packed__)) { LookupTxtData, LookupDict, LookupLocal, LookupPerLine } LookupType;

typedef struct __attribute__ ((__packed__)) { // 9 bytes
    TxtWord;           // for TxtData, Local, PerLine - index into buffer. For Dict - word index (note: gcc/clang flag -fms-extensions is needed for this type of anonymous struct use)
    LookupType lookup    : 2;
    uint8_t ctx_specific : 6;
} HistoryWord;

extern void reconstruct_initialize (void);

extern int32_t reconstruct_from_ctx_do (VBlockP vb, DidIType did_i, char sep, bool reconstruct, rom func);
#define reconstruct_from_ctx(vb,did_i,sep,reconstruct) reconstruct_from_ctx_do ((VBlockP)(vb),(did_i),(sep),(reconstruct), __FUNCTION__)

extern void reconstruct_one_snip (VBlockP vb, ContextP ctx, WordIndex word_index, STRp(snip), bool reconstruct);
extern void reconstruct_from_local_sequence (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct);
extern int64_t reconstruct_from_local_int (VBlockP vb, ContextP ctx, char separator /* 0 if none */, bool reconstruct);
extern HasNewValue reconstruct_demultiplex (VBlockP vb, ContextP ctx, STRp(snip), int channel_i, ValueType *new_value, bool reconstruct);

extern ContextP reconstruct_get_other_ctx_from_snip (VBlockP vb, ContextP ctx, pSTRp (snip));

extern ContextP recon_multi_dict_id_get_ctx_first_time (VBlockP vb, ContextP ctx, STRp(snip), unsigned ctx_i);
#define MCTX(ctx_i,snip,snip_len) ((ctx->con_cache.len32 && *B(ContextP, ctx->con_cache, ctx_i)) \
                                        ? *B(ContextP, ctx->con_cache, ctx_i)                    \
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
#define SCTX0(snip) ({ rom snip0 = (snip)-1;                    \
                       uint32_t snip0_len = base64_sizeof (DictId) + 1; \
                       SCTX(snip0); })

//--------------
// Peeking
//--------------
extern ValueType reconstruct_peek (VBlockP vb, ContextP ctx, pSTRp(txt));
extern ValueType reconstruct_peek_do (VBlockP vb, DictId dict_id, pSTRp(txt));
#define reconstruct_peek_(vb, dict_id, txt, txt_len) reconstruct_peek_do ((VBlockP)(vb), (DictId)(dict_id), (txt), (txt_len))
extern int64_t reconstruct_peek_local_int (VBlockP vb, ContextP ctx, int offset);
extern uint32_t reconstruct_peek_repeats (VBlockP vb, ContextP ctx, char repsep);

//--------------
// history stuff
//--------------
extern void reconstruct_store_history (VBlockP vb, ContextP ctx);
extern void reconstruct_copy_dropped_line_history (VBlockP vb);
extern void reconstruct_set_buddy (VBlockP vb);
extern HasNewValue reconstruct_from_buddy (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct, ValueType *new_value);
extern void reconstruct_from_buddy_get_textual_snip (VBlockP vb, ContextP ctx, pSTRp(snip));


typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, DidIType did_i, DictId dict_id);

// gets snip, snip_len from b250 data
#define LOAD_SNIP(did_i) ctx_get_next_snip (VB, CTX(did_i), false, &snip, &snip_len) 
#define PEEK_SNIP(did_i) ctx_peek_next_snip (VB, CTX(did_i), &snip, &snip_len)

#define LOAD_SNIP_FROM_LOCAL(ctx) ( {           \
    uint32_t start = ctx->next_local;           \
    ARRAY (char, data, ctx->local);             \
    uint32_t next_local = ctx->next_local; /* automatic for speed */ \
    uint32_t len = ctx->local.len32;            \
    while (next_local < len && data[next_local]) next_local++; \
    snip = &data[start];                        \
    snip_len = next_local - start;              \
    ctx->next_local = next_local + 1; /* skip the separator */ \
    start;                                      \
} ) 

#define NEXTLOCAL(type, ctx) \
    ({ ASSPIZ ((ctx)->next_local < (ctx)->local.len32, "NEXTLOCAL: not enough data in %s.local", (ctx)->tag_name); \
       *B(type, (ctx)->local, (ctx)->next_local++); })

#define PEEKNEXTLOCAL(type, ctx, offset) \
    ({ ASSPIZ ((ctx)->next_local + (offset) < (ctx)->local.len32, "PEEKNEXTLOCAL: not enough data in %s.local", (ctx)->tag_name); \
       *B(type, (ctx)->local, (ctx)->next_local + (offset)); })

#define RECONSTRUCT(s,len) buf_add (&vb->txt_data, (char*)(s), (len))
#define RECONSTRUCT1(c) BNXTc (vb->txt_data) = c
#define RECONSTRUCT_snip RECONSTRUCT (snip, snip_len)
#define RECONSTRUCT_SEP(s,len,sep) ({ RECONSTRUCT((s), (len)); RECONSTRUCT1 (sep); })
#define RECONSTRUCT_TABBED(s,len) RECONSTRUCT_SEP (s, len, '\t')
#define RECONSTRUCT_BUF(buf) RECONSTRUCT((buf).data,(buf).len)
#define RECONSTRUCT_NEXT(ctx,recon_len) \
    ({ ASSPIZ ((ctx)->next_local + (recon_len) <= (ctx)->local.len, "RECONSTRUCT_NEXT: not enough data in %s.local", (ctx)->tag_name); \
       if (reconstruct) RECONSTRUCT(Bc((ctx)->local, (ctx)->next_local), (recon_len)); \
       (ctx)->next_local += (recon_len); })

#define RECONSTRUCT_NEXT_REV(ctx,recon_len) /* reconstructs reverse string */\
    ({ ASSPIZ ((ctx)->next_local + (recon_len) <= (ctx)->local.len, "RECONSTRUCT_NEXT_REV: not enough data in %s.local", (ctx)->tag_name); \
       if (reconstruct) { \
           str_reverse (BAFTtxt, Bc((ctx)->local, (ctx)->next_local), (recon_len)); \
           vb->txt_data.len += (recon_len); \
       } \
       (ctx)->next_local += (recon_len); })

#define RECONSTRUCT_INT(n) buf_add_int_as_text (&vb->txt_data, (n)) /* note: don't add nul-terminator, as sometimes we expect reconstruction to be exact, and memory after could be already occupied - eg sam_piz_prim_add_QNAME */ 
#define RECONSTRUCT_HEX(n, uppercase) buf_add_hex_as_text (&vb->txt_data, (n), (uppercase)) 

#define RECONSTRUCT_FROM_DICT(did_i,add_tab)    \
    ({ WordIndex wi = LOAD_SNIP (did_i);        \
       RECONSTRUCT_snip;            \
       if (add_tab) RECONSTRUCT1 ('\t');        \
       wi;                               })

// binary reconstructions
#define RECONSTRUCT_BIN8(n)  RECONSTRUCT (&n, 1)
#define RECONSTRUCT_BIN16(n) ({ uint16_t lten = (uint16_t)(n); lten = LTEN16(lten); RECONSTRUCT (&lten, sizeof (uint16_t)); })
#define RECONSTRUCT_BIN32(n) ({ uint32_t lten = (uint32_t)(n); lten = LTEN32(lten); RECONSTRUCT (&lten, sizeof (uint32_t)); })
