// ------------------------------------------------------------------
//   reconstruct.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "writer.h"
#include "file.h"

typedef struct { char s[100]; } PizDisCoords; 
extern PizDisCoords piz_dis_coords (VBlockP vb); // for ASSPIZ

typedef struct { char s[100]; } PizDisQname; 
extern PizDisQname piz_dis_qname (VBlockP vb); // for ASSPIZ

extern void asspiz_text (VBlockP vb, FUNCLINE);

#define ASSPIZ(condition, format, ...) ({ \
    if (!(condition)) { \
        DO_ONCE { /* first thread to fail at this point prints error and starts exit flow, other threads stall */ \
            asspiz_text ((VBlockP)vb, __FUNCLINE);\
            fprintf (stderr, format "\n", __VA_ARGS__); \
            exit_on_error(true); \
        } \
        else stall(); \
    } \
})

#define ASSPIZ0(condition, string) ASSPIZ (condition, string "%s", "")

#define ABORT_PIZ(format, ...) ({ \
    DO_ONCE { /* first thread to fail at this point prints error and starts exit flow, other threads stall */ \
        asspiz_text ((VBlockP)vb, __FUNCLINE);\
        fprintf (stderr, format "\n", __VA_ARGS__); \
        exit_on_error(true); \
    } \
    else stall(); \
})

#define ABORT_PIZ0(string) ABORT_PIZ (string "%s", "")

// goes into ctx->history if not STORE_INT
typedef enum __attribute__ ((__packed__)) { LookupTxtData, LookupDict, LookupLocal, LookupPerLine } LookupType;
#define LOOKUP_TYPE_NAMES { "LookupTxtData", "LookupDict", "LookupLocal", "LookupPerLine" }
extern rom lookup_type_name (LookupType lookup);

typedef struct __attribute__ ((__packed__)) { // 9 bytes
    TxtWord;           // for TxtData, Local, PerLine - index into buffer. For Dict - word index (note: gcc/clang flag -fms-extensions is needed for this type of anonymous struct use)
    LookupType lookup    : 2;
    uint8_t ctx_specific : 6;
} HistoryWord;

extern int32_t reconstruct_from_ctx_do (VBlockP vb, Did did_i, char sep, ReconType reconstruct, rom func);
#define reconstruct_from_ctx(vb,did_i,sep,reconstruct) reconstruct_from_ctx_do ((VBlockP)(vb),(did_i),(sep),(reconstruct), __FUNCTION__)

extern void reconstruct_one_snip (VBlockP vb, ContextP ctx, WordIndex word_index, STRp(snip), ReconType reconstruct);
extern uint32_t reconstruct_from_local_sequence (VBlockP vb, ContextP ctx, uint32_t len, ReconType reconstruct);
extern int64_t reconstruct_from_local_int (VBlockP vb, ContextP ctx, char separator /* 0 if none */, ReconType reconstruct);
extern HasNewValue reconstruct_demultiplex (VBlockP vb, ContextP ctx, STRp(snip), int channel_i, ValueType *new_value, ReconType reconstruct);

extern ContextP reconstruct_get_other_ctx_from_snip (VBlockP vb, ContextP ctx, pSTRp (snip));

extern ContextP recon_multi_dict_id_get_ctx_first_time (VBlockP vb, ContextP ctx, STRp(snip), unsigned ctx_i);
#define MCTX(ctx_i,snip,snip_len) ((ctx->ctx_cache.len32 && *B(ContextP, ctx->ctx_cache, ctx_i)) \
                                        ? *B(ContextP, ctx->ctx_cache, ctx_i)                    \
                                        : recon_multi_dict_id_get_ctx_first_time ((VBlockP)vb, ctx, (snip), (snip_len), (ctx_i)))

// use SCTX if we are certain that ctx can only be one other_dict_id in its snips 
// snip is expected to be : 1-char-code + base64-dict_id + other stuff. snip is modified to be after the dict_id
#define SCTX(snip) ({ ContextP sctx;                                \
                      if (ctx->other_did_i != DID_NONE)  {          \
                          snip       += base64_sizeof (DictId) + 1; \
                          snip##_len -= base64_sizeof (DictId) + 1; \
                          sctx = CTX(ctx->other_did_i);             \
                      }                                             \
                      else                                          \
                          sctx = reconstruct_get_other_ctx_from_snip (VB, ctx, &snip, &snip##_len); \
                      sctx;                                         \
                   })

// a version of SCTX with where just a base64-dict_id
#define SCTX0(snip) ({ rom snip0 = (snip)-1;                        \
                       uint32_t snip0_len = base64_sizeof (DictId) + 1; \
                       SCTX(snip0); })

//--------------
// Peeking
//--------------
extern void recon_stack_initialize (void);
extern void recon_stack_push (VBlockP vb, ContextP ctx);
extern void recon_stack_pop (VBlockP vb, ContextP ctx, bool is_done_peek);

extern ValueType reconstruct_peek (VBlockP vb, ContextP ctx, pSTRp(txt));
extern ValueType reconstruct_peek_by_dict_id (VBlockP vb, DictId dict_id, pSTRp(txt));
#define reconstruct_peek_(vb, dict_id, txt, txt_len) reconstruct_peek_by_dict_id ((VBlockP)(vb), (dict_id), (txt), (txt_len))

extern int64_t reconstruct_peek_local_int (VBlockP vb, ContextP ctx, int offset);

//--------------
// history stuff
//--------------
extern void reconstruct_store_history (VBlockP vb, ContextP ctx);
extern void reconstruct_store_history_rollback_recon (VBlockP vb, ContextP ctx, rom recon_start);
extern void reconstruct_copy_dropped_line_history (VBlockP vb);
extern void reconstruct_to_history (VBlockP vb, ContextP ctx);

typedef bool (*PizReconstructSpecialInfoSubfields) (VBlockP vb, Did did_i, DictId dict_id);

// gets snip, snip_len from b250 data
#define LOAD_SNIP(did_i) ctx_get_next_snip (VB, CTX(did_i), false, &snip, &snip_len) 
#define PEEK_SNIP(did_i) ctx_peek_next_snip (VB, CTX(did_i), &snip, &snip_len)

#define LOAD_SNIP_FROM_LOCAL(ctx) ( {                                   \
    ASSPIZ (ctx->next_local < ctx->local.len32, "%s.local exhausted: next_local=%u len=%u%s", (ctx)->tag_name, (ctx)->next_local, (ctx)->local.len32, !(ctx)->local.len32 ? " (since len=0, perhaps it is not loaded? check IS_SKIP function)" : ""); \
    uint32_t start = ctx->next_local;                                   \
    ARRAY (char, data, ctx->local);                                     \
    uint32_t next_local = ctx->next_local; /* automatic for speed */    \
    uint32_t len = ctx->local.len32;                                    \
    while (next_local < len && data[next_local]) next_local++;          \
    snip = &data[start];                                                \
    snip_len = next_local - start;                                      \
    ctx->next_local = next_local + 1; /* skip the separator */          \
    start;                                                              \
} ) 

#define NEXT_ERRFMT "%s: not enough data in %s.local: next_local=%u + recon_len=%u > local.len=%u"

#define NEXTLOCAL(type, ctx) \
    ({ ASSPIZ (ctx->next_local < ctx->local.len32, "%s.local exhausted: next_local=%u len=%u", (ctx)->tag_name, (ctx)->next_local, (ctx)->local.len32); \
       *B(type, (ctx)->local, (ctx)->next_local++); })

#define PEEKNEXTLOCAL(type, ctx, offset) \
    ({ ASSPIZ ((ctx)->next_local + (offset) < (ctx)->local.len32, "PEEKNEXTLOCAL: %s.local exhausted: next_local=%u len=%u", (ctx)->tag_name, (ctx)->next_local, (ctx)->local.len32); \
       *B(type, (ctx)->local, (ctx)->next_local + (offset)); })

#ifdef DEBUG
#define RECONSTRUCT(s,s_len)                                                    \
    ({ uint32_t new_len = (uint32_t)(s_len); /* copy in case caller uses ++ */  \
       ASSPIZ (vb->txt_data.len + new_len <= vb->txt_data.size, "RECONSTRUCT: txt_data overflow: txt_data.len=%"PRIu64" str_len=%u vb->txt_data.size=%u. vb->txt_data dumped to %s.gz", \
               vb->txt_data.len, new_len, (uint32_t)vb->txt_data.size, txtfile_dump_vb ((VBlockP)(vb), z_name));  /* leave vb->txt_data.len 64b to detect bugs */     \
       memcpy (BAFTtxt, (s), new_len);                                          \
       Ltxt += new_len; })
#else
#define RECONSTRUCT(s,s_len) /* note: for speed, we don't check for overflow here. instead, we check after each line in container_reconstruct() */ \
    ({ uint32_t new_len = (uint32_t)(s_len); /* copy in case caller uses ++ */  \
       memcpy (BAFTtxt, (s), new_len);                                          \
       Ltxt += new_len; })
#endif

#define RECONSTRUCT1(c) ({ vb->txt_data.data[Ltxt++] = (c); }) // we don't use BNXTc bc it is too expensive
#define RECONSTRUCT_snip RECONSTRUCT (snip, snip_len)
#define RECONSTRUCT_SEP(s,len,sep) ({ RECONSTRUCT((s), (len)); RECONSTRUCT1 (sep); })
#define RECONSTRUCT_TABBED(s,len) RECONSTRUCT_SEP (s, len, '\t')
#define RECONSTRUCT_BUF(buf) RECONSTRUCT((buf).data,(buf).len32)
#define RECONSTRUCT_LAST_TXT(ctx) ({ ASSERT_LAST_TXT_VALID(ctx); RECONSTRUCT (last_txtx (vb, (ctx)), (ctx)->last_txt.len); })
#define RECONSTRUCT_NEXT(ctx,recon_len)                 \
    ({ ASSPIZ ((ctx)->next_local + (recon_len) <= (ctx)->local.len32, NEXT_ERRFMT, "RECONSTRUCT_NEXT", (ctx)->tag_name, (ctx)->next_local, (recon_len), (ctx)->local.len32); \
       rom next = Bc((ctx)->local, (ctx)->next_local);  \
       if (reconstruct) RECONSTRUCT(next, (recon_len)); \
       (ctx)->next_local += (recon_len);                \
       next; })

#define RECONSTRUCT_NEXT_REV(ctx,recon_len) /* reconstructs reverse string */       \
    ({ ASSPIZ ((ctx)->next_local + (recon_len) <= (ctx)->local.len32, NEXT_ERRFMT, "RECONSTRUCT_NEXT_REV", (ctx)->tag_name, (ctx)->next_local, (recon_len), (ctx)->local.len32); \
       if (reconstruct) {                                                           \
           str_reverse (BAFTtxt, Bc((ctx)->local, (ctx)->next_local), (recon_len)); \
           Ltxt += (recon_len);                                         \
       }                                                                            \
       (ctx)->next_local += (recon_len); })

#define RECONSTRUCT_INT(n)            buf_add_int_as_text (&vb->txt_data, (n)) /* note: don't add nul-terminator, as sometimes we expect reconstruction to be exact, and memory after could be already occupied - eg sam_piz_prim_add_QNAME */ 
#define RECONSTRUCT_HEX(n, uppercase) buf_add_hex_as_text (&vb->txt_data, (n), (uppercase)) 

#define RECONSTRUCT_FROM_DICT(did_i,add_tab)    \
    ({ WordIndex wi = LOAD_SNIP (did_i);        \
       RECONSTRUCT_snip;                        \
       if (add_tab) RECONSTRUCT1 ('\t');        \
       wi;                               })

// binary reconstructions
#define RECONSTRUCT_BIN8(n)  RECONSTRUCT (&n, 1)
#define RECONSTRUCT_BIN16(n) ({ uint16_t lten = (uint16_t)(n); lten = LTEN16(lten); RECONSTRUCT (&lten, sizeof (uint16_t)); })
#define RECONSTRUCT_BIN32(n) ({ uint32_t lten = (uint32_t)(n); lten = LTEN32(lten); RECONSTRUCT (&lten, sizeof (uint32_t)); })
