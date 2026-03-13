// ------------------------------------------------------------------
//   seg_resizeable.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "buffer.h"
#include "file.h"
#include "local_type.h"

static LocalType lt_order[] = { 0/*uninitialzed*/, LT_UINT8, LT_INT8, LT_UINT16, LT_INT16, LT_UINT32, LT_INT32, LT_INT64 }; // seniority levels
static LocalType lt_order_hex[] = { -1,            LT_hex8,  -1,      LT_hex16,  -1,       LT_hex32,  -1,       LT_hex64 }; // LT_DYN_INT_h
static LocalType lt_order_HEX[] = { -1,            LT_HEX8,  -1,      LT_HEX16,  -1,       LT_HEX32,  -1,       LT_HEX64 }; // LT_DYN_INT_H

rom dyn_int_lt_order_name (uint8_t dyn_lt_order)
{
    return lt_name (lt_order[dyn_lt_order]);
}

// ZIP only
LocalType dyn_int_get_ltype (ContextP ctx)
{
    LocalType actual_ltype;

    switch (ctx->ltype) {
        case LT_DYN_INT   : actual_ltype = lt_order    [ctx->dyn_lt_order]; break;
        case LT_DYN_INT_h : actual_ltype = lt_order_hex[ctx->dyn_lt_order]; break;
        case LT_DYN_INT_H : actual_ltype = lt_order_HEX[ctx->dyn_lt_order]; break;
        default           : actual_ltype = ctx->ltype;
    }

    ASSERT (actual_ltype != (LocalType)-1, "Invalid actual_ltype for %.*s.ltype=%s and dyn_lt_order=%u", 
            MAX_TAG_LEN-1, ctx->tag_name, lt_name(ctx->ltype), ctx->dyn_lt_order);

    return actual_ltype;
}

// ZIP only. note: in PIZ, these are untransposed in eg BGEN_transpose_u32_buf
void dyn_int_transpose (VBlockP vb, ContextP ctx)
{
    uint32_t cols;
    ContextP copy_ctx = CTX(VCF_COPY_SAMPLE);

    // note: caller needs to set local.n_cols to number of columns. or if it remains 0, it is interprets as vcf_num_samples
    if (ctx->local.n_cols) {
        cols = ctx->local.n_cols;

        // we're restricted to 255 columns, because this number goes into uint8_t SectionHeaderCtx.param
        ASSERT (cols >= 0 && cols <= 255, "columns=%u ∉ [1,255] in transposed matrix %s", cols, ctx->tag_name);

        ctx->local_param = true; // SectionHeaderCtx.param carries cols
    }
    else {
        ASSERT (VB_DT(VCF), "%s: cols=0 for ctx=%s", VB_NAME, ctx->tag_name);
        cols = vcf_header_get_num_samples(); // not restricted to 255
    } 

    bool *missing = NULL; // if non-NULL: lines x cols matrix - element is true if sample was copied so its value is missing in local

    // case: VCF with vcf_sample_copys: all samples have data, but some is absent from local due vcf_sample_copy
    if (segconf.vcf_sample_copy && 
        ctx != copy_ctx         && // exclude copy_ctx itself: it is transposed, but not part of the sample
        copy_ctx->local.len     && // some samples were copied (see vcf_copy_sample_seg_finalize)
        vb->lines.len32 * cols == copy_ctx->num_samples_copied + ctx->local.len32) // all samples in all variants have this field
        
        missing = B1ST (bool, copy_ctx->local); // 0 and 1 are pre-defined node indicies (see vcf_copy_sample_seg_initialize)

    // case: local is not a rectangle so not transposable 
    else if (ctx->local.len32 % cols) {
        if (ctx->local.n_cols) ctx->local_param = false; // undo
        return; 
    }

    uint32_t rows = missing ? vb->lines.len32 
                  :           (ctx->local.len32 / cols); // if vcf_sample_copy not used, we allow some rows to not have data at all 

    switch (ctx->ltype) { // note: the casting also correctly converts 0xffffffff to eg 0xff
        #define case_width(n)                                                                       \
        case LT_UINT##n: {                                                                          \
            uint##n##_t *data = (uint##n##_t *)ctx->local.data;                                     \
            ARRAY_alloc (uint##n##_t, trans_full, rows * cols, false, vb->scratch, vb, "scratch");  \
            for (uint32_t r=0; r < rows; r++)                                                       \
                for (uint32_t c=0; c < cols; c++)                                                   \
                    if (!missing || !(*missing++))                                                  \
                        trans_full[c * rows + r] = *data++; /* note: if missing, we set only the elements of scratch which were are available in local (i.e. not copied), leaving the remaining scratch elements uninitialized */ \
            ctx->ltype = missing ? LT_UINT##n##_PTR : LT_UINT##n##_TR;                              \
            break;                                                                                  \
        }

        case_width(8);
        case_width(16);
        case_width(32);

        default: ABORT ("Bad ltype=%s in %s", lt_name (ctx->ltype), ctx->tag_name);
    }

    // case: copy back transposed array: rows X cols elements
    if (!missing)
        buf_copy_do (vb, &ctx->local, &vb->scratch, lt_width(ctx), 0, 0, __FUNCLINE, C_LOCAL); // copy and not move, so we can keep local's memory for next vb

    // case: copy to local only the available data (i.e. not uninitialized scratch elements due to copied samples)
    else {
        missing = B1ST (bool, copy_ctx->local); // re-init
         
        switch (ctx->ltype) { 
            #define case_width_copy(n)                                      \
            case LT_UINT##n##_PTR: {                                        \
                uint##n##_t *data = (uint##n##_t *)ctx->local.data;         \
                uint##n##_t *trans_full = (uint##n##_t *)vb->scratch.data;  \
                for (uint32_t c=0; c < cols; c++)                           \
                    for (uint32_t r=0; r < rows; r++)                       \
                        if (!missing[r * cols + c])                         \
                            *data++ = trans_full[c * rows + r];             \
                ASSERT (BNUM(ctx->local, data) == ctx->local.len32, "bad copy: bnum=%u len=%u", BNUM(ctx->local, data), ctx->local.len32);/*sanity*/\
                break;                                                      \
            }

            case_width_copy(8);
            case_width_copy(16);
            case_width_copy(32);
            default: {} // already tested in previous switch
        }
    }

    buf_free (vb->scratch);
}

typedef void (*Resizer) (VBlockP, ContextP, LocalType, LocalType);
#define resize(src_type, dst_size)                                                                                  \
static inline void resize_##src_type##_to_##dst_size (VBlockP vb, ContextP ctx, LocalType src_lt, LocalType dst_lt) \
{                                                                                                                   \
    BufferP buf = IS_ZIP ? &ctx->local : &ctx->history;                                                             \
    if (IS_ZIP)                                                                                                     \
        buf_alloc (vb, buf, 0, (buf->len + 1), int##dst_size##_t, 0, C_LOCAL);                                      \
    else {                                                                                                          \
        bool initialize = (buf->len == 0);                                                                          \
        buf_alloc_exact (vb, *buf, vb->lines.len, int##dst_size##_t, C_HISTORY);                                    \
        if (initialize) memset (buf->data, 0, buf->len * sizeof (int##dst_size##_t));                               \
        SWAP (vb->line_i, buf->len32); /* resize only until before current line */                                  \
    }                                                                                                               \
                                                                                                                    \
    if (IS_PIZ || !ctx->nothing_char) {                                                                             \
        for_buf_tandem_back (src_type, src, *buf, int##dst_size##_t, dst, *buf)                                     \
            *dst = *src;                                                                                            \
    }                                                                                                               \
    else {                                                                                                          \
        int64_t src_max = lt_max (src_lt);                                                                          \
        int64_t dst_max = lt_max (dst_lt);                                                                          \
        for_buf_tandem_back (src_type, src, *buf, uint##dst_size##_t, dst, *buf)                                    \
            *dst = *src == src_max ? dst_max : *src; /* src/dst_max represent a period */                           \
    }                                                                                                               \
                                                                                                                    \
    if (IS_PIZ) SWAP (vb->line_i, buf->len32); /* restore */                                                        \
}

resize(uint8_t,  8 )
resize(uint8_t,  16)
resize(uint8_t,  32)
resize(uint8_t,  64)
resize(int8_t,   16)
resize(int8_t,   32)
resize(int8_t,   64)
resize(uint16_t, 16)
resize(uint16_t, 32)
resize(uint16_t, 64)
resize(int16_t,  32)
resize(int16_t,  64)
resize(uint32_t, 32)
resize(uint32_t, 64)
resize(int32_t,  64)

// ZIP/PIZ
static void dyn_int_resize (VBlockP vb, ContextP ctx, LocalType old_lt, LocalType new_lt)
{
    // if "resizing" from uint to int of the same size - no need to change local (unless period_as_int)
    // but we CANNOT skip ZIP is nothing_char and we CANNOT skip PIZ if not allocated yet
    if (((IS_ZIP && !ctx->nothing_char) || (IS_PIZ && ctx->history.len)) && 
         ((old_lt == LT_UINT8  && new_lt == LT_INT8 ) ||
          (old_lt == LT_UINT16 && new_lt == LT_INT16) ||
          (old_lt == LT_UINT32 && new_lt == LT_INT32))) return;

    static Resizer resizers[NUM_LTYPES/*src*/][NUM_LTYPES/*dst*/] =
        { [LT_UINT8]  = { [LT_INT8  ... LT_UINT8 ] = resize_uint8_t_to_8,
                          [LT_INT16 ... LT_UINT16] = resize_uint8_t_to_16,
                          [LT_INT32 ... LT_UINT32] = resize_uint8_t_to_32,
                          [LT_INT64]               = resize_uint8_t_to_64  },
          [LT_INT8]   = { [LT_INT16]               = resize_int8_t_to_16,
                          [LT_INT32]               = resize_int8_t_to_32,
                          [LT_INT64]               = resize_int8_t_to_64   },
          [LT_UINT16] = { [LT_INT16]               = resize_uint16_t_to_16,
                          [LT_INT32 ... LT_UINT32] = resize_uint16_t_to_32,
                          [LT_INT64]               = resize_uint16_t_to_64 },
          [LT_INT16]  = { [LT_INT32]               = resize_int16_t_to_32,
                          [LT_INT64]               = resize_int16_t_to_64  },
          [LT_UINT32] = { [LT_INT32]               = resize_uint32_t_to_32,
                          [LT_INT64]               = resize_uint32_t_to_64 },
          [LT_INT32]  = { [LT_INT64]               = resize_int32_t_to_64  } };
                         
    ASSERT (resizers[old_lt][new_lt], "missing resizer %s to %s", lt_name (old_lt), lt_name (new_lt));

    resizers[old_lt][new_lt] (vb, ctx, old_lt, new_lt);
}

// ZIP/PIZ
void dyn_int_init_ctx (VBlockP vb, ContextP ctx, int64_t value)
{
    if (ctx->dyn_lt_order) return; // already initialized

    ctx->dyn_lt_order = 1; // start with LT_UINT8 
    ctx->dyn_int_min = ctx->dyn_int_max = value;

    if (IS_ZIP) {
        // if set, segger *may* call dyn_int_append_nothing_char 
        if (VB_DT(VCF))
            ctx->nothing_char = '.';

        ASSERT (!ctx->local.len, "%s: expecting %s.local to be empty", LN_NAME, ctx->tag_name);

        if (ctx->ltype == LT_SINGLETON) 
            ctx->ltype = LT_DYN_INT;

        ASSERT (IS_LT_DYN (ctx->ltype), "%s: incompatible %s.ltype=%s", LN_NAME, ctx->tag_name, lt_name (ctx->ltype));

        // set STORE_INT 
        if (ctx->flags.store != STORE_INT) {
            ASSERT (ctx->flags.store == STORE_NONE, "%s: Incompatible %s.store_type=%s", LN_NAME, ctx->tag_name, store_type_name (ctx->flags.store));

            ctx->flags.store = STORE_INT;
        }
    }
}

// ZIP only
static inline void dyn_init_alloc (VBlockP vb, ContextP ctx)
{
    #define MIN_LOCAL_ALLOCATION 1024
    #define AT_LEAST(did_i) ROUNDUP64 (MAX_(MIN_LOCAL_ALLOCATION * sizeof(int64_t), ((uint64_t)(((did_i) < MAX_NUM_PREDEFINED) ? segconf.local_per_line[did_i] * (float)(vb->lines.len32) : 0))))

    // non-initial alloc. skip the complicated computation in AT_LEAST
    if (ctx->local.len32)
        buf_alloc (vb, &ctx->local, 1, MIN_LOCAL_ALLOCATION, int64_t, CTX_GROWTH, C_LOCAL); 

    // initial alloc. note: not in dyn_int_init_ctx, bc we initialize the context in seg_initialize, but then sometimes never use it
    else 
        buf_alloc (vb, &ctx->local, 1, AT_LEAST(ctx->did_i), char, CTX_GROWTH, C_LOCAL); // complicated math only on first call
}

// ZIP/PIZ
static LocalType dyn_init_prepare (VBlockP vb, ContextP ctx, int64_t value)
{
    // initialize if needed
    if (!ctx->dyn_lt_order) 
        dyn_int_init_ctx (vb, ctx, value);

    if (IS_ZIP) // note: in PIZ, history is allocated in dyn_int_resize 
        dyn_init_alloc (vb, ctx);

    if      (value < ctx->dyn_int_min) ctx->dyn_int_min = value;
    else if (value > ctx->dyn_int_max) ctx->dyn_int_max = value;

    LocalType dyn_ltype = lt_order[ctx->dyn_lt_order]; 

    uint8_t nothing_char = IS_ZIP ? ctx->nothing_char : 0; // we don't support storing nothing_char in PIZ in history, but no reason not to support if needed in the future 
    bool is_allocated = IS_ZIP || (ctx->history.len > 0);

    // case: the current ltype cannot accomodate the value - we will need to resize
    if (value < lt_min(dyn_ltype) || value > lt_max(dyn_ltype) - (nothing_char != 0) ||
        (IS_PIZ && !ctx->history.len)) {
        // search for the next order up, which this and all previous values can fit into
        bool resized = false;
        for (int i=ctx->dyn_lt_order + is_allocated; i < ARRAY_LEN(lt_order); i++) 
            if (ctx->dyn_int_min >= lt_min (lt_order[i]) && 
                ctx->dyn_int_max <= lt_max (lt_order[i]) - (nothing_char != 0)) { // maximum is reduced by 1 if nothing_char!=0 as nothing_char is represented as max_int

                if (IS_PIZ || ctx->local.len) 
                    dyn_int_resize (vb, ctx, dyn_ltype, lt_order[i]); // note: we will always get here, because at least int64_t always works 

                if (flag.debug_dyn_int) 
                    iprintf ("%s: %s[%u].%s resized from %s to %s due to %"PRId64"\n", 
                            VB_NAME, ctx->tag_name, ctx->did_i, (IS_ZIP ? "local" : "history"),
                            lt_name (dyn_ltype), lt_name (lt_order[i]), value);

                ctx->dyn_lt_order = i;
                dyn_ltype = lt_order[ctx->dyn_lt_order]; // update after resizing
                resized = true;
                break; 
            }

        ASSERT (resized, "value=%"PRId64" in ctx=%s%s cannot be resized", // expected only if nothing_char!=0 and value=LT_INT64.max_int
                value, ctx->tag_name, 
                cond_str (nothing_char, " with nothing_char=", char_to_printable (ctx->nothing_char).s));
    }

    return dyn_ltype;
}

// ZIP: appends ctx->local: requires setting ltype=LT_DYN_INT* in seg_initialize
void dyn_int_append (VBlockP vb, ContextP ctx, int64_t value, unsigned add_bytes)
{
    LocalType dyn_ltype = dyn_init_prepare (vb, ctx, value);

    // add value to local
    switch (dyn_ltype) {
        case LT_INT8  : case LT_UINT8  : BNXT (int8_t,  ctx->local) = (int8_t) value; break; // note: if value >= 128, converting it to int8_t will not change its binary representation 
        case LT_INT16 : case LT_UINT16 : BNXT (int16_t, ctx->local) = (int16_t)value; break;
        case LT_INT32 : case LT_UINT32 : BNXT (int32_t, ctx->local) = (int32_t)value; break;
        case LT_INT64 :                  BNXT (int64_t, ctx->local) =          value; break;            
        default : ABORT ("Unexpected dyn_ltype=%u", dyn_ltype);
    }
    
    if (add_bytes) ctx->txt_len += add_bytes;
    ctx->local_num_words++;
}

// ZIP only
// - this MAY be called for contexts marked that have nothing_char - currently VCF samples contexts
// - reconstruction of the period happens in reconstruct_from_local_int 
void dyn_int_append_nothing_char (VBlockP vb, ContextP ctx, unsigned add_bytes)
{
    // initialize if needed
    if (!ctx->dyn_lt_order) 
        dyn_int_init_ctx (vb, ctx, 0xff);

    LocalType dyn_ltype = lt_order[ctx->dyn_lt_order]; 

    dyn_init_alloc (vb, ctx);

    ASSSEG (ctx->nothing_char, "cannot use for %s because nothing_char=0", ctx->tag_name);

    switch (dyn_ltype) {
        case LT_INT8  : case LT_UINT8  : BNXT8  (ctx->local) = lt_max (dyn_ltype); break; 
        case LT_INT16 : case LT_UINT16 : BNXT16 (ctx->local) = lt_max (dyn_ltype); break;
        case LT_INT32 : case LT_UINT32 : BNXT32 (ctx->local) = lt_max (dyn_ltype); break;
        case LT_INT64 :                  BNXT64 (ctx->local) = lt_max (dyn_ltype); break;
        default : ABORT ("Unexpected dyn_ltype=%u", dyn_ltype);
    }
    
    if (add_bytes) ctx->txt_len += add_bytes;
    ctx->local_num_words++;
}

// PIZ: store value in history
void dyn_int_store_history (VBlockP vb, ContextP ctx, int64_t value)
{
    LocalType dyn_ltype = dyn_init_prepare (vb, ctx, value); // allocate / resize

    // add value to local
    switch (dyn_ltype) {
        case LT_INT8  : case LT_UINT8  : *B(int8_t,  ctx->history, vb->line_i) = (int8_t) value; break; // note: if value >= 128, converting it to int8_t will not change its binary representation 
        case LT_INT16 : case LT_UINT16 : *B(int16_t, ctx->history, vb->line_i) = (int16_t)value; break;
        case LT_INT32 : case LT_UINT32 : *B(int32_t, ctx->history, vb->line_i) = (int32_t)value; break;
        case LT_INT64 :                  *B(int64_t, ctx->history, vb->line_i) =          value; break;            
        default : ABORT ("Unexpected dyn_ltype=%u", dyn_ltype);
    }   
}

// PIZ: get value from history
int64_t piz_get_history (ContextP ctx, uint32_t line_i)
{
    switch (lt_order[ctx->dyn_lt_order]) {
        case LT_INT8   : return *B(int8_t,   ctx->history, line_i);
        case LT_UINT8  : return *B(uint8_t,  ctx->history, line_i);
        case LT_INT16  : return *B(int16_t,  ctx->history, line_i);
        case LT_UINT16 : return *B(uint16_t, ctx->history, line_i);
        case LT_INT32  : return *B(int32_t,  ctx->history, line_i);
        case LT_UINT32 : return *B(uint32_t, ctx->history, line_i);
        case LT_INT64  : return *B(int64_t,  ctx->history, line_i);
        default : ABORT ("Unexpected dyn_lt_order=%u", ctx->dyn_lt_order);
    }
} 
