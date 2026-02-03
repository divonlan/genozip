// ------------------------------------------------------------------
//   lookback.c
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "lookback.h"
#include "context.h"

// A round-robin buffer that can handle 1^num_bits-1 items. An item is inserted *before* the previous item, so
// that the newest item as the lowest index (modulo the size) and when we search for the most recent item, we search
// forward.

#define lookback_size(lb_ctx) (1 << ((lb_ctx)->local.prm8[0] + 10))

#define RR(value, size) (((value) < 0) ? ((value)+(size)) : ((value)>= size) ? ((value)-(size)) : (value))

void lookback_init (VBlockP vb, ContextP lb_ctx, ContextP ctx, StoreType store_type)
{
    if (IS_ZIP) {
        ASSERT (!ctx->is_initialized, "Context %s already initialized", ctx->tag_name);

        ctx->flags.store    = store_type; // tell PIZ store store values, so that the container callback can insert them to the lookback
        ctx->no_drop_b250   = true;       // we cannot have all_the_same, bc we need the b250 section to pass the param (lookback bits)
        ctx->is_initialized = true;
    }

    buf_alloc (vb, &ctx->lookback, 0, lookback_size(lb_ctx) * (store_type == STORE_INDEX ? sizeof (WordIndex) : sizeof (ValueType)), char, 1, "contexts->lookback_buf");
 }

// Seg and PIZ
void lookback_insert (VBlockP vb, Did lb_did_i, Did did_i, bool copy_last_value, ValueType value)
{
    decl_ctx (did_i);
    uint32_t lb_size = lookback_size (CTX(lb_did_i));

    ctx->lookback.newest_index = RR(ctx->lookback.newest_index - 1, lb_size);

    // case: buffer is full, slide gap_index down, thereby discarding the oldest item
    if (ctx->lookback.newest_index == ctx->lookback.gap_index) 
        ctx->lookback.gap_index = RR((int64_t)ctx->lookback.gap_index - 1, lb_size);

    if (copy_last_value)
        value = ctx->last_value;

    if (ctx->flags.store == STORE_INDEX) 
        *B(WordIndex, ctx->lookback, ctx->lookback.newest_index) = (WordIndex)value.i; // insert index
    else              
        *B(ValueType, ctx->lookback, ctx->lookback.newest_index) = value;              // insert value
}

static inline unsigned lookback_len (ContextP ctx, uint32_t lb_size)
{
    if (ctx->lookback.newest_index <= ctx->lookback.gap_index) 
        return ctx->lookback.gap_index - ctx->lookback.newest_index;
    else
        return ctx->lookback.gap_index + lb_size - ctx->lookback.newest_index;
}

const void *lookback_get_do (VBlockP vb, ContextP lb_ctx, ContextP ctx, 
                             unsigned lookback) // 1 means the newest item, 2 is 2nd newest etc
{
    uint32_t lb_size = lookback_size (lb_ctx);

    ASSERT (lookback <= lookback_len (ctx, lb_size), "%s: expecting lookback=%u <= lookback_len=%u for ctx=%s%s lb_size=%u", 
            LN_NAME, lookback, lookback_len(ctx, lb_size), ctx->tag_name, cond_int (VB_DT(VCF), " sample_i=", vb->sample_i), lb_size);
            
    unsigned index = RR(ctx->lookback.newest_index + lookback - 1, lb_size);

    // cases where we segged "SNIP_LOOKBACK" when there is no lookback, to improve compression and knowing that we won't be using this value
    if (lookback == 0 && ctx->flags.lookback0_ok) {
        static const int64_t dummy_value = 0;
        return &dummy_value;
    }

    ASSERT (lookback > 0 && lookback < lb_size, "%s: Expecting lookback=%d in ctx=%s%s to be in the range [1,%u]", 
            LN_NAME, lookback, ctx->tag_name, cond_int (VB_DT(VCF), " sample_i=", vb->sample_i), lb_size-1);

    return (ctx->flags.store == STORE_INDEX) ? (void *)B(WordIndex, ctx->lookback, index) 
                                             : (void *)B(ValueType, ctx->lookback, index);
}

// shift existing lookups after insertion into txt_data
void lookback_shift_txt_index (VBlockP vb, ContextP lb_ctx, ContextP ctx, STRp (insert))
{
    if (!buf_is_alloc (&ctx->lookback)) return;

    uint32_t lb_size = lookback_size (lb_ctx);
    unsigned lb_len = lookback_len (ctx, lb_size);        

    for (unsigned lookback=1; lookback <= lb_len; lookback++) {
        unsigned index = RR(ctx->lookback.newest_index + lookback - 1, lb_size);
        ValueType *value = B(ValueType, ctx->lookback, index);
        
        if (value->index > BNUMtxt (insert)) // this lookback is after the insertion, therefore affected by it
            value->index += insert_len;

        else 
            break; // we reached the lookbacks that are before the insertion, and thus unaffected it 
    }
}

// Seg: check if a string is the same of a back txt at a certain lookback
bool lookback_is_same_txt (VBlockP vb, Did lb_did_i, ContextP ctx, uint32_t lookback, STRp(str))
{
    ContextP lb_ctx = CTX(lb_did_i);
    uint32_t lb_size = lookback_size (lb_ctx);
    if (lookback > lookback_len (ctx, lb_size)) return false; // no lookup available - not enough lookback data yet

    ValueType value = lookback_get_value (vb, CTX(lb_did_i), ctx, lookback);

    return str_issame_(STRa(str), STRtxt(value));
}


// Seg: Returns the next lookup value that contains the WordIndex search_for, or 0 if there isn't one.  
uint32_t lookback_get_next (VBlockP vb, ContextP lb_ctx, ContextP ctx, WordIndex search_for, 
                            int64_t *iterator) // iterator should be initialized to -1 by caller. updates to the first item to be tested next call.
{
    uint32_t lb_size = lookback_size (lb_ctx);

    if (ctx->lookback.newest_index == ctx->lookback.gap_index) return 0; // buffer is empty
    
    if (*iterator == -1) *iterator = ctx->lookback.newest_index;
    uint32_t lookback=0; // initialize to "not found"

    for (; !lookback && *iterator != ctx->lookback.gap_index ; *iterator = RR(*iterator + 1, lb_size))
        if (*B(WordIndex, ctx->lookback, *iterator) == search_for) 
            lookback = (RR(*iterator - ctx->lookback.newest_index + 1, lb_size));

    ASSERTINRANGE (lookback, 0, lb_size);
    return lookback;
}

uint8_t lookback_size_to_local_param (uint32_t size)
{
    return size ? MAX_(0, 31 - __builtin_clz (2*size-1) - 10) : 0; // round up to the next power of 2
}

// ZIP
void lookback_flush (VBlockP vb, ConstMediumContainerP con)
{
    for (unsigned i=1; i < con->nitems_lo; i++)
        if (con->items[i].separator[1] == CI1_LOOKBACK) {
            ContextP ctx = ctx_get_ctx (vb, con->items[i].dict_id);
            ctx->lookback.gap_index = ctx->lookback.newest_index = 0;
        }
}

// PIZ: insert all values of lookbackable items in a container
void lookback_insert_container (VBlockP vb, ConstContainerP con, unsigned num_items, ContextP *item_ctxs)
{        
    if (!item_ctxs[0]->is_initialized) {
        for (unsigned i=1; i < num_items; i++)
            if (con->items[i].separator[1] == CI1_LOOKBACK)
                lookback_init (vb, item_ctxs[0], item_ctxs[i], item_ctxs[i]->flags.store);

        item_ctxs[0]->is_initialized = true;
    }

    for (unsigned i=1; i < num_items; i++) 
        if (con->items[i].separator[1] == CI1_LOOKBACK) {
            if (item_ctxs[i]->flags.store != STORE_LAST_TXT)
                lookback_insert (vb, item_ctxs[0]->did_i, item_ctxs[i]->did_i, true, NO_VALUE); // copy last_value
            else
                lookback_insert (vb, item_ctxs[0]->did_i, item_ctxs[i]->did_i, false, item_ctxs[i]->last_txt);
        }
}
