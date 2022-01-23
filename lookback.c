#include "genozip.h"
#include "lookback.h"
#include "context.h"
#include "buffer.h"

// A round-robin buffer that can handle 1^num_bits-1 items. An item is inserted *before* the previous item, so
// that the newest item as the lowest index (modulo the size) and when we search for the most recent item, we search
// forward.

#define lookback_buf(ctx) ((command == ZIP) ? &ctx->zip_lookback_buf : &ctx->history)
#define lookback_size(lb_ctx) (1 << ((lb_ctx)->local.param + 10))
#define gap_index len // the index of the entry that is not used, one before (i.e. higher) that the oldest entry
#define newest_index param

#define RR(value, size) (((value) < 0) ? ((value)+(size)) : ((value)>= size) ? ((value)-(size)) : (value))

void lookback_init (VBlockP vb, ContextP lb_ctx, ContextP ctx, StoreType store_type)
{
    if (command == ZIP) {
        ctx->flags.store  = store_type; // tell PIZ store store values, so that the container callback can insert them to the lookback
        ctx->no_drop_b250 = true;       // we cannot have all_the_same, bc we need the b250 section to pass the param (lookback bits)
    }

    buf_alloc (vb, lookback_buf(ctx), 0, lookback_size(lb_ctx) * (store_type == STORE_INT ? sizeof (int64_t) : sizeof (WordIndex)), char, 1, "lookback_buf");
}

void lookback_insert (VBlockP vb, DidIType lb_did_i, DidIType did_i, bool copy_last_value, ValueType value, bool is_word_index)
{
    ContextP ctx = CTX(did_i);
    Buffer *buf = lookback_buf(ctx);
    uint32_t lb_size = lookback_size (CTX(lb_did_i));

    buf->newest_index = RR(buf->newest_index - 1, lb_size);

    // case: buffer is full, slide gap_index down, thereby discarding the oldest item
    if (buf->newest_index == buf->gap_index) 
        buf->gap_index = RR((int64_t)buf->gap_index - 1, lb_size);

    if (copy_last_value)
        value = ctx->last_value;

    if (is_word_index) 
        *ENT (WordIndex, *buf, buf->newest_index) = (WordIndex)value.i; // insert index
    else              
        *ENT (ValueType, *buf, buf->newest_index) = value;              // insert value
}

void lookback_insert_txt (VBlockP vb, DidIType lb_did_i, DidIType did_i, STRp(txt))
{ 
    ValueType value = { .txt = { .index =  ENTNUM (vb->txt_data, txt), .len = txt_len } };
    lookback_insert (vb, lb_did_i, did_i, false, value, false);
}

static inline unsigned lookback_len (ContextP ctx, uint32_t lb_size)
{
    Buffer *buf = lookback_buf(ctx);

    if (buf->newest_index <= buf->gap_index) 
        return buf->gap_index - buf->newest_index;
    else
        return buf->gap_index + lb_size - buf->newest_index;
}

const void *lookback_get_do (VBlockP vb, ContextP lb_ctx, ContextP ctx, 
                             unsigned lookback, // 1 means the newest item, 2 is 2nd newest etc
                             bool is_word_index)
{
    uint32_t lb_size = lookback_size (lb_ctx);

    ASSERT (lookback <= lookback_len (ctx, lb_size), "expecting lookback=%u <= lookback_len=%u for ctx=%s vb=%d line_i=%"PRIu64"%s%s lb_size=%u", 
            lookback, lookback_len(ctx, lb_size), ctx->tag_name, vb->vblock_i, vb->line_i, (VB_DT(DT_VCF) ? " sample_i=" : ""), (VB_DT(DT_VCF) ? str_int_s (vb->sample_i).s : ""), lb_size);
            
    Buffer *buf = lookback_buf(ctx);
    unsigned index = RR(buf->newest_index + lookback - 1, lb_size);

    ASSERT (lookback > 0 && lookback < lb_size, "Expecting lookback=%d in ctx=%s vb=%d line_i=%"PRIu64"%s%s to be in the range [1,%u]", 
            lookback, ctx->tag_name, vb->vblock_i, vb->line_i, (VB_DT(DT_VCF) ? " sample_i=" : ""), (VB_DT(DT_VCF) ? str_int_s (vb->sample_i).s : ""), lb_size-1);

    return is_word_index ? (void *)ENT (WordIndex, *buf, index) : (void *)ENT (int64_t, *buf, index);
}

// Seg: check if a string is the same of a back txt at a certain lookback
bool lookback_is_same_txt (VBlockP vb, DidIType lb_did_i, ContextP ctx, uint32_t lookback, STRp(str))
{
    ContextP lb_ctx = CTX(lb_did_i);
    uint32_t lb_size = lookback_size (lb_ctx);
    if (lookback > lookback_len (ctx, lb_size)) return false; // no lookup available - not enough lookback data yet

    ValueType value = lookback_get_value (vb, CTX(lb_did_i), ctx, lookback);

    return str_issame_(STRa(str), ENT(char, vb->txt_data, value.txt.index), value.txt.len);
}


// Returns the next lookup value that contains the WordIndex search_for, or 0 if there isn't one.  
uint32_t lookback_get_next (VBlockP vb, ContextP lb_ctx, ContextP ctx, WordIndex search_for, 
                            int64_t *iterator) // iterator should be initialized to -1 by caller. updates to the first item to be tested next call.
{
    Buffer *buf = lookback_buf(ctx);
    uint32_t lb_size = lookback_size (lb_ctx);

    if (buf->newest_index == buf->gap_index) return 0; // buffer is empty
    
    if (*iterator == -1) *iterator = buf->newest_index;
    uint32_t lookback=0; // initialize to "not found"

    for (; !lookback && *iterator != buf->gap_index ; *iterator = RR(*iterator + 1, lb_size))
        if (*ENT (WordIndex, *buf, *iterator) == search_for) 
            lookback = (RR(*iterator - buf->newest_index + 1, lb_size));

    ASSERT (lookback >= 0 && lookback < lb_size, "Invalid lookback=%d", lookback);
    return lookback;
}

void lookback_flush (VBlockP vb, ContextP ctx)
{
    Buffer *buf = lookback_buf(ctx);
    buf->gap_index = buf->newest_index = 0;
}

uint8_t lookback_size_to_local_param (uint32_t size)
{
    return size ? MAX_(0, 31 - __builtin_clz (2*size-1) - 10) : 0; // round up to the next power of 2
}
