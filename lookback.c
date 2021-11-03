#include "genozip.h"
#include "lookback.h"
#include "context.h"
#include "buffer.h"

// A round-robin buffer that can handle 1^num_bits-1 items. An item is inserted *before* the previous item, so
// that the newest item as the lowest index (modulo the size) and when we search for the most recent item, we search
// forward.

#define ZIP_LOOKBACK_SIZE 1024  // this can be changed up or down, up to PIZ_LOOKBACK_SIZE
#define PIZ_LOOKBACK_SIZE 1024  // this needs to be at least what ZIP_LOOKBACK_SIZE was at any Genozip release

#define lookback_buf(ctx) ((command == ZIP) ? &ctx->zip_lookback_buf : &ctx->history)
#define lookback_size (command == ZIP ? ZIP_LOOKBACK_SIZE : PIZ_LOOKBACK_SIZE)
#define gap_index len // the index of the entry that is not used, one before (i.e. higher) that the oldest entry
#define newest_index param

#define RR(value, size) (((value) < 0) ? ((value)+(size)) : ((value)>= size) ? ((value)-(size)) : (value))

void lookback_init (VBlockP vb, ContextP ctx, StoreType store_type)
{
    if (command == ZIP) {
        ctx->flags.store  = store_type; // tell PIZ store store values, so that the container callback can insert them to the lookback
        ctx->no_drop_b250 = true;       // we cannot have all_the_same, bc we need the b250 section to pass the param (lookback bits)
    }

    buf_alloc (vb, lookback_buf(ctx), 0, lookback_size * (store_type == STORE_INT ? sizeof (int64_t) : sizeof (WordIndex)), char, 1, "lookback_buf");
}

void lookback_insert (VBlockP vb, DidIType did_i, int64_t value, bool is_word_index)
{
    ContextP ctx = CTX(did_i);
    Buffer *buf = lookback_buf(ctx);

    buf->newest_index = RR(buf->newest_index - 1, lookback_size);

    // case: buffer is full, slide gap_index down, thereby discarding the oldest item
    if (buf->newest_index == buf->gap_index) 
        buf->gap_index = RR((int64_t)buf->gap_index - 1, lookback_size);

    if (value == TAKE_LAST_VALUE)
        value = ctx->last_value.i;

    if (is_word_index) 
        *ENT (WordIndex, *buf, buf->newest_index) = (WordIndex)value; // insert index
    else              
        *ENT (int64_t, *buf, buf->newest_index) = value;              // insert value
}

static inline unsigned lookback_len (ContextP ctx)
{
    Buffer *buf = lookback_buf(ctx);

    if (buf->newest_index <= buf->gap_index) 
        return buf->gap_index - buf->newest_index;
    else
        return buf->gap_index + lookback_size - buf->newest_index;
}

const void *lookback_get_do (VBlockP vb, ContextP ctx, 
                             unsigned lookback, // 1 means the newest item, 2 is 2nd newest etc
                             bool is_word_index)
{
    ASSERT (lookback <= lookback_len (ctx), "expecting lookback=%u <= lookback_len=%u for ctx=%s vb=%d line_i=%"PRIu64, 
            lookback, lookback_len(ctx), ctx->tag_name, vb->vblock_i, vb->line_i);
            
    Buffer *buf = lookback_buf(ctx);
    unsigned index = RR(buf->newest_index + lookback - 1, lookback_size);

    ASSERT (lookback > 0 && lookback < lookback_size, "Expecting lookback=%d in ctx=%s to be in the range [1,%u]", 
            lookback, ctx->tag_name, lookback_size-1);

    return is_word_index ? (void *)ENT (WordIndex, *buf, index) : (void *)ENT (int64_t, *buf, index);
}

// Returns the next lookup value that contains the WordIndex search_for, or 0 if there isn't one.  
uint32_t lookback_get_next (VBlockP vb, ContextP ctx, WordIndex search_for, 
                            int64_t *iterator) // iterator should be initialized to -1 by caller. updates to the first item to be tested next call.
{
    Buffer *buf = lookback_buf(ctx);

    if (buf->newest_index == buf->gap_index) return 0; // buffer is empty
    
    if (*iterator == -1) *iterator = buf->newest_index;
    uint32_t lookback=0; // initialize to "not found"

    for (; !lookback && *iterator != buf->gap_index ; *iterator = RR(*iterator + 1, lookback_size))
        if (*ENT (WordIndex, *buf, *iterator) == search_for) 
            lookback = (RR(*iterator - buf->newest_index + 1, lookback_size));

    ASSERT (lookback >= 0 && lookback < lookback_size, "Invalid lookback=%d", lookback);
    return lookback;
}

void lookback_flush (VBlockP vb, ContextP ctx)
{
    Buffer *buf = lookback_buf(ctx);
    buf->gap_index = buf->newest_index = 0;
}
