
// ------------------------------------------------------------------
//   container.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "container.h"
#include "vblock.h"
#include "data_types.h"
#include "base64.h"
#include "seg.h"
#include "piz.h"
#include "dict_id.h"

//----------------------
// Segmentation
//----------------------

WordIndex container_seg_by_ctx (VBlock *vb, Context *ctx, Container *con, 
                                 // prefixes can be one of 3 options:
                                 // 1. NULL
                                 // 2. a "container-wide prefix" that will be reconstructed once, at the beginning of the Container
                                 // 3. a "container-wide prefix" followed by exactly one prefix per item. the per-item prefixes will be
                                 //    displayed once per repeat, before their respective items. in this case, the container-wide prefix
                                 //    may be empty. 
                                 // Each prefix is terminated by a SNIP_CONTAINER character
                                 const char *prefixes, unsigned prefixes_len, // a container-wide prefix (may be empty), followed (or not) by one prefix per item. Each prefixes is terminated by SNIP_CONTAINER.
                                 unsigned add_bytes)
{
    con->repeats = BGEN32 (con->repeats);
    char snip[1 + base64_sizeof(Container) + CONTAINER_MAX_PREFIXES_LEN]; // maximal size
    snip[0] = SNIP_CONTAINER;
    unsigned b64_len = base64_encode ((uint8_t*)con, sizeof_container (*con), &snip[1]);
    con->repeats = BGEN32 (con->repeats); // restore

    if (prefixes_len) memcpy (&snip[1+b64_len], prefixes, prefixes_len);
    uint32_t snip_len = 1 + b64_len + prefixes_len;

    ctx->flags |= CTX_FL_CONTAINER;

    // store in struct cache
    return seg_by_ctx (vb, snip, snip_len, ctx, add_bytes, NULL); 
}

//----------------------
// Reconstruction
//----------------------

static inline void container_reconstruct_prefix (VBlockP vb, const char **prefixes, uint32_t *prefixes_len)
{
    if (! (*prefixes_len)) return; // nothing to do
    
    const char *start = *prefixes;
    while (**prefixes != SNIP_CONTAINER) (*prefixes)++; // prefixes are terminated by SNIP_CONTAINER
    uint32_t len = (unsigned)((*prefixes) - start);

    RECONSTRUCT (start, len);

    (*prefixes)++; // skip SNIP_CONTAINER seperator
    (*prefixes_len) -= len + 1;
}

static inline void container_reconstruct_do (VBlock *vb, DictId dict_id, const Container *con, const char *prefixes, uint32_t prefixes_len)
{
    TimeSpecType profiler_timer={0}; 
    if (flag_show_time && (con->flags & CONTAINER_TOPLEVEL)) 
        clock_gettime (CLOCK_REALTIME, &profiler_timer);

    // container wide prefix - it will be missing if Container has no prefixes, or empty if it has only items prefixes
    container_reconstruct_prefix (vb, &prefixes, &prefixes_len); 

    ASSERT (DTP (container_filter) || !(con->flags & CONTAINER_FILTER_REPEATS) || !(con->flags & CONTAINER_FILTER_ITEMS), 
            "Error: data_type=%s doesn't support container_filter", dt_name (vb->data_type));

    for (uint32_t rep_i=0; rep_i < con->repeats; rep_i++) {

        // case this is the top-level snip
        if (con->flags & CONTAINER_TOPLEVEL) {
            vb->line_i = vb->first_line + rep_i;
            vb->line_start = vb->txt_data.len;
            vb->dont_show_curr_line = false; 
        }
    
        if ((con->flags & CONTAINER_FILTER_REPEATS) && !(DT_FUNC (vb, container_filter) (vb, dict_id, con, rep_i, -1))) continue; // repeat is filtered out

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t item_prefixes_len = prefixes_len;

        for (unsigned i=0; i < con->num_items; i++) {

            if ((con->flags & CONTAINER_FILTER_ITEMS) && !(DT_FUNC (vb, container_filter) (vb, dict_id, con, rep_i, i))) continue; // item is filtered out

            container_reconstruct_prefix (vb, &item_prefixes, &item_prefixes_len); // item prefix (we will have one per item or none at all)

            const ContainerItem *item = &con->items[i];
            int32_t reconstructed_len=0;
            if (item->dict_id.num) {  // not a prefix-only item

                if (flag_show_containers) // show container reconstruction 
                    fprintf (stderr, "Line=%u Repeat=%u %.*s->%s\n", vb->line_i, rep_i, DICT_ID_LEN, dict_id_printable(dict_id).id, vb->contexts[item->did_i].name);
                
                reconstructed_len = piz_reconstruct_from_ctx (vb, item->did_i, item->seperator[0] != CI_DONT_RECONSTRUCT, 0);
            }            
            if (reconstructed_len == -1 && i > 0)  // not WORD_INDEX_MISSING_SF - delete previous item's separator
                vb->txt_data.len -= ((item-1)->seperator[0] != 0) + ((item-1)->seperator[1] != 0);

            // seperator determines what to do after the item
            switch (item->seperator[0]) {
                case CI_MOVE:
                    vb->txt_data.len += item->seperator[1]; break;
                case CI_NUL_TERMINATE:
                    RECONSTRUCT1 (0); break;
                case CI_DONT_RECONSTRUCT:
                    break;
                default: 
                    // note: we emit this item's separator even if this item is missing - next item will delete it if also missing (last item in Samples doesn't have a seperator)
                    if (item->seperator[0]) RECONSTRUCT1 (item->seperator[0]);
                    if (item->seperator[1]) RECONSTRUCT1 (item->seperator[1]);
            }
        }

        if (rep_i+1 < con->repeats || !(con->flags & CONTAINER_DROP_FINAL_REPEAT_SEP)) {
            if (con->repsep[0]) RECONSTRUCT1 (con->repsep[0]);
            if (con->repsep[1]) RECONSTRUCT1 (con->repsep[1]);
        }

        // in top level: after consuming the line's data, if it is not to be outputted - trim txt_data back to start of line
        if ((con->flags & CONTAINER_TOPLEVEL) && vb->dont_show_curr_line) 
            vb->txt_data.len = vb->line_start; 
    }

    if (con->flags & CONTAINER_DROP_FINAL_ITEM_SEP)
        vb->txt_data.len -= (con->items[con->num_items-1].seperator[0] != 0) + 
                            (con->items[con->num_items-1].seperator[1] != 0);

    if (con->flags & CONTAINER_TOPLEVEL)   
        COPY_TIMER (piz_reconstruct_vb);
}

void container_reconstruct (VBlock *vb, Context *ctx, WordIndex word_index, const char *snip, unsigned snip_len)
{
    ASSERT (snip_len <= base64_sizeof(Container), "Error in container_reconstruct: snip_len=%u exceed base64_sizeof(Container)=%u",
            snip_len, base64_sizeof(Container));

    Container con, *con_p=NULL;
    const char *prefixes;
    uint32_t prefixes_len;

    bool cache_exists = buf_is_allocated (&ctx->con_cache);
    uint16_t cache_item_len;

    // if this container exists in the cache - use the cached one
    if (cache_exists && word_index != WORD_INDEX_NONE && ((cache_item_len = *ENT (uint16_t, ctx->con_len, word_index)))) {
        con_p = (Container *)ENT (char, ctx->con_cache, *ENT (uint32_t, ctx->con_index, word_index));
        
        unsigned st_size = sizeof_container (*con_p);
        prefixes = (char *)con_p + st_size; // prefixes are stored after the Container
        prefixes_len = cache_item_len - st_size;
    }

    // case: not cached - decode it and optionally cache it
    if (!con_p) {
        // decode
        unsigned b64_len = snip_len;
        base64_decode (snip, &b64_len, (uint8_t*)&con);
        con.repeats = BGEN32 (con.repeats);

        // get the did_i for each dict_id
        for (uint8_t item_i=0; item_i < con.num_items; item_i++)
            if (con.items[item_i].dict_id.num)  // not a prefix-only item
                con.items[item_i].did_i = mtf_get_existing_did_i (vb, con.items[item_i].dict_id);

        // get prefixes
        unsigned st_size = sizeof_container (con);
        con_p         = &con;
        prefixes     = (b64_len < snip_len) ? &snip[b64_len+1]       : NULL;
        prefixes_len = (b64_len < snip_len) ? snip_len - (b64_len+1) : 0;

        ASSERT (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "Error in container_reconstruct: prefixes_len=%u longer than CONTAINER_MAX_PREFIXES_LEN=%u", 
                prefixes_len, CONTAINER_MAX_PREFIXES_LEN);

        ASSERT (!prefixes_len || prefixes[prefixes_len-1] == SNIP_CONTAINER, "Error in container_reconstruct: prefixes array does end with a SNIP_CONTAINER: %.*s",
                prefixes_len, prefixes);

        // cache it if it is cacheable 
        if (word_index != WORD_INDEX_NONE) {

            ASSERT (st_size + prefixes_len <= 65535, "st_size=%u + prefixes_len=%u too large", st_size, prefixes_len);

            // first encounter with Container for this context - allocate the cache
            if (!cache_exists) {
                buf_alloc (vb, &ctx->con_index, ctx->word_list.len * sizeof (uint32_t), 1, "context->con_index", ctx->did_i);
                buf_alloc (vb, &ctx->con_len,   ctx->word_list.len * sizeof (uint16_t),  1, "context->con_len",   ctx->did_i);
                buf_zero (&ctx->con_len);
            }

            // place Container followed by prefix in the cache
            *ENT (uint32_t, ctx->con_index, word_index) = (uint32_t)ctx->con_cache.len;
            *ENT (uint16_t, ctx->con_len,   word_index) = (uint16_t)(st_size + prefixes_len);

            buf_alloc (vb, &ctx->con_cache, ctx->con_cache.len + st_size + prefixes_len, 2, "context->con_cache", ctx->did_i);
            buf_add (&ctx->con_cache, con_p, st_size);
            if (prefixes_len) buf_add (&ctx->con_cache, prefixes, prefixes_len);
        }
    }

    container_reconstruct_do (vb, ctx->dict_id, con_p, prefixes, prefixes_len); 
}
