
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
#include "reconstruct.h"
#include "dict_id.h"
#include "endianness.h"
#include "file.h"

//----------------------
// Segmentation
//----------------------

WordIndex container_seg_by_ctx (VBlock *vb, Context *ctx, ContainerP con, 
                                // prefixes can be one of 3 options:
                                // 1. NULL
                                // 2. a "container-wide prefix" that will be reconstructed once, at the beginning of the Container
                                // 3. a "container-wide prefix" followed by exactly one prefix per item. the per-item prefixes will be
                                //    displayed once per repeat, before their respective items. in this case, the container-wide prefix
                                //    may be empty. 
                                // Each prefix is terminated by a CON_PREFIX_SEP character
                                const char *prefixes, unsigned prefixes_len, // a container-wide prefix (may be empty), followed (or not) by one prefix per item. Each prefixes is terminated by CON_PREFIX_SEP.
                                unsigned add_bytes)
{
    ctx->no_stons = true; // we need the word index to for container caching

    // con=NULL means MISSING Container (see container_reconstruct_do)
    if (!con) return seg_by_ctx (vb, NULL, 0, ctx, 0); 

    ASSERTE (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "prefixes_len=%u is beyond maximum of %u in ctx=%s",
             prefixes_len, CONTAINER_MAX_PREFIXES_LEN, ctx->name);

    con->repeats = BGEN24 (con->repeats); 
    unsigned con_size = con_sizeof (*con);
    char snip[1 + base64_size(con_size) + prefixes_len]; 
    snip[0] = SNIP_CONTAINER;
    
    unsigned b64_len = base64_encode ((uint8_t*)con, con_size, &snip[1]);
    ASSERTE (b64_len <= base64_size(con_size), "b64_len=%u larger than base64_size(%u)=%u in ctx=%s",
             b64_len, con_size, base64_size(con_size), ctx->name);
    
    con->repeats = BGEN24 (con->repeats); // restore

    if (prefixes_len) memcpy (&snip[1+b64_len], prefixes, prefixes_len);
    uint32_t snip_len = 1 + b64_len + prefixes_len;

    if (flag.show_containers) { 
        fprintf (info_stream, "VB=%u Line=%u Ctx=%u:%s Repeats=%u RepSep=%u,%u Items=", 
                 vb->vblock_i, vb->line_i, ctx->did_i, ctx->name, con->repeats, con->repsep[0], con->repsep[1]);
        for (unsigned i=0; i < con_nitems (*con); i++) {
            ContainerItem *item = &con->items[i];
            if (item->dict_id.num) 
                fprintf (info_stream, "%u:%s ", ctx_get_existing_ctx (vb, item->dict_id)->did_i, dis_dict_id (item->dict_id).s); 
        }
        iprint0 ("\n");
    }

    // store in struct cache
    return seg_by_ctx (vb, snip, snip_len, ctx, add_bytes); 
}

//----------------------
// Reconstruction
//----------------------

static inline void container_reconstruct_prefix (VBlockP vb, ConstContainerP con, const char **prefixes, uint32_t *prefixes_len)
{
    ASSERTE (*prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "prefixes_len=%u is too big", *prefixes_len);

    if (! (*prefixes_len)) return; // nothing to do
    
    const char *start = *prefixes;
    while (**prefixes != CON_PREFIX_SEP && **prefixes != CON_PREFIX_SEP_SHOW_REPEATS) (*prefixes)++; // prefixes are terminated by CON_PREFIX_SEP
    uint32_t len = (unsigned)((*prefixes) - start);

    if (len) RECONSTRUCT (start, len);

    // if the seperator is CON_PREFIX_SEP_SHOW_REPEATS, output the number of repeats. This is for BAM 'B' array 'count' field.
    if (**prefixes == CON_PREFIX_SEP_SHOW_REPEATS)
        RECONSTRUCT_BIN32 (con->repeats);

    (*prefixes)++; // skip seperator
    (*prefixes_len) -= len + 1;
}

static inline LastValueType container_reconstruct_do (VBlock *vb, Context *ctx, ConstContainerP con, 
                                                      const char *prefixes, uint32_t prefixes_len)
{
    #define IS_CI_SET(flag) (CI_ITEM_HAS_FLAG (item) && ((uint8_t)item->seperator[0] & ~(uint8_t)0x80 & flag))

    TimeSpecType profiler_timer = {}; 
    if (flag.show_time && con->is_toplevel) 
        clock_gettime (CLOCK_REALTIME, &profiler_timer);
    
    int32_t last_non_filtered_item_i = -1;

    // container wide prefix - it will be missing if Container has no prefixes, or empty if it has only items prefixes
    container_reconstruct_prefix (vb, con, &prefixes, &prefixes_len); 

    ASSERTE (DTP (container_filter) || (!con->filter_repeats && !con->filter_items), 
             "data_type=%s doesn't support container_filter, despite being specified in the Container", dt_name (vb->data_type));

    uint32_t num_items = con_nitems(*con);
    Context *item_ctxs[num_items];
    
    // we can cache did_i up to 254. dues historical reasons the field is only 8 bit. that's enough in most cases anyway.
    for (unsigned i=0; i < num_items; i++) 
        item_ctxs[i] = !con->items[i].dict_id.num      ? NULL 
                     : con->items[i].did_i_small < 255 ? &vb->contexts[con->items[i].did_i_small]
                     :                                   ctx_get_existing_ctx (vb, con->items[i].dict_id);

    // for containers, new_value is the some of all its items, all repeats last_value (either int or float)
    LastValueType new_value = {};

    for (uint32_t rep_i=0; rep_i < con->repeats; rep_i++) {

        // case this is the top-level snip
        if (con->is_toplevel) {
            vb->line_i = vb->first_line + rep_i;
            vb->line_start = vb->txt_data.len;

            // show (or not) the line based on our downsampling rate
            vb->dont_show_curr_line = flag.downsample && (vb->line_i % flag.downsample); 
        }
    
        if (con->filter_repeats && !(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, -1))) continue; // repeat is filtered out

        char *rep_reconstruction_start = AFTERENT (char, vb->txt_data);

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t item_prefixes_len = prefixes_len;

        last_non_filtered_item_i = -1;
        for (unsigned i=0; i < num_items; i++) {
            const ContainerItem *item = &con->items[i];

            if (con->filter_items && !(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, i))) continue; // item is filtered out

            last_non_filtered_item_i = i;

            if (flag.show_containers && item_ctxs[i]) // show container reconstruction 
                iprintf ("VB=%u Line=%u Repeat=%u %s->%s txt_data.len=%"PRIu64" (0x%04"PRIx64") (BEFORE)\n", 
                         vb->vblock_i, vb->line_i, rep_i, dis_dict_id (ctx->dict_id).s, item_ctxs[i]->name, 
                         vb->vb_position_txt_file + vb->txt_data.len, vb->vb_position_txt_file + vb->txt_data.len);

            container_reconstruct_prefix (vb, con, &item_prefixes, &item_prefixes_len); // item prefix (we will have one per item or none at all)

            int32_t reconstructed_len=0;
            if (item->dict_id.num) {  // not a prefix-only or translator-only item
                char *reconstruction_start = AFTERENT (char, vb->txt_data);
                bool reconstruct = !flag.trans_containers ||  // not translating Or...
                                   !IS_CI_SET (CI_TRANS_NOR); // no prohibition on reconstructing when translating

                reconstructed_len = reconstruct_from_ctx (vb, item_ctxs[i]->did_i, 0, reconstruct);

                // sum up items' values if needed
                if (ctx->flags.store == STORE_INT)
                    new_value.i += item_ctxs[i]->last_value.i;
                
                else if (ctx->flags.store == STORE_FLOAT)
                    new_value.f += item_ctxs[i]->last_value.f;

                // if we're reconstructing to a translated format (eg SAM2BAM) - re-reconstruct this item
                // using the designated "translator" function, if one is available
                if (flag.trans_containers && item->translator) 
                    DT_FUNC(vb, translator)[item->translator](vb, item_ctxs[i], reconstruction_start, reconstructed_len);  
            }            

            // case: WORD_INDEX_MISSING_SF - delete previous item's separator if it has one (used by SAM_OPTIONAL - sam_seg_optional_all)
            if (reconstructed_len == -1 && i > 0 && !con->keep_empty_item_sep && !CI_ITEM_HAS_FLAG(item-1))  
                vb->txt_data.len -= ((item-1)->seperator[0] != 0) + ((item-1)->seperator[1] != 0);

            // seperator / flags determines what to do after the item            
            if (flag.trans_containers && IS_CI_SET (CI_TRANS_NUL))
                RECONSTRUCT1 (0);

            else if (!flag.trans_containers && IS_CI_SET (CI_NATIVE_NEXT))
                RECONSTRUCT1 (item->seperator[1]);

            else if (!CI_ITEM_HAS_FLAG(item)) { // seperator, not flag
                if (item->seperator[0]) RECONSTRUCT1 ((char)item->seperator[0]);
                if (item->seperator[1]) RECONSTRUCT1 ((char)item->seperator[1]);
            }

            // after all reconstruction and translation is done - move if needed
            if (flag.trans_containers && IS_CI_SET (CI_TRANS_MOVE))
                vb->txt_data.len += (uint8_t)item->seperator[1];
        }

        if (rep_i+1 < con->repeats || !con->drop_final_repeat_sep) {
            if (con->repsep[0]) RECONSTRUCT1 (con->repsep[0]);
            if (con->repsep[1]) RECONSTRUCT1 (con->repsep[1]);
        }

        // call callback if needed now that repeat reconstruction is done
        if (con->callback)
            DT_FUNC(vb, container_cb)(vb, ctx->dict_id, rep_i, rep_reconstruction_start, AFTERENT (char, vb->txt_data) - rep_reconstruction_start);

        // in top level: after consuming the line's data, if it is not to be outputted - trim txt_data back to start of line
        if (con->is_toplevel && vb->dont_show_curr_line) 
            vb->txt_data.len = vb->line_start; 
    }

    // remove final seperator, if we need to
    if (con->drop_final_item_sep && last_non_filtered_item_i >= 0) {
        const ContainerItem *item = &con->items[last_non_filtered_item_i]; // last_non_filtered_item_i is the last item that survived the filter, of the last repeat

        vb->txt_data.len -= CI_ITEM_HAS_FLAG(item) ? (flag.trans_containers ? 0 : !!IS_CI_SET (CI_NATIVE_NEXT))
                                                   : (!!item->seperator[0] + !!item->seperator[1]);
    }
     
    if (con->is_toplevel)   
        COPY_TIMER (reconstruct_vb);

    return new_value;
}

LastValueType container_reconstruct (VBlock *vb, Context *ctx, WordIndex word_index, const char *snip, unsigned snip_len)
{
    Container con, *con_p=NULL;
    const char *prefixes;
    uint32_t prefixes_len;

    bool cache_exists = buf_is_allocated (&ctx->con_cache);
    uint16_t cache_item_len;

    // if this container exists in the cache - use the cached one
    if (cache_exists && word_index != WORD_INDEX_NONE && ((cache_item_len = *ENT (uint16_t, ctx->con_len, word_index)))) {
        con_p = (ContainerP)ENT (char, ctx->con_cache, *ENT (uint32_t, ctx->con_index, word_index));
        
        unsigned st_size = con_sizeof (*con_p);
        prefixes = (char *)con_p + st_size; // prefixes are stored after the Container
        prefixes_len = cache_item_len - st_size;
    }

    // case: not cached - decode it and optionally cache it
    if (!con_p) {
        // decode
        unsigned b64_len = snip_len; // maximum length of b64 - it will be shorter if snip includes prefixes too
        base64_decode (snip, &b64_len, (uint8_t*)&con);
        con.repeats = BGEN24 (con.repeats);

        ASSERTE (con_nitems (con) <= MAX_SUBFIELDS, "A container of %s has %u items which is beyond MAX_SUBFIELDS=%u. Please upgrade to latest version of genozip to access this file.",
                 ctx->name, con_nitems (con), MAX_SUBFIELDS);

        // get the did_i for each dict_id - unfortunately we can only store did_i up to 254 (changing this would be a change in the file format)
        for (uint32_t item_i=0; item_i < con_nitems (con); item_i++)
            if (con.items[item_i].dict_id.num) { // not a prefix-only item
                DidIType did_i = ctx_get_existing_did_i (vb, con.items[item_i].dict_id);
                ASSERTE (did_i != DID_I_NONE, "analyzing a %s container: unable to find did_i for item %s",
                         ctx->name, dis_dict_id (con.items[item_i].dict_id).s);

                con.items[item_i].did_i_small = (did_i <= 254 ? (uint8_t)did_i : 255);
            }
            else 
                con.items[item_i].did_i_small = 255;
 
        // get prefixes
        unsigned st_size = con_sizeof (con);
        con_p         = &con;
        prefixes     = (b64_len < snip_len) ? &snip[b64_len+1]       : NULL;
        prefixes_len = (b64_len < snip_len) ? snip_len - (b64_len+1) : 0;

        ASSERTE (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "ctx=%s: prefixes_len=%u longer than CONTAINER_MAX_PREFIXES_LEN=%u", 
                 ctx->name, prefixes_len, CONTAINER_MAX_PREFIXES_LEN);

        ASSERTE (!prefixes_len || prefixes[prefixes_len-1] == CON_PREFIX_SEP, 
                 "ctx=%s: prefixes array does end with a CON_PREFIX_SEP: %.*s", ctx->name, prefixes_len, prefixes);

        // condition testing needed only for backward compatability with 8.0.x: prior to version 8.1 containers didn't have NO_STONS 
        // so singletons were in ctx.local without a word_index)
        if (word_index != WORD_INDEX_NONE) {

            ASSERTE (st_size + prefixes_len <= 65535, "ctx=%s: st_size=%u + prefixes_len=%u too large", 
                     ctx->name, st_size, prefixes_len);

            // first encounter with Container for this context - allocate the cache
            if (!cache_exists) {
                buf_alloc_more (vb, &ctx->con_index, 0, ctx->word_list.len, uint32_t, 1, "context->con_index");
                buf_alloc_more (vb, &ctx->con_len,   0, ctx->word_list.len, uint16_t, 1, "context->con_len");
                buf_zero (&ctx->con_len); // zero entire buffer, not just part added
            }

            // place Container followed by prefix in the cache
            *ENT (uint32_t, ctx->con_index, word_index) = (uint32_t)ctx->con_cache.len;

            buf_alloc (vb, &ctx->con_cache, ctx->con_cache.len + st_size + prefixes_len + CONTAINER_MAX_SELF_TRANS_CHANGE, 2, "context->con_cache");
            
            char *cached_con = AFTERENT (char, ctx->con_cache);
            buf_add (&ctx->con_cache, con_p, st_size);
            if (prefixes_len) buf_add (&ctx->con_cache, prefixes, prefixes_len);

            con_p    = (ContainerP)cached_con; // update so we reconstruct from translated cached item
            prefixes = cached_con + st_size;

            // if item[0] is a translator-only item, use it to translate the Container itself (used by SAM_OPTIONAL)
            ContainerItem *item0 = &con.items[0];
            if (flag.trans_containers && !item0->dict_id.num && item0->translator) {
                int32_t prefixes_len_change = DT_FUNC(vb, translator)[item0->translator](vb, ctx, cached_con, st_size + prefixes_len);  
                ASSERTE (prefixes_len_change <= CONTAINER_MAX_SELF_TRANS_CHANGE, 
                         "ctx=%s: prefixes_len_change=%d exceeds range maximum %u", 
                         ctx->name, prefixes_len_change, CONTAINER_MAX_SELF_TRANS_CHANGE);
                
                prefixes_len       += prefixes_len_change;
                ctx->con_cache.len += prefixes_len_change;
            }

            // finally, record the length (possibly updated by the translator)
            *ENT (uint16_t, ctx->con_len, word_index) = (uint16_t)(st_size + prefixes_len);
        }
    }

    return container_reconstruct_do (vb, ctx, con_p, prefixes, prefixes_len); 
}

// display
void container_display (ConstContainerP con)
{
    uint32_t num_items = con_nitems (*con);

    fprintf (info_stream, "repeats=%u\nnum_items=%u\ndrop_final_item_sep=%u\ndrop_final_repeat_sep=%u\n"
                          "filter_repeats=%u\nfilter_items=%u\nis_toplevel=%u\nrepsep={ %u %u }\n",
             con->repeats, num_items, con->drop_final_item_sep, con->drop_final_repeat_sep, 
             con->filter_repeats, con->filter_items, con->is_toplevel, con->repsep[0], con->repsep[1]);
    
    for (unsigned i=0; i < num_items; i++)
        fprintf (info_stream, "item %u: dict_id=%s seperator={ %u %u } translator=%u\n",
                 i, dis_dict_id (con->items[i].dict_id).s,  
                 con->items[i].seperator[0], con->items[i].seperator[1], con->items[i].translator);
    
    fflush (info_stream);
}

// Translators reconstructing last_value as a little endian binary
#define SET_n(type,mn,mx) type n = (type)ctx->last_value.i; \
                           ASSINP (ctx->last_value.i>=(int64_t)(mn) && ctx->last_value.i<=(int64_t)(mx), "Error: Failed to convert %s to %s because of bad data in line %u of the %s file: value of %s=%"PRId64" is out of range [%"PRId64"-%"PRId64"]",\
                                   dt_name (z_file->data_type), dt_name (flag.out_dt), vb->line_i, dt_name (z_file->data_type), \
                                   ctx->name, (int64_t)(ctx->last_value.i), (int64_t)(mn), (int64_t)(mx))

TRANSLATOR_FUNC (container_translate_I8)       { SET_n (int8_t,   INT8_MIN,  INT8_MAX  );              RECONSTRUCT1(n);       return 0; }
TRANSLATOR_FUNC (container_translate_U8)       { SET_n (uint8_t,  0,         UINT8_MAX );              RECONSTRUCT1((char)n); return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_I16) { SET_n (int16_t,  INT16_MIN, INT16_MAX ); n=LTEN16(n); RECONSTRUCT(&n,2);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_U16) { SET_n (uint16_t, 0,         UINT16_MAX); n=LTEN16(n); RECONSTRUCT(&n,2);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_I32) { SET_n (int32_t,  INT32_MIN, INT32_MAX ); n=LTEN32(n); RECONSTRUCT(&n,4);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_U32) { SET_n (uint32_t, 0,         UINT32_MAX); n=LTEN32(n); RECONSTRUCT(&n,4);     return 0; }
