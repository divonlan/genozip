
// ------------------------------------------------------------------
//   container.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
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

void container_prepare_snip (ContainerP con, const char *prefixes, unsigned prefixes_len, 
                             char *snip, unsigned *snip_len) // in / out
{
    ASSERT (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "prefixes_len=%u is beyond maximum of %u",
            prefixes_len, CONTAINER_MAX_PREFIXES_LEN);

    // make sure we have enough memory
    unsigned con_size = con_sizeof (*con);
    unsigned snip_len_needed = 1 + base64_size(con_size) + prefixes_len;
    ASSERT (*snip_len >= snip_len_needed, "snip_len=%u too short, need at least %u", *snip_len, snip_len_needed);

    con->repeats = BGEN24 (con->repeats); 
    snip[0] = SNIP_CONTAINER;
    
    unsigned b64_len = base64_encode ((uint8_t*)con, con_size, &snip[1]);
    ASSERT (b64_len <= base64_size(con_size), "b64_len=%u larger than base64_size(%u)=%u",
            b64_len, con_size, base64_size(con_size));
    
    con->repeats = BGEN24 (con->repeats); // restore

    if (prefixes_len) memcpy (&snip[1+b64_len], prefixes, prefixes_len);
    *snip_len = 1 + b64_len + prefixes_len;
}

WordIndex container_seg_by_ctx (VBlock *vb, Context *ctx, ContainerP con, 
                                // prefixes may be NULL or may be contain a container-wide prefix and per-item prefixes.
                                // The "container-wide prefix" is reconstructed once, at the beginning of the Container.
                                // The per-item prefixes are reconstructed once per repeat, before their respective items. 
                                // - prefixes (if not NULL) starts with a CON_PREFIX_SEP
                                // - then the container-wide prefix terminated by CON_PREFIX_SEP 
                                // - then one prefix per item terminated by CON_PREFIX_SEP
                                // all prefixes may be empty; empty prefixes at the end of the prefix string maybe omitted.
                                const char *prefixes, unsigned prefixes_len, 
                                unsigned add_bytes)
{
    ctx->no_stons = true; // we need the word index to for container caching

    // con=NULL means MISSING Container (see container_reconstruct_do)
    if (!con) return seg_by_ctx (vb, NULL, 0, ctx, 0); 

    unsigned snip_len = 1 + base64_size(con_sizeof (*con)) + prefixes_len;
    char snip[snip_len]; 
    container_prepare_snip (con, prefixes, prefixes_len, snip, &snip_len);

    if (flag.show_containers) { 
        iprintf ("VB=%u Line=%u Ctx=%u:%s Repeats=%u RepSep=%u,%u Items=", 
                 vb->vblock_i, vb->line_i, ctx->did_i, ctx->name, con->repeats, con->repsep[0], con->repsep[1]);
        for (unsigned i=0; i < con_nitems (*con); i++) {
            ContainerItem *item = &con->items[i];
            if (item->dict_id.num) 
                iprintf ("%u:%s ", ctx->did_i, dis_dict_id (item->dict_id).s); 
        }
        iprint0 ("\n");
    }

    // store in struct cache
    return seg_by_ctx (vb, snip, snip_len, ctx, add_bytes); 
}

//----------------------
// Reconstruction
//----------------------

static inline void container_reconstruct_prefix (VBlockP vb, ConstContainerP con, const char **prefixes, uint32_t *prefixes_len,
                                                 bool skip, DictId con_dict_id, DictId item_dict_id)
{
    ASSERT (*prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "prefixes_len=%u is too big", *prefixes_len);

    if (! (*prefixes_len)) return; // nothing to do
    
    const char *start = *prefixes;
    while (**prefixes != CON_PREFIX_SEP && **prefixes != CON_PREFIX_SEP_SHOW_REPEATS) (*prefixes)++; // prefixes are terminated by CON_PREFIX_SEP
    uint32_t len = (unsigned)((*prefixes) - start);

    if (!skip) {
        if (len) {
            RECONSTRUCT (start, len);

            if (flag.show_containers)
                iprintf ("Prefix(%s,%s)=\"%.*s\"\n", dis_dict_id (con_dict_id).s, dis_dict_id (item_dict_id).s, len, start);
        }

        // if the seperator is CON_PREFIX_SEP_SHOW_REPEATS, output the number of repeats. This is for BAM 'B' array 'count' field.
        if (**prefixes == CON_PREFIX_SEP_SHOW_REPEATS)
            RECONSTRUCT_BIN32 (con->repeats);
    }

    (*prefixes)++; // skip seperator
    (*prefixes_len) -= len + 1;
}


CONTAINER_FILTER_FUNC (container_no_filter)
{
    ABORT ("File %s requires a filter for dict_id=%s item=%u. Please upgrade to the latest version of genozip",
           z_name, dis_dict_id (dict_id).s, item);

    return true;    
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
    container_reconstruct_prefix (vb, con, &prefixes, &prefixes_len, false, ctx->dict_id, DICT_ID_NONE); 

    ASSERT (DTP (container_filter) || (!con->filter_repeats && !con->filter_items), 
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
            vb->line_i = vb->first_line + rep_i; // 1-based line from the begginging for the file, including the header
            vb->line_start = vb->txt_data.len;

            vb->dont_show_curr_line = false; // initialize for this line
        }

        if (con->filter_repeats && !(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, -1, NULL))) 
            continue; // repeat is filtered out

        char *rep_reconstruction_start = AFTERENT (char, vb->txt_data);

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t item_prefixes_len = prefixes_len;

        last_non_filtered_item_i = -1;
        for (unsigned i=0; i < num_items; i++) {
            const ContainerItem *item = &con->items[i];
            bool reconstruct = true;

            // an item filter may filter items in two ways:
            // - returns true - item is filter out, and data is not consumed
            // - sets reconstruct=false - data is consumed, but item is not reconstructed
            if (con->filter_items && 
                (!(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, i, &reconstruct) || !reconstruct))) {
                
                container_reconstruct_prefix (vb, con, &item_prefixes, &item_prefixes_len, true, DICT_ID_NONE, DICT_ID_NONE); // skip prefix            
                continue; // item is filtered out
            }

            last_non_filtered_item_i = i;

            if (flag.show_containers && item_ctxs[i]) // show container reconstruction 
                iprintf ("VB=%u Line=%u Repeat=%u %s->%s txt_data.len=%"PRIu64" (0x%04"PRIx64") (BEFORE)\n", 
                         vb->vblock_i, vb->line_i, rep_i, dis_dict_id (ctx->dict_id).s, item_ctxs[i]->name, 
                         vb->vb_position_txt_file + vb->txt_data.len, vb->vb_position_txt_file + vb->txt_data.len);

/*BRKPOINT*/container_reconstruct_prefix (vb, con, &item_prefixes, &item_prefixes_len, false, ctx->dict_id, con->items[i].dict_id); // item prefix (we will have one per item or none at all)

            int32_t reconstructed_len=0;
            if (item->dict_id.num) {  // not a prefix-only or translator-only item
                char *reconstruction_start = AFTERENT (char, vb->txt_data);
                reconstruct = reconstruct && !flag.show_sex && !flag.show_coverage && !flag.idxstats && // no reconstruction with these flags
                                  (  !flag.trans_containers ||   // not translating OR... 
                                     !IS_CI_SET (CI_TRANS_NOR)); // no prohibition on reconstructing when translating

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
        } // items loop

        // remove final seperator, if we need to (introduced v12)
        if (con->drop_final_item_sep && last_non_filtered_item_i >= 0) {
            const ContainerItem *item = &con->items[last_non_filtered_item_i]; // last_non_filtered_item_i is the last item that survived the filter, of the last repeat

            vb->txt_data.len -= CI_ITEM_HAS_FLAG(item) ? (flag.trans_containers ? 0 : !!IS_CI_SET (CI_NATIVE_NEXT))
                                                    : (!!item->seperator[0] + !!item->seperator[1]);
        }

        // reconstruct repeats separator, if neeeded
        if (rep_i+1 < con->repeats || !con->drop_final_repeat_sep) {
            if (con->repsep[0]) RECONSTRUCT1 (con->repsep[0]);
            if (con->repsep[1]) RECONSTRUCT1 (con->repsep[1]);
        }

        // call callback if needed now that repeat reconstruction is done
        if (con->callback) 
            DT_FUNC(vb, container_cb)(vb, ctx->dict_id, rep_i, rep_reconstruction_start, AFTERENT (char, vb->txt_data) - rep_reconstruction_start);

        // in top level: after consuming the line's data, if it is not to be outputted - drop it
        if (con->is_toplevel && vb->dont_show_curr_line) {
            ASSERT0 (flag.may_drop_lines, "Lines cannot be dropped because flag.may_drop_lines=false");
            vb->txt_data.len = vb->line_start;
            
            // only add "dropped line" (\b\n) for sorter, if its a beginning of a textual line (in case of --sequential, it wont
            // be for lines 2+ of the FASTA sequence)
            if (!vb->txt_data.len || *LASTENT (char, vb->txt_data) == '\n') {
                static const char *drop_lines[] = { "", "\b\n", "\b\n\b\n", "\b\n\b\n\b\n", "\b\n\b\n\b\n\b\n" };
                unsigned ht = DTPT (line_height);
                ASSERT (ht < sizeof (drop_lines) / sizeof (drop_lines[0]), "ht=%u is too large", ht);
                RECONSTRUCT (drop_lines[ht], ht*2); // tell sorter_piz_writer to drop line
            }
        }
    } // repeats loop

    // remove final seperator, only of final item (since v12, only remained used (by mistake) by SAM arrays - TODO: obsolete this. see sam_seg_array_field)
    if (con->drop_final_item_sep_of_final_repeat && last_non_filtered_item_i >= 0) {
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

        ASSERT (con_nitems (con) <= MAX_SUBFIELDS, "A container of %s has %u items which is beyond MAX_SUBFIELDS=%u. Please upgrade to latest version of genozip to access this file.",
                ctx->name, con_nitems (con), MAX_SUBFIELDS);

        // get the did_i for each dict_id - unfortunately we can only store did_i up to 254 (changing this would be a change in the file format)
        for (uint32_t item_i=0; item_i < con_nitems (con); item_i++)
            if (con.items[item_i].dict_id.num) { // not a prefix-only item
                DidIType did_i = ctx_get_existing_did_i (vb, con.items[item_i].dict_id);
                ASSERT (did_i != DID_I_NONE, "analyzing a %s container: unable to find did_i for item %s",
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

        ASSERT (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "ctx=%s: prefixes_len=%u longer than CONTAINER_MAX_PREFIXES_LEN=%u", 
                ctx->name, prefixes_len, CONTAINER_MAX_PREFIXES_LEN);

        ASSERT (!prefixes_len || prefixes[prefixes_len-1] == CON_PREFIX_SEP, 
                "ctx=%s: prefixes array does end with a CON_PREFIX_SEP: %.*s", ctx->name, prefixes_len, prefixes);

        // condition testing needed only for backward compatability with 8.0.x: prior to version 8.1 containers didn't have NO_STONS 
        // so singletons were in ctx.local without a word_index)
        if (word_index != WORD_INDEX_NONE) {

            ASSERT (st_size + prefixes_len <= 65535, "ctx=%s: st_size=%u + prefixes_len=%u too large", 
                    ctx->name, st_size, prefixes_len);

            // first encounter with Container for this context - allocate the cache
            if (!cache_exists) {
                buf_alloc (vb, &ctx->con_index, 0, ctx->word_list.len, uint32_t, 1, "context->con_index");
                buf_alloc (vb, &ctx->con_len,   0, ctx->word_list.len, uint16_t, 1, "context->con_len");
                buf_zero (&ctx->con_len); // zero entire buffer, not just part added
            }

            // place Container followed by prefix in the cache
            *ENT (uint32_t, ctx->con_index, word_index) = (uint32_t)ctx->con_cache.len;

            buf_alloc_old (vb, &ctx->con_cache, ctx->con_cache.len + st_size + prefixes_len + CONTAINER_MAX_SELF_TRANS_CHANGE, 2, "context->con_cache");
            
            char *cached_con = AFTERENT (char, ctx->con_cache);
            buf_add (&ctx->con_cache, con_p, st_size);
            if (prefixes_len) buf_add (&ctx->con_cache, prefixes, prefixes_len);

            con_p    = (ContainerP)cached_con; // update so we reconstruct from translated cached item
            prefixes = cached_con + st_size;

            // if item[0] is a translator-only item, use it to translate the Container itself (used by SAM_OPTIONAL)
            ContainerItem *item0 = &con.items[0];
            if (flag.trans_containers && !item0->dict_id.num && item0->translator) {
                int32_t prefixes_len_change = DT_FUNC(vb, translator)[item0->translator](vb, ctx, cached_con, st_size + prefixes_len);  
                ASSERT (prefixes_len_change <= CONTAINER_MAX_SELF_TRANS_CHANGE, 
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

    iprintf ("repeats=%u\nnum_items=%u\ndrop_final_item_sep=%u\ndrop_final_repeat_sep=%u\n"
                          "filter_repeats=%u\nfilter_items=%u\nis_toplevel=%u\nrepsep={ %u %u }\n",
             con->repeats, num_items, con->drop_final_item_sep, con->drop_final_repeat_sep, 
             con->filter_repeats, con->filter_items, con->is_toplevel, con->repsep[0], con->repsep[1]);
    
    for (unsigned i=0; i < num_items; i++)
        iprintf ("item %u: dict_id=%s seperator={ %u %u } translator=%u\n",
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
