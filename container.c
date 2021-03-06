
// ------------------------------------------------------------------
//   container.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#if defined __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif
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
#include "regions.h"
#include "piz.h"

//----------------------
// Segmentation
//----------------------

void container_prepare_snip (ConstContainerP con, const char *prefixes, unsigned prefixes_len, 
                             char *snip, unsigned *snip_len) // in / out
{
    ASSERT (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "prefixes_len=%u is beyond maximum of %u",
            prefixes_len, CONTAINER_MAX_PREFIXES_LEN);

    // make sure we have enough memory
    unsigned con_size = con_sizeof (*con);
    unsigned snip_len_needed = 1 + base64_size(con_size) + prefixes_len;
    ASSERT (*snip_len >= snip_len_needed, "snip_len=%u too short, need at least %u", *snip_len, snip_len_needed);

    ((Container*)con)->repeats = BGEN24 (con->repeats); 
    snip[0] = SNIP_CONTAINER;
    
    unsigned b64_len = base64_encode ((uint8_t*)con, con_size, &snip[1]);
    ASSERT (b64_len <= base64_size(con_size), "b64_len=%u larger than base64_size(%u)=%u",
            b64_len, con_size, base64_size(con_size));
    
    ((Container*)con)->repeats = BGEN24 (con->repeats); // restore - honoring our "const" contract

    if (prefixes_len) memcpy (&snip[1+b64_len], prefixes, prefixes_len);
    *snip_len = 1 + b64_len + prefixes_len;
}

WordIndex container_seg_by_ctx_do (VBlock *vb, Context *ctx, ConstContainerP con, 
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
        iprintf ("VB=%u Line=%"PRIu64" Ctx=%u:%s Repeats=%u RepSep=%u,%u Items=", 
                 vb->vblock_i, vb->line_i, ctx->did_i, ctx->name, con->repeats, con->repsep[0], con->repsep[1]);
        for (unsigned i=0; i < con_nitems (*con); i++) {
            const ContainerItem *item = &con->items[i];
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
                     : con->items[i].did_i_small < 255 ? CTX(con->items[i].did_i_small)
                     :                                   ctx_get_existing_ctx (vb, con->items[i].dict_id);

    // for containers, new_value is the some of all its items, all repeats last_value (either int or float)
    LastValueType new_value = {};

    if (flag.show_containers) // show container reconstruction 
        iprintf ("VB=%u Container: %s repeats=%u items=%u\n", 
                 vb->vblock_i, dis_dict_id (ctx->dict_id).s, con->repeats, con_nitems(*con));

    if (con->is_toplevel) 
        buf_alloc (vb, &vb->lines, 0, con->repeats+1, char *, 1.1, "lines"); // note: lines.len was set in piz_read_one_vb
    
    for (uint32_t rep_i=0; rep_i < con->repeats; rep_i++) {

        // case this is the top-level snip
        if (con->is_toplevel) {
            vb->line_i = vb->first_line + rep_i; // 1-based line from the begginging for the file, including the header
            
            *ENT (char *, vb->lines, rep_i) = AFTERENT (char, vb->txt_data); 
            vb->line_start = vb->txt_data.len;

            vb->drop_curr_line = NULL; // initialize for this line
        }

        if (con->filter_repeats && !(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, -1, NULL))) 
            continue; // repeat is filtered out

        char *rep_reconstruction_start = AFTERENT (char, vb->txt_data);

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t item_prefixes_len = prefixes_len;

        last_non_filtered_item_i = -1;
        unsigned num_preceding_seps = 0;
        for (unsigned i=0; i < num_items; i++) {
            const ContainerItem *item = &con->items[i];
            Context *item_ctx = item_ctxs[i];
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

            if (flag.show_containers && item_ctx) // show container reconstruction 
                iprintf ("VB=%u Line=%"PRIu64" Repeat=%u %s->%s trans_id=%u txt_data.len=%"PRIu64" (0x%04"PRIx64") (BEFORE)\n", 
                         vb->vblock_i, vb->line_i, rep_i, dis_dict_id (ctx->dict_id).s, item_ctx->name,
                         vb->translation.trans_containers ? item->translator : 0, 
                         vb->vb_position_txt_file + vb->txt_data.len, vb->vb_position_txt_file + vb->txt_data.len);

/*BRKPOINT*/container_reconstruct_prefix (vb, con, &item_prefixes, &item_prefixes_len, false, ctx->dict_id, con->items[i].dict_id); // item prefix (we will have one per item or none at all)

            int32_t recon_len=0;
            if (item->dict_id.num) {  // not a prefix-only or translator-only item
                char *reconstruction_start = AFTERENT (char, vb->txt_data);
                reconstruct = reconstruct && !flag.collect_coverage && // no reconstruction with these flags
                                  (  !vb->translation.trans_containers ||   // not translating OR... 
                                     !IS_CI_SET (CI_TRANS_NOR)); // no prohibition on reconstructing when translating

                recon_len = reconstruct_from_ctx (vb, item_ctx->did_i, 0, reconstruct); // -1 if WORD_INDEX_MISSING_SF

                // sum up items' values if needed
                if (ctx->flags.store == STORE_INT)
                    new_value.i += item_ctx->last_value.i;
                
                else if (ctx->flags.store == STORE_FLOAT)
                    new_value.f += item_ctx->last_value.f;

                // if we're reconstructing to a translated format (eg SAM2BAM) - re-reconstruct this item
                // using the designated "translator" function, if one is available
                if (vb->translation.trans_containers && item->translator && recon_len != -1 &&
                    !(flag.missing_contexts_allowed && !item_ctx->dict.len && !item_ctx->local.len)) // skip if missing contexts are allowed, and this context it missing 
                    DT_FUNC(vb, translator)[item->translator](vb, item_ctx, reconstruction_start, recon_len, false);  
            }            

            // case: WORD_INDEX_MISSING_SF - delete previous item's separator if it has one (used by SAM_OPTIONAL - sam_seg_optional_all)
            if (recon_len == -1 && i > 0 && !con->keep_empty_item_sep && !CI_ITEM_HAS_FLAG(item-1))  
                vb->txt_data.len -= num_preceding_seps;

            // reconstruct separator(s) as needed
            if (reconstruct) {
                if (vb->translation.trans_containers && IS_CI_SET (CI_TRANS_NUL)) {
                    RECONSTRUCT1 (0);
                    num_preceding_seps = 1;
                }

                else if (!vb->translation.trans_containers && IS_CI_SET (CI_NATIVE_NEXT)) {
                    RECONSTRUCT1 (item->seperator[1]);
                    num_preceding_seps = 1;
                }

                else if (!CI_ITEM_HAS_FLAG(item)) { // seperator, not flag
                    if (item->seperator[0]) RECONSTRUCT1 ((char)item->seperator[0]);
                    if (item->seperator[1]) RECONSTRUCT1 ((char)item->seperator[1]);
                    num_preceding_seps = (item->seperator[0] != 0) + (item->seperator[1] != 0);
                }

                else
                    num_preceding_seps = 0;
            }

            // after all reconstruction and translation is done - move if needed
            if (vb->translation.trans_containers && IS_CI_SET (CI_TRANS_MOVE))
                vb->txt_data.len += (uint8_t)item->seperator[1];

            // test for exclusion of the line due to --regions
            if (flag.regions && con->is_toplevel && item_ctx->did_i == DTF(test_regions) && 
                !regions_is_site_included (vb))
                vb->drop_curr_line = "regions";
        } // items loop

        // remove final seperator, if we need to (introduced v12)
        if (con->drop_final_item_sep && last_non_filtered_item_i >= 0) {
            const ContainerItem *item = &con->items[last_non_filtered_item_i]; // last_non_filtered_item_i is the last item that survived the filter, of the last repeat

            vb->txt_data.len -= CI_ITEM_HAS_FLAG(item) ? (vb->translation.trans_containers ? 0 : !!IS_CI_SET (CI_NATIVE_NEXT))
                                                    : (!!item->seperator[0] + !!item->seperator[1]);
        }

        // reconstruct repeats separator, if neeeded
        if (rep_i+1 < con->repeats || !con->drop_final_repeat_sep) {
            if (con->repsep[0]) RECONSTRUCT1 (con->repsep[0]);
            if (con->repsep[1]) RECONSTRUCT1 (con->repsep[1]);
        }

        // call callback if needed now that repeat reconstruction is done (always callback for top level)
        if (con->callback || (con->is_toplevel && DTP (container_cb)))
            DT_FUNC(vb, container_cb)(vb, ctx->dict_id, con->is_toplevel, rep_i, con->repeats, rep_reconstruction_start, 
                    AFTERENT (char, vb->txt_data) - rep_reconstruction_start);

        // in top level: after consuming the line's data, if it is not to be outputted - drop it
        if (con->is_toplevel) {

            // filter by --lines
            if (flag.lines_first >= 0) {
                int64_t recon_line_i = vb->first_line + (int64_t)rep_i - 1LL; // 0-based (note: can't use vb->line_i as it is textual lines rather than reps)
                if (recon_line_i < flag.lines_first) 
                    vb->drop_curr_line = "lines";

                else if (recon_line_i > flag.lines_last) {
                    vb->drop_curr_line = "lines";
                    vb->txt_data.len = vb->line_start;
                    
                     // skip to the end - no need to reconstruct any further lines
                    for (rep_i=rep_i+1; rep_i < con->repeats; rep_i++)
                        *ENT (char *, vb->lines, rep_i) = AFTERENT (char, vb->txt_data); 
                }
            }

            // filter by --tail
            if (!vb->drop_curr_line && flag.tail && vb->vblock_i == txt_file->tail_1st_vb &&
                rep_i < txt_file->tail_1st_line_1st_vb)
                vb->drop_curr_line = "tail";

            // filter by --grep (but not for FASTQ or FASTA - they implement their own logic)
            if (!vb->drop_curr_line && flag.grep && txt_file->data_type != DT_FASTA && txt_file->data_type != DT_FASTQ
                && !piz_grep_match (rep_reconstruction_start, AFTERENT (char, vb->txt_data)))
                vb->drop_curr_line = "grep";
                /*
                //SAFE_NUL (&rep_reconstruction_start[AFTERENT (char, vb->txt_data) - rep_reconstruction_start]);
                SAFE_NUL (AFTERENT (char, vb->txt_data));

                if (!flag.grepw && strstr (rep_reconstruction_start, flag.grep))
                    vb->drop_curr_line = "grep";

                // case: --grepw - grep whole word
                else if (flag.grepw) { 
                    const char *s = rep_reconstruction_start;
                    while (s <= AFTERENT (char, vb->txt_data) - flag.grep_len) {
                        if (!(s = strstr (s, flag.grep))) break;

                        char before = (s == rep_reconstruction_start ? ' ' : s[-1]);
                        char after  = s[flag.grep_len];
                    
                        if (!IS_LETTER(before) && !IS_DIGIT(before) &&
                            !IS_LETTER(after) && !IS_DIGIT(after)) {
                            
                            vb->drop_curr_line = "grep";
                            break;
                        }

                        s += flag.grep_len;
                    }
                }
                            
                SAFE_RESTORE;
            }*/

            if (vb->drop_curr_line) {
                ASSERT (flag.maybe_vb_modified_by_reconstructor, "Attempting drop_curr_line=\"%s\", but lines cannot be dropped because flag.maybe_vb_modified_by_reconstructor=false. This is bug in the code. vb_i=%u line_i=%"PRIu64, 
                        vb->drop_curr_line, vb->vblock_i, vb->line_i);

                vb->txt_data.len = vb->line_start;
            }
            else
                vb->num_nondrop_lines++;

            if (flag.show_containers && vb->drop_curr_line) // show container reconstruction 
                iprintf ("VB=%u Line=%"PRIu64" dropped due to \"%s\"\n", vb->vblock_i, vb->line_i, vb->drop_curr_line);
        }
    } // repeats loop

    if (con->is_toplevel) {
        // final line_start entry (it has a total of num_lines+1 entires) - done last, after possibly dropping line
        *AFTERENT (char *, vb->lines) = AFTERENT (char, vb->txt_data); // we allocated vb->lines.len+1 entries in lines

        // sanity checks
        ASSERT (vb->lines.len == con->repeats, "Expected vb->lines.len=%"PRIu64" == con->repeats=%u", vb->lines.len, con->repeats);
        ASSERT (flag.data_modified || flag_loading_auxiliary || vb->num_nondrop_lines == vb->recon_num_lines, 
                "Expected vb->num_nondrop_lines=%u == vb->recon_num_lines=%u", vb->num_nondrop_lines, vb->recon_num_lines);
    }

    // remove final seperator, only of final item (since v12, only remained used (by mistake) by SAM arrays - TODO: obsolete this. see sam_seg_array_field)
    if (con->drop_final_item_sep_of_final_repeat && last_non_filtered_item_i >= 0) {
        const ContainerItem *item = &con->items[last_non_filtered_item_i]; // last_non_filtered_item_i is the last item that survived the filter, of the last repeat

        vb->txt_data.len -= CI_ITEM_HAS_FLAG(item) ? (vb->translation.trans_containers ? 0 : !!IS_CI_SET (CI_NATIVE_NEXT))
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

    bool cache_exists = buf_is_alloc (&ctx->con_cache);
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

        ASSERT (con_nitems (con) <= MAX_FIELDS, "A container of %s has %u items which is beyond MAX_FIELDS=%u. Please upgrade to latest version of genozip to access this file.",
                ctx->name, con_nitems (con), MAX_FIELDS);

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

            buf_alloc (vb, &ctx->con_cache, st_size + prefixes_len + CONTAINER_MAX_SELF_TRANS_CHANGE, 0, char, 2, "context->con_cache");
            
            char *cached_con = AFTERENT (char, ctx->con_cache);
            buf_add (&ctx->con_cache, con_p, st_size);
            if (prefixes_len) buf_add (&ctx->con_cache, prefixes, prefixes_len);

            con_p    = (ContainerP)cached_con; // update so we reconstruct from translated cached item
            prefixes = cached_con + st_size;

            // if item[0] is a translator-only item, use it to translate the Container itself (used by SAM_OPTIONAL)
            ContainerItem *item0 = &con.items[0];
            if (vb->translation.trans_containers && !item0->dict_id.num && item0->translator) {
                int32_t prefixes_len_change = (DT_FUNC(vb, translator)[item0->translator](vb, ctx, cached_con, st_size + prefixes_len, false));  
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
                           ASSINP (ctx->last_value.i>=(int64_t)(mn) && ctx->last_value.i<=(int64_t)(mx), "Error: Failed to convert %s to %s because of bad data in line %"PRIu64" of the %s file: value of %s=%"PRId64" is out of range [%"PRId64"-%"PRId64"]",\
                                   dt_name (z_file->data_type), dt_name (flag.out_dt), vb->line_i, dt_name (z_file->data_type), \
                                   ctx->name, (int64_t)(ctx->last_value.i), (int64_t)(mn), (int64_t)(mx))

TRANSLATOR_FUNC (container_translate_I8)       { SET_n (int8_t,   INT8_MIN,  INT8_MAX  );              RECONSTRUCT1(n);       return 0; }
TRANSLATOR_FUNC (container_translate_U8)       { SET_n (uint8_t,  0,         UINT8_MAX );              RECONSTRUCT1((char)n); return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_I16) { SET_n (int16_t,  INT16_MIN, INT16_MAX ); n=LTEN16(n); RECONSTRUCT(&n,2);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_U16) { SET_n (uint16_t, 0,         UINT16_MAX); n=LTEN16(n); RECONSTRUCT(&n,2);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_I32) { SET_n (int32_t,  INT32_MIN, INT32_MAX ); n=LTEN32(n); RECONSTRUCT(&n,4);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_U32) { SET_n (uint32_t, 0,         UINT32_MAX); n=LTEN32(n); RECONSTRUCT(&n,4);     return 0; }
