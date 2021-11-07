
// ------------------------------------------------------------------
//   container.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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

void container_prepare_snip (ConstContainerP con, STRp(prefixes), 
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

// WARNING: don't pass in con an address of a static container, it is not thread safe! container_prepare_snip temporarily BGEN fields of the container.
WordIndex container_seg_do (VBlock *vb, Context *ctx, ConstContainerP con, 
                            // prefixes may be NULL or may be contain a container-wide prefix and per-item prefixes.
                            // The "container-wide prefix" is reconstructed once, at the beginning of the Container.
                            // The per-item prefixes are reconstructed once per repeat, before their respective items. 
                            // - prefixes (if not NULL) starts with a CON_PX_SEP
                            // - then the container-wide prefix terminated by CON_PX_SEP 
                            // - then one prefix per item terminated by CON_PX_SEP
                            // all prefixes may be empty; empty prefixes at the end of the prefix string maybe omitted.
                            const char *prefixes, unsigned prefixes_len, 
                            const char *ren_prefixes, unsigned ren_prefixes_len, 
                            unsigned add_bytes,
                            bool *is_new) // optional out
{
    ctx->no_stons = true; // we need the word index to for container caching

    // con=NULL means MISSING Container (see container_reconstruct_do)
    if (!con) return seg_by_ctx (VB, NULL, 0, ctx, 0); 
    
    unsigned con_b64_size = base64_size(con_sizeof (*con));
    unsigned con1_snip_len = 1 + con_b64_size + prefixes_len;
    unsigned con2_snip_len = 1 + con_b64_size + ren_prefixes_len;

    unsigned snip_len = ren_prefixes_len ? 2 + con1_snip_len + con2_snip_len : con1_snip_len;
    char snip[snip_len];

    // case: we have renamed prefixes due to dual-coordinates with tag renaming
    if (ren_prefixes_len) {
        snip[0] = SNIP_DUAL;
        container_prepare_snip (con, prefixes, prefixes_len, &snip[1], &con1_snip_len);
        snip[1 + con1_snip_len] = SNIP_DUAL;
        container_prepare_snip (con, ren_prefixes, ren_prefixes_len, &snip[2+con1_snip_len], &con2_snip_len);
        snip_len = 2 + con1_snip_len + con2_snip_len;
    }
    else 
        container_prepare_snip (con, prefixes, prefixes_len, snip, &snip_len);

    if (flag.show_containers) { 
        iprintf ("VB=%u Line=%"PRIu64" Ctx=%u:%s Repeats=%u RepSep=%u,%u Items=", 
                 vb->vblock_i, vb->line_i, ctx->did_i, ctx->tag_name, con->repeats, con->repsep[0], con->repsep[1]);
        for (unsigned i=0; i < con_nitems (*con); i++) {
            const ContainerItem *item = &con->items[i];
            if (item->dict_id.num) 
                iprintf ("%u:%s ", ctx->did_i, dis_dict_id (item->dict_id).s); 
        }
        iprint0 ("\n");
    }

    return seg_by_ctx_ex (vb, snip, snip_len, ctx, add_bytes, is_new); 
}

//----------------------
// Reconstruction
//----------------------

static inline uint32_t container_reconstruct_prefix (VBlockP vb, ConstContainerP con, pSTRp(prefixes),
                                                     bool skip, DictId con_dict_id, DictId item_dict_id)
{
    ASSERT (*prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "prefixes_len=%u is too big", *prefixes_len);

    if (! (*prefixes_len)) return 0; // nothing to do
    
    const char *start = *prefixes;
    while (**prefixes != CON_PX_SEP && **prefixes != CON_PX_SEP_SHOW_REPEATS) (*prefixes)++; // prefixes are terminated by CON_PX_SEP
    uint32_t len = (unsigned)((*prefixes) - start);

    if (!skip) {
        if (len) {
            RECONSTRUCT (start, len);

            if (flag.show_containers)
                iprintf ("Prefix(%s,%s)=\"%.*s\"\n", dis_dict_id (con_dict_id).s, dis_dict_id (item_dict_id).s, len, start);
        }

        // if the separator is CON_PX_SEP_SHOW_REPEATS, output the number of repeats. This is for BAM 'B' array 'count' field.
        if (**prefixes == CON_PX_SEP_SHOW_REPEATS)
            RECONSTRUCT_BIN32 (con->repeats);
    }

    (*prefixes)++; // skip separator
    (*prefixes_len) -= len + 1;

    return len;
}


CONTAINER_FILTER_FUNC (container_no_filter)
{
    ABORT ("File %s requires a filter for dict_id=%s item=%u. Please upgrade to the latest version of genozip",
           z_name, dis_dict_id (dict_id).s, item);

    return true;    
}

#define IS_CI_SET(flag) (CI_ITEM_HAS_FLAG (item) && ((uint8_t)item->separator[0] & ~(uint8_t)0x80 & flag))

// reconstructs the item separator, returning the number of characters reconstructed
static inline unsigned container_reconstruct_item_seperator (VBlockP vb, const ContainerItem *item, char *reconstruction_start, bool translating)
{
    if (!item->separator[0])
        return 0;

    else if (translating && IS_CI_SET (CI_TRANS_NUL)) {
        RECONSTRUCT1 (0);
        return 1;
    }

    else if (!translating && IS_CI_SET (CI_NATIVE_NEXT)) {
        RECONSTRUCT1 (item->separator[1]);
        return 1;
    }

    else if (item->separator[0] == CI_FIXED_0_PAD) {
        int recon_len = AFTERENT (char, vb->txt_data) - reconstruction_start;
        int pad = (int)item->separator[1] - recon_len;
        if (pad > 0) {
            memmove (reconstruction_start + pad, reconstruction_start, recon_len);
            memset (reconstruction_start, '0', pad);
            vb->txt_data.len += pad;
        }
        return 0; // no seperators to delete if needed
    }

    else if (item->separator[0] < '\t') {
        ASSPIZ (false, "Unrecognized special seperator %u. Please upgrade to latest version of Genozip", item->separator[0]); // a seperator from the future
        return 0;
    }

    else if (!CI_ITEM_HAS_FLAG(item)) { // separator, not flag        unsigned sep_0_valid = item->separator[0] >= '\t'; // 1->8 reserved for special separators ; 0=no seperator
        unsigned sep_1_valid = item->separator[1] >= '\t';
        
        RECONSTRUCT1 ((char)item->separator[0]);
        if (sep_1_valid) RECONSTRUCT1 ((char)item->separator[1]);
        
        return 1 + sep_1_valid;
    }

    else
        return 0;
}

ValueType container_reconstruct (VBlockP vb, ContextP ctx, ConstContainerP con, STRp(prefixes))
{
    TimeSpecType profiler_timer = {}; 
    if (flag.show_time && con->is_toplevel) 
        clock_gettime (CLOCK_REALTIME, &profiler_timer);
    
    int32_t last_non_filtered_item_i = -1;

    bool translating = vb->translation.trans_containers && !con->no_translation;

    // container wide prefix - it will be missing if Container has no prefixes, or empty if it has only items prefixes
    container_reconstruct_prefix (vb, con, pSTRa(prefixes), false, ctx->dict_id, DICT_ID_NONE); 

    ASSERT (DTP (container_filter) || (!con->filter_repeats && !con->filter_items), 
            "data_type=%s doesn't support container_filter, despite being specified in the Container. Please upgrade to the latest version of Genozip.", 
            dt_name (vb->data_type));

    uint32_t num_items = con_nitems(*con);
    Context *item_ctxs[num_items];
    
    // we can cache did_i up to 254. dues historical reasons the field is only 8 bit. that's enough in most cases anyway.
    for (unsigned i=0; i < num_items; i++) 
        item_ctxs[i] = !con->items[i].dict_id.num      ? NULL 
                     : con->items[i].did_i_small < 255 ? CTX(con->items[i].did_i_small)
                     :                                   ECTX (con->items[i].dict_id);

    // for containers, new_value is the sum of all its items, all repeats last_value (either int or float)
    ValueType new_value = {};

    if (flag.show_containers) // show container reconstruction 
        iprintf ("VB=%u Container: %s repeats=%u items=%u filter_items=%u filter_repeats=%u callback=%u\n", 
                 vb->vblock_i, dis_dict_id (ctx->dict_id).s, con->repeats, con_nitems(*con), con->filter_items, con->filter_repeats, con->callback);

    if (con->is_toplevel) {
        buf_alloc (vb, &vb->lines, 0, con->repeats+1, char *, 1.1, "lines"); // note: lines.len was set in piz_read_one_vb
        buf_alloc_bitarr (vb, &vb->is_dropped, con->repeats, "is_dropped");
        buf_zero (&vb->is_dropped);
    }

    for (uint32_t rep_i=0; rep_i < con->repeats; rep_i++) {

        // case this is the top-level snip: initialize line
        if (con->is_toplevel) {
            vb->line_i         = vb->first_line + rep_i; // 1-based line from the begginging for the file, including the txt header            
            vb->sample_i       = 0;
            vb->line_start     = vb->txt_data.len;
            vb->drop_curr_line = NULL;    
            vb->buddy_line_i   = NO_BUDDY;
            *ENT (char *, vb->lines, rep_i) = AFTERENT (char, vb->txt_data); 
        }

        if (con->filter_repeats && !(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, -1, NULL))) 
            continue; // repeat is filtered out

        char *rep_reconstruction_start = AFTERENT (char, vb->txt_data);

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t remaining_prefix_len = prefixes_len; 

        last_non_filtered_item_i = -1;
        unsigned num_preceding_seps = 0;

        if (flag.show_containers) // show container reconstruction 
            iprintf ("VB=%u Line=%"PRIu64" Repeat=%u LastRepeat=%u %s\n", vb->vblock_i, vb->line_i, rep_i, con->repeats-1, ctx->tag_name);

        for (unsigned i=0; i < num_items; i++) {
            const ContainerItem *item = &con->items[i];
            Context *item_ctx = item_ctxs[i];
            bool reconstruct = !flag.genocat_no_reconstruct;
            bool trans_nor = translating && IS_CI_SET (CI_TRANS_NOR); // check for prohibition on reconstructing when translating

            // an item filter may filter items in two ways:
            // - returns true - item is filter out, and data is not consumed
            // - sets reconstruct=false - data is consumed, but item is not reconstructed
            if (con->filter_items && (!(DT_FUNC (vb, container_filter) (vb, ctx->dict_id, con, rep_i, i, &reconstruct)))) {
                
                container_reconstruct_prefix (vb, con, &item_prefixes, &remaining_prefix_len, true, DICT_ID_NONE, DICT_ID_NONE); // skip prefix            
                continue; // item is filtered out
            }

            last_non_filtered_item_i = i;

            if (flag.show_containers && item_ctx && (!flag.dict_id_show_containers.num || dict_id_typeless (item_ctx->dict_id).num == flag.dict_id_show_containers.num)) // show container reconstruction 
                iprintf ("VB=%u Line=%"PRIu64" Repeat=%u %s->%s trans_id=%u txt_data.len=%"PRIu64" (0x%04"PRIx64") reconstruct_prefix=%d reconstruct_value=%d%s", 
                         vb->vblock_i, vb->line_i, rep_i, ctx->tag_name, item_ctx->tag_name,
                         translating ? item->translator : 0, 
                         vb->vb_position_txt_file + vb->txt_data.len, vb->vb_position_txt_file + vb->txt_data.len,
                         reconstruct, reconstruct && !trans_nor, (flag.dict_id_show_containers.num ? " : " : "\n"));

/*BRKPOINT*/uint32_t item_prefix_len = 
                reconstruct ? container_reconstruct_prefix (vb, con, &item_prefixes, &remaining_prefix_len, false, ctx->dict_id, con->items[i].dict_id) : 0; // item prefix (we will have one per item or none at all)

            int32_t recon_len=0;
            char *reconstruction_start = AFTERENT (char, vb->txt_data);

            if (item->dict_id.num) {  // not a prefix-only or translator-only item
                reconstruct &= !trans_nor; // check for prohibition on reconstructing when translating
                reconstruct &= (item->separator[0] != CI_INVISIBLE); // check if this item should never be reconstructed

                recon_len = reconstruct_from_ctx (vb, item_ctx->did_i, 0, reconstruct); // -1 if WORD_INDEX_MISSING

                // sum up items' values if needed
                if (ctx->flags.store == STORE_INT)
                    new_value.i += item_ctx->last_value.i;
                
                else if (ctx->flags.store == STORE_FLOAT)
                    new_value.f += item_ctx->last_value.f;

                // case: reconstructing to a translated format (eg SAM2BAM) - modify the reconstruction ("translate") this item
                if (translating && item->translator && recon_len != -1 &&
                    !(flag.missing_contexts_allowed && !item_ctx->dict.len && !item_ctx->local.len)) // skip if missing contexts are allowed, and this context it missing 
                    DT_FUNC(vb, translator)[item->translator](vb, item_ctx, reconstruction_start, recon_len, item_prefix_len, false);  
            }            

            // case: WORD_INDEX_MISSING - delete previous item's separator if it has one (used by SAM_OPTIONAL - sam_seg_optional_all)
            if (recon_len == -1 && i > 0 && !con->keep_empty_item_sep && !CI_ITEM_HAS_FLAG(item-1))  
                vb->txt_data.len -= num_preceding_seps;

            // reconstruct separator(s) as needed
            if (reconstruct) 
                num_preceding_seps = container_reconstruct_item_seperator (vb, item, reconstruction_start, translating);

            // after all reconstruction and translation is done - move if needed
            if (translating && IS_CI_SET (CI_TRANS_MOVE))
                vb->txt_data.len += (uint8_t)item->separator[1];

            if (flag.show_containers && item_ctx && dict_id_typeless (item_ctx->dict_id).num == flag.dict_id_show_containers.num)
                iprintf ("\"%.*s\"\n", (int)(AFTERENT(char, vb->txt_data)-reconstruction_start), reconstruction_start);
                
        } // items loop

        // remove final separator, if we need to (introduced v12)
        if (con->drop_final_item_sep && last_non_filtered_item_i >= 0) {
            const ContainerItem *item = &con->items[last_non_filtered_item_i]; // last_non_filtered_item_i is the last item that survived the filter, of the last repeat

            vb->txt_data.len -= CI_ITEM_HAS_FLAG(item) ? (translating ? 0 : !!IS_CI_SET (CI_NATIVE_NEXT))
                                                       : (!!item->separator[0] + !!item->separator[1]);
        }

        // reconstruct repeats separator, if neeeded
        if (rep_i+1 < con->repeats || !con->drop_final_repeat_sep) {
            if (con->repsep[0]) RECONSTRUCT1 (con->repsep[0]);
            if (con->repsep[1]) RECONSTRUCT1 (con->repsep[1]);
        }

        // call callback if needed now that repeat reconstruction is done (always callback for top level)
        if (con->callback || (con->is_toplevel && DTP (container_cb)))
            DT_FUNC(vb, container_cb)(vb, ctx->dict_id, con->is_toplevel, rep_i, con->repeats, rep_reconstruction_start, 
                    AFTERENT (char, vb->txt_data) - rep_reconstruction_start, prefixes, prefixes_len);

        // in top level: after consuming the line's data, if it is not to be outputted - drop it
        if (con->is_toplevel) {

            // filter by --regions
            if (!vb->drop_curr_line && flag.regions && !regions_is_site_included (vb)) 
                vb->drop_curr_line = "regions";

            // filter by --lines
            if (flag.lines_first >= 0) {
                int64_t recon_line_i = vb->first_line + (int64_t)rep_i - 1LL; // 0-based (note: can't use vb->line_i as it is textual lines rather than reps)
                if (recon_line_i < flag.lines_first) 
                    vb->drop_curr_line = "lines";

                else if (recon_line_i > flag.lines_last) {
                    vb->drop_curr_line = "lines";
                    vb->txt_data.len = vb->line_start;
                    
                     // skip to the end - no need to reconstruct any further lines
                    for (rep_i=rep_i+1; rep_i < con->repeats; rep_i++) { 
                        bit_array_set (buf_get_bitarray (&vb->is_dropped), rep_i);
                        *ENT (char *, vb->lines, rep_i) = AFTERENT (char, vb->txt_data); 
                    }
                }
            } 

            // filter by --tail
            if (!vb->drop_curr_line && flag.tail && vb->vblock_i == txt_file->tail_1st_vb &&
                rep_i < txt_file->tail_1st_line_1st_vb)
                vb->drop_curr_line = "tail";

            // filter by --grep (but not for FASTQ or FASTA - they implement their own logic)
            if (!vb->drop_curr_line && flag.grep && txt_file->data_type != DT_FASTA /*&& txt_file->data_type != DT_FASTQ*/
                && !piz_grep_match (rep_reconstruction_start, AFTERENT (char, vb->txt_data)))
                vb->drop_curr_line = "grep";

            if (!vb->drop_curr_line && TXT_DT(DT_FASTQ)) {

                // in FASTQ --header-only - remove the 3 non-header lines only after --grep
                if (flag.header_only_fast) 
                    vb->txt_data.len = ENTNUM (vb->txt_data, strchr (rep_reconstruction_start, '\n') + 1);

                else if (flag.seq_only) {} // TO DO 

                else if (flag.qual_only) {} // TO DO 
            }

            if (vb->drop_curr_line) {
                ASSERT (flag.maybe_vb_modified_by_reconstructor, "Attempting drop_curr_line=\"%s\", but lines cannot be dropped because flag.maybe_vb_modified_by_reconstructor=false. This is bug in the code. vb_i=%u line_i=%"PRIu64, 
                        vb->drop_curr_line, vb->vblock_i, vb->line_i);

                // remove reconstructed text (to allow writer_flush_vb()), except if interleave, as we might un-drop the line in 
                if (!flag.interleave) vb->txt_data.len = vb->line_start;

                bit_array_set (buf_get_bitarray (&vb->is_dropped), rep_i);
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

    // remove final separator, only of final item (since v12, only remained used (by mistake) by SAM arrays - TODO: obsolete this. see sam_seg_array_field)
    if (con->drop_final_item_sep_of_final_repeat && last_non_filtered_item_i >= 0) {
        const ContainerItem *item = &con->items[last_non_filtered_item_i]; // last_non_filtered_item_i is the last item that survived the filter, of the last repeat

        vb->txt_data.len -= CI_ITEM_HAS_FLAG(item) ? (translating ? 0 : !!IS_CI_SET (CI_NATIVE_NEXT))
                                                   : (!!item->separator[0] + !!item->separator[1]);
    }
     
    if (con->is_toplevel)   
        COPY_TIMER (reconstruct_vb);

    return new_value;
}

ContainerP container_retrieve (VBlockP vb, ContextP ctx, WordIndex word_index, STRp(snip),
                               pSTRp(out_prefixes))
{
    Container con, *con_p=NULL;
    STR(prefixes);

    bool cache_exists = buf_is_alloc (&ctx->con_cache);
    uint16_t cache_item_len;

    // if this container exists in the cache - use the cached one
    if (cache_exists && word_index != WORD_INDEX_NONE && ((cache_item_len = *ENT (uint16_t, ctx->con_len, word_index)))) {
        con_p = (ContainerP)ENT (char, ctx->con_cache, *ENT (uint32_t, ctx->con_index, word_index));
        
        unsigned st_size = con_sizeof (*con_p);
        prefixes = (char *)con_p + st_size; // prefixes are stored after the Container
        prefixes_len = cache_item_len - st_size;
    }

    // case: not cached - decode and cache it
    if (!con_p) {
        // decode
        unsigned b64_len = snip_len; // maximum length of b64 - it will be shorter if snip includes prefixes too
        base64_decode (snip, &b64_len, (uint8_t*)&con);
        con.repeats = BGEN24 (con.repeats);

        ASSERT (con_nitems (con) <= MAX_FIELDS, "A container of %s has %u items which is beyond MAX_FIELDS=%u. Please upgrade to latest version of genozip to access this file.",
                ctx->tag_name, con_nitems (con), MAX_FIELDS);

        // get the did_i for each dict_id - unfortunately we can only store did_i up to 254 (changing this would be a change in the file format)
        for (uint32_t item_i=0; item_i < con_nitems (con); item_i++)
            if (con.items[item_i].dict_id.num) { // not a prefix-only item
                DidIType did_i = ctx_get_existing_did_i (vb, con.items[item_i].dict_id);
                ASSERT (did_i != DID_I_NONE, "analyzing a %s container in vb=%u: unable to find did_i for item %s/%s",
                        ctx->tag_name, vb->vblock_i, dtype_name_z(con.items[item_i].dict_id), dis_dict_id (con.items[item_i].dict_id).s);

                con.items[item_i].did_i_small = (did_i <= 254 ? (uint8_t)did_i : 255);
            }
            else 
                con.items[item_i].did_i_small = 255;
 
        // get prefixes
        unsigned st_size = con_sizeof (con);
        con_p        = &con;
        prefixes     = (b64_len < snip_len) ? &snip[b64_len+1]       : NULL;
        prefixes_len = (b64_len < snip_len) ? snip_len - (b64_len+1) : 0;

        ASSERT (prefixes_len <= CONTAINER_MAX_PREFIXES_LEN, "ctx=%s: prefixes_len=%u longer than CONTAINER_MAX_PREFIXES_LEN=%u", 
                ctx->tag_name, prefixes_len, CONTAINER_MAX_PREFIXES_LEN);

        ASSERT (!prefixes_len || prefixes[prefixes_len-1] == CON_PX_SEP, 
                "ctx=%s: prefixes array does end with a CON_PX_SEP: %.*s", ctx->tag_name, prefixes_len, prefixes);

        ASSERT (st_size + prefixes_len <= 65535, "ctx=%s: st_size=%u + prefixes_len=%u too large", 
                ctx->tag_name, st_size, prefixes_len);

        // first encounter with Container for this context - allocate the cache index
        if (!cache_exists) {
            buf_alloc (vb, &ctx->con_index,    0, ctx->word_list.len, uint32_t, 1, "contexts->con_index");
            buf_alloc_zero (vb, &ctx->con_len, 0, ctx->word_list.len, uint16_t, 1, "contexts->con_len");
        }

        // case: add container to cache index - only if it is not a singleton (i.e. has word_index). 
        // note: singleton containers only occur in old files compressed with v8 (since v9 no_stons is set in container_seg_do)
        if (word_index != WORD_INDEX_NONE) 
            *ENT (uint32_t, ctx->con_index, word_index) = (uint32_t)ctx->con_cache.len;

        // place Container followed by prefix in the cache (even if its a singleton)
        buf_alloc (vb, &ctx->con_cache, st_size + prefixes_len + CONTAINER_MAX_SELF_TRANS_CHANGE, 0, char, 2, "contexts->con_cache");
        
        char *cached_con = AFTERENT (char, ctx->con_cache);
        buf_add (&ctx->con_cache, con_p, st_size);
        if (prefixes_len) buf_add (&ctx->con_cache, prefixes, prefixes_len);

        con_p    = (ContainerP)cached_con; // update so we reconstruct from translated cached item
        prefixes = cached_con + st_size;

        // if item[0] is a translator-only item, use it to translate the Container itself (used by SAM_OPTIONAL)
        ContainerItem *item0 = &con.items[0];
        
        if (vb->translation.trans_containers && !item0->dict_id.num && item0->translator) {
            int32_t prefixes_len_change = (DT_FUNC(vb, translator)[item0->translator](vb, ctx, cached_con, st_size + prefixes_len, 0, false));  
            ASSERT (prefixes_len_change <= CONTAINER_MAX_SELF_TRANS_CHANGE, 
                    "ctx=%s: prefixes_len_change=%d exceeds range maximum %u", 
                    ctx->tag_name, prefixes_len_change, CONTAINER_MAX_SELF_TRANS_CHANGE);
            
            prefixes_len       += prefixes_len_change;
            ctx->con_cache.len += prefixes_len_change;
        }

        // finally, record the length (possibly updated by the translator) - but not for singletons
        if (word_index != WORD_INDEX_NONE) 
            *ENT (uint16_t, ctx->con_len, word_index) = (uint16_t)(st_size + prefixes_len);
    }

    if (out_prefixes) {
        *out_prefixes = prefixes;
        *out_prefixes_len = prefixes_len;
    }

    return con_p;
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
        iprintf ("item %u: dict_id=%s separator={ %u %u } translator=%u\n",
                 i, dis_dict_id (con->items[i].dict_id).s,  
                 con->items[i].separator[0], con->items[i].separator[1], con->items[i].translator);
    
    fflush (info_stream);
}

// Translators reconstructing last_value as a little endian binary
#define SET_n(type,mn,mx) type n = (type)ctx->last_value.i; \
                           ASSPIZ (ctx->last_value.i>=(int64_t)(mn) && ctx->last_value.i<=(int64_t)(mx), "Error: Failed to convert %s to %s because of bad data in line %"PRIu64" of the %s file: value of %s=%"PRId64" is out of range [%"PRId64"-%"PRId64"]",\
                                   dt_name (z_file->data_type), dt_name (flag.out_dt), vb->line_i, dt_name (z_file->data_type), \
                                   ctx->tag_name, ctx->last_value.i, (int64_t)(mn), (int64_t)(mx))

TRANSLATOR_FUNC (container_translate_I8)       { SET_n (int8_t,   INT8_MIN,  INT8_MAX  );              RECONSTRUCT1(n);       return 0; }
TRANSLATOR_FUNC (container_translate_U8)       { SET_n (uint8_t,  0,         UINT8_MAX );              RECONSTRUCT1((char)n); return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_I16) { SET_n (int16_t,  INT16_MIN, INT16_MAX ); n=LTEN16(n); RECONSTRUCT(&n,2);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_U16) { SET_n (uint16_t, 0,         UINT16_MAX); n=LTEN16(n); RECONSTRUCT(&n,2);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_I32) { SET_n (int32_t,  INT32_MIN, INT32_MAX ); n=LTEN32(n); RECONSTRUCT(&n,4);     return 0; }
TRANSLATOR_FUNC (container_translate_LTEN_U32) { SET_n (uint32_t, 0,         UINT32_MAX); n=LTEN32(n); RECONSTRUCT(&n,4);     return 0; }
