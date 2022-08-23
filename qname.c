// ------------------------------------------------------------------
//   qname.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "qname.h"
#include "tokenizer.h"
#include "buffer.h"
#include "vblock.h"
#include "segconf.h"
#include "strings.h"
#include "seg.h"
#include "sam.h"
#include "fastq.h"
#include "stats.h"
#include "qname_flavors.h"
#include "file.h"

static STRl(copy_qname, 50);

static inline unsigned qname_get_px_str_len (QnameFlavorStruct *qfs)
{
    unsigned i=0; 
    while (qfs->px_strs[i]) i++;
    return i;
}

static void qname_remove_skip (SmallContainerP con_no_skip, uint32_t *prefix_no_skip_len, 
                               ConstSmallContainerP con, STRp (con_prefix)) // out
{
    // all container data, except for items
    memcpy (con_no_skip, con, sizeof (MiniContainer) - sizeof (ContainerItem));
    con_no_skip->nitems_lo = 0;

    // copy items, except SKIP items    
    for (int item_i=0; item_i < con->nitems_lo; item_i++)
        if (con->items[item_i].separator[0] != CI0_SKIP) 
            con_no_skip->items[con_no_skip->nitems_lo++] = con->items[item_i];

    // copy prefix, except last item if it is a SKIP (we only support SKIP for final prefix items)
    rom px_skip = con_prefix ? memchr (con_prefix, CI0_SKIP, con_prefix_len) : 0;
    *prefix_no_skip_len = px_skip ? (px_skip - con_prefix) : *prefix_no_skip_len;
}

static void qname_genarate_qfs_with_mate (QnameFlavorStruct *qfs)
{
    // find mate item (can't be item 0 ; usually, but not always, its the last item)
    int mate_item_i;
    for (mate_item_i=qfs->con.nitems_lo-1; mate_item_i >= 1; mate_item_i--)
        if (qfs->con.items[mate_item_i].dict_id.num == _SAM_QmNAME) break;

    ASSERT (mate_item_i, "Can't find mate item for qfs=%s", qfs->name);

    // replace CI0_SKIP with fixed-len of length 1
    qfs->con.items[mate_item_i].separator[0] = CI0_FIXED_0_PAD;
    qfs->con.items[mate_item_i].separator[1] = 1;
    
    // name
    strcpy (&qfs->name[strlen(qfs->name)], "/");
    
    // examples
    for (int i=0; i < QFS_MAX_EXAMPLES; i++) {
        unsigned example_len = strlen(qfs->example[i]);
        if (!example_len) break;

        char *after = &qfs->example[i][example_len];
        if (!qfs->fq_only)
            strcpy (after, "/1");
        else {
            str_split (qfs->example[i], example_len, 0, ' ', str, false);
            char *slash = (char *)strs[qfs->fq_only-1] + str_lens[qfs->fq_only-1];
            memmove (slash+2, slash, after - slash);
            slash[0] = '/';
            slash[1] = '1';
        }
    }

    // case 1: previous item is not fixed - move its separator (possibly 0) to the mate item and make it '/'
    if (qfs->con.items[mate_item_i-1].separator[0] != CI0_FIXED_0_PAD) {
        qfs->con.items[mate_item_i].separator[0] = qfs->con.items[mate_item_i-1].separator[0]; // 0 or ' '
        qfs->con.items[mate_item_i].separator[1] = 0;
        qfs->con.items[mate_item_i-1].separator[0] = '/';
    }

    // case 2: previous item is fixed - add / as a prefix. 
    else { // eg con_roche_454
        // qfs must already have prefix item for the mate item - initially empty string
        ASSERT (qname_get_px_str_len (qfs) > mate_item_i, "QnameFlavor=%s: expecting prefix to exist for mate_item_i=%u",
                qfs->name, mate_item_i);

        qfs->px_strs[mate_item_i] = "/";
    }

    // add new item as numeric item (at beginning of array - easier, and the order doesn't matter)
    if (qfs->con.items[mate_item_i].separator[0] == CI0_FIXED_0_PAD) {
        memmove (&qfs->numeric_items[1], qfs->numeric_items, sizeof(qfs->numeric_items) - sizeof(qfs->numeric_items[0])); // make room
        qfs->numeric_items[0] = mate_item_i;
    }

    qfs->num_seps++; // we added / to length of prefixes/seperator

    if (qfs->fixed_len) // note: qfs can have fixed_len even if items aren't fixed (eg IonTorrent)
        qfs->fixed_len += 2;
}

// note: we need to re-initialize in each file, lest the data type has changed
void qname_zip_initialize (Did qname_did_i)
{
    ContextP qname_zctx = ZCTX(qname_did_i);
    ContextP line3_zctx = ZCTX(FASTQ_LINE3);

    // we need to prepare the containers only once and they can serve all files, as the containers
    // are data-type independent (since all data types use the dict_id Q?NAME for qnames)
    DO_ONCE {
        for (QnameFlavorStruct *qfs=&qf[0]; qfs < &qf[NUM_QFs + NUM_QF_L3s]; qfs++) {

            // case: this qfs version of the next qfs, just with room for /1 /2 - generate it now
            if (!qfs->name[0]) { 
                *qfs = *(qfs+1);
                qfs->con = *qfs->con_template;
                qname_genarate_qfs_with_mate (qfs);
            }
            else
                qfs->con = *qfs->con_template; // create container

            qfs->qname2 = -1; // initialize pessimistically
            ContextP next_zctx = (qfs < &qf[NUM_QFs] ? qname_zctx : line3_zctx) + 1;
            for (int item_i=0; item_i < qfs->con.nitems_lo; item_i++) {
                uint64_t dnum = qfs->con.items[item_i].dict_id.num;
                // set qname2
                if (dnum == _FASTQ_QNAME2) {
                    ASSERT0 (qfs->fq_only && qfs->tech==TECH_UNKNOWN, "Bad embedded qname2 definition"); // a con with FASTQ_QNAME2 must have these properties
                    qfs->qname2 = item_i;
                }
                // verify dict_id (in some cases)
                else if (dnum != _SAM_QmNAME && dnum != _FASTQ_COPY_Q && (Z_DT(DT_FASTQ) || qfs < &qf[NUM_QFs])) {
                    ASSERT (next_zctx->dict_id.num == dnum, "Expecting item #%u of %s to have to be %s", item_i, qfs->name, next_zctx->tag_name);
                    next_zctx++;
                }
            }

            // case: this QF needs a container prefix - generate it
            if (qfs->px_strs[0]) {
                qfs->con_prefix[qfs->con_prefix_len++] = CON_PX_SEP; 
                qfs->con_prefix[qfs->con_prefix_len++] = CON_PX_SEP; // no container-wide prefix

                for (rom *px = qfs->px_strs; *px ; px++) {
                    unsigned len = strlen (*px);
                    memcpy (&qfs->con_prefix[qfs->con_prefix_len], *px, len);
                    qfs->con_prefix_len += len;
                    qfs->con_prefix[qfs->con_prefix_len++] = CON_PX_SEP; 
                }
            } 

            // remove SKIP item, if one exists, from container and prefix
            SmallContainer con_no_skip = qfs->con; // copy
            uint32_t prefix_no_skip_len=qfs->con_prefix_len;
            qname_remove_skip (&con_no_skip, &prefix_no_skip_len, &qfs->con, STRa (qfs->con_prefix));
            
            // prepare container snip for this QF
            qfs->con_snip_len = sizeof (qfs->con_snip);            
            container_prepare_snip ((Container*)&con_no_skip, qfs->con_prefix, prefix_no_skip_len, qfs->con_snip, &qfs->con_snip_len);

            // prepare snip of identical container, except dict_ids are of QNAME2
            Did next_qname2_did_i = FASTQ_Q0NAME2;
            for (unsigned item_i=0; item_i < con_no_skip.nitems_lo; item_i++) 
                if (con_no_skip.items[item_i].dict_id.num != _SAM_QmNAME && con_no_skip.items[item_i].dict_id.num != _FASTQ_QNAME2) {
                    con_no_skip.items[item_i].dict_id = dt_fields[DT_FASTQ].predefined[next_qname2_did_i].dict_id; // only FASTQ has QNAME2 dict_ids, this is not necessarily the current data type
                    next_qname2_did_i++;
                }

            qfs->con_snip2_len = sizeof (qfs->con_snip2);            
            container_prepare_snip ((Container*)&con_no_skip, qfs->con_prefix, prefix_no_skip_len, qfs->con_snip2, &qfs->con_snip2_len);

            // in qname_flavors.h, we keep lists in the form of index lists, for maintenanbility. now
            // we convert them to a bitmap for ease of segging.
            for (unsigned i=0; qfs->integer_items[i] != -1; i++) {
                ASSERT (qfs->con.items[qfs->integer_items[i]].separator[0] != CI0_FIXED_0_PAD,
                        "since qfs=%s item=%u has seperator[0] == CI0_FIXED_0_PAD, expecting it to be numeric, not int", qfs->name, qfs->integer_items[i]);
                qfs->is_int[qfs->integer_items[i]] = true;
            }

            for (unsigned i=0; qfs->numeric_items[i] != -1; i++) {
                ASSERT (qfs->con.items[qfs->numeric_items[i]].separator[0] == CI0_FIXED_0_PAD,
                        "since qfs=%s item=%u has seperator[0] != CI0_FIXED_0_PAD, expecting it to be int, not numeric", qfs->name, qfs->numeric_items[i]);
                qfs->is_numeric[qfs->numeric_items[i]] = true;
            }

            for (unsigned i=0; qfs->hex_items[i] != -1; i++)
                qfs->is_hex[qfs->hex_items[i]] = true;
            
            for (unsigned i=0; qfs->in_local[i] != -1; i++)
                qfs->is_in_local[qfs->in_local[i]] = true;
        }
    }

    // prepare copy_qname snip - every time, as container_did_i changes between data types
    seg_prepare_snip_other (SNIP_COPY, qname_zctx->dict_id, false, 0, copy_qname); // QNAME dict_id is the same for SAM, FASTQ, KRAKEN
}

static void qname_seg_initialize_do (VBlockP vb, QnameFlavor qfs, QnameFlavor qfs2, Did qname_did_i, Did st_did_i, unsigned num_qname_items, int pairing)
{
    if (!qfs) return;

    // consolidate all items under QNAME in --stats
    ContextP qname_ctxs[num_qname_items+1];
    for (int i=0; i < num_qname_items+1; i++) {
        qname_ctxs[i] = CTX(qname_did_i + i);

        if (pairing) {
            qname_ctxs[i]->no_stons = true; // prevent singletons, so pair_1 and pair_2 are comparable based on b250 only
            
            if (pairing == PAIR_READ_2)
                ctx_create_node (vb, qname_did_i + i, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_mate_lookup }, 2); // required by ctx_convert_generated_b250_to_mate_lookup
        }
    }

    ctx_consolidate_stats_(vb, st_did_i, num_qname_items+1, qname_ctxs); 

    // set STORE_INT as appropriate
    int ctx_delta = (qname_did_i == FASTQ_QNAME2 ? (FASTQ_QNAME2 - FASTQ_DESC) : 0); // con.items[] are expressed relative to QNAME, for QNAME2, we need to shift this much
    if (qfs->ordered_item1   != -1) (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->ordered_item1].dict_id))    ->flags.store = STORE_INT; 
    if (qfs->ordered_item2   != -1) (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->ordered_item2].dict_id))    ->flags.store = STORE_INT; 
    if (qfs->range_end_item1 != -1) (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->range_end_item1-1].dict_id))->flags.store = STORE_INT; 
    if (qfs->range_end_item2 != -1) (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->range_end_item2-1].dict_id))->flags.store = STORE_INT; 
    if (qfs->seq_len_item    != -1) (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->seq_len_item].dict_id))     ->flags.store = STORE_INT; 

    // no singletons for in_local items (but not if int or numeric, as these go to seg_integer_or_not which interpret no_stons differently)
    for (int i=0; qfs->in_local[i] != -1; i++)
        if (!qfs->is_int[qfs->in_local[i]] && !qfs->is_numeric[qfs->in_local[i]])
            (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->in_local[i]].dict_id))->no_stons = true;

    // hex items
    for (int i=0; qfs->hex_items[i] != -1; i++)
        if (qfs->is_in_local[qfs->hex_items[i]] && (qfs->is_int[qfs->hex_items[i]] || qfs->is_numeric[qfs->hex_items[i]]))
            (ctx_delta + ctx_get_ctx (vb, qfs->con.items[qfs->hex_items[i]].dict_id))->ltype = LT_DYN_INT_h;

    if (qfs->qname2 != -1 && qfs2)
        qname_seg_initialize_do (vb, qfs2, NULL, FASTQ_QNAME2, st_did_i, num_qname_items-1/*-1 bc no mate*/, pairing);
}

void qname_seg_initialize (VBlockP vb, Did qname_did_i) 
{   
    qname_seg_initialize_do (vb, segconf.qname_flavor, segconf.qname_flavor2, qname_did_i, qname_did_i, MAX_QNAME_ITEMS, flag.pair); 

    if (segconf.line3_flavor) // true for FASTQ with SRA line3, but not SRA qname
        qname_seg_initialize_do (vb, segconf.line3_flavor, 0, FASTQ_LINE3, FASTQ_LINE3, MAX_LINE3_ITEMS, 0); 
}

// note: we run this function only in discovery, not in segging, because it is quite expesive - checking all numerics.
// returns 0 if qname is indeed the flavor, or error code if not
int qname_test_flavor (STRp(qname), QnameFlavor qfs,
                       pSTRp (qname2)) // out
{
    if (!qname_len) return -1;

    // test fixed length, if applicable
    if (qfs->fixed_len && qname_len != qfs->fixed_len) return -2;

    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item);
    if (!n_items) return -4; // failed to split according to container - qname is not of this flavor

    // items cannot include whitespace or non-printable chars
    for (uint32_t item_i=0; item_i < n_items; item_i++)
        if (!str_is_no_ws (STRi(item, item_i))) return -6;

    // check that all the items expected to be numeric (leading zeros ok) are indeed so
    for (const int *item_i = qfs->numeric_items; *item_i != -1; item_i++)
        if (!qfs->is_hex[*item_i] && !str_is_numeric (STRi(item, *item_i))) return -7;

    // check that all the items expected to be lower case hexadecimal are indeed so
    for (const int *item_i = qfs->hex_items; *item_i != -1; item_i++)
        if (!str_is_hexlo (STRi(item, *item_i))) return -8;

    // check that all the items expected to be integer (no leading zeros) are indeed so
    for (const int *item_i = qfs->integer_items; *item_i != -1; item_i++)
        if (!qfs->is_hex[*item_i] && !str_is_int (STRi(item, *item_i))) return -9;

    // return embedded qname2, eg "82a6ce63-eb2d-4812-ad19-136092a95f3d" in "@ERR3278978.1 82a6ce63-eb2d-4812-ad19-136092a95f3d/1"
    if (qname2 && qfs->qname2 != -1) {
        *qname2 = items[qfs->qname2];
        *qname2_len = item_lens[qfs->qname2];
    }

    return 0; // yes, qname is of this flavor
}

// called for the first line in segconf.running
void qname_segconf_discover_flavor (VBlockP vb, Did qname_did_i, STRp(qname))
{
    segconf.qname_flavor = 0; // unknown
    STR0(qname2);

    for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) 
        if ((VB_DT(DT_FASTQ) || !qfs->fq_only) && !qname_test_flavor (STRa(qname), qfs, pSTRa(qname2))) {
            segconf.qname_flavor = qfs;
            segconf.tech = qfs->tech;
            
            if (qfs->seq_len_item != -1) 
                segconf.qname_seq_len_dict_id = qfs->con.items[qfs->seq_len_item].dict_id;
            
            // if it has an embedded qname2 - find the true tech from qname2 (only possible for FASTQ, in SAM we can't find the tech from NCBI qname)
            if (qname2) 
                for (QnameFlavor qf2=&qf[0]; qf2 < &qf[NUM_QFs]; qf2++) 
                    if (!qf2->fq_only && !qname_test_flavor (STRa (qname2), qf2, 0, 0)) { // qname2 cannot be "fq_only"
                        segconf.qname_flavor2 = qf2;
                        segconf.tech = qf2->tech;
                        break;
                    }

            qname_seg_initialize (vb, qname_did_i); // so the rest of segconf.running can seg fast using the discovered container
            break;
        }

    // when optimizing qname with --optimize_DESC - capture the correct TECH ^ but set flavor to Genozip-opt
    if (flag.optimize_DESC) {
        segconf.qname_flavor = &qf[NUM_QFs-1]; // Genozip-opt is last
        segconf.qname_flavor2 = NULL;
    }

    if (flag.debug || flag.debug_qname) // unit test
        for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {

            for (int i=0; i < QFS_MAX_EXAMPLES; i++) {
                unsigned example_len = strlen(qfs->example[i]);
                if (!example_len) break;

                if (flag.debug_qname) iprintf ("Testing qname flavor=\"%s\" example=\"%s\"\n", qfs->name, qfs->example[i]);
                int res = qname_test_flavor (qfs->example[i], example_len, qfs, 0, 0); 
                ASSERT (!res, "Failed to identify qname \"%s\" as %s: res=%d", qfs->example[i], qfs->name, res);
            }
        }
}

// called for the first line in segconf.running
bool qname_segconf_discover_fastq_line3_sra_flavor (VBlockP vb, STRp(line3))
{
    segconf.line3_flavor = 0; // unknown
    STR0(line3_2);

    for (QnameFlavor qfs=&qf[FIRST_QF_L3]; qfs < &qf[NUM_QFs + NUM_QF_L3s]; qfs++) 
        if (!qname_test_flavor (STRa(line3), qfs, pSTRa(line3_2))) {

            segconf.line3_flavor = qfs;
            qname_seg_initialize_do (vb, qfs, 0, FASTQ_LINE3, FASTQ_LINE3, MAX_LINE3_ITEMS, 0); 
            return true;
        }

    return false;
}

// attempt to seg according to the qf - return true if successful
bool qname_seg_qf (VBlockP vb, ContextP qname_ctx, QnameFlavor qfs, STRp(qname), bool use_qname2,
                   unsigned add_additional_bytes)
{
    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item);
    if (!n_items) return false; // failed to split according to container - qname is not of this flavor

    // seg
    int64_t value, prev_value;

    // seg container
    if (!use_qname2)
        seg_by_ctx (vb, STRa(qfs->con_snip),  qname_ctx, qfs->num_seps + add_additional_bytes); // account for container separators, prefixes and caller-requested add_additional_bytes 
    else
        seg_by_ctx (vb, STRa(qfs->con_snip2), qname_ctx, qfs->num_seps + add_additional_bytes); 

    int encountered_mName=0, encountered_QNAME2=0;
    for (unsigned item_i=0; item_i < qfs->con.nitems_lo; item_i++) {

        // get item ctx - Q?NAME did_i are immediately following QNAME, and QmNAME si the last.
        const ContainerItem *item = &qfs->con.items[item_i];
        value=-1;

        // calculate context - a bit messy, but faster than looking up by dict_id
        ContextP item_ctx = (item->dict_id.num == _SAM_QmNAME)   ? (qname_ctx + MAX_QNAME_ITEMS) 
                          : (item->dict_id.num == _FASTQ_QNAME2) ? CTX(FASTQ_QNAME2)
                          : (item->dict_id.num == _FASTQ_COPY_Q) ? CTX(FASTQ_COPY_Q)
                          :                                        (qname_ctx + 1+item_i - encountered_mName - encountered_QNAME2);
                          
        encountered_mName  |= item->dict_id.num == _SAM_QmNAME;
        encountered_QNAME2 |= item->dict_id.num == _FASTQ_QNAME2 || item->dict_id.num == _FASTQ_COPY_Q;

        // case: this is the file is sorted by qname - delta against previous
        if ((item_i == qfs->ordered_item1 || item_i == qfs->ordered_item2) && 
            !segconf.is_sorted &&
            ( (!qfs->is_hex[item_i] && str_get_int_dec (STRi(item, item_i), (uint64_t*)&value)) ||
              ( qfs->is_hex[item_i] && str_get_int_hex (STRi(item, item_i), true, false, (uint64_t*)&value)))) // lower-case hex

            seg_self_delta (vb, item_ctx, value, (qfs->is_hex[item_i] ? 'x' : 0), item_lens[item_i]);
        
        // case: end-of-range item, seg as delta vs previous item which is start-of-range
        else if ((qfs->range_end_item1 == item_i || qfs->range_end_item2 == item_i) &&
                 !flag.pair && // we don't use delta against other if we're pairing - too complicated
                 str_get_int (STRi(item, item_i), &value) && 
                 str_get_int (STRi(item, item_i-1), &prev_value)) {
            SNIP(32);
            seg_prepare_snip_other (SNIP_OTHER_DELTA, qfs->con.items[item_i-1].dict_id, true, value - prev_value, snip);
            seg_by_ctx (vb, STRa(snip), item_ctx, item_lens[item_i]);      
        }

        else if (qfs->is_in_local[item_i] && !flag.pair) { // note: we can't store in local if pairing        
            if (qfs->is_int[item_i]) 
                seg_integer_or_not (vb, item_ctx, STRi(item, item_i), item_lens[item_i]);

            else if (qfs->is_numeric[item_i])
                seg_numeric_or_not (vb, item_ctx, STRi(item, item_i), item->separator[1], item_lens[item_i]);

            else 
                seg_add_to_local_text (vb, item_ctx, STRi(item, item_i), true, item_lens[item_i]);
        }

        // TO DO - add field xor_diff allowing specifying xor_diff method (if we find a case where its useful)
        // else if (xor_diff[item_i] && !flag.pair) {
        //    seg_diff (vb, item_ctx, STRi(item, item_i), item_ctx->flags.all_the_same, item_lens[item_i]);
        //     seg_set_last_txt (vb, item_ctx, STRi(item, item_i));
        // }
        
        else if (item_i == qfs->qname2 && segconf.qname_flavor2) 
            qname_seg_qf (vb, CTX (FASTQ_QNAME2), segconf.qname_flavor2, STRi(item, item_i), true, 0);

        else if (item->dict_id.num == _FASTQ_COPY_Q) {
            // for now, this field is marked as an ordered so we never reach here. We will need a SPECIAL to
            // copy from DESC, because we need to potentially remove the /1 /2
            ABORT0 ("Copy from DESC not supported yet"); 
        }

        else if (item->separator[0] == CI0_SKIP)
            {} // no segging a skipped item

        // case: textual item
        else
            seg_by_ctx (vb, STRi(item, item_i), item_ctx, item_lens[item_i]); 

        // case: this item is qname_seq_len - set last_value by the beneficial field (CIGAR in SAM, ? in FASTQ)
        if (item_i == qfs->seq_len_item && (value >= 0 || str_get_int (STRi(item, item_i), &value))) 
            ctx_set_last_value (vb, item_ctx, value);  
    }

    return true;
}

void qname_seg (VBlockP vb, Context *qname_ctx, STRp (qname), unsigned add_additional_bytes)  // account for characters in addition to the field
{
    START_TIMER;
    
    // copy if identical to previous (> 50% of lines in collated) - small improvement in compression and compression time
    // no need in is_sorted, as already handled in sam_seg_QNAME with buddy 
    if (!segconf.is_sorted && vb->line_i && !flag.optimize_DESC && is_same_last_txt (vb, qname_ctx, STRa(qname))) {
        seg_by_ctx (vb, STRa(copy_qname), qname_ctx, qname_len + add_additional_bytes);
        goto done;
    }

    bool success = false;
    if (segconf.qname_flavor)
        success = qname_seg_qf (VB, qname_ctx, segconf.qname_flavor, STRa(qname), false, add_additional_bytes);

    if (!success) 
        tokenizer_seg (VB, qname_ctx, STRa(qname), 
                       VB_DT(DT_FASTQ) ? sep_with_space : sep_without_space, 
                       add_additional_bytes);

done:
    seg_set_last_txt (vb, qname_ctx, STRa(qname)); // needed also for kraken in sam
    COPY_TIMER (qname_seg);
}

rom qf_name (QnameFlavor qf)
{
    return qf ? qf->name : "";
}
