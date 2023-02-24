// ------------------------------------------------------------------
//   qname.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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

static inline Did did_by_q (QType q) 
{
    ASSERT (q >= 0 && q < NUM_QTYPES, "Invalid q=%d", q);

    return (Did[]){ FASTQ_QNAME, FASTQ_QNAME2, FASTQ_LINE3 }[q]; // note: SAM, KRAKEN and FASTQ have the same dids for QNAMEs
}

static inline unsigned qname_get_px_str_len (QnameFlavorStruct *qfs)
{
    unsigned i=0; 
    while (qfs->px_strs[i]) i++;
    return i;
}

// remove an unused mate item (i.e. one that still have a CI0_SKIP) from the container
static void qname_remove_skipped_mate (SmallContainerP con_no_skip, uint32_t *prefix_no_skip_len, 
                                       ConstSmallContainerP con, STRp (con_prefix)) // out
{
    // remove last item, if it is an unused QmNAME (i.e. with a CI0_SKIP) 
    *con_no_skip = *con;
    if (con->items[con->nitems_lo-1].separator[0] == CI0_SKIP) {
        con_no_skip->nitems_lo--;

        // copy prefix, except last item if it is a SKIP 
        rom px_skip = con_prefix ? memchr (con_prefix, CI0_SKIP, con_prefix_len) : 0;
        *prefix_no_skip_len = px_skip ? (px_skip - con_prefix) : *prefix_no_skip_len;
    }
}

static void qname_genarate_qfs_with_mate (QnameFlavorStruct *qfs)
{
    // mate item is always the final item
    int mate_item_i = qfs->con.nitems_lo-1;

    ASSERT (qfs->con.items[mate_item_i].dict_id.num == _SAM_QmNAME, "Expecting the last item of qfs=%s to be a mate item", qfs->name);
    ASSERT (qfs->con.items[mate_item_i].separator[0] == CI0_SKIP, "Expecting QmNAME of %s to have a CI0_SKIP separator", qfs->name);

    // replace CI0_SKIP with fixed-len of length 1
    qfs->con.items[mate_item_i].separator[0] = CI0_FIXED_0_PAD;
    qfs->con.items[mate_item_i].separator[1] = 1;
    
    // name
    strcpy (&qfs->name[strlen(qfs->name)], "/");
    
    // examples
    for (int i=0; i < QFS_MAX_EXAMPLES; i++) {
        unsigned example_len = strlen(qfs->example[i]);
        if (!example_len) break;

        strcpy (&qfs->example[i][example_len], "/1");
    }

    // case 1: previous item is not fixed - move its separator (possibly 0) to the mate item and make it '/'
    if (qfs->con.items[mate_item_i-1].separator[0] != CI0_FIXED_0_PAD) {
        qfs->con.items[mate_item_i].separator[0] = qfs->con.items[mate_item_i-1].separator[0]; // 0 or ' '
        qfs->con.items[mate_item_i].separator[1] = 0;
        qfs->con.items[mate_item_i-1].separator[0] = '/';
    }

    // case 2: previous item is fixed - add / as a prefix. 
    else { // eg con_roche_454
        // qfs must already a PX_MATE prefix item for the mate item, if previous item is fixed
        ASSERT (qname_get_px_str_len (qfs) > mate_item_i, "QnameFlavor=%s: expecting PX_MATE prefix item to exist for mate_item_i=%u",
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

// we need to prepare the containers only once and they can serve data types, as the containers
// are data-type independent (since all data types use the dict_id for QNAMEs)
void qname_zip_initialize (void)
{
    DO_ONCE {
        for (QnameFlavorStruct *qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {

            // case: this qfs version of the next qfs, just with room for /1 /2 - generate it now
            if (!qfs->name[0]) { 
                *qfs = *(qfs+1);
                qfs->con = *qfs->con_template;
                qname_genarate_qfs_with_mate (qfs);
            }
            else
                qfs->con = *qfs->con_template; // create container

            // verify that the fields are consecutive Q*NAME contexts, and verify QmNAME requirements
            ContextP expected_zctx = ZCTX(SAM_Q0NAME);
            for (int item_i=0; item_i < qfs->con.nitems_lo; item_i++) {
                uint64_t dnum = qfs->con.items[item_i].dict_id.num;

                // verify QmNAME requirements
                if (dnum == _SAM_QmNAME) 
                    ASSERT (item_i == qfs->con.nitems_lo - 1, "Item=%d QmNAME of %s is not the last item in the container", item_i, qfs->name);

                // veriry Q*NAME requirements
                else { 
                    ASSERT (dnum, "dnum=0 for item_i=%u in qfs=\"%s\"", item_i, qfs->name);

                    ASSERT (dnum == expected_zctx->dict_id.num, "Expecting item #%u of %s to be %s", item_i, qfs->name, expected_zctx->tag_name);

                    ASSERT (expected_zctx->did_i < SAM_Q0NAME + MAX_QNAME_ITEMS,
                            "Unexpected dict_id=%s in container of flavor %s. Perhaps MAX_QNAME_ITEMS=%u needs updating?", 
                            dis_dict_id(qfs->con.items[item_i].dict_id).s, qfs->name, MAX_QNAME_ITEMS);

                    expected_zctx++;
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

            // remove QmNAME from the container if not mated 
            SmallContainer con_no_skip = qfs->con; // copy
            uint32_t prefix_no_skip_len=qfs->con_prefix_len;
            qname_remove_skipped_mate (&con_no_skip, &prefix_no_skip_len, &qfs->con, STRa (qfs->con_prefix));
            
            // prepare container snips - with correct dict_ids for each QType
            for (QType q=QNAME1; q < NUM_QTYPES; q++) {
                Did first_did_i = did_by_q (q);
                
                for (unsigned item_i=0; item_i < con_no_skip.nitems_lo; item_i++) 
                    if (con_no_skip.items[item_i].dict_id.num != _SAM_QmNAME) 
                        con_no_skip.items[item_i].dict_id = dt_fields[DT_FASTQ].predefined[first_did_i + item_i + 1].dict_id; 
                    else
                        con_no_skip.items[item_i].dict_id = dt_fields[DT_FASTQ].predefined[first_did_i + MAX_QNAME_ITEMS].dict_id; 

                // prepare container snip for this QF
                qfs->con_snip_lens[q] = sizeof (qfs->con_snips[0]);            
                container_prepare_snip ((Container*)&con_no_skip, qfs->con_prefix, prefix_no_skip_len, qSTRi(qfs->con_snip,q));
            }

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

                ASSERT (!qfs->is_int[qfs->numeric_items[i]], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an integer_item and numeric_item",
                        qfs->name, qfs->numeric_items[i]);
            }

            for (unsigned i=0; qfs->hex_items[i] != -1; i++) {
                qfs->is_hex[qfs->hex_items[i]] = true;

                ASSERT (!qfs->is_int[qfs->hex_items[i]], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an hex_item and integer_item",
                        qfs->name, qfs->hex_items[i]);
            }

            for (unsigned i=0; qfs->in_local[i] != -1; i++)
                qfs->is_in_local[qfs->in_local[i]] = true;

            ASSERT (qfs->ordered_item1 == -1 || qfs->is_int[qfs->ordered_item1] || qfs->is_numeric[qfs->ordered_item1] || qfs->is_hex[qfs->ordered_item1], 
                    "Error in definition of QNAME flavor=%s: item=%u is one of \"ordered_item\" - expecting it to be is_integer or is_numeric or is_hex", qfs->name, qfs->ordered_item1);
        }

        seg_prepare_snip_other (SNIP_COPY, (DictId)_SAM_QNAME, false, 0, copy_qname); // QNAME dict_id is the same for SAM, FASTQ, KRAKEN
    }
}

void qname_seg_initialize (VBlockP vb, QType q, int pairing)
{
    QnameFlavor qfs = segconf.qname_flavor[q];
    if (!qfs) return;

    Did did_q = did_by_q (q);

    // initialize for pairing if needed
    if (pairing) 
        for (int i=0; i < MAX_QNAME_ITEMS + 1; i++) 
            fastq_initialize_ctx_for_pairing (vb, did_q + i);

    // consolidate all items under did_q in --stats
    ctx_consolidate_statsN (vb, did_q, did_q + 1, MAX_QNAME_ITEMS); 

    // set STORE_INT as appropriate
    int ctx_delta = did_q - SAM_QNAME;
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
}

// note: we run this function only in discovery, not in segging, because it is quite expensive - checking all numerics.
// returns 0 if qname is indeed the flavor, or error code if not

QnameTestResult qname_test_flavor (STRp(qname), QType q, QnameFlavor qfs, bool quiet) // out
{
    // return embedded qname2, eg "82a6ce63-eb2d-4812-ad19-136092a95f3d" in "@ERR3278978.1 82a6ce63-eb2d-4812-ad19-136092a95f3d/1"
    if (q != QANY && qfs->only_q != QANY && qfs->only_q != q
        && !(qfs->only_q == Q1or3 && (q == QNAME1 || QLINE3))) return QTR_WRONG_Q; 

    if (q==QNAME2 && segconf.tech != qfs->qname1_tech) 
        return QTR_TECH_MISMATCH;

    if (!qname_len) return QTR_QNAME_LEN_0;

    // test fixed length, if applicable
    if (qfs->fixed_len && qname_len != qfs->fixed_len) return QTR_FIXED_LEN_MISMATCH;

    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item, quiet ? NULL : qfs->name);
    if (!n_items) return QTR_CONTAINER_MISMATCH; // failed to split according to container - qname is not of this flavor

    // items cannot include whitespace or non-printable chars
    for (uint32_t item_i=0; item_i < n_items; item_i++)
        if (!str_is_no_ws (STRi(item, item_i))) return QTR_BAD_CHARS;

    // check that all the items expected to be numeric (leading zeros ok) are indeed so
    for (const int *item_i = qfs->numeric_items; *item_i != -1; item_i++)
        if (!qfs->is_hex[*item_i] && !str_is_numeric (STRi(item, *item_i))) return QTR_BAD_NUMERIC;

    // check that all the items expected to be lower case hexadecimal are indeed so
    for (const int *item_i = qfs->hex_items; *item_i != -1; item_i++)
        if (!str_is_hexlo (STRi(item, *item_i))) return QTR_BAD_HEX;

    // check that all the items expected to be integer (no leading zeros) are indeed so
    for (const int *item_i = qfs->integer_items; *item_i != -1; item_i++)
        if (!qfs->is_hex[*item_i] && !str_is_int (STRi(item, *item_i))) return QTR_BAD_INTEGER;

    return QTR_SUCCESS; // yes, qname is of this flavor
}

// called for the first line in segconf.running
void qname_segconf_discover_flavor (VBlockP vb, QType q, STRp(qname))
{
    ASSERTNOTZERO (segconf.running, "");

    static rom reasons[] = QTR_NAME;

    segconf.qname_flavor[q] = NULL; // unknown

    for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {
        QnameTestResult reason = QTR_SUCCESS;
        if (!(reason = qname_test_flavor (STRa(qname), q, qfs, true))) {
            segconf.qname_flavor[q] = qfs;
            segconf.tech = qfs->tech; // note: if this is QNAME2, we update tech according to QNAME2 (instead of NCBI)
            
            if (qfs->seq_len_item != -1) 
                segconf.seq_len_dict_id = qfs->con.items[qfs->seq_len_item].dict_id;
            
            qname_seg_initialize (vb, q, 0); // so the rest of segconf.running can seg fast using the discovered container

            if (flag.debug_qname) 
                iprintf ("%.*s is DISCOVERED as %s - for %s\n", STRf(qname), qfs->name, qtype_name(q));
            
            break;
        }

        else if (flag.debug_qname && reason) 
            iprintf ("%.*s is not flavor \"%s\". Reason: %s\n", STRf(qname), qfs->name, reasons[reason]);
    }

    if (flag.debug_qname && !segconf.qname_flavor[q])
        iprintf ("%.*s - flavor is NOT DISCOVERED - for %s\n", STRf(qname), qtype_name(q));

    // when optimizing qname with --optimize_DESC - capture the correct TECH ^ but set flavor to Genozip-opt
    if (flag.optimize_DESC) 
        segconf.qname_flavor[q] = (q==QNAME1) ? &qf[NUM_QFs-1] : NULL; // relying on Genozip-opt being last

    DO_ONCE // unit test
    if (flag.debug || flag.debug_qname)
        for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {

            for (int i=0; i < QFS_MAX_EXAMPLES; i++) {
                if (!qfs->example[i][0]) break;

                if (flag.debug_qname) iprintf ("Unit-testing qname flavor=\"%s\" example=\"%s\"\n", qfs->name, qfs->example[i]);
                QnameTestResult reason = qname_test_flavor (qfs->example[i], strlen(qfs->example[i]), QANY, qfs, !flag.debug_qname); 
                ASSERT (reason == QTR_SUCCESS, "Failed to identify qname \"%s\" as %s: reason=%s", qfs->example[i], qfs->name, reasons[reason]);
            }
        }
}

// attempt to seg according to the qf - return true if successful
bool qname_seg_qf (VBlockP vb, QType q, STRp(qname), unsigned add_additional_bytes)
{
    QnameFlavor qfs = segconf.qname_flavor[q];
    ASSERTNOTNULL (qfs);
    
    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item, NULL);
    if (!n_items) return false; // failed to split according to container - qname is not of this flavor

    // seg
    int64_t value, prev_value;
    ContextP qname_ctx = CTX(did_by_q(q));
    
    // seg container
    seg_by_ctx (vb, STRi(qfs->con_snip, q),  qname_ctx, qfs->num_seps + add_additional_bytes); // account for container separators, prefixes and caller-requested add_additional_bytes 

    for (unsigned item_i=0; item_i < qfs->con.nitems_lo; item_i++) {

        // get item ctx - Q?NAME did_i are immediately following QNAME, and QmNAME si the last.
        const ContainerItem *item = &qfs->con.items[item_i];
        value=-1;

        // calculate context - a bit messy, but faster than looking up by dict_id (relying on Q?NAME being sequential in the container)
        ContextP item_ctx = qname_ctx + ((item->dict_id.num == _SAM_QmNAME) ? MAX_QNAME_ITEMS : (1+item_i)); // note: QmName, if exists, is always the last item
                        
        // case: this is the file is collated by qname - delta against previous
        if ((item_i == qfs->ordered_item1 || item_i == qfs->ordered_item2) && 
            (segconf.is_collated || VB_DT(FASTQ) || VB_DT(KRAKEN) || qfs->id == QF_GENOZIP_OPT) &&
            (ctx_has_value_in_prev_line_(vb, item_ctx) || vb->line_i==0) &&
            !(item_lens[item_i] >= 2 && items[item_i][0] == '0') && // can't yet handle reproducing leading zeros with a delta
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
                seg_add_to_local_text (vb, item_ctx, STRi(item, item_i), LOOKUP_SIMPLE, item_lens[item_i]);
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

void qname_seg (VBlockP vb, QType q, STRp (qname), unsigned add_additional_bytes)  // account for characters in addition to the field
{
    START_TIMER;

    QnameFlavor flavor = segconf.qname_flavor[q];
    ContextP qname_ctx = CTX(did_by_q(q)); 

    // copy if identical to previous (> 50% of lines in collated SAM/BAM) - small improvement in compression and compression time
    // no need in is_sorted, as already handled in sam_seg_QNAME with buddy 
    if ((VB_DT(SAM) || VB_DT(BAM)) && !segconf.is_sorted && vb->line_i && is_same_last_txt (vb, qname_ctx, STRa(qname))) {
        seg_by_ctx (vb, STRa(copy_qname), qname_ctx, qname_len + add_additional_bytes);
        goto done;
    }

    bool success = false;
    if (flavor)
        success = qname_seg_qf (vb, q, STRa(qname), add_additional_bytes);

    if (!success) 
        tokenizer_seg (VB, qname_ctx, STRa(qname), 
                       VB_DT(FASTQ) ? sep_with_space : sep_without_space, 
                       add_additional_bytes);

    if (segconf.running) {
        // collect the first 6 qnames / q2names, if flavor is unknown    
        if (vb->line_i < 6 && !flavor) // unrecognized flavor
            memcpy (segconf.unknown_flavor_qnames[q][vb->line_i], qname, MIN_(qname_len, UNK_QNANE_LEN));
    }

done:
    seg_set_last_txt (vb, qname_ctx, STRa(qname)); // needed also for kraken in sam
    COPY_TIMER (qname_seg);
}

rom segconf_qf_name (QType q)
{
    return segconf.qname_flavor[q] ? segconf.qname_flavor[q]->name : "N/A";
}

QnameFlavorId segconf_qf_id (QType q)
{
    return segconf.qname_flavor[q]->id;
}

rom qtype_name (QType q)
{
    static rom qtype_names[] = QTYPE_NAME;
    
    return (q >= 0 && q < ARRAY_LEN (qtype_names)) ? qtype_names[q] : "Invalid qtype";
}
