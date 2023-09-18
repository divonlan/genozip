// ------------------------------------------------------------------
//   qname.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "libdeflate/libdeflate.h"
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
#include "sam.h"
#include "codec.h"

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

static void qname_generate_qfs_with_mate (QnameFlavorStruct *qfs)
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
    qfs->is_mated = true;

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
                qname_generate_qfs_with_mate (qfs);
            }
            else
                qfs->con = *qfs->con_template; // create container

            // generate barcode_item2 if we have a barcode with a '+' separator
            if (qfs->barcode_item != -1 && qfs->barcode_item != qfs->con.nitems_lo - 1 &&
                qfs->con.items[qfs->barcode_item].separator[0] == '+')
                qfs->barcode_item2 = qfs->barcode_item + 1;
            else
                qfs->barcode_item2 = -1;

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
                qfs->is_integer[qfs->integer_items[i]] = true;

                ASSERT (qfs->barcode_item != qfs->integer_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an integer_item and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item2 != qfs->integer_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an integer_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);
            }

            for (unsigned i=0; qfs->numeric_items[i] != -1; i++) {
                ASSERT (qfs->con.items[qfs->numeric_items[i]].separator[0] == CI0_FIXED_0_PAD,
                        "since qfs=%s item=%u has seperator[0] != CI0_FIXED_0_PAD, expecting it to be int, not numeric", qfs->name, qfs->numeric_items[i]);
                qfs->is_numeric[qfs->numeric_items[i]] = true;

                ASSERT (!qfs->is_integer[qfs->numeric_items[i]], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an integer_item and numeric_item",
                        qfs->name, qfs->numeric_items[i]);

                ASSERT (qfs->barcode_item != qfs->numeric_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an numeric_item and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item2 != qfs->numeric_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an numeric_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);
            }

            for (unsigned i=0; qfs->hex_items[i] != -1; i++) {
                qfs->is_hex[qfs->hex_items[i]] = true;

                ASSERT (qfs->is_integer[qfs->hex_items[i]] || qfs->is_numeric[qfs->hex_items[i]], 
                        "Error in definition of QNAME flavor=%s: item=%u a hex_item must also be either an integer_item or a numeric_item",
                        qfs->name, qfs->hex_items[i]);

                ASSERT (qfs->barcode_item != qfs->hex_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an hex_item and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item2 != qfs->hex_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an hex_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);                        
            }

            for (unsigned i=0; qfs->in_local[i] != -1; i++) {
                qfs->is_in_local[qfs->in_local[i]] = true;

                ASSERT (qfs->barcode_item != qfs->in_local[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both in_local and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item2 != qfs->in_local[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both in_local and barcode_item2",
                        qfs->name, qfs->barcode_item2);                        
            }
            
            if (qfs->barcode_item != -1) {
                ASSERT (qfs->barcode_item != qfs->range_end_item1 && qfs->barcode_item != qfs->range_end_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both a range_item and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item != qfs->ordered_item1 && qfs->barcode_item != qfs->ordered_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an ordered_item and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item != qfs->seq_len_item, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both a seq_len_item and barcode_item",
                        qfs->name, qfs->barcode_item);
            }

            if (qfs->barcode_item2 != -1) {
                ASSERT (qfs->barcode_item2 != qfs->range_end_item1 && qfs->barcode_item2 != qfs->range_end_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both a range_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);

                ASSERT (qfs->barcode_item2 != qfs->ordered_item1 && qfs->barcode_item2 != qfs->ordered_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an ordered_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);

                ASSERT (qfs->barcode_item2 != qfs->seq_len_item, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both a seq_len_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);
            }

            ASSERT (qfs->ordered_item1 == -1 || qfs->is_integer[qfs->ordered_item1] || qfs->is_numeric[qfs->ordered_item1] || qfs->is_hex[qfs->ordered_item1], 
                    "Error in definition of QNAME flavor=%s: item=%u is one of \"ordered_item\" - expecting it to be is_integer or is_numeric or is_hex", qfs->name, qfs->ordered_item1);

            if (qfs->callback_item != -1) {
                ASSERT (qf_callbacks[qfs->id], "QNAME flavor=%s defines a callback, but the callback is not listed in qf_callbacks", qfs->name);

                ASSERT (qfs->callback_item != qfs->range_end_item1 && qfs->callback_item != qfs->range_end_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both a range_item and callback_item",
                        qfs->name, qfs->callback_item);

                ASSERT (qfs->callback_item != qfs->ordered_item1 && qfs->callback_item != qfs->ordered_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an ordered_item and callback_item",
                        qfs->name, qfs->callback_item);

                ASSERT (!qfs->is_in_local[qfs->callback_item], 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an is_in_local and callback_item",
                        qfs->name, qfs->callback_item);
            }
        }

        seg_prepare_snip_other (SNIP_COPY, (DictId)_SAM_QNAME, false, 0, copy_qname); // QNAME dict_id is the same for SAM, FASTQ, KRAKEN
    }
}

void qname_seg_initialize (VBlockP vb, QType q, Did st_did_i)
{
    #define ctx_by_item(item) ((did_q - SAM_QNAME) + ctx_get_ctx (vb, qfs->con.items[item].dict_id))
    QnameFlavor qfs = segconf.qname_flavor[q];
    if (!qfs) return;

    Did did_q = did_by_q (q);

    // consolidate all items under did_q in --stats
    if (st_did_i != DID_NONE)
        ctx_consolidate_statsN (vb, st_did_i, did_q, MAX_QNAME_ITEMS + 1); 

    // set STORE_INT as appropriate
    if (qfs->ordered_item1   != -1) ctx_by_item (qfs->ordered_item1)    ->flags.store = STORE_INT; 
    if (qfs->ordered_item2   != -1) ctx_by_item (qfs->ordered_item2)    ->flags.store = STORE_INT; 
    if (qfs->range_end_item1 != -1) ctx_by_item (qfs->range_end_item1-1)->flags.store = STORE_INT; 
    if (qfs->range_end_item2 != -1) ctx_by_item (qfs->range_end_item2-1)->flags.store = STORE_INT; 
    if (qfs->seq_len_item    != -1) ctx_by_item (qfs->seq_len_item)     ->flags.store = STORE_INT; 

    for (int i=0; qfs->in_local[i] != -1; i++) {
        #define IS(what) qfs->is_##what[qfs->in_local[i]]

        if (IS(integer) || IS(numeric))
            ctx_by_item (qfs->in_local[i])->ltype = IS(hex) ? LT_DYN_INT_h : LT_DYN_INT; // required by seg_integer_or_not / seg_numeric_or_not

        else // textual
            ctx_by_item (qfs->in_local[i])->no_stons = true; // required by seg_add_to_local_text
    }

    // a barcode that is not in local, is strictly in dict (eg molecular identifiers might be sparsely distributed in the file - we don't want the first VB encountering to mark as singleton)
    if (qfs->barcode_item != -1) 
        ctx_by_item (qfs->barcode_item)->no_stons = true;

    if (qfs->barcode_item2 != -1)
        ctx_by_item (qfs->barcode_item2)->no_stons = true;

    if (qfs->seq_len_item != -1 && (VB_DT(SAM) || VB_DT(BAM))) 
        ctx_by_item (qfs->seq_len_item)->flags.store_per_line = true; // consumed by sam_cigar_special_CIGAR

    // when pairing, we cannot have singletons, bc a singleton in R1, when appearing in R2 will not
    // be a singleton (since not the first appearance), causing both b250 and local to differ between the R's. 
    if (flag.pair)
        for (int i=0; i < qfs->con.nitems_lo; i++) 
            ctx_by_item (i)->no_stons = true;
}

// note: we run this function only in discovery, not in segging, because it is quite expensive - checking all numerics.
// returns 0 if qname is indeed the flavor, or error code if not

QnameTestResult qname_test_flavor (STRp(qname), QType q, QnameFlavor qfs, bool quiet) // out
{
    // return embedded qname2, eg "82a6ce63-eb2d-4812-ad19-136092a95f3d" in "@ERR3278978.1 82a6ce63-eb2d-4812-ad19-136092a95f3d/1"
    if (q != QANY && qfs->only_q != QANY && qfs->only_q != q
        && !(qfs->only_q == Q1or3   && (q == QNAME1 || q == QLINE3))
        && !(qfs->only_q == QSAM    && (Z_DT(BAM) || Z_DT(SAM)))
        && !(qfs->only_q == Q2orSAM && (q == QNAME2 || Z_DT(BAM) || Z_DT(SAM)))) return QTR_WRONG_Q; 

    if (q==QNAME2 && (segconf.tech != qfs->qname1_tech && segconf.tech != TECH_UNKNOWN))
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
    ASSERTNOTZERO (segconf.running);

    static rom reasons[] = QTR_NAME;

    segconf.qname_flavor[q] = NULL; // unknown

    for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {
        QnameTestResult reason = QTR_SUCCESS;
        if (!(reason = qname_test_flavor (STRa(qname), q, qfs, true))) {
            segconf.qname_flavor[q] = qfs;

            static QnameCNN char_to_cnn[256] = CHAR_TO_CNN;
            segconf.flav_prop[q] = (QnameFlavorProp){ .id          = qfs->id, 
                                                      .has_seq_len = qfs->seq_len_item != -1,
                                                      .has_R       = qfs->is_mated, 
                                                      .cnn         = char_to_cnn[(int)qfs->cut_to_canonize] };
            
            if (q == QNAME1 || (q == QNAME2 && segconf.tech == TECH_NCBI))
                segconf.tech = qfs->tech; // note: if this is QNAME2, we update tech according to QNAME2 (instead of NCBI)
            
            ASSERT (!qfs->cut_to_canonize || segconf.flav_prop[q].cnn, "flavor=%s has cut_to_canonize='%c', but it is missing in CHAR_TO_CNN", qfs->name, qfs->cut_to_canonize);

            ASSERT (qname_len <= SAM_MAX_QNAME_LEN, "qname=\"%.*s\" is too long (len=%u)", STRf(qname), qname_len);
            memcpy (segconf.qname_line0[q], qname, qname_len);

            if (qfs->seq_len_item != -1) 
                segconf.seq_len_dict_id = qfs->con.items[qfs->seq_len_item].dict_id;
            
            qname_seg_initialize (vb, q, DID_NONE); // so the rest of segconf.running can seg fast using the discovered container

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

    // set up dict id alias. need to do explicitly, because not predefined
    if (segconf.qname_flavor[q]->barcode_item2 != -1) {
        Did did_i_bc1 = did_by_q (q) + 1 + segconf.qname_flavor[q]->barcode_item;
        Did did_i_bc2 = did_by_q (q) + 1 + segconf.qname_flavor[q]->barcode_item2;
        ZCTX(did_i_bc2)->dict_did_i = did_i_bc1;
    }

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

// attempts to modify flavor to accommodate both qname of line_i==0 and current qname. return true if flavor has been modified.
static bool qname_segconf_rediscover_flavor (VBlockP vb, QType q, STRp(qname))
{
    ASSERTNOTZERO (segconf.running); // one of the reasons this cannot run outside of segconf is that segconf_seq_len_dict_id must be the same for all lines

    if (segconf.qname_flavor_rediscovered[q]) return false; // we already rediscovered once, not going to do it again

    if (flag.debug_qname) iprintf ("%s: Rediscovering %s\n", LN_NAME, qtype_name (q));

    QnameFlavor save_flavor = segconf.qname_flavor[q];
    qname_segconf_discover_flavor (vb, q, STRa(qname));

    if (!segconf.qname_flavor[q]) goto fail; // qname has unknown flavor

    // test if qname of line_i=0 can agree on the new flavor
    QnameTestResult reason = qname_test_flavor (segconf.qname_line0[q], strlen(segconf.qname_line0[q]), QANY, segconf.qname_flavor[q], !flag.debug_qname); 
    if (reason != QTR_SUCCESS) goto fail; // qname of line_i=0 doesn't agree with new flavor, we will stick to the old flavor

    // we found a new flavor which is agreeable to both this qname and line_i=0's qname
    segconf.qname_flavor_rediscovered[q] = true;
    return true;

fail:
    segconf.qname_flavor[q] = save_flavor;
    return false;
}

// attempt to seg according to the qf - return true if successful
static bool qname_seg_qf (VBlockP vb, QType q, STRp(qname), unsigned add_additional_bytes)
{
    QnameFlavor qfs = segconf.qname_flavor[q];
    ASSERTNOTNULL (qfs);
    
    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item, NULL);
    if (!n_items) return false; // failed to split according to container - qname is not of this flavor

    // seg
    int64_t value, prev_value;
    ContextP qname_ctx = CTX(did_by_q(q));
    
    // seg container
    seg_by_ctx (vb, STRi(qfs->con_snip, q), qname_ctx, qfs->num_seps + add_additional_bytes); // account for container separators, prefixes and caller-requested add_additional_bytes 

    bool sorted_by_qname = VB_DT(FASTQ) || VB_DT(KRAKEN) ||   
                           segconf.is_collated || qfs->sam_qname_sorted || 
                           (!segconf.is_sorted && segconf.is_long_reads);

    for (unsigned item_i=0; item_i < qfs->con.nitems_lo; item_i++) {

        // get item ctx - Q?NAME did_i are immediately following QNAME, and QmNAME si the last.
        const ContainerItem *item = &qfs->con.items[item_i];
        rom str = items[item_i];
        uint32_t str_len = item_lens[item_i];
        value = -1;

        // if not integer/numeric as expected, send to rediscovery. 
        // However, if we can no longer rediscover, try to accommodate rather than failing
        if (segconf.running && !segconf.qname_flavor_rediscovered[q]) {
            if (qfs->is_integer[item_i] && !str_is_int (STRa(str))) return false;
            if (qfs->is_numeric[item_i] && !str_is_numeric (STRa(str))) return false;
        }

        // calculate context - a bit messy, but faster than looking up by dict_id (relying on Q?NAME being sequential in the container)
        ContextP item_ctx = qname_ctx + ((item->dict_id.num == _SAM_QmNAME) ? MAX_QNAME_ITEMS : (1+item_i)); // note: QmName, if exists, is always the last item
                        
        // case: this is the file is sorted by qname - delta against previous
        if (sorted_by_qname && 
            (item_i == qfs->ordered_item1 || item_i == qfs->ordered_item2) &&
            (ctx_has_value_in_prev_line_(vb, item_ctx) || vb->line_i==0) &&
            ( (!qfs->is_hex[item_i] && str_get_int_dec (STRa(str), (uint64_t*)&value)) || ( qfs->is_hex[item_i] && str_get_int_hex (STRa(str), true, false, (uint64_t*)&value))) && // lower-case hex
            (ABS(value - item_ctx->last_value.i) < MAX_TOKENIZER_DETLA)) 
            seg_self_delta (vb, item_ctx, value, 
                            (qfs->is_hex[item_i]?'x' : qfs->is_numeric[item_i]?'d' : 0), 
                            (qfs->is_hex[item_i] || qfs->is_numeric[item_i]) ? str_len : 0,
                            str_len);

        // case: end-of-range item, seg as delta vs previous item which is start-of-range
        else if ((qfs->range_end_item1 == item_i || qfs->range_end_item2 == item_i) &&
                 str_get_int (STRa(str), &value) && 
                 str_get_int (STRi(item, item_i-1), &prev_value)) {
            SNIP(32);
            seg_prepare_snip_other (SNIP_OTHER_DELTA, qfs->con.items[item_i-1].dict_id, true, value - prev_value, snip);
            seg_by_ctx (vb, STRa(snip), item_ctx, str_len);      
        }

        else if (qfs->is_in_local[item_i]) { 
            if (qfs->is_integer[item_i]) 
                seg_integer_or_not (vb, item_ctx, STRa(str), str_len); // hex or not

            else if (qfs->is_numeric[item_i])
                seg_numeric_or_not (vb, item_ctx, STRa(str), str_len); // hex or not

            else 
                seg_add_to_local_text (vb, item_ctx, STRa(str), LOOKUP_SIMPLE, str_len);
        }

        else if (qfs->callback_item == item_i)
            qf_callbacks[qfs->id](vb, item_ctx, STRa(str));

        else if (item->separator[0] == CI0_SKIP)
            {} // no segging a skipped item

        // case: textual item
        else 
            seg_by_ctx (vb, STRa(str), item_ctx, str_len); 

        // case: this item is qname_seq_len - set last_value by the beneficial field (CIGAR in SAM, ? in FASTQ)
        if (item_i == qfs->seq_len_item && (value >= 0 || str_get_int (STRa(str), &value))) 
            ctx_set_last_value (vb, item_ctx, value);  
    }

    if (segconf.running && segconf.flav_prop[q].has_R && qname[qname_len-1] != '1' && qname[qname_len-1] != '2')
        segconf.flav_prop[q].has_R = false;

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

    // if we're in segconf - check if there is another flavor that matches both this qname and 
    // the qname of line_i=0 on which flavor is based. Eg BGI-varlen can change to BGI-R7 if this qname has a leading 0 in Q4NAME.
    if (!success && segconf.running && qname_segconf_rediscover_flavor (vb, q, STRa(qname)))
        success = qname_seg_qf (vb, q, STRa(qname), add_additional_bytes); // now with new flavor

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

// reduces qname to its canonical form: possibly reduces qname_len to make a qname more likely compareble between SAM/BAM and FASTQ 
void qname_canonize (QType q, rom qname, uint32_t *qname_len)
{
    QnameFlavorProp *f = &segconf.flav_prop[q];

    // mated: "HSQ1004:134:C0D8DACXX:3:1101:1318:114841/2" ⟶ "HSQ1004:134:C0D8DACXX:3:1101:1318:114841"
    // SRA2:  "ERR2708427.177.1" ⟶ "ERR2708427.177"
    if (f->has_R && *qname_len > 2) 
        *qname_len -= 2;

    // eg: "NOVID_3053_FC625AGAAXX:6:1:1069:11483:0,84" ⟶ "NOVID_3053_FC625AGAAXX:6:1:1069:11483"
    else if (f->cnn) {
        static char cnn_to_char[NUM_CNN] = CNN_TO_CHAR;
        rom cut = memrchr (qname + 1, cnn_to_char[f->cnn], *qname_len); // +1 so at least one character survives
        if (cut) *qname_len = cut - qname;
    }
}

uint32_t qname_calc_hash (QType q, STRp(qname), thool is_last, bool canonical, 
                          uint32_t *uncanonical_suffix_len) // optional out
{
    uint32_t save_qame_len = qname_len;
    if (canonical) qname_canonize (q, qSTRa(qname));

    if (uncanonical_suffix_len) 
        *uncanonical_suffix_len = save_qame_len - qname_len;
        
    if (!qname_len) return 0;

    uint8_t data[qname_len];
    for (int i=0; i < qname_len; i++)
        data[i] = ((((uint8_t *)qname)[i]-33) & 0x3f) | ((((uint8_t *)qname)[qname_len-i-1] & 0x3) << 6);

    uint32_t hash = crc32 (0, data, qname_len);
    if (is_last != unknown) hash = (hash & 0xfffffffe) | is_last;

    return hash;
}

rom segconf_qf_name (QType q)
{
    if (q == QSAM)
        return segconf.deep_sam_qname_flavor ? segconf.deep_sam_qname_flavor->name : "N/A";
    else
        return segconf.qname_flavor[q] ? segconf.qname_flavor[q]->name : "N/A";
}

DictIdAlias qname_get_alias (QType q)
{
    int bc2;
    Did first_did_i = did_by_q (q) + 1; // eg Q0NAME 

    if (!segconf.qname_flavor[q] || (bc2 = segconf.qname_flavor[q]->barcode_item2) == -1)
        return (DictIdAlias){};

    else 
        return (DictIdAlias){ .alias_type = ALIAS_DICT,
                              .alias      = dt_fields[DT_FASTQ].predefined[first_did_i + bc2].dict_id,
                              .dst        = dt_fields[DT_FASTQ].predefined[first_did_i + bc2 - 1].dict_id }; 
}

QnameFlavorId segconf_qf_id (QType q)
{
    return segconf.qname_flavor[q] ? segconf.qname_flavor[q]->id : QF_NO_ID;
}

rom qtype_name (QType q)
{
    static rom qtype_names[] = QTYPE_NAME, qtype_neg_names[] = QTYPE_NEG_NAMES;

    return (q >= 0 && q < ARRAY_LEN (qtype_names))     ? qtype_names[q] 
          :(q < 0 && -q < ARRAY_LEN (qtype_neg_names)) ? qtype_neg_names[-q]
          :                                              "Invalid_qtype";
}

