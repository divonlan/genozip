// ------------------------------------------------------------------
//   qname.c
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "libdeflate_1.19/libdeflate.h"
#include "qname.h"
#include "tokenizer.h"
#include "seg.h"
#include "stats.h"
#include "qname_flavors.h"
#include "file.h"
#include "codec.h"
#include "zip_dyn_int.h"

sSTRl(copy_qname, 16);
sSTRl(snip_redirect_to_QNAME2, 16);

QnameFlavor qname_get_optimize_qf (void)
{
    bool is_mated = segconf.qname_flavor[QNAME1] && qf_is_mated(QNAME1);

    return &qf[NUM_QFs-1 - is_mated]; // relying on Genozip-opt being last
}

static inline Did did_by_q (QType q) 
{
    ASSERTINRANGE (q, 0, NUM_QTYPES);

    return (Did[]){ FASTQ_QNAME, FASTQ_QNAME2, FASTQ_LINE3 }[q]; // note: SAM and FASTQ have the same dids for QNAMEs
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
    strcpy (&qfs->name[strlen(qfs->name)], "@");
    
    // examples
    for (int i=0; i < QFS_MAX_EXAMPLES; i++) {
        unsigned example_len = strlen(qfs->example[i]);
        if (!example_len) break;

        strcpy (&qfs->example[i][example_len], "/1");
    }

    // case 1: previous item is not fixed of invisible - move its separator (possibly 0) to the mate item and make it '/'
    uint8_t *prev_sep = qfs->con.items[mate_item_i-1].separator; 
    uint8_t *mate_sep = qfs->con.items[mate_item_i  ].separator; 
    
    if (prev_sep[0] != CI0_FIXED_0_PAD && prev_sep[0] != CI0_VAR_0_PAD) {
        mate_sep[0] = prev_sep[0]; // 0 or ' '
        mate_sep[1] = 0;
        prev_sep[0] = '/';
    }

    // case 2: previous item is fixed - add / as a prefix. 
    else { // eg con_roche_454, con_pacbio_range
        // qfs must already a PX_MATE_FIXED_0_PAD prefix item for the mate item, if previous item is fixed or invisible
        ASSERT (qname_get_px_str_len (qfs) > mate_item_i, "QnameFlavor=%s: expecting PX_MATE prefix item to exist for mate_item_i=%u",
                qfs->name, mate_item_i);

        qfs->px_strs[mate_item_i] = "/";
    }

    // add new item as numeric item (at beginning of array - easier, and the order doesn't matter)
    if (mate_sep[0] == CI0_FIXED_0_PAD || mate_sep[0] == CI0_VAR_0_PAD) {
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

                char item_sep = qfs->con.items[qfs->con.nitems_lo-1].separator[0];
                ASSERT (item_sep == CI0_SKIP, "In flavor=\"%s\", expecting separator of %s to be I_AM_MATE since it is mated, but it is '%c'(%u)", 
                        qfs->name, dis_dict_id (qfs->con.items[qfs->con.nitems_lo-1].dict_id).s, item_sep, item_sep);

                ASSERT (qfs->con.nitems_lo < 2 || 
                        (qfs->con.items[qfs->con.nitems_lo-2].separator[0] != CI0_FIXED_0_PAD && qfs->con.items[qfs->con.nitems_lo-2].separator[0] != CI0_VAR_0_PAD) || 
                        (qfs->px_strs[qfs->con.nitems_lo-1] && qfs->px_strs[qfs->con.nitems_lo-1][0] == CI0_SKIP), 
                        "In flavor=\"%s\", since %s has separator CI0_FIXED_0_PAD or CI0_VAR_0_PAD, prefix[%u] must be PX_MATE_FIXED_0_PAD, but it is \"%s\"", 
                        qfs->name, dis_dict_id (qfs->con.items[qfs->con.nitems_lo-2].dict_id).s, qfs->con.nitems_lo-1,
                        qfs->px_strs[qfs->con.nitems_lo-1] ? qfs->px_strs[qfs->con.nitems_lo-1] : "(NULL)");

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

                // verify Q*NAME requirements
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
                    unsigned len = strnlen (*px, 100);
                    ASSERT (len < 100, "px not nul-terminated in qfs=%s", qfs->name);

                    ASSERT (qfs->con_prefix_len + len + 1 <= MAX_PREFIX_LEN, "Prefix too long in qfs=%s: item=%d qfs->con_prefix_len(after)=%d",
                            qfs->name, (int)(px - qfs->px_strs), qfs->con_prefix_len + len + 1);

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
                
                for_con2 (&con_no_skip)
                    if (item->dict_id.num != _SAM_QmNAME) 
                        item->dict_id = dt_fields[DT_FASTQ].predefined[first_did_i + item_i + 1].dict_id; 
                    else
                        item->dict_id = dt_fields[DT_FASTQ].predefined[first_did_i + MAX_QNAME_ITEMS].dict_id; 

                // prepare container snip for this QF
                qfs->con_snip_lens[q] = sizeof (qfs->con_snips[0]);            
                container_prepare_snip ((Container*)&con_no_skip, qfs->con_prefix, prefix_no_skip_len, qSTRi(qfs->con_snip,q));
            }

            // in qname_flavors.h, we keep lists in the form of index lists, for maintenanbility. now
            // we convert them to a bitmap for ease of segging.
            for (unsigned i=0; qfs->integer_items[i] != -1; i++) {
                ASSERT (qfs->con.items[qfs->integer_items[i]].separator[0] != CI0_FIXED_0_PAD && 
                        qfs->con.items[qfs->integer_items[i]].separator[0] != CI0_VAR_0_PAD,
                        "since qfs=%s item=%u has seperator[0] ∈ { CI0_FIXED_0_PAD, CI0_VAR_0_PAD }, expecting it to be numeric, not int", qfs->name, qfs->integer_items[i]);
                qfs->is_integer[qfs->integer_items[i]] = true;

                ASSERT (qfs->barcode_item != qfs->integer_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an integer_item and barcode_item",
                        qfs->name, qfs->barcode_item);

                ASSERT (qfs->barcode_item2 != qfs->integer_items[i], "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an integer_item and barcode_item2",
                        qfs->name, qfs->barcode_item2);
            }

            for (unsigned i=0; qfs->numeric_items[i] != -1; i++) {
                ASSERT (qfs->con.items[qfs->numeric_items[i]].separator[0] == CI0_FIXED_0_PAD || 
                        qfs->con.items[qfs->numeric_items[i]].separator[0] == CI0_VAR_0_PAD,
                        "since qfs=%s item=%u has seperator[0] ∉ { CI0_FIXED_0_PAD, CI0_VAR_0_PAD }, expecting it to be int, not numeric", qfs->name, qfs->numeric_items[i]);
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

            if (qfs->callback_item2 != -1) {
                ASSERT (qf_callbacks[qfs->id], "QNAME flavor=%s defines a callback, but the callback is not listed in qf_callbacks", qfs->name);

                ASSERT (qfs->callback_item2 != qfs->range_end_item1 && qfs->callback_item2 != qfs->range_end_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both a range_item and callback_item2",
                        qfs->name, qfs->callback_item2);

                ASSERT (qfs->callback_item2 != qfs->ordered_item1 && qfs->callback_item2 != qfs->ordered_item2, 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an ordered_item and callback_item2",
                        qfs->name, qfs->callback_item2);

                ASSERT (!qfs->is_in_local[qfs->callback_item2], 
                        "Error in definition of QNAME flavor=%s: item=%u is invalidly defined as both an is_in_local and callback_item2",
                        qfs->name, qfs->callback_item2);
            }
        }

        seg_prepare_snip_other (SNIP_COPY, _SAM_QNAME,  false, 0, copy_qname); // QNAME dict_id is the same for SAM, FASTQ
        seg_prepare_snip_other (SNIP_COPY, _SAM_Q5NAME, false, 0, copy_q5name);
        seg_prepare_snip_other (SNIP_COPY, _SAM_Q6NAME, false, 0, copy_q6name);

        seg_prepare_snip_other (SNIP_REDIRECTION, _SAM_QNAME2, false, 0, snip_redirect_to_QNAME2);

        tokenizer_zip_initialize();
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

    // set STORE_INT and ltype as appropriate
    if (segconf.sorted_by_qname[q]) {
        #define INIT_ORDERED(field) ({ ContextP ctx = ctx_by_item (qfs->field); ctx->flags.store = STORE_INT; dyn_int_init_ctx (vb, ctx, 0); })
        if (qfs->ordered_item1 != -1) INIT_ORDERED (ordered_item1);
        if (qfs->ordered_item2 != -1) INIT_ORDERED (ordered_item2);
    }

    #define INIT_RANGE(field) ({ ContextP ctx = ctx_by_item (qfs->field); (ctx-1)->flags.store = ctx->flags.store = STORE_INT; dyn_int_init_ctx (vb, ctx, 0); })
    if (qfs->range_end_item1 != -1) INIT_RANGE (range_end_item1);
    if (qfs->range_end_item2 != -1) INIT_RANGE (range_end_item2);

    if (qfs->seq_len_item != -1) ctx_by_item (qfs->seq_len_item)->flags.store = STORE_INT; 

    for (int i=0; qfs->in_local[i] != -1; i++) {
        #define IS(what) qfs->is_##what[qfs->in_local[i]]
        ContextP item_ctx = ctx_by_item (qfs->in_local[i]);

        if (IS(integer) || IS(numeric)) {
            dyn_int_init_ctx (vb, item_ctx, 0);  // required by seg_integer_or_not / seg_numeric_or_not
            if (IS(hex)) item_ctx->ltype = LT_DYN_INT_h;
        }
        
        else // textual
            item_ctx->ltype = LT_STRING;
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

void qname_segconf_finalize (VBlockP vb)
{
    // SAM: if we have a consensus flavor, it must QNAME2, expected by sam_piz_con_item_cb
    if (VB_DT(SAM) || VB_DT(BAM)) {
        if (segconf.flav_prop[QNAME1].is_consensus && segconf.qname_flavor[QNAME2]) {
            SWAP (segconf.qname_flavor[QNAME1], segconf.qname_flavor[QNAME2]);
            SWAP (segconf.flav_prop[QNAME1],    segconf.flav_prop[QNAME2]);
            SWAP (segconf.qname_line0[QNAME1],  segconf.qname_line0[QNAME2]);
        }

        if (flag.deep) 
            for (QType q=QNAME1; q <= QNAME2; q++)
                segconf.deep_sam_qname_flavor[q] = segconf.qname_flavor[q];
    }

    for (QType q=QNAME1; q < NUM_QTYPES; q++) 
        if (segconf.qname_flavor[q])
            segconf.sorted_by_qname[q] = VB_DT(FASTQ) ||
                                         segconf.is_collated || segconf.sam_is_unmapped || segconf.qname_flavor[q]->sam_qname_sorted || 
                                         (!segconf.is_sorted && !segconf.is_paired);
            
    // if only consensus reads exist, or is unknown
    if (segconf.tech == TECH_CONS || segconf.tech == TECH_UNKNOWN)
        segconf.tech = segconf.tech_by_RG; // usually tech by QNAME is more reliable, but if not available, go by the PL field of @RG in the SAM header (which might also be unavailable)
}

// note: we run this function only in discovery, not in segging, because it is quite expensive - checking all numerics.
// returns 0 if qname is indeed the flavor, or error code if not

QnameTestResult qname_test_flavor (STRp(qname), QType q, QnameFlavor qfs, bool quiet) // out
{
    ASSERTNOTNULL (qfs);

    // return embedded qname2, eg "82a6ce63-eb2d-4812-ad19-136092a95f3d" in "@ERR3278978.1 82a6ce63-eb2d-4812-ad19-136092a95f3d/1"
    if (q != QANY && qfs->only_q != QANY && qfs->only_q != q
        && !(qfs->only_q == Q1or3   && (q == QNAME1 || q == QLINE3))
        && !(qfs->only_q == QSAM    && (Z_DT(BAM) || Z_DT(SAM)))
        && !(qfs->only_q == Q2orSAM && (q == QNAME2 || Z_DT(BAM) || Z_DT(SAM)))) 
        return QTR_WRONG_Q; 

    // in FASTQ, QNAME1 match fq_qname1_tech
    if (q==QNAME2 && (Z_DT(FASTQ) || Z_DT(FASTA)) && 
        qfs->fq_qname1_tech != TECH_ANY && segconf.tech != qfs->fq_qname1_tech && segconf.tech != TECH_UNKNOWN)
            return QTR_TECH_MISMATCH;

    if (!qname_len) 
        return QTR_QNAME_LEN_0;

    // test fixed length, if applicable
    if (qfs->fixed_len && qname_len != qfs->fixed_len) 
        return QTR_FIXED_LEN_MISMATCH;

    if (qfs->is_mated && qname[qname_len-1] != '1' && qname[qname_len-1] != '2') 
        return QTR_NO_MATE; // quick check

    str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item, quiet ? NULL : qfs->name);
    if (!n_items) 
        return QTR_CONTAINER_MISMATCH; // failed to split according to container - qname is not of this flavor

    // items cannot include whitespace or non-printable chars
    for (uint32_t item_i=0; item_i < n_items; item_i++)
        if (!str_is_no_ws (STRi(item, item_i))) 
            return QTR_BAD_CHARS;

    // check that all the items expected to be numeric (leading zeros ok) are indeed so
    for (const int *item_i = qfs->numeric_items; *item_i != -1; item_i++)
        if (!qfs->is_hex[*item_i] && !str_is_numeric (STRi(item, *item_i))) 
            return QTR_BAD_NUMERIC;

    // check that all the items expected to be lower case hexadecimal are indeed so
    for (const int *item_i = qfs->hex_items; *item_i != -1; item_i++)
        if (!str_is_hexlo (STRi(item, *item_i))) 
            return QTR_BAD_HEX;

    // check that all the items expected to be integer (no leading zeros) are indeed so
    for (const int *item_i = qfs->integer_items; *item_i != -1; item_i++)
        if (!qfs->is_hex[*item_i] && !str_is_int (STRi(item, *item_i))) 
            return QTR_BAD_INTEGER;

    if (qfs->is_mated && (item_lens[n_items-1] != 1 || (*items[n_items-1] != '1' && *items[n_items-1] != '2')))
        return QTR_BAD_MATE; // mate is expected to be "1" or "2" - better check than NO_MATE - after splitting
    
    if (!qfs->validate_flavor (STRas(item)))
        return QTR_FAILED_VALIDATE_FUNC;

    return QTR_SUCCESS; // yes, qname is of this flavor
}

// called for the first line in segconf.running 
void qname_segconf_discover_flavor (VBlockP vb, QType q, STRp(qname))
{
    static rom reasons[] = QTR_NAME;

    segconf.qname_flavor[q] = NULL; // unknown

    for (QnameFlavor qfs=&qf[0]; qfs < &qf[NUM_QFs]; qfs++) {
        QnameTestResult reason = QTR_SUCCESS;
        if (!(reason = qname_test_flavor (STRa(qname), q, qfs, true))) {

            // case: first discovery: store discovered qname
            if (!segconf.qname_line0[q].s[0]) 
                memcpy (segconf.qname_line0[q].s, qname, MIN_(qname_len, SAM_MAX_QNAME_LEN));

            // case: rediscovery: check that it is also an acceptable flavor for stored qname_line0
            else {
                if (QTR_SUCCESS != qname_test_flavor (segconf.qname_line0[q].s, strlen(segconf.qname_line0[q].s), QANY, qfs, true))
                    continue;
            }

            segconf.qname_flavor[q] = qfs;

            static QnameCNN char_to_cnn[256] = CHAR_TO_CNN;
            segconf.flav_prop[q] = (QnameFlavorProp){ .has_seq_len  = (qfs->seq_len_item != -1),
                                                      .is_mated     = qfs->is_mated, 
                                                      .is_consensus = (qfs->tech == TECH_CONS),
                                                      .cnn          = char_to_cnn[(int)qfs->cut_to_canonize] };

            if (q == QNAME1 || (q == QNAME2 && (segconf.tech == TECH_NCBI || segconf.tech == TECH_CONS)))
                segconf.tech = qfs->tech; // note: if this is QNAME2, we update tech according to QNAME2 (instead of NCBI)
            
            ASSERT (!qfs->cut_to_canonize || segconf.flav_prop[q].cnn, "flavor=%s has cut_to_canonize='%c', but it is missing in CHAR_TO_CNN", qfs->name, qfs->cut_to_canonize);

            if (qfs->seq_len_item != -1) 
                segconf.seq_len_dict_id = qfs->con.items[qfs->seq_len_item].dict_id;
            
            qname_seg_initialize (vb, q, DID_NONE); // so the rest of segconf.running can seg fast using the discovered container

            if (flag.debug_qname) 
                iprintf ("%.*s is DISCOVERED as %s - for %s\n", STRf(qname), qfs->name, qtype_name(q));
            
            break;
        }

        else if (flag.debug_qname && reason) 
            iprintf ("%.*s is not %s flavor \"%s\". Reason: %s\n", STRf(qname), qtype_name(q), qfs->name, reasons[reason]);
    }

    if (!segconf.qname_flavor[q]) {
        segconf.flav_prop[q].is_tokenized = true; // since 15.0.63

        if (flag.debug_qname)         
            iprintf ("%.*s - flavor is NOT DISCOVERED - for %s\n", STRf(qname), qtype_name(q));
    }

    // set up dict id alias. need to do explicitly, because not predefined
    if (segconf.qname_flavor[q] && segconf.qname_flavor[q]->barcode_item2 != -1) {
        Did did_i_bc1 = did_by_q (q) + 1 + segconf.qname_flavor[q]->barcode_item;
        Did did_i_bc2 = did_by_q (q) + 1 + segconf.qname_flavor[q]->barcode_item2;
        ZCTX(did_i_bc2)->dict_did_i = did_i_bc1;
    }

    if (segconf.qname_flavor[q] && segconf.qname_flavor[q]->id == QF_PACBIO_lbl && qname_len > 4 && !memcmp (&qname[qname_len-4], "/ccs", 4))
        segconf.is_pacbio_ccs = true;
        
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
bool qname_segconf_rediscover_flavor (VBlockP vb, QType q, STRp(qname))
{
    ASSERTNOTZERO (segconf_running); // one of the reasons this cannot run outside of segconf is that segconf_seq_len_dict_id must be the same for all lines

    if (segconf.qname_flavor_rediscovered[q]) return false; // we already rediscovered once, not going to do it again

    if (flag.debug_qname) iprintf ("%s: Rediscovering %s\n", LN_NAME, qtype_name (q));

    QnameFlavor save_flavor = segconf.qname_flavor[q];
    qname_segconf_discover_flavor (vb, q, STRa(qname));

    if (!segconf.qname_flavor[q]) goto fail; // qname has unknown flavor

    // we found a new flavor which is agreeable to both this qname and line_i=0's qname
    segconf.qname_flavor_rediscovered[q] = true;
    return true;

fail:
    segconf.qname_flavor[q] = save_flavor;
    return false;
}

// SAM seg
QType qname_sam_get_qtype (STRp(qname))
{
    // if we have a QNAME2 flavor and qname matches it, return QNAME2, otherwise QNAME1
    if (segconf.qname_flavor[QNAME2]) {
        QnameFlavor qfs = segconf.qname_flavor[QNAME2];
        str_split_by_container (qname, qname_len, &qfs->con, qfs->con_prefix, qfs->con_prefix_len, item, NULL);

        return n_items ? QNAME2 : QNAME1;
    }
    else
        return QNAME1;
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

    for (unsigned item_i=0; item_i < qfs->con.nitems_lo; item_i++) {

        // get item ctx - Q?NAME did_i are immediately following QNAME, and QmNAME si the last.
        const ContainerItem *item = &qfs->con.items[item_i];
        rom str = items[item_i];
        uint32_t str_len = item_lens[item_i];
        value = -1;

        // if not integer/numeric as expected, send to rediscovery. 
        // However, if we can no longer rediscover, try to accommodate rather than failing
        if (segconf_running && !segconf.qname_flavor_rediscovered[q]) {
            if (qfs->is_integer[item_i] && !str_is_int (STRa(str))) return false;
            if (qfs->is_numeric[item_i] && !str_is_numeric (STRa(str))) return false;
        }

        // calculate context - a bit messy, but faster than looking up by dict_id (relying on Q?NAME being sequential in the container)
        ContextP item_ctx = qname_ctx + ((item->dict_id.num == _SAM_QmNAME) ? MAX_QNAME_ITEMS : (1+item_i)); // note: QmName, if exists, is always the last item
                        
        // case: this is the file is sorted by qname - delta against previous
        if (segconf.sorted_by_qname[q] && 
            (item_i == qfs->ordered_item1 || item_i == qfs->ordered_item2) &&
            ( (!qfs->is_hex[item_i] && qfs->is_integer[item_i] && str_get_int     (STRa(str), &value)) || 
              (!qfs->is_hex[item_i] && qfs->is_numeric[item_i] && str_get_int_dec (STRa(str), (uint64_t*)&value)) || 
              ( qfs->is_hex[item_i]                            && str_get_int_hex (STRa(str), true, false, (uint64_t*)&value)))) // lower-case hex

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
                seg_add_to_local_string (vb, item_ctx, STRa(str), LOOKUP_SIMPLE, str_len);
        }

        else if (qfs->callback_item == item_i || qfs->callback_item2 == item_i)
            qf_callbacks[qfs->id](vb, item_ctx, STRa(str));

        else if (item->separator[0] == CI0_SKIP)
            {} // no segging a skipped item

        // case: textual item
        else
            seg_by_ctx (vb, STRa(str), item_ctx, str_len); 

        // case: set last_value if needed and not already set. note: some lines do not set last_value -
        // if qname is not parsed and just copied from elsewhere (previous line, buddy, pair, deep, sag...)
        if (item_ctx->flags.store == STORE_INT && !ctx_has_value_in_line_(vb, item_ctx) &&
            (value >= 0 || str_get_int (STRa(str), &value)))
            ctx_set_last_value (vb, item_ctx, value); 

        set_last_txtC (item_ctx, str, str_len); 
    }

    return true;
}

// returns true is redirected to QNAME2
bool qname_seg (VBlockP vb, QType q, STRp (qname), unsigned add_additional_bytes)  // account for characters in addition to the field
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
    // the qname of line_i=0 on which flavor is based. Eg MGI-varlen can change to MGI-R7 if this qname has a leading 0 in Q4NAME.
    if (!success && segconf_running && qname_segconf_rediscover_flavor (vb, q, STRa(qname)))
        success = qname_seg_qf (vb, q, STRa(qname), add_additional_bytes); // now with new flavor

    // in SAM/BAM, we allow two flavours - try the second one (eg the second could be a consensus read name)
    if (!success && flavor && (VB_DT(SAM) || VB_DT(BAM)) && q == QNAME1) {
        qname_seg (vb, QNAME2, STRa(qname), 0);

        seg_by_ctx (vb, STRa(snip_redirect_to_QNAME2), qname_ctx, add_additional_bytes);
        return true;
    } 

    if (!success) { 
        tokenizer_seg (VB, qname_ctx, STRa(qname), 
                       VB_DT(FASTQ) ? sep_with_space : sep_without_space, 
                       add_additional_bytes);

        if (segconf_running) {
            // case SAM/BAM: reset to QNAME1 to avoid a confusing stats display - as we failed to seg by QNAME1 and QNAME2 flavors - but this is the one and only QNAME field
            if (VB_DT(SAM) || VB_DT(BAM)) q = QNAME1; 

            // collect the first 6 qnames / q2names, if flavor is unknown OR not successful in segging by flavor 
            if (segconf.n_1st_flav_qnames[q] < NUM_COLLECTED_WORDS) { // unrecognized flavor
                memcpy (segconf.unk_flav_qnames[q][segconf.n_1st_flav_qnames[q]], qname, MIN_(qname_len, UNK_QNANE_LEN));
        
                if (!segconf.n_1st_flav_qnames[q] || memcmp (segconf.unk_flav_qnames[q][segconf.n_1st_flav_qnames[q]], segconf.unk_flav_qnames[q][segconf.n_1st_flav_qnames[q]-1], UNK_QNANE_LEN))
                    segconf.n_1st_flav_qnames[q]++; // advance iff qname is different than previous line (not always the case if collated)
            }
        }
    }

done:
    seg_set_last_txt (vb, qname_ctx, STRa(qname)); // needed for bi:Z, parent_read_id...
    COPY_TIMER (qname_seg);

    return false;
}

// reduces qname to its canonical form: possibly reduces qname_len to make a qname more likely compareble between SAM/BAM and FASTQ 
void qname_canonize (QType q, rom qname, uint32_t *qname_len, CompIType comp_i/*only use in PIZ*/)
{
    QnameFlavorProp f = IS_ZIP                                 ? segconf.flav_prop[q] 
                      : (Z_DT(SAM) && comp_i <= SAM_COMP_DEPN) ? z_file->flav_prop[SAM_COMP_MAIN][q] // DEPN and PRIM don't have a TxtHeader section
                      :                                          z_file->flav_prop[comp_i][q];

    if (!f.is_tokenized) {
        // mated: "HSQ1004:134:C0D8DACXX:3:1101:1318:114841/2" ⟶ "HSQ1004:134:C0D8DACXX:3:1101:1318:114841"
        // SRA2:  "ERR2708427.177.1" ⟶ "ERR2708427.177"
        if (f.is_mated && *qname_len > 2) 
            *qname_len -= 2;

        // eg: "NOVID_3053_FC625AGAAXX:6:1:1069:11483:0,84" ⟶ "NOVID_3053_FC625AGAAXX:6:1:1069:11483"
        // note: it is possible that QNAME has both the mate and the cnn removed
        if (f.cnn) {
            static char cnn_to_char[NUM_CNN] = CNN_TO_CHAR;
            rom cut = memrchr (qname + 1, cnn_to_char[f.cnn], *qname_len); // +1 so at least one character survives
            if (cut) *qname_len = cut - qname;
        }
    }
    else { // not a recognized flavor - try our best (added 15.0.63)
        SAFE_NUL(&qname[*qname_len]);
        *qname_len = strcspn (qname, " \t\n\r");
        SAFE_RESTORE;

        if (*qname_len > 2 && qname[*qname_len-2] == '/' && 
            (qname[*qname_len-1] == '1' || qname[*qname_len-1] == '2'))
            *qname_len -= 2;
    }
}

uint32_t qname_hash_change_last (uint32_t hash, bool is_last)
{
    return (hash & 0xfffffffe) | is_last;
}

// ZIP: used for syncing DEPN/PRIM in gencomp and SAM/FASTQ in deep, and for scanning a SAM/BAM  
// PIZ: used for qname filters.
uint64_t qname_calc_hash (QType q, CompIType comp_i/*only used in PIZ*/, STRp(qname), thool is_last, bool canonical, CrcType type, 
                          uint32_t *uncanonical_suffix_len) // optional out
{
    uint32_t save_qame_len = qname_len;
    if (canonical) qname_canonize (q, qSTRa(qname), comp_i); // might shorten qname_len, does not modify the qname string

    if (uncanonical_suffix_len) 
        *uncanonical_suffix_len = save_qame_len - qname_len;
        
    if (!qname_len) return 0;

    uint8_t data[qname_len];
    for (int i=0; i < qname_len; i++)
        data[i] = ((((uint8_t *)qname)[i]-33) & 0x3f) | ((((uint8_t *)qname)[qname_len-i-1] & 0x3) << 6);

    uint64_t hash = (type == CRC64) ? crc64 (0, data, qname_len) : crc32 (0, data, qname_len);
    if (is_last != unknown) hash = (hash & 0xfffffffffffffffe) | is_last;

    return hash;
}

rom segconf_qf_name (QType q)
{
    if (q == QSAM)
        return segconf.deep_sam_qname_flavor[0] ? segconf.deep_sam_qname_flavor[0]->name : "N/A";
    else if (q == QSAM2)
        return segconf.deep_sam_qname_flavor[1] ? segconf.deep_sam_qname_flavor[1]->name : "N/A";
    else
        return segconf.qname_flavor[q] ? segconf.qname_flavor[q]->name 
             : flag.skip_segconf       ? "<skipped segconf>"
             :                           "N/A";
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

bool qf_is_mated (QType q)
{
    return segconf.qname_flavor[q]->is_mated;
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
