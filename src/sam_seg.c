// ------------------------------------------------------------------
//   sam_seg.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "reference.h"
#include "refhash.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"
#include "zip.h"
#include "dict_id.h"
#include "codec.h"
#include "container.h"
#include "stats.h"
#include "kraken.h"
#include "segconf.h"
#include "contigs.h"
#include "chrom.h"
#include "qname.h"
#include "lookback.h"
#include "gencomp.h"
#include "tip.h"
#include "deep.h"
#include "arch.h"
#include "threads.h"
#include "license.h"

typedef enum { QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, AUX } SamFields __attribute__((unused)); // quick way to define constants

STRl(taxid_redirection_snip, 100);
STRl(copy_NM_snip, 30);
STRl(copy_GX_snip, 30);
STRl(copy_POS_snip, 30);
STRl(copy_mate_CIGAR_snip, 30);
STRl(copy_mate_MAPQ_snip, 30);
STRl(copy_mate_MQ_snip, 30);
STRl(XA_lookback_snip, 30);
STRl(TX_lookback_snip, 30);
STRl(AN_lookback_snip, 30);
sSTRl(copy_mate_RNEXT_snip, 30);
sSTRl(copy_mate_RNAME_snip, 30);
sSTRl(copy_saggy_RNAME_snip, 30);
STRl(copy_mate_YS_snip, 30);
STRl(copy_mate_AS_snip, 30);
STRl(copy_mate_PNEXT_snip, 100);
STRl(copy_saggy_PNEXT_snip, 100);
STRl(copy_mate_POS_snip, 100);
STRl(copy_mate_ms_snip, 30);
STRl(copy_mate_nM_snip, 30);
STRl(copy_buddy_NH_snip, 30);
STRl(redirect_to_CR_X_snip, 30);
STRl(redirect_to_GR_X_snip, 30);
STRl(redirect_to_GY_X_snip, 30);
STRl(redirect_to_RX_X_snip, 30);
char copy_buddy_Z_snips[NUM_MATED_Z_TAGS][30]; unsigned copy_buddy_Z_snip_lens[NUM_MATED_Z_TAGS];

WordIndex xa_lookback_strand_word_index = WORD_INDEX_NONE, xa_lookback_rname_word_index = WORD_INDEX_NONE;

Did buddied_Z_dids[NUM_MATED_Z_TAGS] = MATED_Z_DIDs;

// called by zfile_compress_genozip_header to set FlagsGenozipHeader.dt_specific
bool sam_zip_dts_flag (int dts)
{
    return (dts==1) ? IS_REF_INTERNAL : flag.deep;
}

// ----------------------
// Seg stuff
// ----------------------

// detect if a generic file is actually a SAM
bool is_sam (STRp(header), bool *need_more)
{
    #define H(x) !memcmp (header, "@" #x "\t", 4)

    return header_len >= 4 && (H(HD) || H(CO) || H(SQ) || H(RG) || H(PG)); 
}

void sam_zip_free_end_of_z (void)
{
    sam_header_finalize();
}

// main thread, called for each component. called before segconf.
void sam_zip_initialize (void)
{
    bool has_hdr_contigs = sam_hdr_contigs && sam_hdr_contigs->contigs.len;

    // Copy header contigs to RNAME upon first component. This is in the order of the header, and
    // encodes ref_id based on header order. Note: BAM always has header contigs (might be 0 of them, 
    // for an unaligned file), while SAM is allowed to be header-less (but then requires an external reference)
    if (z_file->num_txts_so_far == 1 && has_hdr_contigs) 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, sam_hdr_contigs);

    seg_prepare_snip_other (SNIP_REDIRECTION, _SAM_TAXID, false, 0, taxid_redirection_snip);

    qname_zip_initialize();

    seg_prepare_snip_other (SNIP_COPY, _SAM_POS, false, 0, copy_POS_snip);
    seg_prepare_snip_other (SNIP_COPY, _OPTION_NM_i, false, 0, copy_NM_snip);
    seg_prepare_snip_other (SNIP_COPY, _OPTION_GX_Z, false, 0, copy_GX_snip);
    seg_prepare_snip_other (SNIP_LOOKBACK, _OPTION_XA_LOOKBACK, false, 0, XA_lookback_snip);
    seg_prepare_snip_other (SNIP_LOOKBACK, _OPTION_TX_LOOKBACK, false, 0, TX_lookback_snip);
    seg_prepare_snip_other (SNIP_LOOKBACK, _OPTION_AN_LOOKBACK, false, 0, AN_lookback_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _OPTION_CR_Z_X, false, 0, redirect_to_CR_X_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _OPTION_GR_Z_X, false, 0, redirect_to_GR_X_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _OPTION_GY_Z_X, false, 0, redirect_to_GY_X_snip);
    seg_prepare_snip_other (SNIP_REDIRECTION, _OPTION_RX_Z_X, false, 0, redirect_to_RX_X_snip);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_PNEXT, copy_mate_PNEXT_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_PNEXT, copy_saggy_PNEXT_snip, '0' + BUDDY_SAGGY);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_POS, copy_mate_POS_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_RNAME, copy_mate_RNAME_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_RNAME, copy_saggy_RNAME_snip, '0' + BUDDY_SAGGY);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_RNEXT, copy_mate_RNEXT_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _SAM_MAPQ, copy_mate_MAPQ_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_MQ_i, copy_mate_MQ_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_ms_i, copy_mate_ms_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_nM_i, copy_mate_nM_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_YS_i, copy_mate_YS_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_AS_i, copy_mate_AS_snip, '0' + BUDDY_MATE);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_NH_i, copy_buddy_NH_snip, '0' + BUDDY_EITHER);
    seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY_CIGAR, _SAM_CIGAR, copy_mate_CIGAR_snip, '0' + BUDDY_MATE);

    // we seg into mux channels, but we copy from the parent
    for (MatedZFields f = 0; f < NUM_MATED_Z_TAGS; f++)
    {
        copy_buddy_Z_snips[f][0] = SNIP_SPECIAL;
        copy_buddy_Z_snip_lens[f] = sizeof (copy_buddy_Z_snips[f]) - 1;
        seg_prepare_snip_other_do (SAM_SPECIAL_COPY_BUDDY, ZCTX(buddied_Z_dids[f])->dict_id,
                                  true, 0, '0' + BUDDY_EITHER, &copy_buddy_Z_snips[f][1], &copy_buddy_Z_snip_lens[f]);
        copy_buddy_Z_snip_lens[f]++;
    }

    if (flag.deep)
        license_show_deep_notice();
}

// called after each txt file (including after the SAM component in Deep)
void sam_zip_finalize (bool is_last_user_txt_file)
{
    threads_log_by_vb (evb, "main_thread", "sam_zip_finalize", 0);

    if (!flag.let_OS_cleanup_on_exit) {

        if ((IS_REF_INTERNAL || IS_REF_EXT_STORE) && !flag.deep)
            ref_destroy_reference (gref);

        if (segconf.sag_type) 
            gencomp_destroy();

        if (flag.deep) 
            sam_deep_zip_finalize();
    }
}

// called by main thread after reading txt of one vb into vb->txt_data
void sam_zip_init_vb (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    vb->chrom_node_index = WORD_INDEX_NONE;
    vb->XG_inc_S = unknown;

    // note: we test for sorted and not collated, because we want non-sorted long read files (which are collated)
    // to seg depn against same-VB prim (i.e. not gencomp) - as the depn lines will follow the prim line
    vb->check_for_gc = !segconf.running && segconf.sag_type && !segconf.abort_gencomp && IS_MAIN(vb);

    sam_sag_zip_init_vb (vb);
}

// called compute thread after compress, order of VBs is arbitrary
void sam_zip_after_compress (VBlockP vb)
{
    // Only the MAIN component produces gencomp lines, however we are absorbing VBs in order, so out-of-band VBs
    // need to be sent too, just to keep the order
    if (segconf.sag_type && (IS_MAIN(vb) || IS_PRIM(vb)))
        gencomp_absorb_vb_gencomp_lines (vb);
}

// called by main thread, as VBs complete (might be out-of-order)
void sam_zip_after_compute (VBlockP vb)
{
    if (IS_MAIN(vb))
        sam_zip_gc_after_compute_main (VB_SAM);

    else if (IS_PRIM(vb))
        gencomp_sam_prim_vb_has_been_ingested (vb);
}

// main thread: writing data-type specific fields to genozip header
void sam_zip_genozip_header (SectionHeaderGenozipHeaderP header)
{
    header->sam.segconf_seq_len         = BGEN32(segconf.sam_seq_len);   // v14
    header->sam.segconf_seq_len_cm      = segconf.seq_len_to_cm;         // v14
    header->sam.segconf_ms_type         = segconf.sam_ms_type;           // v14
    header->sam.segconf_has_MD_or_NM    = segconf.has_MD_or_NM;          // v14
    header->sam.segconf_bisulfite       = segconf.sam_bisulfite;         // v14
    header->sam.segconf_is_paired       = segconf.is_paired;             // v14
    header->sam.segconf_sag_type        = segconf.sag_type;              // v14
    header->sam.segconf_sag_has_AS      = segconf.sag_has_AS;            // v14
    header->sam.segconf_SA_HtoS         = (segconf.SA_HtoS == yes);      // v14
    header->sam.segconf_is_sorted       = segconf.is_sorted;             // v14
    header->sam.segconf_is_collated     = segconf.is_collated;           // v14
    header->sam.segconf_pysam_qual      = segconf.pysam_qual;            // v14
    header->sam.segconf_cellranger      = segconf.has_cellranger;        // v14
    header->sam.segconf_seq_len_dict_id = segconf.seq_len_dict_id;       // v14
    header->sam.segconf_MD_NM_by_un     = segconf.MD_NM_by_unconverted;  // v14
    header->sam.segconf_predict_meth    = segconf.sam_predict_meth_call; // v14
    header->sam.segconf_deep_qname1     = (segconf.deep_qtype == QNAME1);// v15
    header->sam.segconf_deep_qname2     = (segconf.deep_qtype == QNAME2);// v15
    header->sam.segconf_deep_no_qual    = segconf.deep_no_qual;          // v15
}

// initialize SA and OA
static void sam_seg_0X_initialize (VBlockP vb, Did strand_did_i)
{
    // create strand nodes (nodes will be deleted in sam_seg_finalize if not used)
    ctx_create_node (vb, strand_did_i, cSTR("-"));
    ctx_create_node (vb, strand_did_i, cSTR("+"));
    CTX(strand_did_i)->no_vb1_sort = true; // keep them in this ^ order
}

static void sam_seg_QNAME_initialize (VBlockSAMP vb)
{
    CTX(SAM_QNAME)->no_stons = true;             // no singletons, bc sam_piz_special_SET_BUDDY uses PEEK_SNIP
    CTX(SAM_QNAME)->flags.store_per_line = true; // 12.0.41

    qname_seg_initialize (VB, QNAME1, SAM_QNAME); 

    if (segconf.running)
        segconf.qname_flavor[0] = 0; // unknown

    // initial allocations based on segconf data
    else {
        vb->qname_hash.prm8[0] = MIN_(20, MAX_(14, 32 - __builtin_clz (vb->lines.len32 * 5))); // between 14 and 20 bits - tested - no additional compression benefit beyond 20 bits
        buf_alloc_255(vb, &vb->qname_hash, 0, (1ULL << vb->qname_hash.prm8[0]), int32_t, 1, "qname_hash");
    }

    // all-the-same for SAM_BUDDY
    seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SET_BUDDY}, 2, SAM_BUDDY, 0);
}

void sam_seg_initialize (VBlockP vb_)
{
    START_TIMER;
    VBlockSAMP vb = (VBlockSAMP)vb_;

    // case: this is a DEPN VB where data needs to be re-read from txt file - reread data now into txt_data
    if (vb->reread_prescription.len)
        gencomp_reread_lines_as_prescribed (VB);

    // all numeric fields need STORE_INT / STORE_FLOAT to be reconstructable to BAM (possibly already set)
    // via the translators set in the SAM_TOP2BAM Container
    #define T(cond, did_i) ((cond) ? (did_i) : DID_NONE)
    
    ctx_set_store (VB, STORE_INT, SAM_TLEN, SAM_MAPQ, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_GPOS,
                   OPTION_NM_i, OPTION_AS_i, OPTION_MQ_i, OPTION_XS_i, OPTION_XM_i, OPTION_mc_i, OPTION_ms_i, OPTION_Z5_i,
                   OPTION_tx_i, OPTION_YS_i, OPTION_XC_i, OPTION_AM_i, OPTION_SM_i, OPTION_X0_i, OPTION_X1_i, OPTION_CP_i,
                   OPTION_OP_i, OPTION_NH_i, OPTION_HI_i, OPTION_UQ_i, OPTION_cm_i, 
                   OPTION_SA_POS, OPTION_OA_POS,
                   T(TECH(PACBIO), OPTION_qs_i), T(TECH(PACBIO), OPTION_qe_i),
                   T(MP(BLASR), OPTION_XS_i), T(MP(BLASR), OPTION_XE_i), T(MP(BLASR), OPTION_XQ_i), T(MP(BLASR), OPTION_XL_i), T(MP(BLASR), OPTION_FI_i),
                   T(MP(NGMLR), OPTION_QS_i), T(MP(NGMLR), OPTION_QE_i), T(MP(NGMLR), OPTION_XR_i),
                   T(segconf.is_minimap2, OPTION_s1_i), OPTION_ZS_i, OPTION_nM_i,
                   T(kraken_is_loaded, SAM_TAXID),
                   DID_EOL);

    ctx_set_store (VB, STORE_INDEX, SAM_RNAME, SAM_RNEXT, // when reconstructing BAM, we output the word_index instead of the string
                   DID_EOL);

    // when reconstructing these contexts against another context (DELTA_OTHER or XOR_DIFF) the other may be before or after
    ctx_set_same_line (VB, OPTION_AS_i, OPTION_s1_i, // AS may be DELTA_OTHER vs ms:i ; s1 vs AS ; XS vs AS
                       T(segconf.sam_has_BWA_XS_i, OPTION_XS_i),
                       OPTION_RX_Z, OPTION_CR_Z, // whe UR and CR are reconstruted as XOR_DIFF against UB and CB respectively, we search for the values on the same line (befor or after)
                       DID_EOL);

    // don't store singletons in local. note: automatically implied if ltype!=LT_TEXT is set
    ctx_set_no_stons (VB,
                      SAM_RNAME, SAM_RNEXT, // BAM reconstruction needs RNAME, RNEXT word indices. also needed for random access.
                      OPTION_MD_Z,
                      OPTION_BI_Z, OPTION_BD_Z,                                               // we can't use local for singletons in BD or BI as next_local is used by sam_piz_special_BD_BI to point into BD_BI
                      OPTION_SA_CIGAR, OPTION_XA_CIGAR, OPTION_OA_CIGAR, OPTION_OC_Z,         // we can't use local for singletons bc sam_seg_other_CIGAR manually offloads CIGARs to local
                      OPTION_QX_Z,
                      SAM_POS, SAM_PNEXT, OPTION_mc_i, OPTION_OP_i, OPTION_Z5_i, OPTION_CP_i, // required by seg_pos_field
                      T(MP(BSSEEKER2), OPTION_XM_Z), T(MP(BSSEEKER2), OPTION_XG_Z),
                      T(MP(BSBOLT), OPTION_XB_Z),
                      T(MP(BISMARK), OPTION_XM_Z),
                      OPTION_CY_Z, OPTION_QT_Z, OPTION_CB_Z, // barcode fallback segging - add to local text
                      T(kraken_is_loaded, SAM_TAXID),
                      DID_EOL);

    // also implies no_stons
    ctx_set_ltype (VB, LT_UINT8, SAM_MAPQ, OPTION_SA_MAPQ, OPTION_SM_i, OPTION_AM_i, OPTION_OA_MAPQ, // MAPQ (and hence other fields carrying mapping quality) is uint8_t by BAM specification
                   OPTION_RX_Z, OPTION_CR_Z, OPTION_GR_Z, OPTION_GY_Z,                               // local has xor_diff data
                   DID_EOL);

    // initialize these to LT_SEQUENCE, the qual-type ones might be changed later to LT_CODEC (eg domq, longr)
    ctx_set_ltype (VB, LT_SEQUENCE, SAM_QUAL, SAM_QUAL_FLANK, OPTION_BD_BI, OPTION_iq_sq_dq, 
                   OPTION_QX_Z, OPTION_CY_ARR, OPTION_QT_ARR,
                   OPTION_CR_Z_X, OPTION_RX_Z_X, OPTION_2R_Z, OPTION_TR_Z, DID_EOL);

    // set ltype=LT_DYN_INT to allow use of seg_integer
    ctx_set_ltype (VB, LT_DYN_INT, SAM_BUDDY, OPTION_HI_i, OPTION_NM_i, OPTION_NH_i, OPTION_XM_i, OPTION_X1_i,
                   OPTION_AS_i, OPTION_XS_i, OPTION_ZS_i, OPTION_cm_i, OPTION_ms_i, OPTION_nM_i, OPTION_UQ_i,
                   T(segconf.has_TLEN_non_zero, SAM_TLEN), // note: we don't set if !has_TLEN_non_zero, bc values are stored in b250 and may require singletons
                   T(TECH(PACBIO), OPTION_qs_i), T(TECH(PACBIO), OPTION_qe_i),
                   T(MP(BLASR), OPTION_XS_i), T(MP(BLASR), OPTION_XE_i), T(MP(BLASR), OPTION_XQ_i), T(MP(BLASR), OPTION_XL_i), T(MP(BLASR), OPTION_FI_i),
                   T(MP(NGMLR), OPTION_QS_i), T(MP(NGMLR), OPTION_QE_i), T(MP(NGMLR), OPTION_XR_i),
                   DID_EOL);

    ctx_set_ltype (VB, LT_UINT32, OPTION_CP_i, DID_EOL);

    if (segconf.is_collated)
        CTX(SAM_POS)->flags.store_delta = true; // since v12.0.41

    // we may use mates (other than for QNAME) if not is_long_reads (meaning: no mates in this file) and not DEPN components (bc we seg against PRIM)
    if (segconf.is_paired && !IS_DEPN(vb))
        ctx_set_store_per_line (VB, SAM_RNAME, SAM_RNEXT, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_MAPQ, SAM_CIGAR,
                                OPTION_MQ_i, OPTION_MC_Z, OPTION_SM_i, DID_EOL);

    // case: some lines may be segged against a in-VB saggy line
    if (IS_MAIN(vb)) // 14.0.0
        ctx_set_store_per_line (VB, SAM_RNAME, SAM_RNEXT, SAM_PNEXT, SAM_POS, SAM_CIGAR, SAM_MAPQ, SAM_FLAG,
                                OPTION_SA_Z, OPTION_NM_i, DID_EOL);

    ctx_set_store_per_line (VB, OPTION_NH_i, T(segconf.is_paired && segconf.sam_multi_RG, OPTION_RG_Z), DID_EOL);

    if (segconf.is_paired)
        for (MatedZFields f = 1; f < NUM_MATED_Z_TAGS; f++) // note: f starts from 1, bc RG is set above with T()
            ctx_set_store_per_line (VB, buddied_Z_dids[f], DID_EOL);

    if (kraken_is_loaded)
        CTX(SAM_TAXID)->counts_section = true;

    if (TECH(PACBIO)) {
        ctx_set_no_stons (VB, OPTION_dt_Z, OPTION_mq_Z, OPTION_st_Z, 
                          OPTION_dq_Z, OPTION_sq_Z, OPTION_iq_Z, DID_EOL); // Following BD_BI logic: we can't use local for singletons as next_local is used to point into iq_sq_dq

        if (segconf.use_pacbio_iqsqdq) {
            ctx_set_ltype (VB, LT_SEQUENCE, OPTION_dq_Z, OPTION_iq_Z, OPTION_sq_Z, DID_EOL);

            ctx_consolidate_stats (VB, OPTION_iq_sq_dq, OPTION_dq_Z, OPTION_iq_Z, OPTION_sq_Z, DID_EOL);
        }
    }

    // in --stats, consolidate stats
    ctx_consolidate_stats (VB, SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, SAM_SEQMIS_A, SAM_SEQMIS_C, SAM_SEQMIS_G, SAM_SEQMIS_T, DID_EOL);
    ctx_consolidate_stats (VB, SAM_QUAL, SAM_DOMQRUNS, SAM_QUALMPLX, SAM_DIVRQUAL, SAM_QUALSA, SAM_QUAL_PACBIO_DIFF,
                           SAM_QUAL_FLANK, SAM_QUAL_FLANK_DOMQRUNS, SAM_QUAL_FLANK_QUALMPLX, SAM_QUAL_FLANK_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_OQ_Z, OPTION_OQ_DOMQRUNS, OPTION_OQ_QUALMPLX, OPTION_OQ_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_TQ_Z, OPTION_TQ_DOMQRUNS, OPTION_TQ_QUALMPLX, OPTION_TQ_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_2Y_Z, OPTION_2Y_DOMQRUNS, OPTION_2Y_QUALMPLX, OPTION_2Y_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_QX_Z, OPTION_QX_DOMQRUNS, OPTION_QX_QUALMPLX, OPTION_QX_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_CY_Z, OPTION_CY_ARR, OPTION_CY_DOMQRUNS, OPTION_CY_QUALMPLX, OPTION_CY_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_QT_Z, OPTION_QT_ARR, OPTION_QT_DOMQRUNS, OPTION_QT_QUALMPLX, OPTION_QT_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_RX_Z, OPTION_RX_Z_X, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_U2_Z, OPTION_U2_DOMQRUNS, OPTION_U2_QUALMPLX, OPTION_U2_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_E2_Z, OPTION_2NONREF, OPTION_N2ONREFX, OPTION_2GPOS, OPTION_S2TRAND, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_OA_Z, OPTION_OA_RNAME, OPTION_OA_POS, OPTION_OA_STRAND, OPTION_OA_CIGAR, OPTION_OA_MAPQ, OPTION_OA_NM, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_XA_Z, OPTION_XA_RNAME, OPTION_XA_POS, OPTION_XA_STRAND, OPTION_XA_CIGAR, OPTION_XA_NM, OPTION_XA_STRAND_POS, OPTION_XA_LOOKBACK, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_CB_Z, OPTION_CB_ARR, OPTION_CB_SUFFIX, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_CR_Z, OPTION_CR_Z_X, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_GR_Z, OPTION_GR_Z_X, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_GY_Z, OPTION_GY_Z_X, DID_EOL);
    ctx_consolidate_stats (VB, SAM_QNAME, SAM_BUDDY, SAM_QNAMESA, SAM_FQ_AUX, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_BD_BI, OPTION_BI_Z, OPTION_BD_Z, DID_EOL);

    if (segconf.has[OPTION_HI_i] && !segconf.has[OPTION_SA_Z])
        ctx_consolidate_stats (VB, OPTION_HI_i, OPTION_SA_Z, DID_EOL);
    else
        ctx_consolidate_stats (VB, OPTION_SA_Z, OPTION_SA_RNAME, OPTION_SA_POS, OPTION_SA_STRAND, OPTION_SA_CIGAR, OPTION_SA_MAPQ, OPTION_SA_NM, OPTION_SA_MAIN, DID_EOL);

    codec_acgt_comp_init (VB, SAM_NONREF);

    sam_seg_QNAME_initialize (vb);
    sam_seg_QUAL_initialize (vb);
    sam_seg_SEQ_initialize (vb);
    sam_seg_cigar_initialize (vb);
    sam_seg_gc_initialize (vb);
    sam_seg_0X_initialize (VB, OPTION_SA_STRAND);
    sam_seg_0X_initialize (VB, OPTION_OA_STRAND);

    if (segconf.has_cellranger) {
        sam_seg_TX_AN_initialize (vb, OPTION_TX_Z);
        sam_seg_TX_AN_initialize (vb, OPTION_AN_Z);
    }

    ctx_set_store (VB, STORE_INDEX, OPTION_XA_Z, DID_EOL); // for containers this stores repeats - used by sam_piz_special_X1->container_peek_repeats

    if (segconf.sam_has_BWA_XS_i) // XS:i is as defined some aligners
        seg_mux_init (VB, CTX(OPTION_XS_i), 4, SAM_SPECIAL_BWA_XS, false, (MultiplexerP)&vb->mux_XS, "0123");

    else if (MP(HISAT2)) // ZS:i is like BWA's XS:i
        seg_mux_init (VB, CTX(OPTION_ZS_i), 4, SAM_SPECIAL_BWA_XS, false, (MultiplexerP)&vb->mux_XS, "0123");

    if (sam_has_bowtie2_YS_i()) {
        ctx_set_store_per_line (VB, OPTION_AS_i, OPTION_YS_i, DID_EOL);
        seg_mux_init (VB, CTX(OPTION_YS_i), 2, SAM_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mux_YS, "01");
    }

    if (MP(STAR) && segconf.is_paired && !IS_DEPN(vb)) {
        seg_mux_init (VB, CTX(OPTION_nM_i), 2, SAM_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mux_nM, "01");
        ctx_set_store_per_line (VB, OPTION_AS_i, OPTION_nM_i, DID_EOL);
    }

    if (segconf.sam_ms_type == ms_BIOBAMBAM) {
        seg_mux_init (VB, CTX(OPTION_ms_i), 2, SAM_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mux_ms, "01");

        CTX(OPTION_ms_i)->flags.spl_custom = true;  // custom store-per-line - SPECIAL will handle the storing
        CTX(OPTION_ms_i)->flags.store = STORE_INT;  // since v14 - store QUAL_score for mate ms:i (in v13 it was stored in QUAL)
    }
    
    seg_mux_init (VB, CTX(SAM_FLAG), 2, SAM_SPECIAL_DEMUX_BY_BUDDY, false, (MultiplexerP)&vb->mux_FLAG, "01");
    seg_mux_init (VB, CTX(SAM_POS), 3, SAM_SPECIAL_DEMUX_BY_MATE_PRIM, true, (MultiplexerP)&vb->mux_POS, "012");
    seg_mux_init (VB, CTX(SAM_PNEXT), 4, SAM_SPECIAL_PNEXT, true, (MultiplexerP)&vb->mux_PNEXT, "0123");
    seg_mux_init (VB, CTX(SAM_MAPQ), 3, SAM_SPECIAL_DEMUX_BY_MATE_PRIM, false, (MultiplexerP)&vb->mux_MAPQ, "012");
    seg_mux_init (VB, CTX(OPTION_MQ_i), 2, SAM_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mux_MQ, "01");
    seg_mux_init (VB, CTX(OPTION_MC_Z), 2, SAM_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mux_MC, "01");
    seg_mux_init (VB, CTX(OPTION_AS_i), 2, SAM_SPECIAL_DEMUX_BY_MATE, false, (MultiplexerP)&vb->mux_AS, "01");
    seg_mux_init (VB, CTX(OPTION_NH_i), 3, SAM_SPECIAL_DEMUX_BY_BUDDY_MAP, false, (MultiplexerP)&vb->mux_NH, "012");

    for (MatedZFields f = 0; f < NUM_MATED_Z_TAGS; f++)
        seg_mux_init (VB, CTX(buddied_Z_dids[f]), 2, SAM_SPECIAL_DEMUX_BY_BUDDY, false, (MultiplexerP)&vb->mux_mated_z_fields[f], "01");
        
    // get counts of qnames in the VB, that will allow us to leave lines in th  e VB for saggy instead of gencomp
    // note: for SAG_BY_SA in --best, we send all prim/depn to gencomp (see bug 629)
    if (vb->check_for_gc && !flag.force_gencomp && (/*IS_SAG_SA ||*/ IS_SAG_NH || IS_SAG_SOLO)) // SA still doesn't work well with saggy - not even QUAL
        scan_index_qnames_seg (vb);

    // note: SA_HtoS might be already set if we segged an line containing SA:Z against a prim saggy
    if (IS_DEPN(vb) && IS_SAG_SA && segconf.SA_HtoS == unknown)
        segconf.SA_HtoS = segconf.depn_CIGAR_can_have_H && // set when segging MAIN component
                          !segconf.SA_CIGAR_can_have_H;    // set when segging PRIM component

    COPY_TIMER(seg_initialize);
#undef T
}

static void sam_seg_toplevel (VBlockP vb)
{
    // in PRIM, SA group loader preprocessor reconstructs QNAME / SQBITMAP / QUAL to load SA Groups, while recon
    // reconstructs QNAMESA, SEQSA, QUALSA to copy from SA Groups. In contrast, in DEPN, copy from SA Groups
    // occurs (if needed) when reconstructing QNAME / SQBITMAP / QUAL.
    uint64_t qname_dict_id = (IS_PRIM(vb) ? _SAM_QNAMESA : _SAM_QNAME);
    uint64_t qual_dict_id  = (IS_PRIM(vb) ? _SAM_QUALSA  : _SAM_QUAL );

    uint8_t deep_cb = flag.deep && (vb->comp_i != SAM_COMP_DEPN) ? CI1_ITEM_CB : 0; // note: DEPN alignments don't participate in Deep

    // top level snip - reconstruction as SAM
    SmallContainer top_level_sam = {
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 14,
        .items = { { .dict_id = { _SAM_BUDDY    }                                 },
                   { .dict_id = { qname_dict_id }, .separator = "\t"              },
                   { .dict_id = { _SAM_FLAG     }, .separator = "\t"              },
                   { .dict_id = { _SAM_RNAME    }, .separator = "\t"              },
                   { .dict_id = { _SAM_POS      }, .separator = "\t"              },
                   { .dict_id = { _SAM_MAPQ     }, .separator = "\t"              },
                   { .dict_id = { _SAM_CIGAR    }, .separator = "\t"              },
                   { .dict_id = { _SAM_RNEXT    }, .separator = "\t"              },
                   { .dict_id = { _SAM_PNEXT    }, .separator = "\t"              },
                   { .dict_id = { _SAM_TLEN     }, .separator = "\t"              },
                   { .dict_id = { _SAM_SQBITMAP }, .separator = { '\t', deep_cb } },
                   { .dict_id = { qual_dict_id  }, .separator = { '\t', deep_cb } },
                   { .dict_id = { _SAM_AUX      }                                 },
                   { .dict_id = { _SAM_EOL      }                                 }}};
                   
    container_seg (vb, CTX(SAM_TOPLEVEL), (ContainerP)&top_level_sam, 0, 0, 0);

    // Deep: top level snip - reconstruction as SAM when genocat --R1 --R2 --R --fq - reconstruct only what's needed for FASTQ 
    SmallContainer top_level_deep_fq = {
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 11,
        .items = { { .dict_id = { _SAM_BUDDY    }                                 },
                   { .dict_id = { qname_dict_id }, .separator = "\t"              },
                   { .dict_id = { _SAM_FLAG     }, .separator = { CI0_INVISIBLE } },
                   { .dict_id = { _SAM_RNAME    }, .separator = { CI0_INVISIBLE } }, // needed for reconstructing seq
                   { .dict_id = { _SAM_POS      }, .separator = { CI0_INVISIBLE } }, // needed for reconstructing seq
                   { .dict_id = { _SAM_CIGAR    }, .separator = { CI0_INVISIBLE } }, // needed for reconstructing seq
                   { .dict_id = { _SAM_RNEXT    }, .separator = { CI0_INVISIBLE } }, // needed for reconstructing RNAME
                   { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_INVISIBLE } }, // needed for reconstructing POS
                   { .dict_id = { _SAM_SQBITMAP }, .separator = { '\t', deep_cb } },
                   { .dict_id = { qual_dict_id  }, .separator = { '\t', deep_cb } },
                   { .dict_id = { _SAM_FQ_AUX   }, .separator = { CI0_INVISIBLE } }}}; // consumes OPTION_MC_Z (needed for mate CIGAR), OPTION_SA_Z (need for saggy RNAME, POS, CIGAR, NM, MAPQ)
                   
    container_seg (vb, CTX(SAM_TOP2NONE), (ContainerP)&top_level_deep_fq, 0, 0, 0);

    // top level snip - reconstruction as BAM
    // strategy: we start by reconstructing the variable-length fields first (after a prefix that sets them in place)
    // - read_name, cigar, seq and qual - and then go back and fill in the fixed-location fields
    // Translation (a feature of Container): items reconstruct their data and then call a translation function to translate it to the desired format
    SmallContainer top_level_bam = {
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 15,
        .items        = { { .dict_id = { _SAM_BUDDY    }                                                                       },
                          { .dict_id = { _SAM_RNAME    }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                          { .dict_id = { _SAM_POS      }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 1 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                          { .dict_id = { _SAM_MAPQ     }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_U8       }, // Translate - textual to binary number
                          { .dict_id = { _SAM_BAM_BIN  }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 2 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                          { .dict_id = { _SAM_FLAG     }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 4 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                          { .dict_id = { _SAM_RNEXT    }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                          { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 4 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                          { .dict_id = { qname_dict_id }, .separator = { CI0_TRANS_NUL                     },                  }, 
                          { .dict_id = { _SAM_CIGAR    }, .separator = ""                                                      }, // handle in special reconstructor - translate textual to BAM CIGAR format + reconstruct l_read_name, n_cigar_op, l_seq
                          { .dict_id = { _SAM_TLEN     }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_TLEN     }, // must be after CIGAR bc sam_piz_special_TLEN_old needs vb->seq_num
                          { .dict_id = { _SAM_SQBITMAP }, .separator = "",                                    SAM2BAM_SEQ      }, // Translate - textual format to BAM format ; if --deep, store in deep_ents
                          { .dict_id = { qual_dict_id  }, .separator = "",                                    SAM2BAM_QUAL     }, // Translate - textual format to BAM format ; set block_size ; if --deep, store in deep_ents
                          { .dict_id = { _SAM_AUX      }, .separator = { CI0_TRANS_NOR                     }                   }, // up to v11, this had the SAM2BAM_AUX translator
                        }
    };

    // no container wide-prefix, skip l_name with a 4-character prefix
    static const char bam_line_prefix[] = {CON_PX_SEP,                      // has prefix
                                           CON_PX_SEP,                      // end of (empty) container-wide prefix
                                           ' ', ' ', ' ', ' ', CON_PX_SEP}; // first item prefix - 4 spaces (place holder for block_size)

    container_seg (vb, CTX(SAM_TOP2BAM), (ContainerP)&top_level_bam, bam_line_prefix, sizeof (bam_line_prefix),
                  IS_BAM_ZIP ? sizeof (uint32_t) * vb->lines.len : 0); // if BAM, account for block_size

    // top level snip - reconstruction as FASTQ
    SmallContainer top_level_fastq = {
        .repeats     = vb->lines.len32,
        .is_toplevel = true,
        .callback    = true, // drop supplementary alignments and alignments without QUAL data
        .nitems_lo   = 11,
        .items       = { { .dict_id = { _SAM_BUDDY    }                                                  },
                         { .dict_id = { qname_dict_id }, .separator = "\n"                               },
                         { .dict_id = { _SAM_RNAME    }, .separator = { CI0_INVISIBLE }                  },// needed for reconstructing seq
                         { .dict_id = { _SAM_POS      }, .separator = { CI0_INVISIBLE }                  },// needed for reconstructing seq
                         { .dict_id = { _SAM_RNEXT    }, .separator = { CI0_INVISIBLE }                  },// needed for reconstructing PNEXT
                         { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_INVISIBLE }                  },// needed for reconstructing POS (in case of BUDDY)
                         { .dict_id = { _SAM_FLAG     }, .separator = { CI0_INVISIBLE }, .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_CIGAR    }, .separator = { CI0_INVISIBLE }                  },// needed for reconstructing seq and also of SA_CIGAR (need vb->hard_clips[] for de-squanking)
                         { .dict_id = { _SAM_FQ_AUX   }, .separator = { CI0_INVISIBLE }                  },// consumes OPTION_MC_Z (needed for mate CIGAR), OPTION_SA_Z (need for saggy RNAME, POS, CIGAR, NM, MAPQ)
                         { .dict_id = { _SAM_SQBITMAP }, .separator = { CI0_TRANS_ALWAYS, '\n' }, .translator = SAM2FASTQ_SEQ  },
                         { .dict_id = { qual_dict_id  }, .separator = { CI0_TRANS_ALWAYS, '\n' }, .translator = SAM2FASTQ_QUAL }, // also moves fastq "line" to R2 (paired file) if needed
                      }
    };

    // add a '@' to the description line, use a prefix to add the + line
    static const char fastq_line_prefix[] = {CON_PX_SEP, CON_PX_SEP, '@', CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP,
                                             CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, '+', '\n', CON_PX_SEP};

    container_seg (vb, CTX(SAM_TOP2FQ), (ContainerP)&top_level_fastq, fastq_line_prefix, sizeof (fastq_line_prefix), 0);

    // top level snip - reconstruction as FASTQ "extended" - with all the SAM fields in the description line
    SmallContainer top_level_fastq_ext = {
        .repeats     = vb->lines.len32,
        .is_toplevel = true,
        .callback    = true, // drop non-primary chimeric reads and reads without QUAL data
        .nitems_lo   = 13,
        .items       = { { .dict_id = { _SAM_BUDDY    }                                                  },
                         { .dict_id = { qname_dict_id }, .separator = "\t"                               },
                         { .dict_id = { _SAM_FLAG     }, .separator = "\t", .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_RNAME    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_POS      }, .separator = "\t"                               },
                         { .dict_id = { _SAM_MAPQ     }, .separator = "\t"                               },
                         { .dict_id = { _SAM_CIGAR    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_RNEXT    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_PNEXT    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_TLEN     }, .separator = "\t"                               },
                         { .dict_id = { _SAM_AUX      }, .separator = "\n"                               },
                         { .dict_id = { _SAM_SQBITMAP }, .separator = { CI0_TRANS_ALWAYS, '\n' }, .translator = SAM2FASTQ_SEQ  },
                         { .dict_id = { qual_dict_id  }, .separator = { CI0_TRANS_ALWAYS, '\n' }, .translator = SAM2FASTQ_QUAL }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line
    static const char fastq_ext_line_prefix[] = {
        CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, 
        '@', CON_PX_SEP, 
        'F', 'L', 'A', 'G', ':', CON_PX_SEP,
        'R', 'N', 'A', 'M', 'E', ':', CON_PX_SEP,
        'P', 'O', 'S', ':', CON_PX_SEP,
        'M', 'A', 'P', 'Q', ':', CON_PX_SEP,
        'C', 'I', 'G', 'A', 'R', ':', CON_PX_SEP,
        'R', 'N', 'E', 'X', 'T', ':', CON_PX_SEP,
        'P', 'N', 'E', 'X', 'T', ':', CON_PX_SEP,
        'T', 'L', 'E', 'N', ':', CON_PX_SEP,
        CON_PX_SEP,             // AUX
        CON_PX_SEP,             // SEQ
        '+', '\n', CON_PX_SEP}; // QUAL

    container_seg (vb, CTX(SAM_TOP2FQEX), (ContainerP)&top_level_fastq_ext,
                  fastq_ext_line_prefix, sizeof (fastq_ext_line_prefix), 0);
}

static uint32_t num_lines_at_max_len (VBlockSAMP vb)
{
    uint32_t count = 0;
    for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++)
        if (DATA_LINE(line_i)->SEQ.len == vb->longest_seq_len)
            count++;

    return count;
}

// finalize Seg configuration parameters
static void sam_seg_finalize_segconf (VBlockSAMP vb)
{
    segconf.longest_seq_len = vb->longest_seq_len;
    segconf.is_long_reads   = segconf_is_long_reads();
    segconf.sam_multi_RG    = CTX(OPTION_RG_Z)->nodes.len32 >= 2;
    segconf.sam_cigar_len   = 1 + ((segconf.sam_cigar_len - 1) / vb->lines.len32);                   // set to the average CIGAR len (rounded up)

    if (num_lines_at_max_len(vb) > vb->lines.len32 / 2 &&  // more than half the lines are at exactly maximal length
        (vb->lines.len32 > 100 || txt_file->is_eof)    &&  // enough lines to be reasonably convinced that this is not by chance
        !segconf.is_long_reads)        // TO DO: trimming long-read qual in FASTQ with --deep would mess up LONGR codec, we need to sort this out
        
        segconf.sam_cropped_at = vb->longest_seq_len; // possibily the FASTQ reads were cropped to be all equal length

    segconf.sam_seq_len = (uint32_t)(0.5 + (double)segconf.sam_seq_len / (double)vb->lines.len32); // set average seq_len - rounded to the nearest

    segconf.seq_len_to_cm /= vb->lines.len32;
    if (segconf.seq_len_to_cm > 255)
        segconf.seq_len_to_cm = 0;

    // in some STAR transcriptome files, there are no @PG header lines. it is important that mapper is set to STAR so we can use liberal saggy_line_i, as SA groups are large
    if (MP(UNKNOWN) && segconf.has[OPTION_NH_i] && segconf.has[OPTION_HI_i] && !segconf.has[OPTION_SA_Z])
        segconf.sam_mapper = MP_STAR;

    if (MP(STAR)) segconf.sag_has_AS = true;

    // evidence of STARsolo or cellranger (this is also detected in the SAM header)
    if ((MP(STAR) || MP(UNKNOWN)) &&
        (((segconf.has[OPTION_UB_Z] > 0) + (segconf.has[OPTION_UR_Z] > 0) + (segconf.has[OPTION_UY_Z] > 0) >= 2) ||
         (segconf.has[OPTION_gn_Z] && segconf.has[OPTION_gx_Z]) ||
         (segconf.has[OPTION_GN_Z] && segconf.has[OPTION_GX_Z]) ||
         (segconf.has[OPTION_2R_Z] && segconf.has[OPTION_2Y_Z]) ||
         (segconf.has[OPTION_GR_Z] && segconf.has[OPTION_GY_Z]) ||
         segconf.has[OPTION_TX_Z]                               ||
         segconf.has[OPTION_AN_Z])) {

        segconf.sam_mapper = MP_STAR;
        segconf.star_solo  = true;     // must be before sam_set_sag_type
        segconf.has_cellranger = true; // STAR Solo mimics cellranger's tags
    }

    // we set AS_is_2ref_consumed, meaning AS tends to be near 2 X ref_consumed, if at least half of the lines say so
    segconf.AS_is_2ref_consumed = (segconf.AS_is_2ref_consumed > vb->lines.len32 / 2);

    // possibly reset sam_bisulfite, first set in sam_header_zip_inspect_PG_lines 
    segconf.sam_bisulfite = MP(BISMARK) || MP(BSSEEKER2) || (MP(DRAGEN) && segconf.has[OPTION_XM_Z]) ||
                            MP(BSBOLT) || (MP(GEM3) && segconf.has[OPTION_XB_A]);

    ASSINP (!segconf.sam_bisulfite        // not a bisulfite file
         || flag.force                    // --force overrides
         || IS_REF_EXTERNAL || IS_REF_EXT_STORE   // reference is provided
         || flag.zip_no_z_file,           // we're not creating a compressed format
            "Compressing bisulfite file %s requires using --reference. Override with --force.", dt_name (vb->data_type));

    segconf.sam_predict_meth_call = segconf.sam_bisulfite          &&
                                    !IS_REF_INTERNAL && // bug 648
                                    (MP(BISMARK) || MP(DRAGEN) || MP(BSBOLT)); // have methylation call tags

    // in bisulfate data, we still calculate MD:Z and NM:i vs unconverted reference
    segconf.MD_NM_by_unconverted = MP(BISMARK) || MP(DRAGEN) || MP(BSBOLT) || MP(GEM3) || MP(BSSEEKER2);

    // if we have @HD-SO "coordinate" or "queryname", then we take that as definitive. Otherwise, we go by our segconf sampling.
    if (sam_hd_so == HD_SO_COORDINATE) {
        segconf.is_sorted = true;
        segconf.is_collated = false;
    }

    else if (sam_hd_so == HD_SO_QUERYNAME || // file is sorted alphabetically by QNAME
             (sam_hd_so == HD_SO_UNSORTED && sam_hd_go == HD_GO_QUERY)) { // same QNAMEs are grouped together, but there is no order between groups
        segconf.is_collated = true;
        segconf.is_sorted = false;
    }

    else {
        // case: if we haven't found any pair of consecutive lines with the same CHROM and non-descreasing POS, this is not a sorted file, despite no evidence of "not sorted". eg could be unique CHROMs.
        if (!segconf.evidence_of_sorted)
            segconf.is_sorted = false;
        
        // we have at leaest one pair of lines with the same QNAME, and the file is not sorted
        if (!segconf.is_sorted && segconf.evidence_of_collated)
            segconf.is_collated = true;
    }

    segconf.has_barcodes = segconf.has[OPTION_CB_Z] || segconf.has[OPTION_CR_Z] || segconf.has[OPTION_CY_Z] || segconf.has[OPTION_BX_Z] || 
                           segconf.has[OPTION_RX_Z] || segconf.has[OPTION_QX_Z] || segconf.has[OPTION_BC_Z] || segconf.has[OPTION_QT_Z];

    // second is_paired test: PNEXT is not 0 for all lines (first test is in sam_seg_FLAG)
    if (segconf.has[SAM_PNEXT])
        segconf.is_paired = true;

    segconf.has_MD_or_NM = segconf.has[OPTION_MD_Z] || segconf.has[OPTION_NM_i] ||
                           (MP(STAR) && segconf.has[OPTION_nM_i] && !segconf.is_paired); // we use the NM method to seg nM in this case

    // note: in case of both is_biobambam2_sort and minimap2, it is likely biobambam's tag
    if (segconf.is_biobambam2_sort && segconf.is_paired) 
        segconf.sam_ms_type = ms_BIOBAMBAM;

    else if (segconf.is_minimap2)
        segconf.sam_ms_type = ms_MINIMAP2; 

    // if we have no conclusive evidence of either minimap2 or biobambam, this is likely samtools, which has the same logic as biobambam
    else if (segconf.is_sorted)
        segconf.sam_ms_type = ms_BIOBAMBAM;

    // update context tag names if this file has UB/UR/UY which are aliased to BX/RX/QX
    if (segconf.has[OPTION_UB_Z] || segconf.has[OPTION_UR_Z]) {
        strcpy (ZCTX(OPTION_BX_Z)->tag_name, "UB:Z");
        strcpy (ZCTX(OPTION_RX_Z)->tag_name, "UR:Z");
        strcpy (ZCTX(OPTION_QX_Z)->tag_name, "UY:Z");
    }

    // allow aligner if unmapped file (usually only enabled in best) if we have unmapped reads in segconf indicating a file enriched
    // in unmapped reads (normally, in sorted BAMs unmapped reads are at the end of the file)
    if (segconf.num_mapped < vb->lines.len32 && !flag.aligner_available && IS_REF_LOADED_ZIP) {
        flag.aligner_available = true;
        refhash_load_standalone();
    }

    // with REF_EXTERNAL and unaligned data, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (z_file->num_txts_so_far == 1 && (flag.aligner_available || !sam_hdr_contigs) && IS_REF_LOADED_ZIP)
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, ref_get_ctgs (gref));

    if (TECH(PACBIO) && segconf.has[OPTION_dq_Z] == vb->lines.len32 && segconf.has[OPTION_iq_Z] == vb->lines.len32 && segconf.has[OPTION_sq_Z] == vb->lines.len32)
        segconf.use_pacbio_iqsqdq = true;

    // Aplogize to user for not doing a good job with PacBio subreads files
    if ((MP(BAZ2BAM) || TECH(PACBIO)) && 
            (segconf.has[OPTION_ip_B_C] || segconf.has[OPTION_pw_B_C] || segconf.has[OPTION_fi_B_C] || 
             segconf.has[OPTION_fp_B_C] || segconf.has[OPTION_ri_B_C] || segconf.has[OPTION_rp_B_C])) {
        TEMP_FLAG(quiet, flag.explicit_quiet); // note: quiet is set to true in segconf
        WARN0 ("FYI: Genozip currently doesn't do a very good job at compressing PacBio kinetic BAM files. This is because we haven't figured out yet a good method to compress kinetic data - the ip:B, pw:B, fi:B, fp:B, ri:B, rp:B fields. Sorry!\n");
        RESTORE_FLAG(quiet);
    }

    if (flag.reference == REF_INTERNAL && !txt_file->redirected && (!segconf.sam_is_unmapped || !segconf.is_long_reads))
        TIP ("Compressing a %s file using a reference file can reduce the size by 7%%-30%%%s.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             dt_name (txt_file->data_type), arch_get_argv0(), MP(UNKNOWN) ? " (even for unaligned files)" : "", txt_file->name);
    
    ASSERT (!flag.debug_or_test || !segconf.has[OPTION_SA_Z] || segconf.sam_has_SA_Z || MP(LONGRANGER) || MP(UNKNOWN),
            "%s produces SA:Z, expecting segconf.sam_has_SA_Z to be set", segconf_sam_mapper_name()); // should be set in sam_header_zip_inspect_PG_lines
}

void sam_seg_finalize (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (vb->lines.len) {
        CTX(SAM_SQBITMAP)->local_always = true; // We always include the SQBITMAP local section, except if no lines
        bits_truncate ((BitsP)&CTX(SAM_SQBITMAP)->local, CTX(SAM_SQBITMAP)->next_local); // remove unused bits due to MAPPING_PERFECT
    }

    // assign the QUAL codec
    if (vb->has_qual && CTX(SAM_QUAL)->local.len32) 
        codec_assign_best_qual_codec (VB, SAM_QUAL, sam_zip_qual, false, true);

    if (CTX(OPTION_OQ_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_OQ_Z, sam_zip_OQ, false, true);

    if (CTX(OPTION_TQ_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_TQ_Z, sam_zip_TQ, true, false);

    if (CTX(OPTION_CY_ARR)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_CY_ARR, NULL, true, false);

    if (CTX(OPTION_QX_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_QX_Z, sam_zip_QX, true, false);

    if (CTX(OPTION_2Y_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_2Y_Z, sam_zip_2Y, true, true);

    if (CTX(OPTION_QT_ARR)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_QT_Z, NULL, true, false);

    if (CTX(OPTION_U2_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_U2_Z, sam_zip_U2, true, true);

    // determine if sam_piz_sam2bam_SEQ ought to store vb->textual_seq
    CTX(SAM_SQBITMAP)->flags.no_textual_seq = CTX(SAM_QUAL)->lcodec != CODEC_LONGR && segconf.sam_mapper != MP_BSSEEKER2;

    if (flag.biopsy_line.line_i == NO_LINE) // no --biopsy-line
        sam_seg_toplevel (VB);

    // finalize Seg configuration parameters
    if (segconf.running)
        sam_seg_finalize_segconf (vb);

    // get rid of some contexts if we ended up not have using them
    static const Did delete_me[] = {
        OPTION_XA_LOOKBACK, OPTION_XA_STRAND, OPTION_OA_STRAND, OPTION_SA_STRAND, // nodes added in sam_seg_0X_initialize. note: we use SA_STRAND in PRIM lines even if there is no SA:Z field
        OPTION_TX_NEGATIVE, OPTION_TX_LOOKBACK, OPTION_TX_SAM_POS,
        OPTION_AN_NEGATIVE, OPTION_AN_LOOKBACK, OPTION_AN_SAM_POS };  // nodes added in sam_seg_TX_AN_initialize

    for (int i = 0; i < ARRAY_LEN(delete_me); i++)
        if (!CTX(delete_me[i])->b250.len && !CTX(delete_me[i])->local.len)
            buflist_free_ctx (VB, CTX(delete_me[i]));

    // Primary VB - ingest VB data into z_file->sa_*
    if (IS_PRIM(vb))
        sam_zip_prim_ingest_vb (vb);

    // collect stats
    if (!segconf.running) {
        __atomic_add_fetch(&z_file->mate_line_count,  (uint64_t)vb->mate_line_count,  __ATOMIC_RELAXED);
        __atomic_add_fetch(&z_file->saggy_near_count, (uint64_t)vb->saggy_near_count, __ATOMIC_RELAXED);
        __atomic_add_fetch(&z_file->prim_far_count,   (uint64_t)vb->prim_far_count,   __ATOMIC_RELAXED);
    }

    if (!CTX(SAM_QUAL)->local.len)
        CTX(SAM_QUAL)->ltype = LT_TEXT; // if no QUAL in this VB, switch back to LT_TEXT to enable all_the_same

    // VB1: if we've not found depn lines in the VB, abort gencomp (likely the depn lines were filtered out)
    if (vb->vblock_i == 1 && !vb->seg_found_depn_line && !flag.force_gencomp)
        segconf.abort_gencomp = true;
}

// main thread: called after all VBs, before compressing global sections
void sam_zip_after_vbs (void)
{
    // case: header-only file, completely get rid of RNAME
    if (!sections_count_sections (SEC_VB_HEADER)) 
        ZCTX(SAM_RNAME)->dict.len = ZCTX(SAM_RNAME)->nodes.len = 0;
    
    // shorten unused words in dictionary strings to "" (dict pre-populated in sam_zip_initialize)
    if (!IS_REF_INTERNAL) // TO DO: this doesn't work for REF_INTERNAL for example with test.transcriptome.bam
        ctx_shorten_unused_dict_words (SAM_RNAME);

    ctx_shorten_unused_dict_words (OPTION_XA_STRAND); // remove gem3 bi-sulfite words if not used
}

bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return
        // typically small
        dict_id.num == _SAM_TOPLEVEL     ||
        dict_id.num == _SAM_TOP2BAM      ||
        dict_id.num == _SAM_TOP2FQ       ||
        dict_id.num == _SAM_TOP2FQEX     ||
        dict_id.num == _SAM_FLAG         ||
        dict_id.num == _SAM_MAPQ         ||
        dict_id.num == _SAM_QNAME        ||
        dict_id.num == _SAM_RNAME        ||
        dict_id.num == _SAM_RNEXT        ||
        dict_id.num == _SAM_TLEN         ||
        dict_id.num == _SAM_AUX          ||
        dict_id.num == _SAM_EOL          ||
        dict_id.num == _SAM_TAXID        ||
        dict_id.num == _SAM_BUDDY        ||

        // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
        dict_id.num == _OPTION_AM_i      ||
        dict_id.num == _OPTION_AS_i      ||
        dict_id.num == _OPTION_CC_Z      ||
        dict_id.num == _OPTION_CM_i      ||
        dict_id.num == _OPTION_LB_Z      ||
        dict_id.num == _OPTION_FI_i      ||
        dict_id.num == _OPTION_H0_i      ||
        dict_id.num == _OPTION_H1_i      ||
        dict_id.num == _OPTION_H2_i      ||
        dict_id.num == _OPTION_HI_i      ||
        dict_id.num == _OPTION_MQ_i      ||
        dict_id.num == _OPTION_NH_i      ||
        dict_id.num == _OPTION_NM_i      ||
        dict_id.num == _OPTION_OA_Z      ||
        dict_id.num == _OPTION_OA_RNAME  ||
        dict_id.num == _OPTION_OA_NM     ||
        dict_id.num == _OPTION_OA_STRAND ||
        dict_id.num == _OPTION_OA_MAPQ   ||
        dict_id.num == _OPTION_OA_CIGAR  ||
        dict_id.num == _OPTION_OC_Z      ||
        dict_id.num == _OPTION_PG_Z      ||
        dict_id.num == _OPTION_PQ_i      ||
        dict_id.num == _OPTION_PU_Z      ||
        dict_id.num == _OPTION_RG_Z      ||
        dict_id.num == _OPTION_RG_Z      ||
        dict_id.num == _OPTION_SA_Z      ||
        dict_id.num == _OPTION_SA_RNAME  ||
        dict_id.num == _OPTION_SA_NM     ||
        dict_id.num == _OPTION_SA_STRAND ||
        dict_id.num == _OPTION_SA_MAPQ   ||
        dict_id.num == _OPTION_SA_CIGAR  ||
        dict_id.num == _OPTION_SM_i      ||
        dict_id.num == _OPTION_TC_i      ||
        dict_id.num == _OPTION_UQ_i      ||

        // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
        dict_id.num == _OPTION_X0_i      ||
        dict_id.num == _OPTION_X1_i      ||
        dict_id.num == _OPTION_XA_Z      ||
        dict_id.num == _OPTION_XA_RNAME  ||
        dict_id.num == _OPTION_XA_NM     ||
        dict_id.num == _OPTION_XA_STRAND ||
        dict_id.num == _OPTION_XA_CIGAR  ||
        dict_id.num == _OPTION_XN_i      ||
        dict_id.num == _OPTION_XM_i      ||
        dict_id.num == _OPTION_XO_i      ||
        dict_id.num == _OPTION_XG_i      ||
        dict_id.num == _OPTION_XS_i      ||
        dict_id.num == _OPTION_XE_i      ||
        0;
}

bool sam_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id)
{
    return
        // typically big
        dict_id.num == _OPTION_BX_Z      ||
        dict_id.num == _OPTION_TX_GENE   ||
        dict_id.num == _OPTION_AN_GENE;
}

static void sam_get_one_aux (VBlockSAMP vb, int16_t idx,
                             rom *tag, char *type, char *array_subtype, pSTRp (value)) // out
{
    rom aux = vb->auxs[idx];
    uint32_t aux_len = vb->aux_lens[idx];

    ASSSEG (aux_len >= 6 && aux[2] == ':' && aux[4] == ':', "invalid optional field format: %.*s", STRf (aux));

    *tag = aux;
    *type = aux[3];

    if (*type == 'B') {
        *array_subtype = aux[5];
        *value = aux + 7;
        *value_len = aux_len - 7;
    }

    else {
        *array_subtype = 0;
        *value = aux + 5;
        *value_len = aux_len - 5;
    }

    ASSSEG0(*value_len < 10000000, "Invalid array field format"); // check that *value_len didn't underflow beneath 0
}

void sam_seg_idx_aux (VBlockSAMP vb)
{
    if ((int32_t)vb->n_auxs < 1)
        return; // this line has no AUX fields (possibly negative)

    bool is_bam = IS_BAM_ZIP;

    for (int16_t f = 0; f < vb->n_auxs; f++) {
        ASSSEG (vb->aux_lens[f] > (is_bam ? 3 : 5), "Invalid auxilliary field format. AUX tag #%u (0-based) has a length of %u, expecting %u: \"%.*s\"",
                f, vb->aux_lens[f], (is_bam ? 3 : 5), MIN_(16, vb->aux_lens[f]), vb->auxs[f]);

        char c1 = vb->auxs[f][0];
        char c2 = vb->auxs[f][1];
        char c3 = is_bam ? sam_seg_bam_type_to_sam_type (vb->auxs[f][2]) : vb->auxs[f][3];

        #define AUXval(c1,c2,c3) (((uint32_t)(c1)) << 16 | ((uint32_t)(c2)) << 8 | ((uint32_t)(c3)))

        #define TEST_AUX(name, c1, c2, c3)        \
            case AUXval(c1,c2,c3): vb->idx_##name = f; break;

        switch (AUXval(c1, c2, c3)) {
            TEST_AUX(NM_i, 'N', 'M', 'i');
            TEST_AUX(NH_i, 'N', 'H', 'i');
            TEST_AUX(MD_Z, 'M', 'D', 'Z');
            TEST_AUX(SA_Z, 'S', 'A', 'Z');
            TEST_AUX(HI_i, 'H', 'I', 'i');
            TEST_AUX(IH_i, 'I', 'H', 'i');
            TEST_AUX(XG_Z, 'X', 'G', 'Z');
            TEST_AUX(XM_Z, 'X', 'M', 'Z');
            TEST_AUX(X0_i, 'X', '0', 'i');
            TEST_AUX(X1_i, 'X', '1', 'i');
            TEST_AUX(XA_Z, 'X', 'A', 'Z');
            TEST_AUX(AS_i, 'A', 'S', 'i');
            TEST_AUX(CC_Z, 'C', 'C', 'Z');
            TEST_AUX(CP_i, 'C', 'P', 'i');
            TEST_AUX(ms_i, 'm', 's', 'i');
            TEST_AUX(SM_i, 'S', 'M', 'i');
            TEST_AUX(UB_Z, 'U', 'B', 'Z');
            TEST_AUX(BX_Z, 'B', 'X', 'Z');
            TEST_AUX(CB_Z, 'C', 'B', 'Z');
            TEST_AUX(CR_Z, 'C', 'R', 'Z');
            TEST_AUX(CY_Z, 'C', 'Y', 'Z');
            TEST_AUX(XO_Z, 'X', 'O', 'Z');
            TEST_AUX(YS_Z, 'Y', 'S', 'Z');
            TEST_AUX(XB_A, 'X', 'B', 'A');
            TEST_AUX(XB_Z, 'X', 'B', 'Z');
            TEST_AUX(dq_Z, 'd', 'q', 'Z');
            TEST_AUX(iq_Z, 'i', 'q', 'Z');
            TEST_AUX(sq_Z, 's', 'q', 'Z');
            default: {}
        }
    }
}

// returns aux field if it exists or NULL if it doesn't
// In BAM, with a numeric field, result is returned in numeric, otherwise in the return value + value_len
static rom sam_seg_get_aux (VBlockSAMP vb, int16_t idx, uint32_t *value_len, ValueType *numeric, bool is_bam)
{
    STR0(my_value);
    rom tag;
    char sam_type = 0, bam_type = 0, array_subtype = 0;

    if (is_bam)
        bam_get_one_aux (vb, idx, &tag, &bam_type, &array_subtype, pSTRa (my_value), numeric);
    else
        sam_get_one_aux (vb, idx, &tag, &sam_type, &array_subtype, pSTRa (my_value));

    if (value_len)
        *value_len = my_value_len;

    return (is_bam && !my_value) ? &vb->auxs[idx][3] : my_value; // in BAM, if NULL, value is in numeric
}

// returns the length of the field, or 0 if it is not found
// note: BAM allows aux ints up to 0xffffffff, but this function supports only up to 0x7fffffff. Don't
// use for fields that might have a larger value
uint32_t sam_seg_get_aux_int (VBlockSAMP vb, int16_t idx,
                              int32_t *number, // modified only if integer is parsed
                              bool is_bam,
                              int32_t min_value, int32_t max_value, FailType soft_fail)
{
    ASSERT (min_value <= max_value, "%s: expecting idx=%u expecting min_value=%d <= max_value=%d", LN_NAME, idx, min_value, max_value);

    ValueType numeric;
    uint32_t value_len;
    rom value = sam_seg_get_aux (vb, idx, &value_len, &numeric, is_bam);

    if (value && is_bam) {
        if (numeric.i < min_value || numeric.i > max_value) goto out_of_range;

        *number = numeric.i;
        return value_len;
    }

    else if (value && str_get_int (STRa (value), &numeric.i)) {
        if (numeric.i < min_value || numeric.i > max_value) goto out_of_range;

        *number = numeric.i;
        return value_len;
    }

    // this line doesn't have this field, or (SAM only) the field is not a valid integer
    if (soft_fail) return 0;
    ABORT("%s: no valid value found for %.2s", LN_NAME, vb->auxs[idx]);

out_of_range:
    if (soft_fail) return 0;
    ABORT("%s: value of %.2s=%" PRId64 " is out of range [%d,%d]", LN_NAME, vb->auxs[idx], numeric.i, min_value, max_value);
}

void sam_seg_get_aux_Z(VBlockSAMP vb, int16_t idx, pSTRp (snip), bool is_bam)
{
    *snip = vb->auxs[idx] + (is_bam ? 3 : 5);
    *snip_len = vb->aux_lens[idx] - (is_bam ? 4 : 5); // bam: remove the count of the \0 that is including in aux_lens[idx] in BAM
}

char sam_seg_get_aux_A (VBlockSAMP vb, int16_t idx, bool is_bam)
{
    return *(vb->auxs[idx] + (is_bam ? 3 : 5));
}

// note: we test RNAME but not POS. POS exceeding contig LN is handled in ref_seg_get_locked_range_loaded
void sam_seg_verify_RNAME (VBlockSAMP vb)
{
    if (segconf.running)
        return;

    if (IS_REF_INTERNAL && (!sam_hdr_contigs /* SQ-less SAM */ || !sam_hdr_contigs->contigs.len /* SQ-less BAM */)) return;
    if (IS_ASTERISK(vb->chrom_name)) return; // unaligned

    // Per SAM spec, if the header has any SQ lines, all RNAMEs must be listed in it (the header contigs appear first in CTX(RNAME), see sam_zip_initialize
    if (sam_hdr_contigs && sam_hdr_contigs->contigs.len)
        ASSSEG (vb->chrom_node_index < sam_hdr_contigs->contigs.len, "RNAME \"%.*s\" does not have an SQ record in the header", STRf (vb->chrom_name));

    // case an external reference and no SQ lines: if an RNAME is missing in the reference file, we compressed the alignment as "unaligned"
    else {  // headerless SAM
        WordIndex ref_index = chrom_2ref_seg_get (gref, VB, vb->chrom_node_index); // possibly an alt contig
        if (ref_index == WORD_INDEX_NONE) {
            WARN_ONCE("FYI: RNAME \"%.*s\" (and possibly others) is missing in the reference file. This might impact the compression ratio.",
                      vb->chrom_name_len, vb->chrom_name);
            return; // the sequence will be segged as unaligned
        }
    }
}

static void sam_seg_get_kraken (VBlockSAMP vb, bool *has_kraken,
                                char *taxid_str, rom *tag, char *sam_type,
                                pSTRp (value), ValueType *numeric, // out - value for SAM, numeric for BAM
                                bool is_bam)
{
    *tag = "tx"; // genozip introduced tag (=taxid)
    *sam_type = 'i';
    *value = taxid_str;
    *has_kraken = false;
    *value_len = kraken_seg_taxid_do (VB, SAM_TAXID, last_txt (VB, SAM_QNAME), vb->last_txt_len (SAM_QNAME),
                                      taxid_str, true);
    *numeric = CTX(SAM_TAXID)->last_value;

    vb->recon_size += is_bam ? 7 : (*value_len + 6); // txt modified
}

void sam_seg_aux_all (VBlockSAMP vb, ZipDataLineSAM *dl)
{
    START_TIMER;

    const bool is_bam = IS_BAM_ZIP;
    Container con = { .repeats = 1, .filter_repeats = true /* v14 */};
    char prefixes[(vb->n_auxs + 1) * 6 + 3];                // each name is 5 characters per SAM specification, eg "MC:Z:" followed by CON_PX_SEP ; +3 for the initial CON_PX_SEP. +1 for kraken
    prefixes[0] = prefixes[1] = prefixes[2] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix followed by separator for empty prefix for translator-only item[0]
    unsigned prefixes_len = 3;

    // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM
    con.items[con_nitems (con)] = (ContainerItem){ .translator = SAM2BAM_AUX_SELF };
    con_inc_nitems (con);

    bool has_kraken = kraken_is_loaded;

    for (int16_t idx = 0; idx < vb->n_auxs + has_kraken; idx++) {

        STR0(value);
        ValueType numeric = {};
        rom tag;
        char sam_type = 0, bam_type = 0, array_subtype = 0;
        char taxid_str[20];

        if (idx == vb->n_auxs)
            sam_seg_get_kraken (vb, &has_kraken, taxid_str, &tag, &bam_type, pSTRa (value), &numeric, is_bam);
        else if (is_bam)
            bam_get_one_aux (vb, idx, &tag, &bam_type, &array_subtype, pSTRa (value), &numeric);
        else
            sam_get_one_aux (vb, idx, &tag, &sam_type, &array_subtype, pSTRa (value));

        if (sam_type == 'i') {
            if (idx == vb->idx_NM_i) // we already converted NM to integer, no need to do it again
                numeric.i = dl->NM;
            else
                ASSERT(str_get_int (STRa (value), &numeric.i), "%s: Expecting integer value for auxiliary field %c%c but found \"%.*s\"",
                       LN_NAME, tag[0], tag[1], STRf (value));
            value = 0;
        }

        if (!bam_type) bam_type = sam_seg_sam_type_to_bam_type (sam_type, numeric.i);
        if (!sam_type) sam_type = sam_seg_bam_type_to_sam_type (bam_type);

        ASSERT (bam_type, "%s: value %" PRId64 " of field %c%c is out of range of the BAM specification: [%d-%u]",
                LN_NAME, numeric.i, tag[0], tag[1], -0x80000000, 0x7fffffff);

        con.items[con_nitems (con)] = (ContainerItem){
            .dict_id    = sam_seg_aux_field (vb, dl, is_bam, tag, bam_type, array_subtype, STRa (value), numeric),
            .translator = aux_field_translator((uint8_t)bam_type), // how to transform the field if reconstructing to BAM
            .separator  = {aux_sep_by_type[is_bam][(uint8_t)bam_type], '\t'},
        };

        con_inc_nitems (con);

        ASSSEG (con_nitems (con) <= MAX_FIELDS, "too many optional fields, limit is %u", MAX_FIELDS);

        // note: in the optional field prefix (unlike array type), all integer types become 'i'.
        char prefix[6] = {tag[0], tag[1], ':', sam_type, ':', CON_PX_SEP};
        memcpy(&prefixes[prefixes_len], prefix, 6);
        prefixes_len += 6;
    }

    uint32_t num_items = con_nitems (con);
    if (num_items > 1) {                                                                  // we have Aux fields, not just the translator item
        if (con.items[num_items - 1].separator[0] & 0x80)                                 // is a flag
            con.items[num_items - 1].separator[0] &= ~(CI0_NATIVE_NEXT & ~(uint8_t)0x80); // last Aux field has no tab
        con.items[num_items - 1].separator[1] = 0;
        container_seg (vb, CTX(SAM_AUX), &con, prefixes, prefixes_len, (is_bam ? 3 : 5) * (num_items - 1)); // account for : SAM: "MX:i:" BAM: "MXi"
    }

    else
        // NULL means MISSING Container item (of the toplevel container) - will cause container_reconstruct_do of
        // the toplevel container to delete of previous separator (\t)
        container_seg (vb, CTX(SAM_AUX), 0, 0, 0, 0);

    COPY_TIMER(sam_seg_aux_all);
}   

static inline bool has_same_qname (VBlockSAMP vb, STRp (qname), LineIType buddy_line_i, bool last_char_flipped)
{
    if (buddy_line_i == NO_LINE) return false;

    TxtWord buddy_q = DATA_LINE(buddy_line_i)->QNAME;

    if (!last_char_flipped) 
        return str_issame_(STRa(qname), STRtxtw (buddy_q));

    else {   
        char r_a = qname[qname_len-1];
        char r_b = *Btxt(buddy_q.index + buddy_q.len - 1);

        return str_issame_(STRa(qname) - 1, STRtxtw (buddy_q) - 1) &&
               ((r_a=='1' && r_b=='2') || (r_a=='2' && r_b=='1'));
    }
}

#define LINE_BY_HASH(hash) *B(LineIType, vb->qname_hash, (hash))

// seg mate as buddy and return true if this line has one
static inline BuddyType sam_seg_mate (VBlockSAMP vb, SamFlags f, STRp (qname), uint32_t my_hash, bool *insert_to_hash)
{
    if (sam_is_depn (f) || !segconf.is_paired) return false;

    uint32_t mate_hash = qname_calc_hash (QNAME1, STRa(qname), !f.is_last, true, NULL) & MAXB(vb->qname_hash.prm8[0]);
    LineIType candidate = LINE_BY_HASH(mate_hash);
    SamFlags mate_f = DATA_LINE(candidate)->FLAG;

    // case: mate is found
    if (has_same_qname (vb, STRa (qname), candidate, segconf.flav_prop[QNAME1].has_R) && 
        !sam_is_depn (mate_f) && 
        mate_f.is_last != f.is_last) {
        
        vb->mate_line_i = candidate;
        vb->mate_line_count++; // for stats

        seg_add_to_local_resizable (VB, CTX(SAM_BUDDY), vb->line_i - candidate, 0); // add buddy (delta) >= 1 .
        return BUDDY_MATE;
    }

    // case: mate not found
    else {
        // if we haven't found our mate - store this line in the hash table so our mate can find us (don't store
        // unneccessarily to reduce hash contention)
        *insert_to_hash = true;
        return BUDDY_NONE;
    }
}

// seg saggy (a previous line that's member of the same sag) as buddy and return true if this line has one
static inline BuddyType sam_seg_saggy (VBlockSAMP vb, SamFlags f, STRp (qname), uint32_t my_hash, bool *insert_to_hash)
{
    if (!IS_MAIN(vb)) return false;

    LineIType candidate = LINE_BY_HASH(my_hash);
    SamFlags saggy_f = DATA_LINE(candidate)->FLAG;

    // case: we found another member of the same sag (i.e. same qname, same is_last)
    if (has_same_qname (vb, STRa (qname), candidate, false) && saggy_f.is_last == f.is_last) { // the "prim" line against which we are segging cannot have hard clips
        vb->saggy_line_i = candidate;
        vb->saggy_is_prim = !sam_is_depn (saggy_f);
        vb->saggy_near_count++;

        // replace value in hash table with current line, so future lines from the same seg have a smaller buddy delta,
        // except if the line is a prim - we prefer the saggy to be a prim because there are many methods that require it
        if (!vb->saggy_is_prim)
            *insert_to_hash = true;

        seg_add_to_local_resizable (VB, CTX(SAM_BUDDY), vb->line_i - candidate, 0); // add buddy (delta) >= 1 .
        return BUDDY_SAGGY;
    }

    // case: we didn't find another member of our sag in this VB
    else {
        *insert_to_hash = true;
        return BUDDY_NONE;
    }
}

static inline void sam_seg_QNAME_segconf (VBlockSAMP vb, ContextP ctx, STRp (qname))
{
    if (segconf.is_collated) { // still collated, try to find evidence to the contrary
        bool is_new = !is_same_last_txt (VB, ctx, STRa(qname));
        if (is_new && ctx->last_is_new) segconf.is_collated = false; // two new QNAMEs in a row = not collated

        // case: at least on pair of consecutive lines has the same QNAME. if the file is not sorted, we will set it as collated in sam_seg_finalize_segconf
        if (!is_new) segconf.evidence_of_collated = true;

        ctx->last_is_new = is_new;
    }

    if (vb->line_i == 0) {
        qname_segconf_discover_flavor (VB, QNAME1, STRa (qname));

        if (flag.deep) 
            memcpy (segconf.deep_1st_qname, qname, MIN_(qname_len, sizeof (segconf.deep_1st_qname)-1));
    }
}

void sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp (qname), unsigned add_additional_bytes)
{
    decl_ctx (SAM_QNAME);

    // in segconf, identify if this file is collated (each QNAME appears in two or more consecutive lines)
    if (segconf.running) {
        sam_seg_QNAME_segconf (vb, ctx, STRa (qname));
        goto normal_seg;
    }

    uint32_t qname_hash = qname_calc_hash (QNAME1, qname, qname_len, dl->FLAG.is_last, true, NULL); // note: canonical=true as we use the same hash for find a mate and a saggy
    uint32_t my_hash = qname_hash & MAXB(vb->qname_hash.prm8[0]);
    bool insert_to_hash = false;

    if (flag.deep || flag.show_deep == 2) 
        sam_deep_set_QNAME_hash (vb, dl, STRa(qname));

    BuddyType bt = sam_seg_mate  (vb, dl->FLAG, STRa (qname), my_hash, &insert_to_hash) | // bitwise or
                   sam_seg_saggy (vb, dl->FLAG, STRa (qname), my_hash, &insert_to_hash);

    if (bt && flag.show_buddy) {
        if (bt == BUDDY_EITHER)
            iprintf("%s: %.*s mate_line_i=%u saggy_line_i=%u\n", LN_NAME, STRf (qname), vb->mate_line_i, vb->saggy_line_i);
        else if (bt == BUDDY_MATE)
            iprintf("%s: %.*s mate_line_i=%u\n", LN_NAME, STRf (qname), vb->mate_line_i);
        else if (bt == BUDDY_SAGGY)
            iprintf("%s: %.*s saggy_line_i=%u\n", LN_NAME, STRf (qname), vb->saggy_line_i);
    }

    if (insert_to_hash)
        LINE_BY_HASH(my_hash) = vb->line_i;

    // note: we add the buddy_type to qname rather than to SAM_BUDDY, as QNAME carries part of the buddy entropy anyway (i.e. has buddy vs hasn't)

    // case: if QNAME is the same as buddy's - seg against buddy
    // note: in PRIM segging against buddy here is used for loading sag, SAM_QNAMESA is used for reconstruction
    if (bt)
        seg_by_ctx (VB, (char[]){SNIP_SPECIAL, SAM_SPECIAL_COPY_BUDDY, '0' + bt}, 3, ctx, qname_len + add_additional_bytes); // seg QNAME as copy-from-buddy

    // case: DEPN with SA Group: seg against SA Group (unless already segged against buddy)
    else if (IS_DEPN(vb) && vb->sag)
        sam_seg_against_sa_group (vb, ctx, qname_len + add_additional_bytes);

    else normal_seg:
        qname_seg (VB, QNAME1, STRa (qname), add_additional_bytes); // note: for PRIM component, this will be consumed with loading SA

    // case: PRIM: additional seg against SA Group - store in SAM_QNAMESA - Reconstruct will take from here in PRIM per Toplevel container
    if (IS_PRIM(vb))
        seg_by_did (VB, (char[]){SNIP_SPECIAL, SAM_SPECIAL_PRIM_QNAME}, 2, SAM_QNAMESA, 0); // consumed when reconstructing PRIM vb
}

WordIndex sam_seg_RNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp (chrom),
                         bool against_sa_group_ok, // if true, vb->chrom_node_index must already be set
                         unsigned add_bytes)
{
    bool normal_seg = false;

    if (segconf.running)
        goto normal_seg;

    WordIndex node_index;

    // case: PRIM or DEPN vb - seg against SA group with alignments
    // Note: in DEPN, rname already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    if (against_sa_group_ok && sam_seg_has_sag_by_SA (vb)) {
        sam_seg_against_sa_group (vb, CTX(SAM_RNAME), add_bytes);

        if (IS_PRIM(vb)) {
            // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
            seg_by_did (VB, STRa (chrom), OPTION_SA_RNAME, 0);

            // count RNAME field contribution to OPTION_SA_RNAME, so sam_stats_reallocate can allocate the z_data between RNAME and SA:Z
            CTX(OPTION_SA_RNAME)->counts.count += add_bytes;
        }

        STRset (vb->chrom_name, chrom);
        random_access_update_chrom (VB, 0, vb->chrom_node_index, STRa (chrom));

        node_index = vb->chrom_node_index;
    }

    // seg RNAME against mate's RNEXT. limit to only if RNAME/RNEXT are different (if the same, then
    // likely all the lines between this line and the mate have the same RNAME)
    // note: for now, this only works in BAM because dl->RNAME/RNEXT are set already. To do: support in SAM.
    else if (IS_BAM_ZIP && sam_has_mate &&
             dl->RNAME != dl->RNEXT &&
             dl->RNAME == DATA_LINE(vb->mate_line_i)->RNEXT &&
             DATA_LINE(vb->mate_line_i)->RNAME != DATA_LINE(vb->mate_line_i)->RNEXT) { // protect from pathological case (observed in test.longranger-wgs.bam) where mate_RNAME==mate_RNEXT==my_RNAME but my_RNEXT=*

        seg_by_did (VB, STRa (copy_mate_RNEXT_snip), SAM_RNAME, add_bytes); // copy POS from earlier-line mate PNEXT
        STRset (vb->chrom_name, chrom);
        random_access_update_chrom (VB, 0, dl->RNAME, STRa (chrom));
        node_index = dl->RNAME;
    }

    else normal_seg: {
        bool is_new;
        node_index = chrom_seg_ex (VB, SAM_RNAME, STRa (chrom), 0, NULL, add_bytes, !IS_BAM_ZIP, &is_new);
        normal_seg = true;

        // don't allow adding chroms to a BAM file or a SAM that has SQ lines in the header, but we do allow to add to a headerless SAM.
        ASSSEG (!is_new || !sam_hdr_contigs || segconf.running || IS_EQUAL_SIGN(chrom) || IS_ASTERISK(chrom),
                "contig '%.*s' appears in file, but is missing in the %s header", STRf (chrom), dt_name (vb->data_type));
    }

    // protect rname from removal by ctx_shorten_unused_dict_words if we didn't seg normally
    if (!normal_seg && node_index != WORD_INDEX_NONE) // note: node_index==-1 when RNAME="*"
        ctx_protect_from_removal (VB, CTX(SAM_RNAME), node_index);

    return node_index;
}

WordIndex sam_seg_RNEXT (VBlockSAMP vb, ZipDataLineSAM *dl, STRp (chrom), unsigned add_bytes)
{
    bool normal_seg = false;

    if (segconf.running) {
        if (chrom_len > 1 || *chrom != '*') segconf.has[SAM_RNEXT] = true;
        goto normal_seg;
    }

    WordIndex node_index = IS_BAM_ZIP ? dl->RNEXT : WORD_INDEX_NONE; // already set in case of BAM

    // RNEXT was not detected in segconf - we seg plainly, expecting that is most cases all lines will be "*"
    if (!segconf.has[SAM_RNEXT]) goto normal_seg;

    // case: seg against mate's RNAME. limit to only if RNAME/RNEXT are different (if the same, then
    // likely all the lines between this line and the mate have the same RNAME)
    // note: for now, this only works in BAM because dl->RNAME/RNEXT are set already. To do: support in SAM.
    else if (IS_BAM_ZIP && sam_has_mate &&
             dl->RNAME != dl->RNEXT && dl->RNEXT == DATA_LINE(vb->mate_line_i)->RNAME)
        seg_by_did (VB, STRa (copy_mate_RNAME_snip), SAM_RNEXT, add_bytes); // copy POS from earlier-line mate PNEXT

    // case: seg RNEXT against prim's RNAME. This happens when RNEXT is the same as prim's RNEXT, but prim's
    // RNEXT is the same as prim's RNAME, so PRIM's RNEXT is segged as "=".
    else if (IS_BAM_ZIP && sam_has_saggy &&
             dl->RNAME != dl->RNEXT && dl->RNEXT == DATA_LINE(vb->saggy_line_i)->RNAME)
        seg_by_did (VB, STRa (copy_saggy_RNAME_snip), SAM_RNEXT, add_bytes);

    // case: seg RNEXT against prim's RNEXT
    else if (IS_BAM_ZIP && sam_has_saggy &&
             dl->RNAME != dl->RNEXT && dl->RNEXT == DATA_LINE(vb->saggy_line_i)->RNEXT)
        seg_by_did (VB, (char[]){SNIP_SPECIAL, SAM_SPECIAL_COPY_BUDDY, '0' + BUDDY_SAGGY}, 3, SAM_RNEXT, add_bytes);

    else normal_seg: {
        bool is_new;
        node_index = chrom_seg_ex (VB, SAM_RNEXT, STRa (chrom), 0, NULL, add_bytes, !IS_BAM_ZIP, &is_new);
        normal_seg = true;

        // don't allow adding chroms to a BAM file or a SAM that has SQ lines in the header, but we do allow to add to a headerless SAM.
        ASSSEG (!is_new || !sam_hdr_contigs || segconf.running || IS_EQUAL_SIGN(chrom) || IS_ASTERISK(chrom),
               "contig '%.*s' appears in file, but is missing in the %s header", STRf (chrom), dt_name (vb->data_type));
    }

    // protect rnext from removal by ctx_shorten_unused_dict_words if we didn't seg normally
    if (!normal_seg && node_index != WORD_INDEX_NONE) // note: node_index==-1 when RNEXT="*"
        ctx_protect_from_removal (VB, CTX(SAM_RNEXT), node_index);

    return node_index;
}

bool sam_seg_test_biopsy_line (VBlockP vb, STRp (line))
{
    if (segconf.running) return false; // we need to let segconf run normally, so we get the correct VB size

    if (flag.biopsy_line.line_i == vb->line_i && flag.biopsy_line.vb_i == vb->vblock_i) {
        PutLineFn fn = file_put_line (VB, STRa (line), "Line biopsy:");

        if (TXT_DT(BAM)) WARN("Tip: You can view the dumped BAM line with:\n   genozip --show-bam %s", fn.s);
        exit_ok;
    }

    vb->recon_size -= line_len;
    return true;
}

void sam_seg_init_bisulfite (VBlockSAMP vb, ZipDataLineSAM *dl)
{
    // the converted reference to which this read was mapped (C->T conversion or G->A conversion)
    // note: we calculate it always to avoid needless adding entropy in the snip
    vb->bisulfite_strand =  !segconf.sam_bisulfite    ? 0 
                            : IS_REF_INTERNAL         ? 0 // bug 648
                            : MP(BISMARK)   && has(XG_Z) ? sam_seg_get_aux_A (vb, vb->idx_XG_Z, IS_BAM_ZIP)
                            : MP(DRAGEN)    && has(XG_Z) ? sam_seg_get_aux_A (vb, vb->idx_XG_Z, IS_BAM_ZIP)
                            : MP(BSSEEKER2) && has(XO_Z) ? "CG"[sam_seg_get_aux_A (vb, vb->idx_XO_Z, IS_BAM_ZIP) == '-']
                            : MP(BSBOLT)    && has(YS_Z) ? "CG"[sam_seg_get_aux_A (vb, vb->idx_YS_Z, IS_BAM_ZIP) == 'C']
                            : MP(GEM3)      && has(XB_A) ? sam_seg_get_aux_A (vb, vb->idx_XB_A, IS_BAM_ZIP)
                            :                           0;

    // enter the converted bases into the reference, in case of REF_INTERNAL 
    if (segconf.running || !vb->bisulfite_strand) return;

    else if (MP(BSSEEKER2))
        sam_seg_bsseeker2_XG_Z_analyze (vb, dl, STRauxZ(XG_Z, false), dl->POS);

    else if (MP(BISMARK) || MP(DRAGEN))
        sam_seg_bismark_XM_Z_analyze (vb, dl);
        
    else if (MP(BSBOLT))
        sam_seg_bsbolt_XB_Z_analyze (vb, dl);
}

rom sam_seg_txt_line (VBlockP vb_, rom next_line, uint32_t remaining_txt_len, bool *has_13)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    sam_reset_line (VB);

    ASSERT (!((remaining_txt_len==1 && *next_line=='\n') || (remaining_txt_len==2 && next_line[1]=='\n')),
            "%s: Bad SAM file - empty line detected", LN_NAME);

    WordIndex prev_line_chrom = vb->chrom_node_index;

    ZipDataLineSAM *dl = DATA_LINE(vb->line_i);

    str_split_by_tab (next_line, remaining_txt_len, MAX_FIELDS + AUX, has_13, false, true); // also advances next_line to next line
    
    ASSSEG (n_flds >= 11, "%s: Bad SAM file: alignment expected to have at least 11 fields, but found only %u", LN_NAME, n_flds);

    vb->n_auxs = n_flds - AUX;
    vb->auxs = &flds[AUX]; // note: pointers to data on the stack
    vb->aux_lens = &fld_lens[AUX];

    // support non-compliant SAM by e.g. ngmlr v0.2.7 - lines end with \t\n
    bool terminal_tab = vb->n_auxs && vb->aux_lens[vb->n_auxs-1] == 0;
    if (terminal_tab) vb->n_auxs--;
    ASSSEG0 (!terminal_tab || ! *has_13, "line ends with \\t\\r\\n: this is not currently supported by Genozip");

    sam_seg_idx_aux (vb);

    dl->SEQ.index = BNUMtxt (flds[SEQ]); // SEQ.len to be determined by sam_cigar_analyze
    dl->QNAME = TXTWORDi (fld, QNAME);
    dl->QUAL  = TXTWORDi (fld, QUAL);

    if (fld_lens[QUAL] == 1 && flds[QUAL][0] == '*')
        vb->qual_missing = dl->no_qual = true;

    // lazy way to get vb->chrom* (inc. if --match-chrom), rollback later if seg RNAMEaa is not needed
    if (vb->check_for_gc || !IS_MAIN(vb))
        seg_create_rollback_point (VB, NULL, 1, SAM_RNAME);

    int32_t rname_shrinkage = vb->recon_size;
    dl->RNAME = sam_seg_RNAME (vb, dl, STRfld (RNAME), false, fld_lens[RNAME] + 1);
    rname_shrinkage -= vb->recon_size; // number of characters RNAME was reduced by, due to --match-chrom-to-reference

    ASSSEG (str_get_int_range32 (STRfld (POS), 0, MAX_POS_SAM, &dl->POS),
            "Invalid POS \"%.*s\": expecting an integer [0,%d]", STRfi (fld, POS), (int)MAX_POS_SAM);

    ASSSEG (str_get_int_range8 (STRfld (MAPQ), 0, 255, &dl->MAPQ),
            "Invalid MAPQ \"%.*s\": expecting an integer [0,255]", STRfi (fld, MAPQ));

    ASSSEG (str_get_int_range16 (STRfld (FLAG), 0, SAM_MAX_FLAG, &dl->FLAG.value), 
            "Invalid FLAG field: \"%.*s\"", STRfi (fld, FLAG));

    // analyze (but don't seg yet) CIGAR
    uint32_t seq_consumed;
    sam_cigar_analyze (vb, STRfld (CIGAR), false, &seq_consumed);
    dl->SEQ.len = seq_consumed; // do it this way to avoid compiler warning

    ASSSEG (dl->SEQ.len == fld_lens[SEQ] || (fld_lens[CIGAR] == 1 && *flds[CIGAR] == '*') || flds[SEQ][0] == '*', 
            "seq_len implied by CIGAR=%.*s is %u, but actual SEQ length is %u, SEQ=%.*s",
            STRfi (fld, CIGAR), dl->SEQ.len, fld_lens[SEQ], STRfi (fld, SEQ));

    // if this is a secondary / supplamentary read (aka Dependent) or a read that has an associated sec/sup read (aka Primary) - move
    // the line to the appropriate component and skip it here (no segging done yet)
    if (vb->check_for_gc && sam_seg_is_gc_line (vb, dl, flds[0], next_line - flds[0], false)) {
        vb->debug_line_hash_skip = true;
        goto rollback_and_done;
    }

    // case: biopsy (only arrives here in MAIN VBs if gencomp) 
    //       or biopsy_line: we just needed to pass sam_seg_is_gc_line and we're done
    if ((flag.biopsy && !segconf.running) ||
        (flag.biopsy_line.line_i != NO_LINE && sam_seg_test_biopsy_line (VB, flds[0], next_line - flds[0])))
        goto rollback_and_done;

    vb->last_cigar = flds[CIGAR];
    SAFE_NUL(&vb->last_cigar[fld_lens[CIGAR]]); // nul-terminate CIGAR string

    if (has(NM_i))
        dl->NM_len = sam_seg_get_aux_int (vb, vb->idx_NM_i, &dl->NM, false, MIN_NM_i, MAX_NM_i, HARD_FAIL) + 1; // +1 for \t or \n

    if (!IS_MAIN(vb)) {
        // set dl->AS needed by sam_seg_prim_add_sag
        if (IS_PRIM(vb) && has(AS_i))
            sam_seg_get_aux_int (vb, vb->idx_AS_i, &dl->AS, false, MIN_AS_i, MAX_AS_i, HARD_FAIL);

        sam_seg_sag_stuff (vb, dl, STRfld (CIGAR), flds[SEQ], false);

        // re-seg rname, against SA group
        seg_rollback (VB);
        vb->recon_size += rname_shrinkage; // grow back, as we will shrink it again when re-segging
        sam_seg_RNAME(vb, dl, STRfld (RNAME), true, fld_lens[RNAME] + 1);
    }

    sam_seg_QNAME (vb, dl, STRfld (QNAME), 1);

    sam_seg_FLAG (vb, dl, fld_lens[FLAG] + 1);

    // note: pos can have a value even if RNAME="*" - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    sam_seg_POS (vb, dl, prev_line_chrom, fld_lens[POS] + 1);

    if (fld_lens[RNAME] != 1 || *flds[RNAME] != '*')
        sam_seg_verify_RNAME (vb);

    sam_seg_MAPQ (vb, dl, fld_lens[MAPQ] + 1);

    dl->RNEXT = sam_seg_RNEXT (vb, dl, STRfld (RNEXT), fld_lens[RNEXT] + 1);

    vb->RNEXT_is_equal = str_is_1chari (fld, RNEXT, '=') || str_issame_(STRfld(RNEXT), STRfld(RNAME));

    sam_seg_PNEXT (vb, dl, STRfld (PNEXT), 0, fld_lens[PNEXT] + 1);

    sam_seg_init_bisulfite (vb, dl);

    // we search forward for MD:Z now, XG:Z as we will need it for SEQ if it exists
    if (has_MD)
        sam_seg_MD_Z_analyze (vb, dl, STRauxZ(MD_Z, false), dl->POS);

    seg_set_last_txt (VB, CTX (SAM_SQBITMAP), STRfld (SEQ));

    // calculate diff vs. reference (denovo or loaded)
    sam_seg_SEQ (vb, dl, STRfld (SEQ), fld_lens[SEQ] + 1);

    ASSSEG (str_is_in_range (flds[QUAL], fld_lens[QUAL], 33, 126), "Invalid QUAL - it contains non-Phred characters: \"%.*s\"",
           STRfi (fld, QUAL));

    ASSSEG (fld_lens[SEQ] == fld_lens[QUAL] || vb->qual_missing, "Expecting QUAL(length=%u) to be of the same length as SEQ(length=%u)",
            fld_lens[QUAL], fld_lens[SEQ]);

    // finally we can seg CIGAR now
    sam_seg_CIGAR (vb, dl, fld_lens[CIGAR], STRfld (SEQ), STRfld (QUAL), fld_lens[CIGAR] + 1 /*\t*/);

    // note: can only be called after sam_seg_CIGAR updates SEQ.len
    sam_seg_QUAL (vb, dl, STRfld (QUAL), fld_lens[QUAL] + 1);

    // add BIN so this file can be reconstructed as BAM
    bam_seg_BIN (vb, dl, 0, false);

    // AUX fields - up to MAX_FIELDS of them
    sam_seg_aux_all (vb, dl);

    if (IS_PRIM(vb)) {
        if (IS_SAG_NH)
            sam_seg_prim_add_sag_NH (vb, dl, dl->NH);
        else if (IS_SAG_CC)
            sam_seg_prim_add_sag_CC (vb, dl, dl->NH);
        else if (IS_SAG_FLAG)
            sam_seg_prim_add_sag (vb, dl, 0, false);
        else if (IS_SAG_SOLO)
            sam_seg_prim_add_sag_SOLO (vb, dl);
    }

    // TLEN - must be after AUX as we might need data from MC:Z
    sam_seg_TLEN (vb, dl, STRfld (TLEN), 0, vb->RNEXT_is_equal);

    if (terminal_tab) seg_by_did (VB, "\t\n", 2, SAM_EOL, 1); // last field accounted for \n
    else              SEG_EOL (SAM_EOL, false);

    SAFE_RESTORE;            // restore \t after CIGAR

    if (dl->SEQ.len > vb->longest_seq_len) vb->longest_seq_len = dl->SEQ.len;

    return next_line;

rollback_and_done:
    memset (dl, 0, sizeof (ZipDataLineSAM));
    seg_rollback (VB); // cancelling segging of RNAME
    vb->recon_size += rname_shrinkage; // grow back, as the change in recon_size due to the change in RNAME, will be accounted for by the PRIM/DEPN VB to which this alignment belongs was sent

    return next_line;
}
