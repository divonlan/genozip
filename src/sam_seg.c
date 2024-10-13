// ------------------------------------------------------------------
//   sam_seg.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "sam_private.h"
#include "refhash.h"
#include "random_access.h"
#include "codec.h"
#include "stats.h"
#include "chrom.h"
#include "qname.h"
#include "lookback.h"
#include "gencomp.h"
#include "tip.h"
#include "deep.h"
#include "arch.h"
#include "threads.h"
#include "zip_dyn_int.h"
#include "huffman.h"

STRl(copy_GX_snip, 30);
STRl(copy_sn_snip, 30);
STRl(copy_POS_snip, 30);
STRl(copy_Q1NAME_int, 30);
STRl(copy_Q2NAME_int, 30);
STRl(copy_Q3NAME_int, 30);
STRl(copy_mate_CIGAR_snip, 30);
STRl(copy_mate_MAPQ_snip, 30);
STRl(copy_mate_MQ_snip, 30);
STRl(copy_mate_PQ_snip, 30);
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
STRl(copy_mate_rb_snip, 30);
STRl(copy_mate_mb_snip, 30);
STRl(copy_mate_GP_snip, 30);
STRl(copy_mate_MP_snip, 30);
STRl(redirect_to_CR_X_snip, 30);
STRl(redirect_to_GR_X_snip, 30);
STRl(redirect_to_GY_X_snip, 30);
STRl(redirect_to_RX_X_snip, 30);
STRl(XC_snip, 48);
STRl_ARRAY(copy_buddy_Z_snip, NUM_MATED_Z_TAGS, 30);

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

static bool is_headerful_sam (STRp(header))
{
    #define H(x) !memcmp (header, "@" #x "\t", 4)

    return header_len >= 4 && (H(HD) || H(CO) || H(SQ) || H(RG) || H(PG)); 
}

static bool is_headerless_sam (STRp(header))
{
    if (!str_is_printable (STRa(header))) return false;

    // test one line to see that the fields make sense
    str_split_by_tab (header, header_len, MAX_FIELDS + AUX, NULL, false, false); // also advances header
    
    if (!header || n_flds < AUX-1 || !sam_is_cigar (STRfld(CIGAR), false) ||
        !str_is_int (STRfld (FLAG)) || !str_is_int (STRfld (POS)) || !str_is_int (STRfld (MAPQ)) || 
        !str_is_int (STRfld (PNEXT)) || !str_is_int (STRfld (TLEN)))
        return false;

    // AUX fields start with xx:x:
    for (uint32_t i=AUX; i < n_flds; i++)
        if (fld_lens[i] < 5 || flds[i][2] != ':' || flds[i][4] != ':') return false;

    return true;
}

// detect if a generic file is actually a SAM
bool is_sam (STRp(header), bool *need_more)
{
    return is_headerful_sam (STRa(header)) || is_headerless_sam (STRa(header));
}

// main thread, called for each component. called after reading txt header, before segconf.
void sam_zip_initialize (void)
{
    bool has_hdr_contigs = sam_hdr_contigs && sam_hdr_contigs->contigs.len;

    // Copy header contigs to RNAME upon first component. This is in the order of the header, and
    // encodes ref_id based on header order. Note: BAM always has header contigs (might be 0 of them, 
    // for an unaligned file), while SAM is allowed to be header-less (but then requires an external reference)
    if (z_file->num_txts_so_far == 1 && has_hdr_contigs) 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, sam_hdr_contigs);

    qname_zip_initialize();

    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _SAM_POS, false, 0, copy_POS_snip);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_Q1NAME, true, 0, copy_Q1NAME_int);   
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_Q2NAME, true, 0, copy_Q2NAME_int);   
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_Q3NAME, true, 0, copy_Q3NAME_int);   
        seg_prepare_snip_other (SNIP_COPY, _OPTION_GX_Z, false, 0, copy_GX_snip);
        seg_prepare_snip_other (SNIP_COPY, _OPTION_sn_B_f, false, 0, copy_sn_snip);
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
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_PQ_i, copy_mate_PQ_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_ms_i, copy_mate_ms_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_nM_i, copy_mate_nM_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_YS_i, copy_mate_YS_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_AS_i, copy_mate_AS_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_NH_i, copy_buddy_NH_snip, '0' + BUDDY_EITHER);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_rb_Z, copy_mate_rb_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_mb_Z, copy_mate_mb_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_GP_i, copy_mate_GP_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY, _OPTION_MP_i, copy_mate_MP_snip, '0' + BUDDY_MATE);
        seg_prepare_snip_special_other (SAM_SPECIAL_COPY_BUDDY_CIGAR, _SAM_CIGAR, copy_mate_CIGAR_snip, '0' + BUDDY_MATE);
        seg_prepare_plus_snip (SAM, 3, ((DictId[]){ {_OPTION_XX_i}, {_OPTION_YY_i}, {_OPTION_XY_i} }), XC_snip);

        // we seg into mux channels, but we copy from the parent
        for (MatedZFields f = 0; f < NUM_MATED_Z_TAGS; f++)
        {
            copy_buddy_Z_snips[f][0] = SNIP_SPECIAL;
            copy_buddy_Z_snip_lens[f] = sizeof (copy_buddy_Z_snips[f]) - 1;
            seg_prepare_snip_other_do (SAM_SPECIAL_COPY_BUDDY, ZCTX(buddied_Z_dids[f])->dict_id,
                                       true, 0, '0' + BUDDY_EITHER, &copy_buddy_Z_snips[f][1], &copy_buddy_Z_snip_lens[f]);
            copy_buddy_Z_snip_lens[f]++;
        }
    }

    sam_MM_zip_initialize();

    if (MP(ULTIMA)) sam_ultima_zip_initialize();
    
    if (MP(STAR)) sam_star_zip_initialize();

    if (segconf.sam_has_abra2) sam_abra2_zip_initialize();
    
    if (flag.deep && txt_file->name)
        strncpy (segconf.sam_deep_filename, txt_file->name, sizeof(segconf.sam_deep_filename)-1);
}

// called after each txt file (including after the SAM component in Deep)
void sam_zip_finalize (bool is_last_user_txt_file)
{
    threads_log_by_vb (evb, "main_thread", "sam_zip_finalize", 0);

    if (flag.analyze_ins)
        sam_zip_report_monochar_inserts();

    if (!flag.let_OS_cleanup_on_exit) {

        if ((IS_REF_INTERNAL || IS_REF_EXT_STORE) && !flag.deep)
            ref_destroy_reference (gref);

        if (segconf.sag_type) 
            gencomp_destroy();

        if (flag.deep) 
            sam_deep_zip_finalize();
    }
}

// main thread: after finishing writing z_file
void sam_zip_end_of_z (void)
{
    sam_header_finalize();
}

// called by main thread after reading txt of one vb into vb->txt_data
void sam_zip_init_vb (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    vb->chrom_node_index = WORD_INDEX_NONE;
    CTX(OPTION_XG_Z)->XG_inc_S = unknown;

    // note: we test for sorted and not collated, because we want non-sorted long read files (which are collated)
    // to seg depn against same-VB prim (i.e. not gencomp) - as the depn lines will follow the prim line
    vb->check_for_gc = !segconf.running && segconf.sag_type && !segconf.no_gc_checking && IS_MAIN(vb);

    vb->vb_plan.prev_comp_i = SAM_COMP_PRIM; // initialize (same as in sam_piz_init_vb)

    sam_sag_zip_init_vb (vb);
}

// called compute thread after compress, order of VBs is arbitrary
void sam_zip_after_compress (VBlockP vb)
{
}

// called by main thread, as VBs complete (might be out-of-order)
void sam_zip_after_compute (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (IS_PRIM(vb))
        gencomp_sam_prim_vb_has_been_ingested (VB);

    // increment stats accumulators
    z_file->num_perfect_matches += vb->num_perfect_matches; 
    z_file->num_aligned         += vb->num_aligned;
    z_file->secondary_count     += vb->secondary_count;
    z_file->supplementary_count += vb->supplementary_count;
    z_file->mate_line_count     += vb->mate_line_count;
    z_file->saggy_near_count    += vb->saggy_near_count;
    z_file->depn_far_count      += vb->depn_far_count;

    if (segconf.sam_is_unmapped && IS_REF_LOADED_ZIP) { 
        DO_ONCE ref_verify_organism (VB);
    }
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
    header->sam.segconf_SA_HtoS         = segconf.SA_HtoS == yes;        // v14
    header->sam.segconf_SA_CIGAR_abb    = segconf.SA_CIGAR_abbreviated == yes; // 15.0.66
    header->sam.segconf_is_sorted       = segconf.is_sorted;             // v14
    header->sam.segconf_is_collated     = segconf.is_collated;           // v14
    header->sam.segconf_pysam_qual      = segconf.pysam_qual;            // v14
    header->sam.segconf_10xGen          = segconf.has_10xGen;            // v14
    header->sam.segconf_seq_len_dict_id = segconf.seq_len_dict_id;       // v14
    header->sam.segconf_MD_NM_by_un     = segconf.MD_NM_by_unconverted;  // v14
    header->sam.segconf_predict_meth    = segconf.sam_predict_meth_call; // v14
    header->sam.segconf_use_ins_ctxs    = segconf.use_insertion_ctxs;    // 15.0.30
    header->sam.segconf_deep_qname1     = (segconf.deep_qtype == QNAME1);// v15
    header->sam.segconf_deep_qname2     = (segconf.deep_qtype == QNAME2);// v15
    header->sam.segconf_deep_no_qual    = segconf.deep_no_qual;          // v15
    header->sam.segconf_MAPQ_use_xq     = segconf.MAPQ_use_xq;           // 15.0.61
    header->sam.segconf_has_MQ          = !!segconf.has[OPTION_MQ_i];    // 15.0.61
    header->sam.segconf_SA_NM_by_X      = segconf.SA_NM_by_CIGAR_X;      // 15.0.66
    header->sam.segconf_CIGAR_has_eqx   = segconf.CIGAR_has_eqx;         // 15.0.68
    if (segconf.deep_N_fq_score && segconf.deep_N_sam_score)
        header->sam.segconf_deep_N_fq_score = segconf.deep_N_fq_score;   // 15.0.39
    
    unsigned est_sam_factor_mult = round (MAX_(segconf.est_sam_factor, 1) * (double)SAM_FACTOR_MULT);
    if (est_sam_factor_mult > 255) {
        WARN ("est_sam_factor_mult=%u, this might affect genocat translation to SAM. %s", est_sam_factor_mult, report_support_if_unexpected());
        est_sam_factor_mult = 0; // fallback to PIZ default
    }
    header->sam.segconf_sam_factor      = est_sam_factor_mult;

    header->sam.conc_writing_vbs        = BGEN16 (z_has_gencomp ? sam_zip_calculate_max_conc_writing_vbs() : 1); // 15.0.64
}

// initialize SA and OA
static void sam_seg_0X_initialize (VBlockP vb, Did strand_did_i)
{
    // create strand nodes (nodes will be deleted in sam_seg_finalize if not used)
    ctx_create_node (vb, strand_did_i, cSTR("-"));
    ctx_create_node (vb, strand_did_i, cSTR("+"));
}

static void sam_seg_QNAME_initialize (VBlockSAMP vb)
{
    decl_ctx (SAM_QNAME);

    ctx->no_stons = true;              // no singletons, because sam_piz_special_SET_BUDDY uses PEEK_SNIP 
    ctx->flags.store_per_line = true;  // 12.0.41

    qname_seg_initialize (VB, QNAME1, SAM_QNAME); 
    qname_seg_initialize (VB, QNAME2, SAM_QNAME); // we support up to two flavors (eg 2nd flavor can be consensus reads) 

    // initial allocations based on segconf data
    ctx->qname_hash.prm8[0] = MIN_(20, MAX_(14, 32 - __builtin_clz (vb->lines.len32 * 5))); // between 14 and 20 bits - tested - no additional compression benefit beyond 20 bits
    buf_alloc_255(vb, &ctx->qname_hash, 0, (1ULL << ctx->qname_hash.prm8[0]), int32_t, 1, "contexts->qname_hash");

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
    
    // note: all numeric fields need STORE_INT / STORE_FLOAT to be reconstructable to BAM (possibly already set)
    // via the translators set in the SAM_TOP2BAM Container
    ctx_set_store (VB, STORE_INT, SAM_TLEN, SAM_MAPQ, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_GPOS,
                   OPTION_NM_i, OPTION_AS_i, OPTION_MQ_i, OPTION_XS_i, OPTION_XM_i, OPTION_mc_i, OPTION_ms_i, OPTION_Z5_i,
                   OPTION_tx_i, OPTION_YS_i, OPTION_XC_i, OPTION_AM_i, OPTION_SM_i, OPTION_X0_i, OPTION_X1_i, OPTION_CP_i,
                   OPTION_OP_i, OPTION_NH_i, OPTION_HI_i, OPTION_UQ_i, OPTION_cm_i, OPTION_PQ_i,
                   OPTION_XX_i, OPTION_YY_i, OPTION_XO_i, OPTION_a3_i, OPTION_CP_i,
                   OPTION_SA_POS, OPTION_OA_POS,
                   T(MP(BLASR), OPTION_XS_i), T(MP(BLASR), OPTION_XE_i), T(MP(BLASR), OPTION_XQ_i), T(MP(BLASR), OPTION_XL_i), T(MP(BLASR), OPTION_FI_i),
                   T(MP(NGMLR), OPTION_QS_i), T(MP(NGMLR), OPTION_QE_i), T(MP(NGMLR), OPTION_XR_i),
                   T(MP(CPU), OPTION_Y0_i), T(MP(CPU), OPTION_Y1_i),
                   T(MP(TMAP), OPTION_XT_i), T(MP(TMAP), OPTION_XM_i),
                   T(segconf.is_minimap2, OPTION_s1_i), OPTION_ZS_i, OPTION_nM_i,
                   DID_EOL);

    ctx_set_store (VB, STORE_INDEX, SAM_RNAME, SAM_RNEXT, // when reconstructing BAM, we output the word_index instead of the string
                   OPTION_XA_Z, // for containers this stores repeats - used by sam_piz_special_X1->container_peek_repeats
                   OPTION_ML_B_C, // used to predict repeats in MM:Z subcontext MNM:Z
                   DID_EOL);

    // when reconstructing these contexts against another context (DELTA_OTHER or XOR_DIFF) the other may be before or after
    ctx_set_same_line (VB, OPTION_AS_i, OPTION_s1_i, OPTION_PQ_i, // AS may be DELTA_OTHER vs ms:i ; s1 vs AS ; XS vs AS ; PQ vs AS
                       T(segconf.sam_has_BWA_XS_i, OPTION_XS_i),
                       T(MP(CPU), OPTION_Y0_i), T(MP(CPU), OPTION_Y1_i),
                       OPTION_RX_Z, OPTION_CR_Z, // whe UR and CR are reconstruted as XOR_DIFF against UB and CB respectively, we search for the values on the same line (befor or after)
                       DID_EOL);

    // don't store singletons in local even if local is not used for anything else
    ctx_set_no_stons (VB,
                      SAM_RNAME, SAM_RNEXT, // BAM reconstruction needs RNAME, RNEXT word indices. also needed for random access.
                      SAM_POS, SAM_PNEXT, OPTION_mc_i, OPTION_OP_i, OPTION_Z5_i, OPTION_CP_i, // required by seg_pos_field
                      T(IS_PRIM(vb), OPTION_SA_CIGAR), // index needed by sam_load_groups_add_aln_cigar
                      DID_EOL);

    ctx_set_ltype (VB, LT_STRING, OPTION_MD_Z,
                   T(MP(BSBOLT), OPTION_XB_Z),
                   T(MP(BISMARK), OPTION_XM_Z),
                   T(MP(BSSEEKER2), OPTION_XM_Z), T(MP(BSSEEKER2), OPTION_XG_Z),
                   T(MP(NOVOALIGN), OPTION_YH_Z), T(MP(NOVOALIGN), OPTION_YQ_Z),
                   OPTION_CY_Z, OPTION_QT_Z, OPTION_CB_Z, // barcode fallback segging - added to local as string
                   SAM_CIGAR, OPTION_MC_Z, OPTION_OC_Z, OPTION_SA_CIGAR, OPTION_OA_CIGAR, OPTION_XA_CIGAR, OPTION_TX_CIGAR, OPTION_AN_CIGAR, // cigar might be stored as a string if too complex or squanked (see sam_seg_CIGAR, sam_seg_other_CIGAR, sam_cigar_seg_MC_Z, squank_seg)
                   DID_EOL);

    ctx_set_ltype (VB, LT_UINT8, SAM_MAPQ, OPTION_SA_MAPQ, OPTION_SM_i, OPTION_AM_i, OPTION_OA_MAPQ, // MAPQ (and hence other fields carrying mapping quality) is uint8_t by BAM specification
                   OPTION_RX_Z, OPTION_CR_Z, OPTION_GR_Z, OPTION_GY_Z,                               // local has xor_diff data
                   DID_EOL);

    // initialize these to LT_BLOB, the qual-type ones might be changed later to LT_CODEC (eg domq, longr)
    ctx_set_ltype (VB, LT_BLOB, SAM_QUAL, SAM_QUAL_FLANK, OPTION_BD_BI, OPTION_BQ_Z, OPTION_iq_sq_dq,
                   OPTION_QX_Z, OPTION_CY_ARR, OPTION_QT_ARR, OPTION_CR_Z_X, OPTION_BX_Z,
                   OPTION_BI_Z, OPTION_BD_Z,
                   OPTION_CR_Z_X, OPTION_RX_Z_X, OPTION_2R_Z, OPTION_TR_Z, DID_EOL);

    // set ltype=LT_DYN_INT to allow use of seg_integer
    ctx_set_dyn_int (VB, SAM_BUDDY, OPTION_HI_i, OPTION_NM_i, OPTION_NH_i, OPTION_XM_i, OPTION_X1_i, OPTION_CP_i,
                     OPTION_AS_i, OPTION_XS_i, OPTION_ZS_i, OPTION_cm_i, OPTION_ms_i, OPTION_nM_i, OPTION_UQ_i,
                     T(segconf.has_TLEN_non_zero, SAM_TLEN), // note: we don't set if !has_TLEN_non_zero, bc values are stored in b250 and may require singletons
                     T(segconf.sam_has_xcons, OPTION_YY_i),
                     T(segconf.sam_has_xcons, OPTION_XO_i),
                     T(MP(BLASR), OPTION_XS_i), T(MP(BLASR), OPTION_XE_i), T(MP(BLASR), OPTION_XQ_i), T(MP(BLASR), OPTION_XL_i), T(MP(BLASR), OPTION_FI_i),
                     T(MP(NGMLR), OPTION_QS_i), T(MP(NGMLR), OPTION_QE_i), T(MP(NGMLR), OPTION_XR_i),
                     DID_EOL);

    if (segconf.is_collated)
        CTX(SAM_POS)->flags.store_delta = true; // since v12.0.41

    // we may use mates (other than for QNAME) if not is_long_reads (meaning: no mates in this file) and not DEPN components (bc we seg against PRIM)
    if (segconf.is_paired && !IS_DEPN(vb))
        ctx_set_store_per_line (VB, SAM_RNAME, SAM_RNEXT, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_MAPQ,
                                OPTION_MQ_i, OPTION_PQ_i, OPTION_MC_Z, OPTION_SM_i, DID_EOL);

    // case: some lines may be segged against a in-VB saggy line
    if (IS_MAIN(vb)) { // 14.0.0
        ctx_set_store_per_line (VB, SAM_RNAME, SAM_RNEXT, SAM_PNEXT, SAM_POS, SAM_MAPQ, SAM_FLAG,
                                OPTION_SA_Z, OPTION_NM_i, DID_EOL);

        if (vb->check_for_gc)
            buf_alloc (vb, &vb->vb_plan, 0, MAX_(1, vb->lines.len32 / 8), VbPlanItem, 0, "vb_plan"); // inital allocation
    }

    ctx_set_store_per_line (VB, SAM_CIGAR, OPTION_NH_i, T(segconf.is_paired && segconf.sam_multi_RG, OPTION_RG_Z), DID_EOL);
    
    if (segconf.is_paired) {
        for (MatedZFields f = 1; f < NUM_MATED_Z_TAGS; f++) // note: f starts from 1, bc RG is set above with T()
            ctx_set_store_per_line (VB, buddied_Z_dids[f], DID_EOL);

        // exception - we don't use buddying for RX:Z in AGeNT Trimmer (canceling allows RX:Z b250 section elimination)
        if (segconf.has_agent_trimmer)
            CTX(OPTION_RX_Z)->flags.store_per_line = CTX(OPTION_RX_Z)->flags.ctx_specific_flag = false;
    }

    // in --stats, consolidate stats
    ctx_consolidate_stats (VB, SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, 
                           SAM_SEQMIS_A, SAM_SEQMIS_C, SAM_SEQMIS_G, SAM_SEQMIS_T, 
                           SAM_SEQINS_A, SAM_SEQINS_C, SAM_SEQINS_G, SAM_SEQINS_T, DID_EOL);
    ctx_consolidate_stats (VB, SAM_QUAL, SAM_DOMQRUNS, SAM_QUALMPLX, SAM_DIVRQUAL, SAM_QUALSA, SAM_QUAL_PACBIO_DIFF,
                           SAM_QUAL_FLANK, SAM_QUAL_FLANK_DOMQRUNS, SAM_QUAL_FLANK_QUALMPLX, SAM_QUAL_FLANK_DIVRQUAL, 
                           SAM_CQUAL, SAM_CDOMQRUNS, SAM_CQUALMPLX, SAM_CDIVRQUAL, DID_EOL);
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
    ctx_consolidate_stats_(VB, CTX(OPTION_MM_Z), (ContainerP)&segconf.MM_con);

    if (segconf.has[OPTION_HI_i] && !segconf.has[OPTION_SA_Z] && !MP(NOVOALIGN))
        ctx_consolidate_stats (VB, OPTION_HI_i, OPTION_SA_Z, DID_EOL);
    else
        ctx_consolidate_stats (VB, OPTION_SA_Z, OPTION_SA_RNAME, OPTION_SA_POS, OPTION_SA_STRAND, OPTION_SA_CIGAR, OPTION_SA_MAPQ, OPTION_SA_NM, OPTION_SA_MAIN, DID_EOL);

    if (!(segconf.pacbio_subreads && flag.best) && !segconf.multiseq)
        codec_acgt_seg_initialize (VB, SAM_NONREF, true);

    sam_seg_QNAME_initialize (vb);
    sam_seg_QUAL_initialize (vb);
    sam_seg_SEQ_initialize (vb);
    sam_seg_cigar_initialize (vb);
    sam_seg_gc_initialize (vb);
    sam_seg_0X_initialize (VB, OPTION_SA_STRAND);
    sam_seg_0X_initialize (VB, OPTION_OA_STRAND);

    if (MP(ULTIMA) || segconf.running) // note: need also in segconf, so we can identify Ultima parameters in case it is Ultima (no harm if it is not)
        sam_ultima_seg_initialize (vb);
    else
        CTX(OPTION_t0_Z)->no_callback = true; // override Ultima's sam_zip_t0
    
    if (segconf.has_10xGen)        sam_10xGen_seg_initialize (vb);
    if (segconf.has_agent_trimmer) agilent_seg_initialize    (VB);
    if (segconf.sam_has_abra2)     sam_abra2_seg_initialize  (vb);
    if (MP(DRAGEN))                sam_dragen_seg_initialize (vb);
    if (MP(STAR))                  sam_star_seg_initialize   (vb);
    if (TECH(PACBIO))              sam_pacbio_seg_initialize (vb);

    if (segconf.sam_has_BWA_XS_i) // XS:i is as defined some aligners
        seg_mux_init (vb, OPTION_XS_i, SAM_SPECIAL_BWA_XS, false, XS);

    else if (MP(HISAT2)) // ZS:i is like BWA's XS:i
        seg_mux_init (vb, OPTION_ZS_i, SAM_SPECIAL_BWA_XS, false, XS);

    if (segconf.sam_has_XM_i_is_mismatches)
        CTX(OPTION_XM_i)->flags.same_line = true;   // copy NM from same line

    if (segconf.sam_has_bowtie2_YS_i) {
        ctx_set_store_per_line (VB, OPTION_AS_i, OPTION_YS_i, DID_EOL);
        seg_mux_init (vb, OPTION_YS_i, SAM_SPECIAL_DEMUX_BY_MATE, false, YS);
    }

    if (MP(STAR) && segconf.is_paired && !IS_DEPN(vb)) {
        seg_mux_init (vb, OPTION_nM_i, SAM_SPECIAL_DEMUX_BY_MATE, false, nM);
        ctx_set_store_per_line (VB, OPTION_AS_i, OPTION_nM_i, DID_EOL);
    }

    if (segconf.sam_is_nanoseq) {
        ctx_set_store_per_line (VB, OPTION_rb_Z, OPTION_mb_Z, DID_EOL);
        seg_mux_init (vb, OPTION_rb_Z, SAM_SPECIAL_DEMUX_BY_MATE, false, rb);
        seg_mux_init (vb, OPTION_mb_Z, SAM_SPECIAL_DEMUX_BY_MATE, false, mb);
    }

    if (segconf.sam_ms_type == ms_BIOBAMBAM) {
        seg_mux_init (vb, OPTION_ms_i, SAM_SPECIAL_DEMUX_BY_MATE, false, ms);

        CTX(OPTION_ms_i)->flags.spl_custom = true;  // custom store-per-line - SPECIAL will handle the storing
        CTX(OPTION_ms_i)->flags.store = STORE_INT;  // since v14 - store QUAL_score for mate ms:i (in v13 it was stored in QUAL)

    }
    if (segconf.sam_has_xcons) {
        seg_mux_init (vb, OPTION_YY_i, SAM_SPECIAL_DEMUX_BY_XX_0, false, YY);
        seg_mux_init (vb, OPTION_XO_i, SAM_SPECIAL_DEMUX_BY_AS,   false, XO);
    }

    seg_mux_init (vb, SAM_FLAG,    SAM_SPECIAL_DEMUX_BY_BUDDY,     false, FLAG);
    seg_mux_init (vb, SAM_POS,     SAM_SPECIAL_DEMUX_BY_MATE_PRIM, true,  POS);
    seg_mux_init (vb, SAM_PNEXT,   SAM_SPECIAL_PNEXT,              true,  PNEXT);
    seg_mux_init (vb, SAM_MAPQ,    SAM_SPECIAL_DEMUX_MAPQ,         false, MAPQ);
    seg_mux_init (vb, OPTION_MQ_i, SAM_SPECIAL_DEMUX_BY_MATE,      false, MQ);
    seg_mux_init (vb, OPTION_MC_Z, SAM_SPECIAL_DEMUX_BY_MATE,      false, MC);
    seg_mux_init (vb, OPTION_AS_i, SAM_SPECIAL_DEMUX_BY_MATE,      false, AS);
    seg_mux_init (vb, OPTION_PQ_i, SAM_SPECIAL_DEMUX_BY_MATE,      false, PQ);
    seg_mux_init (vb, OPTION_NH_i, SAM_SPECIAL_DEMUX_BY_BUDDY_MAP, false, NH);

    for (MatedZFields f = 0; f < NUM_MATED_Z_TAGS; f++)
        seg_mux_init (vb, buddied_Z_dids[f], SAM_SPECIAL_DEMUX_BY_BUDDY, false, mated_z_fields[f]);
        
    // get counts of qnames in the VB, that will allow us to leave lines in the VB for saggy instead of gencomp
    // note: for SAG_BY_SA in --best, we send all prim/depn to gencomp (see bug 629)
    if (vb->check_for_gc && !flag.force_gencomp && (/*IS_SAG_SA ||*/ IS_SAG_NH || IS_SAG_SOLO)) // SA still doesn't work well with saggy - not even QUAL
        scan_index_qnames_seg (vb);

    // note: SA_HtoS might be already set if we segged an line containing SA:Z against a prim saggy
    if (IS_DEPN(vb) && IS_SAG_SA && segconf.SA_HtoS == unknown)
        segconf.SA_HtoS = segconf.depn_CIGAR_can_have_H && // set when segging MAIN component
                          !segconf.SA_CIGAR_can_have_H;    // set when segging PRIM component

    COPY_TIMER(seg_initialize);
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
}

static uint32_t num_lines_at_max_len (VBlockSAMP vb)
{
    uint32_t count = 0;
    for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++)
        if (DATA_LINE(line_i)->SEQ.len == vb->longest_seq_len)
            count++;

    return count;
}

// called after reading SAM header, and sometimes also after segconf
void sam_segconf_set_by_MP (void)
{
    if (MP(bwa)) segconf.sam_mapper = MP_BWA; // we consider "bwa" to be "BWA"

    // a small subset of biobambam2 programs - only those that generate ms:i / mc:i tags
    segconf.is_biobambam2_sort = stats_is_in_programs ("bamsormadup") || stats_is_in_programs ("bamsort") || stats_is_in_programs ("bamtagconversion");

    segconf.has_bwa_meth = MP(BWA) && stats_is_in_programs ("bwa-meth"); // https://github.com/brentp/bwa-meth

    segconf.has_bqsr = stats_is_in_programs ("ApplyBQSR");

    segconf.has_RSEM = stats_is_in_programs ("RSEM") || stats_is_in_programs ("rsem");

    // note: this file *might* be of bisulfite-treated reads. This variable might be reset in sam_segconf_finalize if it fails additonal conditions 
    segconf.sam_bisulfite     = MP(BISMARK) || MP(BSSEEKER2) || MP(DRAGEN) || MP(BSBOLT) || MP(GEM3) || MP(ULTIMA) || segconf.has_bwa_meth;
    segconf.sam_has_bismark_XM_XG_XR = MP(BISMARK) || MP(DRAGEN) || MP(ULTIMA);

    // in bisulfate data, we still calculate MD:Z and NM:i vs unconverted reference
    segconf.MD_NM_by_unconverted = MP(BISMARK) || MP(DRAGEN) || MP(BSBOLT) || MP(GEM3) || MP(BSSEEKER2);

    segconf.is_bwa            = MP(BWA) || MP(bwa) || MP(BSBOLT) || MP(CPU) || MP(BWA_MEM2) || MP(PARABRICKS); // aligners based on bwa
    segconf.is_minimap2       = MP(MINIMAP2) || MP(WINNOWMAP) || MP(PBMM2);   // aligners based on minimap2
    segconf.is_bowtie2        = MP(BOWTIE2) || MP(HISAT2) || MP(TOPHAT) || MP(BISMARK) || MP(BSSEEKER2); // aligners based on bowtie2

    segconf.sam_has_SA_Z      = segconf.is_bwa || segconf.is_minimap2 || MP(NGMLR) || MP(DRAGEN) || MP(NOVOALIGN) || MP(ULTIMA) || MP(ISAAC) || MP(CRDNA); /*|| MP(LONGRANGER); non-standard SA:Z format (POS is off by 1, main-field NM is missing) */ 
    
    segconf.sam_has_BWA_XA_Z  = (segconf.is_bwa || MP(GEM3) || MP(GEM2SAM) || MP(DELVE) || MP(DRAGEN) || MP(ULTIMA) || MP(CRDNA)) ? yes 
                              : MP(TMAP) || MP(TORRENT_BC)                                                                        ? no 
                              :                                                                                                     segconf.sam_has_BWA_XA_Z; // remain "unknown" or as determined by segconf
                              
    segconf.sam_has_BWA_XS_i  = segconf.is_bwa || MP(TMAP) || MP(GEM3) || (segconf.is_bowtie2 && !MP(HISAT2)) || MP(CPU) || MP(LONGRANGER) || MP(DRAGEN) || MP(CRDNA);
    segconf.sam_has_XM_i_is_mismatches  = segconf.is_bwa || segconf.is_bowtie2 || MP(NOVOALIGN) || MP(DRAGEN);
    segconf.sam_has_BWA_XT_A  = segconf.is_bwa || MP(DRAGEN);
    segconf.sam_has_BWA_XC_i  = segconf.is_bwa || MP(DRAGEN); // might be canceled during segconf (see sam_seg_BWA_XC_i)
    segconf.sam_has_BWA_X01_i = segconf.is_bwa || MP(DRAGEN);

    segconf.sam_has_bowtie2_YS_i = MP(BOWTIE2) || MP(BSSEEKER2) || MP(HISAT2);

    if (!segconf.is_minimap2) // we assume SA:Z CIGAR is abbreviated IFF is_minimap2 (this might be refined in the future) 
        segconf.SA_CIGAR_abbreviated = no; 
}

// finalize Seg configuration parameters
void sam_segconf_finalize (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    
    segconf.longest_seq_len = vb->longest_seq_len;
    segconf.is_long_reads   = segconf_is_long_reads();
    segconf.sam_cigar_len   = 1 + ((segconf.sam_cigar_len - 1) / vb->lines.len32);                   // set to the average CIGAR len (rounded up)
    segconf.est_sam_factor  = (double)segconf.est_segconf_sam_size / (double)Ltxt;

    if (num_lines_at_max_len(vb) > vb->lines.len32 / 2      &&  // more than half the lines are at exactly maximal length
        (vb->lines.len32 > 100 || txt_file->no_more_blocks) &&  // enough lines to be reasonably convinced that this is not by chance
        !segconf.is_long_reads)        // TO DO: trimming long-read qual in FASTQ with --deep would mess up LONGR codec, we need to sort this out
        
        segconf.sam_cropped_at = vb->longest_seq_len; // possibily the FASTQ reads were cropped to be all equal length

    segconf.sam_seq_len = (uint32_t)(0.5 + (double)segconf.sam_seq_len / (double)vb->lines.len32); // set average seq_len - rounded to the nearest

    segconf.seq_len_to_cm /= vb->lines.len32;
    if (segconf.seq_len_to_cm > 255)
        segconf.seq_len_to_cm = 0;

    // if unmapped test if this might is multiseq
    if (segconf.sam_is_unmapped && !flag.fast) 
        segconf_test_multiseq (VB, SAM_NONREF);

    // 1. earlier Ultima files where bwa was used, 2. headerless files - if we have ultima fields treat as MP_ULTIMA
    if (TECH(ULTIMA) && (MP(BWA) || MP(UNKNOWN)) && segconf.has[OPTION_t0_Z] && segconf.has[OPTION_tp_B_c]) {
        segconf.sam_mapper = MP_ULTIMA;
        sam_segconf_set_by_MP();
    }

    // in some STAR transcriptome files, there are no @PG header lines. it is important that mapper is set to STAR so we can use liberal saggy_line_i, as SA groups are large
    if (MP(UNKNOWN) && segconf.has[OPTION_NH_i] && segconf.has[OPTION_HI_i] && !segconf.has[OPTION_SA_Z]) {
        segconf.sam_mapper = MP_STAR;
        sam_segconf_set_by_MP();
    }

    // SA:Zs analyzed sam_segconf_SA_cigar_cb found no evidence of no-abbrivation: since 
    if (segconf.has[OPTION_SA_Z] &&              // we analyzed some SA_CIGARs in sam_segconf_SA_cigar_cb
        segconf.SA_CIGAR_abbreviated == unknown) // unknown means that sam_segconf_set_by_MP determined that SA_CIGARs might be abbreivated and sam_segconf_SA_cigar_cb found no evidence that they are not
        segconf.SA_CIGAR_abbreviated = yes;

    // sequencing technologies that results in lots of insertion errors - will benefit from the use_insertion_ctxs method
    if (TECH(PACBIO) || TECH(NANOPORE) || TECH(ULTIMA))
        segconf.use_insertion_ctxs = true;

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
        segconf.star_solo  = true; // must be before sam_set_sag_type
        segconf.has_10xGen = true; // STAR Solo mimics cellranger's tags
    }

    // set optimizations
    if (flag.optimize) 
        sam_segconf_finalize_optimizations();

    segconf.deep_paired_qname = flag.deep && segconf.is_paired && flag.deep_num_fastqs == 2;
    
    segconf.sam_has_zm_by_Q1NAME = TECH(PACBIO) && segconf_qf_id (QNAME1) != QF_PACBIO_3;

    segconf.pacbio_subreads = TECH(PACBIO) && segconf.sam_is_unmapped && segconf.has[OPTION_pw_B_C] && segconf.has[OPTION_ip_B_C];

    if (segconf.pacbio_subreads)
        ctx_segconf_set_hard_coded_lcodec (OPTION_ip_ARR, CODEC_ARTB); // as good as LZMA on ip:B:C and hugely faster - hard-code to prevent assigning LZMA  

    segconf.sam_use_sn_mux = segconf.pacbio_subreads && segconf.sam_is_unmapped && segconf_qf_id (QNAME1) == QF_PACBIO_rng;

    segconf.MAPQ_use_xq = MP(DRAGEN) && segconf.has[OPTION_xq_i];

    segconf.AS_is_ref_consumed  = (segconf.AS_is_ref_consumed  > vb->lines.len32 / 2);
    segconf.AS_is_2ref_consumed = (segconf.AS_is_2ref_consumed > vb->lines.len32 / 2); // AS tends to be near 2 X ref_consumed, if at least half of the lines say so

    // possibly reset sam_bisulfite, first set in sam_header_zip_inspect_PG_lines->sam_segconf_set_by_MP
    segconf.sam_bisulfite = MP(BISMARK) || MP(BSSEEKER2) || MP(BSBOLT) ||
                           (MP(DRAGEN)   && segconf.has[OPTION_XM_Z])  ||
                           (MP(ULTIMA)   && segconf.has[OPTION_XM_Z])  ||
                           (MP(GEM3)     && segconf.has[OPTION_XB_A])  ||
                           segconf.has_bwa_meth; // note: methylation tags MM:Z and ML:B:C do not imply bisulfite - e.g. PacBio detects methylation differently: https://www.pacb.com/products-and-services/applications/epigenetics/

    segconf.sam_has_bismark_XM_XG_XR &= segconf.sam_bisulfite;
    
    ASSINP (!segconf.sam_bisulfite        // not a bisulfite file
         || segconf.sam_is_unmapped       // we can't use the bisulfite methods if unmapped             
         || flag.force                    // --force overrides
         || IS_REF_EXTERNAL || IS_REF_EXT_STORE   // reference is provided
         || flag.zip_no_z_file,           // we're not creating a compressed format
            "Compressing bisulfite file %s requires using --reference. Override with --force.", dt_name (vb->data_type));

    segconf.sam_predict_meth_call = segconf.sam_bisulfite &&
                                    !IS_REF_INTERNAL      && // bug 648
                                    (MP(BISMARK) || MP(DRAGEN) || MP(BSBOLT)); // have methylation call tags


    if (segconf.flav_prop[QNAME2].is_consensus && segconf.has[OPTION_XX_i] && segconf.has[OPTION_YY_i] && segconf.has[OPTION_XC_i] && segconf.has[OPTION_XO_i]) {
        segconf.sam_has_xcons = true;
        segconf.sam_has_BWA_XC_i = false;
    }
    
    if (segconf.has[OPTION_rb_Z] && segconf.has[OPTION_rb_Z] == segconf.has[OPTION_mb_Z])
        segconf.sam_is_nanoseq = true;

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
        
        // we have at least one pair of lines with the same QNAME, and the file is not sorted
        if (!segconf.is_sorted && segconf.evidence_of_collated)
            segconf.is_collated = true;
    }

    segconf.disable_random_acccess = !segconf.is_sorted;

    segconf.has_barcodes = segconf.has[OPTION_CB_Z] || segconf.has[OPTION_CR_Z] || segconf.has[OPTION_CY_Z] || segconf.has[OPTION_BX_Z] || 
                           segconf.has[OPTION_RX_Z] || segconf.has[OPTION_QX_Z] || segconf.has[OPTION_BC_Z] || segconf.has[OPTION_QT_Z];

    // second is_paired test: PNEXT is not 0 for all lines (first test is in sam_seg_FLAG)
    if (segconf.has[SAM_PNEXT])
        segconf.is_paired = true;

    segconf.has_MD_or_NM = segconf.has[OPTION_MD_Z] || segconf.has[OPTION_NM_i] ||
                           (MP(STAR) && segconf.has[OPTION_nM_i] && !segconf.is_paired); // we use the NM method to seg nM in this case

    // use the number of CIGAR's X bases to determine SA_NM. Note that this excludes I and D included in NM:i, so 
    // if NM:i is present as well, we expect it to be different that SA_NM.
    segconf.SA_NM_by_CIGAR_X = segconf.CIGAR_has_eqx; // TO DO: test explicitly if SA_NM is determined by CIGAR's X bases or by NM:i (bug 1139)  

    // note: in case of both is_biobambam2_sort and minimap2, it is likely biobambam's tag
    if (segconf.has[OPTION_ms_i]) {
        if (segconf.is_biobambam2_sort && segconf.is_paired) 
            segconf.sam_ms_type = ms_BIOBAMBAM;

        else if (segconf.is_minimap2)
            segconf.sam_ms_type = ms_MINIMAP2; 

        // if we have no conclusive evidence of either minimap2 or biobambam, this is likely samtools, which has the same logic as biobambam
        else if (segconf.is_sorted)
            segconf.sam_ms_type = ms_BIOBAMBAM;
    }

    // update context tag names if this file has UB/UR/UY which are aliased to BX/RX/QX
    if (segconf.has[OPTION_UB_Z] || segconf.has[OPTION_UR_Z]) 
        sam_segconf_retag_UBURUY();

    // note: same logic as in fastq_segconf_is_saux
    if (segconf.has[OPTION_ZA_Z] && segconf.has[OPTION_ZB_Z] && segconf.has[OPTION_RX_Z] && segconf.has[OPTION_QX_Z] && segconf.has[OPTION_BC_Z]) {
        segconf.has_agent_trimmer = true;
        stats_add_one_program (_S("AGeNT_Trimmer")); // note: Trimmers adds the fields to the FASTQ file (as SAUX), and copied to BAM with eg bwa mem -C, so no @PG line
    }

    // cases where aligner is available (note: called even if reference is not loaded, so that it errors in segconf_calculate)
    if (  flag.best 
       || flag.deep 
       || ( IS_REF_LOADED_ZIP 
         && !segconf.is_long_reads 
         && (  ( !flag.low_memory && (segconf.num_mapped < vb->lines.len32 || bam_txt_file_is_last_alignment_unmapped())) 
            || ( flag.low_memory && segconf.num_mapped == 0)))) { 

        flag.aligner_available = true;
        if (IS_REF_LOADED_ZIP) refhash_load_standalone();
    }

    // cases where aligner is not available, despite setting in main_load_reference
    else if (segconf.is_long_reads) 
        flag.aligner_available = false;

    // with REF_EXTERNAL and unaligned data, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (z_file->num_txts_so_far == 1 && (flag.aligner_available || !sam_hdr_contigs) && IS_REF_LOADED_ZIP)
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, ref_get_ctgs (gref));

    if (TECH(PACBIO) && segconf.has[OPTION_dq_Z] == vb->lines.len32 && segconf.has[OPTION_iq_Z] == vb->lines.len32 && segconf.has[OPTION_sq_Z] == vb->lines.len32
        && flag.force_qual_codec != CODEC_PACB)
        segconf.use_pacbio_iqsqdq = true;

    // Aplogize to user for not doing a good job with PacBio subreads files
    if (segconf.pacbio_subreads) {
        TEMP_FLAG(quiet, flag.explicit_quiet); // note: quiet is set to true in segconf
        WARN0 ("FYI: Genozip currently doesn't do a very good job at compressing PacBio subreads files. This is because we haven't figured out yet a good method to compress PacBio kinetic data - the ip:B and pw:B fields. Sorry!\n");
        RESTORE_FLAG(quiet);
    }

    if (MP(ULTIMA))
        sam_ultima_finalize_segconf (vb);

    if (flag.reference == REF_INTERNAL && !txt_file->redirected && !flag.seg_only && (segconf.sam_is_unmapped && !segconf.is_long_reads))
        TIP ("Compressing this unmapped %s file using a reference would shrink it by an additional 20%%-50%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n",
             dt_name (txt_file->data_type), arch_get_argv0(), txt_file->name);
    
    ASSERT (!flag.debug_or_test || !segconf.has[OPTION_SA_Z] || segconf.sam_has_SA_Z || MP(LONGRANGER) || MP(UNKNOWN),
            "%s produces SA:Z, expecting segconf.sam_has_SA_Z to be set", segconf_sam_mapper_name()); // should be set in sam_header_zip_inspect_PG_lines

    if (codec_pacb_maybe_used (SAM_QUAL)) 
        codec_pacb_segconf_finalize (VB);
    
    if (codec_longr_maybe_used (VB, SAM_QUAL))
        codec_longr_segconf_calculate_bins (VB, CTX(SAM_QUAL + 1), sam_zip_qual);

    // RG method - QNAME - if successful, we expect all RGs to have been added to dict in sam_header_zip_inspect_RG_lines, and on the SPECIAL segged in this VB 
    if (segconf.sam_multi_RG && CTX(OPTION_RG_Z)->nodes.len32 == 1 && *Bc(CTX(OPTION_RG_Z)->dict, 1) == SAM_SPECIAL_RG_by_QNAME)
        segconf.RG_method = RG_CELLRANGER; // set if all segconf lines agree

    // note: we calculate the smux stdv to be reported in stats, even if SMUX is not used
    codec_smux_calc_stats (VB);

    qname_segconf_finalize (VB);
}

// main thread: this is SAM/BAM's zip_after_segconf callback
void sam_zip_after_segconf (VBlockP vb)
{
    if (segconf.line_len) // not header-only file
        sam_set_sag_type();   

    if (flag.deep)
        buf_alloc_zero (evb, &z_file->vb_num_deep_lines, 0, 10000, uint32_t, 0, "z_file->vb_num_deep_lines");

    if (segconf.sag_type || flag.deep) {
        huffman_produce_compressor (SAM_QNAME, (bool[256]){[32 ... 126] = true}); // converts qname_huff from HuffmanChewer to HuffmanCodes

        // check if huffman compressed QUAL better than ARITH (for in-memory compressing of QUAL in gencomp ZIP/PIZ and deep PIZ)
        if (VB_SAM->has_qual) 
            sam_qual_produce_huffman_if_better (VB_SAM);

        if (IS_SAG_SOLO)
            sam_produce_solo_huffmans (VB_SAM);
    }

    if (IS_SAG_SA || flag.deep)
        huffman_produce_compressor (SAM_CIGAR, (bool[256]){['0'...'9']=1, ['M']=1, ['I']=1, ['D']=1, ['N']=1, ['S']=1, ['H']=1, ['P']=1, ['=']=1, ['X']=1, ['*']=1 });
    else
        buf_destroy (CTX(SAM_CIGAR)->huffman); // not needed
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
        codec_assign_best_qual_codec (VB, SAM_QUAL, sam_zip_qual, false, true, &vb->codec_requires_seq);

    if (vb->has_qual && CTX(SAM_CQUAL)->local.len32) 
        codec_assign_best_qual_codec (VB, SAM_CQUAL, sam_zip_cqual, false, true, &vb->codec_requires_seq);

    if (CTX(OPTION_OQ_Z)->local.len32 && !codec_oq_comp_init (VB)) 
        codec_assign_best_qual_codec (VB, OPTION_OQ_Z, sam_zip_OQ, false, true, &vb->codec_requires_seq);
    
    if (CTX(OPTION_TQ_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_TQ_Z, sam_zip_TQ, true, false, &vb->codec_requires_seq);

    if (CTX(OPTION_CY_ARR)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_CY_ARR, NULL, true, false, &vb->codec_requires_seq);

    if (CTX(OPTION_QX_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_QX_Z, segconf.has_agent_trimmer ? NULL : sam_zip_QX, true, false, &vb->codec_requires_seq);

    if (CTX(OPTION_2Y_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_2Y_Z, sam_zip_2Y, true, true, &vb->codec_requires_seq);

    if (CTX(OPTION_QT_ARR)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_QT_Z, NULL, true, false, &vb->codec_requires_seq);

    if (CTX(OPTION_U2_Z)->local.len32)
        codec_assign_best_qual_codec (VB, OPTION_U2_Z, sam_zip_U2, true, true, &vb->codec_requires_seq);

    // determine if sam_piz_sam2bam_SEQ ought to store vb->textual_seq (saves piz time if not)
    CTX(SAM_SQBITMAP)->flags.no_textual_seq = !vb->codec_requires_seq            && // a QUAL-like field uses a codec that requires SEQ
                                              segconf.sam_mapper != MP_BSSEEKER2 &&
                                              !segconf.has[OPTION_tp_B_c];

    if (flag.no_biopsy_line) 
        sam_seg_toplevel (VB);

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

    // VB1: if we've not found depn lines in the VB, no need to check for gencomp lines in future VBs 
    // note that it is possible that some VBs had gencomp lines and already absorbed
    if (vb->vblock_i == 1 && !vb->seg_found_depn_line && !flag.force_gencomp)
        segconf.no_gc_checking = true;

    // note: absorb even if no gc checking, since we need to increment counter
    if (segconf.sag_type && IS_MAIN(vb))
        gencomp_absorb_vb_gencomp_lines (VB);
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

    // case: no PRIM or DEPN lines - this file will have just one (main) component (but don't change if BIND_DEEP)
    if (flag.bind == BIND_SAM && !gencomp_have_any_lines_absorbed())
        flag.bind = BIND_NONE; 
}

bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return
        // typically small
        dict_id.num == _SAM_TOPLEVEL     ||
        dict_id.num == _SAM_TOP2BAM      ||
        dict_id.num == _SAM_FLAG         ||
        dict_id.num == _SAM_MAPQ         ||
        dict_id.num == _SAM_QNAME        ||
        dict_id.num == _SAM_RNAME        ||
        dict_id.num == _SAM_RNEXT        ||
        dict_id.num == _SAM_TLEN         ||
        dict_id.num == _SAM_AUX          ||
        dict_id.num == _SAM_EOL          ||
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
        dict_id.num == _OPTION_XA_RNAME  || // not _OPTION_XA_Z!  it can be large if non-bwa format (which can happen even if using bwa)
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

    ASSSEG ((aux_len >= 6 || (aux_len == 5 && (aux[3]=='Z' || aux[3]=='H'))) && aux[2] == ':' && aux[4] == ':', "invalid optional field format: %.*s", STRf (aux));

    *tag = aux;
    *type = aux[3];

    if (*type == 'B') {
        *array_subtype = aux[5];
        *value = (aux_len >= 7) ? (aux + 7) : (aux + 6);
        *value_len = (aux_len >= 7) ? (aux_len - 7) : 0;
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
        ASSSEG (vb->aux_lens[f] > (is_bam ? 3 : 5) || 
                (!is_bam && vb->aux_lens[f] == 5 && (vb->auxs[f][3] == 'Z' || vb->auxs[f][3] == 'H')), // in SAM, Z, H can be data-less, in BAM, they are nul-terminated
                "Invalid auxilliary field format. AUX tag #%u (0-based) has a length of %u, expecting %u: \"%.*s\"",
                f, vb->aux_lens[f], (is_bam ? 3 : 5), MIN_(16, vb->aux_lens[f]), vb->auxs[f]);

        char c1 = vb->auxs[f][0];
        char c2 = vb->auxs[f][1];
        char c3 = is_bam ? sam_seg_bam_type_to_sam_type (vb->auxs[f][2]) : vb->auxs[f][3];

        #define AUXval(c1,c2,c3) (((uint32_t)(c1)) << 16 | ((uint32_t)(c2)) << 8 | ((uint32_t)(c3)))

        #define TEST_AUX(name, c1, c2, c3) \
                case AUXval(c1,c2,c3): vb->idx_##name = f; break;

        #define TEST_AUX_B(name, c1, c2, c3, array_subtype) \
                case AUXval(c1,c2,c3): if (vb->auxs[f][is_bam ? 3 : 5] == (array_subtype)) vb->idx_##name = f; break;

        switch (AUXval(c1, c2, c3)) {
            TEST_AUX(NM_i, 'N', 'M', 'i');
            TEST_AUX(NH_i, 'N', 'H', 'i');
            TEST_AUX(MD_Z, 'M', 'D', 'Z');
            TEST_AUX(SA_Z, 'S', 'A', 'Z');
            TEST_AUX(HI_i, 'H', 'I', 'i');
            TEST_AUX(IH_i, 'I', 'H', 'i');
            TEST_AUX(XG_Z, 'X', 'G', 'Z');
            TEST_AUX(XM_Z, 'X', 'M', 'Z');
            TEST_AUX(XM_i, 'X', 'M', 'i');
            TEST_AUX(X0_i, 'X', '0', 'i');
            TEST_AUX(X1_i, 'X', '1', 'i');
            TEST_AUX(XA_Z, 'X', 'A', 'Z');
            TEST_AUX(AS_i, 'A', 'S', 'i');
            TEST_AUX(XS_i, 'X', 'S', 'i');
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
            TEST_AUX(ZA_Z, 'Z', 'A', 'Z');
            TEST_AUX(ZB_Z, 'Z', 'B', 'Z');
            TEST_AUX(pr_i, 'p', 'r', 'i');
            TEST_AUX(qs_i, 'q', 's', 'i');
            TEST_AUX(ws_i, 'w', 's', 'i');
            TEST_AUX(xq_i, 'x', 'q', 'i');

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

    if (idx == -1) {
        if (soft_fail) return 0;
        ABORT("%s: idx==-1, perhaps field is not in this line?", LN_NAME);
    } 

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
    ABORT("%s: value of %.2s=%" PRId64 " ∉ [%d,%d]", LN_NAME, vb->auxs[idx], numeric.i, min_value, max_value);
}

// returns the length of the field, or 0 if it is not found
uint32_t sam_seg_get_aux_float (VBlockSAMP vb, int16_t idx,
                                double *number, // modified only if float is parsed
                                bool is_bam, FailType soft_fail)
{
    if (idx == -1) {
        if (soft_fail) return 0;
        ABORT("%s: idx==-1, perhaps field is not in this line?", LN_NAME);
    } 

    ValueType numeric;
    uint32_t value_len;
    char *after;

    rom value = sam_seg_get_aux (vb, idx, &value_len, &numeric, is_bam);

    if (value && 
        (  is_bam || // note: use strtod as in reconstruct_one_snip to avoid discrepencies in float parsing
           ((numeric.f = strtod (value, &after)) && after == value + value_len))) {
        *number = numeric.f;
        return value_len;
    }

    // this line doesn't have this field, or (SAM only) the field is not a valid float
    else if (soft_fail) 
        return 0;
    
    else
        ABORT("%s: no valid float found for %.2s", LN_NAME, vb->auxs[idx]);
}

void sam_seg_get_aux_Z(VBlockSAMP vb, int16_t idx, pSTRp (snip), bool is_bam)
{
    ASSSEG (idx >= 0, "%s: idx==-1, perhaps field is not in this line?", LN_NAME);

    *snip = vb->auxs[idx] + (is_bam ? 3 : 5);
    *snip_len = vb->aux_lens[idx] - (is_bam ? 4 : 5); // bam: remove the count of the \0 that is including in aux_lens[idx] in BAM
}

char sam_seg_get_aux_A (VBlockSAMP vb, int16_t idx, bool is_bam)
{
    ASSSEG (idx >= 0, "%s: idx==-1, perhaps field is not in this line?", LN_NAME);

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
            if (IS_REF_LOADED_ZIP) 
                WARN_ONCE("FYI: RNAME \"%.*s\" (and possibly others) is missing in the reference file. This might impact the compression ratio.",
                          vb->chrom_name_len, vb->chrom_name);
            return; // the sequence will be segged as unaligned
        }
    }
}

void sam_seg_aux_all (VBlockSAMP vb, ZipDataLineSAMP dl)
{
    START_TIMER;

    const bool is_bam = IS_BAM_ZIP;
    Container con = { .repeats = 1, .filter_repeats = true /* v14 */};
    char prefixes[vb->n_auxs * 6 + 3];                // each name is 5 characters per SAM specification, eg "MC:Z:" followed by CON_PX_SEP ; +3 for the initial CON_PX_SEP. 
    prefixes[0] = prefixes[1] = prefixes[2] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix followed by separator for empty prefix for translator-only item[0]
    unsigned prefixes_len = 3;

    // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM
    con.items[con_nitems (con)] = (ContainerItem){ .translator = SAM2BAM_AUX_SELF };
    con_inc_nitems (con);

    for (int16_t idx = 0; idx < vb->n_auxs; idx++) {

        STR0(value);
        ValueType numeric = {};
        rom tag;
        char sam_type = 0, bam_type = 0, array_subtype = 0;

        if (is_bam)
            bam_get_one_aux (vb, idx, &tag, &bam_type, &array_subtype, pSTRa (value), &numeric);
        else
            sam_get_one_aux (vb, idx, &tag, &sam_type, &array_subtype, pSTRa (value));

        if (sam_type == 'i') {
            if (idx == vb->idx_NM_i) // we already converted NM to integer, no need to do it again
                numeric.i = dl->NM;
            else
                ASSSEG (str_get_int (STRa (value), &numeric.i), "%s: Expecting integer value for auxiliary field %c%c but found \"%.*s\"",
                        LN_NAME, tag[0], tag[1], STRf (value));
            value = 0;
        }

        if (!bam_type) bam_type = sam_seg_sam_type_to_bam_type (sam_type, numeric.i);
        if (!sam_type) sam_type = sam_seg_bam_type_to_sam_type (bam_type);

        if (sam_type == 'Z')
            for (int i=0; i < value_len; i++)
                ASSSEG (value[i] >= ' ' && value[i] <= '~', // SAM specification section 1.5
                        "Invalid character in %c%c:Z field [index=%u]: ASCII %u", tag[0], tag[1], i, (uint8_t)value[i]);

        else if (sam_type == 'H')
            ASSSEG (str_is_hexup (STRa(value)), "Invalid character in %c%c:H field", tag[0], tag[1]);

        ASSERT (bam_type, "%s: value %" PRId64 " of field %c%c is out of range of the BAM specification: [%d-%u]",
                LN_NAME, numeric.i, tag[0], tag[1], -0x80000000, 0x7fffffff);

        con.items[con_nitems (con)] = (ContainerItem){
            .dict_id    = sam_seg_aux_field (vb, dl, is_bam, tag, bam_type, array_subtype, STRa(value), numeric, idx),
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
        // NULL means MISSING Container item (of the toplevel container) - will cause container_reconstruct of
        // the toplevel container to delete of previous separator (\t)
        container_seg (vb, CTX(SAM_AUX), 0, 0, 0, 0);

    COPY_TIMER(sam_seg_aux_all);
}   

static inline bool has_same_qname (VBlockSAMP vb, STRp (qname), LineIType buddy_line_i, bool last_char_flipped)
{
    if (buddy_line_i == NO_LINE) return false;

    TxtWord buddy_q = DATA_LINE(buddy_line_i)->QNAME;

    if (!last_char_flipped) 
        return str_issame_(STRa(qname), STRtxt (buddy_q));

    else {   
        char r_a = qname[qname_len-1];
        char r_b = *Btxt(buddy_q.index + buddy_q.len - 1);

        return str_issame_(STRa(qname) - 1, STRtxt (buddy_q) - 1) &&
               ((r_a=='1' && r_b=='2') || (r_a=='2' && r_b=='1'));
    }
}

#define LINE_BY_HASH(hash) *B(LineIType, CTX(SAM_QNAME)->qname_hash, (hash))

// seg mate as buddy and return true if this line has one
static inline BuddyType sam_seg_mate (VBlockSAMP vb, SamFlags f, STRp (qname), uint32_t my_hash, bool *insert_to_hash)
{
    if (sam_is_depn (f) || !segconf.is_paired) return BUDDY_NONE;

    uint32_t mate_hash = qname_calc_hash (QNAME1, COMP_NONE, STRa(qname), !f.is_last, true, CRC32, NULL) & MAXB(CTX(SAM_QNAME)->qname_hash.prm8[0]);
    LineIType candidate = LINE_BY_HASH(mate_hash);
    SamFlags *mate_f = &DATA_LINE(candidate)->FLAG; // invalid pointer if no mate

    // case: mate is found
    if (has_same_qname (vb, STRa (qname), candidate, segconf.flav_prop[QNAME1].is_mated) && 
        !sam_is_depn (*mate_f) && 
        mate_f->is_last != f.is_last) {
        
        vb->mate_line_i = candidate;
        vb->mate_line_count++; // for stats

        dyn_int_append (VB, CTX(SAM_BUDDY), vb->line_i - candidate, 0); // add buddy (delta) >= 1 .
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
    if (!IS_MAIN(vb)) return BUDDY_NONE;

    LineIType candidate = LINE_BY_HASH(my_hash);
    SamFlags *saggy_f = (candidate >= 0) ? &DATA_LINE(candidate)->FLAG : NULL; // invalid pointer if no saggy

    // case: we found another member of the same sag (i.e. same qname, same is_last)
    if (has_same_qname (vb, STRa (qname), candidate, false) && saggy_f->is_last == f.is_last) { // the "prim" line against which we are segging cannot have hard clips

        vb->saggy_is_prim = !sam_is_depn (*saggy_f);

        // case: two alignments with the same QNAME, FLAG.is_first, and both are primary: usually a signof defective BAM 
        if (vb->saggy_is_prim && !sam_is_depn (f)) {
            WARN_ONCE ("WARNING: Non-compliant %s file: Found two alignments with the same QNAME, FLAG.IS_FIRST, and both are non-secondary, non-supplementary: QNAME=\"%.*s\"",
                       dt_name (vb->data_type), STRf(qname));
            goto no_saggy;
        }

        vb->saggy_line_i = candidate;
        vb->saggy_near_count++;

        // replace value in hash table with current line, so future lines from the same seg have a smaller buddy delta,
        // except if the line is a prim - we prefer the saggy to be a prim because there are many methods that require it
        if (!vb->saggy_is_prim)
            *insert_to_hash = true;

        dyn_int_append (VB, CTX(SAM_BUDDY), vb->line_i - candidate, 0); // add buddy (delta) >= 1 .
        return BUDDY_SAGGY;
    }

    // case: we didn't find another member of our sag in this VB
    else no_saggy: {
        *insert_to_hash = true;
        return BUDDY_NONE;
    }
}

static inline void sam_seg_QNAME_segconf (VBlockSAMP vb, ContextP ctx, STRp (qname))
{
    if (segconf.is_collated) { // still collated, try to find evidence to the contrary
        bool is_new = !is_same_last_txt (VB, ctx, STRa(qname));
        if (is_new && ctx->last_is_new) segconf.is_collated = false; // two new QNAMEs in a row = not collated

        // case: at least on pair of consecutive lines has the same QNAME. if the file is not sorted, we will set it as collated in sam_segconf_finalize
        if (!is_new) segconf.evidence_of_collated = true;

        ctx->last_is_new = is_new;
    }

    if (vb->line_i == 0) {
        qname_segconf_discover_flavor (VB, QNAME1, STRa (qname));

        huffman_start_chewing (SAM_QNAME, STRa(qname), 31);
    }

    // build huffman compressor based on segconf qname character probabilities, which is then used 
    // to compress qnames in SAGs (gencomp) in ZIP/PIZ and qnames in Deep in PIZ
    huffman_chew_one_sample (SAM_QNAME, 
                             STRa(qname), // lucky first qname will be the string to which we compare all other qnames hoping for an as long as possible identical prefix, as well as XOR the rest of the string
                             true);       // skip if qname is same as previous line, bc: 1. in Sag we won't store the qname again and 2. in Deep we don't store secondary/supplementary (not perfect stats for mates in sorted files through, as we skip some and some not - that's ok)
}

void sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAMP dl, STRp (qname), unsigned add_additional_bytes)
{
    decl_ctx (SAM_QNAME);

    // in segconf, identify if this file is collated (each QNAME appears in two or more consecutive lines)
    if (segconf.running) {
        sam_seg_QNAME_segconf (vb, ctx, STRa (qname));
        goto normal_seg;
    }

    QType q = qname_sam_get_qtype (STRa(qname)); // QNAME2 if we have QNAME2 and qname matches, else QNAME1

    uint32_t qname_hash = qname_calc_hash (q, COMP_NONE, STRa(qname), dl->FLAG.is_last, true, CRC32, NULL); // note: canonical=true as we use the same hash for find a mate and a saggy
    uint32_t my_hash    = qname_hash & MAXB(CTX(SAM_QNAME)->qname_hash.prm8[0]);
    bool insert_to_hash = false;

    if (flag.deep || flag.show_deep == SHOW_DEEP_ONE_HASH) 
        sam_deep_set_QNAME_hash (vb, dl, q, STRa(qname));

    BuddyType bt = sam_seg_mate  (vb, dl->FLAG, STRa (qname), my_hash, &insert_to_hash) | // bitwise or
                   sam_seg_saggy (vb, dl->FLAG, STRa (qname), my_hash, &insert_to_hash);

    dl->mate_line_i = vb->mate_line_i;

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
    if (bt) {
        seg_by_ctx (VB, (char[]){SNIP_SPECIAL, SAM_SPECIAL_COPY_BUDDY, '0' + bt}, 3, ctx, qname_len + add_additional_bytes); // seg QNAME as copy-from-buddy
        seg_set_last_txt (VB, ctx, STRa(qname));

        // mirror logic of cases that sam_piz_special_COPY_BUDDY sets buddy_copied_exactly
        if ((bt & BUDDY_MATE) && !segconf.flav_prop[QNAME1].is_mated)
            dl->qname_mate_copied_exactly = true;  
    }

    // case: DEPN with SA Group: seg against SA Group (unless already segged against buddy)
    else if (IS_DEPN(vb) && vb->sag) {
        sam_seg_against_sa_group (vb, ctx, qname_len + add_additional_bytes);
        seg_set_last_txt (VB, ctx, STRa(qname)); 
    }

    else normal_seg: {
        bool is_qname2 = qname_seg (VB, QNAME1/*must always start from QNAME1*/, STRa (qname), add_additional_bytes); // note: for PRIM component, this will be consumed with loading SA
        dl->is_consensus = is_qname2 && !IS_PRIM(vb)/*bug 949*/ && segconf.flav_prop[QNAME2].is_consensus;
    }

    // case: PRIM: additional seg against SA Group - store in SAM_QNAMESA - Reconstruct will take from here in PRIM per Toplevel container
    if (IS_PRIM(vb))
        seg_by_did (VB, (char[]){SNIP_SPECIAL, SAM_SPECIAL_PRIM_QNAME}, 2, SAM_QNAMESA, 0); // consumed when reconstructing PRIM vb
}

WordIndex sam_seg_RNAME (VBlockSAMP vb, ZipDataLineSAMP dl, STRp (chrom),
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
        random_access_update_chrom (VB, vb->chrom_node_index, STRa (chrom));

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
        random_access_update_chrom (VB, dl->RNAME, STRa (chrom));
        node_index = dl->RNAME;
    }

    else normal_seg: {
        bool is_new;
        node_index = chrom_seg_ex (VB, SAM_RNAME, STRa (chrom), 0, add_bytes, &is_new);
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

WordIndex sam_seg_RNEXT (VBlockSAMP vb, ZipDataLineSAMP dl, STRp (chrom), unsigned add_bytes)
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
        node_index = chrom_seg_ex (VB, SAM_RNEXT, STRa (chrom), 0, add_bytes, &is_new);
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

void sam_seg_init_bisulfite (VBlockSAMP vb, ZipDataLineSAMP dl)
{
    // the converted reference to which this read was mapped (C->T conversion or G->A conversion)
    // note: we calculate it always to avoid needless adding entropy in the snip
    vb->bisulfite_strand =  !segconf.sam_bisulfite    ? 0 
                            : IS_REF_INTERNAL         ? 0  // bug 648
                            : MP(BISMARK)   && has(XG_Z) ? sam_seg_get_aux_A (vb, vb->idx_XG_Z, IS_BAM_ZIP)
                            : MP(DRAGEN)    && has(XG_Z) ? sam_seg_get_aux_A (vb, vb->idx_XG_Z, IS_BAM_ZIP)
                            : MP(BSSEEKER2) && has(XO_Z) ? "CG"[sam_seg_get_aux_A (vb, vb->idx_XO_Z, IS_BAM_ZIP) == '-']
                            : MP(BSBOLT)    && has(YS_Z) ? "CG"[sam_seg_get_aux_A (vb, vb->idx_YS_Z, IS_BAM_ZIP) == 'C']
                            : MP(GEM3)      && has(XB_A) ? sam_seg_get_aux_A (vb, vb->idx_XB_A, IS_BAM_ZIP)
                            :                           0; // including MM/ML tags

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

    ZipDataLineSAMP dl = DATA_LINE(vb->line_i);

    str_split_by_tab (next_line, remaining_txt_len, MAX_FIELDS + AUX, has_13, false, true); // also advances next_line to next line
    
    ASSSEG (n_flds >= 11, "%s: (sam_seg_txt_line) Bad SAM file: alignment expected to have at least 11 fields, but found only %u", LN_NAME, n_flds);

    vb->n_auxs = n_flds - AUX;
    vb->auxs = &flds[AUX]; // note: pointers to data on the stack
    vb->aux_lens = &fld_lens[AUX];

    // support non-compliant SAM by e.g. ngmlr v0.2.7 - lines end with \t\n
    bool terminal_tab = vb->n_auxs && vb->aux_lens[vb->n_auxs-1] == 0;
    if (terminal_tab) vb->n_auxs--;
    ASSSEG0 (!terminal_tab || ! *has_13, "line ends with \\t\\r\\n: this is not currently supported by Genozip");

    sam_seg_idx_aux (vb);

    dl->SEQ.index = BNUMtxt (flds[SEQ]); // SEQ.len to be determined by sam_cigar_analyze
    vb->textual_seq_str = flds[SEQ]; 

    dl->QNAME = TXTWORDi (fld, QNAME);
    dl->QUAL  = TXTWORDi (fld, QUAL);

    if (fld_lens[QUAL] == 1 && flds[QUAL][0] == '*')
        vb->qual_missing = dl->no_qual = true;

    // lazy way to get vb->chrom*, rollback later if seg RNAME is not needed
    if (vb->check_for_gc || !IS_MAIN(vb))
        seg_create_rollback_point (VB, NULL, 1, SAM_RNAME);

    dl->RNAME = sam_seg_RNAME (vb, dl, STRfld (RNAME), false, fld_lens[RNAME] + 1);

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
        (flag.has_biopsy_line && sam_seg_test_biopsy_line (VB, flds[0], next_line - flds[0])))
        goto rollback_and_done;

    vb->last_cigar = flds[CIGAR];
    SAFE_NUL(&vb->last_cigar[fld_lens[CIGAR]]); // nul-terminate CIGAR string

    if (has(NM_i))
        dl->NM_len = sam_seg_get_aux_int (vb, vb->idx_NM_i, &dl->NM, false, MIN_NM_i, MAX_NM_i, HARD_FAIL) + 1; // +1 for \t or \n

    // set dl->AS needed by sam_seg_prim_add_sag (in PRIM) and several fields that delta against it
    if (has(AS_i))
        sam_seg_get_aux_int (vb, vb->idx_AS_i, &dl->AS, false, MIN_AS_i, MAX_AS_i, HARD_FAIL);

    if (!IS_MAIN(vb)) {
        sam_seg_sag_stuff (vb, dl, STRfld (CIGAR), flds[SEQ], false);

        // re-seg rname, against SA group
        seg_rollback (VB);
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

    // we analyze MD:Z now (if it exists) as we will need it for SEQ
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
        if      (IS_SAG_NH)   sam_seg_prim_add_sag_NH (vb, dl, dl->NH);
        else if (IS_SAG_CC)   sam_seg_prim_add_sag_CC (vb, dl, dl->NH);
        else if (IS_SAG_FLAG) sam_seg_prim_add_sag (vb, dl, 0, false);
        else if (IS_SAG_SOLO) sam_seg_prim_add_sag_SOLO (vb, dl);
    }

    // TLEN - must be after AUX as we might need data from MC:Z
    sam_seg_TLEN (vb, dl, STRfld (TLEN), 0, vb->RNEXT_is_equal);

    if (terminal_tab) seg_by_did (VB, "\t\n", 2, SAM_EOL, 1); // last field accounted for \n
    else              SEG_EOL (SAM_EOL, false);

    SAFE_RESTORE;            // restore \t after CIGAR

    MAXIMIZE (vb->longest_seq_len, dl->SEQ.len);

    return next_line;

rollback_and_done:
    memset (dl, 0, sizeof (ZipDataLineSAM));
    seg_rollback (VB); // cancelling segging of RNAME

    return next_line;
}
