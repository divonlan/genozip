// ------------------------------------------------------------------
//   sam_seg.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "sam_private.h"
#include "reference.h"
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
#include "libdeflate/libdeflate.h"

// Buddy & lookback parameters
typedef enum { QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, AUX } SamFields __attribute__((unused)); // quick way to define constants

char taxid_redirection_snip[100], xa_strand_pos_snip[100], copy_AS_snip[30], copy_NM_snip[30], copy_CR_snip[30],
     copy_buddy_CIGAR_snip[30], 
     copy_buddy_MAPQ_snip[30], copy_buddy_MQ_snip[30], XA_lookback_snip[30], copy_buddy_RNEXT_snip[30], copy_buddy_RNAME_snip[30],
     copy_buddy_YS_snip[30], copy_buddy_Z_snips[NUM_BUDDIED_Z_FIELDS][30], copy_buddy_AS_snip[30], copy_buddy_PNEXT_snip[100], copy_buddy_POS_snip[100],
     copy_buddy_NH_snip[30];
unsigned taxid_redirection_snip_len, xa_strand_pos_snip_len, copy_AS_snip_len, copy_NM_snip_len, copy_CR_snip_len, copy_buddy_CIGAR_snip_len,
     copy_buddy_MAPQ_snip_len, copy_buddy_MQ_snip_len, XA_lookback_snip_len, copy_buddy_RNEXT_snip_len, copy_buddy_RNAME_snip_len,
     copy_buddy_YS_snip_len, copy_buddy_Z_snip_lens[NUM_BUDDIED_Z_FIELDS], copy_buddy_AS_snip_len, copy_buddy_PNEXT_snip_len, copy_buddy_POS_snip_len,
     copy_buddy_NH_snip_len; 
WordIndex xa_lookback_strand_word_index = WORD_INDEX_NONE, xa_lookback_rname_word_index = WORD_INDEX_NONE;

Did buddied_Z_dids[NUM_BUDDIED_Z_FIELDS] = BUDDIED_Z_DIDs;

// called by zfile_compress_genozip_header to set FlagsGenozipHeader.dt_specific
bool sam_zip_dts_flag (void)
{
    return flag.reference == REF_INTERNAL;
}

// ----------------------
// Seg stuff
// ----------------------

void sam_zip_free_end_of_z (void)
{
    sam_header_finalize(); 
}

// main thread, called for each component. called before segconf.
void sam_zip_initialize (void)
{
    bool has_hdr_contigs = sam_hdr_contigs && sam_hdr_contigs->contigs.len;
    
    // Copy header contigs to RNAME and RNEXT upon first component. This is in the order of the
    // header, as required by BAM (it encodes ref_id based on header order). Note, subsequent
    // bound files are required to have the same contigs or at least a contiguous subset starting a contig 0.
    // BAM always has header contigs (might be 0 of them, for an unaligned file), while SAM is allowed to be header-less
    if (z_file->num_txts_so_far == 1 && has_hdr_contigs) { 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, sam_hdr_contigs); 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNEXT, sam_hdr_contigs);
    }

    // with REF_EXTERNAL and unaligned data, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (z_file->num_txts_so_far == 1 && flag.aligner_available && flag.reference == REF_EXTERNAL) 
        ctx_populate_zf_ctx_from_contigs (gref, SAM_RNAME, ref_get_ctgs (gref)); 

    seg_prepare_snip_other (SNIP_REDIRECTION, _SAM_TAXID, false, 0, taxid_redirection_snip);

    static SmallContainer xa_strand_pos_con = { .repeats=1, .nitems_lo=2, .items = { { { _OPTION_XA_STRAND } }, { { _OPTION_XA_POS } } } };
    xa_strand_pos_snip_len = sizeof (xa_strand_pos_snip);
    container_prepare_snip ((ConstContainerP)&xa_strand_pos_con, 0, 0, xa_strand_pos_snip, &xa_strand_pos_snip_len); 

    qname_zip_initialize (SAM_QNAME);
    sam_zip_QUAL_initialize();
        
    seg_prepare_snip_other (SNIP_COPY,       _OPTION_AS_i, false, 0, copy_AS_snip);
    seg_prepare_snip_other (SNIP_COPY,       _OPTION_NM_i, false, 0, copy_NM_snip);
    seg_prepare_snip_other (SNIP_COPY,       _OPTION_CR_Z, false, 0, copy_CR_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_PNEXT,   false, 0, copy_buddy_PNEXT_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_POS,     false, 0, copy_buddy_POS_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_RNAME,   false, 0, copy_buddy_RNAME_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_RNEXT,   false, 0, copy_buddy_RNEXT_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _SAM_MAPQ,    false, 0, copy_buddy_MAPQ_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_MQ_i, false, 0, copy_buddy_MQ_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_YS_i, false, 0, copy_buddy_YS_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_AS_i, false, 0, copy_buddy_AS_snip);
    seg_prepare_snip_other (SNIP_COPY_BUDDY, _OPTION_NH_i, false, 0, copy_buddy_NH_snip);
    seg_prepare_snip_other (SNIP_LOOKBACK,   _OPTION_XA_LOOKBACK, false, 0, XA_lookback_snip);
    seg_prepare_snip_special_other (SAM_SPECIAL_RECON_BUDDY_CIGAR, copy_buddy_CIGAR_snip, _SAM_CIGAR);

    for (BuddiedZFields f=0; f < NUM_BUDDIED_Z_FIELDS; f++) {
        copy_buddy_Z_snip_lens[f] = sizeof (copy_buddy_Z_snips[f]);
        seg_prepare_snip_other_do (SNIP_COPY_BUDDY, ZCTX(buddied_Z_dids[f])->dict_id, false, 0, 0, copy_buddy_Z_snips[f], &copy_buddy_Z_snip_lens[f]);
    }

    strcpy (ZCTX(OPTION_UR_Z)->tag_name, "UR_UB:Z"); // for stats - as UB:Z is an alias of UR:Z - both these data go into UR:Z
}

// called after each file
void sam_zip_finalize (void)
{
    gencomp_destroy();
}

// called by main thread after reading txt of one vb into vb->txt_data
void sam_zip_init_vb (VBlockP vb)
{
    vb->chrom_node_index = NODE_INDEX_NONE;    

    // note: we test for sorted and not collated, because we want non-sorted long read files (which are collated)
    // to seg depn against same-VB prim (i.e. not gencomp) - as the depn lines will follow the prim line
    VB_SAM->check_for_gc = (!segconf.running && segconf.sag_type && sam_is_main_vb);

    sam_sag_zip_init_vb (vb);
}

// called compute thread after compress, order of VBs is arbitrary
void sam_zip_after_compress (VBlockP vb)
{
    // Only the MAIN component produces gencomp lines, however we are processing VBs in order, so out-of-band VBs
    // need to be sent too, just to advance the serializing mutex
    if (segconf.sag_type && (sam_is_main_vb || sam_is_prim_vb))
        gencomp_absorb_vb_gencomp_lines (vb);
}

// called main thread, in order of VBs
void sam_zip_after_compute (VBlockP vb)
{
    if (vb->comp_i == SAM_COMP_MAIN) 
        sam_zip_gc_after_compute_main (VB_SAM);

    else if (vb->comp_i == SAM_COMP_PRIM) 
        gencomp_sam_prim_vb_has_been_ingested (vb);
}

// main thread: writing data-type specific fields to genozip header
void sam_zip_genozip_header (SectionHeaderGenozipHeader *header)
{
    header->sam.segconf_seq_len         = BGEN32 (segconf.sam_seq_len);  // added v14
    header->sam.segconf_seq_len_cm      = segconf.seq_len_to_cm;         // v14
    header->sam.segconf_ms_type         = segconf.sam_ms_type;           // v14
    header->sam.segconf_has_MD_or_NM    = segconf.has_MD_or_NM;          // v14
    header->sam.segconf_bisulfite       = segconf.sam_bisulfite;         // v14
    header->sam.segconf_is_paired       = segconf.is_paired;             // v14
    header->sam.segconf_sag_type        = segconf.sag_type;              // v14
    header->sam.segconf_sag_has_AS      = segconf.sag_has_AS;       // v14
    header->sam.segconf_AS_is_2refc     = segconf.AS_is_2ref_consumed;   // v14
    header->sam.segconf_pysam_qual      = segconf.pysam_qual;            // v14
    header->sam.segconf_seq_len_dict_id = segconf.qname_seq_len_dict_id; // v14
}

// initialize SA and OA
static void sam_seg_0X_initialize (VBlockP vb, Did strand_did_i)
{
    // create strand nodes (nodes will be deleted in sam_seg_finalize if not used)
    ctx_create_node (vb, strand_did_i, cSTR("-"));
    ctx_create_node (vb, strand_did_i, cSTR("+"));
    CTX(strand_did_i)->no_vb1_sort = true; // keep them in this ^ order
}

static void sam_seg_qname_initialize (VBlockSAMP vb)
{
    CTX(SAM_QNAME)->no_stons = true;             // no singletons, bc sam_piz_filter uses PEEK_SNIP
    CTX(SAM_QNAME)->flags.store_per_line = true; // 12.0.41 

    // a bug that existed 12.0.41-13.0.1 (bug 367): we stored buddy in machine endianty instead of BGEN32.
    // we use local.prm8[0]=1 to indicate to reconstruct_set_buddy that this bug is now fixed.
    CTX(SAM_BUDDY)->local.prm8[0] = 1;
    CTX(SAM_BUDDY)->local_param = true;
    CTX(SAM_BUDDY)->ltype = LT_DYN_INT;

    qname_seg_initialize (VB, SAM_QNAME);

    if (segconf.running) 
        segconf.qname_flavor = 0; // unknown

    // initial allocations based on segconf data
    else {
        vb->qname_hash.prm8[0] = MIN_(20, MAX_(14, 32 - __builtin_clzl (vb->lines.len32 * 5))); // between 14 and 20 bits - tested - no additional compression benefit beyond 20 bits
        buf_alloc_255 (vb, &vb->qname_hash, 0, (1ULL << vb->qname_hash.prm8[0]), int32_t, 1, "qname_hash");    
    }
}

void sam_seg_initialize (VBlockP vb_)
{
    START_TIMER;
    VBlockSAMP vb = (VBlockSAMP)vb_;

    // all numeric fields need STORE_INT / STORE_FLOAT to be reconstructable to BAM (possibly already set)
    // via the translators set in the SAM_TOP2BAM Container
    #define T(cond, did_i) ((cond) ? (did_i) : DID_NONE)
    ctx_set_store (VB, STORE_INT, 41, SAM_TLEN, SAM_MAPQ, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_GPOS,
                   OPTION_NM_i, OPTION_AS_i, OPTION_MQ_i, OPTION_XS_i, OPTION_XM_i, OPTION_mc_i, OPTION_ms_i, OPTION_Z5_i, 
                   OPTION_tx_i, OPTION_YS_i, OPTION_XC_i, OPTION_AM_i, OPTION_SM_i, OPTION_X0_i, OPTION_X1_i, OPTION_CP_i,
                   OPTION_OP_i, OPTION_NH_i, OPTION_HI_i, OPTION_cm_i, OPTION_SA_POS, OPTION_OA_POS,
                   T(segconf.tech == TECH_PACBIO, OPTION_qs_i), T(segconf.tech == TECH_PACBIO, OPTION_qe_i),
                   T(MP(BLASR), OPTION_XS_i), T(MP(BLASR), OPTION_XE_i), T(MP(BLASR), OPTION_XQ_i), T(MP(BLASR), OPTION_XL_i), T(MP(BLASR), OPTION_FI_i),
                   T(MP(NGMLR), OPTION_QS_i), T(MP(NGMLR), OPTION_QE_i), T(MP(NGMLR), OPTION_XR_i), 
                   T(MP(MINIMAP2) || MP(WINNOWMAP), OPTION_s1_i),
                   T(kraken_is_loaded, SAM_TAXID),
                   DID_EOL);

    ctx_set_store (VB, STORE_INDEX, 3, SAM_RNAME, SAM_RNEXT, DID_EOL); // when reconstructing BAM, we output the word_index instead of the string
    
    ctx_set_delta_peek (VB, 3, OPTION_AS_i, OPTION_s1_i, DID_EOL); // when eg AS is reconstructed as DELTA_OTHER, the other (ms:i) may be before or after

    // don't store singletons in local. note: automatically implied if ltype!=LT_TEXT is set 
    ctx_set_no_stons (VB, 20,
                      SAM_RNAME, SAM_RNEXT,      // BAM reconstruction needs RNAME, RNEXT word indices. also needed for random access.
                      OPTION_MD_Z, 
                      OPTION_BI_Z, OPTION_BD_Z,  // we can't use local for singletons in BD or BI as next_local is used by sam_piz_special_BD_BI to point into BD_BI
                      OPTION_SA_CIGAR, OPTION_XA_CIGAR, OPTION_OA_CIGAR, OPTION_OC_Z, // we can't use local for singletons bc sam_seg_other_CIGAR manually offloads CIGARs to local
                      SAM_POS, SAM_PNEXT, OPTION_mc_i, OPTION_OP_i, OPTION_Z5_i, OPTION_CP_i, // required by seg_pos_field
                      T(MP(BSSEEKER2), OPTION_XM_Z), T(MP(BSSEEKER2), OPTION_XG_Z),
                      T(MP(BSBOLT), OPTION_XB_Z),
                      T(kraken_is_loaded, SAM_TAXID),
                      DID_EOL);

    // MAPQ (and hence other fields carrying mapping quality) is uint8_t by BAM specification
    ctx_set_ltype (VB, LT_UINT8, 6, SAM_MAPQ, OPTION_SA_MAPQ, OPTION_SM_i, OPTION_AM_i, OPTION_OA_MAPQ, DID_EOL);
    
    // initialize these to LT_SEQUENCE, the qual-type ones might be changed later to LT_CODEC (eg domq)
    ctx_set_ltype (VB, LT_SEQUENCE, 4, OPTION_BD_BI, OPTION_UY_Z, OPTION_CY_ARR, DID_EOL);

    // set ltype=LT_DYN_INT to allow use of seg_integer
    ctx_set_ltype (VB, LT_DYN_INT, 10, OPTION_HI_i, OPTION_NM_i, OPTION_NH_i, OPTION_NH_PRIM,  OPTION_XM_i, OPTION_X1_i, OPTION_AS_i,
                   OPTION_cm_i,
                   T(segconf.has_TLEN_non_zero, SAM_TLEN), // note: we don't set if !has_TLEN_non_zero, bc values are stored in b250 and may require singletons
                   DID_EOL);

    ctx_set_ltype (VB, LT_UINT32, 2, OPTION_CP_i, DID_EOL);

    if (segconf.is_collated) 
        CTX(SAM_POS)->flags.store_delta = true; // since v12.0.41
    
    // we may use mates (other than for QNAME) if not is_long_reads (meaning: no mates in this file) and not DEPN components (bc we seg against PRIM)
    if (segconf.is_paired && !sam_is_depn_vb) 
        ctx_set_store_per_line (VB, 11, SAM_RNAME, SAM_RNEXT, SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_MAPQ, SAM_CIGAR,
                                OPTION_MQ_i, OPTION_MC_Z, OPTION_SM_i, DID_EOL);

    // case: some rows may be segged against an in-VB prim line
    if (sam_is_main_vb && !VB_SAM->check_for_gc)  // 14.0.0
        ctx_set_store_per_line (VB, 10, SAM_RNAME, SAM_RNEXT, SAM_PNEXT, SAM_POS, SAM_CIGAR, SAM_MAPQ, SAM_FLAG,
                                OPTION_SA_Z, OPTION_NM_i, DID_EOL);

    ctx_set_store_per_line (VB, 3, OPTION_NH_i, T(segconf.is_paired && segconf.sam_multi_RG, OPTION_RG_Z), DID_EOL); 

    if (segconf.is_paired)
        for (BuddiedZFields f=1; f < NUM_BUDDIED_Z_FIELDS; f++) // note: f starts from 1, bc RG is set above with T()
            ctx_set_store_per_line (VB, 2, buddied_Z_dids[f], DID_EOL); 

    if (kraken_is_loaded) 
        CTX(SAM_TAXID)->counts_section = true; 

    // in --stats, consolidate stats 
    ctx_consolidate_stats (VB, SAM_SQBITMAP, 9, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND, SAM_SEQMIS_A, SAM_SEQMIS_C, SAM_SEQMIS_G, SAM_SEQMIS_T, DID_EOL);
    ctx_consolidate_stats (VB, SAM_QUAL,     5, SAM_DOMQRUNS, SAM_QUALMPLX, SAM_DIVRQUAL, SAM_QUALSA, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_OQ_Z,  4, OPTION_OQ_DOMQRUNS, OPTION_OQ_QUALMPLX, OPTION_OQ_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_UY_Z,  4, OPTION_UY_DOMQRUNS, OPTION_UY_QUALMPLX, OPTION_UY_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_CY_Z,  5, OPTION_CY_ARR, OPTION_CY_DOMQRUNS, OPTION_CY_QUALMPLX, OPTION_CY_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_U2_Z,  4, OPTION_U2_DOMQRUNS, OPTION_U2_QUALMPLX, OPTION_U2_DIVRQUAL, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_E2_Z,  5, OPTION_2NONREF, OPTION_N2ONREFX, OPTION_2GPOS, OPTION_S2TRAND, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_OA_Z,  7, OPTION_OA_RNAME, OPTION_OA_POS, OPTION_OA_STRAND, OPTION_OA_CIGAR, OPTION_OA_MAPQ, OPTION_OA_NM, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_XA_Z,  8, OPTION_XA_RNAME, OPTION_XA_POS, OPTION_XA_STRAND, OPTION_XA_CIGAR, OPTION_XA_NM, OPTION_XA_STRAND_POS, OPTION_XA_LOOKBACK, DID_EOL);
    ctx_consolidate_stats (VB, SAM_AUX,      2, SAM_MC_Z, DID_EOL); // note: this is *not* OPTION_MC_Z
    ctx_consolidate_stats (VB, SAM_QNAME,    3, SAM_BUDDY, SAM_QNAMESA, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_BD_BI, 3, OPTION_BI_Z, OPTION_BD_Z, DID_EOL);
    ctx_consolidate_stats (VB, OPTION_CR_CB, 3, OPTION_CR_Z, OPTION_CB_Z, DID_EOL);

    if (segconf.has[OPTION_NH_i]) { // uses NH:i / HI:i instead of SA
        ctx_consolidate_stats (VB, OPTION_HI_i, 2, OPTION_SA_Z, DID_EOL);
        ctx_consolidate_stats (VB, OPTION_NH_i, 2, OPTION_NH_PRIM, DID_EOL);
    }
    else
        ctx_consolidate_stats (VB, OPTION_SA_Z, 8, OPTION_SA_RNAME, OPTION_SA_POS, OPTION_SA_STRAND, OPTION_SA_CIGAR, OPTION_SA_MAPQ, OPTION_SA_NM, OPTION_SA_MAIN, DID_EOL);

    codec_acgt_comp_init (VB);
    sam_seg_qname_initialize (vb);
    sam_seg_QUAL_initialize (vb);
    sam_seg_SEQ_initialize (vb);
    sam_seg_cigar_initialize (vb);
    sam_seg_gc_initialize (vb);
    sam_seg_0X_initialize (VB, OPTION_SA_STRAND);
    sam_seg_0X_initialize (VB, OPTION_OA_STRAND);

    if (sam_has_BWA_XA_Z())
        sam_seg_BWA_XA_initialize (vb);

    ctx_set_store (VB, STORE_INDEX, 2, OPTION_XA_Z, DID_EOL); // for containers this stores repeats - used by sam_piz_special_X1->container_peek_repeats

    if (sam_has_BWA_XS_i()) // XS:i is as defined some aligners
        seg_mux_init (VB, 4, SAM_SPECIAL_BWA_XS, OPTION_XS_i, OPTION_XS_i, STORE_INT, LT_TEXT, false, (MultiplexerP)&vb->mux_XS, "0123");

    else if (MP(HISAT2)) // ZS:i is like BWA's XS:i
        seg_mux_init (VB, 4, SAM_SPECIAL_BWA_XS, OPTION_ZS_i, OPTION_ZS_i, STORE_INT, LT_TEXT, false, (MultiplexerP)&vb->mux_XS, "0123");

    if (sam_has_bowtie2_YS_i()) {
        ctx_set_store_per_line (VB, 3, OPTION_AS_i, OPTION_YS_i, DID_EOL);
        seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_MATE,  OPTION_YS_i, OPTION_YS_i, STORE_INT,  LT_TEXT, false, (MultiplexerP)&vb->mux_YS, "01");
    }

    seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_BUDDY,     SAM_FLAG,    SAM_FLAG,    STORE_INT,  LT_TEXT,    false, (MultiplexerP)&vb->mux_FLAG, "01");
    seg_mux_init (VB, 3, SAM_SPECIAL_DEMUX_BY_MATE_PRIM, SAM_POS,     SAM_POS,     STORE_INT,  LT_TEXT,    true,  (MultiplexerP)&vb->mux_POS, "012");
    seg_mux_init (VB, 4, SAM_SPECIAL_PNEXT,              SAM_PNEXT,   SAM_PNEXT,   STORE_INT,  LT_TEXT,    true,  (MultiplexerP)&vb->mux_PNEXT, "0123");
    seg_mux_init (VB, 3, SAM_SPECIAL_DEMUX_BY_MATE_PRIM, SAM_MAPQ,    SAM_MAPQ,    STORE_INT,  LT_TEXT,    false, (MultiplexerP)&vb->mux_MAPQ, "012");
    seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_MATE,      OPTION_MQ_i, OPTION_MQ_i, STORE_INT,  LT_TEXT,    false, (MultiplexerP)&vb->mux_MQ, "01");
    seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_MATE,      OPTION_MC_Z, OPTION_MC_Z, STORE_NONE, LT_TEXT,    false, (MultiplexerP)&vb->mux_MC, "01");
    seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_MATE,      OPTION_ms_i, OPTION_ms_i, STORE_INT,  LT_TEXT,    false, (MultiplexerP)&vb->mux_ms, "01");
    seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_MATE,      OPTION_AS_i, OPTION_AS_i, STORE_INT,  LT_DYN_INT, false, (MultiplexerP)&vb->mux_AS, "01");

    for (BuddiedZFields f=0; f < NUM_BUDDIED_Z_FIELDS; f++)
        seg_mux_init (VB, 2, SAM_SPECIAL_DEMUX_BY_BUDDY, buddied_Z_dids[f], buddied_Z_dids[f], STORE_NONE, LT_TEXT, false, (MultiplexerP)&vb->mux_buddied_z_fields[f], "01");
    
    // needed in QUAL-type fields that might have '*' for missing quality
    CTX(SAM_QUAL)->flags.is_qual = CTX(OPTION_OQ_Z)->flags.is_qual = CTX(OPTION_U2_Z)->flags.is_qual = true;

    if (segconf.running) {
        segconf.is_sorted = segconf.is_collated = segconf.MAPQ_has_single_value = segconf.NM_after_MD = true; // initialize optimistically
        segconf.sam_is_unmapped = true;  // we will reset this if finding a line with POS>0
    }
    
    COPY_TIMER (seg_initialize);
    #undef T
}

static void sam_seg_toplevel (VBlockP vb)
{
    // in PRIM, SA group loader preprocessor reconstructs QNAME / SQBITMAP / QUAL to load SA Groups, while recon
    // reconstructs QNAMESA, SEQSA, QUALSA to copy from SA Groups. In contrast, in DEPN, copy from SA Groups 
    // occurs (if needed) when reconstructing QNAME / SQBITMAP / QUAL.
    uint64_t qname_dict_id = (sam_is_prim_vb ? _SAM_QNAMESA : _SAM_QNAME);
    uint64_t seq_dict_id   = _SAM_SQBITMAP; 
    uint64_t qual_dict_id  = (sam_is_prim_vb ? _SAM_QUALSA  : _SAM_QUAL);
    
    // top level snip - reconstruction as SAM
    SmallContainer top_level_sam = { 
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 13,
        .items        = { { .dict_id = { qname_dict_id }, .separator = "\t" },
                          { .dict_id = { _SAM_FLAG     }, .separator = "\t" },
                          { .dict_id = { _SAM_RNAME    }, .separator = "\t" },
                          { .dict_id = { _SAM_POS      }, .separator = "\t" },
                          { .dict_id = { _SAM_MAPQ     }, .separator = "\t" },
                          { .dict_id = { _SAM_CIGAR    }, .separator = "\t" },
                          { .dict_id = { _SAM_RNEXT    }, .separator = "\t" },
                          { .dict_id = { _SAM_PNEXT    }, .separator = "\t" },
                          { .dict_id = { _SAM_TLEN     }, .separator = "\t" },
                          { .dict_id = { seq_dict_id   }, .separator = "\t" },
                          { .dict_id = { qual_dict_id  }, .separator = "\t" },
                          { .dict_id = { _SAM_AUX      },                   },
                          { .dict_id = { _SAM_EOL      },                   } 
                        }
    };
    container_seg (vb, CTX(SAM_TOPLEVEL), (ContainerP)&top_level_sam, 0, 0, 0);

    // top level snip - reconstruction as BAM
    // strategy: we start by reconstructing the variable-length fields first (after a prefix that sets them in place) 
    // - read_name, cigar, seq and qual - and then go back and fill in the fixed-location fields
    // Translation (a feature of Container): items reconstruct their data and then call a translation function to translate it to the desired format
    SmallContainer top_level_bam = { 
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .callback     = true,
        .filter_items = true,
        .nitems_lo    = 14,
        .items        = { { .dict_id = { _SAM_RNAME    }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                          { .dict_id = { _SAM_POS      }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 1 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                          { .dict_id = { _SAM_MAPQ     }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_U8       }, // Translate - textual to binary number
                          { .dict_id = { _SAM_BAM_BIN  }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 2 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                          { .dict_id = { _SAM_FLAG     }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 4 }, SAM2BAM_LTEN_U16 }, // Translate - textual to binary number
                          { .dict_id = { _SAM_RNEXT    }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_RNAME    }, // Translate - output word_index instead of string
                          { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_TRANS_NOR | CI0_TRANS_MOVE, 4 }, SAM2BAM_POS      }, // Translate - output little endian POS-1
                          { .dict_id = { qname_dict_id }, .separator = { CI0_TRANS_NUL                     }                   }, // normal 
                          { .dict_id = { _SAM_CIGAR    }, .separator = ""                                                      }, // handle in special reconstructor - translate textual to BAM CIGAR format + reconstruct l_read_name, n_cigar_op, l_seq
                          { .dict_id = { _SAM_TLEN     }, .separator = { CI0_TRANS_NOR                     }, SAM2BAM_TLEN     }, // must be after CIGAR bc sam_piz_special_TLEN_old needs vb->seq_num
                          { .dict_id = { seq_dict_id   }, .separator = "",                                    SAM2BAM_SEQ      }, // Translate - textual format to BAM format
                          { .dict_id = { qual_dict_id  }, .separator = "",                                    SAM2BAM_QUAL     }, // Translate - textual format to BAM format, set block_size
                          { .dict_id = { _SAM_AUX      }, .separator = { CI0_TRANS_NOR                     }                   }, // up to v11, this had the SAM2BAM_AUX translator
                        }
    };

    // no container wide-prefix, skip l_name with a 4-character prefix
    static const char bam_line_prefix[] = { CON_PX_SEP, // has prefix 
                                            CON_PX_SEP, // end of (empty) container-wide prefix
                                            ' ',' ',' ',' ', CON_PX_SEP }; // first item prefix - 4 spaces (place holder for block_size)

    container_seg (vb, CTX(SAM_TOP2BAM), (ContainerP)&top_level_bam, bam_line_prefix, sizeof(bam_line_prefix), 
                   IS_BAM_ZIP ? sizeof (uint32_t) * vb->lines.len : 0); // if BAM, account for block_size

    // top level snip - reconstruction as FASTQ
    SmallContainer top_level_fastq = { 
        .repeats     = vb->lines.len32,
        .is_toplevel = true,
        .callback    = true,  // drop supplementary alignments and alignments without QUAL data
        .nitems_lo   = 9,
        .items       = { { .dict_id = { qname_dict_id }, .separator = "\n"                  }, 
                         { .dict_id = { _SAM_RNAME    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing seq 
                         { .dict_id = { _SAM_POS      }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing seq
                         { .dict_id = { _SAM_RNEXT    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing PNEXT 
                         { .dict_id = { _SAM_PNEXT    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing POS (in case of BUDDY)
                         { .dict_id = { _SAM_FLAG     }, .separator = { CI0_TRANS_NOR }, .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_CIGAR    }, .separator = { CI0_TRANS_NOR }     }, // needed for reconstructing seq
                         { .dict_id = { _SAM_MC_Z     }, .separator = { CI0_TRANS_NOR }     }, // consumes OPTION_MC_Z if its on this line, might be needed for mate CIGAR
                         { .dict_id = { seq_dict_id   }, .separator = "\n",              .translator = SAM2FASTQ_SEQ  }, 
                         { .dict_id = { qual_dict_id  }, .separator = "\n",              .translator = SAM2FASTQ_QUAL }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line    
    static const char fastq_line_prefix[] = { CON_PX_SEP, CON_PX_SEP, '@', CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, 
                                              CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, CON_PX_SEP, '+', '\n', CON_PX_SEP };

    container_seg (vb, CTX(SAM_TOP2FQ), (ContainerP)&top_level_fastq, fastq_line_prefix, sizeof(fastq_line_prefix), 0);

    // top level snip - reconstruction as FASTQ "extended" - with all the SAM fields in the description line
    SmallContainer top_level_fastq_ext = { 
        .repeats     = vb->lines.len32,
        .is_toplevel = true,
        .callback    = true,  // drop non-primary chimeric reads and reads without QUAL data
        .nitems_lo   = 12,
        .items       = { { .dict_id = { qname_dict_id }, .separator = "\t"                               }, 
                         { .dict_id = { _SAM_FLAG     }, .separator = "\t", .translator = SAM2FASTQ_FLAG }, // need to know if seq is reverse complemented & if it is R2 ; reconstructs "1" for R1 and "2" for R2
                         { .dict_id = { _SAM_RNAME    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_POS      }, .separator = "\t"                               },
                         { .dict_id = { _SAM_MAPQ     }, .separator = "\t"                               },
                         { .dict_id = { _SAM_CIGAR    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_RNEXT    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_PNEXT    }, .separator = "\t"                               },
                         { .dict_id = { _SAM_TLEN     }, .separator = "\t"                               },
                         { .dict_id = { _SAM_AUX      }, .separator = "\n"                               },
                         { .dict_id = { seq_dict_id   }, .separator = "\n", .translator =SAM2FASTQ_SEQ   }, 
                         { .dict_id = { qual_dict_id  }, .separator = "\n", .translator =SAM2FASTQ_QUAL  }, // also moves fastq "line" to R2 (paired file) if needed
                       }
    };

    // add a '@' to the description line, use a prefix to add the + line    
    static const char fastq_ext_line_prefix[] = { 
        CON_PX_SEP, CON_PX_SEP, 
        '@', CON_PX_SEP, 
        'F','L','A','G',':', CON_PX_SEP, 
        'R','N','A','M','E',':', CON_PX_SEP, 
        'P','O','S',':', CON_PX_SEP, 
        'M','A','P','Q',':', CON_PX_SEP, 
        'C','I','G','A','R',':', CON_PX_SEP,
        'R','N','E','X','T',':', CON_PX_SEP, 
        'P','N','E','X','T',':', CON_PX_SEP, 
        'T','L','E','N',':', CON_PX_SEP, 
        CON_PX_SEP,              // optional
        CON_PX_SEP,              // SEQ
        '+', '\n', CON_PX_SEP }; // QUAL

    container_seg (vb, CTX(SAM_TOP2FQEX), (ContainerP)&top_level_fastq_ext, 
                   fastq_ext_line_prefix, sizeof(fastq_ext_line_prefix), 0);
}

// finalize Seg configuration parameters    
static void sam_seg_finalize_segconf (VBlockP vb)
{
    segconf.sam_multi_RG  = CTX(OPTION_RG_Z)->nodes.len32 >= 2; 
    segconf.sam_cigar_len = 1 + ((segconf.sam_cigar_len-1) / vb->lines.len32); // set to the average CIGAR len (rounded up)
    segconf.sam_seq_len   = (uint32_t)(0.5 + (double)segconf.sam_seq_len / (double)vb->lines.len32); // set average seq_len - rounded to the nearest 
    
    segconf.seq_len_to_cm /= vb->lines.len32;
    if (segconf.seq_len_to_cm > 255) segconf.seq_len_to_cm = 0;

    // in some STAR transcriptome files, there are no @PG header lines. it is important that mapper is set to STAR so we can use liberal prim_line_i, as SA groups are large
    if (MP(UNKNOWN) && segconf.has[OPTION_NH_i] && segconf.has[OPTION_HI_i] && !segconf.has[OPTION_SA_Z])
        segconf.sam_mapper = MP_STAR;

    // cases where its beneficial seg AS:i against prim
    if (MP(STAR))
        segconf.sag_has_AS = true; 

    // we set AS_is_2ref_consumed, meaning AS tends to be near 2 X ref_consumed, if at least half of the lines say so
    segconf.AS_is_2ref_consumed = (segconf.AS_is_2ref_consumed > vb->lines.len32 / 2);

    segconf.sam_bisulfite = MP(BISMARK) || MP(BSSEEKER2) ||
                            MP(BSBOLT)  ||(MP(GEM3) && segconf.has[OPTION_XB_A]);

    if (MP(MINIMAP2) || MP(WINNOWMAP))
        segconf.sam_ms_type = ms_MINIMAP2;

    segconf.is_long_reads = segconf_is_long_reads();

    // if we have @HD-SO "coordinate" or "queryname", then we take that as definitive. Otherwise, we go by our segconf sampling.
    if (sam_hd_so == HD_SO_COORDINATE) {
        segconf.is_sorted   = true;
        segconf.is_collated = false;
    }

    else if (sam_hd_so == HD_SO_QUERYNAME) {
        segconf.is_collated = true;
        segconf.is_sorted   = false;
    }

    else {
        // case: if we haven't found any pair of consecutive lines with the same RNAME and non-descreasing POS, this is not a sorted file, despite no evidence of "not sorted". eg could be unique RNAMEs.
        if (!segconf.evidence_of_sorted)
            segconf.is_sorted = false;

        // we have at leaest one pair of lines with the same QNAME, and the file is not sorted
        if (!segconf.is_sorted && segconf.evidence_of_collated)
            segconf.is_collated = true; 
    }

    // second is_paired test: PNEXT is not 0 for all lines (first test is in sam_seg_FLAG)
    if (segconf.has[SAM_PNEXT])
        segconf.is_paired = true;

    if (segconf.is_long_reads) {
        if (MP(MINIMAP2) || MP(WINNOWMAP))
            segconf.sam_ms_type = ms_MINIMAP2; // definitely not biobambam's MateBaseScore as long reads don't have mates
    }
    else {
        if (segconf.has[OPTION_ms_i] && segconf.is_paired)
            segconf.sam_ms_type = ms_BIOBAMBAM;
    }

    // evidence of STARsolo or cellranger (this is also detected in the SAM header)
    if (segconf.has[OPTION_UR_Z] + segconf.has[OPTION_UY_Z] + segconf.has[OPTION_UB_Z] >= 2 ||
        (segconf.has[OPTION_gn_Z] && segconf.has[OPTION_gx_Z]) ||
        (segconf.has[OPTION_GN_Z] && segconf.has[OPTION_GX_Z]))
        segconf.star_solo = true; // must be before sam_segconf_set_sag_type

    sam_segconf_set_sag_type();
}

void sam_seg_finalize (VBlockP vb)
{
    // We always include the SQBITMAP local section, except if no lines
    if (vb->lines.len)
        CTX(SAM_SQBITMAP)->local_always = true;

    // assign the QUAL codec
    if (VB_SAM->has_qual)
        codec_assign_best_qual_codec (vb, SAM_QUAL, sam_zip_qual, VB_SAM->qual_codec_no_longr);

    if (CTX(OPTION_OQ_Z)->local.len32)
        codec_assign_best_qual_codec (vb, OPTION_OQ_Z, sam_zip_OQ, VB_SAM->qual_codec_no_longr);
    
    if (CTX(OPTION_CY_ARR)->local.len32)
        codec_assign_best_qual_codec (vb, OPTION_CY_ARR, NULL, true);
    
    if (CTX(OPTION_UY_Z)->local.len32)
        codec_assign_best_qual_codec (vb, OPTION_UY_Z, sam_zip_UY, true);
    
    if (CTX(OPTION_U2_Z)->local.len32)
        codec_assign_best_qual_codec (vb, OPTION_U2_Z, sam_zip_U2, true);

    // determine if sam_piz_sam2bam_SEQ ought to store vb->textual_seq
    CTX(SAM_SQBITMAP)->flags.no_textual_seq = CTX(SAM_QUAL)->lcodec != CODEC_LONGR
                                           && segconf.sam_mapper != MP_BSSEEKER2;

    if (flag.biopsy_line.line_i == NO_LINE) // no --biopsy-line
        sam_seg_toplevel (vb);

    // finalize Seg configuration parameters    
    if (segconf.running) sam_seg_finalize_segconf (vb);
  
    // get rid of the 0A strand contexts (nodes added in sam_seg_0X_initialize) if we ended up not have using them
    // (note: we use SA_STRAND in PRIM lines even if there is no SA:Z field)
    if (!CTX(OPTION_XA_STRAND)->b250.len) ctx_free_context (CTX(OPTION_XA_STRAND), OPTION_XA_STRAND);
    if (!CTX(OPTION_OA_STRAND)->b250.len) ctx_free_context (CTX(OPTION_OA_STRAND), OPTION_OA_STRAND);
    if (!CTX(OPTION_SA_STRAND)->b250.len) ctx_free_context (CTX(OPTION_SA_STRAND), OPTION_SA_STRAND);

    // case: we have XA but didn't use lookback (because its a non-sorted file) - create an all-the-same XA_LOOKBACK dict
    if (CTX(OPTION_XA_Z)->b250.len && !CTX(OPTION_XA_Z)->local.len)
        ctx_create_node (vb, OPTION_XA_LOOKBACK, "0", 1);

    // Primary VB - ingest VB data into z_file->sa_*
    if (sam_is_prim_vb)
        sam_zip_prim_ingest_vb (VB_SAM);

    // collect stats
    if (!segconf.running) {
        __atomic_add_fetch (&z_file->mate_line_count, (uint64_t)VB_SAM->mate_line_count, __ATOMIC_RELAXED);
        __atomic_add_fetch (&z_file->prim_near_count, (uint64_t)VB_SAM->prim_near_count, __ATOMIC_RELAXED);
        __atomic_add_fetch (&z_file->prim_far_count,  (uint64_t)VB_SAM->prim_far_count,  __ATOMIC_RELAXED);
    }
}

// main thread: called after all VBs, before compressing global sections
void sam_zip_after_vbs (void)
{
    // shorten unused words in dictionary strings to "" (dict pre-populated in sam_zip_initialize)
    if (flag.reference != REF_INTERNAL)  // TO DO: this doesn't work for REF_INTERNAL for example with test.transcriptome.bam
        ctx_shorten_unused_dict_words (SAM_RNAME);

    ctx_shorten_unused_dict_words (SAM_RNEXT);
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
        dict_id.num == _OPTION_XE_i;
}

static void sam_get_one_aux (VBlockSAMP vb, STRp(aux),
                             rom *tag, char *type, char *array_subtype, pSTRp(value)) // out
{
    ASSSEG (aux_len >= 6 && aux[2] == ':' && aux[4] == ':', aux, "invalid optional field format: %.*s", STRf(aux));

    *tag         = aux;
    *type        = aux[3];
    
    if (*type == 'B') {
        *array_subtype = aux[5];
        *value         = aux + 7;
        *value_len     = aux_len - 7;
    }
    else {
        *array_subtype = 0;
        *value         = aux + 5;
        *value_len     = aux_len - 5;
    }

    ASSSEG0 (*value_len < 10000000, aux, "Invalid array field format"); // check that *value_len didn't underflow beneath 0
}

static int64_t sam_get_one_aux_int (VBlockSAMP vb, STRp(aux))
{
    rom tag                  __attribute__((unused));
    char type, array_subtype __attribute__((unused));
    STR (value_str);

    sam_get_one_aux (vb, STRa(aux), &tag, &type, &array_subtype, pSTRa(value_str));

    int64_t value;    
    ASSSEG (str_get_int (STRa(value_str), &value), aux, "%.*s is not an integer", STRf(aux));

    return value;
}


void sam_seg_idx_aux (VBlockSAMP vb, STRps(aux))
{
    if ((int32_t)n_auxs < 1) return; // this line has no AUX fields (possibly negative)

    bool is_bam = IS_BAM_ZIP;

    for (int16_t f=0; f < n_auxs; f++) {
        ASSSEG0 (aux_lens[f] > (is_bam ? 3 : 5), (is_bam ? NULL : auxs[f]), "Invalid auxilliary field format");

        char c1 = auxs[f][0];
        char c2 = auxs[f][1];
        char c3 = is_bam ? sam_seg_bam_type_to_sam_type(auxs[f][2]) : auxs[f][3];

        #define TEST_AUX(name,f1,f2,f3) if (c1==f1 && c2==f2 && c3==f3) vb->idx_##name = f

             TEST_AUX (NM_i, 'N','M','i');
        else TEST_AUX (NH_i, 'N','H','i');
        else TEST_AUX (MD_Z, 'M','D','Z');
        else TEST_AUX (SA_Z, 'S','A','Z');
        else TEST_AUX (HI_i, 'H','I','i');
        else TEST_AUX (XG_Z, 'X','G','Z');
        else TEST_AUX (X0_i, 'X','0','i');
        else TEST_AUX (X1_i, 'X','1','i');
        else TEST_AUX (XA_Z, 'X','A','Z');
        else TEST_AUX (AS_i, 'A','S','i');
        else TEST_AUX (CC_Z, 'C','C','Z');
        else TEST_AUX (CP_i, 'C','P','i');
        else TEST_AUX (ms_i, 'm','s','i');
    }
}

// returns aux field if it exists or NULL if it doesn't
// In BAM, with a numeric field, result is returned in numeric, otherwise in the return value + value_len
static rom sam_seg_get_aux (VBlockSAMP vb, STRp (aux), uint32_t *value_len, ValueType *numeric, bool is_bam)
{
    STR0 (my_value);
    rom tag;
    char sam_type=0, bam_type=0, array_subtype=0;

    if (is_bam) bam_get_one_aux (vb, STRa(aux), &tag, &bam_type, &array_subtype, pSTRa(my_value), numeric);
    else        sam_get_one_aux (vb, STRa(aux), &tag, &sam_type, &array_subtype, pSTRa(my_value));

    if (value_len) *value_len = my_value_len;
    return (is_bam && !my_value) ? &aux[3] : my_value; // in BAM, if NULL, value is in numeric
}

// returns the length of the field, or 0 if it is not found
// note: BAM allows aux ints up to 0xffffffff, but this function supports only up to 0x7fffffff. Don't
// use for fields that might have a larger value
uint32_t sam_seg_get_aux_int (VBlockSAMP vb, STRp (aux), 
                              int32_t *number, // modified only if integer is parsed
                              bool is_bam,
                              int32_t min_value, int32_t max_value, bool soft_fail)
{
    ValueType numeric; 
    uint32_t value_len;
    rom value = sam_seg_get_aux (vb, STRa(aux), &value_len, &numeric, is_bam);
  
    if (value && is_bam) {
        if (numeric.i < min_value || numeric.i > max_value) goto out_of_range;
        *number = numeric.i;
        return value_len;
    }
  
    else if (value && str_get_int (STRa(value), &numeric.i)) {
        if (numeric.i < min_value || numeric.i > max_value) goto out_of_range;
        *number = numeric.i;
        return value_len;
    }
  
    // this line doesn't have this field, or (SAM only) the field is not a valid integer
    if (soft_fail) return 0;
    ABORT_R ("%s: no valid value found for %.2s", LN_NAME, aux);

out_of_range:
    if (soft_fail) return 0;
    ABORT_R ("%s: value of %.2s=%"PRId64" is out of range [%d,%d]", LN_NAME, aux, numeric.i, min_value, max_value);
}

void sam_seg_get_aux_str (VBlockSAMP vb, STRp (aux), pSTRp (snip), bool is_bam)
{
    *snip = aux         + (is_bam ? 3 : 5);
    *snip_len = aux_len - (is_bam ? 3 : 5);
}

// note: we test RNAME but not POS. POS exceeding contig LN is handled in ref_seg_get_locked_range_loaded
void sam_seg_verify_RNAME (VBlockSAMP vb, rom p_into_txt)
{
    if (segconf.running) return;

    if (flag.reference == REF_INTERNAL && (!sam_hdr_contigs /* SQ-less SAM */ || !sam_hdr_contigs->contigs.len /* SQ-less BAM */)) return;
    if (vb->chrom_name_len==1 && *vb->chrom_name=='*') return; // unaligned
    
    if (sam_hdr_contigs) 
        // since this SAM file has a header, all RNAMEs must be listed in it (the header contigs appear first in CTX(RNAME), see sam_zip_initialize
        ASSSEG (vb->chrom_node_index < sam_hdr_contigs->contigs.len, p_into_txt, "RNAME \"%.*s\" does not have an SQ record in the header", STRf(vb->chrom_name));

    else { // headerless SAM
        WordIndex ref_index = chrom_2ref_seg_get (gref, VB, vb->chrom_node_index); // possibly an alt contig
        if (ref_index == WORD_INDEX_NONE) {
            WARN_ONCE ("FYI: RNAME \"%.*s\" (and possibly others) is missing in the reference file. This might impact the compression ratio.", 
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
    *tag        = "tx"; // genozip introduced tag (=taxid)
    *sam_type   = 'i';
    *value      = taxid_str;
    *has_kraken = false;
    *value_len  = kraken_seg_taxid_do (VB, SAM_TAXID, last_txt (VB, SAM_QNAME), vb->last_txt_len (SAM_QNAME),
                                       taxid_str, true);
    *numeric = CTX(SAM_TAXID)->last_value;

    vb->recon_size += is_bam ? 7 : (*value_len + 6); // txt modified
}

void sam_seg_aux_all (VBlockSAMP vb, ZipDataLineSAM *dl, STRps(aux))
{
    START_TIMER;

    const bool is_bam = IS_BAM_ZIP;
    Container con = { .repeats=1 };
    char prefixes[MAX_FIELDS * 6 + 3]; // each name is 5 characters per SAM specification, eg "MC:Z:" followed by CON_PX_SEP ; +3 for the initial CON_PX_SEP
    prefixes[0] = prefixes[1] = prefixes[2] = CON_PX_SEP; // initial CON_PX_SEP follow by separator of empty Container-wide prefix followed by separator for empty prefix for translator-only item[0]
    unsigned prefixes_len=3;

    // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
    con.items[con_nitems(con)] = (ContainerItem){ .translator = SAM2BAM_AUX_SELF }; 
    con_inc_nitems (con);
    
    bool has_kraken = kraken_is_loaded;

    for (uint32_t f=0 ; f < n_auxs + has_kraken; f++) {

        STR0(value);
        ValueType numeric = {};
        rom tag;
        char sam_type=0, bam_type=0, array_subtype=0;
        char taxid_str[20];

        if (f == n_auxs) sam_seg_get_kraken (vb, &has_kraken, taxid_str, &tag, &bam_type,   pSTRa(value), &numeric, is_bam);
        else if (is_bam) bam_get_one_aux (vb, STRi(aux,f), &tag, &bam_type, &array_subtype, pSTRa(value), &numeric);
        else             sam_get_one_aux (vb, STRi(aux,f), &tag, &sam_type, &array_subtype, pSTRa(value));

        if (sam_type == 'i') {
            if (f == vb->idx_NM_i) // we already converted NM to integer, no need to do it again
                numeric.i = dl->NM; 
            else 
                ASSERT (str_get_int (STRa(value), &numeric.i), "%s: Expecting integer value for auxiliary field %c%c but found \"%.*s\"",
                        LN_NAME, tag[0], tag[1], STRf(value));
            value=0;
        }
            
        if (!bam_type) bam_type = sam_seg_sam_type_to_bam_type (sam_type, numeric.i);
        if (!sam_type) sam_type = sam_seg_bam_type_to_sam_type (bam_type);
        
        ASSERT (bam_type, "%s: value %"PRId64" of field %c%c is out of range of the BAM specification: [%d-%u]", 
                LN_NAME, numeric.i, tag[0], tag[1], -0x80000000, 0x7fffffff);

        con.items[con_nitems(con)] = (ContainerItem) {
            .dict_id    = sam_seg_aux_field (vb, dl, is_bam, tag, bam_type, array_subtype, STRa(value), numeric),
            .translator = aux_field_translator ((uint8_t)bam_type), // how to transform the field if reconstructing to BAM
            .separator  = { aux_sep_by_type[is_bam][(uint8_t)bam_type], '\t' },
        };

        con_inc_nitems (con);

        ASSSEG (con_nitems(con) <= MAX_FIELDS, value, "too many optional fields, limit is %u", MAX_FIELDS);

        // note: in the optional field prefix (unlike array type), all integer types become 'i'.
        char prefix[6] = { tag[0], tag[1], ':', sam_type, ':', CON_PX_SEP}; 
        memcpy (&prefixes[prefixes_len], prefix, 6);
        prefixes_len += 6;
    }

    uint32_t num_items = con_nitems(con);
    if (num_items > 1) { // we have Aux fields, not just the translator item
        if (con.items[num_items-1].separator[0] & 0x80) // is a flag
            con.items[num_items-1].separator[0] &= ~(CI0_NATIVE_NEXT & ~(uint8_t)0x80); // last Aux field has no tab
        con.items[num_items-1].separator[1] = 0;
        container_seg (vb, CTX(SAM_AUX), &con, prefixes, prefixes_len, (is_bam ? 3 : 5) * (num_items-1)); // account for : SAM: "MX:i:" BAM: "MXi"
    }
    else
        // NULL means MISSING Container item (of the toplevel container) - will cause container_reconstruct_do of 
        // the toplevel container to delete of previous separator (\t)
        container_seg (vb, CTX(SAM_AUX), 0, 0, 0, 0); 

    COPY_TIMER (sam_seg_aux_all);
}

static inline bool has_same_qname (VBlockSAMP vb, STRp(qname), LineIType buddy_line_i)
{
    return buddy_line_i != NO_LINE &&  
           ({ ZipDataLineSAM *buddy_dl = DATA_LINE (buddy_line_i); 
              buddy_dl->QNAME.len == qname_len && !memcmp (qname, Btxt (buddy_dl->QNAME.index), qname_len); });
}

// seg BUDDY and return true, if this QNAME (and other fields of this alignment) is to be compressed against a buddy 
static inline bool sam_seg_buddy (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qname))
{
    LineIType buddy_line_i = NO_LINE, candidate;

    uint32_t my_hash = QNAME_HASH (qname, qname_len,  dl->FLAG.is_last) & MAXB(vb->qname_hash.prm8[0]);
    #define mate_hash  QNAME_HASH (qname, qname_len, !dl->FLAG.is_last) & MAXB(vb->qname_hash.prm8[0]) // calculated only if needed
    
    #define LINE_BY_HASH(hash) *B(LineIType, vb->qname_hash, (hash))

    // note there are cases where a line can both have mate and be depn. if a line "can_be_mate" we don't
    // allow segging as prim, even if there isn't a mate. this is because of muxes that predict the channel_i
    // based on can_mate(?)
    // 1) R1 Prim  can_mate, can_be_depn, can_be_prin : not mate, not depn because no matching previous QNAME
    // 2) R1 Depn  can_be_depn                        : prim_line because (1) is a matching QNAME
    // 3) R2 Depn  can depn                           : not mate (bc Depn) and not depn (no matching QNAME)
    // 3) R2 Prim  can_mate, can_be_depn, can_be_prin : prim_line (buddy_line=3) and also mate_line (buddy_line=1). prim_line takes precedece.

    bool line_can_mate    = !sam_line_is_depn (dl) && segconf.is_paired; // only in MAIN and PRIM VBs
    bool line_can_be_depn = !line_can_mate && sam_is_main_vb && !vb->check_for_gc; // only in MAIN VB
    bool line_can_be_prim = sam_is_main_vb && !vb->check_for_gc && !dl->has_hard_clips; // only in MAIN VB

    // bool line_can_be_depn = sam_is_main_vb   && !vb->check_for_gc;           // only in MAIN VB
    // bool line_can_be_prim = line_can_be_depn && !dl->has_hard_clips;         // only in MAIN VB
    // bool line_can_mate    = !sam_line_is_depn (dl) && segconf.is_paired; // only in MAIN and PRIM VBs

    // case: depn line in main VB - seg against "primary" (if we in cases that there cannot be a mate)
    // note: the only requirement for a "primary" line here is that it has no hard clips. no consideration of SUPP/SECONARY flags.
    if (line_can_be_depn && has_same_qname (vb, STRa(qname), (candidate = LINE_BY_HASH (my_hash))) &&
             !DATA_LINE(candidate)->has_hard_clips) { // the "prim" line against which we are segging cannot have hard clips
        vb->prim_line_i = buddy_line_i = candidate; 
        vb->prim_near_count++;
    }

    // case: buddy is a mate. note relevant for dependent lines. 
    else if (line_can_mate && has_same_qname (vb, STRa(qname), (candidate = LINE_BY_HASH (mate_hash))) &&
             sam_line_is_prim(DATA_LINE(candidate))) {
        vb->mate_line_i = buddy_line_i = candidate;
        vb->mate_line_count++; // for stats
    }

    // store this line if there's a chance we will seg against it (we don't hash unneeded lines to reduce hash contention)
    // Note: In case of multiple consecutive depns in collated file, we overwrite the hash with each line. This keeps BUDDY=1 and hence more compressible.
    if (line_can_mate || line_can_be_prim)
        LINE_BY_HASH (my_hash) = vb->line_i;

    // case: if QNAME is the same as proposed buddy's - seg against buddy
    if (buddy_line_i != NO_LINE) {
        seg_add_to_local_resizable (VB, CTX(SAM_BUDDY), vb->line_i - buddy_line_i, 0); // add buddy (delta) >= 1 . 
        return true;
    }
    else
        return false;
}

static inline void sam_seg_QNAME_segconf (VBlockSAMP vb, ContextP ctx, STRp(qname))
{
    if (segconf.is_collated) { // still collated, try to find evidence to the contrary
        bool is_new = (qname_len != ctx->last_txt.len || memcmp (qname, last_txtx(vb, ctx), qname_len));
        if (is_new && ctx->last_is_new)
            segconf.is_collated = false; // two new QNAMEs in a row = not collated 
        
        // case: at least on pair of consecutive lines has the same QNAME. if the file is not sorted, we will set it as collated in sam_seg_finalize_segconf
        if (!is_new)
            segconf.evidence_of_collated = true; 

        ctx->last_is_new = is_new;
    }
    
    if (vb->line_i==0)
        qname_segconf_discover_flavor (VB, SAM_QNAME, STRa(qname));
}

void sam_seg_QNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(qname), unsigned add_additional_bytes)
{
    ContextP ctx = CTX(SAM_QNAME);

    // in segconf, identify if this file is collated (each QNAME appears in two or more consecutive lines)
    if (segconf.running) {
        sam_seg_QNAME_segconf (vb, ctx, STRa(qname));
        goto normal_seg;
    }

    // case: if QNAME is the same as buddy's - seg against buddy
    if (sam_seg_buddy (vb, dl, STRa(qname))) 
        // case: Seg against buddy - the second SNIP_COPY_BUDDY means that reconstruction of this snip will also set buddy
        // note: in PRIM segging against mate here is used for loading SA Groups, SAM_QNAMESA is used for reconstruction 
        seg_by_ctx (VB, (char[]){ SNIP_COPY_BUDDY, SNIP_COPY_BUDDY }, 2, ctx, qname_len + add_additional_bytes); // seg QNAME as copy-from-buddy (an extra SNIP_COPY_BUDDY indicates that reconstruct_from_buddy should set buddy_line_i here)

    // case: DEPN with SA Group: seg against SA Group (unless already segged against buddy)
    else if (sam_is_depn_vb && vb->sag) 
        sam_seg_against_sa_group (vb, ctx, qname_len + add_additional_bytes);

    else normal_seg:    
        qname_seg (VB, ctx, STRa(qname), add_additional_bytes); // note: for PRIM component, this will be consumed with loading SA

    // case: PRIM: additional seg against SA Group - store in SAM_QNAMESA - Reconstruct will take from here in PRIM per Toplevel container
    if (sam_is_prim_vb) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PRIM_QNAME }, 2, SAM_QNAMESA, 0); // consumed when reconstructing PRIM vb
}

WordIndex sam_seg_RNAME (VBlockSAMP vb, ZipDataLineSAM *dl, STRp (chrom), 
                         bool against_sa_group_ok, // if true, vb->chrom_node_index must already be set
                         unsigned add_bytes)
{
    bool normal_seg = false;

    if (segconf.running) goto normal_seg;

    WordIndex node_index;

    // case: PRIM or DEPN vb - seg against SA group with alignments
    // Note: in DEPN, rname already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    if (against_sa_group_ok && sam_seg_has_sag_by_SA (vb)) {
        sam_seg_against_sa_group (vb, CTX(SAM_RNAME), add_bytes);

        if (sam_is_prim_vb) {
            // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
            seg_by_did (VB, STRa(chrom), OPTION_SA_RNAME, 0);

            // count RNAME field contribution to OPTION_SA_RNAME, so sam_stats_reallocate can allocate the z_data between RNAME and SA:Z
            CTX(OPTION_SA_RNAME)->counts.count += add_bytes; 
        }

        STRset (vb->chrom_name, chrom);
        random_access_update_chrom (VB, 0, vb->chrom_node_index, STRa(chrom)); 

        node_index = vb->chrom_node_index;
    }

    // seg RNAME against mate's RNEXT. limit to only if RNAME/RNEXT are different (if the same, then 
    // likely all the lines between this line and the mate have the same RNAME) 
    // note: for now, this only works in BAM because dl->RNAME/RNEXT are set already. To do: support in SAM.
    else if (IS_BAM_ZIP && zip_has_mate &&
             dl->RNAME != dl->RNEXT && dl->RNAME == DATA_LINE(vb->mate_line_i)->RNEXT) {
        
        seg_by_did (VB, STRa(copy_buddy_RNEXT_snip), SAM_RNAME, add_bytes); // copy POS from earlier-line mate PNEXT
        STRset (vb->chrom_name, chrom);
        random_access_update_chrom (VB, 0, dl->RNAME, STRa(chrom)); 
        node_index = dl->RNAME;
    }

    else normal_seg: {
        bool is_new;
        node_index = chrom_seg_ex (VB, SAM_RNAME, STRa(chrom), 0, NULL, add_bytes, !IS_BAM_ZIP, &is_new);
        normal_seg = true;

        // don't allow adding chroms to a BAM file or a SAM that has SQ lines in the header, but we do allow to add to a headerless SAM.
        ASSSEG (!is_new || !sam_hdr_contigs || segconf.running || (chrom_len==1 && (*chrom=='*' || *chrom=='=')),
                chrom, "contig '%.*s' appears in file, but is missing in the %s header", STRf(chrom), dt_name (vb->data_type));
    }

    // protect rname from removal by ctx_shorten_unused_dict_words if we didn't seg normally
    if (!normal_seg && node_index != NODE_INDEX_NONE) // note: node_index==-1 when RNAME="*"
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
    if (!segconf.has[SAM_RNEXT]) 
        goto normal_seg;

    // case: seg against mate's RNAME. limit to only if RNAME/RNEXT are different (if the same, then 
    // likely all the lines between this line and the mate have the same RNAME) 
    // note: for now, this only works in BAM because dl->RNAME/RNEXT are set already. To do: support in SAM.
    else if (IS_BAM_ZIP && zip_has_mate &&
             dl->RNAME != dl->RNEXT && dl->RNEXT == DATA_LINE(vb->mate_line_i)->RNAME) 
        seg_by_did (VB, STRa(copy_buddy_RNAME_snip), SAM_RNEXT, add_bytes); // copy POS from earlier-line mate PNEXT

    // case: seg RNEXT against prim's RNAME. This happens when RNEXT is the same as prim's RNEXT, but prim's
    // RNEXT is the same as prim's RNAME, so PRIM's RNEXT is segged as "=".
    else if (IS_BAM_ZIP && zip_has_prim &&
             dl->RNAME != dl->RNEXT && dl->RNEXT == DATA_LINE(vb->prim_line_i)->RNAME) 
        seg_by_did (VB, STRa(copy_buddy_RNAME_snip), SAM_RNEXT, add_bytes);

    // case: seg RNEXT against prim's RNEXT
    else if (IS_BAM_ZIP && zip_has_prim &&
             dl->RNAME != dl->RNEXT && dl->RNEXT == DATA_LINE(vb->prim_line_i)->RNEXT) 
        seg_by_did (VB, (char[]){ SNIP_COPY_BUDDY }, 1, SAM_RNEXT, add_bytes);

    else normal_seg: {
        bool is_new;
        node_index = chrom_seg_ex (VB, SAM_RNEXT, STRa(chrom), 0, NULL, add_bytes, !IS_BAM_ZIP, &is_new);
        normal_seg = true;

        // don't allow adding chroms to a BAM file or a SAM that has SQ lines in the header, but we do allow to add to a headerless SAM.
        ASSSEG (!is_new || !sam_hdr_contigs || segconf.running || (chrom_len==1 && (*chrom=='*' || *chrom=='=')),
                chrom, "contig '%.*s' appears in file, but is missing in the %s header", STRf(chrom), dt_name (vb->data_type));
    }

    // protect rnext from removal by ctx_shorten_unused_dict_words if we didn't seg normally
    if (!normal_seg && node_index != NODE_INDEX_NONE) // note: node_index==-1 when RNEXT="*"
        ctx_protect_from_removal (VB, CTX(SAM_RNEXT), node_index); 

    return node_index;
}

bool sam_seg_test_biopsy_line (VBlockP vb, STRp(line))
{
    if (segconf.running) return false; // we need to let segconf run normally, so we get the correct VB size

    if (flag.biopsy_line.line_i == vb->line_i && flag.biopsy_line.vb_i == vb->vblock_i) {
        PutLineFn fn = file_put_line (VB, STRa(line), "Line biopsy:");
        
        if (TXT_DT(DT_BAM)) WARN ("Tip: You can view the dumped BAM line with:\n   genozip --show-bam %s", fn.s);
        exit_ok();
    }
    vb->recon_size -= line_len;
    return true;
}

rom sam_seg_txt_line (VBlockP vb_, rom next_line, uint32_t remaining_txt_len, bool *has_13)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    sam_reset_line (VB);

    WordIndex prev_line_chrom = vb->chrom_node_index;

    ZipDataLineSAM *dl = DATA_LINE (vb->line_i);
        
    str_split_by_tab (next_line, remaining_txt_len, MAX_FIELDS+AUX, has_13); // also advances next_line to next line
    rom *auxs          = &flds[AUX];
    uint32_t *aux_lens = &fld_lens[AUX];
    uint32_t n_auxs    = n_flds - AUX;

    sam_seg_idx_aux (vb, STRas(aux));

    dl->SEQ.index = BNUMtxt (flds[SEQ]); // SEQ.len to be determined by sam_cigar_analyze
    dl->QNAME = TXTWORDi(fld, QNAME);
    dl->QUAL  = TXTWORDi(fld, QUAL);
    
    if (fld_lens[QUAL]==1 && flds[QUAL][0]=='*') 
        vb->qual_missing = dl->no_qual = true;

    // lazy way to get vb->chrom* (inc. if --match-chrom), rollback later if seg RNAME is not needed
    if (vb->check_for_gc || !sam_is_main_vb)
        seg_create_rollback_point (VB, NULL, 1, SAM_RNAME);

    dl->RNAME = sam_seg_RNAME (vb, dl, STRfld(RNAME), false, fld_lens[RNAME]+1);
 
    ASSSEG (str_get_int_range32 (STRfld(POS), 0, MAX_POS_SAM, &dl->POS), 
            flds[POS], "Invalid POS \"%.*s\": expecting an integer [0,%d]", STRfi(fld,POS), (int)MAX_POS_SAM);

    ASSSEG (str_get_int_range8 (STRfld(MAPQ), 0, 255, &dl->MAPQ), 
            flds[MAPQ], "Invalid MAPQ \"%.*s\": expecting an integer [0,255]", STRfi(fld,MAPQ));

    ASSSEG (str_get_int_range16 (STRfld(FLAG), 0, SAM_MAX_FLAG, &dl->FLAG.value), flds[FLAG], "invalid FLAG field: \"%.*s\"", STRfi(fld,FLAG));

    // analyze (but don't seg yet) CIGAR
    uint32_t seq_len;
    sam_cigar_analyze (vb, STRfld(CIGAR), false, &seq_len);
    dl->SEQ.len = seq_len; // do it this way to avoid compiler warning

    ASSSEG (dl->SEQ.len == fld_lens[SEQ] || (fld_lens[CIGAR] == 1 && *flds[CIGAR] == '*') || flds[SEQ][0] == '*', flds[SEQ], 
            "seq_len implied by CIGAR=%.*s is %u, but actual SEQ length is %u, SEQ=%.*s", 
            STRfi(fld,CIGAR), dl->SEQ.len, fld_lens[SEQ], STRfi(fld,SEQ));

    // if this is a secondary / supplamentary read (aka Dependent) or a read that has an associated sec/sup read (aka Primary) - move
    // the line to the appropriate component and skip it here (no segging done yet)
    if (vb->check_for_gc && sam_seg_is_gc_line (vb, dl, flds[0], next_line - flds[0], STRasi(fld,AUX), false)) {
        vb->debug_line_hash_skip = true;
        goto rollback_and_done;
    }

    // case: to biopsy_line: we just needed to pass sam_seg_is_gc_line and we're done
    if (flag.biopsy_line.line_i != NO_LINE && sam_seg_test_biopsy_line (VB, flds[0], next_line - flds[0])) 
        goto rollback_and_done;  

    vb->last_cigar = flds[CIGAR];
    SAFE_NUL (&vb->last_cigar[fld_lens[CIGAR]]); // nul-terminate CIGAR string

    if (has_NM)
        dl->NM_len = sam_seg_get_aux_int (vb, STRi(aux, vb->idx_NM_i), &dl->NM, false, MIN_NM_i, MAX_NM_i, false) + 1; // +1 for \t or \n

    if (!sam_is_main_vb) {

        // set dl->AS needed by sam_seg_prim_add_sa_group
        if (sam_is_prim_vb && has_AS)
            sam_seg_get_aux_int (vb, STRi(aux, vb->idx_AS_i), &dl->AS, false, MIN_AS_i, MAX_AS_i, false);

        sam_seg_sa_group_stuff (vb, dl , STRasi(fld,AUX), STRfld(CIGAR), flds[SEQ], false);

        // re-seg rname, against SA group
        seg_rollback (VB);
        sam_seg_RNAME (vb, dl, STRfld(RNAME), true, fld_lens[RNAME]+1);
    }

    sam_seg_QNAME (vb, dl, STRfld(QNAME), 1);

    sam_seg_FLAG (vb, dl, fld_lens[FLAG]+1);

    // note: pos can have a value even if RNAME="*" - this happens if a SAM with a RNAME that is not in the header is converted to BAM with samtools
    sam_seg_POS (vb, dl, prev_line_chrom, fld_lens[POS]+1);

    if (fld_lens[RNAME] != 1 || *flds[RNAME] != '*')
        sam_seg_verify_RNAME (vb, flds[RNAME]);

    sam_seg_MAPQ (vb, dl, fld_lens[MAPQ]+1);

    dl->RNEXT = sam_seg_RNEXT (vb, dl, STRfld(RNEXT), fld_lens[RNEXT]+1);

    vb->RNEXT_is_equal = ((fld_lens[RNEXT]==1 && *flds[RNEXT] == '=') || str_issame_(STRfld(RNAME), STRfld(RNEXT)));

    sam_seg_PNEXT (vb, dl, STRfld(PNEXT), 0, fld_lens[PNEXT]+1);

    // we search forward for MD:Z now, XG:Z as we will need it for SEQ if it exists
    if (has_MD)
        sam_seg_MD_Z_analyze (vb, dl, STRauxZ(MD_Z,false), dl->POS);

    if (has_XG)
        sam_seg_bsseeker2_XG_Z_analyze (vb, dl, STRauxZ(XG_Z,false), dl->POS);

    // analyzing X0 - needed for segging XT:A - set to 1 X0==1
    if (has_X0)
        ctx_set_last_value (VB, CTX(OPTION_X0_i), sam_get_one_aux_int (vb, STRi(fld, AUX+vb->idx_X0_i))); 

    // analyzing XA - needed for segging X1:i - set to number of repeats
    if (has_XA && has_X1)
        ctx_set_last_value (VB, CTX(OPTION_XA_Z), (int64_t)str_count_char (flds[AUX+vb->idx_XA_Z]+5, fld_lens[AUX+vb->idx_XA_Z]-5, ';')); 

    if (has_ms)
        ctx_set_last_value (VB, CTX(OPTION_s1_i), sam_get_one_aux_int (vb, STRi(fld, AUX+vb->idx_ms_i))); 

    seg_set_last_txt (VB, CTX(SAM_SQBITMAP), STRfld(SEQ));

    // calculate diff vs. reference (denovo or loaded)
    sam_seg_SEQ (vb, dl, STRfld(SEQ), fld_lens[SEQ]+1);
    
    sam_seg_QUAL (vb, dl, STRfld(QUAL), fld_lens[QUAL] + 1); 
    
    ASSSEG (str_is_in_range (flds[QUAL], fld_lens[QUAL], 33, 126), flds[QUAL], "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", 
            STRfi(fld,QUAL));

    if (fld_lens[SEQ] != fld_lens[QUAL])
        vb->qual_codec_no_longr = true; // we cannot compress QUAL with CODEC_LONGR in this case

    // finally we can seg CIGAR now
    sam_cigar_seg_textual (vb, dl, fld_lens[CIGAR], STRfld(SEQ), STRfld(QUAL));
    
    // add BIN so this file can be reconstructed as BAM
    bam_seg_BIN (vb, dl, 0, false);

    // AUX fields - up to MAX_FIELDS of them
    sam_seg_aux_all (vb, dl, n_flds-AUX, &flds[AUX], &fld_lens[AUX]);

    if (sam_is_prim_vb) {
        if (IS_SAG_FLAG)
            sam_seg_prim_add_sa_group (vb, dl, 0, false);
        
        else if (IS_SAG_SOLO) 
            sam_seg_prim_add_sa_group_SOLO (vb, dl);
    }
    
    // TLEN - must be after AUX as we might need data from MC:Z
    bool is_rname_rnext_same = (fld_lens[RNEXT]==1 && *flds[RNEXT]=='=') || 
                               (fld_lens[RNEXT]==fld_lens[RNAME] && !memcmp (flds[RNEXT], flds[RNAME], fld_lens[RNAME]));

    sam_seg_TLEN (vb, dl, STRfld(TLEN), 0, is_rname_rnext_same);

    SEG_EOL (SAM_EOL, false); /* last field accounted for \n */
    SAFE_RESTORE; // restore \t after CIGAR

    return next_line;

rollback_and_done:
    seg_rollback (VB); // cancelling segging of RNAME
    return next_line;
}
