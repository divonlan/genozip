// ------------------------------------------------------------------
//   qname.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"
#include "segconf.h"
#include "aliases.h"

// NOTE: values are ZIP-only, not part of file format
typedef packed_enum { 
    QF_NO_ID=0, 
    // Illumina flavors
    QF_ILLUM_7, QF_ILLUM_7i, QF_ILLUM_7umi, QF_ILLUM_7_bc, QF_ILLUM_7gs, QF_ILLUM_5i, QF_ILLUM_5, QF_ILLUM_5rng, QF_ILLUM_6,
    QF_ILLUM_X_0bc, QF_ILLUM_X_1bc, QF_ILLUM_X_2bc, QF_ILLUM_S_0bc, QF_ILLUM_S_1bc, QF_ILLUM_S_2bc, QF_ILLUM_7gsFQ, QF_ILLUM_7_2bc,
    QF_ILLUM_7_rbc,
    // Illumina-style FASTQ QNAME2 flavors (also appears in Ultima, Singular...)
    QF_ILLUM_2bc, QF_ILLUM_1bc, QF_ILLUM_0bc, 
    // MGI flavors
    QF_MGI_NEW6, QF_MGI_NEW7, QF_MGI_NEW8, QF_MGI_SAP8, QF_MGI_MFT, QF_MGI_varlen, QF_MGI_r6, QF_MGI_die6, QF_MGI_r7, QF_MGI_r8, QF_MGI_ll7, QF_MGI_cl, QF_MGI_ml, QF_MGI_dl, QF_MGI_rgs8, QF_MGI_rgs8FQ,
    // PacBio flavors 
    QF_PACBIO_3, QF_PACBIO_rng, QF_PACBIO_lbl, QF_PACBIO_pln, QF_ONSO,
    // Nanopore flavors
    QF_NANOPORE, QF_NANOPORE_rng, QF_NANOPORE_ext,
    // Ultima flavors
    QF_ULTIMA_a, QF_ULTIMA_b6_bc, QF_ULTIMA_a_bc, QF_ULTIMA_c, QF_ULTIMA_c_bc, QF_ULTIMA_b6, QF_ULTIMA_d, QF_ULTIMA_d_bc, 
    QF_ULTIMA_b9, QF_ULTIMA_b9_bc, QF_ULTIMA_n, QF_ULTIMA_Z9,
    // Other sequencer flavors
    QF_ION_TORR_3, QF_ROCHE_454, QF_HELICOS, QF_SINGULAR, QF_ELEMENT, QF_ELEMENT_AV, QF_ELEMENT_bc, 
    // NCBI flavors
    QF_SRA_L, QF_SRA2, QF_SRA, QF_SRA_sra, QF_SRA_label,
    // Consensus alignments flavors
    QF_CONSENSUS, QF_XCON, 
    // Other syntheic (i.e. non-sequencer) flavors
    QF_GENOZIP_OPT, QF_INTEGER, QF_HEX_CHR, QF_BAMSURGEON, QF_SEQAN, QF_CLC_GW, QF_STR_INT, QF_Sint, QF_GENERATED,
    NUM_FLAVORS
} QnameFlavorId;

#define FLAVOR(q,flavor) (segconf.qname_flavor[q] && segconf.qname_flavor[q]->id == QF_##flavor)

#define dict_id_is_qname_sf dict_id_is_type_1
#define dict_id_qname_sf dict_id_type_1

extern const char sep_with_space[], sep_without_space[];
extern void qname_seg (VBlockP vb, QType q, STRp(qname), unsigned add_additional_bytes);
extern void qname_segconf_discover_flavor (VBlockP vb, QType q, STRp(qname));
extern bool qname_segconf_rediscover_flavor (VBlockP vb, QType q, STRp(qname));
extern QType qname_sam_get_qtype (STRp(qname));
extern QnameFlavor qname_get_optimize_qf (void);

extern void qname_zip_initialize (void);
extern void qname_seg_initialize (VBlockP vb, QType q, Did st_did_i);
extern void qname_segconf_finalize (VBlockP vb);

typedef enum  { QTR_SUCCESS, QTR_QNAME_LEN_0, QTR_FIXED_LEN_MISMATCH, QTR_WRONG_Q, QTR_CONTAINER_MISMATCH, QTR_BAD_INTEGER, QTR_BAD_CHARS, QTR_BAD_NUMERIC, QTR_BAD_HEX, QTR_TECH_MISMATCH, QTR_NOT_BARCODE, QTR_NOT_BARCODE2, QTR_NO_MATE, QTR_BAD_MATE, QTR_FAILED_VALIDATE_FUNC, NUM_QTRs } QnameTestResult;
#define QTR_NAME { "SUCCESS",   "QNAME_LEN=0",   "FIXED_LEN_MISMATCH",   "WRONG_Q",   "CONTAINER_MISMATCH",   "BAD_INTEGER",   "BAD_CHARS",   "BAD_NUMERIC",   "BAD_HEX",   "TECH_MISMATCH",   "NOT_BARCODE",   "NOT_BARCODE2",   "NO_MATE",   "BAD_MATE",  "FAILED_VALIDATE_FUNC"}
extern QnameTestResult qname_test_flavor (STRp(qname), QType q, QnameFlavor qf, bool quiet);

typedef enum { CRC32, CRC64 } CrcType;
extern uint64_t qname_calc_hash (QType q, CompIType comp_i, STRp(qname), thool is_last, bool canonical, CrcType type, uint32_t *uncanonical_suffix_len);
extern uint32_t qname_hash_change_last (uint32_t hash, bool is_last);

extern void qname_canonize (QType q, rom qname, uint32_t *qname_len, CompIType comp_i);

extern rom segconf_qf_name (QType q);
extern QnameFlavorId segconf_qf_id (QType q);
extern rom qtype_name (QType q);
extern DictIdAlias qname_get_alias (QType q);
extern bool qf_is_mated (QType q);

typedef void (*QnameSegCallback) (VBlockP vb, ContextP ctx, STRp(value));

// flavor-specific callbacks
extern void ultima_c_Q5NAME_cb (VBlockP vb, ContextP ctx, STRp(value));
extern void seg_qname_rng2seq_len_cb (VBlockP vb, ContextP ctx, STRp(value));

