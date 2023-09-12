// ------------------------------------------------------------------
//   qname.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"
#include "segconf.h"

#define dict_id_is_qname_sf dict_id_is_type_1
#define dict_id_qname_sf dict_id_type_1

extern const char sep_with_space[], sep_without_space[];
extern void qname_seg (VBlockP vb, QType q, STRp(qname), unsigned add_additional_bytes);
extern void qname_segconf_discover_flavor (VBlockP vb, QType q, STRp(qname));

extern void qname_zip_initialize (void);
extern void qname_seg_initialize (VBlockP vb, QType q, Did st_did_i);

typedef enum  { QTR_SUCCESS, QTR_QNAME_LEN_0, QTR_FIXED_LEN_MISMATCH, QTR_WRONG_Q, QTR_CONTAINER_MISMATCH, QTR_BAD_INTEGER, QTR_BAD_CHARS, QTR_BAD_NUMERIC, QTR_BAD_HEX, QTR_TECH_MISMATCH, QTR_NOT_BARCODE, NUM_QTRs } QnameTestResult;
#define QTR_NAME { "SUCCESS",   "QNAME_LEN=0",   "FIXED_LEN_MISMATCH",   "WRONG_Q",   "CONTAINER_MISMATCH",   "BAD_INTEGER",   "BAD_CHARS",   "BAD_NUMERIC",   "BAD_HEX",   "TECH_MISMATCH",   "NOT_BARCODE" }
extern QnameTestResult qname_test_flavor (STRp(qname), QType q, QnameFlavor qf, bool quiet);

extern uint32_t qname_calc_hash (QType q, STRp(qname), thool is_last, bool canonical, uint32_t *uncanonical_suffix_len);
extern void qname_canonize (QType q, rom qname, uint32_t *qname_len);

extern rom segconf_qf_name (QType q);
extern QnameFlavorId segconf_qf_id (QType q);
extern rom qtype_name (QType q);

typedef void (*QnameSegCallback) (VBlockP vb, ContextP ctx, STRp(value));

// flavor-specific callbacks
extern void sultima_Q5NAME_cb (VBlockP vb, ContextP ctx, STRp(value));

