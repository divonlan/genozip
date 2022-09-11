// ------------------------------------------------------------------
//   qname.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#define dict_id_is_qname_sf dict_id_is_type_1
#define dict_id_qname_sf dict_id_type_1

extern const char sep_with_space[], sep_without_space[];
extern void qname_seg (VBlockP vb, ContextP qname_ctx, STRp(qname), unsigned add_additional_bytes);
extern bool qname_seg_qf (VBlockP vb, ContextP qname_ctx, QnameFlavor qfs, STRp(qname), bool use_qname2, unsigned add_additional_bytes);
extern int qname_test_flavor (STRp(qname), QnameFlavor qf, pSTRp (qname2));

extern void qname_segconf_discover_flavor (VBlockP vb, Did qname_did_i, STRp(qname));
extern bool qname_segconf_discover_fastq_line3_sra_flavor (VBlockP vb, STRp(line3));

extern void qname_zip_initialize (Did qname_did_i);
extern void qname_seg_initialize (VBlockP vb, Did qname_did_i);

extern rom qf_name (QnameFlavor qf);
