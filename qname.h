// ------------------------------------------------------------------
//   qname.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

#define dict_id_is_qname_sf dict_id_is_type_1
#define dict_id_qname_sf dict_id_type_1

extern const char sep_with_space[], sep_without_space[];
extern void qname_seg (VBlockP vb, ContextP qname_ctx, STRp(qname), unsigned add_additional_bytes);
extern bool qname_is_flavor (STRp(qname), QnameFlavor qf, pSTRp (qname2));

extern void qname_segconf_discover_flavor (VBlockP vb, DidIType qname_did_i, STRp(qname));
extern void qname_zip_initialize (DidIType qname_did_i);
extern void qname_seg_initialize (VBlockP vb, DidIType qname_did_i);

extern const char *qf_name (QnameFlavor qf);
