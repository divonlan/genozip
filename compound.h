// ------------------------------------------------------------------
//   compound.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef COMPOUND_INCLUDED
#define COMPOUND_INCLUDED

#include "genozip.h"

extern const char sep_with_space[], sep_without_space[];
extern void compound_seg (VBlockP vb, ContextP field_ctx, const char *field, unsigned field_len, 
                          const char *is_sep, unsigned nonoptimized_len, unsigned add_additional_bytes);

extern void compound_segconf_test (STRp(qname));
extern void compound_zip_initialize (DictId qname_dict_id);
extern void compound_seg_initialize (VBlockP vb, DidIType qname_did_i);


#endif
