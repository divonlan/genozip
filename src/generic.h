// ------------------------------------------------------------------
//   generic.h
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#pragma GENDICT_PREFIX GNRIC
#pragma GENDICT GNRIC_DATA=DTYPE_FIELD=DATA
#pragma GENDICT GNRIC_TOPLEVEL=DTYPE_FIELD=TOPLEVEL

extern int32_t generic_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern int32_t generic_is_header_done (bool is_eof);
extern void generic_seg_initialize (VBlockP vb);
extern void generic_seg_finalize (VBlockP vb);
extern bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern rom generic_get_magic (void);
extern rom generic_get_ext (void);

#define GNRIC_SPECIAL { generic_piz_TOPLEVEL }
SPECIAL (GNRIC, 0, TOPLEVEL, generic_piz_TOPLEVEL);
#define NUM_GNRIC_SPECIAL 1