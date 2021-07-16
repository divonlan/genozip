// ------------------------------------------------------------------
//   generic.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef GNRIC_INCLUDED
#define GNRIC_INCLUDED

#include "genozip.h"

extern int32_t generic_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern void generic_seg_finalize (VBlockP vb);
extern bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id);

#define GNRIC_SPECIAL { generic_piz_TOPLEVEL }
SPECIAL (GNRIC, 0, TOPLEVEL, generic_piz_TOPLEVEL);
#define NUM_GNRIC_SPECIAL 1

#endif