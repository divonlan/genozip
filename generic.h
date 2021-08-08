// ------------------------------------------------------------------
//   generic.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef GNRIC_INCLUDED
#define GNRIC_INCLUDED

#include "genozip.h"

#define _GNRIC_DATA     DICT_ID_MAKEF_4 ("DATA")
#define _GNRIC_TOPLEVEL DICT_ID_MAKEF_L (TOPLEVEL)

typedef enum { GNRIC_DATA, GNRIC_TOPLEVEL, NUM_GNRIC_FIELDS } GenericFields;
#define GNRIC_MAPPING { V(GNRIC_DATA), V(GNRIC_TOPLEVEL) }

extern int32_t generic_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern void generic_seg_finalize (VBlockP vb);
extern bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id);

#define GNRIC_SPECIAL { generic_piz_TOPLEVEL }
SPECIAL (GNRIC, 0, TOPLEVEL, generic_piz_TOPLEVEL);
#define NUM_GNRIC_SPECIAL 1

#endif