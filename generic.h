// ------------------------------------------------------------------
//   generic.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GNRIC_INCLUDED
#define GNRIC_INCLUDED

#include "genozip.h"

extern void generic_seg_finalize (VBlockP vb);

#define GNRIC_SPECIAL { generic_piz_TOPLEVEL }
SPECIAL (GNRIC, 0, TOPLEVEL, generic_piz_TOPLEVEL);
#define NUM_GNRIC_SPECIAL 1

#endif