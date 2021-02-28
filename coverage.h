// ------------------------------------------------------------------
//   coverage.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

extern void coverage_initialize (VBlockP vb);
extern void coverage_add_one_vb (VBlockP vb);
extern void coverage_show_coverage (void);
extern void coverage_sex_classifier (void);
