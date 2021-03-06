// ------------------------------------------------------------------
//   coverage.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"

// if updating, update cvr_names in coverage.c too
typedef enum { CVR_SOFT_CLIP, CVR_UNMAPPED, CVR_SECONDARY, CVR_FAILED, CVR_DUPLICATE, CVR_CONTIGS, CVR_TOTAL, NUM_COVER_TYPES } CoverTypes;

extern void coverage_initialize (VBlockP vb);
extern void coverage_add_one_vb (VBlockP vb);
extern void coverage_show_coverage (void);
extern void coverage_show_idxstats (void);
extern void coverage_sex_classifier (bool is_first_z_file);
