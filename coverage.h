// ------------------------------------------------------------------
//   coverage.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted,
//   under penalties specified in the license.

#include "genozip.h"

// if updating, update cvr_names in coverage.c too
typedef enum { CVR_OTHER_CONTIGS, CVR_ALL_CONTIGS, CVR_SOFT_CLIP, CVR_UNMAPPED, CVR_SECONDARY, CVR_SUPPLEMENTARY, CVR_FAILED, CVR_DUPLICATE, CVR_TOTAL, NUM_COVER_TYPES } CoverTypes;

extern void coverage_initialize (VBlockP vb);
extern void coverage_add_one_vb (VBlockP vb);
extern void coverage_show_coverage (void);
extern void coverage_show_idxstats (void);
extern void coverage_sex_classifier (bool is_first_z_file);
