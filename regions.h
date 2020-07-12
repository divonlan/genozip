// ------------------------------------------------------------------
//   regions.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REGIONS_INCLUDED
#define REGIONS_INCLUDED

#include "genozip.h"

extern void regions_add (const char *reg_str);
extern void regions_make_chregs (void);
extern void regions_transform_negative_to_positive_complement(void);
extern bool regions_get_ra_intersection (uint32_t chrom_node_index, int64_t min_pos, int64_t max_pos, char *intersection_one_ra);
extern unsigned regions_max_num_chregs(void);
extern void regions_display(const char *title);
extern bool regions_is_site_included (uint32_t chrom_word_index, int64_t pos);
extern bool regions_is_range_included (int32_t chrom, int64_t start_pos, int64_t end_pos, bool completely_included);
#define regions_is_ra_included(ra) regions_is_range_included(ra->chrom_index, ra->min_pos, ra->max_pos, false)

#endif
