// ------------------------------------------------------------------
//   regions.h
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REGIONS_INCLUDED
#define REGIONS_INCLUDED

#include "genozip.h"

extern void regions_add (const char *reg_str);
extern void regions_add_by_file (const char *regions_filename);
extern void regions_make_chregs (ContextP chrom_ctx);
extern void regions_transform_negative_to_positive_complement(void);
extern bool regions_get_ra_intersection (WordIndex chrom_node_index, PosType min_pos, PosType max_pos);

extern unsigned regions_get_num_range_intersections (WordIndex chrom_word_index);
extern bool regions_get_range_intersection (WordIndex chrom_word_index, PosType min_pos, PosType max_pos, unsigned intersect_i, PosType *intersect_min_pos, PosType *intersect_max_pos, bool *revcomp);

extern unsigned regions_max_num_chregs(void);
extern void regions_display(const char *title);
extern bool regions_is_site_included (VBlockP vb);
extern bool regions_is_range_included (WordIndex chrom, PosType start_pos, PosType end_pos, bool completely_included);
#define regions_is_ra_included(ra) regions_is_range_included(ra->chrom_index, ra->min_pos, ra->max_pos, false)

#endif
