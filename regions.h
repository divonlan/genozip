// ------------------------------------------------------------------
//   regions.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REGIONS_INCLUDED
#define REGIONS_INCLUDED

#include "genozip.h"

extern void regions_add (const char *reg_str);
extern void regions_make_chregs(void);
extern void regions_transform_negative_to_positive_complement(void);
extern bool regions_get_ra_intersection (uint32_t chrom_node_index, uint32_t min_pos, uint32_t max_pos, char *intersection_one_ra);
extern unsigned regions_max_num_chregs(void);
extern void regions_display(const char *title);
extern bool regions_is_site_included (uint32_t chrom_word_index, uint32_t pos);

#endif
