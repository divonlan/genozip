// ------------------------------------------------------------------
//   regions.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REGIONS_INCLUDED
#define REGIONS_INCLUDED

#include "genozip.h"

extern void regions_add (const char *reg_str);

extern void regions_get_chrom_index(void);

extern bool regions_is_region_included_in_requested_regions (uint32_t chrom_node_index, uint32_t min_pos, uint32_t max_pos);

#endif
