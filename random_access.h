// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef RANDOM_ACCESS_INCLUDED
#define RANDOM_ACCESS_INCLUDED

#include <stdint.h>
#include "genozip.h"

extern int32_t random_access_get_last_chrom_node_index (VariantBlockP vb);
extern void random_access_new_entry (VariantBlockP vb, uint32_t vb_line_i, int32_t chrom_node_index);
extern void random_access_update_last_entry (VariantBlockP vb, int32_t this_pos);
extern void random_access_merge_in_vb (VariantBlockP vb);
extern void BGEN_random_access (BufferP ra_buf);
extern unsigned random_access_sizeof_entry();
extern void random_access_show (ConstBufferP ra_buf, bool is_big_endian);

#endif