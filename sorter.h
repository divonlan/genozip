// ------------------------------------------------------------------
//   sorter.h
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SORTER_INCLUDED
#define SORTER_INCLUDED

#include "genozip.h"

// ZIP side
extern void sorter_zip_merge_vb (VBlockP vb);
extern void sorter_compress_recon_plan (void);

// PIZ side
extern void sorter_piz_get_reconstruction_plan (uint32_t component_i);

#endif

