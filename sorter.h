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
extern void sorter_piz_start_writing (uint32_t component_i);
extern void sorter_piz_handover_data (VBlockP vb);
extern void sorter_piz_finish_writing (void);
extern unsigned sorter_piz_get_pair (uint32_t vb_i, uint32_t *pair_vb_i);
extern void sorter_remove_liftover_rejects (void);        // called if no --luft
extern void sorter_move_liftover_rejects_to_front (void); // called if --luft

#endif

