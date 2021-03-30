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
extern void sorter_piz_create_plan (void);
extern void sorter_piz_start_writing (uint32_t txt_file_i);
extern void sorter_piz_handover_data (VBlockP vb);
extern void sorter_piz_handover_txtheader (uint32_t component_i);
extern void sorter_piz_finish_writing (void);
extern unsigned sorter_piz_get_pair (uint32_t vb_i, uint32_t *pair_vb_i);

// PIZ - getting and setting VB / txtheader properties
extern bool sorter_is_txtheader_in_plan (uint32_t component_i);
extern bool sorter_is_component_skipped (uint32_t component_i);
extern bool sorter_is_vb_in_plan (uint32_t vb_i);
extern bool sorter_is_liftover_rejects (uint32_t component_i);
extern void sorter_filter_out_vb (uint32_t vb_i);
extern BufferP sorter_piz_get_region_X_ra_matrix (uint32_t vb_i);
extern void sorter_piz_get_txt_file_info (uint32_t txt_file_i, uint32_t *first_comp_i, uint32_t *num_comps, ConstSecLiEntP *start_sl);

#endif

