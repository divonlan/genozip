// ------------------------------------------------------------------
//   writer.h
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include "genozip.h"

// Reconstructin plan 
extern void writer_create_plan (void);
extern void writer_remove_vb_from_recon_plan (uint32_t vb_i);

// Writer thread
extern void writer_start_writing (uint32_t txt_file_i);
extern void writer_handover_data (VBlockP vb);
extern void writer_handover_txtheader (uint32_t component_i);
extern void writer_finish_writing (bool is_last_txt_file);

// VBlock and Txt Header properties
extern unsigned writer_get_pair (uint32_t vb_i, uint32_t *pair_vb_i);
extern bool writer_is_txtheader_in_plan (uint32_t component_i);
extern bool writer_is_component_no_read (uint32_t component_i);
extern bool writer_is_vb_no_read (uint32_t vb_i);
extern bool writer_is_liftover_rejects (uint32_t component_i);
extern void writer_get_txt_file_info (uint32_t txt_file_i, uint32_t *first_comp_i, uint32_t *num_comps, ConstSecLiEntP *start_sl);

#endif

