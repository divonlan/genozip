// ------------------------------------------------------------------
//   writer.h
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// Reconstruction plan 
extern void writer_create_plan (void);

// Writer thread
extern void writer_handover_data (VBlockP *vb_p);
extern void writer_handover_txtheader (VBlockP *comp_vb_p, uint32_t component_i);
extern void writer_finish_writing (bool is_last_txt_file);

// VBlock and Txt Header properties
extern unsigned writer_get_pair (uint32_t vb_i, uint32_t *pair_vb_i);
extern bool writer_is_txtheader_in_plan (uint32_t component_i);
extern bool writer_is_component_no_read (uint32_t component_i);
extern bool writer_is_vb_no_read (uint32_t vb_i);
extern bool writer_is_vb_full_vb (uint32_t vb_i);
extern void writer_get_txt_file_info (uint32_t *first_comp_i, uint32_t *num_comps, Section *start_sl);
