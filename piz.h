// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"

#define piz_is_skip_section(vb,st,dict_id) (vb->data_type != DT_NONE     && (DT_FUNC (vb, is_skip_secetion) ((VBlockP)(vb), (st), (dict_id))))
#define piz_is_skip_sectionz(st,dict_id)   (z_file->data_type != DT_NONE && DTPZ(is_skip_secetion) && DTPZ(is_skip_secetion)(NULL, (st), (dict_id)))

extern void piz_one_file (uint32_t component_i, bool is_first_z_file, bool is_last_file);
extern uint32_t piz_uncompress_all_ctxs (VBlockP vb, uint32_t pair_vb_i);
extern bool piz_test_grep (VBlockP vb);

#endif

