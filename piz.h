// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"
#include "dispatcher.h"

#define piz_is_skip_section(vb,st,dict_id) (vb->data_type != DT_NONE     && (DT_FUNC (vb, is_skip_secetion) ((VBlockP)(vb), (st), (dict_id))))
#define piz_is_skip_sectionz(st,dict_id)   (z_file->data_type != DT_NONE && DTPZ(is_skip_secetion) && DTPZ(is_skip_secetion)(NULL, (st), (dict_id)))

extern Dispatcher piz_z_file_initialize (bool is_last_z_file);
extern DataType piz_read_global_area (void);
extern void piz_one_txt_file (Dispatcher dispatcher, int txt_file_i, bool is_last_txt_file, bool is_first_z_file);
extern uint32_t piz_uncompress_all_ctxs (VBlockP vb, uint32_t pair_vb_i);
extern bool piz_test_grep (VBlockP vb);

#endif


