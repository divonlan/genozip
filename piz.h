// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"
#include "dispatcher.h"

extern bool piz_default_skip_section (VBlockP vb, SectionType st, DictId dict_id);
#define piz_is_skip_section(vb,st,dict_id) (vb->data_type     != DT_NONE && (piz_default_skip_section ((VBlockP)(vb), (st), (dict_id)) || (DT_FUNC (vb, is_skip_section)((VBlockP)(vb), (st), (dict_id)))))
#define piz_is_skip_sectionz(st,dict_id)   (z_file->data_type != DT_NONE && (piz_default_skip_section (NULL, (st), (dict_id)) || (DTPZ(is_skip_section) && DTPZ(is_skip_section)(NULL, (st), (dict_id)))))

extern Dispatcher piz_z_file_initialize (bool is_last_z_file);
extern DataType piz_read_global_area (Reference ref);
extern bool piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file);
extern uint32_t piz_uncompress_all_ctxs (VBlockP vb, uint32_t pair_vb_i);

extern bool piz_grep_match (const char *start, const char *after);
extern bool piz_test_grep (VBlockP vb);

#define ASSPIZ(condition, format, ...)       do { if (!(condition)) { progress_newline; fprintf (stderr, "Error in %s:%u vb_i=%u line_i=%"PRIu64": ", __FUNCTION__, __LINE__, vb->vblock_i, vb->line_i); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSPIZ0(condition, string)           do { if (!(condition)) { progress_newline; fprintf (stderr, "Error in %s:%u vb_i=%u line_i=%"PRIu64": %s\n", __FUNCTION__, __LINE__, vb->vblock_i, vb->line_i, string); exit_on_error(true); }} while(0)

#endif


