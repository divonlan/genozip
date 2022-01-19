// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "reconstruct.h"

extern bool piz_default_skip_section (VBlockP vb, SectionType st, DictId dict_id);
#define piz_is_skip_section(vb,st,dict_id) (vb->data_type     != DT_NONE && (piz_default_skip_section ((VBlockP)(vb), (st), (dict_id)) || (DT_FUNC (vb, is_skip_section)((VBlockP)(vb), (st), (dict_id)))))
#define piz_is_skip_sectionz(st,dict_id)   (z_file->data_type != DT_NONE && (piz_default_skip_section (NULL, (st), (dict_id)) || (DTPZ(is_skip_section) && DTPZ(is_skip_section)(NULL, (st), (dict_id)))))

extern Dispatcher piz_z_file_initialize (void);
extern DataType piz_read_global_area (Reference ref);
extern bool piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file, bool is_last_z_file);
extern uint32_t piz_uncompress_all_ctxs (VBlockP vb, uint32_t pair_vb_i);

extern bool piz_grep_match (const char *start, const char *after);
extern bool piz_test_grep (VBlockP vb);

typedef struct { char s[100]; } PizDisCoords; 
extern PizDisCoords piz_dis_coords (VBlockP vb); // for ASSPIZ

#define ASSPIZ(condition, format, ...) do { if (!(condition)) { progress_newline(); fprintf (stderr, "Error in %s:%u vb_i=%u line_in_file=%"PRIu64" line_i=%"PRId64"%s%s%s: ", __FUNCTION__, __LINE__, vb->vblock_i, vb->line_i, vb->line_i ? (vb->line_i - vb->first_line) : 0, (Z_DT(DT_VCF) ? " sample_i=" : ""), (Z_DT(DT_VCF) ? str_int_s (vb->sample_i).s : ""), piz_dis_coords(vb).s); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSPIZ0(condition, string)     do { if (!(condition)) { progress_newline(); fprintf (stderr, "Error in %s:%u vb_i=%u line_in_file=%"PRIu64" line_i=%"PRId64"%s: %s\n", __FUNCTION__, __LINE__, vb->vblock_i, vb->line_i, vb->line_i ? (vb->line_i - vb->first_line) : 0, piz_dis_coords(vb).s, string); exit_on_error(true); }} while(0)

extern bool piz_digest_failed;
