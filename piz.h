// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "reconstruct.h"
#include "sam.h"
#include "writer.h"

#define PIZ_TASK_NAME "piz"

extern bool piz_default_skip_section (SectionType st, DictId dict_id);
#define piz_is_skip_section(st,comp_i,dict_id,preprocessing) \
    (z_file->data_type != DT_NONE && (piz_default_skip_section ((st), (dict_id)) || \
    (DTPZ(is_skip_section) && DTPZ(is_skip_section)((st), (comp_i), (dict_id), (preprocessing)))))
#define piz_is_skip_undicted_section(st) (z_file->data_type != DT_NONE && DTPZ(is_skip_section) && DTPZ(is_skip_section)((st), COMP_NONE, DICT_ID_NONE, false))

extern Dispatcher piz_z_file_initialize (void);
extern DataType piz_read_global_area (Reference ref);
extern bool piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file, bool is_last_z_file, CompIType unbind_comp_i);
extern void piz_read_all_ctxs (VBlockP vb, Section *sec, bool is_pair_data);
extern uint32_t piz_uncompress_all_ctxs (VBlockP vb);
extern bool piz_read_one_vb (VBlockP vb, bool for_reconstruction);

extern bool piz_grep_match (rom start, rom after);

typedef struct { char s[100]; } PizDisCoords; 
extern PizDisCoords piz_dis_coords (VBlockP vb); // for ASSPIZ

typedef struct { char s[100]; } PizDisQname; 
extern PizDisQname piz_dis_qname (VBlockP vb); // for ASSPIZ

#define ASSPIZ(condition, format, ...) do { if (!(condition)) { progress_newline(); fprintf (stderr, "%s %s: Error in %s:%u line_in_file(1-based)=%"PRIu64" %s%s%s%s: ", str_time().s, LN_NAME, __FUNCLINE, writer_get_txt_line_i ((VBlockP)(vb)), (Z_DT(VCF) ? " sample_i=" : ""), (Z_DT(VCF) ? str_int_s (vb->sample_i).s : ""), piz_dis_coords((VBlockP)(vb)).s, piz_dis_qname((VBlockP)(vb)).s); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); exit_on_error(true); }} while(0)
#define ASSPIZ0(condition, string) ASSPIZ (condition, string "%s", "")
#define ASSISLOADED(ctx) ASSPIZ((ctx)->is_loaded, "%s is not loaded", ctx->tag_name)

extern bool piz_digest_failed;
