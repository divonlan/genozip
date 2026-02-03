// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "reconstruct.h"
#include "file.h"

extern bool piz_default_skip_section (SectionType st, DictId dict_id);

#define piz_is_skip_section(dt,st,comp_i,dict_id,f,preprocessing) \
    (vb->data_type != DT_NONE && (piz_default_skip_section ((st), (dict_id)) || \
    (dt_props[dt].is_skip_section && dt_props[dt].is_skip_section ((st), (comp_i), (dict_id), (f), (preprocessing)))))

#define piz_is_skip_undicted_section(st) (z_file->data_type != DT_NONE && DTPZ(is_skip_section) && DTPZ(is_skip_section)((st), COMP_NONE, DICT_ID_NONE, 0, false))

extern Dispatcher piz_z_file_initialize (void);
extern DataType piz_read_global_area (void);
extern void piz_one_txt_file (Dispatcher dispatcher, bool is_first_z_file, bool is_last_z_file, CompIType first_comp_i, CompIType last_comp_i, bool allow_skip_cleaning);
extern void piz_read_all_ctxs (VBlockP vb, Section *sec, bool is_pair_data);
typedef enum { PUR_RECON, PUR_FASTA_WRITER_INIT, PUR_FASTQ_READ_R1, PUR_SAM_LOAD_SAG } PizUncompressReason;
extern void piz_uncompress_all_ctxs (VBlockP vb, PizUncompressReason reason);
extern SectionHeaderVbHeader piz_read_vb_header (VBlockP vb);
extern bool piz_read_one_vb (VBlockP vb, bool for_reconstruction);
extern void piz_set_main_dispatcher (Dispatcher dispatcher);
extern void piz_allow_out_of_order (void);

extern bool piz_grep_match (rom start, rom after);

extern TRANSLATOR_FUNC (piz_obsolete_translator);

#define ASSISLOADED(ctx) ASSPIZ((ctx)->is_loaded, "%s is not loaded", ctx->tag_name)
