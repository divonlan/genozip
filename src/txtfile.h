// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"

extern rom txtfile_piz_get_filename (rom orig_name, rom prefix, bool is_orig_name_genozip);
extern rom txtfile_dump_vb (VBlockP vb, rom base_name);
extern void txtfile_read_header (bool is_first_txt);

#define TXTFILE_READ_VB_PADDING 16 // txtfile_read_vblock ensure this quantity of bytes at the end of vb.txt_data are unused
extern uint64_t txtfile_max_memory_per_vb (void);
extern void txtfile_read_vblock (VBlockP vb);
extern int64_t txtfile_get_seggable_size (void);

typedef bool (*TxtIteratorCallback)(rom line, unsigned line_len, void *cb_param1, void *cb_param2, unsigned cb_param3);
extern char *txtfile_foreach_line (BufferP txt_header, bool reverse, TxtIteratorCallback callback, void *cb_param1, void *cb_param2, unsigned cb_param3, int64_t *line_len);

// callbacks
extern int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern int32_t def_is_header_done (bool is_eof);

extern DataType txtfile_get_file_dt (rom filename);
