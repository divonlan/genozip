// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef TXTFILE_INCLUDED
#define TXTFILE_INCLUDED

#include "genozip.h"
#include "digest.h"

extern const char *txtfile_piz_get_filename (const char *orig_name, const char *prefix, bool is_orig_name_genozip);

typedef bool (*TxtFileTestFunc)(const char *, int);
extern bool txtfile_test_data (char first_char, unsigned num_lines_to_test, double success_threashold, TxtFileTestFunc test_func);

extern int64_t txtfile_estimate_txt_data_size (VBlockP vb);

extern const char *txtfile_dump_filename (VBlockP vb, const char *base_name, const char *ext);
extern const char *txtfile_dump_vb (VBlockP vb, const char *base_name);

extern Digest txtfile_read_header (bool is_first_txt);

#define TXTFILE_READ_VB_PADDING 16 // txtfile_read_vblock ensure this quanity of bytes at the end of vb.txt_data are unused
extern void txtfile_read_vblock (VBlockP vb, bool force_uncompress);

typedef bool (*TxtIteratorCallback)(const char *line, unsigned line_len, void *cb_param1, void *cb_param2, unsigned cb_param3);
extern char *txtfile_foreach_line (BufferP txt_header, bool reverse, TxtIteratorCallback callback, void *cb_param1, void *cb_param2, unsigned cb_param3, int64_t *line_len);

// callbacks
extern int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern int32_t def_is_header_done (bool is_eof);

extern DataType txtfile_get_file_dt (const char *filename);

#endif