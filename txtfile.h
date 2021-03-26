// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef TXTFILE_INCLUDED
#define TXTFILE_INCLUDED

#include "genozip.h"
#include "digest.h"

extern const char *txtfile_piz_get_filename (const char *orig_name, const char *prefix, bool is_orig_name_genozip);

typedef bool (*TxtFileTestFunc)(const char *, int);
extern bool txtfile_test_data (char first_char, unsigned num_lines_to_test, double success_threashold, TxtFileTestFunc test_func);

extern int64_t txtfile_estimate_txt_data_size (VBlockP vb);
extern void txtfile_write_one_vblock (ConstBufferP txt_data, ConstBufferP bgzf_blocks, BufferP compressed, 
                                      uint32_t vb_i, uint32_t vb_data_size, uint32_t num_lines, uint32_t first_line);

extern void txtfile_write_4_lines (ConstBufferP txt_data, const char *line_start[5], unsigned pair);

extern const char *txtfile_dump_filename (VBlockP vb, const char *base_name, const char *ext);
extern const char *txtfile_dump_vb (VBlockP vb, const char *base_name);

extern bool txtfile_header_to_genozip (uint32_t *vcf_line_i);
extern void txtfile_genozip_to_txt_header (ConstSectionListEntryP sl, Digest *digest);

uint32_t txtfile_get_bound_headers_len(void); // for stats

extern void txtfile_read_vblock (VBlockP vb, bool force_uncompress);
extern void txtfile_write_to_disk (ConstBufferP buf, const char *data, unsigned len);
extern void txtfile_flush_if_stdout (void);

typedef bool (*TxtIteratorCallback)(const char *line, unsigned line_len, void *cb_param1, void *cb_param2, unsigned cb_param3);
extern char *txtfile_foreach_line (BufferP txt_header, bool reverse, TxtIteratorCallback callback, void *cb_param1, void *cb_param2, unsigned cb_param3, int64_t *line_len);

// callbacks
extern int32_t def_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern int32_t def_is_header_done (bool is_eof);

extern DataType txtfile_get_file_dt (const char *filename);

#endif