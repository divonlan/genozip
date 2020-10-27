// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef TXTFILE_INCLUDED
#define TXTFILE_INCLUDED

#include "genozip.h"
#include "md5.h"

extern Md5Hash txtfile_read_header (bool is_first_txt, bool header_required, char first_char);

typedef bool (*TxtFileTestFunc)(const char *, int);
extern bool txtfile_test_data (char first_char, unsigned num_lines_to_test, double success_threashold, TxtFileTestFunc test_func);

extern void txtfile_estimate_txt_data_size (VBlockP vb);
extern void txtfile_write_one_vblock (VBlockP vb);

extern bool txtfile_header_to_genozip (uint32_t *vcf_line_i);
extern bool txtfile_genozip_to_txt_header (ConstSectionListEntryP sl, Md5Hash *digest);

extern void txtfile_header_initialize(void);
uint32_t txtfile_get_bound_headers_len(void); // for stats

extern uint32_t txtfile_read_block (char *data, uint32_t max_bytes);

extern uint32_t txtfile_unconsumed (VBlockP vb);
extern void txtfile_read_vblock (VBlockP vb);
extern void txtfile_update_md5 (const char *data, uint32_t len, bool is_2ndplus_txt_header);

#endif