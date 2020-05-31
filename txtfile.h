// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef TXTFILE_INCLUDED
#define TXTFILE_INCLUDED

#include "genozip.h"
#include "md5.h"

extern void txtfile_read_header (bool is_first_txt, bool header_required, char first_char);
extern void txtfile_read_vblock (VBlockP vb);
extern unsigned txtfile_write_to_disk (ConstBufferP buf);
extern void txtfile_estimate_txt_data_size (VBlockP vb);
extern void txtfile_write_one_vblock (VBlockP vb);

extern bool txtfile_header_to_genozip (uint32_t *vcf_line_i);
extern bool txtfile_genozip_to_txt_header (Md5Hash *digest);

extern void txtfile_header_initialize(void);
uint32_t txtfile_get_last_header_len(void); // for stats

#endif