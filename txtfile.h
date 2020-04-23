// ------------------------------------------------------------------
//   txtfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef TXTFILE_INCLUDED
#define TXTFILE_INCLUDED

#include "genozip.h"

extern void txtfile_read_header (bool is_first_txt, bool header_required, char first_char);
extern void txtfile_read_vblock (VBlockP vb);
extern unsigned txtfile_write_to_disk (ConstBufferP buf);
extern void txtfile_estimate_txt_data_size (VBlockP vb);
extern void txtfile_write_one_vblock (VBlockP vb);

#endif