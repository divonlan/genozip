// ------------------------------------------------------------------
//   txtheader.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"

//----------
// ZIP stuff
//----------

extern bool txtheader_zip_read_and_compress (int64_t *txt_header_offset, CompIType comp_i);

//----------
// PIZ stuff
//----------

extern void txtheader_piz_read_and_reconstruct (Section txt_header_sec);

