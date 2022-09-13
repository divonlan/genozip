// ------------------------------------------------------------------
//   txtheader.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"

//----------
// ZIP stuff
//----------

extern int64_t txtheader_zip_read_and_compress (int64_t *txt_header_offset, CompIType comp_i);

//----------
// PIZ stuff
//----------

extern void txtheader_piz_read_and_reconstruct (Section txt_header_sec);

