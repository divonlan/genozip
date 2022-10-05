// ------------------------------------------------------------------
//   txtheader.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"

//----------
// ZIP stuff
//----------

extern void txtheader_compress (BufferP vcf_header_text, uint64_t unmodified_txt_header_len, Digest header_md5, bool is_first_vcf, CompIType comp_i);
extern int64_t txtheader_zip_read_and_compress (int64_t *txt_header_offset, CompIType comp_i);

//----------
// PIZ stuff
//----------

extern void txtheader_piz_read_and_reconstruct (Section txt_header_sec);

