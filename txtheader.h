// ------------------------------------------------------------------
//   txtheader.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"

//----------
// ZIP stuff
//----------

extern bool txtheader_zip_read_and_compress (uint64_t *txt_header_size);

//----------
// PIZ stuff
//----------

extern Coords txtheader_piz_read_and_reconstruct (uint32_t component_i, Section txt_header_sl);

