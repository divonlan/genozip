// ------------------------------------------------------------------
//   compresssor.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "codec.h"

extern uint32_t comp_compress (VBlockP vb, BufferP z_data, SectionHeaderP header, 
                               const char *uncompressed_data, // option 1 - compress contiguous data
                               LocalGetLineCB callback); // option 2 - compress data one line at a time

extern void comp_uncompress (VBlockP vb, Codec codec, Codec sub_codec, uint8_t param,
                             STRp(compressed_data), Buffer *uncompressed_data, uint64_t uncompressed_len);
