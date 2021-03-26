// ------------------------------------------------------------------
//   compresssor.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef COMPRESSOR_INCLUDED
#define COMPRESSOR_INCLUDED

#include "genozip.h"
#include "codec.h"

extern uint32_t comp_compress (VBlockP vb, BufferP z_data, SectionHeaderP header, 
                               const char *uncompressed_data, // option 1 - compress contiguous data
                               LocalGetLineCB callback); // option 2 - compress data one line at a time

extern void comp_uncompress (VBlockP vb, Codec codec, Codec sub_codec, uint8_t param,
                             const char *compressed_data, uint32_t compressed_data_len,
                             Buffer *uncompressed_data, uint64_t uncompressed_len);

#endif
