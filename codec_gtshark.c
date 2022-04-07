// ------------------------------------------------------------------
//   codec_gtshark.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// support for compressing with --gtshark was available until v10 and for decompresssing until v11. it is now retired.

#include "genozip.h"

void codec_gtshark_uncompress (VBlockP vb_, Codec codec, uint8_t param,
                               rom compressed, uint32_t compressed_len,
                               BufferP uncompressed_buf, uint64_t uncompressed_len,
                               Codec sub_codec, rom name)
{
    ABORT0 ("Support for the gtshark codec has been discontinued. Please use genozip v11 to decompress VCF files compressed with --gtshark");
}
