// ------------------------------------------------------------------
//   codec_gtshark.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// support for compressing with --gtshark was available until v10 and for decompresssing until v11. it is now retired.

#include "genozip.h"

void codec_gtshark_uncompress (VBlockP vb_, Codec codec, uint8_t param,
                               const char *compressed, uint32_t compressed_len,
                               BufferP uncompressed_buf, uint64_t uncompressed_len,
                               Codec sub_codec)
{
    ABORT0 ("Support for the gtshark codec has been discontinued. Please use genozip v11 to decompress VCF files compressed with --gtshark");
}
