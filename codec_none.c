// ------------------------------------------------------------------
//   codec_none.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"

bool codec_none_compress (VBlockP vb, SectionHeader *header,
                          rom uncompressed,   // option 1 - compress contiguous data
                          uint32_t *uncompressed_len, // in/out (might be modified by complex codecs)
                          LocalGetLineCB callback,    // option 2 - compress data one line at a tim
                          char *compressed, uint32_t *compressed_len /* in/out */, 
                          bool soft_fail, rom name)
{
    if (*compressed_len < *uncompressed_len && soft_fail) 
        return false;
    
    ASSERT (*compressed_len >= *uncompressed_len, "\"%s\": expecting compressed_len=%u >= uncompressed_len=%u for ", name, *compressed_len, *uncompressed_len);

    if (callback) {
        char *next = compressed;
        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {
            char *start1=0;
            uint32_t len1=0;        
            
            callback (vb, line_i, &start1, &len1, *uncompressed_len - (next - compressed), NULL);

            if (start1 && len1) { memcpy (next, start1, len1); next += len1; }
        }
    }
    else
        memcpy (compressed, uncompressed, *uncompressed_len);

    *compressed_len = *uncompressed_len;

    return true;
}

void codec_none_uncompress (VBlockP vb, Codec codec, uint8_t param,
                           rom compressed, uint32_t compressed_len,
                           Buffer *uncompressed_buf, uint64_t uncompressed_len, 
                           Codec unused, rom name)
{
    memcpy (uncompressed_buf->data, compressed, compressed_len);
}

uint32_t codec_none_est_size (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len;
}
