// ------------------------------------------------------------------
//   codec_none.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"

bool codec_none_compress (VBlock *vb, SectionHeader *header,
                         const char *uncompressed, uint32_t *uncompressed_len, // option 1 - compress contiguous data
                         LocalGetLineCB callback,                        // option 2 - compress data one line at a tim
                         char *compressed, uint32_t *compressed_len /* in/out */, 
                         bool soft_fail)
{
    if (*compressed_len < *uncompressed_len && soft_fail) return false;
    ASSERT0 (*compressed_len >= *uncompressed_len, "Error in codec_none_compress: compressed_len too small");

    if (callback) {
        char *next = compressed;
        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {
            char *start1=0;
            uint32_t len1=0;        
            
            callback (vb, line_i, &start1, &len1, *uncompressed_len - (next - compressed));

            if (start1 && len1) { memcpy (next, start1, len1); next += len1; }
        }
    }
    else
        memcpy (compressed, uncompressed, *uncompressed_len);

    *compressed_len = *uncompressed_len;

    return true;
}

void codec_none_uncompress (VBlock *vb, Codec codec,
                           const char *compressed, uint32_t compressed_len,
                           Buffer *uncompressed_buf, uint64_t uncompressed_len, 
                           Codec unused)
{
    memcpy (uncompressed_buf->data, compressed, compressed_len);
}

uint32_t codec_none_est_size (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len;
}
