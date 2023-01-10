// ------------------------------------------------------------------
//   codec_none.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"

COMPRESS (codec_none_compress)
{
    if (*compressed_len < *uncompressed_len && soft_fail) 
        return false;
    
    ASSERT (*compressed_len >= *uncompressed_len, "\"%s\": expecting compressed_len=%u >= uncompressed_len=%u for ", name, *compressed_len, *uncompressed_len);

    if (get_line_cb) {
        char *next = compressed;
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
            char *start1=0;
            uint32_t len1=0;        
            
            get_line_cb (vb, ctx, line_i, &start1, &len1, *uncompressed_len - (next - compressed), NULL);

            if (start1 && len1) { memcpy (next, start1, len1); next += len1; }
        }
    }
    else
        memcpy (compressed, uncompressed, *uncompressed_len);

    *compressed_len = *uncompressed_len;

    return true;
}

UNCOMPRESS (codec_none_uncompress)
{
    memcpy (uncompressed_buf->data, compressed, compressed_len);
}

uint32_t codec_none_est_size (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len;
}
