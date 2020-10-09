// ------------------------------------------------------------------
//   comp_none.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "comp_private.h"
#include "vblock.h"
#include "buffer.h"

bool comp_none_compress (VBlock *vb, Codec codec,
                         const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                         LocalGetLineCallback callback,                        // option 2 - compress data one line at a tim
                         char *compressed, uint32_t *compressed_len /* in/out */, 
                         bool soft_fail)
{
    if (*compressed_len < uncompressed_len && soft_fail) return false;
    ASSERT0 (*compressed_len >= uncompressed_len, "Error in comp_none_compress: compressed_len too small");

    if (callback) {
        char *next = compressed;
        for (unsigned line_i=0; line_i < vb->lines.len; line_i++) {
            char *start1=0, *start2=0;
            uint32_t len1=0, len2=0;        
            callback (vb, line_i, &start1, &len1, &start2, &len2);
            memcpy (next, start1, len1);
            next += len1;
            memcpy (next, start2, len2);
            next += len2;
        }
    }
    else
        memcpy (compressed, uncompressed, uncompressed_len);

    *compressed_len = uncompressed_len;

    return true;
}

void comp_none_uncompress (VBlock *vb, 
                           const char *compressed, uint32_t compressed_len,
                           char *uncompressed_data, uint64_t uncompressed_len)
{
    memcpy (uncompressed_data, compressed, compressed_len);
}

uint32_t comp_none_est_size (uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len;
}
