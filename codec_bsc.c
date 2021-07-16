// ------------------------------------------------------------------
//   codec_none.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "bsc/libbsc.h"

static const char *codec_bsc_errstr (int err)
{
    switch (err) {
        case LIBBSC_NO_ERROR:              return "LIBBSC_NO_ERROR";
        case LIBBSC_BAD_PARAMETER:         return "LIBBSC_BAD_PARAMETER";
        case LIBBSC_NOT_ENOUGH_MEMORY:     return "LIBBSC_NOT_ENOUGH_MEMORY";
        case LIBBSC_NOT_COMPRESSIBLE:      return "LIBBSC_NOT_COMPRESSIBLE";
        case LIBBSC_NOT_SUPPORTED:         return "LIBBSC_NOT_SUPPORTED";
        case LIBBSC_UNEXPECTED_EOB:        return "LIBBSC_UNEXPECTED_EOB";
        case LIBBSC_DATA_CORRUPT:          return "LIBBSC_DATA_CORRUPT";
        case LIBBSC_GPU_ERROR:             return "LIBBSC_GPU_ERROR";
        case LIBBSC_GPU_NOT_SUPPORTED:     return "LIBBSC_GPU_NOT_SUPPORTED";
        case LIBBSC_GPU_NOT_ENOUGH_MEMORY: return "LIBBSC_GPU_NOT_ENOUGH_MEMORY";
        default:                           return "Unknown BSC error";
    }
}

bool codec_bsc_compress (VBlock *vb, SectionHeader *header,
                         const char *uncompressed,      // option 1 - compress contiguous data
                         uint32_t *uncompressed_len, 
                         LocalGetLineCB callback, // option 2 - compress data one line at a time
                         char *compressed, uint32_t *compressed_len /* in/out */, 
                         bool soft_fail)
{
    START_TIMER;

    ASSERT (*compressed_len >= *uncompressed_len + LIBBSC_HEADER_SIZE, "compressed_len too small: compress_len=%u < uncompressed_len=%u + LIBBSC_HEADER_SIZE=%u",
            *compressed_len, *uncompressed_len, LIBBSC_HEADER_SIZE);

    // libbsc doesn't allow piecemiel compression, so we need to copy all the data in case of callback
    if (callback) {

        // copy data to vb->compressed
        ASSERTNOTINUSE (vb->compressed);
        buf_alloc (vb, &vb->compressed, 0, *uncompressed_len, char, 1.2, "compressed");

        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {
            char *start1=0;
            uint32_t len1=0;        
            
            // note: get what we need, might be less than what's available if calling from zip_assign_best_codec
            callback (vb, line_i, &start1, &len1, *uncompressed_len - vb->compressed.len); 

            if (start1 && len1) buf_add (&vb->compressed, start1, len1);
        }

        uncompressed = vb->compressed.data;
    }

    int ret = bsc_compress (vb, (const uint8_t *)uncompressed, (uint8_t *)compressed, *uncompressed_len,
                            16,  // lzp hash size
                            128, // lzp min size
                            LIBBSC_BLOCKSORTER_BWT, // block sorter 
                            flag.fast ? LIBBSC_CODER_QLFC_STATIC : LIBBSC_CODER_QLFC_ADAPTIVE, // coder
                            LIBBSC_FEATURE_FASTMODE); // flags ("features")

    if (ret == LIBBSC_NOT_COMPRESSIBLE)
        ret = bsc_store (vb, (const uint8_t *)uncompressed, (uint8_t *)compressed, *uncompressed_len, LIBBSC_FEATURE_FASTMODE);

    ASSERT (ret >= LIBBSC_NO_ERROR, "bsc_compress or bsc_store returned %s", codec_bsc_errstr (ret));

    if (callback) buf_free (&vb->compressed);

    *compressed_len = ret;

    COPY_TIMER (compressor_bsc); // higher level codecs are accounted for in their codec code

    return true;
}

void codec_bsc_uncompress (VBlock *vb, Codec codec, uint8_t param,
                          const char *compressed, uint32_t compressed_len,
                          Buffer *uncompressed_buf, uint64_t uncompressed_len, 
                          Codec unused)
{
    START_TIMER;
    int ret = bsc_decompress (vb, (const uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, uncompressed_len, LIBBSC_FEATURE_FASTMODE);

    ASSERT (ret >= LIBBSC_NO_ERROR, "%s", codec_bsc_errstr (ret));    

    COPY_TIMER (compressor_bsc);
}

uint32_t codec_bsc_est_size (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len + LIBBSC_HEADER_SIZE; // as required by libbsc
}

static void *codec_bsc_malloc (void *vb, size_t size)
{
    void *mem =  codec_alloc ((VBlock *)vb, size, 1.5); 
    return mem;
}

void codec_bsc_initialize (void)
{
    bsc_init_full (LIBBSC_FEATURE_FASTMODE, codec_bsc_malloc, codec_free);
}

