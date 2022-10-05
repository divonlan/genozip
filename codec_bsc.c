// ------------------------------------------------------------------
//   codec_none.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "bsc/libbsc.h"

static rom codec_bsc_errstr (int err)
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

COMPRESS (codec_bsc_compress)
{
    START_TIMER;

    ASSERT (*compressed_len >= *uncompressed_len + LIBBSC_HEADER_SIZE, "%s: compressed_len too small: \"%s\": compress_len=%u < uncompressed_len=%u + LIBBSC_HEADER_SIZE=%u. ctx=%s",
            VB_NAME, name, *compressed_len, *uncompressed_len, LIBBSC_HEADER_SIZE, TAG_NAME);

    // libbsc doesn't allow piecemiel compression, so we need to copy all the data in case of callback
    if (get_line_cb) {

        // copy data to vb->scratch
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc (vb, &vb->scratch, 0, *uncompressed_len, char, 1.2, "scratch");

        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
            char *start1=0;
            uint32_t len1=0;        
            
            // note: get what we need, might be less than what's available if calling from codec_assign_best_codec
            get_line_cb (vb, ctx, line_i, &start1, &len1, *uncompressed_len - vb->scratch.len, NULL); 

            if (start1 && len1) buf_add (&vb->scratch, start1, len1);
        }

        uncompressed = vb->scratch.data;
    }

    int ret = bsc_compress (vb, (bytes )uncompressed, (uint8_t *)compressed, *uncompressed_len,
                            16,  // lzp hash size
                            128, // lzp min size
                            LIBBSC_BLOCKSORTER_BWT, // block sorter 
                            flag.fast ? LIBBSC_CODER_QLFC_STATIC : LIBBSC_CODER_QLFC_ADAPTIVE, // coder
                            LIBBSC_FEATURE_FASTMODE); // flags ("features")

    if (ret == LIBBSC_NOT_COMPRESSIBLE) {
        ret = bsc_store (vb, (bytes )uncompressed, (uint8_t *)compressed, *uncompressed_len, LIBBSC_FEATURE_FASTMODE);

        ASSERT (ret >= LIBBSC_NO_ERROR, "%s: \"%s\": bsc_store returned %s. ctx=%s uncompressed_len=%u", 
                VB_NAME, name, codec_bsc_errstr (ret), TAG_NAME, *uncompressed_len);
    }

    ASSERT (ret >= LIBBSC_NO_ERROR, "%s: \"%s\": bsc_compress returned %s. ctx=%s uncompressed_len=%u", 
            VB_NAME, name, codec_bsc_errstr (ret), TAG_NAME, *uncompressed_len);

    if (get_line_cb) buf_free (vb->scratch);

    *compressed_len = ret;

    COPY_TIMER_COMPRESS (compressor_bsc); // higher level codecs are accounted for in their codec code

    return true;
}

UNCOMPRESS (codec_bsc_uncompress) 
{ 
    START_TIMER;
    int ret = bsc_decompress (vb, (bytes )compressed, compressed_len, (uint8_t *)uncompressed_buf->data, uncompressed_len, LIBBSC_FEATURE_FASTMODE);

    ASSERT (ret >= LIBBSC_NO_ERROR, "%s: \"%s\": %s", VB_NAME, name, codec_bsc_errstr (ret));    

    COPY_TIMER (compressor_bsc);
}

uint32_t codec_bsc_est_size (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len + LIBBSC_HEADER_SIZE; // as required by libbsc
}

static void *codec_bsc_malloc (void *vb, size_t size, FUNCLINE)
{
    void *mem =  codec_alloc_do (VB, size, 1.5, func, code_line); 
    return mem;
}

void codec_bsc_initialize (void)
{
    bsc_init_full (LIBBSC_FEATURE_FASTMODE, codec_bsc_malloc, codec_free_do);
}

