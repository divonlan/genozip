// ------------------------------------------------------------------
//   comp_none.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "comp_private.h"
#include "vblock.h"
#include "buffer.h"
#include "bsc/libbsc.h"

static const char *comp_bsc_errstr (int err)
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

bool comp_bsc_compress (VBlock *vb, Codec codec,
                        const char *uncompressed,      // option 1 - compress contiguous data
                        uint32_t uncompressed_len, 
                        LocalGetLineCallback callback, // option 2 - compress data one line at a time
                        char *compressed, uint32_t *compressed_len /* in/out */, 
                        bool soft_fail)
{
    ASSERT0 (*compressed_len >= uncompressed_len + LIBBSC_HEADER_SIZE, "Error in comp_bsc_compress: compressed_len too small");

    // libbsc doesn't allow piecemiel compression, so we need to copy all the data in case of callback
    if (callback) {

        // copy data to vb->compressed
        ASSERT0 (!vb->compressed.len, "Error in comp_bsc_compress: expecting vb->compress to be free, but its not");
        buf_alloc (vb, &vb->compressed, uncompressed_len, 1.2, "compressed", 0);

        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {
            char *start1=0, *start2=0;
            uint32_t len1=0, len2=0;        
            callback (vb, line_i, &start1, &len1, &start2, &len2);

            if (start1 && len1) buf_add (&vb->compressed, start1, len1);
            if (start2 && len2) buf_add (&vb->compressed, start2, len2);
        }

        uncompressed = vb->compressed.data;
    }

    int ret = bsc_compress (vb, (const uint8_t *)uncompressed, (uint8_t *)compressed, uncompressed_len,
                            16,  // lzp hash size
                            128, // lzp min size
                            LIBBSC_BLOCKSORTER_BWT, // block sorter 
                            flag_fast ? LIBBSC_CODER_QLFC_STATIC : LIBBSC_CODER_QLFC_ADAPTIVE, // coder
                            LIBBSC_FEATURE_FASTMODE); // flags ("features")

    if (ret == LIBBSC_NOT_COMPRESSIBLE)
        ret = bsc_store (vb, (const uint8_t *)uncompressed, (uint8_t *)compressed, uncompressed_len, LIBBSC_FEATURE_FASTMODE);

    ASSERT (ret >= LIBBSC_NO_ERROR, "Error in comp_bsc_compress: %s", comp_bsc_errstr (ret));

    if (callback) buf_free (&vb->compressed);

    *compressed_len = ret;

    return true;
}

void comp_bsc_uncompress (VBlock *vb, 
                          const char *compressed, uint32_t compressed_len,
                          char *uncompressed_data, uint64_t uncompressed_len)
{
    int ret = bsc_decompress (vb, (const uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_data, uncompressed_len, LIBBSC_FEATURE_FASTMODE);

    ASSERT (ret >= LIBBSC_NO_ERROR, "Error in comp_bsc_uncompress: %s", comp_bsc_errstr (ret));    
}

uint32_t comp_bsc_est_size (uint64_t uncompressed_len)
{
    return (uint32_t)uncompressed_len + LIBBSC_HEADER_SIZE; // as required by libbsc
}

static void *comp_bsc_malloc (void *vb, size_t size)
{
    void *mem =  comp_alloc ((VBlock *)vb, size, 1.5); 
    return mem;
}

static void comp_bsc_free (void *vb, void *addr)
{
    comp_free ((VBlock *)vb, addr);
}

void comp_bsc_initialize (void)
{
    bsc_init_full (LIBBSC_FEATURE_FASTMODE, comp_bsc_malloc, NULL, comp_bsc_free);
}

