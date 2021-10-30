// ------------------------------------------------------------------
//   codec_htscodecs.c
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// See: https://github.com/samtools/htscodecs

#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/arith_dynamic.h"
#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"

#define ORDER_8b  0x01 // entropy-level-1
#define ORDER_32b 0x19 // X_STRIPE | X_NOSZ | entropy-level-1

//----------------------------
// Get maximum compressed size
//----------------------------

// compressed size, as allocated in the original fqzcomp code
uint32_t codec_arith_est_size_8b (Codec codec, uint64_t uncompressed_len)
{
    return arith_compress_bound (uncompressed_len, ORDER_8b);
}

// compressed size, as allocated in the original fqzcomp code
uint32_t codec_arith_est_size_32b (Codec codec, uint64_t uncompressed_len)
{
    return arith_compress_bound (uncompressed_len, ORDER_32b);
}

// compressed size, as allocated in the original fqzcomp code
uint32_t codec_rans_est_size_8b (Codec codec, uint64_t uncompressed_len)
{
    return rans_compress_bound_4x16(uncompressed_len, ORDER_8b);
}
// compressed size, as allocated in the original fqzcomp code
uint32_t codec_rans_est_size_32b (Codec codec, uint64_t uncompressed_len)
{
    return rans_compress_bound_4x16(uncompressed_len, ORDER_32b);
}

//------------------------
// Compress
//------------------------

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
static bool codec_hts_compress (VBlock *vb, 
                                const char *uncompressed,       // option 1 - compress contiguous data
                                uint32_t *uncompressed_len, 
                                LocalGetLineCB callback,  // option 2 - compress data one line at a time
                                char *compressed, uint32_t *compressed_len /* in/out */, 
                                unsigned char *(*func)(unsigned char *in,  unsigned int in_size, unsigned char *out, unsigned int *out_size, int order),                                
                                int order)
{
    START_TIMER;

    if (callback) {
        // lengths, followed by FLAG, followed by in data
        buf_alloc (vb, &vb->codec_bufs[0], 0, *uncompressed_len, char, 1, "codec_bufs[0]");
        
        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {
            char *line; uint32_t line_len;
            callback (vb, line_i, pSTRa(line), *uncompressed_len - vb->codec_bufs[0].len);
            
            buf_add (&vb->codec_bufs[0], line, line_len);
        }
        ASSERT (vb->codec_bufs[0].len == *uncompressed_len, "Expecting in_so_far=%u == uncompressed_len=%u", (unsigned)vb->codec_bufs[0].len, *uncompressed_len);
        
        uncompressed = vb->codec_bufs[0].data;
    }

    func ((uint8_t*)uncompressed, *uncompressed_len, (uint8_t*)compressed, compressed_len, order);

    if (callback) buf_free (&vb->codec_bufs[0]);

    COPY_TIMER (compressor_rans); // higher level codecs are accounted for in their codec code

    return true;
}

bool codec_rans_compress_8b (VBlock *vb, SectionHeader *header,
                             const char *uncompressed,       // option 1 - compress contiguous data
                             uint32_t *uncompressed_len, 
                             LocalGetLineCB callback,  // option 2 - compress data one line at a time
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail)
{
    return codec_hts_compress (vb, uncompressed, uncompressed_len, callback, compressed, compressed_len, 
                               rans_compress_to_4x16, ORDER_8b);
}

bool codec_rans_compress_32b (VBlock *vb, SectionHeader *header,
                              const char *uncompressed,       // option 1 - compress contiguous data
                              uint32_t *uncompressed_len, 
                              LocalGetLineCB callback,  // option 2 - compress data one line at a time
                              char *compressed, uint32_t *compressed_len /* in/out */, 
                              bool soft_fail)
{
    return codec_hts_compress (vb, uncompressed, uncompressed_len, callback, compressed, compressed_len, 
                               rans_compress_to_4x16, ORDER_32b); // X_STRIPE | X_NOSZ | entropy-level-1
}

bool codec_arith_compress_8b (VBlock *vb, SectionHeader *header,
                             const char *uncompressed,       // option 1 - compress contiguous data
                             uint32_t *uncompressed_len, 
                             LocalGetLineCB callback,  // option 2 - compress data one line at a time
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail)
{
    return codec_hts_compress (vb, uncompressed, uncompressed_len, callback, compressed, compressed_len, 
                               arith_compress_to, ORDER_8b); // entropy-level-1
}

bool codec_arith_compress_32b (VBlock *vb, SectionHeader *header,
                               const char *uncompressed,       // option 1 - compress contiguous data
                               uint32_t *uncompressed_len, 
                               LocalGetLineCB callback,  // option 2 - compress data one line at a time
                               char *compressed, uint32_t *compressed_len /* in/out */, 
                               bool soft_fail)
{
    return codec_hts_compress (vb, uncompressed, uncompressed_len, callback, compressed, compressed_len, 
                               arith_compress_to, ORDER_32b); // X_STRIPE | X_NOSZ | entropy-level-1
}

//------------------------
// Uncompress
//------------------------

void codec_rans_uncompress (VBlock *vb, Codec codec, uint8_t param,
                           const char *compressed, uint32_t compressed_len,
                           Buffer *uncompressed_buf, uint64_t uncompressed_len, 
                           Codec unused)
{
    START_TIMER;

    unsigned out_len = (unsigned)uncompressed_len;
    rans_uncompress_to_4x16 ((uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, &out_len);

    ASSERT (out_len == uncompressed_len, "Expecting out_len=%u == uncompressed_len=%"PRIu64, out_len, uncompressed_len);

    COPY_TIMER (compressor_rans);
}

void codec_arith_uncompress (VBlock *vb, Codec codec, uint8_t param,
                             const char *compressed, uint32_t compressed_len,
                             Buffer *uncompressed_buf, uint64_t uncompressed_len, 
                             Codec unused)
{
    START_TIMER;

    unsigned out_len = (unsigned)uncompressed_len;
    arith_uncompress_to ((uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, &out_len);

    ASSERT (out_len == uncompressed_len, "Expecting out_len=%u == uncompressed_len=%"PRIu64, out_len, uncompressed_len);

    COPY_TIMER (compressor_arith);
}
