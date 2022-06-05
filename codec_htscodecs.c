// ------------------------------------------------------------------
//   codec_htscodecs.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// See: https://github.com/samtools/htscodecs

#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/arith_dynamic.h"
#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"

#define ORDER_B 0x01 // entropy-level-1
#define ORDER_W 0x19 // X_STRIPE | X_NOSZ | entropy-level-1
#define ORDER_b 0x81 // X_PACK | entropy-level-1
#define ORDER_w 0x99 // X_PACK| X_STRIPE | X_NOSZ | entropy-level-1

//----------------------------
// Get maximum compressed size
//----------------------------

uint32_t codec_ARTB_est_size (Codec codec, uint64_t uncompressed_len) { return arith_compress_bound    (uncompressed_len, ORDER_B); }
uint32_t codec_ARTW_est_size (Codec codec, uint64_t uncompressed_len) { return arith_compress_bound    (uncompressed_len, ORDER_W); }
uint32_t codec_ARTb_est_size (Codec codec, uint64_t uncompressed_len) { return arith_compress_bound    (uncompressed_len, ORDER_b); }
uint32_t codec_ARTw_est_size (Codec codec, uint64_t uncompressed_len) { return arith_compress_bound    (uncompressed_len, ORDER_w); }
uint32_t codec_RANB_est_size (Codec codec, uint64_t uncompressed_len) { return rans_compress_bound_4x16(uncompressed_len, ORDER_B); }
uint32_t codec_RANW_est_size (Codec codec, uint64_t uncompressed_len) { return rans_compress_bound_4x16(uncompressed_len, ORDER_W); }
uint32_t codec_RANb_est_size (Codec codec, uint64_t uncompressed_len) { return rans_compress_bound_4x16(uncompressed_len, ORDER_b); }
uint32_t codec_RANw_est_size (Codec codec, uint64_t uncompressed_len) { return rans_compress_bound_4x16(uncompressed_len, ORDER_w); }

//------------------------
// Compress
//------------------------

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
static bool codec_hts_compress (VBlockP vb, 
                                rom uncompressed,           // option 1 - compress contiguous data
                                uint32_t *uncompressed_len, 
                                LocalGetLineCB get_line_cb, // option 2 - compress data one line at a time
                                STRe(compressed),           // in/out 
                                unsigned char *(*func)(VBlockP vb, unsigned char *in,  unsigned int in_size, unsigned char *out, unsigned int *out_size, int order),                                
                                int order, bool soft_fail, rom name)
{
    START_TIMER;

    if (get_line_cb) {
        uncompressed = codec_alloc (vb, *uncompressed_len, 1);
        
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
            char *line; uint32_t line_len;
            get_line_cb (vb, line_i, pSTRa(line), *uncompressed_len - vb->codec_bufs[0].len, NULL);
            
            buf_add (&vb->codec_bufs[0], line, line_len);
        }
        ASSERT (vb->codec_bufs[0].len == *uncompressed_len, "\"%s\": Expecting in_so_far=%u == uncompressed_len=%u", name, vb->codec_bufs[0].len32, *uncompressed_len);
    }

    bool ret = !!func (vb, (uint8_t*)uncompressed, *uncompressed_len, (uint8_t*)compressed, compressed_len, order);

    if (func == rans_compress_to_4x16) COPY_TIMER_COMPRESS (compressor_rans);
    else                               COPY_TIMER_COMPRESS (compressor_arith); // higher level codecs are accounted for in their codec code

    return ret;
}

#define COMPRESS_FUNC_TEMPLATE(func,codec,order)                                                             \
COMPRESS (codec_##codec##_compress)                                                     \
{                                                                                                            \
    int ret;                                                                                                 \
    ASSERT ((ret = codec_hts_compress (vb, uncompressed, uncompressed_len, get_line_cb, compressed, compressed_len,     \
                                       func, order, soft_fail, name)) || soft_fail,                          \
            "Failed " #func " uncompressed_len=%u compressed_len=%u", *uncompressed_len, *compressed_len);   \
    return ret;                                                                                              \
}

COMPRESS_FUNC_TEMPLATE(rans_compress_to_4x16, RANB, ORDER_B)
COMPRESS_FUNC_TEMPLATE(rans_compress_to_4x16, RANb, ORDER_b)
COMPRESS_FUNC_TEMPLATE(rans_compress_to_4x16, RANW, ORDER_W)
COMPRESS_FUNC_TEMPLATE(rans_compress_to_4x16, RANw, ORDER_w)
COMPRESS_FUNC_TEMPLATE(arith_compress_to,     ARTB, ORDER_B)
COMPRESS_FUNC_TEMPLATE(arith_compress_to,     ARTb, ORDER_b)
COMPRESS_FUNC_TEMPLATE(arith_compress_to,     ARTW, ORDER_W)
COMPRESS_FUNC_TEMPLATE(arith_compress_to,     ARTw, ORDER_w)

//------------------------
// Uncompress
//------------------------

UNCOMPRESS (codec_rans_uncompress)
{
    START_TIMER;
    ASSERTNOTZERO (uncompressed_len, name);
    ASSERTNOTZERO (compressed_len, name);

    unsigned out_len = (unsigned)uncompressed_len;
    ASSERT (rans_uncompress_to_4x16 (vb, (uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, &out_len),
            "Failed rans_uncompress_to_4x16: \"%s\" compressed_len=%u uncompressed_len=%u (expected: %"PRIu64")", 
            name, compressed_len, out_len, uncompressed_len);

    ASSERT (out_len == uncompressed_len, "Expecting, for \"%s\" out_len=%u == uncompressed_len=%"PRIu64, name, out_len, uncompressed_len);

    COPY_TIMER (compressor_rans);
}

UNCOMPRESS (codec_arith_uncompress)
{
    START_TIMER;
    ASSERTNOTZERO (uncompressed_len, name);
    ASSERTNOTZERO (compressed_len, name);

    unsigned out_len = (unsigned)uncompressed_len;
    ASSERT (arith_uncompress_to (vb, (uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, &out_len),
            "Failed arith_uncompress_to: \"%s\" compressed_len=%u uncompressed_len=%"PRIu64, name, compressed_len, uncompressed_len);

    ASSERT (out_len == uncompressed_len, "For \"%s\": expecting out_len=%u == uncompressed_len=%"PRIu64, name, out_len, uncompressed_len);

    COPY_TIMER (compressor_arith);
}
