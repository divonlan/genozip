// ------------------------------------------------------------------
//   codec_htscodecs.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// See: https://github.com/samtools/htscodecs

#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/arith_dynamic.h"
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

uint32_t codec_ARTB_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + arith_compress_bound    (uncompressed_len, ORDER_B); } // +1 KB - see bug 1131
uint32_t codec_ARTW_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + arith_compress_bound    (uncompressed_len, ORDER_W); }
uint32_t codec_ARTb_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + arith_compress_bound    (uncompressed_len, ORDER_b); }
uint32_t codec_ARTw_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + arith_compress_bound    (uncompressed_len, ORDER_w); }
uint32_t codec_RANB_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + rans_compress_bound_4x16(uncompressed_len, ORDER_B); }
uint32_t codec_RANW_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + rans_compress_bound_4x16(uncompressed_len, ORDER_W); }
uint32_t codec_RANb_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + rans_compress_bound_4x16(uncompressed_len, ORDER_b); }
uint32_t codec_RANw_est_size (Codec codec, uint64_t uncompressed_len) { return 1 KB + rans_compress_bound_4x16(uncompressed_len, ORDER_w); }

//------------------------
// Compress
//------------------------

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
static bool codec_hts_compress (VBlockP vb, ContextP ctx,
                                rom uncompressed,           // option 1 - compress contiguous data
                                uint32_t *uncompressed_len, // function does NOT modify this value
                                LocalGetLineCB get_line_cb, // option 2 - compress data one line at a time
                                qSTRp(compressed),           // in/out 
                                uint8_t *(*func)(VBlockP vb, uint8_t *in, unsigned in_size, uint8_t *out, unsigned *out_size, int order),                                
                                int order, FailType soft_fail, rom name)
{
    START_TIMER;
    unsigned buf_i = 0;
    
    if (get_line_cb) {
        uncompressed = codec_alloc_do (vb, *uncompressed_len, 1, &buf_i, __FUNCLINE);
        
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {
            STRw (line);
            get_line_cb (vb, ctx, line_i, pSTRa(line), *uncompressed_len - vb->codec_bufs[buf_i].len, NULL);
            
            if (line_len)
                buf_add (&vb->codec_bufs[buf_i], STRa(line));
        }

        ASSERT (vb->codec_bufs[buf_i].len == *uncompressed_len, "%s: \"%s\": Expecting total_qual_len_from_callbacks=%u == uncompressed_len=%u ctx=%s", 
                VB_NAME, name, vb->codec_bufs[buf_i].len32, *uncompressed_len, TAG_NAME);
    }

    bool ret = !!func (vb, (uint8_t*)uncompressed, *uncompressed_len, (uint8_t*)compressed, compressed_len, order);

    if (func == rans_compress_to_4x16) COPY_TIMER_COMPRESS (compressor_rans);
    else                               COPY_TIMER_COMPRESS (compressor_arith); // higher level codecs are accounted for in their codec code

    if (get_line_cb) 
        buf_free (vb->codec_bufs[buf_i]); // don't rely on caller to free, because this function is called in many places outside of normal section compression

    return ret;
}

#define COMPRESS_FUNC_TEMPLATE(func,codec,order)                                                                                            \
COMPRESS (codec_##codec##_compress)                                                                                                         \
{                                                                                                                                           \
    int ret;                                                                                                                                \
    ASSERT ((ret = codec_hts_compress (vb, ctx, uncompressed, uncompressed_len, get_line_cb, compressed, compressed_len,                    \
                                       func, order, soft_fail, name)) || soft_fail,                                                         \
            "%s: Failed " #func " uncompressed_len=%u compressed_len=%u. ctx=%s", VB_NAME, *uncompressed_len, *compressed_len, TAG_NAME);   \
    return ret;                                                                                                                             \
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
    ASSERTNOTZEROn (uncompressed_len, name);
    ASSERTNOTZEROn (compressed_len, name);

    unsigned out_len = (unsigned)uncompressed_len;
    ASSERT (rans_uncompress_to_4x16 (vb, (uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, &out_len),
            "%s: Failed rans_uncompress_to_4x16: \"%s\" compressed_len=%u uncompressed_len=%u (expected: %"PRIu64")", 
            VB_NAME, name, compressed_len, out_len, uncompressed_len);

    ASSERT (out_len == uncompressed_len, "%s: Expecting, for \"%s\" out_len=%u == uncompressed_len=%"PRIu64, VB_NAME, name, out_len, uncompressed_len);

    COPY_TIMER (compressor_rans);
}

UNCOMPRESS (codec_arith_uncompress)
{
    START_TIMER;
    ASSERTNOTZEROn (uncompressed_len, name);
    ASSERTNOTZEROn (compressed_len, name);

    unsigned out_len = (unsigned)uncompressed_len;
    ASSERT (arith_uncompress_to (vb, (uint8_t *)compressed, compressed_len, (uint8_t *)uncompressed_buf->data, &out_len),
            "%s: Failed arith_uncompress_to: \"%s\" compressed_len=%u uncompressed_len=%"PRIu64, VB_NAME, name, compressed_len, uncompressed_len);

    ASSERT (out_len == uncompressed_len, "%s: For \"%s\": expecting out_len=%u == uncompressed_len=%"PRIu64, VB_NAME, name, out_len, uncompressed_len);

    COPY_TIMER (compressor_arith);
}
