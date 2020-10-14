// ------------------------------------------------------------------
//   codec.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "codec.h"
#include "vblock.h"
#include "strings.h"

// --------------------------------------
// memory functions that serve the codecs
// --------------------------------------

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
void *codec_alloc (VBlock *vb, int size, double grow_at_least_factor)
{
    // get the next buffer - allocations are always in the same order in bzlib and lzma -
    // so subsequent VBs will allocate roughly the same amount of memory for each buffer
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (!buf_is_allocated (&vb->codec_bufs[i])) {
            buf_alloc (vb, &vb->codec_bufs[i], size, grow_at_least_factor, "codec_bufs", i);
            //printf ("codec_alloc: %u bytes buf=%u\n", size, i);
            return vb->codec_bufs[i].data;
        }

    ABORT ("Error: codec_alloc could not find a free buffer. vb_i=%d", vb->vblock_i);
    return 0; // squash compiler warning
}

void codec_free (VBlock *vb, void *addr)
{
    if (!addr) return; // already freed

    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        if (vb->codec_bufs[i].data == addr) {
            buf_free (&vb->codec_bufs[i]);
            //printf ("codec_free: buf=%u\n", i);
            return;
        }

    char addr_str[POINTER_STR_LEN];
    ABORT ("Error: codec_free failed to find buffer to free. vb_i=%d addr=%s", 
           vb->vblock_i, str_pointer (addr, addr_str));
}

void codec_free_all (VBlock *vb)
{
    for (unsigned i=0; i < NUM_CODEC_BUFS ; i++) 
        buf_free (&vb->codec_bufs[i]);
}

static bool codec_compress_error (VBlock *vb, Codec *codec, const char *uncompressed, uint32_t *uncompressed_len, LocalGetLineCB callback,
                                 char *compressed, uint32_t *compressed_len, bool soft_fail) 
{
    ABORT ("Error in comp_compress: Unsupported codec: %s", codec_name (*codec));
    return false;
}


static void codec_uncompress_error (VBlock *vb, Codec codec, 
                                   const char *compressed, uint32_t compressed_len,
                                   Buffer *uncompressed_buf, uint64_t uncompressed_len,
                                   Codec sub_codec)
{
    ABORT ("Error in comp_uncompress: Unsupported codec: %s", codec_name (codec));
}

static uint32_t codec_est_size_default (Codec codec, uint64_t uncompressed_len)
{
    return (uint32_t)MAX (uncompressed_len / 2, 500);
}

static uint32_t codec_sc1_est_size (Codec codec, uint64_t uncompressed_len)
{
    Codec sub_codec1 = codec_args[codec].sub_codec1;
    return codec_args[sub_codec1].est_size (sub_codec1, uncompressed_len);
}

const char *codec_name (Codec codec)
{
    return (codec >=0 && codec < NUM_CODECS) ? codec_args[codec].name : "BAD!";    
}

void codec_initialize (void)
{
    codec_bsc_initialize();
}

CodecArgs codec_args[NUM_CODECS] = CODEC_ARGS;
