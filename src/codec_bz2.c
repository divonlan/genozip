// ------------------------------------------------------------------
//   comp_bz2.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "bzlib/bzlib.h"
#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"

static void *codec_bz2_alloc (void *vb_, int items, int size, FUNCLINE)
{
    return codec_alloc_do ((VBlockP )vb_, (uint64_t)items * (uint64_t)size, 1, func, code_line); // all bzlib buffers are constant in size between subsequent compressions
}

static rom BZ2_errstr (int err)
{
    switch (err) {
        case BZ_OK:               return "BZ_OK";
        case BZ_RUN_OK:           return "BZ_RUN_OK";
        case BZ_FLUSH_OK:         return "BZ_FLUSH_OK";
        case BZ_FINISH_OK:        return "BZ_FINISH_OK";
        case BZ_STREAM_END:       return "BZ_STREAM_END";
        case BZ_SEQUENCE_ERROR:   return "BZ_SEQUENCE_ERROR";
        case BZ_PARAM_ERROR:      return "BZ_PARAM_ERROR";
        case BZ_MEM_ERROR:        return "BZ_MEM_ERROR";
        case BZ_DATA_ERROR:       return "BZ_DATA_ERROR";
        case BZ_DATA_ERROR_MAGIC: return "BZ_DATA_ERROR_MAGIC";
        case BZ_IO_ERROR:         return "BZ_IO_ERROR";
        case BZ_UNEXPECTED_EOF:   return "BZ_UNEXPECTED_EOF";
        case BZ_OUTBUFF_FULL:     return "BZ_OUTBUFF_FULL";
        case BZ_CONFIG_ERROR:     return "BZ_CONFIG_ERROR";
        default:                  return "Unknown BZ2 error";
    }
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
COMPRESS (codec_bz2_compress)
{
    // good manual: http://linux.math.tifr.res.in/manuals/html/manual_3.html
    START_TIMER;
    
    bz_stream strm;
    memset (&strm, 0, sizeof (strm)); // safety

    strm.bzalloc = codec_bz2_alloc;
    strm.bzfree  = codec_free_do;
    strm.opaque  = vb; // just passed to malloc/free
    
    int init_ret = BZ2_bzCompressInit (&strm, flag.fast ? 1 : 9, 0, 30); // we optimize for size (normally) or speed (if user selected --fast)
    ASSERT (init_ret == BZ_OK, "%s: \"%s\": BZ2_bzCompressInit failed for ctx=%s: %s", VB_NAME, name, TAG_NAME, BZ2_errstr(init_ret));

    strm.next_out  = compressed;
    strm.avail_out = *compressed_len;
    bool success = true; // optimistic intialization
    int ret;

    // option 1 - compress contiguous data
    if (uncompressed) {
        strm.next_in  = (char*)uncompressed;
        strm.avail_in = *uncompressed_len;

        ret = BZ2_bzCompress (&strm, BZ_FINISH);
        if (soft_fail && ret == BZ_FINISH_OK)
            success = false; // data_compressed_len too small
        else 
            ASSERT (ret == BZ_STREAM_END, "%s: \"%s\": BZ2_bzCompress failed for ctx=%s: %s", VB_NAME, name, TAG_NAME, BZ2_errstr (ret));
    }
    
    // option 2 - compress data one line at a time
    else if (get_line_cb) {

        uint32_t in_so_far = 0;
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {

            ASSERT (!strm.avail_in, "%s: \"%s\": expecting strm.avail_in to be 0, but it is %u. ctx=%s", VB_NAME, name, strm.avail_in, TAG_NAME);

            // initialize to 0
            strm.next_in=0;
            strm.avail_in=0;

            get_line_cb (vb, ctx, line_i, &strm.next_in, &strm.avail_in, *uncompressed_len - in_so_far, NULL);
            in_so_far += strm.avail_in;

            bool is_last_line = (line_i == vb->lines.len - 1);
            if (!strm.avail_in && !is_last_line) continue; // this line has no SEQ data - move to next line (this happens eg in FASTA)

            bool final = is_last_line;

            ret = BZ2_bzCompress (&strm, final ? BZ_FINISH : BZ_RUN);

            if (soft_fail && ((!final && !strm.avail_out) || (final && ret != BZ_STREAM_END))) {
                success = false; // data_compressed_len too small
                break;
            }
            else 
                ASSERT (ret == (final ? BZ_STREAM_END : BZ_RUN_OK), "%s: \"%s\": BZ2_bzCompress failed for ctx=%s: %s", 
                        VB_NAME, name, TAG_NAME, BZ2_errstr (ret));
        }
    }
    else 
        ABORT ("%s: \"%s\": neither src_data nor callback is provided", VB_NAME, name);
    
    ret = BZ2_bzCompressEnd (&strm);
    ASSERT (ret == BZ_OK, "%s: \"%s\": BZ2_bzCompressEnd failed for ctx=%s: %s", VB_NAME, name, TAG_NAME, BZ2_errstr (ret));

    *compressed_len -= strm.avail_out;

    COPY_TIMER_COMPRESS (compressor_bz2); // higher level codecs are accounted for in their codec code

    return success;
}

UNCOMPRESS (codec_bz2_uncompress)
{
    START_TIMER;

    bz_stream strm;
    strm.bzalloc = codec_bz2_alloc;
    strm.bzfree  = codec_free_do;
    strm.opaque  = vb; // just passed to malloc/free

    int ret = BZ2_bzDecompressInit (&strm, 0, 0);
    ASSERT (ret == BZ_OK, "%s: \"%s\": BZ2_bzDecompressInit failed", VB_NAME, name);

    strm.next_in   = (char *)compressed;
    strm.avail_in  = compressed_len;
    strm.next_out  = uncompressed_buf->data;
    strm.avail_out = uncompressed_len;

    ret = BZ2_bzDecompress (&strm);
    ASSERT (ret == BZ_STREAM_END || ret == BZ_OK, "%s: \"%s\": BZ2_bzDecompress failed: %s, avail_in=%d, avail_out=%d", 
            VB_NAME, name, BZ2_errstr(ret), strm.avail_in, strm.avail_out);

    BZ2_bzDecompressEnd (&strm);

    COPY_TIMER (compressor_bz2);
}
