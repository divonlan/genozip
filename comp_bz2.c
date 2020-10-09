// ------------------------------------------------------------------
//   comp_bz2.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <bzlib.h>
#include "genozip.h"
#include "comp_private.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"

static void *comp_bzalloc (void *vb_, int items, int size)
{
    return comp_alloc ((VBlock *)vb_, items * size, 1); // all bzlib buffers are constant in size between subsequent compressions
}

static void comp_bzfree (void *vb_, void *addr)
{
    comp_free ((VBlock *)vb_, addr);
}

static const char *BZ2_errstr (int err)
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

// a hacky addition to bzlib - this should go into bzlib.c
uint64_t BZ2_consumed (void *bz_file)
{
    // note that this struct is not aligned to 32/64 bit. we need to trust that the bzlib compiler
    // and its options produce a similar alignment to ours...
    typedef struct {
        void *a;
        char  b[5000];
        int32_t c;
        uint8_t d;
        bz_stream strm;
    } bzFile;

    bz_stream *strm = &((bzFile*)bz_file)->strm;

    uint64_t total_in = ((uint64_t)(strm->total_in_hi32) << 32) |
                        ((uint64_t) strm->total_in_lo32);

    return total_in - strm->avail_in; // don't include unconsumed data
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool comp_bzlib_compress (VBlock *vb, Codec codec,
                          const char *uncompressed,       // option 1 - compress contiguous data
                          uint32_t uncompressed_len, 
                          LocalGetLineCallback callback,  // option 2 - compress data one line at a tim
                          char *compressed, uint32_t *compressed_len /* in/out */, 
                          bool soft_fail)
{
    // good manual: http://linux.math.tifr.res.in/manuals/html/manual_3.html
    START_TIMER;
    
    bz_stream strm;
    memset (&strm, 0, sizeof (strm)); // safety

    strm.bzalloc = comp_bzalloc;
    strm.bzfree  = comp_bzfree;
    strm.opaque  = vb; // just passed to malloc/free
    
    int init_ret = BZ2_bzCompressInit (&strm, flag_fast ? 1 : 9, 0, 30); // we optimize for size (normally) or speed (if user selected --fast)
    ASSERT (init_ret == BZ_OK, "Error: BZ2_bzCompressInit failed: %s", BZ2_errstr(init_ret));

    strm.next_out  = compressed;
    strm.avail_out = *compressed_len;
    bool success = true; // optimistic intialization
    int ret;

    // option 1 - compress contiguous data
    if (uncompressed) {
        strm.next_in   = (char*)uncompressed;
        strm.avail_in  = uncompressed_len;

        ret = BZ2_bzCompress (&strm, BZ_FINISH);
        if (soft_fail && ret == BZ_FINISH_OK)
            success = false; // data_compressed_len too small
        else 
            ASSERT (ret == BZ_STREAM_END, "Error: BZ2_bzCompress failed: %s", BZ2_errstr (ret));
    }
    
    // option 2 - compress data one line at a time
    else if (callback) {

        for (unsigned line_i=0; line_i < vb->lines.len; line_i++) {

            ASSERT (!strm.avail_in, "Error in comp_bzlib_compress: expecting strm.avail_in to be 0, but it is %u", strm.avail_in);

            // initialize to 0
            char *next_in_2=0;     strm.next_in=0;
            uint32_t avail_in_2=0; strm.avail_in=0;

            callback (vb, line_i, &strm.next_in, &strm.avail_in, &next_in_2, &avail_in_2);

            if (!strm.avail_in && !avail_in_2) continue; // this line has no SEQ data - move to next line (this happens eg in FASTA)

            bool final = (line_i == vb->lines.len - 1) && !avail_in_2;

            ret = BZ2_bzCompress (&strm, final ? BZ_FINISH : BZ_RUN);

            if (soft_fail && ((!final && !strm.avail_out) || (final && ret != BZ_STREAM_END))) {
                success = false; // data_compressed_len too small
                break;
            }
            else ASSERT (ret == (final ? BZ_STREAM_END : BZ_RUN_OK), 
                         "Error: BZ2_bzCompress failed: %s", BZ2_errstr (ret));

            // now the second part, if there is one
            if (avail_in_2) {
                final = (line_i == vb->lines.len - 1);

                strm.next_in  = next_in_2;
                strm.avail_in = avail_in_2;

                ret = BZ2_bzCompress (&strm, final ? BZ_FINISH : BZ_RUN);

                if (soft_fail && ret == BZ_FINISH_OK) { // TO DO - what is the condition for out of output space in BZ_RUN?
                    success = false; // data_compressed_len too small
                    break;
                }
                else ASSERT (ret == (final ? BZ_STREAM_END : BZ_RUN_OK), 
                             "Error: BZ2_bzCompress failed: %s", BZ2_errstr (ret));
            }
        }
    }
    else 
        ABORT0 ("Error in comp_bzlib_compress: neither src_data nor callback is provided");
    
    ret = BZ2_bzCompressEnd (&strm);
    ASSERT (ret == BZ_OK, "Error: BZ2_bzCompressEnd failed: %s", BZ2_errstr (ret));

    *compressed_len -= strm.avail_out;

    COPY_TIMER(vb->profile.compressor_bz2);

    return success;
}

void comp_bzlib_uncompress (VBlock *vb, 
                            const char *compressed, uint32_t compressed_len,
                            char *uncompressed_data, uint64_t uncompressed_len)
{
    bz_stream strm;
    strm.bzalloc = comp_bzalloc;
    strm.bzfree  = comp_bzfree;
    strm.opaque  = vb; // just passed to malloc/free

    int ret = BZ2_bzDecompressInit (&strm, 0, 0);
    ASSERT0 (ret == BZ_OK, "Error: BZ2_bzDecompressInit failed");

    strm.next_in   = (char *)compressed;
    strm.avail_in  = compressed_len;
    strm.next_out  = uncompressed_data;
    strm.avail_out = uncompressed_len;

    ret = BZ2_bzDecompress (&strm);
    ASSERT (ret == BZ_STREAM_END || ret == BZ_OK, "Error: BZ2_bzDecompress failed: %s, avail_in=%d, avail_out=%d", BZ2_errstr(ret), strm.avail_in, strm.avail_out);

    BZ2_bzDecompressEnd (&strm);
}
