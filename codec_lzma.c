// ------------------------------------------------------------------
//   comp_bz2.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "lzma/LzmaEnc.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"

void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size)
{
    return codec_alloc ((VBlock *)alloc_stuff->vb, size, 1.15); // lzma 5th buffer (the largest) may vary in size between subsequent compressions
}

void lzma_free (ISzAllocPtr alloc_stuff, void *addr)
{
    codec_free ((VBlock *)alloc_stuff->vb, addr);
}

const char *lzma_errstr (SRes res) 
{
    static const char *lzma_errors[] = { // from lzma/7zTypes.h
        "SZ_OK", "SZ_ERROR_DATA", "SZ_ERROR_MEM", "SZ_ERROR_CRC", "SZ_ERROR_UNSUPPORTED", "SZ_ERROR_PARAM", 
        "SZ_ERROR_INPUT_EOF", "SZ_ERROR_OUTPUT_EOF", "SZ_ERROR_READ", "SZ_ERROR_WRITE", "SZ_ERROR_PROGRESS",
        "SZ_ERROR_FAIL", "SZ_ERROR_THREAD", "Unknown lzma error", "Unknown lzma error", "Unknown lzma error",
        "SZ_ERROR_ARCHIVE", "SZ_ERROR_NO_ARCHIVE" };
    
    return lzma_errors[(unsigned)res <= 17 ? res : 13];
}

const char *lzma_status (ELzmaStatus status)
{
    static const char *lzma_statuses[] = { // from lzma/LzmaDec.h
        "LZMA_STATUS_NOT_SPECIFIED", "LZMA_STATUS_FINISHED_WITH_MARK", "LZMA_STATUS_NOT_FINISHED",                /* stream was not finished */
        "LZMA_STATUS_NEEDS_MORE_INPUT", "LZMA_STATUS_MAYBE_FINISHED_WITHOUT_MARK" };

    return ((unsigned)status <= 4) ? lzma_statuses[status] : "Unrecognized lzma status";
}


static SRes codec_lzma_data_in_callback (const ISeqInStream *p, void *buf, size_t *size)
{
    ISeqInStream *instream = (ISeqInStream *)p; // discard the const

    // case: we're done serving all the data
    if (!instream->avail_in) {
        *size = 0; // we're done
        return SZ_OK;
    }

    // get next line if we have no data - keep on calling back until there is a line with data 
    // (not all lines must have seq data - for example, in FASTA they don't or SAM optional fields BI/BD/E2/U2 might not appear on every line)
    while (instream->line_i < ((VBlockP)instream->vb)->lines.len && 
           !instream->avail_in_1 && !instream->avail_in_2) {

        // initialize to 0
        instream->next_in_1  = instream->next_in_2  = 0;
        instream->avail_in_1 = instream->avail_in_2 = 0;

        instream->callback (instream->vb, instream->line_i, 
                            &instream->next_in_1, &instream->avail_in_1,
                            instream->avail_in);

        instream->line_i++;
    }

    ASSERT (instream->avail_in_1 <= instream->avail_in, "Expecting avail_in_1=%u <= avail_in=%u",
            instream->avail_in_1, instream->avail_in);
            
    uint32_t bytes_served_1 = MIN (instream->avail_in_1, *size);
    if (bytes_served_1) {
        memcpy (buf, instream->next_in_1, bytes_served_1);
        instream->next_in_1  += bytes_served_1;
        instream->avail_in_1 -= bytes_served_1;
    }

    *size = bytes_served_1;
    instream->avail_in -= bytes_served_1;

    return SZ_OK;
}

static size_t codec_lzma_data_out_callback (const ISeqOutStream *p, const void *buf, size_t size)
{
    ISeqOutStream *outstream = (ISeqOutStream *)p; // discard the const

    uint32_t bytes_written = MIN (outstream->avail_out, (uint32_t)size);
    memcpy (outstream->next_out, buf, bytes_written);

    outstream->avail_out -= bytes_written;
    outstream->next_out  += bytes_written;

    return (size_t)bytes_written;
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool codec_lzma_compress (VBlock *vb, SectionHeader *header,
                         const char *uncompressed,    // option 1 - compress contiguous data
                         uint32_t *uncompressed_len,
                         LocalGetLineCB callback,     // option 2 - compress data one line at a time
                         char *compressed, uint32_t *compressed_len /* in/out */, 
                         bool soft_fail)
{
    START_TIMER;
    ISzAlloc alloc_stuff = { .Alloc = lzma_alloc, .Free = lzma_free, .vb = vb};

    // for documentation on these parameters, see lzma/LzmaLib.h
    CLzmaEncProps props;
    LzmaEncProps_Init (&props);
    props.level        = 5;    // Level 5 consumes < 200MB ; level 7 consumes up to 350MB per VB. negligible difference between level 5,7,9 (< 0.1% file size)
    props.fb           = 273;  // a bit better compression with no noticable impact on memory or speed
    props.writeEndMark = true; // add an "end of compression" mark - better error detection during decompress

    CLzmaEncHandle lzma_handle = LzmaEnc_Create (&alloc_stuff);
    ASSERT0 (lzma_handle, "Error: LzmaEnc_Create failed");

    SRes res = LzmaEnc_SetProps (lzma_handle, &props);
    ASSERT (res == SZ_OK, "Error: LzmaEnc_SetProps failed: %s", lzma_errstr (res));
    
    // write encoding properties as first 5 bytes of compressed stream
    SizeT props_size = LZMA_PROPS_SIZE; // per documentation in LzmaLib.h
    res = LzmaEnc_WriteProperties (lzma_handle, (uint8_t*)compressed, &props_size);
    ASSERT (res == SZ_OK && props_size==LZMA_PROPS_SIZE, "Error: LzmaEnc_WriteProperties failed: %s", lzma_errstr (res));

    bool success = true;

    // option 1 - compress contiguous data
    if (uncompressed) {

        SizeT data_compressed_len64 = (SizeT)*compressed_len - LZMA_PROPS_SIZE;
        res = LzmaEnc_MemEncode (lzma_handle, 
                                (uint8_t *)compressed + LZMA_PROPS_SIZE, &data_compressed_len64, 
                                (uint8_t *)uncompressed, *uncompressed_len,
                                true, NULL, &alloc_stuff, &alloc_stuff);
        
        *compressed_len = (uint32_t)data_compressed_len64 + LZMA_PROPS_SIZE;
    }
    // option 2 - compress data one line at a time
    else if (callback) {

        ISeqInStream instream =   { .Read          = codec_lzma_data_in_callback, 
                                    .vb            = vb,
                                    .line_i        = 0,
                                    .avail_in      = *uncompressed_len,
                                    .next_in_1     = NULL,
                                    .avail_in_1    = 0,
                                    .next_in_2     = NULL,
                                    .avail_in_2    = 0,
                                    .callback      = callback };
                                  
        ISeqOutStream outstream = { .Write        = codec_lzma_data_out_callback,
                                    .next_out     = compressed + LZMA_PROPS_SIZE,
                                    .avail_out    = *compressed_len - LZMA_PROPS_SIZE};
        
        res = LzmaEnc_Encode (lzma_handle, &outstream, &instream, NULL, &alloc_stuff, &alloc_stuff);        

        *compressed_len -= outstream.avail_out; 
    }

    if (soft_fail && ((callback && res == SZ_ERROR_WRITE) || (!callback && res == SZ_ERROR_OUTPUT_EOF)))  // data_compressed_len is too small
        success = false;
    else
        ASSERT (res == SZ_OK, "Error: LzmaEnc_MemEncode failed: %s", lzma_errstr (res));

    LzmaEnc_Destroy (lzma_handle, &alloc_stuff, &alloc_stuff);

    COPY_TIMER (compressor_lzma); // higher level codecs are accounted for in their codec code

    return success;
}

void codec_lzma_uncompress (VBlock *vb, Codec codec,
                           const char *compressed, uint32_t compressed_len,
                           Buffer *uncompressed_buf, uint64_t uncompressed_len, 
                           Codec unused)
{
    ISzAlloc alloc_stuff = { .Alloc = lzma_alloc, .Free = lzma_free, .vb = vb};
    ELzmaStatus status;

    SizeT compressed_len64 = (uint64_t)compressed_len - LZMA_PROPS_SIZE; // first 5 bytes in compressed stream are the encoding properties
    
    SRes ret = LzmaDecode ((uint8_t *)uncompressed_buf->data, &uncompressed_len, 
                           (uint8_t *)compressed + LZMA_PROPS_SIZE, &compressed_len64, 
                           (uint8_t *)compressed, LZMA_PROPS_SIZE, 
                           LZMA_FINISH_END, &status, &alloc_stuff);

    ASSERT (ret == SZ_OK && status == LZMA_STATUS_FINISHED_WITH_MARK, 
            "Error in codec_lzma_uncompress: LzmaDecode failed: ret=%s status=%s", lzma_errstr (ret), lzma_status (status)); 
}
