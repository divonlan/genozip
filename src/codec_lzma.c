// ------------------------------------------------------------------
//   comp_bz2.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "lzma/LzmaEnc.h"
#include "lzma/LzmaDec.h"
#include "compressor.h"
#include "vblock.h"
#include "segconf.h"

static rom lzma_errstr (SRes res) 
{
    static rom lzma_errors[] = { // from lzma/7zTypes.h
        "SZ_OK", "SZ_ERROR_DATA", "SZ_ERROR_MEM", "SZ_ERROR_CRC", "SZ_ERROR_UNSUPPORTED", "SZ_ERROR_PARAM", 
        "SZ_ERROR_INPUT_EOF", "SZ_ERROR_OUTPUT_EOF", "SZ_ERROR_READ", "SZ_ERROR_WRITE", "SZ_ERROR_PROGRESS",
        "SZ_ERROR_FAIL", "SZ_ERROR_THREAD", "Unknown lzma error", "Unknown lzma error", "Unknown lzma error",
        "SZ_ERROR_ARCHIVE", "SZ_ERROR_NO_ARCHIVE" };
    
    return lzma_errors[(unsigned)res <= 17 ? res : 13];
}

static rom lzma_status (ELzmaStatus status)
{
    static rom lzma_statuses[] = { // from lzma/LzmaDec.h
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

        #pragma GCC diagnostic push 
        #pragma GCC diagnostic ignored "-Wpragmas"         // avoid warning if "-Wuse-after-free" is not defined in this version of gcc
        #pragma GCC diagnostic ignored "-Wunknown-warning-option" // same
        #pragma GCC diagnostic ignored "-Wdeprecated-non-prototype"
        instream->callback (instream->vb, instream->ctx, instream->line_i, 
                            &instream->next_in_1, &instream->avail_in_1,
                            instream->avail_in, NULL);
        #pragma GCC diagnostic pop

        instream->line_i++;
    }

    ASSERT (instream->avail_in_1 <= instream->avail_in, "Expecting avail_in_1=%u <= avail_in=%u. ctx=%s",
            instream->avail_in_1, instream->avail_in, (instream->ctx ? ((ContextP)instream->ctx)->tag_name : "NoContext"));
            
    uint32_t bytes_served_1 = MIN_(instream->avail_in_1, *size);
    if (bytes_served_1) {
        instream->next_in_1 = mempcpy (buf, instream->next_in_1, bytes_served_1);
        instream->avail_in_1 -= bytes_served_1;
    }

    *size = bytes_served_1;
    instream->avail_in -= bytes_served_1;

    return SZ_OK;
}

static size_t codec_lzma_data_out_callback (const ISeqOutStream *p, const void *buf, size_t size)
{
    ISeqOutStream *outstream = (ISeqOutStream *)p; // discard the const

    uint32_t bytes_written = MIN_(outstream->avail_out, (uint32_t)size);
    
    outstream->next_out = mempcpy (outstream->next_out, buf, bytes_written);
    outstream->avail_out -= bytes_written;

    return (size_t)bytes_written;
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
COMPRESS (codec_lzma_compress)
{
    START_TIMER;

    // for documentation on these parameters, see lzma/LzmaLib.h
    CLzmaEncProps props;
    LzmaEncProps_Init (&props);
    props.level        = 5;    // Without setting dictSize, Level 5 consumes < 200MB ; level 7 consumes up to 350MB per VB. negligible difference between level 5,7,9 (< 0.1% file size)
    props.fb           = 273;  // a bit better compression with no noticable impact on memory or speed
    props.writeEndMark = true; // add an "end of compression" mark - better error detection during decompress
    props.dictSize     = MIN_(*uncompressed_len, segconf.vb_size);

    char lzma_handle[LzmaEnc_LzmaHandleSize()];
    LzmaEnc_Create (lzma_handle, vb, ctx);

    SRes res = LzmaEnc_SetProps (lzma_handle, &props);
    ASSERT (res == SZ_OK, "%s: \"%s\": LzmaEnc_SetProps failed for ctx=%s: %s", VB_NAME, name, TAG_NAME, lzma_errstr (res));
    
    // write encoding properties as first 5 bytes of compressed stream
    SizeT props_size = LZMA_PROPS_SIZE; // per documentation in LzmaLib.h
    res = LzmaEnc_WriteProperties (lzma_handle, (uint8_t*)compressed, &props_size);
    ASSERT (res == SZ_OK && props_size==LZMA_PROPS_SIZE, "%s: \"%s\": LzmaEnc_WriteProperties failed for ctx=%s: %s", 
            VB_NAME, name, TAG_NAME, lzma_errstr (res));

    bool success = true;

    // option 1 - compress contiguous data
    if (uncompressed) {

        SizeT data_compressed_len64 = (SizeT)*compressed_len - LZMA_PROPS_SIZE;
        res = LzmaEnc_MemEncode (lzma_handle, 
                                (uint8_t *)compressed + LZMA_PROPS_SIZE, &data_compressed_len64, 
                                (uint8_t *)uncompressed, *uncompressed_len, true);
        
        *compressed_len = (uint32_t)data_compressed_len64 + LZMA_PROPS_SIZE;
    }
    // option 2 - compress data one line at a time
    else if (get_line_cb) {

        ISeqInStream instream =   { .Read          = codec_lzma_data_in_callback, 
                                    .vb            = vb,
                                    .line_i        = 0,
                                    .avail_in      = *uncompressed_len,
                                    .next_in_1     = NULL,
                                    .avail_in_1    = 0,
                                    .next_in_2     = NULL,
                                    .avail_in_2    = 0,
                                    .callback      = get_line_cb };
                                  
        ISeqOutStream outstream = { .Write        = codec_lzma_data_out_callback,
                                    .next_out     = compressed + LZMA_PROPS_SIZE,
                                    .avail_out    = *compressed_len - LZMA_PROPS_SIZE};
        
        res = LzmaEnc_Encode (lzma_handle, &outstream, &instream);        

        *compressed_len -= outstream.avail_out; 
    }

    if (soft_fail && (res == SZ_ERROR_WRITE || (!get_line_cb && res == SZ_ERROR_OUTPUT_EOF)))  // data_compressed_len is too small
        success = false;
    else
        ASSERT (res == SZ_OK, "%s: \"%s\": LzmaEnc_MemEncode failed for ctx=%s: %s", VB_NAME, name, TAG_NAME, lzma_errstr (res));

    LzmaEnc_Destroy (lzma_handle);

    COPY_TIMER_COMPRESS (compressor_lzma); // higher level codecs are accounted for in their codec code

    return success;
}

UNCOMPRESS (codec_lzma_uncompress)
{
    START_TIMER;
    
    ELzmaStatus status;

    SizeT compressed_len64 = (uint64_t)compressed_len - LZMA_PROPS_SIZE; // first 5 bytes in compressed stream are the encoding properties
    
    SRes ret = LzmaDecode (vb, (uint8_t *)uncompressed_buf->data, &uncompressed_len, 
                           (uint8_t *)compressed + LZMA_PROPS_SIZE, &compressed_len64, 
                           (uint8_t *)compressed, LZMA_PROPS_SIZE, 
                           LZMA_FINISH_END, &status);

    ASSERT (ret == SZ_OK && status == LZMA_STATUS_FINISHED_WITH_MARK, 
            "%s: \"%s\": LzmaDecode failed: ret=%s status=%s", VB_NAME, name, lzma_errstr (ret), lzma_status (status)); 

    COPY_TIMER (compressor_lzma);
}
