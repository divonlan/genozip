// ------------------------------------------------------------------
//   compressor.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <bzlib.h>
#include "lzma/LzmaEnc.h"
#include "lzma/LzmaDec.h"
#include "lzma/7zTypes.h"
#include "genozip.h"
#include "compressor.h"
#include "vblock.h"
#include "buffer.h"
#include "endianness.h"
#include "crypt.h"
#include "zfile.h"
#include "file.h"
#include "strings.h"

// -----------------------------------------------------
// memory functions that serve the compression libraries
// -----------------------------------------------------

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
static void *comp_alloc (VBlock *vb, int size, double grow_at_least_factor)
{
    // get the next buffer - allocations are always in the same order in bzlib and lzma -
    // so subsequent VBs will allocate roughly the same amount of memory for each buffer
    for (unsigned i=0; i < NUM_COMPRESS_BUFS ; i++) 
        if (!buf_is_allocated (&vb->compress_bufs[i])) {
            buf_alloc (vb, &vb->compress_bufs[i], size, grow_at_least_factor, "compress_bufs", i);
            return vb->compress_bufs[i].data;
        }

    ABORT ("Error: comp_alloc could not find a free buffer. vb_i=%d", vb->vblock_i);
    return 0; // squash compiler warning
}

static void comp_free (VBlock *vb, void *addr)
{
    if (!addr) return; // already freed

    for (unsigned i=0; i < NUM_COMPRESS_BUFS ; i++) 
        if (vb->compress_bufs[i].data == addr) {
            buf_free (&vb->compress_bufs[i]);
            return;
        }

    char addr_str[POINTER_STR_LEN];
    ABORT ("Error: comp_free failed to find buffer to free. vb_i=%d addr=%s", 
           vb->vblock_i, str_pointer (addr, addr_str));
}

static void comp_free_all (VBlock *vb)
{
    for (unsigned i=0; i < NUM_COMPRESS_BUFS ; i++) 
        buf_free (&vb->compress_bufs[i]);
}

static void *comp_bzalloc (void *vb_, int items, int size)
{
    return comp_alloc ((VBlock *)vb_, items * size, 1); // all bzlib buffers are constant in size between subsequent compressions
}

static void comp_bzfree (void *vb_, void *addr)
{
    comp_free ((VBlock *)vb_, addr);
}

static void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size)
{
//printf ("vb=%u size=%d\n", ((VBlock *)alloc_stuff->vb)->vblock_i, (int)size);
    return comp_alloc ((VBlock *)alloc_stuff->vb, size, 1.15); // lzma 5th buffer (the largest) may vary in size between subsequent compressions
}

static void lzma_free (ISzAllocPtr alloc_stuff, void *addr)
{
    comp_free ((VBlock *)alloc_stuff->vb, addr);
}

// -----------------------------------------------------
// bzlib stuff
// -----------------------------------------------------

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
bool comp_compress_bzlib (VBlock *vb, 
                          const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                          CompGetLineCallback callback,                        // option 2 - compress data one line at a tim
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

            ASSERT (!strm.avail_in, "Error in comp_compress_bzlib: expecting strm.avail_in to be 0, but it is %u", strm.avail_in);

            char *next_in_2;
            uint32_t avail_in_2;
            callback (vb, line_i, &strm.next_in, &strm.avail_in, &next_in_2, &avail_in_2);

            if (!strm.avail_in && !strm.avail_out) continue; // this line has no SEQ data - move to next line (this happens eg in FASTA)

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
        ABORT0 ("Error in comp_compress_bzlib: neither src_data nor callback is provided");
    
    ret = BZ2_bzCompressEnd (&strm);
    ASSERT (ret == BZ_OK, "Error: BZ2_bzCompressEnd failed: %s", BZ2_errstr (ret));

    *compressed_len -= strm.avail_out;

    COPY_TIMER(vb->profile.compressor);

    return success;
}

// -----------------------------------------------------
// lzma stuff
// -----------------------------------------------------

static const char *lzma_errstr (SRes res) 
{
    static const char *lzma_errors[] = { // from lzma/7zTypes.h
        "SZ_OK", "SZ_ERROR_DATA", "SZ_ERROR_MEM", "SZ_ERROR_CRC", "SZ_ERROR_UNSUPPORTED", "SZ_ERROR_PARAM", 
        "SZ_ERROR_INPUT_EOF", "SZ_ERROR_OUTPUT_EOF", "SZ_ERROR_READ", "SZ_ERROR_WRITE", "SZ_ERROR_PROGRESS",
        "SZ_ERROR_FAIL", "SZ_ERROR_THREAD", "Unknown lzma error", "Unknown lzma error", "Unknown lzma error",
        "SZ_ERROR_ARCHIVE", "SZ_ERROR_NO_ARCHIVE" };
    
    return lzma_errors[(unsigned)res <= 17 ? res : 13];
}

static const char *lzma_status (ELzmaStatus status)
{
    static const char *lzma_statuses[] = { // from lzma/LzmaDec.h
        "LZMA_STATUS_NOT_SPECIFIED", "LZMA_STATUS_FINISHED_WITH_MARK", "LZMA_STATUS_NOT_FINISHED",                /* stream was not finished */
        "LZMA_STATUS_NEEDS_MORE_INPUT", "LZMA_STATUS_MAYBE_FINISHED_WITHOUT_MARK" };

    return ((unsigned)status <= 4) ? lzma_statuses[status] : "Unrecognized lzma status";
}

static SRes comp_lzma_data_in_callback (const ISeqInStream *p, void *buf, size_t *size)
{
    ISeqInStream *instream = (ISeqInStream *)p; // discard the const

    // case: we're done serving all the data
    if (!instream->avail_in) {
        *size = 0;
        return SZ_OK;
    }

    // get next line if we have no data - keep on calling back until there is a line with data 
    // (not all lines must have seq data - for example, in FASTA they don't)
    while (instream->line_i < ((VBlockP)instream->vb)->lines.len && 
           !instream->avail_in_1 && !instream->avail_in_2) {

        instream->callback (instream->vb, instream->line_i, 
                            &instream->next_in_1, &instream->avail_in_1,
                            &instream->next_in_2, &instream->avail_in_2);

        instream->line_i++;
    }

    ASSERT (instream->avail_in_1 + instream->avail_in_2 <= instream->avail_in, 
            "Expecting avail_in_1=%u + avail_in_2=%u <= avail_in=%u but avail_in_1+avail_in_2=%u",
            instream->avail_in_1, instream->avail_in_2, instream->avail_in, instream->avail_in_1+instream->avail_in_2);
            
    uint32_t bytes_served_1 = MIN (instream->avail_in_1, *size);
    if (bytes_served_1) {
        memcpy (buf, instream->next_in_1, bytes_served_1);
        instream->next_in_1  += bytes_served_1;
        instream->avail_in_1 -= bytes_served_1;
    }

    uint32_t bytes_served_2 = MIN (instream->avail_in_2, *size - bytes_served_1);
    if (bytes_served_2) {    
        memcpy (buf + bytes_served_1, instream->next_in_2, bytes_served_2);
        instream->next_in_2  += bytes_served_2;
        instream->avail_in_2 -= bytes_served_2;
    }

    *size = bytes_served_1 + bytes_served_2;
    instream->avail_in -= bytes_served_1 + bytes_served_2;

    return SZ_OK;
}

static size_t comp_lzma_data_out_callback (const ISeqOutStream *p, const void *buf, size_t size)
{
    ISeqOutStream *outstream = (ISeqOutStream *)p; // discard the const

    uint32_t bytes_written = MIN (outstream->avail_out, (uint32_t)size);
    memcpy (outstream->next_out, buf, bytes_written);

    outstream->avail_out -= bytes_written;
    outstream->next_out  += bytes_written;

    return (size_t)bytes_written;
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool comp_compress_lzma (VBlock *vb, 
                         const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                         CompGetLineCallback callback,                        // option 2 - compress data one line at a tim
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
        uint64_t data_compressed_len64 = (uint64_t)*compressed_len - LZMA_PROPS_SIZE;
        res = LzmaEnc_MemEncode (lzma_handle, 
                                (uint8_t *)compressed + LZMA_PROPS_SIZE, &data_compressed_len64, 
                                (uint8_t *)uncompressed, uncompressed_len,
                                true, NULL, &alloc_stuff, &alloc_stuff);
        
        *compressed_len = (uint32_t)data_compressed_len64 + LZMA_PROPS_SIZE;
    }
    // option 2 - compress data one line at a time
    else if (callback) {

        ISeqInStream instream =   { .Read          = comp_lzma_data_in_callback, 
                                    .vb            = vb,
                                    .line_i        = 0,
                                    .avail_in      = uncompressed_len,
                                    .next_in_1     = NULL,
                                    .avail_in_1    = 0,
                                    .next_in_2     = NULL,
                                    .avail_in_2    = 0,
                                    .callback      = callback };
                                  
        ISeqOutStream outstream = { .Write        = comp_lzma_data_out_callback,
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

    COPY_TIMER(vb->profile.compressor);

    return success;
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool comp_compress_none (VBlock *vb, 
                         const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                         CompGetLineCallback callback,                        // option 2 - compress data one line at a tim
                         char *compressed, uint32_t *compressed_len /* in/out */, 
                         bool soft_fail)
{
    ASSERT0 (uncompressed, "Error in comp_compress_none: only option 1 supported in the mean time");

    if (*compressed_len < uncompressed_len && soft_fail) return false;
    ASSERT0 (*compressed_len >= uncompressed_len, "Error in comp_compress_none: compressed_len too small");
    
    memcpy (compressed, uncompressed, uncompressed_len);

    *compressed_len = uncompressed_len;

    return true;
}

bool comp_error (VBlock *vb, const char *uncompressed, uint32_t uncompressed_len, CompGetLineCallback callback,
                 char *compressed, uint32_t *compressed_len, bool soft_fail) 
{
    ABORT0 ("Error in comp_compress: Unsupported section compression algorithm");
    return false;
}

#define MIN_LEN_FOR_COMPRESSION 30 // less that this size its not worth compressing

// compresses data - either a conitguous block or one line at a time. If both are NULL that there is no data to compress.
void comp_compress (VBlock *vb, Buffer *z_data, bool is_z_file_buf,
                    SectionHeader *header, 
                    const char *uncompressed_data, // option 1 - compress contiguous data
                    CompGetLineCallback callback)  // option 2 - compress data one line at a time
{ 
    ASSERT0 (!uncompressed_data || !callback, "Error in comp_compress: expecting either uncompressed_data or callback but not both");

    // if the user requested --fast - we always use BZLIB, never LZMA
    if (flag_fast && header->sec_compression_alg == COMP_LZMA)
        header->sec_compression_alg = COMP_BZ2;

    static Compressor compressors[NUM_COMPRESSION_ALGS] = { 
        comp_compress_none, comp_error, comp_compress_bzlib, comp_error, comp_error, comp_error, comp_error, comp_compress_lzma, comp_error };

    ASSERT (header->sec_compression_alg < NUM_COMPRESSION_ALGS, "Error in comp_compress: unsupported section compressor=%u", header->sec_compression_alg);

    unsigned compressed_offset     = BGEN32 (header->compressed_offset);
    unsigned data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    unsigned data_compressed_len   = 0;
    unsigned data_encrypted_len=0, data_padding=0, header_padding=0;

    ASSERT0 (!data_uncompressed_len || uncompressed_data || callback, "Error in comp_compress: data_uncompressed_len!=0 but neither uncompressed_data nor callback are provided");

    bool is_encrypted = false;
    unsigned encryption_padding_reserve = 0;

    if (header->section_type != SEC_GENOZIP_HEADER) { // genozip header is never encrypted
        is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
        encryption_padding_reserve = crypt_max_padding_len(); // padding for the body
    }

    // if there's no data to compress, or its too small, don't compress
    if (data_uncompressed_len < MIN_LEN_FOR_COMPRESSION && !callback) 
        header->sec_compression_alg = COMP_PLN;

    uint32_t est_compressed_len = 
        (header->sec_compression_alg != COMP_PLN) ? MAX (data_uncompressed_len / 2, 500) : data_uncompressed_len;

    // allocate what we think will be enough memory. usually this alloc does nothing, as the memory we pre-allocate for z_data is sufficient
    // note: its ok for other threads to allocate evb data because we have a special mutex in buffer protecting the 
    // evb buffer list
    buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + est_compressed_len + encryption_padding_reserve, 1.5, z_data->name, z_data->param);

    // compress the data, if we have it...
    if (data_uncompressed_len) {
        
        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the alloc, because z_data is pre-allocated

        bool success = 
            compressors[header->sec_compression_alg](vb, uncompressed_data, data_uncompressed_len,
                                                     callback,  
                                                     &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                                     true);
        comp_free_all (vb); // just in case

        // if output buffer is too small, increase it, and try again
        if (!success) {
            buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + data_uncompressed_len  + encryption_padding_reserve + 50 /* > BZ_N_OVERSHOOT */, 1,
                       z_data->name ? z_data->name : "z_data", z_data->param);
            
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;

            compressors[header->sec_compression_alg](vb, 
                                                     uncompressed_data, data_uncompressed_len,
                                                     callback,  
                                                     &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                                     false);

            comp_free_all (vb); // just in case
        }
        
        // get encryption related lengths
        if (is_encrypted) {
            data_encrypted_len = data_compressed_len;
            crypt_get_encrypted_len (&data_encrypted_len, &data_padding); // both are set to 0 if no encryption
        }
    }

    // finalize & copy header
    header->compressed_offset   = BGEN32 (compressed_offset);  // updated value to include header padding
    header->data_compressed_len = BGEN32 (data_compressed_len);   
    header->data_encrypted_len  = BGEN32 (data_encrypted_len); 
    memcpy (&z_data->data[z_data->len], header, compressed_offset);

    // encrypt if needed - header & body separately
    unsigned total_z_len;
    if (is_encrypted) {
        // create good padding (just padding with 0 exposes a cryptographic volnurability)
        crypt_pad ((uint8_t*)&z_data->data[z_data->len], compressed_offset, header_padding);
        
        if (data_uncompressed_len) 
            crypt_pad ((uint8_t*)&z_data->data[z_data->len + compressed_offset], data_encrypted_len, data_padding);

        // encrypt the header - we use vb_i and section_i to generate a different AES key for each section
        uint32_t vb_i  = BGEN32 (header->vblock_i);

        // note: for SEC_VB_HEADER we will encrypt at the end of calculating this VB in zfile_update_compressed_vb_header() 
        // and we will then update z_data in memory prior to writing the encrypted data to disk
        if (header->section_type != SEC_VB_HEADER || header->vblock_i == 0 /* terminator vb header */)
            crypt_do (vb, (uint8_t*)&z_data->data[z_data->len], compressed_offset, vb_i, header->section_type, true);

        // encrypt the data body 
        if (data_uncompressed_len) 
            crypt_do (vb, (uint8_t*)&z_data->data[z_data->len + compressed_offset], data_encrypted_len, vb_i, header->section_type, false);
        
        total_z_len = compressed_offset + data_encrypted_len;
    }
    else
        total_z_len = compressed_offset + data_compressed_len;

    // add section to the list - except for genozip header which we already added in zfile_compress_genozip_header()
    if (header->section_type != SEC_GENOZIP_HEADER)
        sections_add_to_list (vb, header);

    z_data->len += total_z_len;

    if (flag_show_headers) zfile_show_header (header, vb->vblock_i ? vb : NULL); // store and print upon about for vb sections, and print immediately for non-vb sections
}

void comp_uncompress (VBlock *vb, CompressionAlg alg, 
                      const char *compressed, uint32_t compressed_len,
                      Buffer *uncompressed)
{
    switch (alg) {

    case COMP_BZ2: {
        bz_stream strm;
        strm.bzalloc = comp_bzalloc;
        strm.bzfree  = comp_bzfree;
        strm.opaque  = vb; // just passed to malloc/free

        int ret = BZ2_bzDecompressInit (&strm, 0, 0);
        ASSERT0 (ret == BZ_OK, "Error: BZ2_bzDecompressInit failed");

        strm.next_in   = (char *)compressed;
        strm.avail_in  = compressed_len;
        strm.next_out  = uncompressed->data;
        strm.avail_out = uncompressed->len;

        ret = BZ2_bzDecompress (&strm);
        ASSERT (ret == BZ_STREAM_END || ret == BZ_OK, "Error: BZ2_bzDecompress failed: %s, avail_in=%d, avail_out=%d", BZ2_errstr(ret), strm.avail_in, strm.avail_out);

        BZ2_bzDecompressEnd (&strm);
        break;
    }
    case COMP_LZMA: {
        ISzAlloc alloc_stuff = { .Alloc = lzma_alloc, .Free = lzma_free, .vb = vb};
        ELzmaStatus status;
        //uint64_t uncompressed_len64 = uncompressed->len;
        uint64_t compressed_len64 = (uint64_t)compressed_len - LZMA_PROPS_SIZE; // first 5 bytes in compressed stream are the encoding properties
        
        SRes ret = LzmaDecode ((uint8_t *)uncompressed->data, &uncompressed->len, 
                               (uint8_t *)compressed + LZMA_PROPS_SIZE, &compressed_len64, 
                               (uint8_t *)compressed, LZMA_PROPS_SIZE, 
                               LZMA_FINISH_END, &status, &alloc_stuff);

        ASSERT (ret == SZ_OK && status == LZMA_STATUS_FINISHED_WITH_MARK, 
                "Error: LzmaDecode failed: ret=%s status=%s", lzma_errstr (ret), lzma_status (status)); 

        break;
    }
    case COMP_PLN:
        memcpy (uncompressed->data, compressed, compressed_len);
        break;

    default:
        ABORT ("Error in comp_uncompress: invalid compression algorithm %u", alg);
    }

    comp_free_all (vb); // just in case
}
