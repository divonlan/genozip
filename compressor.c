// ------------------------------------------------------------------
//   compressor.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <bzlib.h>
#include "lzma/LzmaEnc.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "compressor.h"
#include "vblock.h"
#include "buffer.h"
#include "endianness.h"
#include "crypt.h"
#include "zfile.h"

#define BZLIB_BLOCKSIZE100K 9 /* maximum mem allocation for bzlib */

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
uint64_t BZ2_consumed (BZFILE* b)
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

    bz_stream *strm = &((bzFile*)b)->strm;

    uint64_t total_in = ((uint64_t)(strm->total_in_hi32) << 32) |
                        ((uint64_t) strm->total_in_lo32);

    return total_in - strm->avail_in; // don't include unconsumed data
}

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
static void *comp_bzalloc (void *vb_, int items, int size)
{
    VBlock *vb = (VBlock *)vb_;
    Buffer *use_buf = NULL;
    
    // search for a free buffer - preferring one of the size requested 
   for (unsigned i=0; i < NUM_COMPRESS_BUFS ; i++) {
        if (!buf_is_allocated (&vb->compress_bufs[i])) {
            use_buf = &vb->compress_bufs[i];
            if (use_buf->size == (unsigned)(items * size)) break;
        }
    }
    ASSERT0 (use_buf, "Error: comp_bzalloc could not find a free buffer");

    buf_alloc (vb, use_buf, items * size, 1, "compress_bufs", 0);

    return use_buf->data;
}

static void comp_bzfree (void *vb_, void *addr)
{
    VBlock *vb = (VBlock *)vb_;
    //Buffer *bufs = vb ? vb->compress_bufs : compress_bufs;
    Buffer *bufs = vb->compress_bufs;

    unsigned i; for (i=0; i < NUM_COMPRESS_BUFS ; i++) 
        if (bufs[i].data == addr) {
            buf_free (&bufs[i]);
            break;
        }

    ASSERT0 (i < NUM_COMPRESS_BUFS, "Error: comp_bzfree failed to find buffer to free");
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

    int init_ret = BZ2_bzCompressInit (&strm, BZLIB_BLOCKSIZE100K, 0, 30);
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

        for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {

            ASSERT (!strm.avail_in, "Error in comp_compress_bzlib: expecting strm.avail_in to be 0, but it is %u", strm.avail_in);

            callback (vb, line_i, &strm.next_in, &strm.avail_in);

            if (line_i < vb->num_lines - 1) { // not last line
//printf("BEFORE: line_i=%u compress_ret=%d avail_in=%u avail_out=%u\n", line_i, compress_ret, strm.avail_in, strm.avail_out);                
                ret = BZ2_bzCompress (&strm, BZ_RUN);
//printf("AFTER:  line_i=%u compress_ret=%d avail_in=%u avail_out=%u\n", line_i, compress_ret, strm.avail_in, strm.avail_out); // DEBUG
                if (soft_fail && !strm.avail_out /* ret == BZ_FINISH_OK */) {
                    success = false; // data_compressed_len too small
                    break;
                }
                else 
                    ASSERT (ret == BZ_RUN_OK, "Error: BZ2_bzCompress failed: %s", BZ2_errstr (ret));
            }
            else { // last line
                ret = BZ2_bzCompress (&strm, BZ_FINISH);
                if (soft_fail && ret == BZ_FINISH_OK)
                    success = false; // data_compressed_len too small
                else 
                    ASSERT (ret == BZ_STREAM_END, "Error: BZ2_bzCompress failed: %s", BZ2_errstr (ret));
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

void *lzma_alloc (ISzAllocPtr p, size_t size)
{
  fprintf(stderr, "\nAlloc %10u bytes", (unsigned)size);
  return malloc (size);
}

void lzma_free (ISzAllocPtr p, void *address)
{
  // if (address != 0) fprintf(stderr, "\nFree; count = %10d", g_allocCount);
  free (address);
}

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

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
bool comp_compress_lzma (VBlock *vb, 
                         const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                         CompGetLineCallback callback,                        // option 2 - compress data one line at a tim
                         char *compressed, uint32_t *compressed_len /* in/out */, 
                         bool soft_fail)
{
    START_TIMER;
    static const ISzAlloc memfuncs = { lzma_alloc, lzma_free };

    // for documentation on these parameters, see lzma/LzmaLib.h
    CLzmaEncProps props;
    LzmaEncProps_Init (&props);
    props.level = 7; // level 7 consumes up to 350MB. level 9 can be up to 540MB and is only 0.05% smaller file in file tested
    props.fb    = 273;

    CLzmaEncHandle lzma_handle = LzmaEnc_Create (&memfuncs);
    ASSERT0 (lzma_handle, "Error: LzmaEnc_Create failed");

    SRes res = LzmaEnc_SetProps (lzma_handle, &props);
    ASSERT (res == SZ_OK, "Error: LzmaEnc_SetProps failed: %s", lzma_errstr (res));
    
    // write encoding properties as first 5 bytes of compressed stream
    SizeT props_size = LZMA_PROPS_SIZE; // per documentation in LzmaLib.h
    res = LzmaEnc_WriteProperties (lzma_handle, (uint8_t*)compressed, &props_size);
    ASSERT (res == SZ_OK && props_size==LZMA_PROPS_SIZE, "Error: LzmaEnc_WriteProperties failed: %s", lzma_errstr (res));

    bool success = true;

    uint64_t data_compressed_len64 = (uint64_t)*compressed_len - LZMA_PROPS_SIZE;
    res = LzmaEnc_MemEncode (lzma_handle, 
                             (uint8_t *)compressed + LZMA_PROPS_SIZE, &data_compressed_len64, 
                             (uint8_t *)uncompressed, uncompressed_len,
                             true, NULL, &memfuncs, &memfuncs);
    
    *compressed_len = (uint32_t)data_compressed_len64 + LZMA_PROPS_SIZE;

    if (res == SZ_ERROR_OUTPUT_EOF && soft_fail)  // data_compressed_len is too small
        success = false;
    else
        ASSERT (res == SZ_OK, "Error: LzmaEnc_MemEncode failed: %s", lzma_errstr (res));

    LzmaEnc_Destroy (lzma_handle, &memfuncs, &memfuncs);

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
    ASSERT0 (uncompressed, "comp_compress_none: only option 1 supported in the mean time");

    if (*compressed_len < uncompressed_len && soft_fail) return false;
    ASSERT0 (*compressed_len >= uncompressed_len, "Error in comp_compress_none: compressed_len too small");
    
    memcpy (compressed, uncompressed, uncompressed_len);

    *compressed_len = uncompressed_len;

    return true;
}

// compresses data - either a conitguous block or one line at a time. If both are NULL that there is no data to compress.
void comp_compress (VBlock *vb, Buffer *z_data, bool is_z_file_buf,
                    SectionHeader *header, 
                    const char *uncompressed_data, // option 1 - compress contiguous data
                    CompGetLineCallback callback)  // option 2 - compress data one line at a time
{ 
    static Compressor compressors[NUM_COMPRESSOR_ALGS] = { comp_compress_none, comp_compress_bzlib, comp_compress_lzma };

    unsigned compressed_offset     = BGEN32 (header->compressed_offset);
    unsigned data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    unsigned data_compressed_len   = 0;
    unsigned data_encrypted_len=0, data_padding=0, header_padding=0;

    bool is_encrypted = false;
    unsigned encryption_padding_reserve = 0;

    if (header->section_type != SEC_GENOZIP_HEADER) { // genozip header is never encrypted
        is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
        encryption_padding_reserve = crypt_max_padding_len(); // padding for the body
    }

    // if there's no data to compress, set compression to NONE
    if (!data_uncompressed_len) header->data_compression_alg = COMPRESS_NONE;
       
    uint32_t est_compressed_len = 
        (header->data_compression_alg != COMPRESS_NONE) ? MAX (data_uncompressed_len / 2, 500) : data_uncompressed_len;

    // allocate what we think will be enough memory. usually this alloc does nothing, as the memory we pre-allocate for z_data is sufficient
    // note: its ok for other threads to allocate evb data because we have a special mutex in buffer protecting the 
    // evb buffer list
    buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + est_compressed_len + encryption_padding_reserve, 1.5, z_data->name, z_data->param);

    // compress the data, if we have it...
    if (data_uncompressed_len) {

        ASSERT ((unsigned)header->data_compression_alg < NUM_COMPRESSOR_ALGS, 
                "Error in comp_compress: unrecognized data_compression_alg=%d", header->data_compression_alg);
        
        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the alloc, because z_data is pre-allocated

        bool success = 
            compressors[header->data_compression_alg](vb, 
                                                      uncompressed_data, data_uncompressed_len,
                                                      callback,  
                                                      &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                                      true);

        // if output buffer is too small, increase it, and try again
        if (!success) {
            buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + data_uncompressed_len  + encryption_padding_reserve + 50 /* > BZ_N_OVERSHOOT */, 1,
                       z_data->name ? z_data->name : "z_data", z_data->param);
            
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;

            compressors[header->data_compression_alg](vb, 
                                                      uncompressed_data, data_uncompressed_len,
                                                      callback,  
                                                      &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                                      false);
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

        // note: for SEC_VCF_VB_HEADER we will encrypt at the end of calculating this VB when the index data is
        // known, and we will then update z_data in memory prior to writing the encrypted data to disk
        if (header->section_type != SEC_VCF_VB_HEADER || header->vblock_i == 0 /* terminator vb header */)
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

    // do calculations for --show-content and --show-sections options
    vb->z_section_bytes[header->section_type] += total_z_len;
    vb->z_num_sections [header->section_type] ++;
    
    if (flag_show_headers) zfile_show_header (header, vb->vblock_i ? vb : NULL); // store and print upon about for vb sections, and print immediately for non-vb sections
}

void comp_uncompress (VBlock *vb, CompressorAlg alg, 
                      const char *compressed, uint32_t compressed_len,
                      Buffer *uncompressed)
{
    switch (alg) {

    case COMPRESS_BZLIB: {
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
    case COMPRESS_LZMA: {
        static const ISzAlloc memfuncs = { lzma_alloc, lzma_free };
        ELzmaStatus status;
        uint64_t uncompressed_len64 = (uint64_t)uncompressed->len;
        uint64_t compressed_len64 = (uint64_t)compressed_len - LZMA_PROPS_SIZE; // first 5 bytes in compressed stream are the encoding properties
        
        SRes ret = LzmaDecode ((uint8_t *)uncompressed->data, &uncompressed_len64, 
                               (uint8_t *)compressed + LZMA_PROPS_SIZE, &compressed_len64, 
                               (uint8_t *)compressed, LZMA_PROPS_SIZE, 
                               LZMA_FINISH_END, &status, &memfuncs);

        ASSERT (ret == SZ_OK && status == LZMA_STATUS_FINISHED_WITH_MARK, 
                "Error: LzmaDecode failed: ret=%s status=%s", lzma_errstr (ret), lzma_status (status));

        break;
    }
    case COMPRESS_NONE:
        memcpy (uncompressed->data, compressed, compressed_len);
        break;

    default:
        ABORT ("Error in comp_uncompress: invalid compression algorithm %u", alg);
    }
}
