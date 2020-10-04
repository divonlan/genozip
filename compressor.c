// ------------------------------------------------------------------
//   compressor.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "comp_private.h"
#include "vblock.h"
#include "buffer.h"
#include "endianness.h"
#include "crypt.h"
#include "zfile.h"
#include "strings.h"

// -----------------------------------------------------
// memory functions that serve the compression libraries
// -----------------------------------------------------

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
void *comp_alloc (VBlock *vb, int size, double grow_at_least_factor)
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

void comp_free (VBlock *vb, void *addr)
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

void comp_free_all (VBlock *vb)
{
    for (unsigned i=0; i < NUM_COMPRESS_BUFS ; i++) 
        buf_free (&vb->compress_bufs[i]);
}

// CODEC_NONE (no compression)
static bool comp_compress_none (VBlock *vb, Codec codec,
                                const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                                LocalGetLineCallback callback,                        // option 2 - compress data one line at a tim
                                char *compressed, uint32_t *compressed_len /* in/out */, 
                                bool soft_fail)
{
    if (*compressed_len < uncompressed_len && soft_fail) return false;
    ASSERT0 (*compressed_len >= uncompressed_len, "Error in comp_compress_none: compressed_len too small");

    if (callback) {
        char *next = compressed;
        for (unsigned line_i=0; line_i < vb->lines.len; line_i++) {
            char *start1=0, *start2=0;
            uint32_t len1=0, len2=0;        
            callback (vb, line_i, &start1, &len1, &start2, &len2);
            memcpy (next, start1, len1);
            next += len1;
            memcpy (next, start2, len2);
            next += len2;
        }
    }
    else
        memcpy (compressed, uncompressed, uncompressed_len);

    *compressed_len = uncompressed_len;

    return true;
}

static void comp_uncompress_none (VBlock *vb, 
                                  const char *compressed, uint32_t compressed_len,
                                  char *uncompressed_data, uint64_t uncompressed_len)
{
    memcpy (uncompressed_data, compressed, compressed_len);
}

static bool comp_compress_error (VBlock *vb, Codec codec, const char *uncompressed, uint32_t uncompressed_len, LocalGetLineCallback callback,
                                 char *compressed, uint32_t *compressed_len, bool soft_fail) 
{
    ABORT0 ("Error in comp_compress: Unsupported section compression codecorithm");
    return false;
}


static void comp_uncompress_error (VBlock *vb,
                                   const char *compressed, uint32_t compressed_len,
                                   char *uncompressed_data, uint64_t uncompressed_len)
{
    ABORT0 ("Error in comp_uncompress: Unsupported section compression codecorithm");
}

#define MIN_LEN_FOR_COMPRESSION 90 // less that this size, and compressed size is typically larger than uncompressed size

// compresses data - either a contiguous block or one line at a time. If both are NULL that there is no data to compress.
void comp_compress (VBlock *vb, Buffer *z_data, bool is_z_file_buf,
                    SectionHeader *header, 
                    const char *uncompressed_data,  // option 1 - compress contiguous data
                    LocalGetLineCallback callback)  // option 2 - compress data one line at a time
{ 
    ASSERT0 (!uncompressed_data || !callback, "Error in comp_compress: expecting either uncompressed_data or callback but not both");

    ASSERT0 (BGEN32 (header->magic) == GENOZIP_MAGIC, "Error in comp_compress: corrupt header - bad magic");

    // if the user requested --fast - we always use BZLIB, never LZMA
    if (flag_fast && header->codec == CODEC_LZMA)
        header->codec = CODEC_BZ2;

#   define NA  comp_compress_error
    static Compressor compressors[NUM_CODECS] = { 
        comp_compress_none, NA, comp_compress_bzlib,  comp_compress_lzma, NA, NA, NA, NA, NA, NA, 
        comp_compress_lzma /* acgt */, comp_compress_non_acgt, 
        comp_compress_ht, comp_compress_lzma /* ht_index */, 
        NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA  };
#   undef NA

    ASSERT (header->codec < NUM_CODECS, "Error in comp_compress: unsupported section compressor=%u", header->codec);

    unsigned compressed_offset     = BGEN32 (header->compressed_offset);
    unsigned data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    unsigned data_compressed_len   = 0;
    unsigned data_encrypted_len=0, data_padding=0, header_padding=0;

    ASSERT0 (!data_uncompressed_len || uncompressed_data || callback, "Error in comp_compress: data_uncompressed_len!=0 but neither uncompressed_data nor callback are provided");

    bool is_encrypted = false;
    unsigned encryption_padding_reserve = 0;

    if (header->section_type != SEC_GENOZIP_HEADER &&  // genozip header is never encrypted
        !(header->section_type == SEC_REFERENCE && flag_reference == REF_EXT_STORE)) { // external reference copied over is never encrypted
        is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
        encryption_padding_reserve = crypt_max_padding_len(); // padding for the body
    }

    // if there's no data to compress, or its too small, don't compress (except for HT, as it generates an index too)
    if (data_uncompressed_len < MIN_LEN_FOR_COMPRESSION && header->codec != CODEC_HT) 
        header->codec = CODEC_NONE;

    uint32_t est_compressed_len = 
        (header->codec != CODEC_NONE) ? MAX (data_uncompressed_len / 2, 500) : data_uncompressed_len;

    // allocate what we think will be enough memory. usually this alloc does nothing, as the memory we pre-allocate for z_data is sufficient
    // note: its ok for other threads to allocate evb data because we have a special mutex in buffer protecting the 
    // evb buffer list
    buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + est_compressed_len + encryption_padding_reserve, 1.5, z_data->name, z_data->param);

    // compress the data, if we have it...
    if (data_uncompressed_len) {
        
        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the alloc, because z_data is pre-allocated

        bool success = 
            compressors[header->codec](vb, header->codec, uncompressed_data, data_uncompressed_len,
                                                     callback,  
                                                     &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                                     true);
        comp_free_all (vb); // just in case

        // if output buffer is too small, increase it, and try again
        if (!success) {
            buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + data_uncompressed_len  + encryption_padding_reserve + 50 /* > BZ_N_OVERSHOOT */, 1,
                       z_data->name ? z_data->name : "z_data", z_data->param);
            
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;

            compressors[header->codec](vb, header->codec,
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

        // encrypt the header - we use vb_i, section_type and is_header to generate a different AES key for each section
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
    uint64_t offset=0;
    if (header->section_type != SEC_GENOZIP_HEADER)
        offset = sections_add_to_list (vb, header);

    z_data->len += total_z_len;

    if (flag_show_headers) 
        zfile_show_header (header, vb->vblock_i ? vb : NULL, offset, 'W'); // store and print upon about for vb sections, and print immediately for non-vb sections
}

void comp_uncompress (VBlock *vb, Codec codec, 
                      const char *compressed, uint32_t compressed_len,
                      char *uncompressed_data, uint64_t uncompressed_len)
{
    ASSERT0 (compressed_len, "Error in comp_uncompress: compressed_len=0");

#   define NA  comp_uncompress_error
    static Uncompressor uncompressors[NUM_CODECS] = { 
        comp_uncompress_none, NA, comp_uncompress_bzlib,  comp_uncompress_lzma, NA, NA, NA, NA, NA, NA, 
        comp_uncompress_acgt, comp_uncompress_non_acgt, comp_uncompress_bzlib /* ht */, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA  };
#   undef NA

    uncompressors[codec] (vb, compressed, compressed_len, uncompressed_data, uncompressed_len);

    comp_free_all (vb); // just in case
}

const char *codec_name (Codec codec)
{
    static const char *comp_names[NUM_CODECS] = CODEC_NAMES;

    if (codec >=0 && codec < NUM_CODECS) 
        return comp_names[codec];

    else
        return "BAD!";    
}

/*
static Codec vcf_zip_get_best_gt_compressor (VBlock *vb, Buffer *test_data)
{
    static Codec best_gt_data_compressor = CODEC_UNKNOWN;
    static Buffer compressed = EMPTY_BUFFER; // no thread issues as protected my mutex

    // get best compression codec for gt data - lzma or bzlib - their performance varies considerably with
    // the type of data - with either winning by a big margin
    pthread_mutex_lock (&best_gt_data_compressor_mutex);    

    if (best_gt_data_compressor != CODEC_UNKNOWN) goto finish; // answer already known

    #define TEST_BLOCK_SIZE 100000
    buf_alloc (vb, &compressed, TEST_BLOCK_SIZE+1000, 1, "compressed_data_test", 0);

    uint32_t uncompressed_len = MIN (test_data->len, TEST_BLOCK_SIZE);

    uint32_t bzlib_comp_len = compressed.size;
    comp_compress_bzlib (vb, CODEC_BZ2, test_data->data, uncompressed_len, NULL, compressed.data, &bzlib_comp_len, false);
    
    uint32_t lzma_comp_len = compressed.size;
    comp_compress_lzma (vb, CODEC_LZMA, test_data->data, uncompressed_len, NULL, compressed.data, &lzma_comp_len, false);
    
    if      (bzlib_comp_len < uncompressed_len && bzlib_comp_len < lzma_comp_len) best_gt_data_compressor = CODEC_BZ2;
    else if (lzma_comp_len  < uncompressed_len && lzma_comp_len < bzlib_comp_len) best_gt_data_compressor = CODEC_LZMA;
    else                                                                          best_gt_data_compressor = CODEC_NONE;

    buf_free (&compressed);

finish:
    pthread_mutex_unlock (&best_gt_data_compressor_mutex);
    return best_gt_data_compressor;
}
*/