// ------------------------------------------------------------------
//   compressor.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "endianness.h"
#include "crypt.h"
#include "zfile.h"
#include "strings.h"
#include "compressor.h"

// compresses data - either a contiguous block or one line at a time. If both are NULL that there is no data to compress.
// returns data_compressed_len
uint32_t comp_compress (VBlock *vb, Buffer *z_data, bool is_z_file_buf,
                        SectionHeader *header, 
                        const char *uncompressed_data,  // option 1 - compress contiguous data
                        LocalGetLineCB callback)  // option 2 - compress data one line at a time
{ 
    ASSERT0 (!uncompressed_data || !callback, "Error in comp_compress: expecting either uncompressed_data or callback but not both");

    ASSERT0 (BGEN32 (header->magic) == GENOZIP_MAGIC, "Error in comp_compress: corrupt header - bad magic");

    ASSERT (header->codec < NUM_CODECS, "Error in comp_compress: unsupported section compressor=%u", header->codec);

    unsigned compressed_offset     = BGEN32 (header->compressed_offset);
    unsigned data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    unsigned data_compressed_len   = 0;
    unsigned data_encrypted_len=0, data_padding=0, header_padding=0;

    ASSERT0 (!data_uncompressed_len || uncompressed_data || callback, "Error in comp_compress: data_uncompressed_len!=0 but neither uncompressed_data nor callback are provided");

    bool is_encrypted = false;
    unsigned encryption_padding_reserve = 0;

    if (header->section_type != SEC_GENOZIP_HEADER &&  // genozip header is never encrypted
        !(header->section_type == SEC_REFERENCE && flag.reference == REF_EXT_STORE)) { // external reference copied over is never encrypted
        is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
        encryption_padding_reserve = crypt_max_padding_len(); // padding for the body
    }

    // if there's no data to compress, or its too small, and its a simple codec - don't compress
    if (data_uncompressed_len < MIN_LEN_FOR_COMPRESSION && 
        codec_args[header->codec].is_simple) 
        header->codec = CODEC_NONE;

    // use codec's compress function, but if its marked as USE_SUBCODEC, then use sub_codec instead
    Codec est_size_codec = codec_args[header->codec].est_size ? header->codec : codec_args[header->codec].sub_codec;

    uint32_t est_compressed_len = codec_args[est_size_codec].est_size (header->codec, data_uncompressed_len); 

    // allocate what we think will be enough memory. usually this alloc does nothing, as the memory we pre-allocate for z_data is sufficient
    // note: its ok for other threads to allocate evb data because we have a special mutex in buffer protecting the 
    // evb buffer list
    buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + est_compressed_len + encryption_padding_reserve, 1.5, z_data->name);

    // use codec's compress function, but if its marked as USE_SUBCODEC, then use sub_codec instead
    Codec comp_codec = codec_args[header->codec].compress ? header->codec : codec_args[header->codec].sub_codec;

    // compress the data, if we have it...
    if (data_uncompressed_len) {
        
        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the alloc, because z_data is pre-allocated

        bool success = 
            codec_args[comp_codec].compress (vb, header, uncompressed_data, &data_uncompressed_len,
                                             callback,  
                                             &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                             true);
        codec_free_all (vb); // just in case

        // if output buffer is too small, increase it, and try again
        if (!success) {
            buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + data_uncompressed_len  + encryption_padding_reserve + 50 /* > BZ_N_OVERSHOOT, LIBBSC_HEADER_SIZE */, 1,
                       z_data->name ? z_data->name : "z_data");
            
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;
            data_uncompressed_len = BGEN32 (header->data_uncompressed_len); // reset

            codec_args[comp_codec].compress (vb, header,
                                             uncompressed_data, &data_uncompressed_len,
                                             callback,  
                                             &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                             false);

            codec_free_all (vb); // just in case
        }
    
        // update uncompressed length - complex codecs (like domqual) might change it
        header->data_uncompressed_len = BGEN32 (data_uncompressed_len);
        
        // get encryption related lengths
        if (is_encrypted) {
            data_encrypted_len = data_compressed_len;
            crypt_get_encrypted_len (&data_encrypted_len, &data_padding); // both are set to 0 if no encryption
        }
    }

    // if there is no compressed data, then no need to write this section (this can happen when
    // novel codecs create other sections, eg CODEC_GTSHARK)
    if (data_uncompressed_len && !data_compressed_len) goto done;

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

    if (flag.show_headers) 
        zfile_show_header (header, vb->vblock_i ? vb : NULL, offset, 'W'); // store and print upon about for vb sections, and print immediately for non-vb sections

done:
    return data_compressed_len;
}

void comp_uncompress (VBlock *vb, Codec codec, Codec sub_codec,
                      const char *compressed, uint32_t compressed_len,
                      Buffer *uncompressed_data, uint64_t uncompressed_len)
{
    ASSERT0 (compressed_len, "Error in comp_uncompress: compressed_len=0");

    if (codec && codec_args[codec].uncompress) // might be UNKNOWN (eg GT_X_ALLELES) or not have an uncompressor (eg: HT, DOMQ don't have a special uncompressor) and only a sub-codec available
        codec_args[codec].uncompress (vb, codec, compressed, compressed_len, uncompressed_data, uncompressed_len, sub_codec);

    // if this is (1) a simple codec (including CODEC_UNKNOWN) that has a sub-codec or
    // (2) or no codec uncompressor - the sub-codec now run it now
    // note: for non-simple codecs, the sub-codec is used in the codec decompression
    if (sub_codec && (codec_args[codec].is_simple || !codec_args[codec].uncompress))
        codec_args[sub_codec].uncompress (vb, sub_codec, compressed, compressed_len, uncompressed_data, uncompressed_len, CODEC_UNKNOWN);
}

void comp_unit_test (SectionHeader *header)
{
    uint32_t size = 1000000;
    //int size = 3380084;
    char *data = MALLOC (size);
    for (int i=0; i < size; i++) data[i] = 'A' + (i%26);

    //FILE *fp = fopen ("bugq.bz2", "rb");
    //ASSERT0 (fread (data, 1, size, fp) == size, "read failed");

    uint32_t comp_len = codec_args[header->codec].est_size (header->codec, size);
    char *comp = MALLOC (comp_len);
    
    codec_args[header->codec].compress (evb, header, data, &size, 0, comp, &comp_len, false);

    Buffer uncomp = EMPTY_BUFFER;
    buf_alloc (evb, &uncomp, size, 1, "uncomp");
    codec_args[header->codec].uncompress (evb, header->codec, comp, comp_len, &uncomp, size, CODEC_NONE);

    printf ("Unit test %s!\n", memcmp (data, uncomp.data, size) ? "failed" : "succeeded");

    exit(0);
}
