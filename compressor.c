// ------------------------------------------------------------------
//   compressor.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
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
#include "zfile.h"

// compresses data - either a contiguous block or one line at a time. If both are NULL that there is no data to compress.
// returns data_compressed_len
uint32_t comp_compress (VBlock *vb, Buffer *z_data,
                        SectionHeader *header, 
                        const char *uncompressed_data,  // option 1 - compress contiguous data
                        LocalGetLineCB callback)        // option 2 - compress data one line at a time
{ 
    ASSERT0 (!uncompressed_data || !callback, "expecting either uncompressed_data or callback but not both");

    ASSERT0 (BGEN32 (header->magic) == GENOZIP_MAGIC, "corrupt header - bad magic");

    ASSERT (header->codec < NUM_CODECS, "unsupported section compressor=%u", header->codec);

    unsigned compressed_offset     = BGEN32 (header->compressed_offset);
    unsigned data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    unsigned data_compressed_len   = 0;
    unsigned data_encrypted_len=0, data_padding=0, header_padding=0;

    ASSERT0 (!data_uncompressed_len || uncompressed_data || callback, "data_uncompressed_len!=0 but neither uncompressed_data nor callback are provided");

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
    buf_alloc_old (vb, z_data, z_data->len + compressed_offset + est_compressed_len + encryption_padding_reserve, 1.5, 
               z_data->name ? z_data->name : "z_data");

    // use codec's compress function, but if its marked as USE_SUBCODEC, then use sub_codec instead
    Codec comp_codec = (codec_args[header->codec].compress != USE_SUBCODEC) ? header->codec : codec_args[header->codec].sub_codec;

    // compress the data, if we have it...
    if (data_uncompressed_len) {
        
        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the alloc, because z_data is pre-allocated

        codec_verify_free_all (vb, "compressor", comp_codec);

        vb->codec_using_codec_bufs = comp_codec;
        bool success = 
            codec_args[comp_codec].compress (vb, header, uncompressed_data, &data_uncompressed_len,
                                             callback,  
                                             &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                             true);
        codec_free_all (vb); // just in case

        // if output buffer is too small, increase it, and try again
        if (!success) {
            buf_alloc_old (vb, z_data, z_data->len + compressed_offset + data_uncompressed_len * 1.5 + encryption_padding_reserve + 50 /* > BZ_N_OVERSHOOT, LIBBSC_HEADER_SIZE */, 1,
                       z_data->name ? z_data->name : "z_data");
            
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;
            data_uncompressed_len = BGEN32 (header->data_uncompressed_len); // reset

            codec_verify_free_all (vb, "compressor", comp_codec);

            codec_args[comp_codec].compress (vb, header,
                                             uncompressed_data, &data_uncompressed_len,
                                             callback,  
                                             &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                             false);

            codec_free_all (vb); // just in case
        }
        vb->codec_using_codec_bufs = CODEC_UNKNOWN;
    
        // update uncompressed length - complex codecs (like domqual, pbwt) might change it
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

    // add section to the section list
    if (header->section_type != SEC_GENOZIP_HEADER && // the section list is part of the genozip header section currently being compressed, so can't add to it
        header->section_type != SEC_NONE /* from codec_assign_best_codec */)
        sections_add_to_list (vb, header);

    if (flag.show_headers &&
        header->section_type != SEC_NONE) // not only testing from codec_assign_best_codec
        zfile_show_header (header, vb->vblock_i ? vb : NULL, z_data->len, 'W'); // store and print upon about for vb sections, and print immediately for non-vb sections

    z_data->len += total_z_len;

    // if we're compressing a global buffer in the main thread, we can write it immeidately
    if (vb == evb && header->section_type != SEC_NONE) 
        zfile_output_processed_vb (vb);

done:
    return data_compressed_len;
}

void comp_uncompress (VBlock *vb, Codec codec, Codec sub_codec, uint8_t param,
                      const char *compressed, uint32_t compressed_len,
                      Buffer *uncompressed_data, uint64_t uncompressed_len)
{
    ASSERT0 (compressed_len, "compressed_len=0");

    // if this is (1) a simple codec (including CODEC_UNKNOWN) that has a sub-codec or
    // (2) or no codec uncompressor - the sub-codec now run it now
    // for non-simple codecs, the subcodec is handled within the main codec decompressor. run_subcodec is true if we need to handle it here.
    bool run_subcodec = sub_codec && (codec_args[codec].is_simple || !codec_args[codec].uncompress);

    if (codec && codec_args[codec].uncompress) { // might be UNKNOWN (eg GT_X_ALLELES) or not have an uncompressor (eg: HT, DOMQ don't have a special uncompressor) and only a sub-codec available

        codec_verify_free_all (vb, "uncompressor", codec);

        vb->codec_using_codec_bufs = codec;
        codec_args[codec].uncompress (vb, codec, param, compressed, compressed_len, uncompressed_data, uncompressed_len, sub_codec);
        codec_free_all (vb); // just in case
        vb->codec_using_codec_bufs = CODEC_UNKNOWN;

        // case: uncompressed data is now going to be the compressed data for the sub-codec
        if (run_subcodec) {
            // move the intermediary data after primary codec decompression to vb->compressed
            buf_copy (vb, &vb->compressed, uncompressed_data, char, 0,0, "compressed");
            compressed     = vb->compressed.data;
            compressed_len = vb->compressed.len;
            
            // free uncompressed_data to accept the uncompressed data after subcodec decompression
            buf_free (uncompressed_data);
            uncompressed_len = 0; // the header doesn't say the uncompressed len after the subcodec - the subcodec should know this using its own means
        }
    }

    if (run_subcodec) { // might run with or without first running the primary codec
        codec_verify_free_all (vb, "uncompressor", sub_codec);
        vb->codec_using_codec_bufs = sub_codec;
        codec_args[sub_codec].uncompress (vb, sub_codec, param, compressed, compressed_len, uncompressed_data, uncompressed_len, CODEC_UNKNOWN);
        codec_free_all (vb); // just in case
        vb->codec_using_codec_bufs = CODEC_UNKNOWN;
    }
}

