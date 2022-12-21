// ------------------------------------------------------------------
//   compressor.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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
#include "libdeflate/libdeflate.h"

// compresses data - either a contiguous block or one line at a time. If both are NULL that there is no data to compress.
// returns data_compressed_len
uint32_t comp_compress (VBlockP vb, 
                        ContextP ctx,            // note: NULL if not compressing context
                        BufferP z_data,
                        SectionHeaderP header, 
                        rom uncompressed_data,   // option 1 - compress contiguous data
                        LocalGetLineCB callback, // option 2 - compress data one line at a time
                        rom name)                // tag name for b250/local, otherwise section name
{ 
    ASSERT (!uncompressed_data || !callback, "\"%s\": expecting either uncompressed_data or callback but not both", name);

    ASSERT (BGEN32 (header->magic) == GENOZIP_MAGIC, "\"%s\": corrupt header - bad magic", name);

    ASSERT (header->codec < NUM_CODECS, "\"%s\": unsupported section compressor=%u", name, header->codec);
             
    uint32_t compressed_offset     = BGEN32 (header->compressed_offset);
    uint32_t data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    uint32_t data_compressed_len   = 0;
    uint32_t data_encrypted_len=0, data_padding=0, header_padding=0;

    ASSERT (!data_uncompressed_len || uncompressed_data || callback, "\"%s\": data_uncompressed_len!=0 but neither uncompressed_data nor callback are provided", name);

    ASSERTW (data_uncompressed_len < (1<<30), "%s: Excessive uncompressed_data_len=%u: %s. Please report to support@genozip.com", 
             VB_NAME, data_uncompressed_len, name); // compressing a buffer over 1GB is likely an indication of not handling some edge case well

    bool is_encrypted = false;
    uint32_t encryption_padding_reserve = 0;

    if (header->section_type != SEC_GENOZIP_HEADER &&  // genozip header is never encrypted
        !(header->section_type == SEC_REFERENCE && IS_REF_EXT_STORE)) { // external reference copied over is never encrypted
        is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
        encryption_padding_reserve = crypt_max_padding_len(); // padding for the body
    }

    // if there's no data to compress, or its too small, and its a simple codec - don't compress
    if (data_uncompressed_len < MIN_LEN_FOR_COMPRESSION && 
        codec_args[header->codec].is_simple) 
        header->codec = CODEC_NONE;

    // use codec's compress function, but if its marked as USE_SUBCODEC, then use sub_codec instead
    Codec comp_codec = (codec_args[header->codec].compress != USE_SUBCODEC) ? header->codec : header->sub_codec;

    uint32_t est_compressed_len = codec_args[comp_codec].est_size (header->codec, data_uncompressed_len); 

    // allocate what we think will be enough memory. usually this alloc does nothing, as the memory we pre-allocate for z_data is sufficient
    buf_alloc (vb, z_data, compressed_offset + est_compressed_len + encryption_padding_reserve, 0, char, 1.5, 
               z_data->name ? z_data->name : "z_data");

    // compress the data, if we have it...
    if (data_uncompressed_len) {

        // flag.verify-codec: add adler32 of the uncompressed data to the Header, instead of magic
        if (flag.verify_codec && uncompressed_data && data_uncompressed_len && !in_assign_codec)  // TODO: support the callback option
            header->uncomp_adler32 = adler32 (1, uncompressed_data, data_uncompressed_len);

        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the alloc, because z_data is pre-allocated

        codec_verify_free_all (vb, "compressor", comp_codec);

        vb->codec_using_codec_bufs = comp_codec;

        bool success = 
            codec_args[comp_codec].compress (vb, ctx, header, uncompressed_data, &data_uncompressed_len,
                                             callback,  
                                             BAFTc (*z_data) + compressed_offset, &data_compressed_len,
                                             true, name);
        codec_free_all (vb); // just in case

        // if output buffer is too small, increase it, and try again
        if (!success) {

            uint32_t at_least = header->sub_codec ? codec_args[header->sub_codec].est_size (header->sub_codec, data_uncompressed_len)
                                                  : MAX_(data_uncompressed_len * 1.5 + 500000,
                                                         codec_args[header->codec].est_size (header->sub_codec, data_uncompressed_len));

            buf_alloc (vb, z_data, compressed_offset + at_least, 0, char, 1, z_data->name ? z_data->name : "z_data");
            
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;
            data_uncompressed_len = BGEN32 (header->data_uncompressed_len); // reset

            codec_verify_free_all (vb, "compressor", comp_codec);

            codec_args[comp_codec].compress (vb, ctx, header,
                                             uncompressed_data, &data_uncompressed_len,
                                             callback,  
                                             &z_data->data[z_data->len + compressed_offset], &data_compressed_len,
                                             false, name);

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
    memcpy (BAFT(uint8_t, *z_data), header, compressed_offset);

    // encrypt if needed - header & body separately
    unsigned total_z_len;
    if (is_encrypted) { // create good padding (just padding with 0 exposes a cryptographic vulnerability)
        crypt_pad (BAFT(uint8_t, *z_data), compressed_offset, header_padding);
        
        if (data_uncompressed_len) 
            crypt_pad (BAFT(uint8_t, *z_data) + compressed_offset, data_encrypted_len, data_padding);

        // encrypt the header - we use vb_i, section_type and is_header to generate a different AES key for each section
        VBIType vb_i = BGEN32 (header->vblock_i);

        // note: for SEC_VB_HEADER we will encrypt at the end of calculating this VB in zfile_update_compressed_vb_header() 
        // and we will then update z_data in memory prior to writing the encrypted data to disk
        if (header->section_type != SEC_VB_HEADER || header->vblock_i == 0 /* terminator vb header */)
            crypt_do (vb, BAFT(uint8_t, *z_data), compressed_offset, vb_i, header->section_type, true);

        // encrypt the data body 
        if (data_uncompressed_len) 
            crypt_do (vb, BAFT(uint8_t, *z_data) + compressed_offset, data_encrypted_len, vb_i, header->section_type, false);
        
        total_z_len = compressed_offset + data_encrypted_len;
    }
    else
        total_z_len = compressed_offset + data_compressed_len;

    // add section to the section list
    if (header->section_type != SEC_GENOZIP_HEADER && // the section list is part of the genozip header section currently being compressed, so can't add to it
        !in_assign_codec)
        sections_add_to_list (vb, header);

    if (flag.show_headers && !in_assign_codec)
        sections_show_header (header, vb->vblock_i ? vb : NULL, z_data->len, 'W'); // store and print upon about for vb sections, and print immediately for non-vb sections

    z_data->len += total_z_len;

    // if we're compressing a global buffer in the main thread, we can write it immeidately
    if (vb == evb && !in_assign_codec) 
        zfile_output_processed_vb (vb);

done:
    return data_compressed_len;
}

// compress primary context of a complex codec, after codec code as prepared the data in ctx->local. The other contexts 
// of the complex codec are marked with DEP_L* and will be compressed in the normal ctx->local compression loop
bool comp_compress_complex_codec (VBlockP vb, ContextP ctx, SectionHeaderP header, bool is_2nd_try,
                                  uint32_t *uncompressed_len, STRe (compressed), rom name)
{
    if (!is_2nd_try) {
        Codec save_lcodec = ctx->lcodec;
        ctx->lcodec = CODEC_UNKNOWN;
        header->sub_codec = codec_assign_best_codec (vb, ctx, &ctx->local, SEC_LOCAL); // provide BufferP to override callback
        if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small
        ctx->lcodec = save_lcodec;
    }

    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = ctx->local.len32;

    // make sure we have enough memory
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, ctx->local.len);
    if (*compressed_len < min_required_compressed_len) {
        if (!is_2nd_try) return false; // call me again with more memory
        ABORT ("%s: Compressing %s with codec %s + sub-codec %s need %u bytes, but allocated only %u", 
                VB_NAME, ctx->tag_name, codec_name(ctx->lcodec), codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);
    }

    return compress (vb, ctx, header, ctx->local.data, uncompressed_len, NULL, compressed, compressed_len, false, name);
}

void comp_uncompress (VBlockP vb, Codec codec, Codec sub_codec, uint8_t param, STRp(compressed),
                      BufferP uncompressed_data, uint64_t uncompressed_len, rom name)
{
    ASSERTNOTZERO (compressed_len, name);
    ASSERTNOTZERO (uncompressed_len, name);
    ASSERT (codec != sub_codec, "vb=%s \"%s\": Expectedly, codec=%s == sub_codec=%s", VB_NAME, name, codec_name(codec), codec_name(sub_codec));

    // if this is (1) a simple codec (including CODEC_UNKNOWN) that has a sub-codec or
    // (2) or no codec uncompressor - run the sub-codec uncompressor directly
    // for non-simple codecs, the subcodec is handled within the main codec decompressor. run_subcodec is true if we need to handle it here.
    bool run_subcodec = sub_codec && (codec_args[codec].is_simple || !codec_args[codec].uncompress);

    if (codec && codec_args[codec].uncompress) { // might be UNKNOWN (eg GT_X_ALLELES) or not have an uncompressor (eg: HT, DOMQ don't have a special uncompressor) and only a sub-codec available

        codec_verify_free_all (vb, "uncompressor", codec);

        vb->codec_using_codec_bufs = codec;
        codec_args[codec].uncompress (vb, codec, param, STRa(compressed), uncompressed_data, uncompressed_len, sub_codec, name);
        codec_free_all (vb); // just in case
        vb->codec_using_codec_bufs = CODEC_UNKNOWN;

        // case: uncompressed data is now going to be the compressed data for the sub-codec
        if (run_subcodec) {
            // move the intermediary data after primary codec decompression to vb->scratch
            buf_copy (vb, &vb->scratch, uncompressed_data, char, 0,0, "scratch");
            compressed     = vb->scratch.data;
            compressed_len = vb->scratch.len;
            
            // free uncompressed_data to accept the uncompressed data after subcodec decompression
            buf_free (*uncompressed_data);
            uncompressed_len = 0; // the header doesn't say the uncompressed len after the subcodec - the subcodec should know this using its own means
        }
    }

    if (run_subcodec) { // might run with or without first running the primary codec
        codec_verify_free_all (vb, "uncompressor", sub_codec);
        vb->codec_using_codec_bufs = sub_codec;
        codec_args[sub_codec].uncompress (vb, sub_codec, param, STRa(compressed), uncompressed_data, uncompressed_len, CODEC_UNKNOWN, name);
        codec_free_all (vb); // just in case
        vb->codec_using_codec_bufs = CODEC_UNKNOWN;
    }
}

