// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#include <time.h>
#include <math.h>
#include <bzlib.h>
#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "squeeze.h"
#include "vcf_header.h"
#include "crypt.h"
#include "vb.h"
#include "move_to_front.h"
#include "file.h"
#include "endianness.h"
#include "version.h"
#include "sections.h"
#include "gtshark.h"

#define BZLIB_BLOCKSIZE100K 9 /* maximum mem allocation for bzlib */

static Buffer compress_bufs[NUM_COMPRESS_BUFS] = {EMPTY_BUFFER, EMPTY_BUFFER, EMPTY_BUFFER, EMPTY_BUFFER}; // used by I/O thread when compressing out of context of a VB (i.e. VCF header)

static const char *password_test_string = "WhenIThinkBackOnAllTheCrapIlearntInHighschool";

// memory management for bzlib - tesing shows that compress allocates 4 times, and decompress 2 times. Allocations are the same set of sizes
// every call to compress/decompress with the same parameters, independent on the contents or size of the compressed/decompressed data.
static void *zfile_bzalloc (void *vb_, int items, int size)
{
    VariantBlock *vb = (VariantBlock *)vb_;
    Buffer *bufs = vb ? vb->compress_bufs : compress_bufs;

    Buffer *use_buf = NULL;
    
    // search for a free buffer - preferring one of the size requested 
   for (unsigned i=0; i < NUM_COMPRESS_BUFS ; i++) {
        if (!buf_is_allocated (&bufs[i])) {
            use_buf = &bufs[i];
            if (use_buf->size == (unsigned)(items * size)) break;
        }
    }
    ASSERT0 (use_buf, "Error: zfile_bzalloc could not find a free buffer");

    buf_alloc (vb, use_buf, items * size, 1, "compress_bufs", 0);

    return use_buf->data;
}

static void zfile_bzfree (void *vb_, void *addr)
{
    VariantBlock *vb = (VariantBlock *)vb_;
    Buffer *bufs = vb ? vb->compress_bufs : compress_bufs;

    unsigned i; for (i=0; i < NUM_COMPRESS_BUFS ; i++) 
        if (bufs[i].data == addr) {
            buf_free (&bufs[i]);
            break;
        }

    ASSERT0 (i < NUM_COMPRESS_BUFS, "Error: zfile_bzfree failed to find buffer to free");
}

static int zfile_compress_do (VariantBlock *vb, Buffer *z_data, const char *data,
                              unsigned compressed_offset, unsigned data_uncompressed_len, 
                              unsigned *data_compressed_len /* in/out */, bool ok_is_ok)
{
    START_TIMER;

    bz_stream strm;
    strm.bzalloc = zfile_bzalloc;
    strm.bzfree  = zfile_bzfree;
    strm.opaque  = vb; // just passed to malloc/free

    int ret = BZ2_bzCompressInit (&strm, BZLIB_BLOCKSIZE100K, 0, 30);
    ASSERT (ret == BZ_OK, "Error: BZ2_bzCompressInit failed: %s", BZ2_errstr(ret));

    strm.next_in   = (char*)data;
    strm.next_out  = &z_data->data[z_data->len + compressed_offset];
    strm.avail_in  = data_uncompressed_len;
    strm.avail_out = *data_compressed_len;

    ret = BZ2_bzCompress (&strm, BZ_FINISH);
    ASSERT ((ok_is_ok && ret == BZ_FINISH_OK) || ret == BZ_STREAM_END, "Error: BZ2_bzCompress failed: %s", BZ2_errstr(ret));

    BZ2_bzCompressEnd (&strm);

    *data_compressed_len -= strm.avail_out;

    if (vb) COPY_TIMER(vb->profile.compressor);

    return ret;
}

static void zfile_show_header (const SectionHeader *header, VariantBlock *vb /* optional if output to buffer */)
{
    DictIdType dict_id = {0};

    if (section_type_is_dictionary (header->section_type)) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    if (section_type_is_b250       (header->section_type)) dict_id = ((SectionHeaderBase250    *)header)->dict_id;

    char str[1000];

    sprintf (str, "%-23s %*.*s vb_i=%-3u sec_i=%-2u comp_offset=%-6u uncomp_len=%-6u comp_len=%-6u enc_len=%-6u magic=%8.8x\n",
             st_name(header->section_type), -DICT_ID_LEN, DICT_ID_LEN, dict_id.num ? dict_id_printable (dict_id).id : dict_id.id,
             BGEN32 (header->variant_block_i), BGEN16 (header->section_i), 
             BGEN32 (header->compressed_offset), BGEN32 (header->data_uncompressed_len),
             BGEN32 (header->data_compressed_len), BGEN32 (header->data_encrypted_len), BGEN32 (header->magic));

    if (vb) {
        unsigned len = strlen (str);
        buf_alloc (vb, &vb->show_headers_buf, vb->show_headers_buf.len + len + 1, 2, "show_headers_buf", 0);
        strcpy (&vb->show_headers_buf.data[vb->show_headers_buf.len], str);
        vb->show_headers_buf.len += len;
    }
    else 
        fprintf (stderr, "%s", str);
}

static void zfile_compress (VariantBlock *vb, Buffer *z_data, bool is_z_file_buf,
                            SectionHeader *header, const char *data /* NULL if section has no data */) 
{ 
    unsigned compressed_offset     = BGEN32 (header->compressed_offset);
    unsigned data_uncompressed_len = BGEN32 (header->data_uncompressed_len);
    unsigned data_compressed_len   = 0;
    unsigned data_encrypted_len=0, data_padding=0, header_padding=0;

    // for DB_DB and DB_GT we don't compress as they are already compressed
    bool need_compression = (header->section_type != SEC_HT_GTSHARK_DB_DB && header->section_type != SEC_HT_GTSHARK_DB_GT);

    bool is_encrypted = false;
    unsigned encryption_padding_reserve = 0;

    if (header->section_type != SEC_GENOZIP_HEADER) { // genozip header is never encrypted
        is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
        encryption_padding_reserve = crypt_max_padding_len(); // padding for the body
    }

    // allocate what we think will be enough memory. usually this realloc does nothing, as the
    // memory we pre-allocate for z_data is sufficient
    if (need_compression)
        buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + MAX (data_uncompressed_len / 2, 500) + encryption_padding_reserve, // usually compression is a lot better than 2X so this should be enough
                   1.5, z_data->name, z_data->param); // name and param need to be set in advance
    else 
        buf_alloc (vb, z_data, z_data->len + compressed_offset + data_uncompressed_len, 1.5, z_data->name, z_data->param);

    // compress the data, if we have it...
    if (data_uncompressed_len) {
     
        if (need_compression) {

            // compress data
            data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the realloc, because z_data is pre-allocated

            int ret = zfile_compress_do (vb, z_data, data, compressed_offset, data_uncompressed_len, &data_compressed_len, true);

            // if output buffer is too small, increase it, and try again
            if (ret == BZ_FINISH_OK) {
                buf_alloc (is_z_file_buf ? evb : vb, z_data, z_data->len + compressed_offset + data_uncompressed_len + 50 /* > BZ_N_OVERSHOOT */, 1,
                        z_data->name ? z_data->name : "z_data", z_data->param);
                data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;

                zfile_compress_do (vb, z_data, data, compressed_offset, data_uncompressed_len, &data_compressed_len, false);
            }
        }
        else { // no compression
            data_compressed_len = data_uncompressed_len;
            memcpy (&z_data->data[z_data->len + compressed_offset], data, data_uncompressed_len);
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
        uint32_t vb_i  = BGEN32 (header->variant_block_i);

        // note: for SEC_VB_HEADER we will encrypt at the end of calculating this VB when the index data is
        // known, and we will then update z_data in memory prior to writing the encrypted data to disk
        if (header->section_type != SEC_VB_HEADER || header->variant_block_i == 0 /* terminator vb header */)
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
    
    if (flag_show_headers) zfile_show_header (header, vb->variant_block_i ? vb : NULL); // store and print upon about for vb sections, and print immediately for non-vb sections
}

static void zfile_show_b250_section (void *section_header_p, Buffer *b250_data)
{
    SectionHeaderBase250 *header = (SectionHeaderBase250 *)section_header_p;

    if (!flag_show_b250 && dict_id_printable (header->dict_id).num != dict_id_show_one_b250.num) return;

    if (flag_show_b250 && header->h.section_type == SEC_CHROM_B250)
        fprintf (stderr, "Base-250 data for VB %u (result of '--show-b250'):\n", BGEN32 (header->h.variant_block_i));
        
    fprintf (stderr, "  %*.*s: ", -DICT_ID_LEN-1, DICT_ID_LEN, dict_id_printable (header->dict_id).id);

    const uint8_t *data = (const uint8_t *)b250_data->data;
    for (unsigned i=0; i < BGEN32 (header->num_b250_items); i++) {
        uint32_t word_index = base250_decode (&data);
        switch (word_index) {
            case WORD_INDEX_ONE_UP     : fprintf (stderr, "ONE_UP "); break;
            case WORD_INDEX_EMPTY_SF   : fprintf (stderr, "EMPTY "); break;
            case WORD_INDEX_MISSING_SF : fprintf (stderr, "MISSING "); break;
            default: fprintf (stderr, "%u ", word_index);
        }
    }
    fprintf (stderr, "\n");
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_read_one_section
void zfile_uncompress_section (VariantBlock *vb,
                               void *section_header_p,
                               Buffer *uncompressed_data,
                               const char *uncompressed_data_buf_name,
                               SectionType expected_section_type) 
{
    START_TIMER;

    // for DB_DB and DB_GT are not compressed by us as they are already compressed by gtshark
    bool need_decompression = (expected_section_type != SEC_HT_GTSHARK_DB_DB && expected_section_type != SEC_HT_GTSHARK_DB_GT);

    SectionHeader *section_header = (SectionHeader *)section_header_p;
    
    uint32_t compressed_offset     = BGEN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = BGEN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = BGEN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = BGEN32 (section_header->data_uncompressed_len);
    uint32_t variant_block_i       = BGEN32 (section_header->variant_block_i);
    uint16_t section_i             = BGEN16 (section_header->section_i);

    // sanity checks
    ASSERT (section_header->section_type == expected_section_type, "Error: expecting section type %s but seeing %s", st_name(expected_section_type), st_name(section_header->section_type));
    
    bool expecting_vb_i = !section_type_is_dictionary (expected_section_type) && expected_section_type != SEC_VCF_HEADER;
    ASSERT (variant_block_i == vb->variant_block_i || !expecting_vb_i, // dictionaries are uncompressed by the I/O thread with pseduo_vb (vb_i=0) 
             "Error: bad variant_block_i: in file=%u in vb=%u", variant_block_i, vb->variant_block_i);

    // decrypt data (in-place) if needed
    if (data_encrypted_len) {
        if (z_file->genozip_version > 1)
            crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, variant_block_i, section_header->section_type, false);
        else
            v1_crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, variant_block_i, section_i);
    }

    if (data_uncompressed_len > 0) { // FORMAT, for example, can be missing in a sample-less file

        buf_alloc (vb, uncompressed_data, data_uncompressed_len + 1, 1.1, uncompressed_data_buf_name, 0); // +1 for \0
            
        uncompressed_data->len = data_uncompressed_len;

        if (need_decompression) {
            bz_stream strm;
            strm.bzalloc = zfile_bzalloc;
            strm.bzfree  = zfile_bzfree;
            strm.opaque  = vb; // just passed to malloc/free

            int ret = BZ2_bzDecompressInit (&strm, 0, 0);
            ASSERT0 (ret == BZ_OK, "Error: BZ2_bzDecompressInit failed");

            strm.next_in   = (char*)section_header + compressed_offset;
            strm.next_out  = uncompressed_data->data;
            strm.avail_in  = data_compressed_len;
            strm.avail_out = uncompressed_data->len;

            ret = BZ2_bzDecompress (&strm);
            ASSERT (ret == BZ_STREAM_END || ret == BZ_OK, "Error: BZ2_bzDecompress failed: %s, avail_in=%d, avail_out=%d", BZ2_errstr(ret), strm.avail_in, strm.avail_out);

            BZ2_bzDecompressEnd (&strm);
        }
        else 
            memcpy (uncompressed_data->data, (char*)section_header + compressed_offset, data_compressed_len);
    }

    if (flag_show_b250 && section_type_is_b250 (expected_section_type)) 
        zfile_show_b250_section (section_header_p, uncompressed_data);
    
    if (vb) COPY_TIMER(vb->profile.zfile_uncompress_section);
}

// generate the metadata string eg "user@domain.com 2019-11-16 16:48:19"
static void zfile_get_metadata(char *metadata)
{
    time_t timer;
    char time_buf[20];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", tm_info);

    const char *user = getenv ("USER");
    const char *host = getenv ("HOSTNAME");
    sprintf(metadata, "%.19s %.20s%s%.30s", time_buf,  // 71 chars + the string termintor \0 = 72
            user ? user : "",
            user && host ? "@" : "",
            host ? host : "");

    // sanity
    ASSERT0 (strlen (metadata) < FILE_METADATA_LEN, "Error: metadata too long");
}

// ZIP
void zfile_write_vcf_header (Buffer *vcf_header_text, bool is_first_vcf)
{
    SectionHeaderVCFHeader vcf_header;
    memset (&vcf_header, 0, sizeof(vcf_header)); // safety

    vcf_header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    vcf_header.h.section_type          = SEC_VCF_HEADER;
    vcf_header.h.data_uncompressed_len = BGEN32 (vcf_header_text->len);
    vcf_header.h.compressed_offset     = BGEN32 (sizeof (SectionHeaderVCFHeader));
    vcf_header.num_samples             = BGEN32 (global_num_samples);
    vcf_header.vcf_data_size           = BGEN64 (vcf_file->vcf_data_size_single) /* 0 if gzipped - will be updated later. */; 
    vcf_header.num_lines               = NUM_LINES_UNKNOWN; 
    
    switch (vcf_file->type) {
        case VCF     : vcf_header.compression_type = COMPRESSION_TYPE_NONE  ; break;
        case VCF_GZ  :
        case VCF_BGZ : vcf_header.compression_type = COMPRESSION_TYPE_GZIP  ; break;
        case VCF_BZ2 : vcf_header.compression_type = COMPRESSION_TYPE_BZIP2 ; break;
        default      : ABORT ("Error: invalid vcf_file->type=%u", vcf_file->type);
    }

    file_basename (vcf_file->name, false, "(stdin)", vcf_header.vcf_filename, VCF_FILENAME_LEN);
    vcf_header.vcf_filename[strlen(vcf_header.vcf_filename)- (strlen(file_exts[vcf_file->type])-4)] = '\0'; // remove the .gz/.bgz/.bz2
    
    static Buffer vcf_header_buf = EMPTY_BUFFER;

    buf_alloc ((VariantBlock*)evb, &vcf_header_buf, vcf_header_text->len / 3, // generous guess of compressed size
               1, "vcf_header_buf", 0); 

    zfile_compress ((VariantBlock*)evb, &vcf_header_buf, true, (SectionHeader*)&vcf_header, vcf_header_text->data);

    file_write (z_file, vcf_header_buf.data, vcf_header_buf.len);

    z_file->disk_so_far      += vcf_header_buf.len;   // length of GENOZIP data writen to disk
    z_file->vcf_data_so_far  += vcf_header_text->len; // length of the original VCF header

    buf_free (&vcf_header_buf); 

    // copy it to z_file - we might need to update it at the very end in zfile_update_vcf_header_section_header()
    if (is_first_vcf)
        memcpy (&z_file->vcf_header_first, &vcf_header, sizeof (vcf_header));

    memcpy (&z_file->vcf_header_single, &vcf_header, sizeof (vcf_header));
}

void zfile_compress_vb_header (VariantBlock *vb)
{
    uint32_t my_squeeze_len = squeeze_len (vb->num_haplotypes_per_line);
    uint32_t sizeof_header = sizeof (SectionHeaderVbHeader);

    SectionHeaderVbHeader vb_header;
    memset (&vb_header, 0, sizeof(SectionHeaderVbHeader)); // safety
    
    vb_header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    vb_header.h.section_type          = SEC_VB_HEADER;
    vb_header.h.data_uncompressed_len = BGEN32 (my_squeeze_len);
    vb_header.h.compressed_offset     = BGEN32 (sizeof_header);
    vb_header.h.variant_block_i       = BGEN32 (vb->variant_block_i);
    vb_header.h.section_i             = BGEN16 (vb->z_next_header_i++); // always 0
    vb_header.first_line              = BGEN32 (vb->first_line);
    vb_header.num_lines               = BGEN32 (vb->num_lines);
    vb_header.phase_type              = (char)vb->phase_type; 
    vb_header.has_genotype_data       = vb->has_genotype_data;
    vb_header.num_samples             = BGEN32 (global_num_samples);
    vb_header.num_haplotypes_per_line = BGEN32 (vb->num_haplotypes_per_line);
    vb_header.num_samples_per_block   = BGEN32 (vb->num_samples_per_block);
    vb_header.num_sample_blocks       = BGEN32 (vb->num_sample_blocks);
    vb_header.ploidy                  = BGEN16 (vb->ploidy);
    vb_header.vb_data_size            = BGEN32 (vb->vb_data_size);
    vb_header.max_gt_line_len         = BGEN32 (vb->max_gt_line_len);
    vb_header.is_gtshark              = flag_gtshark;

    // create squeezed index - IF we have haplotype data AND more than one haplotype per line (i.e. my_squeeze_len > 0)
    if (my_squeeze_len) {
        buf_alloc (vb, &vb->haplotype_permutation_index_squeezed, my_squeeze_len, 1, "haplotype_permutation_index_squeezed", 0);

        uint16_t haplotype_index_checksum;
        squeeze (vb, (uint8_t *)vb->haplotype_permutation_index_squeezed.data, &haplotype_index_checksum,
                (unsigned *)vb->haplotype_permutation_index.data, 
                vb->num_haplotypes_per_line);

        vb_header.haplotype_index_checksum = BGEN16(haplotype_index_checksum);
    }
    
    // compress section into z_data - to be eventually written to disk by the I/O thread
    vb->z_data.name = "z_data"; // this is the first allocation of z_data - zfile_compress requires that it is pre-named
    vb->z_data.param = vb->variant_block_i;
    zfile_compress (vb, &vb->z_data, false, (SectionHeader*)&vb_header, 
                    my_squeeze_len ? vb->haplotype_permutation_index_squeezed.data : NULL);

    // account for haplotype_index as part of the haplotype_data for compression statistics
    // even though we store it in the variant_data header
    vb->z_section_bytes[SEC_HAPLOTYPE_DATA] += my_squeeze_len;
    vb->z_section_bytes[SEC_VB_HEADER]   -= my_squeeze_len;
}

// ZIP only: updating of the already compressed variant data section, after completion of all other sections
// note: this updates the z_data in memory (not on disk)
void zfile_update_compressed_vb_header (VariantBlock *vb, uint32_t vcf_first_line_i)
{
    uint32_t pos=0;

    SectionHeaderVbHeader *vb_header = (SectionHeaderVbHeader *)&vb->z_data.data[pos];
    vb_header->z_data_bytes = BGEN32 ((uint32_t)vb->z_data.len);
    vb_header->first_line   = BGEN32 (vcf_first_line_i);

    if (flag_show_vblocks) {
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u uncomprssed=%u compressed=%u num_sample_blocks=%u ht_per_line=%u max_gt_line_len=%u\n",
                 vb->variant_block_i, vcf_first_line_i, BGEN32 (vb_header->num_lines), 
                 BGEN32 (vb_header->vb_data_size), BGEN32 (vb_header->z_data_bytes), 
                 BGEN32 (vb_header->num_sample_blocks), BGEN32 (vb_header->num_haplotypes_per_line), 
                 BGEN32 (vb_header->max_gt_line_len));
    }
    // now we can finally encrypt the header - if needed
    if (vb_header->h.data_encrypted_len)  // non-zero if encrypted
        crypt_do (vb, (uint8_t*)vb_header, BGEN32 (vb_header->h.compressed_offset),
                  BGEN32 (vb_header->h.variant_block_i), vb_header->h.section_type, true);
}

void zfile_compress_dictionary_data (VariantBlock *vb, MtfContext *ctx, 
                                     uint32_t num_words, const char *data, uint32_t num_chars)
{
//printf ("Start compress dict vb_i=%u did_i=%u\n", vb->variant_block_i, ctx->did_i);
    START_TIMER;

    ASSERT (section_type_is_dictionary(ctx->dict_section_type),
            "Error: dict_section_type=%s is not a dictionary section", st_name(ctx->dict_section_type));

    SectionHeaderDictionary header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = ctx->dict_section_type;
    header.h.data_uncompressed_len = BGEN32 (num_chars);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderDictionary));
    header.h.variant_block_i       = BGEN32 (vb->variant_block_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.num_snips               = BGEN32 (num_words);
    header.dict_id                 = ctx->dict_id;

    if (flag_show_dict)
        fprintf (stderr, "%.*s (vb_i=%u, %s, did=%u, num_snips=%u):\t%.*s\n", DICT_ID_LEN, 
                dict_id_printable (ctx->dict_id).id, vb->variant_block_i, st_name(ctx->dict_section_type), 
                ctx->did_i, num_words, num_chars, data);

    if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num)
        fprintf (stderr, "%.*s\t", num_chars, data);

    z_file->dict_data.name = "z_file->dict_data"; // zfile_compress requires that it is set in advance
    zfile_compress (vb, &z_file->dict_data, true, (SectionHeader*)&header, data);

    COPY_TIMER (vb->profile.zfile_compress_dictionary_data)    
//printf ("End compress dict vb_i=%u did_i=%u\n", vb->variant_block_i, ctx->did_i);
}

void zfile_compress_b250_data (VariantBlock *vb, MtfContext *ctx)
{
    SectionHeaderBase250 header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.h.section_type          = ctx->b250_section_type;
    header.h.data_uncompressed_len = BGEN32 (ctx->b250.len);
    header.h.compressed_offset     = BGEN32 (sizeof(SectionHeaderBase250));
    header.h.variant_block_i       = BGEN32 (vb->variant_block_i);
    header.h.section_i             = BGEN16 (vb->z_next_header_i++);
    header.dict_id                 = ctx->dict_id;
    header.num_b250_items          = BGEN32 (ctx->mtf_i.len);
            
    zfile_compress (vb, &vb->z_data, false, (SectionHeader*)&header, ctx->b250.data);
}

void zfile_compress_haplotype_data_gtshark (VariantBlock *vb, const Buffer *haplotype_sections_data, unsigned sb_i)
{
    gtshark_compress_haplotype_data (vb, haplotype_sections_data, sb_i); // populates vb->gtshark_*
    
    vb->gtshark_exceptions_line_i.len *= sizeof (uint32_t);
    vb->gtshark_exceptions_ht_i.len *= sizeof (uint16_t);

    zfile_compress_section_data (vb, SEC_HT_GTSHARK_X_LINE, &vb->gtshark_exceptions_line_i);
    zfile_compress_section_data (vb, SEC_HT_GTSHARK_X_HTI, &vb->gtshark_exceptions_ht_i);
    zfile_compress_section_data (vb, SEC_HT_GTSHARK_X_ALLELE, &vb->gtshark_exceptions_allele);
    zfile_compress_section_data (vb, SEC_HT_GTSHARK_DB_DB, &vb->gtshark_db_db_data);
    zfile_compress_section_data (vb, SEC_HT_GTSHARK_DB_GT, &vb->gtshark_db_gt_data);

    // free buffers - they will be needed by the next section
    buf_free (&vb->gtshark_exceptions_line_i);
    buf_free (&vb->gtshark_exceptions_ht_i);
    buf_free (&vb->gtshark_exceptions_allele);
    buf_free (&vb->gtshark_db_db_data);
    buf_free (&vb->gtshark_db_gt_data);
}

void zfile_compress_section_data (VariantBlock *vb, SectionType section_type, Buffer *section_data)
{
    SectionHeader header;
    memset (&header, 0, sizeof(header)); // safety

    header.magic                 = BGEN32 (GENOZIP_MAGIC);
    header.section_type          = section_type;
    header.data_uncompressed_len = BGEN32 (section_data->len);
    header.compressed_offset     = BGEN32 (sizeof(header));
    header.variant_block_i       = BGEN32 (vb->variant_block_i);
    header.section_i             = BGEN16 (vb->z_next_header_i++);
    
    vb->z_data.name  = "z_data"; // zfile_compress requires that these are pre-set
    vb->z_data.param = vb->variant_block_i;
    zfile_compress (vb, &vb->z_data, false, (SectionHeader*)&header, section_data->data);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return true if data was read as requested, false if file has reached EOF and error otherwise
static void *zfile_read_from_disk (VariantBlock *vb, Buffer *buf, unsigned len, bool fail_quietly_if_not_enough_data)
{
    START_TIMER;

    ASSERT0 (len, "Error: in zfile_read_from_disk, len is 0");

    void *start = (void *)&buf->data[buf->len];

    unsigned memcpyied = 0;
    unsigned len_save  = len;

    while (len) {

        if (z_file->next_read == z_file->last_read) { // nothing left in read_buffer - replenish it from disk

            if (z_file->last_read != READ_BUFFER_SIZE && !memcpyied) return NULL; // EOF reached last time, nothing more to read

            if (z_file->last_read != READ_BUFFER_SIZE && memcpyied && fail_quietly_if_not_enough_data) return NULL; // partial read = quiet error as requested
            
            ASSERT (z_file->last_read == READ_BUFFER_SIZE, 
                    "Error: end-of-file while reading %s, read %u bytes, but wanted to read %u bytes", 
                    file_printname (z_file), memcpyied, len_save);

            z_file->last_read = fread (z_file->read_buffer, 1, READ_BUFFER_SIZE, (FILE *)z_file->file);
            z_file->next_read = 0;
            z_file->disk_so_far += z_file->last_read;

            ASSERT (z_file->last_read || !memcpyied, "Error: data requested could not be read bytes_so_far=%"PRIu64"", z_file->disk_so_far);
            if (!z_file->last_read) {
                COPY_TIMER (vb->profile.read);
                return NULL; // file is exhausted - nothing read
            }
        }

        unsigned memcpy_len = MIN (len, z_file->last_read - z_file->next_read);

        buf_add (buf, z_file->read_buffer + z_file->next_read, memcpy_len);
        len             -= memcpy_len;
        z_file->next_read += memcpy_len;

        memcpyied += memcpy_len;
    }

    COPY_TIMER (vb->profile.read);

    return start;
}

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file
int zfile_read_one_section (VariantBlock *vb, 
                            uint32_t original_vb_i, // the variant_block_i used for compressing. this is part of the encryption key. dictionaries are compressed by the compute thread/vb, but uncompressed by the I/O thread (vb=0)
                            Buffer *data, const char *buf_name, /* buffer to append */
                            unsigned header_size, SectionType expected_sec_type)
{
    // note: the first section is always read by zfile_read_one_section(). if it is a v1, or
    // if it is not readable (perhaps v1 encrypted) - we assume its a v1 VCF header
    // in v2+ the first section is always an unencrypted SectionHeaderGenozipHeader
    ASSERT0 (z_file->genozip_version != 1, "Error: zfile_read_one_section cannot read v1 data");

    bool is_encrypted = (expected_sec_type != SEC_GENOZIP_HEADER) &&
                         crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    unsigned header_offset = data->len;
    buf_alloc (vb, data, header_offset + header_size, 2, buf_name, 1);

    SectionHeader *header = zfile_read_from_disk (vb, data, header_size, false); // note: header in file can be shorter than header_size if its an earlier version

    // case: we're done! no more concatenated files
    if (!header && expected_sec_type == SEC_VCF_HEADER) return EOF; 

    ASSERT (header, "Error: Failed to read data from file %s while expecting section type %s: %s", 
            file_printname (z_file), st_name(expected_sec_type), strerror (errno));
    
    if (flag_show_headers) zfile_show_header (header, NULL);

    bool is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC;

    // decrypt header (note: except for SEC_GENOZIP_HEADER - this header is never encrypted)
    if (is_encrypted) {
        ASSERT (BGEN32 (header->magic) != GENOZIP_MAGIC, 
                "Error: password provided, but file  %s is not encrypted", file_printname (z_file));

        crypt_do (vb, (uint8_t*)header, header_size, original_vb_i, expected_sec_type, true); // negative section_i for a header
    
        is_magical = BGEN32 (header->magic) == GENOZIP_MAGIC; // update after decryption
    }

    ASSERT (is_magical, "Error: corrupt data (magic is wrong) when attempting to read section %s of variant_block_i=%u in file %s", 
            st_name (expected_sec_type), vb->variant_block_i, file_printname (z_file));

    unsigned compressed_offset   = BGEN32 (header->compressed_offset);
    ASSERT (compressed_offset, "Error: header.compressed_offset is 0 when reading section_type=%s", st_name(expected_sec_type));

    unsigned data_compressed_len = BGEN32 (header->data_compressed_len);
    unsigned data_encrypted_len  = BGEN32 (header->data_encrypted_len);

    unsigned data_len = MAX (data_compressed_len, data_encrypted_len);
    
    // check that we received the section type we expect, 
    ASSERT (expected_sec_type == header->section_type,
            "Error: Unexpected section type: expecting %s, found %s", st_name(expected_sec_type), st_name(header->section_type));

    ASSERT (compressed_offset == header_size,
            "Error: invalid header - expecting compressed_offset to be %u but found %u. section_type=%s", 
            header_size, compressed_offset, st_name(header->section_type));

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "zfile_read_one_section", 2);
    header = (SectionHeader *)&data->data[header_offset]; // update after realloc

    // read section data
    if (data_len) {
        ASSERT (zfile_read_from_disk (vb, data, data_len, false), 
                "Error: failed to read section data, section_type=%s: %s", st_name(header->section_type), strerror (errno));
    }

    return header_offset;
}

void zfile_read_all_dictionaries (uint32_t last_vb_i /* 0 means all VBs */)
{
    bool first=true;
    SectionListEntry *sl_ent;
    while (sections_get_next_dictionary (&sl_ent)) {

        if (last_vb_i && sl_ent->variant_block_i > last_vb_i) break;

        // upon encountering the first dictionary, move the cursor to the dictionaries, and reset the read buffers
        if (first) {
            file_seek (z_file, sl_ent->offset, SEEK_SET, false);
            first = false;
        }

        zfile_read_one_section (evb, sl_ent->variant_block_i,
                                &evb->z_data, "z_data", sizeof(SectionHeaderDictionary), sl_ent->section_type);    

        // update dictionaries in z_file->mtf_ctx with dictionary data 
        mtf_integrate_dictionary_fragment (evb, evb->z_data.data);

        buf_free (&evb->z_data);
    }

    if (flag_show_dict || dict_id_show_one_dict.num) 
        for (uint32_t did_i=0; did_i < z_file->num_dict_ids; did_i++) {
            MtfContext *ctx = &z_file->mtf_ctx[did_i];

            if (dict_id_printable (ctx->dict_id).num == dict_id_show_one_dict.num) 
                fprintf (stderr, "%.*s\t", ctx->dict.len, ctx->dict.data);
            
            if (flag_show_dict)
                fprintf (stderr, "%.*s (%s, did=%u, num_snips=%u):\t%.*s\n", DICT_ID_LEN, 
                        dict_id_printable (ctx->dict_id).id, st_name(ctx->dict_section_type), 
                                        did_i, ctx->word_list.len, ctx->dict.len, ctx->dict.data);
        }
}

void zfile_read_one_vb (VariantBlock *vb)
{ 
    START_TIMER;

    // The VB is read from disk here, in the I/O thread, and is decompressed in piz_uncompress_all_sections() in the 
    // Compute thread, with the exception of dictionaries that are handled here - this VBs dictionary fragments are
    // integrated into the global dictionaries.
    // Order of sections in a V2 VB:
    // 1. SEC_VB_HEADER - its data is the haplotype index
    // 2. SEC_INFO_SUBFIELD_B250 - Fields 1-9 b250 data
    // 3. SEC_INFO_SUBFIELD_B250 - All INFO subfield data
    // 4. All sample data: up 3 sections per sample block:
    //    4a. SEC_GENOTYPE_DATA - genotype data
    //    4b. SEC_PHASE_DATA - phase data
    //    4c. SEC_HAPLOTYPE_DATA or SEC_HAPLOTYPE_GTSHARK - haplotype data

    int vb_header_offset = zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data",
                                                   sizeof(SectionHeaderVbHeader), SEC_VB_HEADER);

    // note - use a macro and not a variable bc vb_header changes when z_data gets realloced as we read more data
    #define vb_header ((SectionHeaderVbHeader *)&vb->z_data.data[vb_header_offset])

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading variant_block_i=%u", vb->variant_block_i);

    // overlay all dictionaries (not just those that have fragments in this variant block) to the vb
    mtf_overlay_dictionaries_to_vb (vb);

    // read the the data sections (fields, info sub fields, genotype, phase, haplotype)

    buf_alloc (vb, &vb->z_section_headers, 1000 /* arbitrary initial value */ * sizeof(char*), 0, "z_section_headers", 1);
    
    ((unsigned *)vb->z_section_headers.data)[0] = vb_header_offset; // variant data header is at index 0

    unsigned section_i=1;

    // read the 8 fields (CHROM to FORMAT)    
    for (VcfFields f=CHROM; f <= FORMAT; f++) {
        ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
        zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeaderBase250), SEC_CHROM_B250 + f*2);   
    }

    // read the info subfield sections into memory (if any)
    vb->num_info_subfields = sections_count_info_b250s (vb->variant_block_i); // also used later in piz_uncompress_all_sections()
    for (unsigned sf_i=0; sf_i < vb->num_info_subfields; sf_i++) {
        ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
        zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeaderBase250), SEC_INFO_SUBFIELD_B250);    
    }

    if (flag_drop_genotypes) return; // if we have --drop-genotypes, no need to read the sample data
    
    // read the sample data
    uint32_t num_sample_blocks = BGEN32 (vb_header->num_sample_blocks);
    for (unsigned sb_i=0; sb_i < num_sample_blocks; sb_i++) {

        // make sure we have enough space for the section pointers
        buf_alloc (vb, &vb->z_section_headers, sizeof (unsigned) * (section_i + 3), 2, "z_section_headers", 2);

        if (vb_header->has_genotype_data) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_GENOTYPE_DATA);
        }

        if (vb_header->phase_type == PHASE_MIXED_PHASED) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_PHASE_DATA);
        }

        if (vb_header->num_haplotypes_per_line != 0 && !vb_header->is_gtshark) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HAPLOTYPE_DATA);    
        }

        if (vb_header->num_haplotypes_per_line != 0 && vb_header->is_gtshark) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HT_GTSHARK_X_LINE);    

            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HT_GTSHARK_X_HTI);    

            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HT_GTSHARK_X_ALLELE);    

            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HT_GTSHARK_DB_DB);    

            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, vb->variant_block_i, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HT_GTSHARK_DB_GT);    
        }
    }
    
    COPY_TIMER (vb->profile.zfile_read_one_vb);

    #undef vb_header
}

// returns the read data IF this is an invalid v2+ file, and hence might be v1
int16_t zfile_read_genozip_header (Md5Hash *digest) // out
{
    // read the footer from the end of the file
    file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, false);

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", file_printname (z_file));
    if (!ret) return EOF;
    
    // case: this is not a valid genozip v2+ file... maybe its v1?
    if (BGEN32 (footer.magic) != GENOZIP_MAGIC) {
        file_seek (z_file, 0, SEEK_SET, false);
        return MAYBE_V1;
    }

    // read genozip header
    uint64_t genozip_header_offset = BGEN64 (footer.genozip_header_offset);
    file_seek (z_file, genozip_header_offset, SEEK_SET, false);

    // note: for v1, we will use this function only for the very first VCF header (which will tell us this is v1)
    ret = zfile_read_one_section (evb, 0, &evb->z_data, "genozip_header", sizeof(SectionHeaderGenozipHeader), SEC_GENOZIP_HEADER);
    ASSERT0 (ret != EOF, "Error: unexpected EOF when reading genozip header");
    
    SectionHeaderGenozipHeader *header = (SectionHeaderGenozipHeader *)evb->z_data.data;

    int data_type = BGEN16 (header->data_type); 
    ASSERT (data_type == DATA_TYPE_VCF, "Error: unrecognized data_type=%d", data_type);

    ASSERT (header->genozip_version <= GENOZIP_FILE_FORMAT_VERSION, 
            "Error: %s cannot be openned because it was compressed with a newer version of genozip (version %u.x.x) while the version you're running is older (version %s). You might want to consider upgrading genozip to the newest version.",
            file_printname (z_file), header->genozip_version, GENOZIP_CODE_VERSION);

    ASSERT (header->encryption_type != ENCRYPTION_TYPE_NONE || !crypt_have_password(), 
            "Error: password provided, but file %s is not encrypted", file_printname (z_file));

    // get & test password, if file is encrypted
    if (header->encryption_type != ENCRYPTION_TYPE_NONE) {

        if (!crypt_have_password()) crypt_prompt_for_password();

        crypt_do (evb, header->password_test, sizeof(header->password_test), 0, SEC_EOF, true); // decrypt password test

        ASSERT (!memcmp (header->password_test, password_test_string, sizeof(header->password_test)),
                "Error: password is wrong for file %s", file_printname (z_file));
    }

    global_num_samples        = BGEN32 (header->num_samples); // possibly 0, if genozip header was not rewritten. in this case, piz will get it from the first VCF header, but genols will show 0
    z_file->genozip_version    = header->genozip_version;
    z_file->num_vcf_components = BGEN32 (header->num_vcf_components);
    *digest                   = header->md5_hash_concat; 

    zfile_uncompress_section (evb, header, &z_file->section_list_buf, "z_file->section_list_buf", SEC_GENOZIP_HEADER);
    z_file->section_list_buf.len /= sizeof (SectionListEntry); // fix len
    BGEN_sections_list();

    buf_free (&evb->z_data);

    return data_type;
}

void zfile_compress_genozip_header (uint16_t data_type, const Md5Hash *single_component_md5)
{
    SectionHeaderGenozipHeader header;

    // start with just the fields needed by sections_add_to_list
    memset (&header, 0, sizeof(SectionHeaderGenozipHeader)); // safety
    header.h.section_type = SEC_GENOZIP_HEADER;

    // "manually" add the genozip section to the section list - normally it is added in zfile_compress()
    // but in this case the genozip section containing the list will already be ready...
    sections_add_to_list (evb, &header.h);

    bool is_encrypted = crypt_have_password();

    uint32_t num_sections = z_file->section_list_buf.len;

    header.h.magic                      = BGEN32 (GENOZIP_MAGIC);
    header.h.compressed_offset          = BGEN32 (sizeof (SectionHeaderGenozipHeader));
    header.h.data_uncompressed_len      = BGEN32 (z_file->section_list_buf.len * sizeof (SectionListEntry));
    // gtshark is the only difference between version 2 and 3. therefore without gtshark we can leave it at 2.
    header.genozip_version              = GENOZIP_FILE_FORMAT_VERSION;
    header.data_type                    = BGEN16 (data_type);
    header.encryption_type              = is_encrypted ? ENCRYPTION_TYPE_AES256 : ENCRYPTION_TYPE_NONE;
    header.uncompressed_data_size       = BGEN64 (z_file->vcf_data_size_concat);
    header.num_samples                  = BGEN32 (global_num_samples);
    header.num_items_concat             = BGEN64 (z_file->num_lines_concat);
    header.num_sections                 = BGEN32 (num_sections); 
    header.num_vcf_components           = BGEN32 (z_file->num_vcf_components_so_far);

    if (flag_md5) {
        if (flag_concat) {
            md5_finalize (&z_file->md5_ctx_concat, &header.md5_hash_concat);
            if (flag_md5 && z_file->num_vcf_components_so_far > 1) fprintf (stderr, "Concatenated VCF MD5 = %s\n", md5_display (&header.md5_hash_concat, false));
        } 
        else 
            header.md5_hash_concat = *single_component_md5; // if not in concat mode - just copy the md5 of the single file
    }

    zfile_get_metadata (header.created);

    if (is_encrypted) {
        memcpy (header.password_test, password_test_string, sizeof(header.password_test));
        crypt_do (evb, header.password_test, sizeof(header.password_test), 0, SEC_EOF, true);
    }

    if (flag_show_gheader) sections_show_gheader (&header); // before BGEN_sections_list

    Buffer *z_data = &evb->z_data;

    uint64_t genozip_header_offset = z_file->disk_so_far + z_data->len; // capture before zfile_compress that increases len

    BGEN_sections_list();

    // compress section into z_data - to be eventually written to disk by the I/O thread
    z_file->section_list_buf.len *= sizeof (SectionListEntry); // change to counting bytes
    zfile_compress (evb, z_data, true, &header.h, z_file->section_list_buf.data);

    // add a footer to this section - this footer appears AFTER the genozip header data, 
    // facilitating reading the genozip header in reverse from the end of the file
    SectionFooterGenozipHeader footer;
    footer.magic                 = BGEN32 (GENOZIP_MAGIC);
    footer.genozip_header_offset = BGEN64 (genozip_header_offset);

    buf_alloc (evb, z_data, z_data->len + sizeof(SectionFooterGenozipHeader), 1.5, "z_data", 0);
    memcpy (&z_data->data[z_data->len], &footer, sizeof(SectionFooterGenozipHeader));
    z_data->len += sizeof(SectionFooterGenozipHeader);
    evb->z_section_bytes[SEC_GENOZIP_HEADER] += sizeof(SectionFooterGenozipHeader);
}

// reads the the genozip header section's header from a GENOZIP file - used by main_list, returns true if successful
bool zfile_get_genozip_header (uint64_t *uncompressed_data_size,
                               uint32_t *num_samples,
                               uint64_t *num_items_concat,
                               Md5Hash  *md5_hash_concat,
                               char *created, unsigned created_len /* caller allocates space */)
{
    // read the footer from the end of the file
    if (!file_seek (z_file, -sizeof(SectionFooterGenozipHeader), SEEK_END, true))
        return false;

    SectionFooterGenozipHeader footer;
    int ret = fread (&footer, sizeof (footer), 1, (FILE *)z_file->file);
    ASSERTW (ret == 1, "Skipping empty file %s", file_printname (z_file));    
    if (!ret) return false; // empty file / cannot read
    
    // case: this is not a valid genozip v2+ file... maybe its v1?
    if (BGEN32 (footer.magic) != GENOZIP_MAGIC) {
        file_seek (z_file, 0, SEEK_SET, false);
        return v1_vcf_header_get_vcf_header (uncompressed_data_size, num_samples, num_items_concat, 
                                             md5_hash_concat, created, created_len);
    }

    // read genozip header
    uint64_t genozip_header_offset = BGEN64 (footer.genozip_header_offset);
    if (!file_seek (z_file, genozip_header_offset, SEEK_SET, true))
        return false;

    SectionHeaderGenozipHeader header;
    int bytes = fread ((char*)&header, 1, sizeof(SectionHeaderGenozipHeader), (FILE *)z_file->file);
    if (bytes < sizeof(SectionHeaderGenozipHeader)) return false;

    ASSERTW (BGEN32 (header.h.magic) == GENOZIP_MAGIC, "Error reading %s: corrupt data", file_printname (z_file));
    if (BGEN32 (header.h.magic) != GENOZIP_MAGIC) return false;

    *uncompressed_data_size = BGEN64 (header.uncompressed_data_size);
    *num_samples            = BGEN32 (header.num_samples);
    *num_items_concat       = BGEN64 (header.num_items_concat);
    *md5_hash_concat        = header.md5_hash_concat;
    memcpy (created, header.created, MIN (FILE_METADATA_LEN, created_len));

    return true;
}

// updating the VCF bytes of a GENOZIP file. If we're compressing a simple VCF file, we will know
// the bytes upfront, but if we're concatenating or compressing a VCF.GZ, we will need to update it
// when we're done. num_lines can only be known after we're done with this VCF component.
// if we cannot update the header - that's fine, these fields are only used for the progress indicator on --list
bool zfile_update_vcf_header_section_header (off64_t pos_of_current_vcf_header, uint32_t max_lines_per_vb,
                                             Md5Hash *md5 /* out */)
{
    // rewind to beginning of current (latest) vcf header - nothing to do if we can't
    if (!file_seek (z_file, pos_of_current_vcf_header, SEEK_SET, true)) return false;

    unsigned len = crypt_padded_len (sizeof (SectionHeaderVCFHeader));

    // update the header of the single (current) vcf. 
    SectionHeaderVCFHeader *curr_header = &z_file->vcf_header_single;
    curr_header->vcf_data_size    = BGEN64 (z_file->vcf_data_size_single);
    curr_header->num_lines        = BGEN64 (z_file->num_lines_single);
    curr_header->max_lines_per_vb = BGEN32 (max_lines_per_vb);

    md5_finalize (&z_file->md5_ctx_single, &curr_header->md5_hash_single);
    *md5 = curr_header->md5_hash_single;
    if (flag_md5) fprintf (stderr, "MD5 = %s\n", md5_display (&curr_header->md5_hash_single, false));

    if (pos_of_current_vcf_header == 0) 
        z_file->vcf_header_first.md5_hash_single = curr_header->md5_hash_single; // first vcf - update the stored header 

    // encrypt if needed
    if (crypt_have_password()) 
        crypt_do (evb, (uint8_t *)curr_header, len, 0, curr_header->h.section_type, true);

    file_write (z_file, curr_header, len);
    fflush ((FILE*)z_file->file); // its not clear why, but without this fflush the bytes immediately after the first header get corrupted (at least on Windows with gcc)
    
    file_seek (z_file, 0, SEEK_END, false); // return to the end of the file

    return true; // success
}

#define V1_ZFILE // select the zfile functions of v1.c
#include "v1.c"
