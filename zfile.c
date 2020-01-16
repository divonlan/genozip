// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <bzlib.h>
#include "genozip.h"

#define BZLIB_BLOCKSIZE100K 9 /* maximum mem allocation for bzlib */

static Buffer compress_bufs[NUM_COMPRESS_BUFS] = {EMPTY_BUFFER, EMPTY_BUFFER, EMPTY_BUFFER, EMPTY_BUFFER}; // used by I/O thread when compressing out of context of a VB (i.e. VCF header)

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
    ASSERT (ret == BZ_OK, "Error: BZ2_bzCompressInit failed, ret=%d", ret);

    strm.next_in   = (char*)data;
    strm.next_out  = &z_data->data[z_data->len + compressed_offset];
    strm.avail_in  = data_uncompressed_len;
    strm.avail_out = *data_compressed_len;

    ret = BZ2_bzCompress (&strm, BZ_FINISH);
    ASSERT ((ok_is_ok && ret == BZ_FINISH_OK) || ret == BZ_STREAM_END, "Error: BZ2_bzCompress failed, ret=%d", ret);

    BZ2_bzCompressEnd (&strm);

    *data_compressed_len -= strm.avail_out;

    if (vb) COPY_TIMER(vb->profile.compressor);

    return ret;
}

static void zfile_compress (VariantBlock *vb, Buffer *z_data, SectionHeader *header, 
                            const char *data, SectionType section_type) 
{ 
    unsigned compressed_offset     = ENDN32 (header->compressed_offset);
    unsigned data_uncompressed_len = ENDN32 (header->data_uncompressed_len);
    
    unsigned header_padding=0;
    bool is_encrypted = crypt_get_encrypted_len (&compressed_offset, &header_padding); // set to 0 if no encryption
    
    const unsigned encryption_padding_reserve = crypt_max_padding_len(); // padding for the body

    // allocate what we think will be enough memory. usually this realloc does nothing, as the
    // memory we pre-allocate for z_data is sufficient
    buf_alloc (vb, z_data, z_data->len + compressed_offset + MAX (data_uncompressed_len / 10, 500) + encryption_padding_reserve, // usually compression is a lot better than 10X so this should be enough, except for small data
               1.5, "z_data", 0);

    // compress data
    unsigned data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve; // actual memory available - usually more than we asked for in the realloc, because z_data is pre-allocated

    int ret = zfile_compress_do (vb, z_data, data, compressed_offset, data_uncompressed_len, &data_compressed_len, true);

    // if output buffer is too small, increase it, and try again
    if (ret == BZ_FINISH_OK) {

        buf_alloc (vb, z_data, z_data->len + compressed_offset + data_uncompressed_len + 50 /* > BZ_N_OVERSHOOT */, 1, "zfile_compress", 0);
        data_compressed_len = z_data->size - z_data->len - compressed_offset - encryption_padding_reserve;

        zfile_compress_do (vb, z_data, data, compressed_offset, data_uncompressed_len, &data_compressed_len, false);
    }

    // get encryption related lengths
    unsigned data_encrypted_len=0, data_padding=0;
    if (is_encrypted) {
        data_encrypted_len = data_compressed_len;
        crypt_get_encrypted_len (&data_encrypted_len, &data_padding); // both are set to 0 if no encryption
    }

    // finalize & copy header
    header->compressed_offset   = ENDN32 (compressed_offset);  // updated value to include header padding
    header->data_compressed_len = ENDN32 (data_compressed_len);   
    header->data_encrypted_len  = ENDN32 (data_encrypted_len); 
    memcpy (&z_data->data[z_data->len], header, compressed_offset);

    // encrypt if needed - header & body separately
    unsigned total_z_len;
    if (is_encrypted) {
        // create good padding (just padding with 0 exposes a cryptographic volnurability)
        crypt_pad ((uint8_t*)&z_data->data[z_data->len], compressed_offset, header_padding);
        crypt_pad ((uint8_t*)&z_data->data[z_data->len + compressed_offset], data_encrypted_len, data_padding);

        // encrypt the header - except the variant data header that we encrypt later in zfile_update_compressed_variant_data_header()
        // we use section_i to generate a different AES key for each section
        uint32_t vb_i  = ENDN32 (header->variant_block_i);
        uint16_t sec_i = ENDN16 (header->section_i);

        if (section_type != SEC_VARIANT_DATA)
            crypt_do (vb, (uint8_t*)&z_data->data[z_data->len], compressed_offset, vb_i, -1-sec_i); // (use (-1-section_i) - different than header's +section_i)

        // encrypt the data body 
        crypt_do (vb, (uint8_t*)&z_data->data[z_data->len + compressed_offset], data_encrypted_len, vb_i, sec_i);
        
        total_z_len = compressed_offset + data_encrypted_len;
    }
    else
        total_z_len = compressed_offset + data_compressed_len;

    z_data->len += total_z_len;

    // do calculations for --showcontent option
    vb->z_section_bytes  [section_type] += total_z_len;
    vb->vcf_section_bytes[section_type] += data_uncompressed_len;
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
// when we get here, the header is already unencrypted zfile_read_one_section
void zfile_uncompress_section (VariantBlock *vb,
                               void *section_header_p,
                               Buffer *uncompressed_data,
                               SectionType expected_section_type) 
{
    START_TIMER;
    SectionHeader *section_header = (SectionHeader *)section_header_p;
    
    uint32_t compressed_offset     = ENDN32 (section_header->compressed_offset);
    uint32_t data_encrypted_len    = ENDN32 (section_header->data_encrypted_len);
    uint32_t data_compressed_len   = ENDN32 (section_header->data_compressed_len);
    uint32_t data_uncompressed_len = ENDN32 (section_header->data_uncompressed_len);
    uint32_t variant_block_i       = ENDN32 (section_header->variant_block_i);
    uint16_t section_i             = ENDN16 (section_header->section_i);

    // sanity checks
    ASSERT (section_header->section_type == expected_section_type, "Error: expecting section type %u but seeing %u", expected_section_type, section_header->section_type);
    ASSERT (variant_block_i == vb->variant_block_i, "Error: bad variant_block_i: in file=%u in vb=%u", variant_block_i, vb->variant_block_i);

    // decrypt data (in-place) if needed
    if (data_encrypted_len) 
        crypt_do (vb, (uint8_t*)section_header + compressed_offset, data_encrypted_len, variant_block_i, section_i);

    buf_alloc (vb, uncompressed_data, data_uncompressed_len + 1, 1.1, "zfile_uncompress_section", 0); // +1 for \0
        
    uncompressed_data->len = data_uncompressed_len;
    
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
    ASSERT (ret == BZ_STREAM_END || ret == BZ_OK, "Error: BZ2_bzDecompress failed, ret=%d, avail_in=%d, avail_out=%d", ret, strm.avail_in, strm.avail_out);

    BZ2_bzDecompressEnd (&strm);

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

void zfile_write_vcf_header (VariantBlock *vb, Buffer *vcf_header_text)
{
    File *file = vb->z_file;
    
    SectionHeaderVCFHeader vcf_header;
    memset (&vcf_header, 0, sizeof(vcf_header)); // safety

    vcf_header.h.magic                 = ENDN32 (GENOZIP_MAGIC);
    vcf_header.h.section_type          = SEC_VCF_HEADER;
    vcf_header.h.data_uncompressed_len = ENDN32 (vcf_header_text->len);
    vcf_header.h.compressed_offset     = ENDN32 (sizeof (SectionHeaderVCFHeader));
    vcf_header.genozip_version         = GENOZIP_VERSION;
    vcf_header.num_samples             = ENDN32 (global_num_samples);
    vcf_header.vcf_data_size           = ENDN64 (vb->vcf_file->vcf_data_size) /* 0 if gzipped - will be updated later*/; 
    vcf_header.num_lines               = NUM_LINES_UNKNOWN; 
    zfile_get_metadata (vcf_header.created);

    static Buffer vcf_header_buf = EMPTY_BUFFER;

    buf_alloc ((VariantBlock*)vb, &vcf_header_buf, vcf_header_text->len / 3, // generous guess of compressed size
               1, "vcf_header_buf", 0); 

    zfile_compress ((VariantBlock*)vb, &vcf_header_buf, (SectionHeader*)&vcf_header, vcf_header_text->data, SEC_VCF_HEADER);

    file_write (file, vcf_header_buf.data, vcf_header_buf.len);

    file->disk_so_far      += vcf_header_buf.len;   // length of GENOZIP data writen to disk
    file->vcf_data_so_far  += vcf_header_text->len; // length of the original VCF header

    buf_free (&vcf_header_buf); 

    // copy it to z_file - we might need to update it at the very end in zfile_update_vcf_header_section_header()
    memcpy (&vb->z_file->vcf_header, &vcf_header, sizeof (vcf_header));
}

void zfile_compress_variant_data (VariantBlock *vb)
{
    unsigned my_squeeze_len = squeeze_len (vb->num_haplotypes_per_line);
    unsigned sizeof_header = sizeof (SectionHeaderVariantData) + my_squeeze_len;

    buf_alloc (vb, &vb->vardata_header_buf, sizeof_header, 0, "zfile_compress_variant_data", vb->first_line);
    SectionHeaderVariantData *vardata_header = (SectionHeaderVariantData *)vb->vardata_header_buf.data;
    memset (vardata_header, 0, sizeof(SectionHeaderVariantData)); // safety
    
    vardata_header->h.magic                 = ENDN32 (GENOZIP_MAGIC);
    vardata_header->h.section_type          = SEC_VARIANT_DATA;
    vardata_header->h.data_uncompressed_len = ENDN32 (vb->variant_data_section_data.len);
    vardata_header->h.compressed_offset     = ENDN32 (sizeof_header);
    vardata_header->h.variant_block_i       = ENDN32 (vb->variant_block_i);
    vardata_header->h.section_i             = ENDN16 (vb->z_next_header_i++); // always 0
    vardata_header->first_line              = ENDN32 (vb->first_line);
    vardata_header->num_lines               = ENDN32 (vb->num_lines);
    vardata_header->phase_type              = (char)vb->phase_type; 
    vardata_header->has_genotype_data       = vb->has_genotype_data;
    vardata_header->num_samples             = ENDN32 (global_num_samples);
    vardata_header->num_haplotypes_per_line = ENDN32 (vb->num_haplotypes_per_line);
    vardata_header->num_samples_per_block   = ENDN32 (vb->num_samples_per_block);
    vardata_header->num_sample_blocks       = ENDN32 (vb->num_sample_blocks);
    vardata_header->ploidy                  = ENDN16 (vb->ploidy);
    vardata_header->vb_data_size            = ENDN32 (vb->vb_data_size);
    vardata_header->max_gt_line_len         = ENDN32 (vb->max_gt_line_len);
    vardata_header->num_subfields           = ENDN16 (vb->num_subfields);
    vardata_header->min_pos                 = ENDN64 (vb->min_pos);
    vardata_header->max_pos                 = ENDN64 (vb->max_pos);
    vardata_header->is_sorted_by_pos        = vb->is_sorted_by_pos;
    memcpy (vardata_header->chrom, vb->chrom, MAX_CHROM_LEN);

    uint16_t haplotype_index_checksum;
    squeeze (vb, vardata_header->haplotype_index, &haplotype_index_checksum,
             (unsigned *)vb->haplotype_permutation_index.data, 
             vb->num_haplotypes_per_line);

    vardata_header->haplotype_index_checksum = ENDN16(haplotype_index_checksum);

    zfile_compress (vb, &vb->z_data, (SectionHeader*)vardata_header, vb->variant_data_section_data.data, SEC_VARIANT_DATA);

    // account for haplotype_index as part of the haplotype_data for compression statistics
    // even though we store it in the variant_data header
    vb->z_section_bytes[SEC_HAPLOTYPE_DATA] += my_squeeze_len;
    vb->z_section_bytes[SEC_VARIANT_DATA]   -= my_squeeze_len;
    
    buf_free (&vb->vardata_header_buf);
}

// updating of the already compressed variant data section, after completion of all other sections
// note: this updates the z_data in memory (not on disk)
void zfile_update_compressed_variant_data_header (VariantBlock *vb,
                                                  unsigned pos, // index in vb->z_data
                                                  unsigned num_dictionary_sections)
{
    SectionHeaderVariantData *vardata_header = (SectionHeaderVariantData *)&vb->z_data.data[pos];
    vardata_header->num_dictionary_sections = ENDN16 ((uint16_t)num_dictionary_sections);
    vardata_header->z_data_bytes            = ENDN32 ((uint32_t)vb->z_data.len);

    // now we can finally encrypt the header - if needed
    if (vardata_header->h.data_encrypted_len)  // non-zero if encrypted
        crypt_do (vb, (uint8_t*)vardata_header, ENDN32 (vardata_header->h.compressed_offset),
                  ENDN32 (vardata_header->h.variant_block_i), -1-ENDN16 (vardata_header->h.section_i)); // negative section_i for a header
}

void zfile_compress_dictionary_data (VariantBlock *vb, SubfieldIdType subfield, 
                                     uint32_t num_words, const char *data, uint32_t num_chars)
{
    SectionHeaderDictionary header;
    memset (&header, 0, sizeof(header)); // safety

    header.h.magic                 = ENDN32 (GENOZIP_MAGIC);
    header.h.section_type          = SEC_DICTIONARY;
    header.h.data_uncompressed_len = ENDN32 (num_chars);
    header.h.compressed_offset     = ENDN32 (sizeof(SectionHeaderDictionary));
    header.h.variant_block_i       = ENDN32 (vb->variant_block_i);
    header.h.section_i             = ENDN16 (vb->z_next_header_i++);
    header.num_snips               = ENDN32 (num_words);
    memcpy (header.subfield_id, subfield.id, SUBFIELD_ID_LEN);

    zfile_compress (vb, &vb->z_data, (SectionHeader*)&header, data, SEC_DICTIONARY);
}

void zfile_compress_section_data (VariantBlock *vb, SectionType section_type, Buffer *section_data)
{
    SectionHeader header;
    memset (&header, 0, sizeof(header)); // safety

    header.magic                 = ENDN32 (GENOZIP_MAGIC);
    header.section_type          = section_type;
    header.data_uncompressed_len = ENDN32 (section_data->len);
    header.compressed_offset     = ENDN32 (sizeof(header));
    header.variant_block_i       = ENDN32 (vb->variant_block_i);
    header.section_i             = ENDN16 (vb->z_next_header_i++);
    
    zfile_compress (vb, &vb->z_data, (SectionHeader*)&header, section_data->data, section_type);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return true if data was read as requested, false if file has reached EOF and error otherwise
static void *zfile_read_from_disk (VariantBlock *vb, Buffer *buf, unsigned len)
{
    START_TIMER;

    ASSERT0 (len, "Error: in zfile_read_from_disk, len is 0");

    void *start = (void *)&buf->data[buf->len];

    bool memcpyied = false;
    File *file = vb->z_file; // for code readability

    while (len) {

        if (file->next_read == file->last_read) { // nothing left in read_buffer - replenish it from disk

            if (file->last_read != READ_BUFFER_SIZE && !memcpyied) return false; // EOF reached last time, nothing more to read

            ASSERT (file->last_read == READ_BUFFER_SIZE, "Error: end-of-file of input file, read %"PRIu64" bytes", file->disk_so_far)
            file->last_read = fread (file->read_buffer, 1, READ_BUFFER_SIZE, (FILE *)file->file);
            file->next_read = 0;
            file->disk_so_far += file->last_read;

            ASSERT (file->last_read || !memcpyied, "Error: data requested could not be read bytes_so_far=%"PRIu64"", file->disk_so_far);
            if (!file->last_read) {
                COPY_TIMER (vb->profile.read);
                return NULL; // file is exhausted - nothing read
            }
        }

        unsigned memcpy_len = MIN (len, file->last_read - file->next_read);

        buf_add (buf, file->read_buffer + file->next_read, memcpy_len);
        len             -= memcpy_len;
        file->next_read += memcpy_len;

        memcpyied = true;
    }

    COPY_TIMER (vb->profile.read);

    return start;
}

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file
int zfile_read_one_section (VariantBlock *vb,
                            Buffer *data, const char *buf_name, /* buffer to append */
                            unsigned header_size, SectionType section_type,
                            bool allow_eof)
{
    bool is_encrypted = crypt_get_encrypted_len (&header_size, NULL); // update header size if encrypted
    
    unsigned header_offset = data->len;
    buf_alloc (vb, data, header_offset + header_size, 2, buf_name, 1);

    SectionHeader *header = zfile_read_from_disk (vb, data, header_size);
    ASSERT0 (header || allow_eof, "Failed to read header");
    
    if (!header) return EOF; 

    // decrypt header
    if (is_encrypted)
        crypt_do (vb, (uint8_t*)header, header_size, vb->variant_block_i, --vb->z_next_header_i); // negative section_i for a header

    bool is_magical = ENDN32 (header->magic) == GENOZIP_MAGIC;
    if (!is_magical && !is_encrypted && section_type == SEC_VCF_HEADER) {

        // file appears to be encrypted but user hasn't provided a password    
        crypt_prompt_for_password();

        unsigned padding;
        is_encrypted = crypt_get_encrypted_len (&header_size, &padding); // update header size if encrypted
            
        if (padding) {
            char *header_extra_bytes = zfile_read_from_disk (vb, data, padding);
            ASSERT0 (header_extra_bytes, "Failed to read header padding");
        }

        crypt_do (vb, (uint8_t*)header, header_size, vb->variant_block_i, --vb->z_next_header_i);
        is_magical = ENDN32 (header->magic) == GENOZIP_MAGIC;
    }

    if (is_encrypted && section_type == SEC_VCF_HEADER) {
        ASSERT0 (is_magical, "Error: decryption failed - password is wrong"); // mostly likely its because of a wrong password
    } else {
        ASSERT (is_magical, "Error: corrupt data (magic is wrong) when reading file %s", vb->vcf_file->name ? vb->vcf_file->name : "(stdin)");
    }

    unsigned compressed_offset   = ENDN32 (header->compressed_offset);
    ASSERT (compressed_offset, "Error: header.compressed_offset is 0 when reading section_type=%u", section_type);

    unsigned data_compressed_len = ENDN32 (header->data_compressed_len);
    ASSERT (data_compressed_len, "Error: header.data_compressed_len is 0 when reading section_type=%u", section_type);

    unsigned data_encrypted_len  = ENDN32 (header->data_encrypted_len);

    unsigned data_len = MAX (data_compressed_len, data_encrypted_len);

    // if we expected a new VB and got a VCF header - that's ok - we allow concatenated GENOZIP files
    bool skip_this_vcf_header = (header->section_type == SEC_VCF_HEADER && section_type == SEC_VARIANT_DATA);

    ASSERT0 (section_type != SEC_VCF_HEADER || header->section_type == section_type || skip_this_vcf_header,  
             "Error: Input file is not a genozip file");

    // check that we received the section type we expect, 
    ASSERT (header->section_type == section_type || skip_this_vcf_header,  
            "Error: when reading file: expecting section type %u, but seeing %u", section_type, header->section_type);

    ASSERT (compressed_offset == crypt_padded_len (header_size) || section_type == SEC_VARIANT_DATA, // for variant data, we also have the permutation index
            "Error: invalid header - expecting compressed_offset to be %u but found %u", header_size, compressed_offset);

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_len, 2, "zfile_read_one_section", 2);

    // read the rest of the header (eg variant header has the haplotype index after the header)
    // if skip_this_vcf_header, VCF header is bigger than Variant Data header that was read
    if (section_type == SEC_VARIANT_DATA || skip_this_vcf_header) {

        int bytes_left_over = compressed_offset - header_size;
        ASSERT (bytes_left_over >= 0, "Error: expected bytes_left_over=%d to be >=0", bytes_left_over)

        if (bytes_left_over) { // there will be an Index only if this VCF has samples
            uint8_t *left_over_data = zfile_read_from_disk (vb, data, bytes_left_over);
            ASSERT (left_over_data, "Failed to read variant header left over. bytes_left_over=%u", bytes_left_over);

            if (is_encrypted) {
                ASSERT (bytes_left_over == crypt_padded_len(bytes_left_over), "Error: bad length of bytes_left_over=%u", bytes_left_over); // we expected it to be aligned, bc the total encryption block is aligned and header_size is aligned 
                
                // for the haplotype index - it is part of the header so we just continue the encryption stream
                if (section_type == SEC_VARIANT_DATA)
                    crypt_continue (vb, left_over_data, bytes_left_over);

                // for reading the VCF header's data - its a new encryption of the section's body
                else
                    crypt_do (vb, left_over_data, bytes_left_over, 0, 0);
            }
        }
    }

    // read section data. note: if this is a encrypted section, the remaining data might be up to 15
    // bytes longer than data_compressed_len, due to padding to encryption block boundary (16 bytes)
    ASSERT (zfile_read_from_disk (vb, data, data_len), 
            "Error: failed to read section data, section_type=%u", section_type);

    // if we need to skip this section, return the next section to the user instead
    // note: even though skipped, we read the section data into z_data, or else our progress stats get messed up
    if (skip_this_vcf_header)
        return zfile_read_one_section (vb, data, buf_name, header_size, section_type, allow_eof);

    return header_offset;
}

bool zfile_read_one_vb (VariantBlock *vb)
{ 
    START_TIMER;

    int vardata_header_offset = zfile_read_one_section (vb, &vb->z_data, "z_data",
                                                        sizeof(SectionHeaderVariantData), SEC_VARIANT_DATA, true);
    if (vardata_header_offset == EOF) {

        // update size - in case they were not known (pipe, gzip etc)
        if (!vb->z_file->disk_size) 
            vb->z_file->disk_size = vb->z_file->disk_so_far;
            
        vb->z_file->eof = true;

        COPY_TIMER (vb->profile.zfile_read_one_vb);
        return false; // end of file
    }

    // note - copying values here z_data.data can get reallocated each call to zfile_read_one_section
    SectionHeaderVariantData *vardata_header = (SectionHeaderVariantData *)&vb->z_data.data[vardata_header_offset];
    unsigned num_sample_blocks       = ENDN32 (vardata_header->num_sample_blocks);
    bool has_genotype_data           = vardata_header->has_genotype_data;
    PhaseType phase_type             = (PhaseType)vardata_header->phase_type;
    unsigned num_haplotypes_per_line = ENDN32 (vardata_header->num_haplotypes_per_line);
    unsigned num_dictionary_sections = ENDN16 (vardata_header->num_dictionary_sections);

    // dictionaries are processed right here by the dispatcher thread - the compute
    // thread only access the dictionaries on the z_file->mtf_ctx
    if (num_dictionary_sections) {

        unsigned start_dictionary_sections = vb->z_data.len;

        // read all sections into memory
        for (unsigned sf_i=0; sf_i < num_dictionary_sections; sf_i++) {

            unsigned start_i = vb->z_data.len; // vb->z_data.len is updated next, by zfile_read_one_section()
            zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeaderDictionary), SEC_DICTIONARY, false);    

            // update dictionaries in z_file->mtf_ctx with dictionary data from this VB
            mtf_integrate_dictionary_fragment (vb, &vb->z_data.data[start_i]);
        }

        vb->z_data.len = start_dictionary_sections; // shrink z_data back
    }

    // overlay all available dictionaries (not just those that have fragments in this variant block) to the vb
    mtf_overlay_dictionaries_to_vb (vb);

    // read the other sections

    buf_alloc (vb, &vb->z_section_headers, 1000 /* arbitrary initial value */ * sizeof(char*), 0, "z_section_headers", 1);
    
    ((unsigned *)vb->z_section_headers.data)[0] = vardata_header_offset; // variant data header is at index 0

    unsigned section_i=1;

    for (unsigned sb_i=0; sb_i < num_sample_blocks; sb_i++) {

        // make sure we have enough space for the section pointers
        buf_alloc (vb, &vb->z_section_headers, sizeof (unsigned) * (section_i + 3), 2, "z_section_headers", 2);

        if (has_genotype_data) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_GENOTYPE_DATA, false);
        }

        if (phase_type == PHASE_MIXED_PHASED) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_PHASE_DATA, false);
        }

        if (num_haplotypes_per_line) {
            ((unsigned *)vb->z_section_headers.data)[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, &vb->z_data, "z_data", sizeof(SectionHeader), SEC_HAPLOTYPE_DATA, false);    
        }
    }
    
    COPY_TIMER (vb->profile.zfile_read_one_vb);

    return true; 
}

// updating the VCF bytes of a GENOZIP file. If we're compressing a simple VCF file, we will know
// the bytes upfront, but if we're concatenating or compressing a VCF.GZ, we will need to update it
// when we're done. num_lines can only be known after we're done.
// if we cannot update the header - that's fine, these fields are only used for the progress indicator on --list
void zfile_update_vcf_header_section_header (VariantBlock *vb)
{
    SectionHeaderVCFHeader *header = &vb->z_file->vcf_header;
    unsigned len = crypt_padded_len (sizeof (SectionHeaderVCFHeader));
    
    // update the header contents
    header->vcf_data_size = ENDN64 (vb->z_file->vcf_concat_data_size);
    header->num_lines     = ENDN64 (vb->z_file->num_lines);
    if (flag_md5) md5_finalize (&vb->z_file->md5_ctx, &header->md5_hash);

    if (!flag_quiet) fprintf (stderr, "MD5 = %s\n", md5_display (&header->md5_hash, false));

    // encrypt if needed
    if (crypt_have_password()) crypt_do (vb, (uint8_t *)header, len, 0, -1); // 0,-1 are the VCF header's header

    // write back to disk
    if (fseek ((FILE *)vb->z_file->file, 0, SEEK_SET)) return; // rewind to the start
    file_write (vb->z_file, header, len);
}

