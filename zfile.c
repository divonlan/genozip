// ------------------------------------------------------------------
//   zfile.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the file LICENSE.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <inttypes.h>
#include <bzlib.h>

#include "vczip.h"

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
            if (use_buf->size == items * size) break;
        }
    }
    ASSERT (use_buf, "Error: zfile_bzalloc could not find a free buffer", "");

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

    ASSERT (i < NUM_COMPRESS_BUFS, "Error: zfile_bzfree failed to find buffer to free", "");
}

static void zfile_compress(VariantBlock *vb, Buffer *z_data, SectionHeader *header, char *data, SectionType section_type) 
{ 
    unsigned compressed_offset = ENDN32 (header->compressed_offset);
    unsigned data_uncompressed_len = ENDN32 (header->data_uncompressed_len);

    // allocate what we think will be enough memory. usually this realloc does nothing, as the
    // memory we pre-allocate for z_data is sufficient
    buf_alloc (vb, z_data, z_data->len + compressed_offset + MAX (data_uncompressed_len / 10, 500), // usually compression is a lot better than 10X so this should be enough, except for small data
               1.5, "z_data", 0);

    // compress data
    unsigned data_compressed_len = z_data->size - z_data->len - compressed_offset; // actual memory available - usually more than we asked for in the realloc, because z_data is pre-allocated
    header->data_compressed_len = ENDN32 (data_compressed_len); 

    START_TIMER;

    bz_stream strm;
    strm.bzalloc = zfile_bzalloc;
    strm.bzfree  = zfile_bzfree;
    strm.opaque  = vb; // just passed to malloc/free

    int ret = BZ2_bzCompressInit (&strm, BZLIB_BLOCKSIZE100K, 0, 30);
    ASSERT (ret == BZ_OK, "Error: BZ2_bzCompressInit failed, ret=%d", ret);

    strm.next_in   = data;
    strm.next_out  = &z_data->data[z_data->len + compressed_offset];
    strm.avail_in  = data_uncompressed_len;
    strm.avail_out = data_compressed_len;

    ret = BZ2_bzCompress (&strm, BZ_FINISH);
    ASSERT (ret == BZ_FINISH_OK || ret == BZ_STREAM_END, "Error: BZ2_bzCompress failed, ret=%d", ret);

    BZ2_bzCompressEnd (&strm);

    // if output buffer is too small, increase it, and continue compressing where we left off
    if (ret == BZ_FINISH_OK) {

        buf_alloc (vb, z_data, z_data->len + compressed_offset + data_uncompressed_len + 50 /* > BZ_N_OVERSHOOT */, 1, "zfile_compress", 0);
        data_compressed_len = z_data->size - z_data->len - compressed_offset;

        strm.next_in   = data;
        strm.next_out  = &z_data->data[z_data->len + compressed_offset];
        strm.avail_in  = data_uncompressed_len;
        strm.avail_out = data_compressed_len;
        
        int ret = BZ2_bzCompressInit (&strm, BZLIB_BLOCKSIZE100K, 0, 30);
        ASSERT (ret == BZ_OK, "Error: BZ2_bzCompressInit failed, ret=%d", ret);

        ret = BZ2_bzCompress (&strm, BZ_FINISH);
        ASSERT (ret == BZ_STREAM_END,  "Error: second BZ2_bzCompress failed, ret=%d", ret);
    
        BZ2_bzCompressEnd (&strm);
    }
    if (vb) COPY_TIMER(vb->profile.compressor);

    data_compressed_len -= strm.avail_out;
    header->data_compressed_len = ENDN32 (data_compressed_len);   

    // copy header
    memcpy (&z_data->data[z_data->len], header, compressed_offset);

    z_data->len += compressed_offset + data_compressed_len;

    // do calculations for --showcontent option
    vb->z_section_bytes[section_type] += compressed_offset + data_compressed_len;
    vb->vcf_section_bytes[section_type] += data_uncompressed_len;
}

// uncompressed a block and adds a \0 at its end. Returns the length of the uncompressed block, without the \0.
void zfile_uncompress_section(VariantBlock *vb,
                              const void *section_header_p,
                              Buffer *uncompressed_data,
                              SectionType expected_section_type) 
{
    START_TIMER;
    const SectionHeader *section_header = section_header_p;
    
    unsigned compressed_offset     = ENDN32 (section_header->compressed_offset);
    unsigned data_compressed_len   = ENDN32 (section_header->data_compressed_len);
    unsigned data_uncompressed_len = ENDN32 (section_header->data_uncompressed_len);

    ASSERT (expected_section_type != SEC_VCF_HEADER ||
            ENDN32(((SectionHeaderVCFHeader *)section_header)->magic) == VCZIP_MAGIC,
            "Error: Input file is not a vcz file", "");

    ASSERT (section_header->section_type == expected_section_type, "Error: expecting section type %u but seeing %u", expected_section_type, section_header->section_type);

    buf_alloc (vb, uncompressed_data, data_uncompressed_len + 1, 1.1, "zfile_uncompress_section", 0); // +1 for \0
        
    uncompressed_data->len = data_uncompressed_len;
    
    bz_stream strm;
    strm.bzalloc = zfile_bzalloc;
    strm.bzfree  = zfile_bzfree;
    strm.opaque  = vb; // just passed to malloc/free
    {
        START_TIMER; // seperate block for another time

        int ret = BZ2_bzDecompressInit (&strm, 0, 0);
        ASSERT (ret == BZ_OK, "Error: BZ2_bzDecompressInit failed", "");

        strm.next_in   = (char*)section_header + compressed_offset;
        strm.next_out  = uncompressed_data->data;
        strm.avail_in  = data_compressed_len;
        strm.avail_out = uncompressed_data->len;

        ret = BZ2_bzDecompress (&strm);
        ASSERT (ret == BZ_STREAM_END || ret == BZ_OK, "Error: BZ2_bzDecompress failed, ret=%d, avail_in=%d, avail_out=%d", ret, strm.avail_in, strm.avail_out);

        BZ2_bzDecompressEnd (&strm);

        if (vb) COPY_TIMER(vb->profile.compressor);
    } 

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

    char *user = getenv ("USER")     ? getenv ("USER")     : "unknown";
    char *host = getenv ("HOSTNAME") ? getenv ("HOSTNAME") : "unknown";
    sprintf(metadata, "%.19s %.20s@%.30s", time_buf, user, host); // 71 chars + the string termintor \0 = 72
}

void zfile_write_vcf_header (VariantBlock *vb, Buffer *vcf_header_text)
{
    File *file = vb->z_file;
    
    SectionHeaderVCFHeader vcf_header;
    vcf_header.h.section_type          = SEC_VCF_HEADER;
    vcf_header.h.data_uncompressed_len = ENDN32 (vcf_header_text->len);
    vcf_header.h.compressed_offset     = ENDN32 (sizeof (SectionHeaderVCFHeader));
    vcf_header.magic                   = ENDN32 (VCZIP_MAGIC);
    vcf_header.vczip_version           = VCZIP_VERSION;
    vcf_header.compression_alg         = COMPRESSION_ALG_BZLIB;
    vcf_header.num_samples             = ENDN32 (global_num_samples);

    // placeholder values for these - we will update them later in zfile_update_vcf_header_section_header
    vcf_header.vcf_data_size           = ENDN64 (vb->vcf_file->vcf_data_size) /* 0 if gzipped */; 
    vcf_header.num_lines               = 0; // not yet known

    zfile_get_metadata (vcf_header.created);
    zfile_get_metadata (vcf_header.modified);

    static Buffer vcf_header_buf = EMPTY_BUFFER;

    buf_alloc ((VariantBlock*)vb, &vcf_header_buf, vcf_header_text->len / 3, // generous guess of compressed size
               1, "vcf_header_buf", 0); 

    zfile_compress ((VariantBlock*)vb, &vcf_header_buf, (SectionHeader*)&vcf_header, vcf_header_text->data, SEC_VCF_HEADER);

    fwrite (vcf_header_buf.data, 1, vcf_header_buf.len, file->file);

    file->disk_so_far      += vcf_header_buf.len;  // length of VCZ data writen to disk
    file->vcf_data_so_far  += vcf_header_text->len; // length of the original VCF header

    buf_free (&vcf_header_buf); 
}

void zfile_compress_variant_data (VariantBlock *vb)
{
    unsigned sizeof_header = sizeof (SectionHeaderVariantData) + squeeze_len (vb->num_haplotypes_per_line);

    buf_alloc (vb, &vb->vardata_header_buf, sizeof_header, 0, "zfile_compress_variant_data", vb->first_line);
    SectionHeaderVariantData *vardata_header = (SectionHeaderVariantData *)vb->vardata_header_buf.data;

    vardata_header->h.section_type          = SEC_VARIANT_DATA;
    vardata_header->h.data_uncompressed_len = ENDN32 (vb->variant_data_section_data.len);
    vardata_header->h.compressed_offset     = ENDN32 (sizeof_header);
    vardata_header->first_line              = ENDN32 (vb->first_line);
    vardata_header->num_lines               = ENDN16 (vb->num_lines);
    vardata_header->phase_type              = (char)vb->phase_type; 
    vardata_header->has_genotype_data       = vb->has_genotype_data;
    vardata_header->num_samples             = ENDN32 (global_num_samples);
    vardata_header->num_haplotypes_per_line = ENDN32 (vb->num_haplotypes_per_line);
    vardata_header->num_samples_per_block   = ENDN32 (vb->num_samples_per_block);
    vardata_header->num_sample_blocks       = ENDN32 (vb->num_sample_blocks);
    vardata_header->ploidy                  = ENDN16 (vb->ploidy);
    vardata_header->vcf_data_size           = ENDN32 (vb->vcf_data_size);

    unsigned short haplotype_index_checksum;
    squeeze (vb, vardata_header->haplotype_index, &haplotype_index_checksum,
             (unsigned *)vb->haplotype_permutation_index.data, 
             vb->num_haplotypes_per_line);

    vardata_header->haplotype_index_checksum = ENDN16(haplotype_index_checksum);

    zfile_compress (vb, &vb->z_data, (SectionHeader*)vardata_header, vb->variant_data_section_data.data, SEC_VARIANT_DATA);

    buf_free (&vb->vardata_header_buf);
}

void zfile_compress_section_data (VariantBlock *vb, unsigned section_type, Buffer *section_data)
{
    SectionHeader header;
    header.section_type          = section_type;
    header.data_uncompressed_len = ENDN32 (section_data->len);
    header.compressed_offset     = ENDN32 (sizeof(header));
    
    zfile_compress (vb, &vb->z_data, (SectionHeader*)&header, section_data->data, section_type);
}

// reads exactly the length required, error otherwise. manages read buffers to optimize I/O performance.
// this doesn't make a big difference for SSD, but makes a huge difference for HD
// return true if data was read as requested, false if file has reached EOF and error otherwise
static bool zfile_read_from_disk (VariantBlock *vb, char *out_str, unsigned len)
{
    ASSERT (len, "Error: in zfile_read_from_disk, len is 0", "");

    bool memcpyied = false;
    File *file = vb->z_file; // for code readability

    while (len) {

        if (file->next_read == file->last_read) { // nothing left in read_buffer - replenish it from disk

            if (file->last_read != READ_BUFFER_SIZE && !memcpyied) return false; // EOF reached last time, nothing more to read

            ASSERT (file->last_read == READ_BUFFER_SIZE, "Error: end-of-file of input file, read %"PRIu64" bytes", file->disk_so_far)
            file->last_read = fread (file->read_buffer, 1, READ_BUFFER_SIZE, file->file);
            file->next_read = 0;
            file->disk_so_far += file->last_read;

            ASSERT (file->last_read || !memcpyied, "Error: data requested could not be read bytes_so_far=%"PRIu64"", file->disk_so_far);
            if (!file->last_read) return false; // file is exhausted - nothing read
        }

        unsigned memcpy_len = MIN (len, file->last_read - file->next_read);

        if (out_str) memcpy (out_str, file->read_buffer + file->next_read, memcpy_len);

        len             -= memcpy_len;
        out_str         += memcpy_len;
        file->next_read += memcpy_len;

        memcpyied = true;
    }
    return true;
}

// read section header - called from the I/O thread, but for a specific VB
// returns offset of header within data, EOF if end of file
int zfile_read_one_section (VariantBlock *vb,
                            Buffer *data, /* buffer to append */
                            unsigned header_size, unsigned section_type,
                            bool allow_eof)
{
    unsigned header_offset = data->len;
    buf_alloc (vb, data, header_offset + header_size, 2, "zfile_read_one_section", 1);

    bool success = zfile_read_from_disk (vb, &data->data[header_offset], header_size);
    ASSERT (success || allow_eof, "Failed to read header", "");
    
    if (!success) return EOF; 

    data->len += header_size;

    SectionHeader *header = (SectionHeader *)&data->data[header_offset];
    unsigned compressed_offset   = ENDN32 (header->compressed_offset);
    unsigned data_compressed_len = ENDN32 (header->data_compressed_len);
    
    // if we expected a new VB and got a VCF header - that's ok - we allow concatenated VCZ files
    bool skip_this_vcf_header = (header->section_type == SEC_VCF_HEADER && section_type == SEC_VARIANT_DATA);

    // check that we received the section type we expect, 
    ASSERT (header->section_type == section_type || skip_this_vcf_header,  
            "Error: when reading file: expecting section type %u, but seeing %u", section_type, header->section_type);

    ASSERT (compressed_offset == header_size || section_type == SEC_VARIANT_DATA, // for variant data, we also have the permutation index
            "Error: invalid header - expecting compressed_offset to be %u but found %u", header_size, compressed_offset);

    // allocate more memory for the rest of the header + data (note: after this realloc, header pointer is no longer valid)
    buf_alloc (vb, data, header_offset + compressed_offset + data_compressed_len, 1, "zfile_read_one_section", 2);

    // read the rest of the header
    // if skip_this_vcf_header, VCF header is bigger than Variant Data header that was read
    if (section_type == SEC_VARIANT_DATA || skip_this_vcf_header) {

        int bytes_left_over = compressed_offset - header_size;
        ASSERT (bytes_left_over >= 0, "Error: expected bytes_left_over=%d to be >=0", bytes_left_over)

        if (bytes_left_over) { // there will be an Index only if this VCF has samples
            success = zfile_read_from_disk (vb, &data->data[data->len], bytes_left_over);
            ASSERT (success, "Failed to read variant header left over. bytes_left_over=%u", bytes_left_over);
        }
        data->len += bytes_left_over;
    }

    // read section data
    success = zfile_read_from_disk (vb, &data->data[data->len], data_compressed_len);
    ASSERT (success, "Error: failed to read section data, section_type=%u", section_type);

    data->len += data_compressed_len;

    // if we need to skip this section, return the next section to the user instead
    // note: even though skipped, we read the section data into z_data, or else our progress stats get messed up
    if (skip_this_vcf_header)
        return zfile_read_one_section (vb, data, header_size, section_type, allow_eof);

    return header_offset;
}

bool zfile_read_one_vb (VariantBlock *vb)
{ 
    START_TIMER;

    int vardata_header_offset = zfile_read_one_section (vb, &vb->z_data, 
                                                       sizeof(SectionHeaderVariantData), SEC_VARIANT_DATA, true);
    if (vardata_header_offset == EOF) {

        // update size - in case they were not known (pipe, gzip etc)
        vb->z_file->disk_size = vb->z_file->disk_so_far;
        vb->z_file->eof = true;

        COPY_TIMER (vb->profile.read);
        return false; // end of file
    }

    // note - copying values here z_data.data can get reallocated each call to zfile_read_one_section
    SectionHeaderVariantData *vardata_header = (SectionHeaderVariantData *)&vb->z_data.data[vardata_header_offset];
    unsigned num_sample_blocks       = ENDN32 (vardata_header->num_sample_blocks);
    bool has_genotype_data           = vardata_header->has_genotype_data;
    PhaseType phase_type             = vardata_header->phase_type;
    unsigned num_haplotypes_per_line = ENDN32 (vardata_header->num_haplotypes_per_line);

    vb->vcf_data_size = ENDN32 (vardata_header->vcf_data_size);

    buf_alloc (vb, &vb->z_section_headers, 1000 /* arbitrary initial value */ * sizeof(char*), 0, "z_section_headers", 1);
    
    ((unsigned *)vb->z_section_headers.data)[0] = vardata_header_offset; // variant data header is at index 0

    unsigned section_i=1;

    for (unsigned sb_i=0; sb_i < num_sample_blocks; sb_i++) {

        // make sure we have enough space for the section pointers
        buf_alloc (vb, &vb->z_section_headers, sizeof (unsigned) * (section_i + 3), 2, "z_section_headers", 2);
        unsigned *section_headers = (unsigned *)vb->z_section_headers.data; // update after each realloc

        if (has_genotype_data) {
            section_headers[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, &vb->z_data, sizeof(SectionHeader), SEC_GENOTYPE_DATA, false);
        }

        if (phase_type == PHASE_MIXED_PHASED) {
            section_headers[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, &vb->z_data, sizeof(SectionHeader), SEC_PHASE_DATA, false);
        }

        if (num_haplotypes_per_line) {
            section_headers[section_i++] = vb->z_data.len;
            zfile_read_one_section (vb, &vb->z_data, sizeof(SectionHeader), SEC_HAPLOTYPE_DATA, false);    
        }
    }

    COPY_TIMER (vb->profile.read);

    return true;
}

// updating the VCF bytes of a VCZ file. If we're compressing a simple VCF file, we will know
// the bytes upfront, but if we're concatenating or compressing a VCF.GZ, we will need to update it
// when we're done
void zfile_update_vcf_header_section_header (File *z_file, 
                                             long long vcf_data_size,
                                             long long vcf_num_lines)
{
    SectionHeaderVCFHeader vcf_header; // just for measuring the offset (must be static or it will be optimized away)

    ASSERT((char*)&vcf_header.num_lines - (char*)&vcf_header.vcf_data_size == sizeof (long long), "Error: looks like SectionHeaderVCFHeader changed", "");

    int ret = fseek (z_file->file, (char*)&vcf_header.vcf_data_size - (char*)&vcf_header, SEEK_SET);
    if (!ret) {
        long long header_vcf_data_size = ENDN64(vcf_data_size);
        fwrite (&header_vcf_data_size, sizeof (long long), 1, z_file->file); 

        long long header_num_lines = ENDN64(vcf_num_lines);
        fwrite (&header_num_lines, sizeof (long long), 1, z_file->file); 

    } else
        ASSERTW (false, "Warning: failed to update the vcf header section", "");
}

