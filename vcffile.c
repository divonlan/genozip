// ------------------------------------------------------------------
//   vcffile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
 
#include "profiler.h"

#ifdef __APPLE__
#define off64_t __int64_t // needed for for conda mac - otherwise zlib.h throws compilation errors
#endif
#define Z_LARGE64
#include <errno.h>
#include <bzlib.h>

#include "genozip.h"
#include "vcffile.h"
#include "vb.h"
#include "file.h"
#include "zlib/zlib.h"

static void vcffile_update_md5 (const char *data, uint32_t len, bool is_2ndplus_vcf_header)
{
    if (flag_md5) {
        if (flag_concat && !is_2ndplus_vcf_header)
            md5_update (&z_file->md5_ctx_concat, data, len);
        
        md5_update (&z_file->md5_ctx_single, data, len);
    }
}

// peformms a single I/O read operation - returns number of bytes read 
static uint32_t vcffile_read_block (char *data)
{
    START_TIMER;

    int32_t bytes_read=0;

    if (file_is_plain_vcf (vcf_file)) {
        
        bytes_read = read (fileno((FILE *)vcf_file->file), data, READ_BUFFER_SIZE);
        ASSERT (bytes_read >= 0, "Error: read failed from %s: %s", file_printname(vcf_file), strerror(errno));

        vcf_file->disk_so_far += (int64_t)bytes_read;

#ifdef _WIN32
        // in Windows using Powershell, the first 3 characters on an stdin pipe are BOM: 0xEF,0xBB,0xBF https://en.wikipedia.org/wiki/Byte_order_mark
        // these charactes are not in 7-bit ASCII, so highly unlikely to be present natrually in a VCF file
        if (vcf_file->type == STDIN && 
            vcf_file->disk_so_far == (int64_t)bytes_read &&  // start of file
            bytes_read >= 3  && 
            (uint8_t)data[0] == 0xEF && 
            (uint8_t)data[1] == 0xBB && 
            (uint8_t)data[2] == 0xBF) {

            // Bomb the BOM
            bytes_read -= 3;
            memcpy (data, data + 3, bytes_read);
            vcf_file->disk_so_far -= 3;
        }
#endif
    }
    else if (vcf_file->type == VCF_GZ || vcf_file->type == VCF_BGZ) { 
        bytes_read = gzfread (data, 1, READ_BUFFER_SIZE, (gzFile)vcf_file->file);
        
        if (bytes_read)
            vcf_file->disk_so_far = gzconsumed64 ((gzFile)vcf_file->file); 
    }
    else if (vcf_file->type == VCF_BZ2) { 
        bytes_read = BZ2_bzread ((BZFILE *)vcf_file->file, data, READ_BUFFER_SIZE);

        if (bytes_read)
            vcf_file->disk_so_far = BZ2_consumed ((BZFILE *)vcf_file->file); 
    } 
    else {
        ABORT0 ("Invalid file type");
    }
    
    COPY_TIMER (evb->profile.read);

    return bytes_read;
}


// ZIP: returns the number of lines read 
void vcffile_read_vcf_header (bool is_first_vcf) 
{
    START_TIMER;

    int32_t bytes_read;
    char prev_char=0;

    // read data from the file until either 1. EOF is reached 2. end of vcf header is reached
    while (1) { 

        // enlarge if needed        
        if (!evb->vcf_data.data || evb->vcf_data.size - evb->vcf_data.len < READ_BUFFER_SIZE) 
            buf_alloc (evb, &evb->vcf_data, evb->vcf_data.size + READ_BUFFER_SIZE, 1.2, "vcf_data", 0);    

        bytes_read = vcffile_read_block (&evb->vcf_data.data[evb->vcf_data.len]);

        if (!bytes_read) { // EOF
            ASSERT (!evb->vcf_data.len || evb->vcf_data.data[evb->vcf_data.len-1] == '\n', 
                    "Error: invalid VCF file %s while reading VCF header - expecting it to end with a newline", file_printname (vcf_file));
            goto finish;
        }

        const char *this_read = &evb->vcf_data.data[evb->vcf_data.len];

        ASSERT (evb->vcf_data.len || this_read[0] == '#',
                "Error: %s is missing a VCF header - expecting first character in file to be #", file_printname (vcf_file));

        // case VB header: check stop condition - a line beginning with a non-#
        for (int i=0; i < bytes_read; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
            if (this_read[i] == '\n') 
                evb->num_lines++;   

            if (prev_char == '\n' && this_read[i] != '#') {  

                uint32_t vcf_header_len = evb->vcf_data.len + i;
                evb->vcf_data.len += bytes_read; // increase all the way - just for buf_copy

                // the excess data is for the next vb to read 
                buf_copy (evb, &vcf_file->vcf_unconsumed_data, &evb->vcf_data, 1, vcf_header_len,
                          bytes_read - i, "vcf_file->vcf_unconsumed_data", 0);

                vcf_file->vcf_data_so_far_single += i; 
                evb->vcf_data.len = vcf_header_len;

                goto finish;
            }
            prev_char = this_read[i];
        }

        evb->vcf_data.len += bytes_read;
        vcf_file->vcf_data_so_far_single += bytes_read;
    }

finish:        
    // md5 header - with logic related to is_first_vcf
    vcffile_update_md5 (evb->vcf_data.data, evb->vcf_data.len, !is_first_vcf);

    // md5 vcf_unconsumed_data - always
    vcffile_update_md5 (vcf_file->vcf_unconsumed_data.data, vcf_file->vcf_unconsumed_data.len, false);

    COPY_TIMER (evb->profile.vcffile_read_vcf_header);
}

// ZIP
void vcffile_read_variant_block (VariantBlock *vb) 
{
    START_TIMER;

    uint64_t pos_before = file_tell (vcf_file);

    buf_alloc (vb, &vb->vcf_data, global_max_memory_per_vb, 1, "vcf_data", vb->variant_block_i);    

    // start with using the unconsumed data from the previous VB (note: copy & free and not move! so we can reuse vcf_data next vb)
    if (buf_is_allocated (&vcf_file->vcf_unconsumed_data)) {
        buf_copy (vb, &vb->vcf_data, &vcf_file->vcf_unconsumed_data, 0 ,0 ,0, "vcf_data", vb->variant_block_i);
        buf_free (&vcf_file->vcf_unconsumed_data);
    }

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    while (vb->vcf_data.len <= global_max_memory_per_vb - READ_BUFFER_SIZE) {  // make sure there's at least READ_BUFFER_SIZE space available

        uint32_t bytes_one_read = vcffile_read_block (&vb->vcf_data.data[vb->vcf_data.len]);

        if (!bytes_one_read) { // EOF - we're expecting to have consumed all lines when reaching EOF (this will happen if the last line ends with newline as expected)
            ASSERT (!vb->vcf_data.len || vb->vcf_data.data[vb->vcf_data.len-1] == '\n', "Error: invalid VCF file %s - expecting it to end with a newline", file_printname (vcf_file));
            break;
        }

        // note: we md_udpate after every block, rather on the complete data (vb or vcf header) when its done
        // because this way the OS read buffers / disk cache get pre-filled in parallel to our md5
        // Note: we md5 everything we read - even unconsumed data
        vcffile_update_md5 (&vb->vcf_data.data[vb->vcf_data.len], bytes_one_read, false);

        vb->vcf_data.len += bytes_one_read;
    }

    // drop the final partial line which we will move to the next vb
    for (int32_t i=vb->vcf_data.len-1; i >= 0; i--) {

        if (vb->vcf_data.data[i] == '\n') {
            // case: still have some unconsumed data, that we wish  to pass to the next vb
            uint32_t unconsumed_len = vb->vcf_data.len-1 - i;
            if (unconsumed_len) {

                // the unconcusmed data is for the next vb to read 
                buf_copy (evb, &vcf_file->vcf_unconsumed_data, &vb->vcf_data, 1, // evb, because dst buffer belongs to File
                          vb->vcf_data.len - unconsumed_len, unconsumed_len, "vcf_file->vcf_unconsumed_data", vb->variant_block_i);

                vb->vcf_data.len -= unconsumed_len;
            }
            break;
        }
    }

    vcf_file->vcf_data_so_far_single += vb->vcf_data.len;
    vb->vb_data_size = vb->vcf_data.len; // initial value. it may change if --optimize is used.
    vb->vb_data_read_size = file_tell (vcf_file) - pos_before; // plain or gz/bz2 compressed bytes read

    COPY_TIMER (vb->profile.vcffile_read_variant_block);
}

// PIZ
unsigned vcffile_write_to_disk (const Buffer *buf)
{
    unsigned len = buf->len;
    char *next = buf->data;

    if (!flag_test) {
        while (len) {
            unsigned bytes_written = file_write (vcf_file, next, len);
            len  -= bytes_written;
            next += bytes_written;
        }
    }

    if (flag_md5) md5_update (&vcf_file->md5_ctx_concat, buf->data, buf->len);

    vcf_file->vcf_data_so_far_single += buf->len;
    vcf_file->disk_so_far            += buf->len;

    return buf->len;
}

void vcffile_write_one_variant_block (VariantBlock *vb)
{
    START_TIMER;

    unsigned size_written_this_vb = 0;

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {
        Buffer *line = &vb->data_lines.piz[line_i].line;

        if (line->len) // if this line is not filtered out
            size_written_this_vb += vcffile_write_to_disk (line);
    }

    char s1[20], s2[20];
    ASSERTW (size_written_this_vb == vb->vb_data_size || exe_type == EXE_GENOCAT, 
            "Warning: Variant block %u (first_line=%u last_line=%u num_lines=%u) had %s bytes in the original VCF file but %s bytes in the reconstructed file (diff=%d)", 
            vb->variant_block_i, vb->first_line, vb->first_line+vb->num_lines-1, vb->num_lines, 
            buf_display_uint (vb->vb_data_size, s1), buf_display_uint (size_written_this_vb, s2), 
            (int32_t)size_written_this_vb - (int32_t)vb->vb_data_size);

    COPY_TIMER (vb->profile.write);
}

// ZIP only - estimate the size of the vcf data in this file. affects the hash table size and the progress indicator.
void vcffile_estimate_vcf_data_size (VariantBlock *vb)
{
    if (!vcf_file->disk_size) return; // we're unable to estimate if the disk size is not known
    
    double ratio=1;

    // if we decomprssed gz/bz2 data directly - we extrapolate from the observed compression ratio
    if (vcf_file->type == VCF_GZ || vcf_file->type == VCF_BGZ || vcf_file->type == VCF_BZ2) 
        ratio = (double)vb->vb_data_size / (double)vb->vb_data_read_size;

    // for compressed files for which we don't have their size (eg streaming from an http server) - we use
    // estimates based on a benchmark compression ratio of files with and without genotype data
    else if (vcf_file->type == BCF || vcf_file->type == BCF_GZ || vcf_file->type == BCF_BGZ)
        // note: .bcf files might be compressed or uncompressed - we have no way of knowing as 
        // "bcftools view" always serves them to us in plain VCF format. These ratios are assuming
        // the bcf is compressed as it normally is.
        ratio = vb->has_genotype_data ? 8.5 : 55;

    else if (vcf_file->type == VCF_XZ)
        ratio = vb->has_genotype_data ? 12.7 : 171;

    else if (vcf_file->type == VCF)
        ratio = 1;

    else ABORT ("Error in file_estimate_vcf_data_size: unspecified file_type=%u", vcf_file->type);

    vcf_file->vcf_data_size_single = vcf_file->disk_size * ratio;
}
