// ------------------------------------------------------------------
//   txtfile.c
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
#include "txtfile.h"
#include "vblock.h"
#include "vblock.h"
#include "vblock.h"
#include "file.h"
#include "compressor.h"
#include "strings.h"
#include "zlib/zlib.h"

static void txtfile_update_md5 (const char *data, uint32_t len, bool is_2ndplus_txt_header)
{
    if (flag_md5) {
        if (flag_concat && !is_2ndplus_txt_header)
            md5_update (&z_file->md5_ctx_concat, data, len);
        
        md5_update (&z_file->md5_ctx_single, data, len);
    }
}

// peformms a single I/O read operation - returns number of bytes read 
static uint32_t txtfile_read_block (char *data, uint32_t max_bytes)
{
    START_TIMER;

    int32_t bytes_read=0;

    if (file_is_plain_or_ext_decompressor (txt_file)) {
        
        bytes_read = read (fileno((FILE *)txt_file->file), data, max_bytes); // -1 if error in libc
        ASSERT (bytes_read >= 0, "Error: read failed from %s: %s", txt_name, strerror(errno));

        // bytes_read=0 and we're using an external decompressor - it is either EOF or
        // there is an error. In any event, the decompressor is done and we can suck in its stderr to inspect it
        if (!bytes_read && file_is_read_via_ext_decompressor_type (txt_file->type)) {
            file_assert_ext_decompressor();
            goto finish; // all is good - just a normal end-of-file
        }
        
        txt_file->disk_so_far += (int64_t)bytes_read;

#ifdef _WIN32
        // in Windows using Powershell, the first 3 characters on an stdin pipe are BOM: 0xEF,0xBB,0xBF https://en.wikipedia.org/wiki/Byte_order_mark
        // these charactes are not in 7-bit ASCII, so highly unlikely to be present natrually in a VCF file
        if (txt_file->redirected && 
            txt_file->disk_so_far == (int64_t)bytes_read &&  // start of file
            bytes_read >= 3  && 
            (uint8_t)data[0] == 0xEF && 
            (uint8_t)data[1] == 0xBB && 
            (uint8_t)data[2] == 0xBF) {

            // Bomb the BOM
            bytes_read -= 3;
            memcpy (data, data + 3, bytes_read);
            txt_file->disk_so_far -= 3;
        }
#endif
    }
    else if (file_is_gz_type (txt_file->type)) {
        bytes_read = gzfread (data, 1, max_bytes, (gzFile)txt_file->file);
        
        if (bytes_read)
            txt_file->disk_so_far = gzconsumed64 ((gzFile)txt_file->file); 
    }
    else if (file_is_bz2_type (txt_file->type)) { 
        bytes_read = BZ2_bzread ((BZFILE *)txt_file->file, data, max_bytes);

        if (bytes_read)
            txt_file->disk_so_far = BZ2_consumed ((BZFILE *)txt_file->file); 
    } 
    else {
        ABORT ("txtfile_read_block: Invalid file type %s", ft_name (txt_file->type));
    }
    
finish:
    COPY_TIMER (evb->profile.read);

    return bytes_read;
}

// ZIP: returns the number of lines read 
void txtfile_read_header (bool is_first_txt, bool header_required,
                          char first_char)  // first character in every header line
{
    START_TIMER;

    int32_t bytes_read;
    char prev_char='\n';

    // read data from the file until either 1. EOF is reached 2. end of txt header is reached
    while (1) { 

        // enlarge if needed        
        if (!evb->txt_data.data || evb->txt_data.size - evb->txt_data.len < READ_BUFFER_SIZE) 
            buf_alloc (evb, &evb->txt_data, evb->txt_data.size + READ_BUFFER_SIZE, 1.2, "txt_data", 0);    

        bytes_read = txtfile_read_block (&evb->txt_data.data[evb->txt_data.len], READ_BUFFER_SIZE);

        if (!bytes_read) { // EOF
            ASSERT (!evb->txt_data.len || evb->txt_data.data[evb->txt_data.len-1] == '\n', 
                    "Error: invalid %s header in %s - expecting it to end with a newline", dt_name (txt_file->data_type), txt_name);
            goto finish;
        }

        const char *this_read = &evb->txt_data.data[evb->txt_data.len];

        ASSERT (!header_required || evb->txt_data.len || this_read[0] == first_char,
                "Error: %s is missing a %s header - expecting first character in file to be %c", 
                txt_name, dt_name (txt_file->data_type), first_char);

        // check stop condition - a line not beginning with a 'first_char'
        for (int i=0; i < bytes_read; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
            if (this_read[i] == '\n') 
                evb->num_lines++;   

            if (prev_char == '\n' && this_read[i] != first_char) {  

                uint32_t txt_header_len = evb->txt_data.len + i;
                evb->txt_data.len += bytes_read; // increase all the way - just for buf_copy

                // the excess data is for the next vb to read 
                buf_copy (evb, &txt_file->unconsumed_txt, &evb->txt_data, 1, txt_header_len,
                          bytes_read - i, "txt_file->unconsumed_txt", 0);

                txt_file->txt_data_so_far_single += i; 
                evb->txt_data.len = txt_header_len;

                goto finish;
            }
            prev_char = this_read[i];
        }

        evb->txt_data.len += bytes_read;
        txt_file->txt_data_so_far_single += bytes_read;
    }

finish:        
    // md5 header - with logic related to is_first
    txtfile_update_md5 (evb->txt_data.data, evb->txt_data.len, !is_first_txt);

    // md5 unconsumed_txt - always
    txtfile_update_md5 (txt_file->unconsumed_txt.data, txt_file->unconsumed_txt.len, false);

    COPY_TIMER (evb->profile.txtfile_read_header);
}

// ZIP
void txtfile_read_vblock (VBlock *vb) 
{
    START_TIMER;

    uint64_t pos_before = file_tell (txt_file);

    buf_alloc (vb, &vb->txt_data, global_max_memory_per_vb, 1, "txt_data", vb->vblock_i);    

    // start with using the unconsumed data from the previous VB (note: copy & free and not move! so we can reuse txt_data next vb)
    if (buf_is_allocated (&txt_file->unconsumed_txt)) {
        buf_copy (vb, &vb->txt_data, &txt_file->unconsumed_txt, 0 ,0 ,0, "txt_data", vb->vblock_i);
        buf_free (&txt_file->unconsumed_txt);
    }

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    while (vb->txt_data.len < global_max_memory_per_vb) {  // make sure there's at least READ_BUFFER_SIZE space available

        uint32_t bytes_one_read = txtfile_read_block (&vb->txt_data.data[vb->txt_data.len], 
                                                      MIN (READ_BUFFER_SIZE, global_max_memory_per_vb - vb->txt_data.len));

        if (!bytes_one_read) { // EOF - we're expecting to have consumed all lines when reaching EOF (this will happen if the last line ends with newline as expected)
            ASSERT (!vb->txt_data.len || vb->txt_data.data[vb->txt_data.len-1] == '\n', "Error: invalid input file %s - expecting it to end with a newline", txt_name);
            break;
        }

        // note: we md_udpate after every block, rather on the complete data (vb or txt header) when its done
        // because this way the OS read buffers / disk cache get pre-filled in parallel to our md5
        // Note: we md5 everything we read - even unconsumed data
        txtfile_update_md5 (&vb->txt_data.data[vb->txt_data.len], bytes_one_read, false);

        vb->txt_data.len += bytes_one_read;
    }

    // drop the final partial line which we will move to the next vb
    for (int32_t i=vb->txt_data.len-1; i >= 0; i--) {

        if (vb->txt_data.data[i] == '\n') {

            // in FASTQ - an "end of line" is one that the next character is @, or it is the end of the file
            if (txt_file->data_type == DT_FASTQ && i < vb->txt_data.len-1 && vb->txt_data.data[i+1] != '@') continue;

            // case: still have some unconsumed data, that we wish  to pass to the next vb
            uint32_t unconsumed_len = vb->txt_data.len-1 - i;
            if (unconsumed_len) {

                // the unconcusmed data is for the next vb to read 
                buf_copy (evb, &txt_file->unconsumed_txt, &vb->txt_data, 1, // evb, because dst buffer belongs to File
                          vb->txt_data.len - unconsumed_len, unconsumed_len, "txt_file->unconsumed_txt", vb->vblock_i);

                vb->txt_data.len -= unconsumed_len;
            }
            break;
        }
    }

    txt_file->txt_data_so_far_single += vb->txt_data.len;
    vb->vb_data_size = vb->txt_data.len; // initial value. it may change if --optimize is used.
    vb->vb_data_read_size = file_tell (txt_file) - pos_before; // plain or gz/bz2 compressed bytes read

    COPY_TIMER (vb->profile.txtfile_read_vblock);
}

// PIZ
unsigned txtfile_write_to_disk (const Buffer *buf)
{
    unsigned len = buf->len;
    char *next = buf->data;

    if (!flag_test) {
        while (len) {
            unsigned bytes_written = file_write (txt_file, next, len);
            len  -= bytes_written;
            next += bytes_written;
        }
    }

    if (flag_md5) md5_update (&txt_file->md5_ctx_concat, buf->data, buf->len);

    txt_file->txt_data_so_far_single += buf->len;
    txt_file->disk_so_far            += buf->len;

    return buf->len;
}

void txtfile_write_one_vblock (VBlockP vb_)
{
    START_TIMER;

    VBlockSAMP vb = (VBlockSAMP)vb_;
    
    txtfile_write_to_disk (&vb->txt_data);

    char s1[20], s2[20];
    ASSERTW (vb->txt_data.len == vb->vb_data_size || exe_type == EXE_GENOCAT, 
            "Warning: vblock %u (num_lines=%u) had %s bytes in the original %s file but %s bytes in the reconstructed file (diff=%d)", 
            vb->vblock_i, vb->num_lines, dt_name (txt_file->data_type),
            str_uint_commas (vb->vb_data_size, s1), str_uint_commas (vb->txt_data.len, s2), 
            (int32_t)vb->txt_data.len - (int32_t)vb->vb_data_size);

    COPY_TIMER (vb->profile.write);
}

// ZIP only - estimate the size of the txt data in this file. affects the hash table size and the progress indicator.
void txtfile_estimate_txt_data_size (VBlock *vb)
{
    uint64_t disk_size = txt_file->disk_size; 

    // case: we don't know the disk file size (because its stdin or a URL where the server doesn't provide the size)
    if (!disk_size) { 
        if (flag_stdin_size) disk_size = flag_stdin_size; // use the user-provided size, if there is one
        else return; // we're unable to estimate if the disk size is not known
    } 
    
    double ratio=1;

    bool is_no_gt_vcf = (txt_file->data_type == DT_VCF && !((VBlockVCFP)vb)->has_genotype_data);

    // if we decomprssed gz/bz2 data directly - we extrapolate from the observed compression ratio
    if (file_is_gz_type (txt_file->type) || file_is_bz2_type (txt_file->type))
        ratio = (double)vb->vb_data_size / (double)vb->vb_data_read_size;

    // for compressed files for which we don't have their size (eg streaming from an http server) - we use
    // estimates based on a benchmark compression ratio of files with and without genotype data
    else if (file_is_bcf_type (txt_file->type))
        // note: .bcf files might be compressed or uncompressed - we have no way of knowing as 
        // "bcftools view" always serves them to us in plain VCF format. These ratios are assuming
        // the bcf is compressed as it normally is.
        ratio = is_no_gt_vcf ? 55 : 8.5;

    else if (file_is_xz_type (txt_file->type)) 
        ratio = is_no_gt_vcf ? 171 : 12.7;

    else if (file_is_bam_type (txt_file->type))
        ratio = 4;

    else if (file_is_plain_type (txt_file->type))
        ratio = 1;

    else ABORT ("Error in file_estimate_txt_data_size: unspecified file_type=%u", txt_file->type);

    txt_file->txt_data_size_single = disk_size * ratio;
}
