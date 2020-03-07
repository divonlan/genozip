// ------------------------------------------------------------------
//   vcffile.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
 
#include "profiler.h"

#ifdef __APPLE__
#define off64_t __int64_t // needed for for conda mac - otherwise zlib.h throws compilation errors
#endif
#define Z_LARGE64
#include <zlib.h>
#include <bzlib.h>

#include "genozip.h"
#include "vcffile.h"
#include "vb.h"
#include "file.h"

// we implement our own "getc" which manages read buffers a lot more efficiently
static inline char vcffile_get_char (VariantBlock *vb)
{
    // read buffer - critically impacts performace. 
    // for hard drives: If the buffer size is too big it drastically slows read time because it forces waits on the disk) 
    // for SSD too large buffers only slightly negative impact performance
    // for both SSD and HD if it is too small slights increases read time - because it consumes more CPU
    // tests on a PC show that 500K yields the near best performance for both HD and SSD (HD works slightly better with 100K)

    File *file = vb->vcf_file; // for code readability

    if (file->next_read == file->last_read) {
        
        START_TIMER;

        file->next_read = 0;

        if (file->type == VCF || file->type == STDIN) {
            file->last_read = fread (file->read_buffer, 1, READ_BUFFER_SIZE, (FILE *)file->file);
            file->disk_so_far += (long long)file->last_read;

            if (file->type == STDIN) {
#ifdef _WIN32
                // in Windows using Powershell, the first 3 characters on an stdin pipe are BOM: 0xEF,0xBB,0xBF https://en.wikipedia.org/wiki/Byte_order_mark
                if (file->last_read >=3 && 
                    (uint8_t)file->read_buffer[0] == 0xEF && 
                    (uint8_t)file->read_buffer[1] == 0xBB && 
                    (uint8_t)file->read_buffer[2] == 0xBF) {

                    file->next_read = 3; // skip the BOM;
                    file->disk_so_far -= 3;
                }
#endif
                file->type = VCF; // we only accept VCF from stdin. TO DO: identify the file type by the magic number (first 2 bytes for gz and bz2)
            }

        }
        else if (file->type == VCF_GZ) { 
            file->last_read = gzfread (file->read_buffer, 1, READ_BUFFER_SIZE, (gzFile)file->file);
            
            if (file->last_read)
                file->disk_so_far = gzoffset64 ((gzFile)file->file); // for compressed files, we update by block read
        }
        else if (file->type == VCF_BZ2) { 
            file->last_read = BZ2_bzread ((BZFILE *)file->file, file->read_buffer, READ_BUFFER_SIZE);

            if (file->last_read)
                file->disk_so_far = BZ2_bzoffset ((BZFILE *)file->file); // for compressed files, we update by block read
        } 
        else {
            ABORT0 ("Invalid file type");
        }

        COPY_TIMER (vb->profile.read);
    }

    char ret_char = file->last_read ? file->read_buffer[file->next_read++] : EOF;

    if (ret_char == EOF && !file->disk_size) 
        file->disk_size = file->disk_so_far; // in case it was not known

    return ret_char;
}

// get the next line in a text file terminated by \n or EOF. 
// Returns the line length (excluding the terminating \0) or 0 if no more lines exist.
// unlike fgets, the line is not limited by any particular length
bool vcffile_get_line(VariantBlock *vb, unsigned line_i_in_file /* 1-based */, bool skip_md5_vcf_header, Buffer *line, const char *buf_name)
{    
    File *file = vb->vcf_file; // for code readability

    double *avg_line_len_so_far = vb ? &file->avg_data_line_len : &file->avg_header_line_len;
    unsigned *lines_so_far      = vb ? &file->data_lines_so_far : &file->header_lines_so_far;

    // note: we start small with only 100, to support single-individual files with global_max_lines_per_vb=128K  
    unsigned buf_len = line->size ? buf_alloc (vb, line, line->size, 1, buf_name, line_i_in_file) // just make buffer allocated, without mallocing new memory
                                  : buf_alloc (vb, line, MAX(100, (unsigned)(*avg_line_len_so_far * 1.2)), 1, buf_name, line_i_in_file); // we reuse the same buffer for every line
    unsigned str_len = 0;

    do {
        char c = vcffile_get_char (vb);

        if (c == EOF) {
            ASSERT0 (!str_len, "Invalid VCF file: Expecting file to end with a newline");

            buf_free(line);
            return false;
        }
        line->data[str_len] = c;
        str_len++;
    
        if (str_len+1 > buf_len) // leave 1 for the \0 to be added 
            buf_len = buf_alloc (vb, line, MAX (str_len * 2, (unsigned)(*avg_line_len_so_far * 2)), 1, buf_name, line_i_in_file);

    } while (line->data[str_len-1] != '\n');
            
    // if this file somehow passed through Windows and had \n replace by \r\n - remove the \r
    bool windows_style_newline = str_len >= 2 && line->data[str_len-2] == '\r';
    
    if (windows_style_newline) {
        line->data[str_len-2] = '\n';
        str_len--;
    }

    line->data[str_len] = '\0'; // terminate string

    *avg_line_len_so_far = (double)((*avg_line_len_so_far * *lines_so_far) + (str_len+1)) / (double)(*lines_so_far+1);
    (*lines_so_far)++;
    
    line->len = str_len;

    file->vcf_data_so_far += str_len + windows_style_newline;

    if (flag_concat && (!skip_md5_vcf_header || line->data[0] != '#'))  // note that we ignore the directive to skip md5 for a concatenated header, if we discover this is actually the first line of the body
        {} //md5_update (&vb->z_file->md5_ctx_concat, line->data, line->len);

    //md5_update (&vb->z_file->md5_ctx_single, line->data, line->len);

    return true;
}

unsigned vcffile_write_to_disk(File *vcf_file, const Buffer *buf)
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

    md5_update (&vcf_file->md5_ctx_concat, buf->data, buf->len);

    vcf_file->vcf_data_so_far += buf->len;
    vcf_file->disk_so_far     += buf->len;

    return buf->len;
}

void vcffile_write_one_variant_block (File *vcf_file, VariantBlock *vb)
{
    START_TIMER;

    unsigned size_written_this_vb = 0;

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {
        Buffer *line = &vb->data_lines[line_i].line;

        if (line->len) // if this line is not filtered out
            size_written_this_vb += vcffile_write_to_disk (vcf_file, line);
    }

    ASSERTW (size_written_this_vb == vb->vb_data_size || exe_type == EXE_GENOCAT, 
            "Warning: Variant block %u (first_line=%u last_line=%u num_lines=%u) had %u bytes in the original VCF file but %u bytes in the reconstructed file", 
            vb->variant_block_i, vb->first_line, vb->first_line+vb->num_lines-1, vb->num_lines, vb->vb_data_size, size_written_this_vb);

    COPY_TIMER (vb->profile.write);
}
