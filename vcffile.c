// ------------------------------------------------------------------
//   vcffile.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the file LICENSE.txt
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <inttypes.h>

#define Z_LARGE64
#include "zlib/zlib.h"

#include "vczip.h"

// we implement our own "getc" which manages read buffers a lot more efficiently
static inline char vcffile_get_char(VariantBlock *vb)
{
    // read buffer - critically impacts performace. 
    // for hard drives: If the buffer size is too big it drastically slows read time because it forces waits on the disk) 
    // for SSD too large buffers only slightly negative impact performance
    // for both SSD and HD if it is too small slights increases read time - because it consumes more CPU
    // tests on a PC show that 500K yields the near best performance for both HD and SSD (HD works slightly better with 100K)

    File *file = vb->vcf_file; // for code readability

    if (file->next_read == file->last_read) {
        
        if (file->type == VCF) {
            file->last_read = fread (file->read_buffer, 1, READ_BUFFER_SIZE, file->file);
            file->disk_so_far += (long long)file->last_read;
        }
        else { // VCF_GZ
            file->last_read   = gzfread (file->read_buffer, 1, READ_BUFFER_SIZE, file->file);
            file->disk_so_far = gzoffset64 (file->file); // for compressed files, we update by block read
        }

        file->next_read = 0;
    }

    char ret_char = file->last_read ? file->read_buffer[file->next_read++] : EOF;

    if (ret_char == EOF) 
        file->disk_size = file->disk_so_far; // in case it was not known

    return ret_char;
}

// get the next line in a text file terminated by \n or EOF. 
// Returns the line length (excluding the terminating \0) or 0 if no more lines exist.
// unlike fgets, the line is not limited by any particular length
bool vcffile_get_line(VariantBlock *vb, unsigned line_i_in_file /* 1-based */, Buffer *line)
{    
    File *file = vb->vcf_file; // for code readability

    double *avg_line_len_so_far = vb ? &file->avg_data_line_len : &file->avg_header_line_len;
    unsigned *lines_so_far      = vb ? &file->data_lines_so_far : &file->header_lines_so_far;

    // allocate only if buffer not allocated already - otherwise try to make do with what we have first
    unsigned buf_len = line->size ? buf_alloc (vb, line, line->size, 1, "vcffile_get_line", line_i_in_file) // just make buffer allocated, without mallocing new memory
                                  : buf_alloc (vb, line, MAX(1000, (unsigned)(*avg_line_len_so_far * 1.2)), 1, "vcffile_get_line", line_i_in_file); // we reuse the same buffer for every line
    unsigned str_len = 0;

    do {
        char c = vcffile_get_char (vb);
        
        if (c == EOF) {
            ASSERT(!str_len, "Invalid VCF file: Expecting file to end with a newline", "");

            file->eof = true;

            buf_free(line);
            return false;
        }
        line->data[str_len] = c;
        str_len++;
    
        if (str_len+1 > buf_len) // leave 1 for the \0 to be added 
            buf_len = buf_alloc (vb, line, MAX (str_len * 2, (unsigned)(*avg_line_len_so_far * 2)), 1, "vcffile_get_line", line_i_in_file);

    } while (line->data[str_len-1] != '\n');
            
    line->data[str_len] = '\0'; // terminate string

    *avg_line_len_so_far = (double)((*avg_line_len_so_far * *lines_so_far) + (str_len+1)) / (double)(*lines_so_far+1);
    (*lines_so_far)++;
    
    line->len = str_len;

    file->vcf_data_so_far += str_len;

    return true;
}

// this function manages a write buffer and writes to the disk only when the buffer is full (or end of data)
// this doesn't make a difference for SSD performance, but it makes a significant different for HD
unsigned vcffile_write_to_disk(File *vcf_file, const Buffer *buf)
{
    unsigned bytes_written = fwrite (buf->data, 1, buf->len, vcf_file->file);
    ASSERT (bytes_written==buf->len, "Error in writing VCF file: expected to write %u but wrote only %u", buf->len, bytes_written);

    vcf_file->vcf_data_so_far += bytes_written;
    vcf_file->disk_so_far     += bytes_written;

    return bytes_written;
}

void vcffile_write_one_variant_block (File *vcf_file, VariantBlock *vb)
{
    START_TIMER;

    unsigned size_written_this_vb = 0;

    for (unsigned line_i=0; line_i < vb->num_lines; line_i++) {
        DataLine *dl = &vb->data_lines[line_i];
        size_written_this_vb += vcffile_write_to_disk (vcf_file, &dl->line);
    }

    ASSERTW (size_written_this_vb == vb->vcf_data_size, 
            "Warning: Variant block %u (first_line=%u last_line=%u num_lines=%u) had %u bytes in the original VCF file but %u bytes in the reconstructed file", 
            vb->variant_block_i, vb->first_line, vb->first_line+vb->num_lines-1, vb->num_lines, vb->vcf_data_size, size_written_this_vb);

    COPY_TIMER (vb->profile.write);
}

void vcffile_compare_pipe_to_file (FILE *from_pipe, File *vcf_file)
{
    const unsigned buf_size = 500000;

    char *data_pipe = calloc (buf_size, 1);
    ASSERT (data_pipe, "Error: Failed to allocate data_pipe", "");

    char *data_file = calloc (buf_size, 1);
    ASSERT (data_file, "Error: Failed to allocate data_file", "");

    unsigned len_file, len_pipe;
    long long total_len=0;
    do {
        len_pipe = fread (data_pipe, 1, buf_size, from_pipe);
        
        if (vcf_file->type == VCF)
            len_file =   fread (data_file, 1, buf_size, vcf_file->file);
        else // VCF_GZ
            len_file = gzfread (data_file, 1, buf_size, vcf_file->file);

        const char *failed_text = "FAILED!!! Please contact vczip@blackpawventures.com to help fix this bug in vczip";

        unsigned min_len = MIN (len_file, len_pipe);
        bool failed = false;
        if (memcmp (data_pipe, data_file, min_len)) {

            for (unsigned i=0; i < min_len ; i++) {
                if (data_pipe[i] != data_file[i]) {
                    printf ("Data differs from pipe in character %"PRIu64" of the file. Showing the buffers starting at this character (100 char max):\n"
                            "***** After ZIP & PIZ *****\n%.*s\n***** ORIGINAL FILE *****\n%.*s\n", total_len + i, MIN (len_pipe, 100), &data_pipe[i], MIN (len_file, 100), &data_file[i]);                     
                    break;
                }
            }
            failed = true;
        } 

        if (len_pipe != len_file) {
            printf ("Length differs - and zip & piz length=%u ; original file length=%u\n", len_pipe, len_file);
            failed = true;
        }

        ASSERT (!failed, "%s", failed_text);

        total_len += len_pipe;
    } while (len_pipe);

    fprintf (stderr, "Success          \b\b\b\b\b\b\b\b\b\b\n");

    free (data_pipe);
    free (data_file);
}
