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
/*
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
bool vcffile_get_line(VariantBlock *vb, unsigned vcf_line_i, // 1-based  
                      bool skip_md5_vcf_header, Buffer *line, const char *buf_name)
{    
    START_TIMER;

    File *file = vb->vcf_file; // for code readability

    double *avg_line_len_so_far = vb ? &file->avg_data_line_len : &file->avg_header_line_len;
    unsigned *lines_so_far      = vb ? &file->data_lines_so_far : &file->header_lines_so_far;

    // note: we start small with only 100, to support single-individual files with global_max_lines_per_vb=128K  
    unsigned buf_len = line->size ? buf_alloc (vb, line, line->size, 1, buf_name, vcf_line_i) // just make buffer allocated, without mallocing new memory
                                  : buf_alloc (vb, line, MAX(100, (unsigned)(*avg_line_len_so_far * 1.2)), 1, buf_name, vcf_line_i); // we reuse the same buffer for every line
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
            buf_len = buf_alloc (vb, line, MAX (str_len * 2, (unsigned)(*avg_line_len_so_far * 2)), 1, buf_name, vcf_line_i);

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

    COPY_TIMER (vb->profile.vcffile_get_line);

    return true;
}
*/

/* cancel list
avg_header_line_len
read_buffer - change to z_read_buffer, 
 -> z_  last_read
*/

// peformms a single I/O read operation - returns number of bytes read 
static uint32_t vcf_read_block (File *file, char *data)
{
    START_TIMER;

    uint32_t bytes_read;

    if (file->type == VCF || file->type == STDIN) {
        
        bytes_read = fread (data, 1, READ_BUFFER_SIZE, (FILE *)file->file);
        file->disk_so_far += (int64_t)bytes_read;

        if (file->type == STDIN) {
#ifdef _WIN32
            // in Windows using Powershell, the first 3 characters on an stdin pipe are BOM: 0xEF,0xBB,0xBF https://en.wikipedia.org/wiki/Byte_order_mark
            if (file->disk_so_far == (int64_t)bytes_read &&  // start of file
                bytes_read >= 3  && 
                (uint8_t)file->read_buffer[0] == 0xEF && 
                (uint8_t)file->read_buffer[1] == 0xBB && 
                (uint8_t)file->read_buffer[2] == 0xBF) {

                // Bomb the BOM
                bytes_read -= 3;
                memcpy (data, data + 3, bytes_read);
                file->disk_so_far -= 3;
            }
#endif
            file->type = VCF; // we only accept VCF from stdin. TO DO: identify the file type by the magic number (first 2 bytes for gz and bz2)
        }
    }
    else if (file->type == VCF_GZ) { 
        bytes_read = gzfread (data, 1, READ_BUFFER_SIZE, (gzFile)file->file);
        
        if (bytes_read)
            file->disk_so_far = gzoffset64 ((gzFile)file->file); // for compressed files, we update by block read
    }
    else if (file->type == VCF_BZ2) { 
        bytes_read = BZ2_bzread ((BZFILE *)file->file, data, READ_BUFFER_SIZE);

        if (bytes_read)
            file->disk_so_far = BZ2_bzoffset ((BZFILE *)file->file); // for compressed files, we update by block read
    } 
    else {
        ABORT0 ("Invalid file type");
    }
    
    COPY_TIMER (evb->profile.read);

    return bytes_read;
}

void vcffile_read_variant_block (VariantBlock *vb) 
{
    START_TIMER;

    static uint32_t *line_starts = NULL;
    static uint32_t *line_lens   = NULL;

    if (!line_starts) {
        line_starts = malloc (global_max_lines_per_vb * sizeof (uint32_t));
        line_lens   = malloc (global_max_lines_per_vb * sizeof (uint32_t));
    }

    File *file = vb->vcf_file;
    uint32_t unconsumed_start=0, unconsumed_len=0; // index and len with vb->vcf_data
    uint32_t windows_newlines = 0;

    if (!vb->data_lines) vb->data_lines = calloc (global_max_lines_per_vb, sizeof (DataLine));

    uint32_t avg_vb_vcf_size_so_far = (vb->variant_block_i > 1) ? (uint32_t)(file->disk_so_far / (vb->variant_block_i-1)) : 0;

    // note: we start small with only 100, to support single-individual files with global_max_lines_per_vb=128K  
    unsigned requested_alloc_size =  MAX (file->vcf_unconsumed_data.len + READ_BUFFER_SIZE, 
                                     MAX (global_max_lines_per_vb * 100, 
                                          avg_vb_vcf_size_so_far));
    
    // correct for tiny files
    if (file->disk_size < (uint64_t)requested_alloc_size) requested_alloc_size = MAX ((uint32_t)file->disk_size, READ_BUFFER_SIZE);

    // read data from the file until either 1. EOF is reached 2. end of block is reached
    while (vb->num_lines < global_max_lines_per_vb) { 

        // enlarge if needed        
        { START_TIMER;
        if (vb->vcf_data.size - vb->vcf_data.len - unconsumed_len < READ_BUFFER_SIZE) 
            buf_alloc (vb, &vb->vcf_data, MAX (requested_alloc_size, vb->vcf_data.size + READ_BUFFER_SIZE), 1.2, "vcf_data", vb->variant_block_i);    
        COPY_TIMER (vb->profile.tmp1); }

        // go ahead and consume the excess data from the previous VB (note: copy & free and not move! so we can reuse vcf_data next vb)
        if (!vb->num_lines && buf_is_allocated (&file->vcf_unconsumed_data)) {
            START_TIMER;
            buf_copy (vb, &vb->vcf_data, &file->vcf_unconsumed_data, 0 ,0 ,0, "vcf_data", vb->variant_block_i);
            buf_free (&file->vcf_unconsumed_data);

            unconsumed_len = vb->vcf_data.len;
            vb->vcf_data.len = 0;
            COPY_TIMER (vb->profile.tmp2);
        }
        else {
            uint32_t bytes_one_read = vcf_read_block (file, &vb->vcf_data.data[vb->vcf_data.len + unconsumed_len]);

            if (!bytes_one_read) { // EOF - we're expecting to have consumed all lines when reaching EOF (this will happen if the last line ends with newline as expected)
                ASSERT (!unconsumed_len, "Error: invalid VCF file %s - expecting it to end with a newline", file_printname (file));
                break;
            }
            unconsumed_len += bytes_one_read;
        }

        char *unconsumed = &vb->vcf_data.data[unconsumed_start];

        // check stop condition - we reached global_max_lines_per_vb lines or EOF
{ START_TIMER;
        for (uint32_t i=0; i < unconsumed_len; i++) {

            // case: Windows-style /r/n line ending - we remove the \r
/*            if (i < unconsumed_len-1 && unconsumed[i] == '\r' && unconsumed[i+1] == '\n') { 
                unconsumed[i] = '\n';
                if (i < unconsumed_len - 2) memcpy (&unconsumed[i+1], &unconsumed[i+2], unconsumed_len - i - 2);
                unconsumed_len--;
                windows_newlines++;
            }
*/
            if (unconsumed[i] == '\n') {

                uint32_t unconsumed_line_len = i+1;

                line_starts[vb->num_lines] = unconsumed_start;
                line_lens  [vb->num_lines] = unconsumed_line_len;
                
                // conceptually, "consume" move bytesby moving them (accountingly) from unconsumped bytes to consumed
                vb->vcf_data.len      += unconsumed_line_len;
                unconsumed_start      += unconsumed_line_len;
                unconsumed_len        -= unconsumed_line_len;
                file->vcf_data_so_far += unconsumed_line_len + windows_newlines;

                vb->num_lines++;

                if (vb->num_lines == global_max_lines_per_vb) break;

                i=0; // update to loop to continue to next unconsumed line
                unconsumed = &vb->vcf_data.data[unconsumed_start];
            }
        }
COPY_TIMER (vb->profile.tmp3);}
    }

    // case: still have some unconsumed data, that we wish  to pass to the next vb
    if (unconsumed_len) {
START_TIMER;

        uint32_t consumed_len = vb->vcf_data.len;
        vb->vcf_data.len += unconsumed_len; // tmp incraese for buf_copy

        // the unconcusmed data is for the next vb to read 
        buf_copy (evb, &file->vcf_unconsumed_data, &vb->vcf_data, 1, // evb, because dst buffer belongs to File
                  consumed_len, unconsumed_len, "vcf_file->vcf_unconsumed_data", 0);

        vb->vcf_data.len = consumed_len; // restore
COPY_TIMER (vb->profile.tmp4);
    }

    // overlay all lines - only now, at the end, after vcf_data is done reallocing
{START_TIMER;
    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {
        DataLine *dl = &vb->data_lines[vb_line_i];
        buf_overlay (&dl->line, &vb->vcf_data, NULL, &line_starts[vb_line_i], "dl->line", dl->vcf_line_i);
        dl->line.len   = line_lens[vb_line_i];
        dl->vcf_line_i = vb->first_line + vb_line_i;  
//printf ("vb_line_i=%u: %.*s\n", vb_line_i, MIN (50,dl->line.len), dl->line.data); // DEBUG
    }
COPY_TIMER (vb->profile.tmp5);}

    file->disk_size = file->disk_so_far; // in case it was not known
    vb->vb_data_size = vb->vcf_data.len; // vb_data_size is redundant in ZIP at least, we can get rid of it one day

    COPY_TIMER (vb->profile.vcffile_read_variant_block);
}

// returns the number of lines read 
void vcffile_read_vcf_header (void) 
{
    START_TIMER;

    File *file = evb->vcf_file;
    uint32_t bytes_read;

    // read data from the file until either 1. EOF is reached 2. end of vcf header is reached
    while (1) { 

        // enlarge if needed        
        if (evb->vcf_data.size - evb->vcf_data.len < READ_BUFFER_SIZE) 
            buf_alloc (evb, &evb->vcf_data, evb->vcf_data.size + READ_BUFFER_SIZE, 1.2, "vcf_data", 0);    

        bytes_read = vcf_read_block (file, &evb->vcf_data.data[evb->vcf_data.len]);

        if (!bytes_read) { // EOF
            ASSERT (!evb->vcf_data.len || evb->vcf_data.data[evb->vcf_data.len-1] == '\n', 
                    "Error: invalid VCF file %s while reading VCF header - expecting it to end with a newline", file_printname (file));
            goto finish;
        }

        const char *this_read = &evb->vcf_data.data[evb->vcf_data.len];

        ASSERT (evb->vcf_data.len || this_read[0] == '#',
                "Error: %s is missing a VCF header - expecting first character in file to be #", file_printname (file));

        // case VB header: check stop condition - a line beginning with a non-#
        for (int i=0; i < bytes_read; i++) { // start from 1 back just in case it is a newline, and end 1 char before bc our test is 2 chars
            if (this_read[i] == '\n') 
                evb->num_lines++;   
                
            if ((i < bytes_read - 1 && this_read[i+1] != '#' && this_read[i] == '\n') ||   // vcf header ended if a line begins with a non-#
                (i==0 && evb->vcf_data.len>0 && this_read[i] != '#' && evb->vcf_data.data[evb->vcf_data.len-1]=='\n')) {

                uint32_t vcf_header_len = evb->vcf_data.len + i + 1;
                evb->vcf_data.len += bytes_read; // increase all the way - just for buf_copy

                // the excess data is for the next vb to read 
                buf_copy (evb, &file->vcf_unconsumed_data, &evb->vcf_data, 1, vcf_header_len,
                          bytes_read - (i+1), "vcf_file->vcf_unconsumed_data", 0);

                file->vcf_data_so_far += i+1; 
                evb->vcf_data.len = vcf_header_len;

                goto finish;
            }
        }

        evb->vcf_data.len += bytes_read;
        file->vcf_data_so_far += bytes_read;
    }

finish:        
    file->disk_size = file->disk_so_far; // in case it was not known

    COPY_TIMER (evb->profile.vcffile_read_vcf_header);
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
