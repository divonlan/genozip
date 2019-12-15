// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the file LICENSE.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <zlib.h>

#include "vczip.h"

File *file_open (const char *filename, FileMode mode, FileType expected_type)
{
    ASSERT (filename, "Error: filename is null", "");

    ASSERT (mode==WRITE || access(filename, F_OK)==0, "%s: cannot open %s", global_cmd, filename);
    
    File *file = calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));

    file->name = filename; // usually comes from the command line - cannot be freed
    unsigned fn_len = strlen ((char*)filename);

    if (expected_type == VCF) {
        if (file_has_ext (file->name, ".vcf")) {
            file->type = VCF;
            file->file = fopen(file->name, mode == READ ? "r" : "wb"); // "wb" so Windows doesn't add ASCII 13
            ASSERT (file->file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno));
        }
        else if (file_has_ext (file->name, ".vcf.gz")) {
            file->type = VCF_GZ;
            file->file = gzopen (file->name, mode == READ ? "rb" : "wb");    
            ASSERT (file, "%s: cannot open file %s", global_cmd, file->name);
        }
        else         
            ASSERT (false, "%s: file: %s - file must have a .vcf or .vcf.gz extension", global_cmd, file->name);
    }
    else { // VCZ
        if (file_has_ext (file->name, ".vcz")) {
            file->type = VCZ;
            file->file = fopen(file->name, mode == READ ? "rb" : "wb"); // "wb" so Windows doesn't add ASCII 13
            
            if (expected_type == VCZ_TEST && !file->file) {
                free (file);
                return NULL;
            }

            ASSERT (file->file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno));
        }
        else if (expected_type == VCZ_TEST)
            return NULL;
        
        else
            ASSERT (false, "%s: file: %s - file must have a .vcz extension", global_cmd, file->name);
    }

    if (mode == READ) {
        // get file size
        struct stat64 st;
        int ret = fstat64(fileno(file->file), &st);

        file->disk_size = ret ? 0 : st.st_size;

        if (file->type == VCF)
            file->vcf_data_size = file->disk_size; 

        // initialize read buffer indices
        file->last_read = file->next_read = READ_BUFFER_SIZE;
    }

    return file;
}

File *file_fdopen (int fd, FileMode mode, FileType type)
{
    File *file = calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));

    file->file = fdopen (fd, mode==READ ? "rb" : "wb");
    ASSERT (file->file, "%s: Failed to file descriptor %u: %s", global_cmd, fd, strerror (errno));

    file->type = type;
    file->last_read = file->next_read = READ_BUFFER_SIZE;

    return file;
}

void file_close (File **file_p)
{
    File *file = *file_p;
    *file_p = NULL;

    if (file->type == VCF_GZ) {
        int ret = gzclose_r(file->file);
        ASSERTW (!ret, "Warning: failed to close vcf.gz file: %s", file->name ? file->name : "");
    }
    else {
        int ret = fclose(file->file);
        ASSERTW (!ret, "Warning: failed to close vcf file %s: %s", file->name ? file->name : "", strerror(errno));
    } 

    // note: we don't free file->name, because it might come from getopt - and should not be freed

    free (file);
}

void file_remove (const char *filename)
{
    int ret = remove (filename); 
    ASSERTW (!ret, "Warning: failed to remove %s: %s", filename, strerror (errno));
}

bool file_has_ext (const char *filename, const char *extension)
{
    if (!filename) return false;

    unsigned ext_len = strlen (extension);
    unsigned fn_len  = strlen (filename);
    
    return fn_len > ext_len && !strncmp (&filename[fn_len-ext_len], extension, ext_len);
}