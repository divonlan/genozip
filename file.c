// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "zlib/zlib.h"
#include "bzlib/bzlib.h"
#include "genozip.h"

File *file_open (const char *filename, FileMode mode, FileType expected_type)
{
    ASSERT0 (filename, "Error: filename is null");

    ASSERT (mode==WRITE || access(filename, F_OK)==0, "%s: cannot open file for reading: %s", global_cmd, filename);
    
    File *file = (File *)calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));

    file->name = filename; // usually comes from the command line - cannot be freed

    if (expected_type == VCF) {
        if (file_has_ext (file->name, ".vcf")) {
            file->type = VCF;
            file->file = fopen(file->name, mode == READ ? "r" : "wb"); // "wb" so Windows doesn't add ASCII 13
        }
        else if (file_has_ext (file->name, ".vcf.gz")) {
            file->type = VCF_GZ;
            file->file = gzopen (file->name, mode == READ ? "rb" : "wb");    
        }
        else if (file_has_ext (file->name, ".vcf.bz2")) {
            file->type = VCF_BZ2;
            file->file = BZ2_bzopen (file->name, mode == READ ? "rb" : "wb");    
        }
        else if (file_has_ext (file->name, ".bcf")) {
            ABORT ("To genozip BCF files, use in conjuction with bcftools, e.g.:\n"
                   "bcftools view -Ov %s | %s -o %.*s.vcf" GENOZIP_EXT, file->name, global_cmd, (int)(strlen(file->name)-4), file->name);
        }
        else if (file_has_ext (file->name, ".bcf.gz")) {
            ABORT ("To genozip BCF files, use in conjuction with bcftools, e.g.:\n"
                   "bcftools view -Ov %s | %s -o %.*s.vcf" GENOZIP_EXT, file->name, global_cmd, (int)(strlen(file->name)-7), file->name);
        }
        else {
            ABORT ("%s: file: %s - file must have a .vcf or .vcf.gz or .vcf.bz2 extension", global_cmd, file->name);
        }

        ASSERT (file->file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 
    }
    else { // GENOZIP
        if (file_has_ext (file->name, ".vcf" GENOZIP_EXT)) {
            file->type = GENOZIP;
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
            ABORT ("%s: file: %s - file must have a .vcf" GENOZIP_EXT " extension", global_cmd, file->name);
    }

    if (mode == READ) {
        // get file size
        struct stat64 st;
        int ret = stat64(file->name, &st);
        ASSERTW (!ret, "Warning: stat64(%s) failed: %s", file->name, strerror(errno));

        file->disk_size = ret ? 0 : st.st_size;

        if (file->type == VCF)
            file->vcf_data_size = file->disk_size; 

        // initialize read buffer indices
        file->last_read = file->next_read = READ_BUFFER_SIZE;
    }

    return file;
}

File *file_fdopen (int fd, FileMode mode, FileType type, bool initialize_mutex)
{
    File *file = (File *)calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));

    file->file = fdopen (fd, mode==READ ? "rb" : "wb");
    ASSERT (file->file, "%s: Failed to file descriptor %u: %s", global_cmd, fd, strerror (errno));

    file->type = type;
    file->last_read = file->next_read = READ_BUFFER_SIZE;

    if (initialize_mutex) {
        unsigned ret = pthread_mutex_init (&file->mutex, NULL);
        file->mutex_initialized = true;
        file->next_variant_i_to_merge = 1;
        ASSERT0 (!ret, "pthread_mutex_init failed");
    }
    
    return file;
}

void file_close (File **file_p)
{
    File *file = *file_p;
    *file_p = NULL;

    if (file->type == VCF_GZ) {
        int ret = gzclose_r((gzFile)file->file);
        ASSERTW (!ret, "Warning: failed to close vcf.gz file: %s", file->name ? file->name : "");
    }
    else if (file->type == VCF_BZ2) {
        BZ2_bzclose((BZFILE *)file->file);
    }
    else {
        int ret = fclose((FILE *)file->file);
        ASSERTW (!ret, "Warning: failed to close vcf file %s: %s", file->name ? file->name : "", strerror(errno));
    } 

    for (unsigned i=0; i < file->num_subfields; i++) 
        mtf_free_context (&file->mtf_ctx[i]);

    if (file->mutex_initialized) 
        pthread_mutex_destroy (&file->mutex);

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