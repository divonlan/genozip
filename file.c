// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#define Z_LARGE64
#ifdef __APPLE__
    #define off64_t __int64_t
#endif
#include <zlib.h>
#include <bzlib.h>

#include "genozip.h"
#include "move_to_front.h"
#include "file.h"

File *file_open (const char *filename, FileMode mode, FileType expected_type)
{
    ASSERT0 (filename, "Error: filename is null");

    bool file_exists = (access (filename, F_OK) == 0);

    ASSERT (mode != READ  || file_exists, "%s: cannot open %s for reading: %s", global_cmd, filename, strerror(errno));
    ASSERT (mode != WRITE || !file_exists || flag_force, "%s: output file %s already exists: you may use --force to overwrite it", global_cmd, filename);

    File *file = (File *)calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));

    // copy filename 
    unsigned fn_size = strlen (filename) + 1; // inc. \0
    file->name = malloc (fn_size);
    memcpy (file->name, filename, fn_size);

    if (expected_type == VCF) {
        if (file_has_ext (file->name, ".vcf")) {
            file->type = VCF;
            file->file = fopen(file->name, mode == READ ? "r" : "wb"); // "wb" so Windows doesn't add ASCII 13
        }
        else if (file_has_ext (file->name, ".vcf.gz")) {
            file->type = VCF_GZ;
            file->file = gzopen64 (file->name, mode == READ ? "rb" : "wb");    
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
            
            if (expected_type == GENOZIP_TEST && !file->file) {
                free (file);
                return NULL;
            }

            ASSERT (file->file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno));
        }
        else if (expected_type == GENOZIP_TEST)
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
            file->vcf_data_size_single = file->vcf_data_size_concat = file->disk_size; 

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

    return file;
}

void file_close (File **file_p, 
                 VariantBlockP pseudo_vb) // optional - used to destroy buffers in the file is closed NOT near the end of the execution, eg when dealing with splitting concatenated files
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

    for (unsigned i=0; i < file->num_dict_ids; i++) 
        mtf_free_context (&file->mtf_ctx[i]);

    if (file->mutex_initialized) 
        pthread_mutex_destroy (&file->mutex);

    if (pseudo_vb) {
        if (file->dict_data.memory) buf_destroy (pseudo_vb, &file->dict_data);
        if (file->ra_buf.memory) buf_destroy (pseudo_vb, &file->ra_buf);
        if (file->section_list_buf.memory) buf_destroy (pseudo_vb, &file->section_list_buf);
        if (file->section_list_dict_buf.memory) buf_destroy (pseudo_vb, &file->section_list_dict_buf);
        if (file->v1_next_vcf_header.memory) buf_destroy (pseudo_vb, &file->v1_next_vcf_header);
    }

    if (file->name) free (file->name);
    
    // note: we don't free file->name, because it might come from getopt - and should not be freed
    free (file);
}

size_t file_write (File *file, const void *data, unsigned len)
{
    size_t bytes_written = fwrite (data, 1, len, (FILE *)file->file);
    ASSERT (bytes_written, "Error: failed to write %u bytes to %s: %s", 
            len, file->name ? file->name : "(stdout)", strerror(errno));

    return bytes_written;
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

// get basename of a filename - we write our own basename for Visual C and Windows compatability
const char *file_basename (const char *filename, bool remove_exe, const char *default_basename,
                           char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */)
{
    if (!filename) filename = default_basename;

    unsigned len = strlen (filename);
    if (remove_exe && file_has_ext (filename, ".exe")) len -= 4; // for Windows

    // get start of basename
    const char *start = filename;
    for (int i=len-1; i >= 0; i--)
        if (filename[i]=='/' || filename[i]=='\\') {
            start = &filename[i+1];
            break;
        }

    len = len - (start-filename);

    if (!basename) 
        basename = (char *)malloc (len + 1); // +1 for \0
    else
        len = MIN (len, basename_size-1);

    sprintf (basename, "%.*s", (int)len, start);

    return basename;
}

// returns true if successful. depending on soft_fail, a failure will either emit an error 
// (and exit) or a warning (and return).
bool file_seek (File *file, int64_t offset, int whence, bool soft_fail)
{
    int ret = fseeko64 ((FILE *)file->file, offset, whence);

    if (soft_fail) {
        if (!flag_stdout) {
            ASSERTW (!ret, errno == EINVAL ? "Error while reading file %s: it is too small%s" 
                                        : "Warning: fseeko failed on file %s: %s", 
                    file_printname (file),  errno == EINVAL ? "" : strerror (errno));
        }
    } 
    else {
        ASSERT (!ret, "Error: fseeko failed on file %s: %s", file_printname (file), strerror (errno));
    }

    // reset the read buffers
    if (!ret) file->next_read = file->last_read = READ_BUFFER_SIZE;

    return !ret;
}