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

// globals
File *z_file   = NULL;
File *vcf_file = NULL;

File *file_open (const char *filename, FileMode mode, FileType expected_type)
{
    ASSERT0 (filename, "Error: filename is null");

    bool file_exists = (access (filename, F_OK) == 0);

    ASSERT (mode != READ  || file_exists, "%s: cannot open %s for reading: %s", global_cmd, filename, strerror(errno));
    ASSERT (mode != WRITE || !file_exists || flag_force || (expected_type==VCF && flag_test), 
            "%s: output file %s already exists: you may use --force to overwrite it", global_cmd, filename);

    File *file = (File *)calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));
    file->mode = mode;
    
    // copy filename 
    unsigned fn_size = strlen (filename) + 1; // inc. \0
    file->name = malloc (fn_size);
    memcpy (file->name, filename, fn_size);

    if (expected_type == VCF) {
        if (file_has_ext (file->name, ".vcf")) {
            file->type = VCF;

            // don't actually open the file if we're just testing in genounzip
            if (flag_test && mode == WRITE) return file;

            file->file = fopen(file->name, mode == READ ? "rb" : "wb"); // "rb"/"wb" so libc on Windows doesn't drop/add '\r' between our code and the disk. we will handle the '\r' explicitly.
        }
        else if (file_has_ext (file->name, ".vcf.gz")) {
            file->type = VCF_GZ;

            // don't actually open the file if we're just testing in genounzip
            if (flag_test && mode == WRITE) return file;

            file->file = gzopen64 (file->name, mode == READ ? "rb" : "wb");    
        }
        else if (file_has_ext (file->name, ".vcf.bz2")) {
            file->type = VCF_BZ2;

            // don't actually open the file if we're just testing in genounzip
            if (flag_test && mode == WRITE) return file;

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
                FREE (file);
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

        file->disk_size = file_get_size (file->name);

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
                 bool cleanup_memory) // optional - used to destroy buffers in the file is closed NOT near the end of the execution, eg when dealing with splitting concatenated files
{
    File *file = *file_p;
    *file_p = NULL;

    if (file->file) {

        if (file->type == VCF_GZ) {
            int ret = gzclose_r((gzFile)file->file);
            ASSERTW (!ret, "Warning: failed to close vcf.gz file: %s", file->name ? file->name : "");
        }
        else if (file->type == VCF_BZ2) {
            BZ2_bzclose((BZFILE *)file->file);
        }
        else {
            int ret = fclose((FILE *)file->file);

            if (ret && !errno) { // this is a telltale sign of a memory overflow
                buf_test_overflows(evb); // failing to close for no reason is a sign of memory issues
                // if its not a buffer - maybe its file->file itself
                fprintf (stderr, "Error: fclose() failed without an error, possible file->file pointer is corrupted\n");
            }

            ASSERTW (!ret, "Warning: failed to close file %s: %s", file->name ? file->name : "", strerror(errno)); // vcf or genozip
        } 
    }

    for (unsigned i=0; i < file->num_dict_ids; i++) 
        mtf_free_context (&file->mtf_ctx[i]);

    if (file->dicts_mutex_initialized) 
        pthread_mutex_destroy (&file->dicts_mutex);

    if (cleanup_memory) {
        buf_destroy (&file->dict_data);
        buf_destroy (&file->ra_buf);
        buf_destroy (&file->section_list_buf);
        buf_destroy (&file->section_list_dict_buf);
        buf_destroy (&file->v1_next_vcf_header);
        buf_destroy (&file->vcf_unconsumed_data);
    }

    if (file->name) FREE (file->name);
    
    FREE (file);
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
bool file_seek (File *file, int64_t offset, 
                int whence, // SEEK_SET, SEEK_CUR or SEEK_END
                bool soft_fail)
{
    // check if we can just move the read buffers rather than seeking
    if (file->mode == READ && file->next_read != file->last_read && whence == SEEK_SET) {
#ifdef __APPLE__
        int64_t move_by = offset - ftello ((FILE *)file->file);
#else
        int64_t move_by = offset - ftello64 ((FILE *)file->file);
#endif

        // case: move is within read buffer already in memory (ftello shows the disk location after read of the last buffer)
        // we just change the buffer pointers rather than discarding the buffer and re-reading
        if (move_by <= 0 && move_by >= -(int64_t)file->last_read) {
            file->next_read = file->last_read + move_by;
            return true;
        }
    }

#ifdef __APPLE__
    int ret = fseeko ((FILE *)file->file, offset, whence);
#else
    int ret = fseeko64 ((FILE *)file->file, offset, whence);
#endif

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

uint64_t file_tell (File *file)
{
#ifdef __APPLE__
    return ftello ((FILE *)file->file);
#else
    return ftello64 ((FILE *)file->file);
#endif
}

uint64_t file_get_size (const char *filename)
{
    struct stat64 st;
    
    int ret = stat64(filename, &st);
    ASSERT (!ret, "Error: failed accessing %s: %s", filename, strerror(errno));
    
    return st.st_size;
}

bool file_is_dir (const char *filename)
{
    struct stat64 st;
    
    int ret = stat64(filename, &st);
    ASSERT (!ret, "Error: failed accessing %s: %s", filename, strerror(errno));
    
    return S_ISDIR (st.st_mode);
}
