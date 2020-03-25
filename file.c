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
#include "stream.h"
#include "url.h"

// globals
File *z_file   = NULL;
File *vcf_file = NULL;

static StreamP input_decompressor  = NULL; // bcftools or xz - only one at a time
static StreamP output_compressor   = NULL; // bgzip

// global pointers - so the can be compared eg "if (mode == READ)"
const char *READ  = "rb";  // use binary mode (b) in read and write so Windows doesn't add \r
const char *WRITE = "wb";

char *file_exts[] = FILE_EXTS;

static FileType file_get_type (const char *filename)
{
    for (FileType ft=UNKNOWN_FILE_TYPE+1; ft < STDIN; ft++)
        if (file_has_ext (filename, file_exts[ft])) return ft;

    return UNKNOWN_FILE_TYPE;
}

File *file_open (const char *filename, FileMode mode, FileType expected_type)
{
    ASSERT0 (filename, "Error: filename is null");

    File *file = (File *)calloc (1, sizeof(File) + (mode == READ ? READ_BUFFER_SIZE : 0));

    bool is_remote = url_is_url (filename);
    bool file_exists;

    // is_remote is only possible in READ mode
    ASSERT (mode != WRITE || !is_remote, "%s: expecting output file %s to be local, not a URL", global_cmd, filename);

    int64_t url_file_size = 0; // will be -1 if the web/ftp site does not provide the file size
    const char *error = NULL;
    if (is_remote) {
        error = url_get_status (filename, &file_exists, &url_file_size); // accessing is expensive - get existance and size in one call
        if (url_file_size >= 0) file->disk_size = (uint64_t)url_file_size;
    }
    else {
        file_exists = (access (filename, F_OK) == 0);
        error = strerror (errno);
        if (file_exists && mode == READ) file->disk_size = file_get_size (filename);
    }

    // return null if genozip input file size is known to be 0, so we can skip it. note: file size of url might be unknown
    if (mode == READ && expected_type == VCF && file_exists && !file->disk_size && !url_file_size) {
        FREE (file);
        return NULL; 
    }

    ASSERT (mode != READ  || file_exists, "%s: cannot open %s for reading: %s", global_cmd, filename, error);
    ASSERT (mode != WRITE || !file_exists || flag_force || (expected_type==VCF && flag_test), 
            "%s: output file %s already exists: you may use --force to overwrite it", global_cmd, filename);

    // copy filename 
    unsigned fn_size = strlen (filename) + 1; // inc. \0
    file->name = malloc (fn_size);
    memcpy (file->name, filename, fn_size);

    file->mode = mode;
    file->type = file_get_type (file->name);

    // sanity check
    if (file->mode == READ  && command == ZIP)   ASSERT (file_is_vcf (file), "%s: input file must have one of the following extensions: " VCF_EXTENSIONS, global_cmd);
    if (file->mode == WRITE && command == ZIP)   ASSERT (file->type == VCF_GENOZIP, "%s: output file must have a " VCF_GENOZIP_ " extension", global_cmd);
    if (file->mode == READ  && command == UNZIP) ASSERT (file->type == VCF_GENOZIP, "%s: input file must have a " VCF_GENOZIP_ " extension", global_cmd); 
    if (file->mode == WRITE && command == UNZIP) ASSERT (file->type == VCF || file->type == VCF_GZ || file->type == VCF_BGZ, 
                                                         "%s: output file must have a .vcf or .vcf.gz or .vcf.bgz extension", global_cmd); 
    if (expected_type == VCF) {

        switch (file->type) {
        case VCF:
            // don't actually open the file if we're just testing in genounzip
            if (flag_test && mode == WRITE) return file;

            file->file = is_remote ? url_open (NULL, file->name) : fopen (file->name, mode);
            break;

        case VCF_GZ:
        case VCF_BGZ:
            if (mode == READ) {
                if (is_remote) { 
                    FILE *url_fp = url_open (NULL, file->name);
                    file->file = gzdopen (fileno(url_fp), mode); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
                }
                else
                    file->file = gzopen64 (file->name, mode); // for local files we decompress ourselves
            }
            else {
                char threads_str[20];
                sprintf (threads_str, "%u", global_max_threads);

                FILE *redirected_stdout_file = NULL;
                if (!flag_stdout) {
                    redirected_stdout_file = fopen (file->name, mode); // bgzip will redirect its output to this file
                    ASSERT (redirected_stdout_file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno));
                }
                output_compressor = stream_create (0, 0, 0, global_max_memory_per_vb, 
                                                   redirected_stdout_file, // output is redirected unless flag_stdout
                                                   0,
                                                   "bgzip", 
                                                   "--stdout", // either to the terminal or redirected to output file
                                                   "--threads", threads_str,
                                                   NULL);
                file->file = stream_to_stream_stdin (output_compressor);
            }
            break;

        case VCF_BZ2:
            if (is_remote) {            
                FILE *url_fp = url_open (NULL, file->name);
                file->file = BZ2_bzdopen (fileno(url_fp), mode); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
            }
            else
                file->file = BZ2_bzopen (file->name, mode);  // for local files we decompress ourselves   
            break;

        case VCF_XZ:
            stream_abort_if_cannot_run ("xz");

            // TO DO (xz and bcftools) - we suck the stderr into a pipe but we never show it - perhaps
            // in some cases we should?

            input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                               is_remote ? file->name : NULL,     // url
                                               "xz",                              // exec_name
                                               is_remote ? SKIP_ARG : file->name, // local file name 
                                               "--threads=8", "--decompress", "--keep", "--stdout", 
                                               flag_quiet ? "--quiet" : SKIP_ARG, 
                                               NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case BCF:
        case BCF_GZ:
        case BCF_BGZ:
            stream_abort_if_cannot_run ("bcftools");

            input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                               is_remote ? file->name : NULL,     // url                                        
                                               "bcftools",                        // exec_name 
                                               "view", "-Ov", "--threads", "8", 
                                               is_remote ? SKIP_ARG : file->name, // local file name 
                                               NULL);
            file->file = stream_from_stream_stdout (input_decompressor);
            break;
        
        default:
            ABORT ("%s: unrecognized file type: %s", global_cmd, file->name);
        }
    }
    
    else if (expected_type == VCF_GENOZIP) 
        file->file = is_remote ? url_open (NULL, file->name) : fopen (file->name, mode);
    
    else ABORT ("Error: invalid expected_type: %u", expected_type);

    ASSERT (file->file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    if (mode == READ) {

        if (file->type == VCF)
            file->vcf_data_size_single = file->disk_size; 

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

        if (file->mode == READ && (file->type == VCF_GZ || file->type == VCF_BGZ)) {
            int ret = gzclose_r((gzFile)file->file);
            ASSERTW (!ret, "Warning: failed to close vcf.gz file: %s", file_printname (file));
        }
        else if (file->mode == READ && file->type == VCF_BZ2) 
            BZ2_bzclose((BZFILE *)file->file);
        
        else if (file->mode == READ && 
                 (file->type == BCF || file->type == BCF_GZ || file->type == BCF_BGZ || file->type == VCF_XZ)) 
            stream_close (&input_decompressor, STREAM_WAIT_FOR_PROCESS);

        else if (file->mode == WRITE && (file->type == VCF_GZ || file->type == VCF_BGZ)) 
            stream_close (&output_compressor, STREAM_WAIT_FOR_PROCESS);

        else 
            FCLOSE (file->file, file_printname (file));
    }

    // free resources if we are NOT near the end of the execution. If we are at the end of the execution
    // it is faster to just let the process die
    if (cleanup_memory) {
            
        if (file->type == VCF_GENOZIP) { // reading or writing a .vcf.genozip (no need to worry about STDIN or STDOUT - they are by definition a single file - so cleaned up when process exits)
            for (unsigned i=0; i < file->num_dict_ids; i++)
                mtf_destroy_context (&file->mtf_ctx[i]);

            if (file->dicts_mutex_initialized) 
                pthread_mutex_destroy (&file->dicts_mutex);

            buf_destroy (&file->dict_data);
            buf_destroy (&file->ra_buf);
            buf_destroy (&file->section_list_buf);
            buf_destroy (&file->section_list_dict_buf);
            buf_destroy (&file->v1_next_vcf_header);
        }

        if (file_is_zip_read(file))
            buf_destroy (&file->vcf_unconsumed_data);

        if (file->name) FREE (file->name);
        
        FREE (file);
    }
}

size_t file_write (File *file, const void *data, unsigned len)
{
    size_t bytes_written = fwrite (data, 1, len, (FILE *)file->file);

    // if we're streaming our genounzip/genocat/genols output to another process and that process has 
    // ended prematurely then exit quietly. in genozip we display an error because this means the resulting
    // .genozip file will be corrupted
    if (!file->name && command != ZIP && errno == EINVAL) exit(0);

    ASSERT (bytes_written, "Error: failed to write %u bytes to %s: %s", 
            len, file->name ? file->name : "(stdout)", strerror(errno));
    return bytes_written;
}

void file_remove (const char *filename, bool fail_quietly)
{
    int ret = remove (filename); 
    ASSERTW (!ret || fail_quietly, "Warning: failed to remove %s: %s", filename, strerror (errno));
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
    if (command == ZIP && (file->type == VCF_GZ || file->type == VCF_BGZ))
        return gzconsumed64 ((gzFile)vcf_file->file); 
    
    if (command == ZIP && file->type == VCF_BZ2)
        return BZ2_consumed ((BZFILE *)vcf_file->file); 

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

void file_get_file (VariantBlockP vb, const char *filename, Buffer *buf, const char *buf_name, unsigned buf_param,
                    bool add_string_terminator)
{
    uint64_t size = file_get_size (filename);

    buf_alloc (vb, buf, size + add_string_terminator, 1, buf_name, buf_param);

    FILE *file = fopen (filename, "rb");
    ASSERT (file, "Error: cannot open %s: %s", filename, strerror (errno));

    size_t bytes_read = fread (buf->data, 1, size, file);
    ASSERT (bytes_read == (size_t)size, "Error reading file %s: %s", filename, strerror (errno));

    buf->len = size;

    if (add_string_terminator) buf->data[size] = 0;

    FCLOSE (file, filename);
}

// used when aborting due to an error. avoid the compressors outputting their own errors after our process is gone
void file_kill_external_compressors (void)
{
    stream_close (&input_decompressor, STREAM_KILL_PROCESS);
    stream_close (&output_compressor, STREAM_KILL_PROCESS);
}
