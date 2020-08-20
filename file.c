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
#include <bzlib.h>
#include "zlib/zlib.h"
#include "genozip.h"
#include "context.h"
#include "file.h"
#include "stream.h"
#include "url.h"
#include "vblock.h"
#include "strings.h"
#include "compressor.h"
#include "mutex.h"

// globals
File *z_file   = NULL;
File *txt_file = NULL;

static StreamP input_decompressor  = NULL; // bcftools or xz, unzip or samtools - only one at a time
static StreamP output_compressor   = NULL; // bgzip, samtools, bcftools

static FileType stdin_type = UNKNOWN_FILE_TYPE; // set by the --input command line option

// global pointers - so the can be compared eg "if (mode == READ)"
const char *READ  = "rb";  // use binary mode (b) in read and write so Windows doesn't add \r
const char *WRITE = "wb";
const char *WRITEREAD = "wb+";

const char *file_exts[] = FILE_EXTS;

static const struct { FileType in; CompressionAlg comp_alg; FileType out; } txt_in_ft_by_dt[NUM_DATATYPES][30] = TXT_IN_FT_BY_DT;
static const FileType txt_out_ft_by_dt[NUM_DATATYPES][20] = TXT_OUT_FT_BY_DT;
static const FileType z_ft_by_dt[NUM_DATATYPES][20] = Z_FT_BY_DT;

// get data type by file type
DataType file_get_data_type (FileType ft, bool is_input)
{
    for (DataType dt=!flag_make_reference; dt < NUM_DATATYPES; dt++) 
        for (unsigned i=0; (is_input ? txt_in_ft_by_dt[dt][i].in : txt_out_ft_by_dt[dt][i]); i++)
            if ((is_input ? txt_in_ft_by_dt[dt][i].in : txt_out_ft_by_dt[dt][i]) == ft) return dt;

    return DT_NONE;
}

const char *file_plain_text_ext_of_dt (DataType dt) 
{
    return file_exts [txt_out_ft_by_dt[dt][0]];
}

// get genozip file type by txt file type
FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft)
{
    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].in == txt_ft) return txt_in_ft_by_dt[dt][i].out;

    return UNKNOWN_FILE_TYPE;
}

// get compression algorithm by txt file type
CompressionAlg file_get_comp_alg_by_txt_ft (DataType dt, FileType txt_ft, FileMode mode)
{
    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].in == txt_ft) return txt_in_ft_by_dt[dt][i].comp_alg;

    return (mode == WRITE ? COMP_NONE : COMP_UNKNOWN);
}

DataType file_get_dt_by_z_ft (FileType z_ft)
{
    for (DataType dt=0; dt < NUM_DATATYPES; dt++)
        for (unsigned i=0; z_ft_by_dt[dt][i]; i++)
            if (z_ft == z_ft_by_dt[dt][i]) 
                return dt;
    
    return DT_NONE;
}

FileType file_get_z_ft_by_dt (DataType dt)
{
    return z_ft_by_dt[dt][0];
}

// possible arguments for --input
static char *file_compressible_extensions(void)
{
    static char s[1000] = {0};

    for (DataType dt=1; dt < NUM_DATATYPES; dt++) { // start from 1, excluding DT_REFERENCE
        sprintf (&s[strlen (s)], "\n%s: ", dt_name (dt));

        for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
            sprintf (&s[strlen(s)], "%s ", &file_exts[txt_in_ft_by_dt[dt][i].in][1]);
    }

    return s;
}

FileType file_get_type (const char *filename, bool enforce_23andme_name_format)
{
    if (!filename) return UNKNOWN_FILE_TYPE;

    // 23andme files have the format "genome_Firstname_Lastname_optionalversion_timestamp.txt" or .zip
    if (enforce_23andme_name_format && file_has_ext (filename, ".txt")) 
        return (strstr (filename, "genome") && strstr (filename, "Full")) ? ME23 : UNKNOWN_FILE_TYPE;
    
    if (enforce_23andme_name_format && file_has_ext (filename, ".zip")) 
        return (strstr (filename, "genome") && strstr (filename, "Full")) ? ME23_ZIP : UNKNOWN_FILE_TYPE;
    
    for (FileType ft=UNKNOWN_FILE_TYPE+1; ft < AFTER_LAST_FILE_TYPE; ft++)
        if (file_has_ext (filename, file_exts[ft])) 
            return ft;

    return UNKNOWN_FILE_TYPE;
}

// returns the filename without the extension eg myfile.1.sam.gz -> myfile.1. 
// if raw_name is given, memory is allocated sufficiently to concatenate a extension. Otherwise, filename is overwritten
void file_get_raw_name_and_type (char *filename, char **raw_name, FileType *out_ft)
{
    unsigned len = strlen (filename);

    if (raw_name) {
        *raw_name = malloc (len + 30);
        memcpy (*raw_name, filename, len);
        (*raw_name)[len] = 0;
    }
    else 
        raw_name = &filename; // overwrite filename

    FileType ft = file_get_type (filename, true);
    if (ft != UNKNOWN_FILE_TYPE) 
        (*raw_name)[len - strlen (file_exts[ft])] = 0;

    if (out_ft) *out_ft = ft;
}

void file_set_input_size (const char *size_str)
{   
    unsigned len = strlen (size_str);

    for (unsigned i=0; i < len; i++) {
        ASSINP (IS_DIGIT (size_str[i]), "%s: expecting the file size in bytes to be a positive integer: %s", global_cmd, size_str);

        flag_stdin_size = flag_stdin_size * 10 + (size_str[i] - '0'); 
    }
}

void file_set_input_type (const char *type_str)
{
    char ext[strlen (type_str) + 2]; // +1 for . +1 for \0
    sprintf (ext, ".%s", type_str);

    str_to_lowercase (ext); // lower-case to allow case-insensitive --input argument (eg vcf or VCF)

    stdin_type = file_get_type (ext, false); // we don't enforce 23andMe name format - any .txt or .zip will be considered ME23

    if (file_get_data_type (stdin_type, true) != DT_NONE) return; // all good 

    // the user's argument is not an accepted input file type - print error message
    ABORT ("%s: --input (or -i) must be ones of these: %s", global_cmd, file_compressible_extensions());
}
    
FileType file_get_stdin_type (void)
{
    return stdin_type;
}

static void file_ask_user_to_confirm_overwrite (const char *filename)
{
    fprintf (stderr, "%s: output file %s already exists: in the future, you may use --force to overwrite\n", global_cmd, filename);
    
    if (!isatty(0) || !isatty(2)) exit_on_error(false); // if we stdin or stderr is redirected - we cannot ask the user an interactive question
    
    // read all chars available on stdin, so that if we're processing multiple files - and we ask this question
    // for a subsequent file later - we don't get left overs of this response
    char read_buf[1000];
    str_query_user ("Do you wish to overwrite it now? (y or [n]) ", read_buf, sizeof(read_buf), str_verify_y_n, "N");

    if (read_buf[0] == 'N') {
        fprintf (stderr, "No worries, I'm stopping here - no damage done!\n");
        exit(0);
    }
}

static void file_redirect_output_to_stream (File *file, char *exec_name, char *stdout_option, char *format_option)
{
    char threads_str[20];
    sprintf (threads_str, "%u", global_max_threads);

    FILE *redirected_stdout_file = NULL;
    if (!flag_stdout) {
        redirected_stdout_file = fopen (file->name, file->mode); // exec_name will redirect its output to this file
        ASSINP (redirected_stdout_file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno));
    }
    char reason[100];
    sprintf (reason, "To output a %s file", file_exts[file->type]);
    output_compressor = stream_create (0, 0, 0, global_max_memory_per_vb, 
                                        redirected_stdout_file, // output is redirected unless flag_stdout
                                        0, reason,
                                        exec_name, 
                                        stdout_option, // either to the terminal or redirected to output file
                                        "--threads", threads_str,
                                        format_option,
                                        NULL);
    file->file = stream_to_stream_stdin (output_compressor);
}

// returns true if successful
bool file_open_txt (File *file)
{
    ASSERT0 (file->mode == READ || file->mode == WRITE, "Error in file_open_txt: only READ and WRITE modes are supported");

    // for READ, set data_type
    if (file->mode == READ) {
        
        // if user provided the type with --input, we use that overriding the type derived from the file name
        if (stdin_type) file->type = stdin_type;

        file->data_type = file_get_data_type (file->type, true);

        if (file->data_type == DT_NONE) { // show meaningful error if file is not a supported type
            if (flag_multiple_files) {
                RETURNW (!file_has_ext (file->name, ".genozip"), true,
                        "Skipping %s - it is already compressed", file_printname(file));

                RETURNW (false, true, "Skipping %s - genozip doesn't know how to compress this file type (use --input to tell it)", 
                        file_printname (file));
            }
            else {
                ASSINP (!file_has_ext (file->name, ".genozip"), 
                        "%s: cannot compress %s because it is already compressed", global_cmd, file_printname(file));

                ABORT ("%s: the type of data in %s cannot be determined by its file name extension.\nPlease use --input (or -i) to specify one of the following types, or provide an input file with an extension matching one of these types.\n\nSupported file types: %s", 
                       global_cmd, file_printname (file),  file_compressible_extensions());
            }
        }

        ASSINP0 (!flag_make_reference || file->data_type == DT_REF, "Error: --make-reference can only be used with FASTA files");
    }
    else { // WRITE - data_type is already set by file_open

        if (file->data_type != DT_NONE && file->data_type < NUM_DATATYPES) {    
            
            // if the file is not a recognized output file
            if (file_get_data_type (file->type, false) == DT_NONE) {

                // case: output file has a .gz extension (eg the user did "genounzip xx.vcf.genozip -o yy.gz")
                if (file_has_ext (file->name, ".gz") &&  
                    file_has_ext (file_exts[txt_out_ft_by_dt[file->data_type][1]], ".gz")) // data type supports .gz txt output
                    
                    file->type = txt_out_ft_by_dt[file->data_type][1]; 
                
                // case: not .gz - use the default plain file format
                else 
                    file->type = txt_out_ft_by_dt[file->data_type][0]; 
            }
        }

        // if we don't know our type based on our own name, consult z_file (this happens when opening 
        // from piz_dispacher)
        else if (file->data_type == DT_NONE && z_file && z_file->data_type != DT_NONE) {
            file->data_type = z_file->data_type;

            #define FORBID_THIS_FLAG(flag,dt) ASSINP (!flag_##flag, "%s: the --" #flag " flag cannot be used with files containing %s data like %s", global_cmd, dt, z_name);
            
            switch (file->data_type) {
                case DT_VCF : 
                    FORBID_THIS_FLAG (bam, "VCF");
                    file->type = (flag_bgzip ? VCF_GZ : (flag_bcf ? BCF : VCF)); 
                    break;
                
                case DT_SAM : 
                    FORBID_THIS_FLAG (bcf, "SAM");
                    FORBID_THIS_FLAG (bgzip, "SAM");
                    file->type = (flag_bam ? BAM : SAM); 
                    break;
                
                case DT_FASTQ : 
                    FORBID_THIS_FLAG (bcf, "FASTQ");
                    FORBID_THIS_FLAG (bam, "FASTQ");
                    file->type = (flag_bgzip ? FASTQ_GZ : FASTQ); 
                    break;
                
                case DT_FASTA : 
                    FORBID_THIS_FLAG (bcf, "FASTA");
                    FORBID_THIS_FLAG (bam, "FASTA");
                    file->type = (flag_bgzip ? FASTA_GZ : FASTA); 
                    break;
                
                case DT_ME23 : 
                    FORBID_THIS_FLAG (bcf, "23andMe");
                    FORBID_THIS_FLAG (bam, "23andMe");
                    FORBID_THIS_FLAG (bgzip, "23andMe");
                    file->type = (flag_bgzip ? FASTA_GZ : FASTA); 
                    break;
                
                default: ABORT ("Error in file_open_txt: unknown data_type=%s", dt_name (file->data_type));
            }
        }
    }

    // open the file, based on the compression algorithm
    file->comp_alg = file_get_comp_alg_by_txt_ft (file->data_type, file->type, file->mode);

    switch (file->comp_alg) { 
        case COMP_NONE:
            // don't actually open the output file if we're just testing in genounzip or PIZing a reference file
            if ((flag_test || flag_reading_reference) && file->mode == WRITE) return true;

            file->file = file->is_remote ? url_open (NULL, file->name) : fopen (file->name, file->mode);
            break;

        case COMP_GZ:
            if (file->mode == READ) {
                if (file->is_remote) { 
                    FILE *url_fp = url_open (NULL, file->name);
                    file->file = gzdopen (fileno(url_fp), file->mode); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
                }
                else
                    file->file = gzopen64 (file->name, file->mode); // for local files we decompress ourselves
            }
            else 
                file_redirect_output_to_stream (file, "bgzip", "--stdout", NULL);
            break;
        
        case COMP_BZ2:
            if (file->is_remote) {            
                FILE *url_fp = url_open (NULL, file->name);
                file->file = BZ2_bzdopen (fileno(url_fp), file->mode); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
            }
            else
                file->file = BZ2_bzopen (file->name, file->mode);  // for local files we decompress ourselves   
            break;

        case COMP_XZ:
            input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,     // url
                                                "To uncompress an .xz file", "xz",       // reason, exec_name
                                                file->is_remote ? SKIP_ARG : file->name, // local file name 
                                                "--threads=8", "--decompress", "--keep", "--stdout", 
                                                flag_quiet ? "--quiet" : SKIP_ARG, 
                                                NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case COMP_ZIP:
            input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,     // url
                                                "To uncompress a .zip file", "unzip",    // reason, exec_name
                                                "-p", // must be before file name
                                                file->is_remote ? SKIP_ARG : file->name, // local file name 
                                                flag_quiet ? "--quiet" : SKIP_ARG, 
                                                NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case COMP_BCF:
        case COMP_BAM: 
        case COMP_CRAM: {
            bool bam =  (file->comp_alg == COMP_BAM);
            bool cram = (file->comp_alg == COMP_CRAM);

            if (file->mode == READ) {
                char reason[100];
                sprintf (reason, "To compress a %s file", file_exts[file->type]);
                input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                                    file->is_remote ? file->name : NULL,         // url                                        
                                                    reason, 
                                                    (bam || cram) ? "samtools" : "bcftools", // exec_name
                                                    "view", 
                                                    "--threads", "8", // in practice, samtools is able to consume 1.3 cores
                                                    (bam || cram) ? "-OSAM" : "-Ov",
                                                    file->is_remote ? SKIP_ARG : file->name,    // local file name 
                                                    (bam || cram) ? "-h" : "--no-version", // BAM: include header
                                                                                           // BCF: do not append version and command line to the header
                                                    cram ? ref_get_cram_ref() : NULL,
                                                    NULL);
                file->file = stream_from_stream_stdout (input_decompressor);
            }
            else { // write
                file_redirect_output_to_stream (file, 
                                                bam ? "samtools" : "bcftools", 
                                                "view", 
                                                bam ? "-OBAM" : "-Ob");            
            }
            break;
        }

        default:
            if (file->mode == WRITE && file->data_type == DT_NONE) {
                // user is trying to genounzip a .genozip file that is not a recognized extension -
                // that's ok - we will discover its type after reading the genozip header in zip_dispatcher
                return true; // let's call this a success anyway
            }
            else {
                ABORT ("%s: unrecognized file type: %s", global_cmd, file->name);
            }
    }

    if (file->mode == READ) 
        file->txt_data_size_single = file->disk_size; 

    return file->file != 0;
}

// we insert all the z_file buffers into the buffer list in advance, to avoid this 
// thread satety issue:
// without our pre-allocation, some of these buffers will be first allocated by a compute threads 
// when the first vb containing a certain did_i is merged in (for the contexts buffers) or
// when the first dictionary is compressed (for dict_data) or ra is merged (for ra_buf). while these operations 
// are done while holding a mutex, so that compute threads don't run over each over, buf_alloc 
// may change buf_lists in evb buffers, while the I/O thread might be doing so concurrently
// resulting in data corruption in evb.buf_list. If evb.buf_list gets corrupted this might result in termination 
// of the execution.
// with these buf_add_to_buffer_list() the buffers will already be in evb's buf_list before any compute thread is run.
static void file_initialize_z_file_data (File *file)
{
    memset (file->dict_id_to_did_i_map, DID_I_NONE, sizeof(file->dict_id_to_did_i_map));

#define INIT(buf) file->contexts[i].buf.name  = #buf; \
                  file->contexts[i].buf.param = i;    \
                  buf_add_to_buffer_list (evb, &file->contexts[i].buf);

    for (unsigned i=0; i < MAX_DICTS; i++) {
        INIT (dict);        
        INIT (b250);
        INIT (mtf);
        INIT (mtf_i);
        INIT (global_hash);
        INIT (ol_dict);
        INIT (ol_mtf);
        INIT (local_hash);
        INIT (word_list);
    }

#undef INIT
#define INIT(buf) file->buf.name = #buf; \
                  buf_add_to_buffer_list (evb, &file->buf);
    INIT (dict_data);
    INIT (ra_buf);
    INIT (ra_min_max_by_chrom);
    INIT (chroms_sorted_index);
    INIT (alt_chrom_map);
    INIT (section_list_buf);
    INIT (section_list_dict_buf);
    INIT (unconsumed_txt);
}
#undef INIT

// returns true if successful
static bool file_open_z (File *file)
{
    // for READ, set data_type
    if (file->mode == READ) {

        if (!file_has_ext (file->name, GENOZIP_EXT)) {
            if (flag_multiple_files) 
                RETURNW (false, true, "Skipping %s - it doesn't have a .genozip extension", file_printname (file))
            else {
                if (flag_reading_reference)
                    ABORT ("%s: with --reference or --REFERENCE, you must specify a genozip reference file (.ref.genozip extension)\n"
                           "Tip: You can create a genozip reference file from a FASTA file with 'genozip --make-reference myfasta.fa'",
                           global_cmd)
                else
                    ABORT ("%s: file %s must have a " GENOZIP_EXT " extension", global_cmd, file_printname (file));
            }
        }

        file->data_type = file_get_dt_by_z_ft (file->type); // if we don't find it (DT_NONE), it will be determined after the genozip header is read
    }
    else { // WRITE or WRITEREAD - data_type is already set by file_open
        ASSINP (file_has_ext (file->name, GENOZIP_EXT), "%s: file %s must have a " GENOZIP_EXT " extension", 
                              global_cmd, file_printname (file));
        // set file->type according to the data type, overriding the previous setting - i.e. if the user
        // uses the --output option, he is unrestricted in the choice of a file name
        file->type = file_get_z_ft_by_txt_in_ft (file->data_type, txt_file->type); 
    }

    ASSERT (!file->is_remote, "Error: it is not possible to access remote genozip files; when attempting to open %s", file->name);
    
    file->file = fopen (file->name, file->mode);

    // initialize read buffer indices
    if (file->mode == READ || file->mode == WRITEREAD) 
        file->z_last_read = file->z_next_read = READ_BUFFER_SIZE;

    file_initialize_z_file_data (file);

    return file->file != 0;
}

File *file_open (const char *filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE or WRITEREAD */)
{
    ASSINP0 (filename, "Error in file_open: filename is null");

    File *file = (File *)calloc (1, sizeof(File) + (((mode == READ || mode == WRITEREAD) && supertype == Z_FILE) ? READ_BUFFER_SIZE : 0));

    file->supertype = supertype;
    file->is_remote = url_is_url (filename);
    bool file_exists;

    // is_remote is only possible in READ mode
    ASSINP (mode == READ || !file->is_remote, "%s: expecting output file %s to be local, not a URL", global_cmd, filename);

    int64_t url_file_size = 0; // will be -1 if the web/ftp site does not provide the file size
    const char *error = NULL;
    if (file->is_remote) {
        error = url_get_status (filename, &file_exists, &url_file_size); // accessing is expensive - get existance and size in one call
        if (url_file_size >= 0) file->disk_size = (uint64_t)url_file_size;
    }
    else {
        file_exists = (access (filename, F_OK) == 0);
        error = strerror (errno);
        if (file_exists && mode == READ) file->disk_size = file_get_size (filename);
    }

    // return null if genozip input file size is known to be 0, so we can skip it. note: file size of url might be unknown
    if (mode == READ && supertype == TXT_FILE && file_exists && !file->disk_size && !url_file_size) {
        FREE (file);
        return NULL; 
    }

    ASSINP (mode != READ || file_exists, "%s: cannot open '%s' for reading: %s", global_cmd, filename, error);

    if ((mode == WRITE || mode == WRITEREAD) && file_exists && !flag_force && !(supertype==TXT_FILE && flag_test))
        file_ask_user_to_confirm_overwrite (filename); // function doesn't return if user responds "no"

    // copy filename 
    unsigned fn_size = strlen (filename) + 1; // inc. \0
    file->name = malloc (fn_size);
    memcpy (file->name, filename, fn_size);

    file->mode = mode;

    if (mode==READ || data_type != DT_NONE) // if its NONE, we will not open now, and try again from piz_dispatch after reading the genozip header
        file->type = file_get_type (file->name, true);

    if (file->mode == WRITE || file->mode == WRITEREAD) 
        file->data_type = data_type; // for WRITE, data_type is set by file_open_*

    bool success=false;
    switch (supertype) {
        case TXT_FILE: success = file_open_txt (file); break;
        case Z_FILE:   success = file_open_z (file);   break;
        default:       ABORT ("Error: invalid supertype: %u", supertype);
    }

    ASSINP (success, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    return file;
}

File *file_open_redirect (FileMode mode, FileSupertype supertype, DataType data_type /* only used for WRITE */)
{
    ASSINP (mode==WRITE || stdin_type != UNKNOWN_FILE_TYPE, 
            "%s: to redirect from standard input use --input (or -i) with one of the supported file types:%s", 
            global_cmd, file_compressible_extensions());

    File *file = (File *)calloc (1, sizeof(File) + ((mode == READ && supertype == Z_FILE) ? READ_BUFFER_SIZE : 0));

    file->file = (mode == READ) ? fdopen (STDIN_FILENO,  "rb")
                                : fdopen (STDOUT_FILENO, "wb");
    ASSINP (file->file, "%s: Failed to redirect %s: %s", global_cmd, (mode==READ ? "stdin" : "stdout"), strerror (errno));

    file->supertype = supertype;
    
    if (mode==READ) {
        file->data_type = file_get_data_type (stdin_type, true);
        file->type = stdin_type; 
    }
    else { // WRITE
        file->data_type = data_type;
        file->type = txt_out_ft_by_dt[data_type][0];
    }

    if (supertype == Z_FILE)
        file_initialize_z_file_data (file);

    file->redirected = true;
    
    return file;
}

void file_close (File **file_p, 
                 bool cleanup_memory) // optional - used to destroy buffers in the file is closed NOT near the end of the execution, eg when dealing with unbinding bound files
{
    File *file = *file_p;
    *file_p = NULL;

    if (!file) return; // nothing to do
    
    if (file->file) {

        if (file->mode == READ && (file->comp_alg == COMP_GZ || file->comp_alg == COMP_BGZ)) {
            int ret = gzclose_r((gzFile)file->file);
            ASSERTW (!ret, "%s: warning: failed to close file: %s", global_cmd, file_printname (file));
        }
        else if (file->mode == READ && file->comp_alg == COMP_BZ2)
            BZ2_bzclose((BZFILE *)file->file);
        
        else if (file->mode == READ && file_is_read_via_ext_decompressor (file)) 
            stream_close (&input_decompressor, STREAM_WAIT_FOR_PROCESS);

        else if (file->mode == WRITE && file_is_written_via_ext_compressor (file))
            stream_close (&output_compressor, STREAM_WAIT_FOR_PROCESS);

        else 
            FCLOSE (file->file, file_printname (file));
    }

    // free resources if we are NOT near the end of the execution. If we are at the end of the execution
    // it is faster to just let the process die
    if (cleanup_memory) {

        mutex_destroy (file->dicts_mutex);
            
        // always destroy all buffers even if unused - for saftey
        for (unsigned i=0; i < MAX_DICTS; i++) // we need to destory all even if unused, because they were initialized in file_initialize_z_file_data
            mtf_destroy_context (&file->contexts[i]);

        buf_destroy (&file->dict_data);
        buf_destroy (&file->ra_buf);
        buf_destroy (&file->ra_min_max_by_chrom);
        buf_destroy (&file->chroms_sorted_index);
        buf_destroy (&file->alt_chrom_map);
        buf_destroy (&file->section_list_buf);
        buf_destroy (&file->section_list_dict_buf);
        buf_destroy (&file->unconsumed_txt);

        FREE (file->name);
        FREE (file->basename);
        FREE (file);
    }
}

size_t file_write (File *file, const void *data, unsigned len)
{
    if (!len) return 0; // nothing to do

    size_t bytes_written = fwrite (data, 1, len, (FILE *)file->file); // use fwrite - let libc manage write buffers for us

    // if we're streaming our genounzip/genocat/genols output to another process and that process has 
    // ended prematurely then exit quietly. in genozip we display an error because this means the resulting
    // .genozip file will be corrupted
    if (!file->name && command != ZIP && errno == EINVAL) exit(0);

    ASSERT (bytes_written, "Error: failed to write %u bytes to %s: %s", 
            len, file->name ? file->name : FILENAME_STDOUT, strerror(errno));
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
    
    return fn_len >= ext_len && !strncmp (&filename[fn_len-ext_len], extension, ext_len);
}

// get basename of a filename - we write our own basename for Visual C and Windows compatibility
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
    ASSERT0 (file->supertype == Z_FILE, "Error: file_seek only works for Z_FILE supertype");

    // check if we can just move the read buffers rather than seeking
    if (file->mode == READ && file->z_next_read != file->z_last_read && whence == SEEK_SET) {
#ifdef __APPLE__
        int64_t move_by = offset - ftello ((FILE *)file->file);
#else
        int64_t move_by = offset - ftello64 ((FILE *)file->file);
#endif

        // case: move is within read buffer already in memory (ftello shows the disk location after read of the last buffer)
        // we just change the buffer pointers rather than discarding the buffer and re-reading
        if (move_by <= 0 && move_by >= -(int64_t)file->z_last_read) {
            file->z_next_read = file->z_last_read + move_by; // only z_file uses next/last_Read
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
    if (!ret || file->mode == WRITEREAD) file->z_next_read = file->z_last_read = READ_BUFFER_SIZE;

    return !ret;
}

uint64_t file_tell (File *file)
{
    if (command == ZIP && file == txt_file && file->comp_alg == COMP_GZ)
        return gzconsumed64 ((gzFile)txt_file->file); 
    
    if (command == ZIP && file == txt_file && file->comp_alg == COMP_BZ2)
        return BZ2_consumed ((BZFILE *)txt_file->file); 

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
    ASSERT (!ret, "Error: failed accessing '%s': %s", filename, strerror(errno));
    
    return st.st_size;
}

bool file_is_dir (const char *filename)
{
    struct stat64 st;
    int ret = stat64(filename, &st); // 0 if successful
    
    return !ret && S_ISDIR (st.st_mode);
}

void file_get_file (VBlockP vb, const char *filename, Buffer *buf, const char *buf_name, unsigned buf_param,
                    bool add_string_terminator)
{
    uint64_t size = file_get_size (filename);

    buf_alloc (vb, buf, size + add_string_terminator, 1, buf_name, buf_param);

    FILE *file = fopen (filename, "rb");
    ASSINP (file, "Error: cannot open %s: %s", filename, strerror (errno));

    size_t bytes_read = fread (buf->data, 1, size, file);
    ASSERT (bytes_read == (size_t)size, "Error reading file %s: %s", filename, strerror (errno));

    buf->len = size;

    if (add_string_terminator) buf->data[size] = 0;

    FCLOSE (file, filename);
}

void file_assert_ext_decompressor (void)
{
    if (!stream_wait_for_exit (input_decompressor)) return; // just a normal EOF - all good!

    // read error from stderr
    #define INPUT_DECOMPRESSOR_RESPSONSE_LEN 4096
    char error_str[INPUT_DECOMPRESSOR_RESPSONSE_LEN];

    FILE *stderr_pipe = stream_from_stream_stderr (input_decompressor);
    int bytes_read = fread (error_str, 1, INPUT_DECOMPRESSOR_RESPSONSE_LEN-1, stderr_pipe);
    error_str[bytes_read] = 0; // string terminator

    ABORT ("%s: failed to read file: %s", global_cmd, error_str);
}

// used when aborting due to an error. avoid the compressors outputting their own errors after our process is gone
void file_kill_external_compressors (void)
{
    stream_close (&input_decompressor, STREAM_KILL_PROCESS);
    stream_close (&output_compressor,  STREAM_KILL_PROCESS);
}

const char *ft_name (FileType ft)
{
    return type_name (ft, &file_exts[ft], sizeof(file_exts)/sizeof(file_exts[0]));
}

const char *file_viewer (const File *file)
{
    static const char *viewer[NUM_COMPRESSION_ALGS] = { 
        /* none */ "cat", 
        /* gz   */ "gunzip -c", 
        /* bz2  */ "bzip2 -d -c", 
        "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", 
        /* bgz  */ "gunzip -c",   
        /* xz   */ "xz -d -c",    
        /* bcf  */ "bcftools view", 
        /* bam  */ "samtools view -h --threads 2",
        /* cram */ "samtools view -h --threads 2",
        /* zip  */ "unzip -p" };  // zip

    return viewer[file->comp_alg];
}

// PIZ: guess original filename from uncompressed txt filename and compression algoritm (allocated memory)
const char *file_guess_original_filename (const File *file)
{
    static const char *comp_alg_exts[NUM_COMPRESSION_ALGS] = COMP_ALG_EXTS; 

    if (file->comp_alg == COMP_NONE) return file->name;

    unsigned len = strlen (file->name) + 10;
    char *org_name = malloc (len);
    strcpy (org_name, file->name);

    // remove existing extension if needed (eg when replacing .sam with .bam)
    if (comp_alg_exts[file->comp_alg][0] == '-') {
        char *last_dot = strrchr (org_name, '.');
        if (last_dot) *last_dot = 0;
    }

    // add new extension
    strcpy (&org_name[strlen(org_name)], &comp_alg_exts[file->comp_alg][1]);

    return org_name;
}

const char *file_plain_ext_by_dt (DataType dt)
{
    FileType plain_ft = txt_in_ft_by_dt[dt][0].in;

    return file_exts[plain_ft];
}
