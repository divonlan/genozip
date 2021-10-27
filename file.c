// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#define Z_LARGE64
#ifdef __APPLE__
    #define off64_t __int64_t
#endif
#include "bzlib/bzlib.h"
#include "zlib/zlib.h"
#include "genozip.h"
#include "context.h"
#include "file.h"
#include "stream.h"
#include "url.h"
#include "vblock.h"
#include "strings.h"
#include "codec.h"
#include "mutex.h"
#include "bgzf.h"
#include "flags.h"
#include "progress.h"
#include "endianness.h"
#include "tar.h"
#include "mutex.h"

// globals
File *z_file   = NULL;
File *txt_file = NULL;

static StreamP input_decompressor  = NULL; // bcftools or xz, unzip or samtools - only one at a time
static StreamP output_compressor   = NULL; // samtools (for cram), bcftools

static FileType stdin_type = UNKNOWN_FILE_TYPE; // set by the --input command line option

// global pointers - so the can be compared eg "if (mode == READ)"
const char *READ  = "rb";  // use binary mode (b) in read and write so Windows doesn't add \r
const char *WRITE = "wb";
const char *WRITEREAD = "wb+"; // only supported for z_file

const char *file_exts[] = FILE_EXTS;

static const struct { FileType in; Codec codec; FileType out; } txt_in_ft_by_dt[NUM_DATATYPES][50] = TXT_IN_FT_BY_DT;
static const FileType txt_out_ft_by_dt[NUM_DATATYPES][20] = TXT_OUT_FT_BY_DT;
static const FileType z_ft_by_dt[NUM_DATATYPES][20] = Z_FT_BY_DT;

static void file_initialize_txt_file_data (File *file);

// get data type by file type
DataType file_get_data_type (FileType ft, bool is_input)
{
    // note: if make-reference, we scan the array from dt=0 (DT_REF), otherwise we ignore DT_REF
    for (DataType dt=!flag.make_reference; dt < NUM_DATATYPES; dt++) 
        for (unsigned i=0; (is_input ? txt_in_ft_by_dt[dt][i].in : txt_out_ft_by_dt[dt][i]); i++)
            if ((is_input ? txt_in_ft_by_dt[dt][i].in : txt_out_ft_by_dt[dt][i]) == ft) 
                return dt;

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

// get codec by txt file type
static Codec file_get_codec_by_txt_ft (DataType dt, FileType txt_ft)
{
    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].in == txt_ft) 
            return txt_in_ft_by_dt[dt][i].codec;

    return CODEC_NONE;
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
char *file_compressible_extensions (bool plain_only)
{
    static char s[1000] = {0};
        
    for (DataType dt=1; dt < NUM_DATATYPES; dt++) { // start from 1, excluding DT_REFERENCE
        
        if (dt == DT_GENERIC || dt == DT_ME23) continue;

        if (plain_only) 
            sprintf (&s[strlen(s)], "%s ", &file_exts[txt_in_ft_by_dt[dt][0].in][1]);
    
        else {
            sprintf (&s[strlen (s)], "\n%-8s: ", dt_name (dt));

            for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
                sprintf (&s[strlen(s)], "%s ", &file_exts[txt_in_ft_by_dt[dt][i].in][1]);
        }
    }

    if (plain_only)
        sprintf (&s[strlen(s)], "23andme generic ");
    else
        sprintf (&s[strlen(s)], "\n23andMe : 23andme 23andme.zip"
                                "\nOther   : generic");
    return s;
}

FileType file_get_type (const char *filename)
{
    if (!filename) return UNKNOWN_FILE_TYPE;

    // 23andme files have the format "genome_Firstname_Lastname_optionalversion_timestamp.txt" or .zip
    if (strstr (filename, "genome") && strstr (filename, "Full")) {
        if (file_has_ext (filename, ".txt")) return ME23;
        if (file_has_ext (filename, ".zip")) return ME23_ZIP;
        if (file_has_ext (filename, ".txt.genozip")) return ME23_GENOZIP;
    }

    for (FileType ft=UNKNOWN_FILE_TYPE+1; ft < AFTER_LAST_FILE_TYPE; ft++) {

        // files that end with .txt/.txt.genozip/.zip are not classified as ME23, we already handled ME23 above
        if (ft == ME23 || ft == ME23_ZIP || ft == ME23_GENOZIP) continue; 

        if (file_has_ext (filename, file_exts[ft])) 
            return ft;
    }

    return UNKNOWN_FILE_TYPE; // this should never happen, because GNRIC_ is "" so it will catch everything
}

// returns the filename without the extension eg myfile.1.sam.gz -> myfile.1. 
// if raw_name is given, memory is allocated sufficiently to concatenate a extension. Otherwise, filename is overwritten
void file_get_raw_name_and_type (char *filename, char **raw_name, FileType *out_ft)
{
    unsigned len = strlen (filename);

    if (raw_name) {
        *raw_name = MALLOC (len + 30);
        memcpy (*raw_name, filename, len);
        (*raw_name)[len] = 0;
    }
    else 
        raw_name = &filename; // overwrite filename

    FileType ft = file_get_type (filename);
    if (ft != UNKNOWN_FILE_TYPE) 
        (*raw_name)[len - strlen (file_exts[ft])] = 0;

    if (out_ft) *out_ft = ft;
}

// rules for pair filename: two filenames need to differ by exactly one character, which is '1' and '2',
// if test_only we validate the rule
// if !test_only, only fn1 needs to be provided, we assume the rule is valid, and we create the output file with "1+2"
char *file_get_fastq_pair_filename (const char *fn1, const char *fn2, bool test_only)
{
    FileType ft1, ft2;
    char *rn1, *rn2;

    file_get_raw_name_and_type ((char *)fn1, &rn1, &ft1);
    file_get_raw_name_and_type ((char *)fn2, &rn2, &ft2);
    
    unsigned len = strlen (rn1);
    if (len != strlen (rn2)) return NULL;

    int df = -1;
    for (unsigned i=0; i < len; i++)
        if (rn1[i] != rn2[i]) {
            if (df >= 0) return NULL; // 2nd differing character
            df = i;
        }

    if (!((rn1[df] == '1' && rn2[df] == '2') || (rn1[df] == '2' && rn2[df] == '1'))) return NULL; // one of them must be '1' and the other '2'

    if (test_only) return (char *)1;

    char *pair_fn = MALLOC (len+20);
    sprintf (pair_fn, "%.*s1+2%s" FASTQ_GENOZIP_, df, rn1, &rn1[df+1]);
    
    return pair_fn;
}

void file_set_input_size (const char *size_str)
{   
    unsigned len = strlen (size_str);

    for (unsigned i=0; i < len; i++) {
        ASSINP (IS_DIGIT (size_str[i]), "expecting the file size in bytes to be a positive integer: %s", size_str);

        flag.stdin_size = flag.stdin_size * 10 + (size_str[i] - '0'); 
    }
}

void file_set_input_type (const char *type_str)
{
    char ext[strlen (type_str) + 2]; // +1 for . +1 for \0
    sprintf (ext, ".%s", type_str);

    str_tolower (ext, ext); // lower-case to allow case-insensitive --input argument (eg vcf or VCF)

    if (!strcmp (ext, ".23andme")) 
        stdin_type = ME23;

    else if (!strcmp (ext, ".23andme.zip")) 
        stdin_type = ME23_ZIP;

    else if (!strcmp (ext, ".generic"))
        stdin_type = GNRIC;

    else {
        stdin_type = file_get_type (ext); // we don't enforce 23andMe name format - any .txt or .zip will be considered ME23
        if (stdin_type == GNRIC) stdin_type = UNKNOWN_FILE_TYPE; // This is an unknown input name - not generic
    }

    if (file_get_data_type (stdin_type, true) != DT_NONE) return; // all good 

    // the user's argument is not an accepted input file type - print error message
    ABORT ("%s: --input (or -i) must be ones of these: %s", global_cmd, file_compressible_extensions(false));
}
    
FileType file_get_stdin_type (void)
{
    return stdin_type;
}

static void file_ask_user_to_confirm_overwrite (const char *filename)
{
    if (!strcmp (filename, "/dev/null")) return; // don't ask for /dev/null
    
    fprintf (stderr, "%s: output file %s already exists: in the future, you may use --force to overwrite\n", global_cmd, filename);
    
    if (!isatty(0) || !isatty(2)) exit_on_error(false); // if we stdin or stderr is redirected - we cannot ask the user an interactive question
    
    // read all chars available on stdin, so that if we're processing multiple files - and we ask this question
    // for a subsequent file later - we don't get left overs of this response
    char read_buf[1000];
    str_query_user ("Do you wish to overwrite it now? (y or [n]) ", read_buf, sizeof(read_buf), str_verify_y_n, "N");

    if (read_buf[0] == 'N') {
        fprintf (stderr, "No worries, I'm stopping here - no damage done!\n");
        exit (EXIT_OK);
    }

    fprintf (stderr, "\n");
}

static void file_redirect_output_to_stream (File *file, char *exec_name, 
                                            const char *stdout_option, const char *format_option_1, const char *format_option_2, const char *format_option_3)
{
    char threads_str[20];
    sprintf (threads_str, "%u", global_max_threads);

    FILE *redirected_stdout_file = NULL;
    if (!flag.to_stdout) {
        redirected_stdout_file = fopen (file->name, file->mode); // exec_name will redirect its output to this file
        ASSINP (redirected_stdout_file, "cannot open file \"%s\": %s", file->name, strerror(errno));
    }
    char reason[100];
    sprintf (reason, "To output a %s file", file_exts[file->type]);
    output_compressor = stream_create (0, 0, 0, DEFAULT_PIPE_SIZE, 
                                       redirected_stdout_file, // output is redirected unless flag.to_stdout
                                       0, false, reason,
                                       exec_name, 
                                       stdout_option, // either to the terminal or redirected to output file
                                       "--threads", threads_str,
                                       format_option_1, format_option_2, format_option_3,
                                       NULL);
    file->file = stream_to_stream_stdin (output_compressor);
}

// starting samtools 1.10, a PG record is added to the SAM header every "samtools view", and an option, --no-PG, 
// is provided to avoid this. See: https://github.com/samtools/samtools/releases/
// returns "--no-PG" if the option exists, or NULL if not
static const char *file_samtools_no_PG (void)
{
    static const char *ret_str[2] = { NULL, "--no-PG" };
    static int has_no_PG = -1; // unknown

    if (has_no_PG >= 0) return ret_str[has_no_PG]; // we already tested

    #define SAMTOOLS_HELP_MAX_LEN 20000
    char samtools_help_text[SAMTOOLS_HELP_MAX_LEN];

    #define MIN_ACCEPTABLE_LEN 100
    int len=0;

    for (unsigned i=1; i < 15 && len < MIN_ACCEPTABLE_LEN; i++) {
        // Tested on samtools 1.11: The normal way to see help is "samtools help view" however it fails if stdin is not the terminal. 
        // Instead, we use samtools view with an invalid option "--junk". This *sometimes* shows the help, and sometimes
        // just shows one line "samtools view:". We overcome this by repeating if the response is not long enough.
        StreamP samtools = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, 0, "To read/write CRAM files",
                                          "samtools", "view", "--junk", NULL);
        usleep (50000 * i); // wait for samtools

        // read both stderr and stdout from samtools
        len  = read (fileno (stream_from_stream_stderr (samtools)), samtools_help_text,       SAMTOOLS_HELP_MAX_LEN-1);
        len += read (fileno (stream_from_stream_stdout (samtools)), &samtools_help_text[len], SAMTOOLS_HELP_MAX_LEN-len-1);

        stream_close (&samtools, STREAM_DONT_WAIT_FOR_PROCESS);
    } 
    samtools_help_text[len] = '\0'; // terminate string (more portable, strnstr and memmem are non-standard)

    ASSERT0 (len >= MIN_ACCEPTABLE_LEN, "no response from \"samtools view --junk\"");

    return ret_str[(has_no_PG = !!strstr (samtools_help_text, "--no-PG"))];
}

// show meaningful error if file is not a supported type and return TRUE if it file should be skipped
static bool file_open_txt_read_test_valid_dt (const File *file)
{ 
    if (file->data_type == DT_NONE) { 

        if (flag.multiple_files || tar_is_tar()) {
            if (file_has_ext (file->name, ".genozip")) {
            
                // case: compressing into a tar file - include .genozip files verbatim
                if (tar_is_tar()) {
                    tar_copy_file (file->name);
                    RETURNW (false, true, "Copied %s to tar file", file_printname(file));
                }    
                else
                    RETURNW (false, true, "Skipping %s - it is already compressed", file_printname(file));

            }

            RETURNW (false, true, "Skipping %s - genozip doesn't know how to compress this file type (use --input to tell it)", 
                    file_printname (file));
        }
        else {
            ASSINP (!file_has_ext (file->name, ".genozip"), 
                    "cannot compress %s because it is already compressed", file_printname(file));

            if (file->redirected)
                ABORT ("%s: to pipe data in, please use --input (or -i) to specify its type, which can be one of the following:\n%s", 
                        global_cmd, file_compressible_extensions (true));
            else
                ABORT ("%s: the type of data in %s cannot be determined by its file name extension.\nPlease use --input (or -i) to specify one of the following types, or provide an input file with an extension matching one of these types.\n\nSupported file types: %s", 
                        global_cmd, file_printname (file),  file_compressible_extensions (false));
        }
    }

    ASSINP0 (!flag.make_reference || file->data_type == DT_REF, "--make-reference can only be used with FASTA files");

    return false; // all good - no need to skip this file
}
static bool file_open_txt_read (File *file)
{
    // if user provided the type with --input, we use that overriding the type derived from the file name
    if (stdin_type) file->type = stdin_type;

    file->data_type = file_get_data_type (file->type, true);

    // show meaningful error if file is not a supported data type
    if (file_open_txt_read_test_valid_dt (file)) return true; // skip this file

    if (file->data_type == DT_GENERIC && !stdin_type) 
        WARN_ONCE ("FYI: genozip can't recognize %s file's type, so it will be compressed as GENERIC. In the future, you may specify the type with \"--input <type>\". To suppress this warning, use \"--input generic\".", file->name);
    
    // open the file, based on the codec
    file->codec = file_get_codec_by_txt_ft (file->data_type, file->type);

    ASSINP0 (!flag.is_windows || !file->redirected || file->codec == CODEC_NONE, 
             "genozip on Windows supports piping in only plain (uncompressed) data");

    switch (file->codec) { 
        case CODEC_GZ:   // we test the first few bytes of the file to differentiate between NONE, GZ and BGZIP
        case CODEC_BGZF: 
        case CODEC_NONE: {
            file->file = file->is_remote                  ? url_open (NULL, file->name)  
                       : file->redirected                 ? fdopen (STDIN_FILENO,  "rb") 
                       : flag.rejects_coord == DC_PRIMARY ? z_file->rejects_file[0]  // already open
                       : flag.rejects_coord == DC_LUFT    ? z_file->rejects_file[1]  // already open
                       :                                    fopen (file->name, READ);
            ASSERT (file->file, "failed to open %s: %s", file->name, strerror (errno));

            if (flag.rejects_coord)
                file_seek (file, 0, SEEK_SET, false); // we finished writing the reject file - now we will read it

            // read the first potential BGZF block to test if this is GZ or BGZF
            uint8_t block[BGZF_MAX_BLOCK_SIZE]; 
            uint32_t block_size;

            // we can't use libc read buffer unfortunately, because when reading the test bgzf block, if it turns out to be GZIP 
            // we will do gzdopen discarding the internal FILE buffer. We also cannot use setvbuf after the file is already partly read (?).
            // TO DO: this doesn't work over a pipe - so we can't pipe-in non-BGZF GZIP files (bug 243)
            setvbuf (file->file, 0, _IONBF, 0); 

            int32_t bgzf_uncompressed_size = bgzf_read_block (file, block, &block_size, true);

            // case: this is indeed a bgzf - we put the still-compressed data in vb->compressed for later consumption
            // is txtfile_read_block_bgzf
            if (bgzf_uncompressed_size > 0) {
                file->codec = CODEC_BGZF;
                
                evb->compressed.param = bgzf_uncompressed_size; // pass uncompressed size in param
                buf_add_more (evb, &evb->compressed, block, block_size, "compressed");
                
                file->bgzf_flags = bgzf_get_compression_level (file->name ? file->name : FILENAME_STDIN, 
                                                               block, block_size, (uint32_t)bgzf_uncompressed_size);
            }

            // for regulars files, we already skipped 0 size files. This can happen in STDIN
            else if (bgzf_uncompressed_size == 0) {

                ASSINP (!flags_pipe_in_process_died(), // only works for Linux
                        "Pipe-in process %s (pid=%u) died without sending any data",
                        flags_pipe_in_process_name(), flags_pipe_in_pid());

                ABORTINP0 ("No input data");
            }

            // case: this is a non-BGZF gzip format - open with zlib and hack back the read bytes 
            // (note: we cannot re-read the bytes from the file as the file might be piped in)
            else if (bgzf_uncompressed_size == BGZF_BLOCK_GZIP_NOT_BGZIP) {
                
                ASSERT0 (!file->redirected, "genozip can't read gzip data from a pipe - piped data must be either plain or in BGZF format - i.e. compressed with bgzip, htslib etc");

                file->codec = CODEC_GZ;
                file->file  = gzdopen (fileno((FILE *)file->file), READ); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
                gzinject (file->file, block, block_size); // a hack - adds a 18 bytes of compressed data to the in stream, which will be consumed next, instead of reading from disk
            } 

            // case: this is not GZIP format at all. treat as a plain file, and put the data read in vb->compressed 
            // for later consumption is txtfile_read_block_plain
            else if (bgzf_uncompressed_size == BGZF_BLOCK_IS_NOT_GZIP) {

                #define BZ2_MAGIC "BZh"
                #define XZ_MAGIC (char[]){ 0xFD, '7', 'z', 'X', 'Z', 0 }

                // we already open the file, so not easy to re-open with BZ2_bzopen as it would require injecting the read data into the BZ2 buffers
                if (block_size >= 3 && !memcmp (block, BZ2_MAGIC, 3)) 
                    ABORTINP0 ("The data seems to be in bz2 format. Please use --input to specify the type (eg: \"genozip --input sam.bz2\"");

                if (block_size >= 6 && !memcmp (block, XZ_MAGIC, 6)) {
                    if (file->redirected) ABORTINP0 ("Compressing piped-in data in xz format is not currently supported");
                    if (file->is_remote) ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\"");
                    ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\"");
                }

                file->codec = CODEC_NONE;
                buf_add_more (evb, &evb->compressed, block, block_size, "compressed");
            }

            // Notes 1. We currently only support plain and BGZF data piped in - see bug 243
            // 2. On Windows, we can't pipe binary files, bc Windows converts \n to \r\n
            // 3. The codec at this point is what the user declared in -i (we haven't tested for gz yet)
#ifdef _WIN32 
            ASSINP0 (!flag.is_windows || !file->redirected || file->codec == CODEC_NONE, 
                    "genozip on Windows supports piping in only plain (uncompressed) data");
#else
            ASSINP (!file->redirected || file->codec == CODEC_NONE || file->codec == CODEC_BGZF, 
                    "genozip only supports piping in data that is either plain (uncompressed) or compressed in BGZF format (typically with .gz extension) (codec=%s)", 
                    codec_name (file->codec));
#endif
            break;
        }
        case CODEC_BZ2:
            if (file->is_remote) {            
                FILE *url_fp = url_open (NULL, file->name);
                file->file = BZ2_bzdopen (fileno(url_fp), READ); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
            }
            else if (file->redirected) 
                file->file = BZ2_bzdopen (STDIN_FILENO, READ);
            
            else
                file->file = BZ2_bzopen (file->name, READ);  // for local files we decompress ourselves   
            break;

        case CODEC_XZ:
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in xz format is not currently supported");

            input_decompressor = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,     // url
                                                file->redirected,
                                                "To uncompress an .xz file", "xz",       // reason, exec_name
                                                file->is_remote ? SKIP_ARG : file->name, // local file name 
                                                "--threads=8", "--decompress", "--keep", "--stdout", 
                                                flag.quiet ? "--quiet" : SKIP_ARG, 
                                                NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case CODEC_ZIP:
            input_decompressor = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,     // url
                                                file->redirected,
                                                "To uncompress a .zip file", "unzip",    // reason, exec_name
                                                "-p", // must be before file name
                                                file->is_remote ? SKIP_ARG : file->name, // local file name 
                                                flag.quiet ? "--quiet" : SKIP_ARG, 
                                                NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case CODEC_BCF: {
            input_decompressor = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,         // url                                        
                                                file->redirected,
                                                "To compress a BCF file", 
                                                "bcftools", "view", "--threads", "8", "-Ov",
                                                file->is_remote ? SKIP_ARG : file->name,    // local file name 
                                                "--no-version", // BCF: do not append version and command line to the header
                                                NULL);
            file->file = stream_from_stream_stdout (input_decompressor);
            break;
        }

        case CODEC_CRAM: {
            input_decompressor = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,      // url                                        
                                                file->redirected,
                                                "To compress a CRAM file", 
                                                "samtools", "view", "-OSAM", "--threads", "8", "-h", // in practice, samtools is able to consume 1.3 cores
                                                file_samtools_no_PG() ? "--no-PG" : "-h", // don't add a PG line to the header (just repeat -h if this is an older samtools without --no-PG - no harm)
                                                file->is_remote ? SKIP_ARG : file->name,  // local file name 
                                                ref_get_cram_ref(gref), NULL);
            file->file = stream_from_stream_stdout (input_decompressor);
            break;
        }

        default:
            ABORT ("%s: invalid filename extension for %s files: %s", global_cmd, dt_name (file->data_type), file->name);
    }

    file_initialize_txt_file_data (file);

    return file->file != 0;
}

// returns true if successful
static bool file_open_txt_write (File *file)
{
    ASSERT (file->data_type > DT_NONE && file->data_type < NUM_DATATYPES ,"invalid data_type=%s (%u)", 
            dt_name (file->data_type), file->data_type);

    // check if type derived from file name is supported, and if not - set it to plain or .gz
    if (file_get_data_type (file->type, false) == DT_NONE) {

        // case: output file has a .gz extension not necessarily related to a file type (eg file.gz rather than file.vcf.gz)
        if ((file_has_ext (file->name, ".gz") || file_has_ext (file->name, ".bgz")) &&  
             file_has_ext (file_exts[txt_out_ft_by_dt[file->data_type][1]], ".gz")) { // data type supports .gz txt output
            
            file->type = txt_out_ft_by_dt[file->data_type][1]; 
            if (flag.bgzf == FLAG_BGZF_BY_ZFILE) flag.bgzf = BGZF_COMP_LEVEL_DEFAULT; // default unless user specified otherwise
        }

        // case: BAM
        else if (file->data_type == DT_BAM) 
            file->type = BAM; 
        
        // case: not .gz and not BAM - use the default plain file format
        else { 
            file->type = txt_out_ft_by_dt[file->data_type][0];  
            ASSINP0 (flag.bgzf == FLAG_BGZF_BY_ZFILE, "using --output in combination with --bgzf, requires the output filename to end with .gz or .bgz");
        }
    }

    // implicit translation flags - eg. "genocat xx.bam.genozip -o xx.fq" implies --fastq
    DataType dt_by_filename = file_get_data_type (file->type, true);
    if ((file->data_type == DT_SAM    && dt_by_filename == DT_BAM  )  ||
        (file->data_type == DT_BAM    && dt_by_filename == DT_SAM  )  ||
        (file->data_type == DT_SAM    && dt_by_filename == DT_FASTQ)  ||
        (file->data_type == DT_BAM    && dt_by_filename == DT_FASTQ)  ||
        (file->data_type == DT_VCF    && dt_by_filename == DT_BCF)    ||
        (file->data_type == DT_FASTA  && dt_by_filename == DT_PHYLIP) ||
        (file->data_type == DT_PHYLIP && dt_by_filename == DT_FASTA)  ||
        (file->data_type == DT_ME23   && dt_by_filename == DT_VCF  ) )
    {
        flag.out_dt = file->data_type = dt_by_filename;
        flags_update_piz_one_file(-1); // update flags accordingly
    }

    if (Z_DT(DT_ME23) && flag.out_dt == DT_VCF)
        ASSINP0 (flag.reference, "--reference must be specified when translating 23andMe to VCF");

    // get the codec    
    file->codec = file_get_codec_by_txt_ft (file->data_type, file->type);
    
    // set to bgzf if the file type implies it
    if (file->codec == CODEC_GZ) {
        file->codec = CODEC_BGZF; // we always write gzip format output with BGZF

        // use default level if source file was not BGZF-compressed or it was of another type (translating)
        if (!z_file->z_flags.bgzf || !flag.reconstruct_as_src)   
            flag.bgzf = BGZF_COMP_LEVEL_DEFAULT; 
    }

    if (file->codec == CODEC_BGZF) {

        // cases we output as plain: 
        // 1. the user overrides with --bgzf=0 (override bgzf set here or before, in flags_update_piz_one_file)
        //    note: a user-specified --bgzf=0 is CODEC_NONE, but a genocat without --bgzf is CODEC_BGZF with level=0 (set in flags_update_piz_one_file)
        // 2. genocat to stdout (terminal or pipe), unless user specified a --bgzf > 0 or its a BAM (binary)
        //    (some BAM downstream tools like GATK expects BGZF)
        if (flag.bgzf == 0 || (flag.bgzf == -1 && flag.to_stdout && flag.out_dt != DT_BAM))
            file->codec = CODEC_NONE;
        
        // set BGZF compression level, in cases we can't or won't compress to the original file's BGZF blocks
        else if (flag.bgzf == FLAG_BGZF_BY_ZFILE) { // user did not use --bgzf to explicitly set the level
            
            // if we're reconstructing as BAM (translated or not) - to stdout (terminal or pipe), use the 
            // fastest BGZF compression if BAM (level) - some downstream tools (eg GATK) expect BAM to always be BGZF-compressed
            if (flag.to_stdout && flag.out_dt == DT_BAM) // genocat of a BAM
                flag.bgzf = 1; 

            // case: translating a BGZF file to another data type (except BAM) - it is not BGZF compressed 
            // (note: if set with --bgzf or implicitly with --output (eg sam.gz) then flag.bgzf>0)
            else if (!flag.reconstruct_as_src && flag.out_dt != DT_BAM)
                file->codec = CODEC_NONE;

            // case: data_modified (inc. non-stdout BAM) we can't recompress using the original BGZF blocks         
            else if (flag.data_modified)             
                flag.bgzf = BGZF_COMP_LEVEL_DEFAULT; 
        }
    }
    
    // don't actually open the output file if we're not going to write to it
    if (flag.no_writer) return true;

    // open the file, based on the codec 
    switch (file->codec) { 
        case CODEC_GZ   :
        case CODEC_BGZF : 
        case CODEC_NONE : file->file = file->redirected ? fdopen (STDOUT_FILENO, "wb") : fopen (file->name, WRITE); break;

        case CODEC_CRAM : file_redirect_output_to_stream (file, "samtools", "view", "-OCRAM", file_samtools_no_PG(), NULL); break;
        
        case CODEC_BCF  : {
            char comp_level[4] = { '-', 'l', '0' + MIN_(flag.bgzf, 9), 0 };
            file_redirect_output_to_stream (file, "bcftools", "view", "-Ob", "/dev/stdin", flag.bgzf >= 0 ? comp_level : NULL); 
            break;
        }

        default         : ABORT ("%s: invalid filename extension for %s files: %s. Please use --output to set another name", global_cmd, dt_name (file->data_type), file->name);
    }

    file_initialize_txt_file_data (file);

    return file->file != 0;
}

static void file_initialize_bufs (File *file)
{
#define INIT(buf) file->buf.name = #buf; \
                  buf_add_to_buffer_list (evb, &file->buf);
    
    INIT (ra_buf);
    INIT (ra_buf_luft);
    INIT (chrom2ref_map);
    INIT (ref2chrom_map);
    INIT (section_list_buf);
    INIT (stats_buf);
    INIT (STATS_buf);
    INIT (bound_txt_names);
    INIT (recon_plan);
    INIT (comp_info);
    INIT (txt_file_info);
    INIT (vb_info[0]);
    INIT (vb_info[1]);
    INIT (line_info[0]);
    INIT (line_info[1]);
    INIT (rejects_report);
    INIT (apriori_tags);
    INIT (unconsumed_txt);
    INIT (bgzf_isizes);
    INIT (coverage);
    INIT (read_count);
    INIT (unmapped_read_count);
}

// we insert all the z_file buffers into the buffer list in advance, to avoid this 
// thread satety issue:
// without our pre-allocation, some of these buffers will be first allocated by a compute threads 
// when the first vb containing a certain did_i is merged in (for the contexts buffers) or
// ra is merged (for ra_buf). while these operations 
// are done while holding a mutex, so that compute threads don't run over each over, buf_alloc 
// may change buf_lists in evb buffers, while the main thread might be doing so concurrently
// resulting in data corruption in evb.buf_list. If evb.buf_list gets corrupted this might result in termination 
// of the execution.
// with these buf_add_to_buffer_list() the buffers will already be in evb's buf_list before any compute thread is run.
static void file_initialize_z_add_to_buf_list (Buffer *buf, const char *func, unsigned line)
{
    buf_add_to_buffer_list (evb, buf);
}

static void file_initialize_z_file_data (File *file)
{
    memset (file->dict_id_to_did_i_map, 0xff, sizeof(file->dict_id_to_did_i_map)); // DID_I_NONE

    for (unsigned i=0; i < MAX_DICTS; i++) 
        ctx_foreach_buffer (&file->contexts[i], true, file_initialize_z_add_to_buf_list);

    file_initialize_bufs (file);
}

static void file_initialize_txt_file_data (File *file)
{
    mutex_initialize (file->recon_plan_mutex[0]);
    mutex_initialize (file->recon_plan_mutex[1]);

    file_initialize_bufs (file);
}

#undef INIT

// returns true if successful
static bool file_open_z (File *file)
{
    ASSINP (!file->is_remote, "it is not possible to access remote genozip files; when attempting to open %s", file->name);
    ASSINP0 (file->name, "it is not possible to redirect genozip files from stdin / stdout");

    // for READ, set data_type
    if (file->mode == READ) {

        ASSINP (!flag.reading_reference || file_has_ext (file->name, REF_GENOZIP_), 
                "You specified file \"%s\", however with --reference or --REFERENCE, you must specify a genozip reference file (%s extension)\n"
                "Tip: You can create a genozip reference file from a FASTA file with 'genozip --make-reference myfasta.fa'",
                file->name, REF_GENOZIP_);

        ASSINP (!flag.reading_chain || file_has_ext (file->name, GENOZIP_EXT), 
                "You specified file \"%s\", however with %s, you must specify a genozip chain file (%s extension)\n"
                "Tip: You can create a genozip chain file from a chain file with eg 'genozip my-chain-file.chain.gz --reference target-coord-ref.ref.genozip'",
                file->name, (command==ZIP) ? "--chain" : "--show-chain", GENOZIP_EXT);

        if (!flag.seg_only || flag_loading_auxiliary) {

            // make sure file is a regular file (not FIFO, directory etc)
            struct stat sb;
            int cause=0, stat_errno=0;
            if (stat (file->name, &sb)) {
                cause = 6; // stat failed
                stat_errno = errno;
            }

            if ((sb.st_mode & S_IFMT) != S_IFREG) cause=7; // not regular file

            if (!cause)
                file->file = fopen (file->name, READ);

            // verify that this is a genozip file 
            // we read the Magic at the end of the file (as the magic at the beginning may be encrypted)
            uint32_t magic;
            if (  cause ||
                  (cause = 1 * !file->file) ||
                  (cause = 2 * !sb.st_size) ||
                  (cause = 3 * !file_seek (file, -(int)sizeof (magic), SEEK_END, true)) || 
                  (cause = 4 * !fread (&magic, sizeof (magic), 1, file->file)) ||
                  (cause = 5 * (BGEN32 (magic) != GENOZIP_MAGIC))) {

                FCLOSE (file->file, file_printname (file));

                if (flag.validate == VLD_REPORT_INVALID) flag.validate = VLD_INVALID_FOUND; 

                const char *cause_str = cause==1 ? strerror (errno)
                                      : cause==2 ? "file is empty"
                                      : cause==3 ? "file_seek failed"
                                      : cause==4 ? "fread failed"
                                      : cause==5 ? "Not a valid genozip file (bad magic)"
                                      : cause==6 ? strerror (stat_errno)
                                      : cause==7 ? "Not a regular file"
                                      :            "no error";

                if (flag.multiple_files) {

                    if (flag.validate == VLD_INVALID_FOUND) // outputs even if --quiet
                        iprintf ("Cannot open %s: %s\n", file_printname (file), cause_str);

                    else if (flag.validate == VLD_NONE) {    // silenced by --quiet
                        static int once=0;
                        WARN ("Skipping %s: %s%s", file_printname (file), cause_str,
                              !(once++) ? " (--quiet to silence this message)" : "");
                    }
                    
                    return true;
                }
                else { // single file
                    if (flag.validate == VLD_REPORT_VALID)
                        exit (EXIT_INVALID_GENOZIP_FILE); // exit quietly - with a return code indicating invalidity
                    else
                        ABORTINP ("Cannot open %s: %s", file_printname (file), cause_str);
                }
            }
            
            // file is valid
            else if (flag.validate == VLD_REPORT_VALID)
                iprintf ("%s\n", file_printname (file)); // print just filename, so a script can use this output
        }

        file->data_type = DT_NONE; // we will get the data type from the genozip header, not by the file name
        file->disk_size_minus_skips = file->disk_size; // will be reduced as we skip sections
    }
    else { // WRITE or WRITEREAD - data_type is already set by file_open
        ASSINP (file_has_ext (file->name, GENOZIP_EXT), 
                "file %s must have a " GENOZIP_EXT " extension", file_printname (file));
        // set file->type according to the data type, overriding the previous setting - i.e. if the user
        // uses the --output option, he is unrestricted in the choice of a file name
        file->type = file_get_z_ft_by_txt_in_ft (file->data_type, txt_file->type); 

        mutex_initialize (file->dicts_mutex);

        if (flag.force) unlink (file->name); // delete file if it already exists (needed in weird cases, eg symlink to non-existing file)

        // if we're writing to a tar file, we get the already-openned tar file
        file->file = tar_is_tar() ? tar_open_file (file->name) 
                                  : fopen (file->name, file->mode);

        if (chain_is_loaded)
            file->z_flags.dual_coords = true;
    }
    
    file_initialize_z_file_data (file);

    return file->file != 0 || flag.seg_only;
}

File *file_open (const char *filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE or WRITEREAD */)
{
    File *file = (File *)CALLOC (sizeof(File));

    file->supertype  = supertype;
    file->is_remote  = filename && url_is_url (filename);
    file->redirected = !filename;
    file->mode       = mode;

    bool is_file_exists = false;

    // is_remote is only possible in READ of a txt_file
    ASSINP ((mode == READ && supertype == TXT_FILE) || !file->is_remote, "expecting file %s to be local, not a URL", filename);

    int64_t url_file_size = 0; // will be -1 if the web/ftp site does not provide the file size
    const char *error = NULL;
    
    if (file->is_remote) {
        error = url_get_status (filename, &is_file_exists, &url_file_size); // accessing is expensive - get existance and size in one call
        if (url_file_size >= 0) file->disk_size = (uint64_t)url_file_size;
    }
    
    else if (!file->redirected) {
        is_file_exists = file_exists (filename);
        error = strerror (errno);

        if (!flag.rejects_coord) {
            if (is_file_exists && mode == READ) 
                file->disk_size = file_get_size (filename);
        }
        else 
            file->disk_size = z_file->rejects_disk_size[flag.rejects_coord-1];
    }

    // return null if genozip input file size is known to be 0, so we can skip it. note: file size of url might be unknown
    if (mode == READ && supertype == TXT_FILE && is_file_exists && !file->disk_size && !url_file_size && 
        !(!file->is_remote && !file->redirected && file_is_fifo (filename))) // a fifo is allowed is be size 0 (as it always is) 
    {
        FREE (file);
        return NULL; 
    }

    if (!file->redirected) {

        ASSINP (mode != READ || is_file_exists, "cannot open \"%s\" for reading: %s", filename, error);
    
        if ((mode == WRITE || mode == WRITEREAD) && 
            is_file_exists && 
            !(!file->is_remote && !file->redirected && file_is_fifo (filename)) && // a fifo can be "overwritten" (that's just normal writing to a fifo)
            !flag.force && 
            !(supertype==TXT_FILE && flag.test) && // not testing piz
            !(supertype==Z_FILE && flag.seg_only)) // not zip with --seg-only
            file_ask_user_to_confirm_overwrite (filename); // function doesn't return if user responds "no"

        // copy filename 
        unsigned fn_size = strlen (filename) + 1; // inc. \0
        file->name = MALLOC (fn_size);
        memcpy (file->name, filename, fn_size);

        if (mode==READ || data_type != DT_NONE) // if its WRITE and DT_NONE, we will not open now, and try again from piz_dispatch after reading the genozip header
            file->type = file_get_type (file->name);

        if (mode == WRITE || mode == WRITEREAD) 
            file->data_type = data_type; // for READ, data_type is set by file_open_*
    }
    else if (mode==READ) {  // stdin

        file->type = stdin_type; 
        file->data_type = file_get_data_type (stdin_type, true);

        Codec codec = file_get_codec_by_txt_ft (file->data_type, file->type);
        if (supertype == TXT_FILE) file->codec = codec;
    }
    else { // stdout
        file->data_type = data_type; 
        file->type = txt_out_ft_by_dt[data_type][flag.bgzf >= 1]; // plain file or .gz file
    }

    file->basename = file_basename (file->name, false, mode==READ ? FILENAME_STDIN : FILENAME_STDOUT, NULL, 0);

    bool success=false;
    switch (supertype) {
        case TXT_FILE: success = (file->mode == READ ? file_open_txt_read : file_open_txt_write) (file); break;
        case Z_FILE:   success = file_open_z (file);   break;
        default:       ABORT ("Error: invalid supertype: %u", supertype);
    }

    ASSINP (success, "cannot open file \"%s\": %s", file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    return file;
}

// index file is it is a disk file of a type that can be indexed
static void file_index_txt (const File *file)
{
    RETURNW (file->name,, "%s: cannot create an index file when output goes to stdout", global_cmd);

    StreamP indexing = NULL;

    switch (file->data_type) {
        case DT_SAM:
        case DT_BAM: 
            RETURNW (file->codec == CODEC_BGZF,, "%s: output file needs to be a .sam.gz or .bam to be indexed", global_cmd); 
            indexing = stream_create (0, 0, 0, 0, 0, 0, 0, "to create an index", "samtools", "index", file->name, NULL); 
            break;
            
        case DT_VCF: 
            RETURNW (file->codec == CODEC_BGZF,, "%s: output file needs to be a .vcf.gz or .bcf to be indexed", global_cmd); 
            RETURNW (vcf_header_get_has_fileformat(),, "%s: file needs to start with ##fileformat=VCF be indexed", global_cmd); 
            indexing = stream_create (0, 0, 0, 0, 0, 0, 0, "to create an index", "bcftools", "index", file->name, NULL); 
            break;

        case DT_FASTQ:
        case DT_FASTA:
            RETURNW (file->codec == CODEC_BGZF || file->codec == CODEC_NONE,, 
                     "%s: To be indexed, the output file cannot be compressed with %s", global_cmd, codec_name (file->codec)); 
            indexing = stream_create (0, 0, 0, 0, 0, 0, 0, "to create an index", "samtools", "faidx", file->name, NULL); 
            break;

        default: break; // we don't know how to create an index for other data types
    }

    if (indexing) {
        progress_new_component (file->basename, "Indexing", false);

        stream_wait_for_exit (indexing);

        progress_finalize_component_time ("Done indexing", DIGEST_NONE);
    }
}

void file_close (File **file_p, 
                 bool index_txt,      // true if we should also index the txt file after it is closed
                 bool cleanup_memory) // true means destroy buffers - use if the closed is NOT near the end of the execution, eg when dealing with unbinding bound files
{
    File *file = *file_p;
    *file_p = NULL;

    if (!file) return; // nothing to do

    // if this is the rejects file, unlink it first to prevent unnecessary flushing when closing
    if (file->supertype == TXT_FILE && flag.rejects_coord) 
        unlink (z_file->rejects_file_name[flag.rejects_coord-1]); // ignore errors (this doesn't work on NTFS)

    if (file->file && file->supertype == TXT_FILE) {

        // finalize a BGZF-compressed reconstructed txt file
        if (file->mode == WRITE && file->codec == CODEC_BGZF)
            bgzf_write_finalize (file);         

        if (file->mode == READ && (file->codec == CODEC_GZ))
            ASSERTW (!gzclose_r((gzFile)file->file), "%s: warning: failed to close file: %s", global_cmd, file_printname (file));

        else if (file->mode == READ && file->codec == CODEC_BZ2)
            BZ2_bzclose((BZFILE *)file->file);
        
        else if (file->mode == READ && file_is_read_via_ext_decompressor (file)) 
            stream_close (&input_decompressor, STREAM_WAIT_FOR_PROCESS);

        else if (file->mode == WRITE && file_is_written_via_ext_compressor (file))
            stream_close (&output_compressor, STREAM_WAIT_FOR_PROCESS);

        // if its stdout - just flush, don't close - we might need it for the next file
        else if (file->mode == WRITE && flag.to_stdout) 
            fflush ((FILE *)file->file);

        else
            FCLOSE (file->file, file_printname (file));
    }

    else if (file->file && file->supertype == Z_FILE) {

        // if there are still DVCF rejects files when close z_file (eg when closing from main_exit()), delete them too
        for (unsigned i=0; i < 2; i++)
            if (file->rejects_file_name[i]) 
                remove (file->rejects_file_name[i]); // ignore errors
    
        if (tar_is_tar() && file->mode != READ)
            tar_close_file (&file->file);
        else
            FCLOSE (file->file, file_printname (file));
    }

    // in case the unlinking didn't work (eg NTFS) - remove the rejects file now that its closed
    if (file->supertype == TXT_FILE && flag.rejects_coord) {
        remove (z_file->rejects_file_name[flag.rejects_coord-1]); // ignore errors
        FREE (z_file->rejects_file_name[flag.rejects_coord-1]);
    }

    // create an index file using samtools, bcftools etc, if applicable
    if (index_txt) file_index_txt (file);

    // free resources if we are NOT near the end of the execution. If we are at the end of the execution
    // it is faster to just let the process die
    if (cleanup_memory) {

        mutex_destroy (file->dicts_mutex);
        mutex_destroy (file->recon_plan_mutex[0]);
        mutex_destroy (file->recon_plan_mutex[1]);
            
        // always destroy all buffers even if unused - for saftey
        for (unsigned i=0; i < MAX_DICTS; i++) // we need to destory all even if unused, because they were initialized in file_initialize_z_file_data
            ctx_destroy_context (&file->contexts[i], i);

        buf_destroy (&file->ra_buf);
        buf_destroy (&file->ra_buf_luft);
        buf_destroy (&file->chrom2ref_map);
        buf_destroy (&file->ref2chrom_map);
        buf_destroy (&file->section_list_buf);
        buf_destroy (&file->unconsumed_txt);
        buf_destroy (&file->bgzf_isizes);
        buf_destroy (&file->coverage);
        buf_destroy (&file->read_count);
        buf_destroy (&file->unmapped_read_count);
        buf_destroy (&file->stats_buf);
        buf_destroy (&file->STATS_buf);
        buf_destroy (&file->bound_txt_names);
        buf_destroy (&file->line_info[0]);
        buf_destroy (&file->line_info[1]);
        buf_destroy (&file->vb_info[0]);
        buf_destroy (&file->vb_info[1]);
        buf_destroy (&file->recon_plan);
        buf_destroy (&file->comp_info);
        buf_destroy (&file->txt_file_info);
        buf_destroy (&file->rejects_report);
        buf_destroy (&file->apriori_tags);

        FREE (file->name);
        FREE (file->basename);
        FREE (file->rejects_file_name[0]);
        FREE (file->rejects_file_name[1]);
        FREE (file);
    }
}

void file_write (File *file, const void *data, unsigned len)
{
    if (!len) return; // nothing to do

    size_t bytes_written = fwrite (data, 1, len, (FILE *)file->file); // use fwrite - let libc manage write buffers for us

    // if we're streaming our genounzip/genocat/genols output to another process and that process has 
    // ended prematurely then exit quietly. In genozip we display an error because this means the resulting
    // .genozip file will be corrupted
    if (!file->name && command != ZIP && errno == EINVAL) exit (EXIT_DOWNSTREAM_LOST);

    // exit quietly if failed to write to stdout - likely downstream consumer (piped executable or terminal) was closed
    if (bytes_written < len && !file->name) exit (EXIT_DOWNSTREAM_LOST);

    // error if failed to write to file
    ASSERT (bytes_written == len, "failed to write %u bytes to %s: %s", len, file->name, strerror(errno));
}

void file_remove (const char *filename, bool fail_quietly)
{
    chmod (filename, S_IRUSR | S_IWUSR); // make sure its +w so we don't get permission denied (ignore errors)

    int ret = remove (filename); 
    ASSERTW (!ret || fail_quietly, "Warning: failed to remove %s: %s", filename, strerror (errno));
}

void file_mkfifo (const char *filename)
{
#ifndef _WIN32
    file_remove (filename, true);
    ASSERT (!mkfifo (filename, 0666), "mkfifo failed for %s: %s", filename, strerror (errno));

#else
    ABORT0 ("file_mkfifo not supported on Windows");
#endif
}

bool file_is_fifo (const char *filename)
{
    if (flag.is_windows) return false; // we don't support FIFOs in Win32 yet

    struct stat st;
    ASSERT (!stat (filename, &st), "stat failed on %s", filename);

    return S_ISFIFO (st.st_mode);
}

bool file_exists (const char *filename)
{
    if (!filename || !filename[0]) return false;

    return !access (filename, F_OK);
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
        basename = (char *)MALLOC (len + 1); // +1 for \0
    else
        len = MIN_(len, basename_size-1);

    sprintf (basename, "%.*s", (int)len, start);

    return basename;
}

// returns true if successful. depending on soft_fail, a failure will either emit an error 
// (and exit) or a warning (and return).
bool file_seek (File *file, int64_t offset, 
                int whence, // SEEK_SET, SEEK_CUR or SEEK_END
                int soft_fail) // 1=warning 2=silent
{
    ASSERTNOTNULL (file);
    ASSERTNOTNULL (file->file);
    
    // in SEEK_SET of a z_file that is being tarred, update the offset to the beginning of the file data in the tar file
    if (file->supertype == Z_FILE && whence == SEEK_SET) {
        offset += tar_file_offset(); // 0 if not using tar

        int64_t old_offset = ftello64 ((FILE *)file->file);
        if (old_offset == offset) return true; // already at the right offset
    }

    int ret = fseeko64 ((FILE *)file->file, offset, whence);

    if (soft_fail) {
        if (!flag.to_stdout && soft_fail==1) {
            ASSERTW (!ret, errno == EINVAL ? "Warning: Error while reading file %s (fseeko64 (whence=%d offset=%"PRId64")): it is too small%s" 
                                           : "Warning: fseeko failed on file %s (whence=%d offset=%"PRId64"): %s", 
                     file_printname (file), whence, offset, errno == EINVAL ? "" : strerror (errno));
        }
    } 
    else
        ASSERT (!ret, "fseeko(offset=%"PRId64" whence=%d) failed on file %s (FILE*=%p remote=%s redirected=%s): %s", 
                       offset, whence, file_printname (file), file->file, TF(file->is_remote), TF(file->redirected), strerror (errno));

    return !ret;
}

int64_t file_tell_do (File *file, bool soft_fail, const char *func, unsigned line)
{
    if (command == ZIP && file->supertype == TXT_FILE && file->codec == CODEC_GZ)
        return gzconsumed64 ((gzFile)file->file); 
    
    if (command == ZIP && file->supertype == TXT_FILE && file->codec == CODEC_BZ2)
        return BZ2_consumed ((BZFILE *)file->file); 
    int64_t offset = ftello64 ((FILE *)file->file);
    ASSERT (offset >= 0 || soft_fail, "called from %s:%u: ftello64 failed for %s (FILE*=%p remote=%s redirected=%s): %s", 
            func, line, file->name, file->file, TF(file->is_remote), TF(file->redirected), strerror (errno));

    if (offset < 0) return -1; // soft fail

    // in in z_file that is being tarred, update the offset to the beginning of the file data in the tar file
    if (file->supertype == Z_FILE)
        offset -= tar_file_offset(); // 0 if not using tar

    return offset;    
}

uint64_t file_get_size (const char *filename)
{
    struct stat64 st;
    
    int ret = stat64(filename, &st);
    ASSERT (!ret, "stat64 failed on '%s': %s", filename, strerror(errno));
    
    return st.st_size;
}

bool file_is_dir (const char *filename)
{
    struct stat64 st;
    int ret = stat64(filename, &st); // 0 if successful
    
    return !ret && S_ISDIR (st.st_mode);
}

// reads an entire file into a buffer. if filename is "-", reads from stdin
void file_get_file (VBlockP vb, const char *filename, Buffer *buf, const char *buf_name,
                    bool add_string_terminator)
{
    #define MAX_STDIN_DATA_SIZE 10000000
    bool is_stdin = !strcmp (filename, "-");

    uint64_t size = is_stdin ? MAX_STDIN_DATA_SIZE : file_get_size (filename);

    buf_alloc (vb, buf, 0, size + add_string_terminator, char, 1, buf_name);

    FILE *file = is_stdin ? stdin : fopen (filename, "rb");
    ASSINP (file, "cannot open \"%s\": %s", filename, strerror (errno));

    buf->len = fread (buf->data, 1, size, file);
    ASSERT (is_stdin || buf->len == size, "Error reading file %s: %s", filename, strerror (errno));
    ASSERT (!is_stdin || buf->len < MAX_STDIN_DATA_SIZE, "Error reading from stdin: too much data for stdin (beyond maximum of %u bytes). Try placing the data in a file instead.", MAX_STDIN_DATA_SIZE-1);

    if (add_string_terminator) 
        *AFTERENT (char, *buf) = 0;

    FCLOSE (file, filename);
}

// writes data to a file and flushes it, returns true if successful

static Mutex put_data_mutex = {};
#define MAX_PUT_DATA_FILES_PER_EXECUTION 10 // maximum files deletable at abort
static const char *put_data_tmp_filenames[MAX_PUT_DATA_FILES_PER_EXECUTION] = {};
static unsigned num_put_files=0;

bool file_put_data (const char *filename, const void *data, uint64_t len, 
                    mode_t mode) // optional - ignored if 0
{
    char *tmp_filename = MALLOC (strlen(filename)+5);
    // we first write to tmp_filename, and after we complete and flush, we rename to the final name
    // this is important, eg for the reference cache files - if a file exists (in its final name) - then it is fully written
    sprintf (tmp_filename, "%s.tmp", filename);

    file_remove (filename, true);
    file_remove (tmp_filename, true);

    FILE *file = fopen (tmp_filename, "wb");
    if (!file) return false;

    // save file name in put_data_tmp_filenames, to be deleted in case of aborting by file_put_data_abort
    mutex_initialize (put_data_mutex); // first caller initializes. not thread safe, but good enough for the purpose.
    mutex_lock (put_data_mutex);
    
    unsigned my_file_i = num_put_files;

    if (num_put_files < MAX_PUT_DATA_FILES_PER_EXECUTION) 
        put_data_tmp_filenames[num_put_files++] = tmp_filename; 

    mutex_unlock (put_data_mutex);

    size_t written = fwrite (data, 1, len, file);
    
    SAVE_VALUE (errno);
    fflush (file);    
    FCLOSE (file, tmp_filename); 
    RESTORE_VALUE (errno); // in cases caller wants to print fwrite error

    if (written != len) {
        put_data_tmp_filenames[my_file_i] = NULL; // no need to lock mutex
        remove (tmp_filename);
        return false;
    }

    // we can't enter if file_put_data_abort is active, or it needs to wait for us before deleting tmp files
    mutex_lock (put_data_mutex);

    int renamed_failed = rename (tmp_filename, filename);

    put_data_tmp_filenames[my_file_i] = NULL; // remove tmp file name from list 
    
    mutex_unlock (put_data_mutex);

    ASSERT (!renamed_failed, "failed to rename %s to %s: %s", tmp_filename, filename, strerror (errno));
    FREE (tmp_filename);

    if (mode) 
        ASSERT (!chmod (filename, mode), "failed to chmod %s: %s", filename, strerror (errno));

    return true;
}

// error handling: unlink files currently being written (the actual writing will terminate when the thread terminates)
void file_put_data_abort (void)
{
    if (!put_data_mutex.initialized) return;

    mutex_lock (put_data_mutex);

    for (unsigned i=0; i < num_put_files; i++) 
        if (put_data_tmp_filenames[i]) {
            unlink (put_data_tmp_filenames[i]); // ignore errors
            remove (put_data_tmp_filenames[i]); // in case unlinked failed - eg NTFS - ignore errors
        }

    // mutex remains locked - no more files can be put after this point
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

// PIZ: guess original filename from uncompressed txt filename and compression algoritm (allocated memory)
const char *file_guess_original_filename (const File *file)
{
    if (file->codec == CODEC_NONE) return file->name;

    unsigned len = strlen (file->name) + 10;
    char *org_name = MALLOC (len);
    strcpy (org_name, file->name);

    const char *ext = codec_args[file->codec].ext;

    // remove existing extension if needed (eg when replacing .sam with .bam)
    if (ext[0] == '-') {
        char *last_dot = strrchr (org_name, '.');
        if (last_dot) *last_dot = 0;
    }

    int org_len = strlen (org_name); 
    if (org_len >= 4 && !strcmp (&org_name[org_len-4], ".bam")) {
        // don't add extension to .bam
    }
    else if (!file_has_ext (org_name, &ext[1])) // add new extension
        strcpy (&org_name[org_len], &ext[1]);

    return org_name;
}

const char *file_plain_ext_by_dt (DataType dt)
{
    FileType plain_ft = txt_in_ft_by_dt[dt][0].in;

    return file_exts[plain_ft];
}

// eg "file.fastq.gz" -> "file.fastq"
void file_remove_codec_ext (char *filename, FileType ft)
{
    const char *codec_ext;
    if (!((codec_ext = strchr (&file_exts[ft][1], '.')))) return; // this file type's extension doesn't have a codec ext (it is, eg ".fastq" not ".fastq.gz")
    
    unsigned codec_ext_len = strlen (codec_ext);
    unsigned fn_len = strlen (filename);

    // make sure filename actually has the codec (eg its not eg "(stdin)")
    if (fn_len > codec_ext_len && !strcmp (&filename[fn_len-codec_ext_len], codec_ext))
        filename[fn_len-codec_ext_len] = 0; // shorten string
}

// change eg "c:\dir\file" to "/c/dir/file" and add full path (allocates memory - caller should free)
char *file_make_unix_filename (char *filename)
{
    char path[PATH_MAX];
    unsigned len = strlen (filename);
    
    if (
#ifdef _WIN32
        (len >= 2 && filename[1] == ':') || // full name - starting with eg C:
#endif
        (len >= 1 && filename[0] == '/') || // full name - starting with /
        !getcwd (path, sizeof (path))) // path too long
        path[0] = 0;  // don't store path

    char *full_fn = malloc (strlen (z_name) + strlen (path) + 2);
    sprintf (full_fn, "%s%s%s", path, *path ? "/" : "", z_name);

#ifdef _WIN32 // convert to Unix-style filename
    len = strlen (full_fn);

    if (len >= 2 && full_fn[1] == ':') {
        full_fn[1] = full_fn[0];
        full_fn[0] = '/';
    }

    for (unsigned i=0; i < len; i++)
        if (full_fn[i] == '\\') full_fn[i] = '/';
#endif

    return full_fn;
}

