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
#include "codec.h"
#include "mutex.h"
#include "bgzf.h"

// globals
File *z_file   = NULL;
File *txt_file = NULL;

static StreamP input_decompressor  = NULL; // bcftools or xz, unzip or samtools - only one at a time
static StreamP output_compressor   = NULL; // bgzip, samtools, bcftools

static FileType stdin_type = UNKNOWN_FILE_TYPE; // set by the --input command line option

// global pointers - so the can be compared eg "if (mode == READ)"
const char *READ  = "rb";  // use binary mode (b) in read and write so Windows doesn't add \r
const char *WRITE = "wb";
const char *WRITEREAD = "wb+"; // only supported for z_file

const char *file_exts[] = FILE_EXTS;

static const struct { FileType in; Codec codec; FileType out; } txt_in_ft_by_dt[NUM_DATATYPES][30] = TXT_IN_FT_BY_DT;
static const FileType txt_out_ft_by_dt[NUM_DATATYPES][20] = TXT_OUT_FT_BY_DT;
static const FileType z_ft_by_dt[NUM_DATATYPES][20] = Z_FT_BY_DT;

// get data type by file type
DataType file_get_data_type (FileType ft, bool is_input)
{
    // note: if make-reference, we scan the array from dt=0 (DT_REF), otherwise we ignore DT_REF
    for (DataType dt=!flag.make_reference; dt < NUM_DATATYPES; dt++) 
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

// get codec by txt file type
static void file_get_codecs_by_txt_ft (DataType dt, FileType txt_ft, FileMode mode, Codec *codec)
{
    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].in == txt_ft) {
            if (codec) *codec         = txt_in_ft_by_dt[dt][i].codec;
            return;
        }

    if (*codec) *codec = (mode == WRITE ? CODEC_NONE : CODEC_UNKNOWN);
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
        *raw_name = MALLOC (len + 30);
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

        flag.stdin_size = flag.stdin_size * 10 + (size_str[i] - '0'); 
    }
}

void file_set_input_type (const char *type_str)
{
    char ext[strlen (type_str) + 2]; // +1 for . +1 for \0
    sprintf (ext, ".%s", type_str);

    str_tolower (ext, ext); // lower-case to allow case-insensitive --input argument (eg vcf or VCF)

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

    fprintf (stderr, "\n");
}

static void file_redirect_output_to_stream (File *file, char *exec_name, 
                                            const char *stdout_option, const char *format_option_1, const char *format_option_2)
{
    char threads_str[20];
    sprintf (threads_str, "%u", global_max_threads);

    FILE *redirected_stdout_file = NULL;
    if (!flag.to_stdout) {
        redirected_stdout_file = fopen (file->name, file->mode); // exec_name will redirect its output to this file
        ASSINP (redirected_stdout_file, "%s: cannot open file %s: %s", global_cmd, file->name, strerror(errno));
    }
    char reason[100];
    sprintf (reason, "To output a %s file", file_exts[file->type]);
    output_compressor = stream_create (0, 0, 0, global_max_memory_per_vb, 
                                       redirected_stdout_file, // output is redirected unless flag.to_stdout
                                       0, reason,
                                       exec_name, 
                                       stdout_option, // either to the terminal or redirected to output file
                                       "--threads", threads_str,
                                       format_option_1, format_option_2,
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
        StreamP samtools = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, "To read/write CRAM files",
                                          "samtools", "view", "--junk", NULL);
        usleep (50000 * i); // wait for samtools

        // read both stderr and stdout from samtools
        len  = read (fileno (stream_from_stream_stderr (samtools)), samtools_help_text,       SAMTOOLS_HELP_MAX_LEN-1);
        len += read (fileno (stream_from_stream_stdout (samtools)), &samtools_help_text[len], SAMTOOLS_HELP_MAX_LEN-len-1);

        stream_close (&samtools, STREAM_DONT_WAIT_FOR_PROCESS);
    } 
    samtools_help_text[len] = '\0'; // terminate string (more portable, strnstr and memmem are non-standard)

    ASSERT0 (len >= MIN_ACCEPTABLE_LEN, "Error in file_samtools_no_PG: no response from \"samtools view --junk\"");

    return ret_str[(has_no_PG = !!strstr (samtools_help_text, "--no-PG"))];
}

// show meaningful error if file is not a supported type and return TRUE if it file should be skipped
static bool file_open_txt_read_test_valid_dt (File *file)
{ 
    if (file->data_type == DT_NONE) { 

        if (flag.multiple_files) {
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

    ASSINP0 (!flag.make_reference || file->data_type == DT_REF, "Error: --make-reference can only be used with FASTA files");

    return false; // all good - no need to skip this file
}

static bool file_open_txt_read (File *file)
{
    // if user provided the type with --input, we use that overriding the type derived from the file name
    if (stdin_type) file->type = stdin_type;

    file->data_type = file_get_data_type (file->type, true);

    // show meaningful error if file is not a supported data type
    if (file_open_txt_read_test_valid_dt (file)) return true; // skip this file

    // open the file, based on the codec
    file_get_codecs_by_txt_ft (file->data_type, file->type, READ, &file->codec);
    
    switch (file->codec) { 
        case CODEC_NONE:
            file->file = file->is_remote ? url_open (NULL, file->name) : fopen (file->name, READ);
            break;

        case CODEC_GZ:
        case CODEC_BGZF: {
            file->file = (FILE *)(file->is_remote ? url_open (NULL, file->name) : fopen (file->name, "rb"));
            ASSERT (file->file, "Error in file_open_txt_read: failed to open %s: %s", file->name, strerror (errno));

            // read the first potential BGZF block to test if this is GZ or BGZF
            uint8_t block[BGZF_MAX_BLOCK_SIZE]; 
            uint32_t block_size;

            int32_t bgzf_uncompressed_size = bgzf_read_block (file, block, &block_size, true);
            
            // case: this is indeed a bgzf - we put the still-compressed data in vb->compressed for later consumption
            // is txtfile_read_block_bgzf
            if (bgzf_uncompressed_size >= 0) {
                file->codec = CODEC_BGZF;
                
                buf_alloc (evb, &evb->compressed, block_size, 1, "compressed"); 
                evb->compressed.param = bgzf_uncompressed_size; // pass uncompressed size in param
                buf_add (&evb->compressed, block, block_size);
                
            }

            // case: this is a non-BGZF gzip format - open with zlib and hack back the read bytes 
            // (note: we cannot re-read the bytes from the file as the file might be piped in)
            else if (bgzf_uncompressed_size == BGZF_BLOCK_GZIP_NOT_BGZIP) {
                file->codec = CODEC_GZ;
                file->file  = gzdopen (fileno((FILE *)file->file), READ); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
                gzinject (file->file, block, block_size); // a hack - adds a 18 bytes of compressed data to the in stream, which will be consumed next, instead of reading from disk
            } 

            // case: this is not GZIP format at all. treat as a plain file, and put the data read in vb->compressed 
            // for later consumption is txtfile_read_block_plain
            else if (bgzf_uncompressed_size == BGZF_BLOCK_IS_NOT_GZIP) {
                file->codec = CODEC_NONE;
                buf_alloc (evb, &evb->compressed, block_size, 1, "compressed"); 
                buf_add (&evb->compressed, block, block_size);
            }

            break;
        }
        case CODEC_BZ2:
            if (file->is_remote) {            
                FILE *url_fp = url_open (NULL, file->name);
                file->file = BZ2_bzdopen (fileno(url_fp), READ); // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
            }
            else
                file->file = BZ2_bzopen (file->name, READ);  // for local files we decompress ourselves   
            break;

        case CODEC_XZ:
            input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,     // url
                                                "To uncompress an .xz file", "xz",       // reason, exec_name
                                                file->is_remote ? SKIP_ARG : file->name, // local file name 
                                                "--threads=8", "--decompress", "--keep", "--stdout", 
                                                flag.quiet ? "--quiet" : SKIP_ARG, 
                                                NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case CODEC_ZIP:
            input_decompressor = stream_create (0, global_max_memory_per_vb, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,     // url
                                                "To uncompress a .zip file", "unzip",    // reason, exec_name
                                                "-p", // must be before file name
                                                file->is_remote ? SKIP_ARG : file->name, // local file name 
                                                flag.quiet ? "--quiet" : SKIP_ARG, 
                                                NULL);            
            file->file = stream_from_stream_stdout (input_decompressor);
            break;

        case CODEC_BCF:
        case CODEC_BAM: 
        case CODEC_CRAM: {
            bool bam =  (file->codec == CODEC_BAM);
            bool cram = (file->codec == CODEC_CRAM);

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
                                                (bam || cram) ? (file_samtools_no_PG() ? "--no-PG" : "-h") : NULL,  // don't add a PG line to the header (just repeat -h if this is an older samtools without --no-PG - no harm)
                                                cram ? ref_get_cram_ref() : NULL,
                                                NULL);
            file->file = stream_from_stream_stdout (input_decompressor);
            break;
        }

        default:
            ABORT ("%s: invalid filename extension for %s files: %s", global_cmd, dt_name (file->data_type), file->name);
    }

    file->txt_data_size_single = file->disk_size; 

    return file->file != 0;
}

// returns true if successful
static bool file_open_txt_write (File *file)
{
    ASSERT (file->data_type > DT_NONE && file->data_type < NUM_DATATYPES ,"Error in file_open_txt_write: invalid data_type=%s (%u)", 
            dt_name (file->data_type), file->data_type);

    // check if type derived from file name is supported, and if not - set it to plain or .gz
    if (file_get_data_type (file->type, false) == DT_NONE) {

        // case: output file has a .gz extension (eg the user did "genounzip xx.vcf.genozip -o yy.gz")
        if ((file_has_ext (file->name, ".gz") || file_has_ext (file->name, ".bgz")) &&  
            file_has_ext (file_exts[txt_out_ft_by_dt[file->data_type][1]], ".gz")) { // data type supports .gz txt output
            
            file->type = txt_out_ft_by_dt[file->data_type][1]; 
            flag.bgzf = true;
        }
        
        // case: not .gz - use the default plain file format
        else { 
            file->type = txt_out_ft_by_dt[file->data_type][0]; 
            ASSINP (!flag.bgzf, "%s: using --output in combination with --bgzip, requires the output filename to end with .gz or .bgz", global_cmd);
        }
    }

    // get the codec    
    file_get_codecs_by_txt_ft (file->data_type, file->type, WRITE, &file->codec);

    // don't actually open the output file if we're just testing in genounzip or PIZing a reference file
    if (flag.test || flag.reading_reference) return true;

    // open the file, based on the codec 
    switch (file->codec) { 
        case CODEC_GZ   :
        case CODEC_BGZF : 
        case CODEC_NONE : file->file = fopen (file->name, WRITE); break;
        //case CODEC_GZ   :
        //case CODEC_BGZF : file_redirect_output_to_stream (file, "bgzip", "--stdout", NULL, NULL); break;
        case CODEC_BCF  : file_redirect_output_to_stream (file, "bcftools", "view", "-Ob", NULL); break;
        case CODEC_BAM  : file_redirect_output_to_stream (file, "samtools", "view", "-OBAM",  file_samtools_no_PG()); break;
        case CODEC_CRAM : file_redirect_output_to_stream (file, "samtools", "view", "-OCRAM", file_samtools_no_PG()); break;

        default         : ABORT ("%s: invalid filename extension for %s files: %s", global_cmd, dt_name (file->data_type), file->name);
    }

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
        INIT (nodes);
        INIT (node_i);
        INIT (global_hash);
        INIT (ol_dict);
        INIT (ol_nodes);
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
    INIT (stats_buf_1);
    INIT (stats_buf_2);
}
#undef INIT

// returns true if successful
static bool file_open_z (File *file)
{
    // for READ, set data_type
    if (file->mode == READ) {

        if (!file_has_ext (file->name, GENOZIP_EXT)) {
            if (flag.multiple_files) 
                RETURNW (false, true, "Skipping %s - it doesn't have a .genozip extension", file_printname (file))
            else {
                if (flag.reading_reference)
                    ABORT ("%s: with --reference or --REFERENCE, you must specify a genozip reference file (.ref.genozip extension)\n"
                           "Tip: You can create a genozip reference file from a FASTA file with 'genozip --make-reference myfasta.fa'",
                           global_cmd)
                else
                    ABORT ("%s: file %s must have a " GENOZIP_EXT " extension", global_cmd, file_printname (file));
            }
        }

        file->data_type = file_get_dt_by_z_ft (file->type); // if we will verify this is correct in zfile_read_genozip_header
    }
    else { // WRITE or WRITEREAD - data_type is already set by file_open
        ASSINP (file_has_ext (file->name, GENOZIP_EXT), "%s: file %s must have a " GENOZIP_EXT " extension", 
                              global_cmd, file_printname (file));
        // set file->type according to the data type, overriding the previous setting - i.e. if the user
        // uses the --output option, he is unrestricted in the choice of a file name
        file->type = file_get_z_ft_by_txt_in_ft (file->data_type, txt_file->type); 
    }

    ASSERT (!file->is_remote, "Error: it is not possible to access remote genozip files; when attempting to open %s", file->name);
    
    if (!flag.test_seg)
        file->file = fopen (file->name, file->mode);

    file_initialize_z_file_data (file);

    return file->file != 0 || flag.test_seg;
}

File *file_open (const char *filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE or WRITEREAD */)
{
    ASSINP0 (filename, "Error in file_open: filename is null");

    //File *file = (File *)CALLOC (sizeof(File) + (((mode == READ || mode == WRITEREAD) && supertype == Z_FILE) ? READ_BUFFER_SIZE : 0));
    File *file = (File *)CALLOC (sizeof(File));

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

    if ((mode == WRITE || mode == WRITEREAD) && file_exists && !flag.force && !(supertype==TXT_FILE && flag.test))
        file_ask_user_to_confirm_overwrite (filename); // function doesn't return if user responds "no"

    // copy filename 
    unsigned fn_size = strlen (filename) + 1; // inc. \0
    file->name = MALLOC (fn_size);
    memcpy (file->name, filename, fn_size);

    file->mode = mode;

    if (mode==READ || data_type != DT_NONE) // if its NONE, we will not open now, and try again from piz_dispatch after reading the genozip header
        file->type = file_get_type (file->name, true);

    if (file->mode == WRITE || file->mode == WRITEREAD) 
        file->data_type = data_type; // for READ, data_type is set by file_open_*

    bool success=false;
    switch (supertype) {
        case TXT_FILE: success = (file->mode == READ ? file_open_txt_read : file_open_txt_write) (file); break;
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

    //File *file = (File *)CALLOC (sizeof(File) + ((mode == READ && supertype == Z_FILE) ? READ_BUFFER_SIZE : 0));
    File *file = (File *)CALLOC (sizeof(File));

    file->file = (mode == READ) ? fdopen (STDIN_FILENO,  "rb")
                                : fdopen (STDOUT_FILENO, "wb");
    ASSINP (file->file, "%s: Failed to redirect %s: %s", global_cmd, (mode==READ ? "stdin" : "stdout"), strerror (errno));

    file->supertype = supertype;
    
    if (mode==READ) {
        file->data_type = file_get_data_type (stdin_type, true);
        file->type  = stdin_type; 
        file->codec = CODEC_NONE;
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

        if (file->mode == WRITE && file->supertype == TXT_FILE && file->codec == CODEC_BGZF)
            ASSERT (fwrite (BGZF_EOF, BGZF_EOF_LEN, 1, (FILE *)file->file), "Failed to write BGZF EOF to %s: %s", file->name, strerror (errno)); 
            // proceed to FCLOSE below     

        if (file->mode == READ && (file->codec == CODEC_GZ))
            ASSERTW (!gzclose_r((gzFile)file->file), "%s: warning: failed to close file: %s", global_cmd, file_printname (file))

        else if (file->mode == READ && file->codec == CODEC_BZ2)
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
            ctx_destroy_context (&file->contexts[i]);

        buf_destroy (&file->dict_data);
        buf_destroy (&file->ra_buf);
        buf_destroy (&file->ra_min_max_by_chrom);
        buf_destroy (&file->chroms_sorted_index);
        buf_destroy (&file->alt_chrom_map);
        buf_destroy (&file->section_list_buf);
        buf_destroy (&file->section_list_dict_buf);
        buf_destroy (&file->unconsumed_txt);
        buf_destroy (&file->bgzf_isizes);
        buf_destroy (&file->stats_buf_1);
        buf_destroy (&file->stats_buf_2);

        FREE (file->name);
        FREE (file->basename);
        FREE (file);
    }
}

void file_write (File *file, const void *data, unsigned len)
{
    if (!len) return; // nothing to do

    size_t bytes_written = fwrite (data, 1, len, (FILE *)file->file); // use fwrite - let libc manage write buffers for us

    // if we're streaming our genounzip/genocat/genols output to another process and that process has 
    // ended prematurely then exit quietly. in genozip we display an error because this means the resulting
    // .genozip file will be corrupted
    if (!file->name && command != ZIP && errno == EINVAL) exit(0);

    // exit quietly if failed to write to stdout - likely downstream consumer (piped executable or terminal) was closed
    if (bytes_written < len && !file->name) exit (0);

    // error if failed to write to file
    ASSERT (bytes_written == len, "Error: failed to write %u bytes to %s: %s", len, file->name, strerror(errno));
}

void file_remove (const char *filename, bool fail_quietly)
{
    int ret = remove (filename); 
    ASSERTW (!ret || fail_quietly, "Warning: failed to remove %s: %s", filename, strerror (errno));
}

void file_mkfifo (const char *filename)
{
#ifndef _WIN32
    file_remove (filename, true);
    ASSERT (!mkfifo (filename, 0666), "Error in file_mkfifo: mkfifo failed for %s: %s", filename, strerror (errno));

#else
    ABORT0 ("file_mkfifo not supported on Windows");
#endif
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
        len = MIN (len, basename_size-1);

    sprintf (basename, "%.*s", (int)len, start);

    return basename;
}

// returns true if successful. depending on soft_fail, a failure will either emit an error 
// (and exit) or a warning (and return).
bool file_seek (File *file, int64_t offset, 
                int whence, // SEEK_SET, SEEK_CUR or SEEK_END
                int soft_fail) // 1=warning 2=silent
{
#ifdef __APPLE__
    int ret = fseeko ((FILE *)file->file, offset, whence);
#else
    int ret = fseeko64 ((FILE *)file->file, offset, whence);
#endif

    if (soft_fail) {
        if (!flag.to_stdout && soft_fail==1) {
            ASSERTW (!ret, errno == EINVAL ? "Warning: Error while reading file %s: it is too small%s" 
                                           : "Warning: fseeko failed on file %s: %s", 
                    file_printname (file),  errno == EINVAL ? "" : strerror (errno));
        }
    } 
    else {
        ASSERT (!ret, "Error: fseeko failed on file %s: %s", file_printname (file), strerror (errno));
    }

    return !ret;
}

uint64_t file_tell (File *file)
{
    if (command == ZIP && file == txt_file && file->codec == CODEC_GZ)
        return gzconsumed64 ((gzFile)txt_file->file); 
    
    if (command == ZIP && file == txt_file && file->codec == CODEC_BZ2)
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

void file_get_file (VBlockP vb, const char *filename, Buffer *buf, const char *buf_name,
                    bool add_string_terminator)
{
    uint64_t size = file_get_size (filename);

    buf_alloc (vb, buf, size + add_string_terminator, 1, buf_name);

    FILE *file = fopen (filename, "rb");
    ASSINP (file, "Error: cannot open %s: %s", filename, strerror (errno));

    size_t bytes_read = fread (buf->data, 1, size, file);
    ASSERT (bytes_read == (size_t)size, "Error reading file %s: %s", filename, strerror (errno));

    buf->len = size;

    if (add_string_terminator) buf->data[size] = 0;

    FCLOSE (file, filename);
}

// writes data to a file, return true if successful
bool file_put_data (const char *filename, void *data, uint64_t len)
{
    FILE *file = fopen (filename, "wb");
    if (!file) return false;
    
    size_t written = fwrite (data, 1, len, file);
    
    SAVE_VALUE (errno);
    FCLOSE (file, filename);
    RESTORE_VALUE (errno); // in cases caller wants to print fwrite error

    return written == len;
}

// writes a buffer to a file, return true if successful
bool file_put_buffer (const char *filename, const Buffer *buf, unsigned buf_word_width)
{
    return file_put_data (filename, buf->data, buf->len * buf_word_width);
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
    else // add new extension
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

    ASSERT (fn_len > codec_ext_len, "Error in file_remove_codec_ext: expecting fn_len=%u > codec_ext_len=%u", fn_len, codec_ext_len)

    filename[fn_len-codec_ext_len] = 0; // shorten string
}

