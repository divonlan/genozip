// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#endif
#define Z_LARGE64
#ifdef __APPLE__
    #define off64_t __int64_t
#endif
#include "bzlib/bzlib.h"
#include "zlib/zlib.h"
#include "file.h"
#include "url.h"
#include "codec.h"
#include "mgzip.h"
#include "progress.h"
#include "tar.h"
#include "writer.h"
#include "filename.h"
#include "huffman.h"

// globals
FileP z_file   = NULL;
FileP txt_file = NULL;

static StreamP input_decompressor = NULL; // bcftools, xz, unzip, samtools or orad - only one at a time
static StreamP output_compressor  = NULL; // samtools (for cram), bcftools

// global pointers - so the can be compared eg "if (mode == READ)"
rom READ      = "rb";  // use binary mode (b) in read and write so Windows doesn't add \r. "b" is accepted but ignored in Linux and MacOS.
rom WRITE     = "wb";
rom WRITEREAD = "wb+"; // only supported for z_file and gencomp disk files

rom file_exts[] = FILE_EXTS;

static const struct { FileType in; Codec codec; FileType out; } txt_in_ft_by_dt[NUM_DATATYPES][50] = TXT_IN_FT_BY_DT;
static const FileType txt_out_ft_by_dt[NUM_DATATYPES][20] = TXT_OUT_FT_BY_DT;
static const FileType z_ft_by_dt[NUM_DATATYPES][20] = Z_FT_BY_DT;

// get data type by file type
DataType file_get_data_type_of_input_file (FileType ft)
{
    // note: if make-reference, we scan the array from dt=0 (DT_REF), otherwise we ignore DT_REF
    for (DataType dt=!flag.make_reference; dt < NUM_DATATYPES; dt++) 
        for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
            if (txt_in_ft_by_dt[dt][i].in == ft) 
                return dt;

    return (ft == UNKNOWN_FILE_TYPE) ? DT_GNRIC : DT_NONE;
}

DataType file_piz_get_dt_of_out_filename (void)
{
    if (!flag.out_filename) return DT_NONE;
    
    FileType ft = file_get_type (flag.out_filename);

    for (DataType dt=1/*skip DT_REF*/; dt < NUM_DATATYPES; dt++) 
        for (unsigned i=0; txt_out_ft_by_dt[dt][i]; i++)
            if (txt_out_ft_by_dt[dt][i] == ft) 
                return dt;

    return DT_NONE;
}

// get genozip file type by txt file type
FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft)
{
    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].in == txt_ft) return txt_in_ft_by_dt[dt][i].out;

    return UNKNOWN_FILE_TYPE;
}

// get codec by txt file type
Codec file_get_codec_by_txt_ft (DataType dt, FileType txt_ft, bool source)
{
    if (source && txt_ft == BAM) return CODEC_BAM; // if !source, it would be CODEC_BGZF
    if (txt_ft == BCF || txt_ft == BCF_GZ || txt_ft == BCF_BGZF) return CODEC_BCF; 

    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].in == txt_ft) 
            return txt_in_ft_by_dt[dt][i].codec;

    return CODEC_NONE;
}

// get codec by txt file type
FileType file_get_txt_ft_by_codec (DataType dt, Codec codec)
{
    for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
        if (txt_in_ft_by_dt[dt][i].codec == codec) 
            return txt_in_ft_by_dt[dt][i].in;

    return UNKNOWN_FILE_TYPE;
}

// get data_type of e.g. "myfile.fastq.genozip" 
// note: DT_SAM is returned for SAM/BAM/CRAM, and DT_VCF for VCF/BCF
DataType file_get_dt_by_z_filename (rom z_filename)
{
    FileType z_ft = file_get_type (z_filename);

    for (DataType dt=0; dt < NUM_DATATYPES; dt++)
        for (unsigned i=0; z_ft_by_dt[dt][i]; i++)
            if (z_ft == z_ft_by_dt[dt][i]) 
                return dt;
    
    return DT_NONE;
}

// returns default file type of the .genozip file of the given dt (e.g. DT_FASTA -> .fasta.genozip) 
FileType file_get_default_z_ft_of_data_type (DataType dt)
{
    return z_ft_by_dt[dt][0];
}

// possible arguments for --input
StrTextLong file_compressible_extensions (bool plain_only)
{
    StrTextLong s;
    int s_len = 0;

    for (DataType dt=1; dt < NUM_DATATYPES; dt++) { // start from 1, excluding DT_REFERENCE
        
        if (dt == DT_GNRIC || dt == DT_ME23 || !txt_in_ft_by_dt[dt][0].in) continue;

        if (plain_only) 
            for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++) {
                Codec codec = txt_in_ft_by_dt[dt][i].codec;
                if (codec != CODEC_BGZF && codec != CODEC_GZ && codec != CODEC_BZ2 && codec != CODEC_XZ && codec != CODEC_ZIP && codec != CODEC_ORA)
                    SNPRINTF (s, "%s ", &file_exts[txt_in_ft_by_dt[dt][i].in][1]);
            }
    
        else {
            SNPRINTF (s, "\n%-8s: ", dt_name (dt));

            for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
                SNPRINTF (s, "%s ", &file_exts[txt_in_ft_by_dt[dt][i].in][1]);
        }
    }

    if (plain_only)
        SNPRINTF0 (s, "23andme generic ");
    else
        SNPRINTF0 (s, "\n23andMe : 23andme 23andme.zip"
                      "\nOther   : generic");
    return s;
}

FileType file_get_type (rom filename)
{
    if (!filename) return UNKNOWN_FILE_TYPE;

    // 23andme files have the format "genome_Firstname_Lastname_optionalversion_timestamp.txt" or .zip
    if (strstr (filename, "genome") && strstr (filename, "Full")) {
        if (filename_has_ext (filename, ".txt")) return ME23;
        if (filename_has_ext (filename, ".zip")) return ME23_ZIP;
        if (filename_has_ext (filename, ".txt.genozip")) return ME23_GENOZIP;
    }

    for (FileType ft=UNKNOWN_FILE_TYPE+1; ft < AFTER_LAST_FILE_TYPE; ft++) {

        // files that end with .txt/.txt.genozip/.zip are not classified as ME23, we already handled ME23 above
        if (ft == ME23 || ft == ME23_ZIP || ft == ME23_GENOZIP) continue; 

        if (filename_has_ext (filename, file_exts[ft])) 
            return ft;
    }

    return UNKNOWN_FILE_TYPE; // this should never happen, because GNRIC_ is "" so it will catch everything
}

static FileType file_get_type_of_generic (rom filename)
{
    if (!filename) return UNKNOWN_FILE_TYPE;

    if (filename_has_ext (filename, ".gz"))  return GNRIC_GZ;
    if (filename_has_ext (filename, ".bz2")) return GNRIC_BZ2;
    if (filename_has_ext (filename, ".xz"))  return GNRIC_XZ;
    if (filename_has_ext (filename, ".zip")) return GNRIC_ZIP;
    else                                     return GNRIC;
}

static FileType file_get_type_force_dt (rom filename, DataType dt)
{
    FileType ft = file_get_type (filename); 

    // if ft cannot be determined by filename (i.e. comes out as generic), get it by data type
    switch (ft) {
        case GNRIC_GZ  : return file_get_txt_ft_by_codec (dt, CODEC_GZ);
        case GNRIC_BZ2 : return file_get_txt_ft_by_codec (dt, CODEC_BZ2);
        case GNRIC_XZ  : return file_get_txt_ft_by_codec (dt, CODEC_XZ);
        case GNRIC_ZIP : return file_get_txt_ft_by_codec (dt, CODEC_ZIP);
        default        : return ft;
    }
}

// returns the filename without the extension eg myfile.1.sam.gz -> myfile.1. 
// if raw_name is given, memory is allocated sufficiently to concatenate a extension. Otherwise, filename is overwritten
uint32_t file_get_raw_name_and_type (rom filename, rom *raw_name, FileType *out_ft)
{
    unsigned len = strlen (filename);

    if (raw_name) {
        *raw_name = MALLOC (len + 30);
        memcpy (*(char **)raw_name, filename, len);
        (*(char **)raw_name)[len] = 0;
    }
    else 
        raw_name = &filename; // overwrite filename

    FileType ft = file_get_type (filename);
    if (ft != UNKNOWN_FILE_TYPE) {
        len -= strlen (file_exts[ft]);
        (*(char **)raw_name)[len] = 0;
    }

    if (out_ft) *out_ft = ft;

    return len;
}

static void file_ask_user_to_confirm_overwrite (rom filename)
{
    if (!strcmp (filename, "/dev/null")) return; // don't ask for /dev/null
    
    fprintf (stderr, "%s: output file %s already exists: in the future, you may use --force to overwrite\n", global_cmd, filename);
    
    if (!isatty(0) || !isatty(2)) exit_on_error(false); // if we stdin or stderr is redirected - we cannot ask the user an interactive question
    
    if (!str_query_user_yn ("Do you wish to overwrite it now?", QDEF_NO)) {
        fprintf (stderr, "No worries, I'm stopping here - no damage done!\n");
        exit (EXIT_OK);
    }

    fprintf (stderr, "\n");
}

static void file_redirect_output_to_stream (FileP file, rom exec_name, 
                                            rom format_option_0, rom format_option_1, rom format_option_2, rom format_option_3)
{
    char threads_str[20];
    snprintf (threads_str, sizeof (threads_str), "%u", global_max_threads);

    FILE *redirected_stdout_file = NULL;
    if (!flag.to_stdout) {
        redirected_stdout_file = fopen (file->name, file->mode); // exec_name will redirect its output to this file
        ASSINP (redirected_stdout_file, "cannot open file \"%s\": %s", file->name, strerror(errno));
    }

    char reason[40];
    snprintf (reason, sizeof (reason), "To output a %s file", file_exts[file->type]);
    output_compressor = stream_create (0, 0, 0, DEFAULT_PIPE_SIZE, 
                                       redirected_stdout_file, // output is redirected unless flag.to_stdout
                                       0, false, reason,
                                       exec_name, 
                                       format_option_0, 
                                       "--threads", threads_str,
                                       format_option_1, format_option_2, format_option_3,
                                       NULL);
    
    file->file = stream_to_stream_stdin (output_compressor);
}

// starting samtools 1.10, a PG record is added to the SAM header every "samtools view", and an option, --no-PG, 
// is provided to avoid this. See: https://github.com/samtools/samtools/releases/
// returns "--no-PG" if the option exists, or NULL if not
static rom file_samtools_no_PG (void)
{
    static rom ret_str[2] = { NULL, "--no-PG" };
    static int has_no_PG = -1; // unknown

    if (has_no_PG >= 0) return ret_str[has_no_PG]; // we already tested

    #define SAMTOOLS_HELP_MAX_LEN 20000
    char samtools_help_text[SAMTOOLS_HELP_MAX_LEN];

    #define MIN_ACCEPTABLE_LEN 100
    int len=0;

    for (unsigned i=1; i < 15 && len < MIN_ACCEPTABLE_LEN; i++) {
        // Tested on samtools 1.11: The normal way to see help is "samtools help view" however it fails if stdin is not the terminal. 
        // Instead, we use samtools view --threads, invalidly without an argument. This *sometimes* shows the help, and sometimes
        // just shows one line "samtools view:". We overcome this by repeating if the response is not long enough.
        StreamP samtools = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 0, 0, "To read/write CRAM files",
                                          "samtools", "view", "--threads", NULL);
        usleep (50000 * i); // wait for samtools

        // read both stderr and stdout from samtools
        len  = read (fileno (stream_from_stream_stderr (samtools)), samtools_help_text,       SAMTOOLS_HELP_MAX_LEN-1);
        len += read (fileno (stream_from_stream_stdout (samtools)), &samtools_help_text[len], SAMTOOLS_HELP_MAX_LEN-len-1);

        stream_close (&samtools, STREAM_DONT_WAIT_FOR_PROCESS);
    } 
    samtools_help_text[len] = '\0'; // terminate string (more portable, strnstr and memmem are non-standard)

    ASSERT0 (len >= MIN_ACCEPTABLE_LEN, "no response from \"samtools view --threads\"");

    return ret_str[(has_no_PG = !!strstr (samtools_help_text, "--no-PG"))];
}

// show meaningful error if file is not a supported type and return TRUE if it file should be skipped
static bool file_open_txt_read_test_valid_dt (ConstFileP file)
{ 
    if (file->data_type == DT_NONE) { 

        if (flag.multiple_files || tar_zip_is_tar()) {
            if (filename_has_ext (file->name, ".genozip")) {
            
                // case: --tar - include .genozip files verbatim
                if (tar_zip_is_tar()) {
                    tar_copy_file (file->name, file->name);
                    RETURNW (false, true, "Copied %s to the tar file", file_printname(file));
                }    
                else
                    RETURNW (false, true, "Skipping %s - it is already compressed", file_printname(file));

            }

            RETURNW (false, true, "Skipping %s - genozip doesn't know how to compress this file type (use --input to tell it)", 
                    file_printname (file));
        }
        else {
            ASSINP (!filename_has_ext (file->name, ".genozip"), 
                    "cannot compress %s because it is already compressed", file_printname(file));

            ABORT0 ("Unexpectedly, data_type==DT_NONE"); // not expecting to ever reach here, bc if file is not recognized, it should have been set to GENERIC
        }
    }

    return false; // all good - no need to skip this file
}

static void file_set_filename (FileP file, rom fn)
{
    // copy filename 
    unsigned fn_size = strlen (fn) + 1; // inc. \0
    file->name = MALLOC (fn_size);
    memcpy (file->name, fn, fn_size);
}

static void file_initialize_txt_file_fields (FileP file)
{
    #define TXT_INIT(buf) ({ buf_set_promiscuous (&file->buf, "txt_file->" #buf); })

    if (IS_ZIP) {
        // initialize evb "promiscuous" buffers - i.e. buffers that can be allocated by any thread
        // promiscuous buffers must be initialized by the main thread, and buffer.c does not verify their integrity.
        TXT_INIT(vb_info[0]);
        TXT_INIT(vb_info[1]);
    }
    else {

    }
}

static void file_open_ext_decompessor (FileP file, rom exec_name, rom subcommand, Codec streamed_codec, bool name_if_not_remote, rom args[7])
{
    char reason[64]; // used for error message if stream_create fails
    snprintf (reason, sizeof(reason), "To compress a %s file", codec_name (file->src_codec));

    input_decompressor = 
        stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0,
                       file->is_remote ? file->name : NULL,      // url                                        
                       file->redirected, reason, exec_name,
                       subcommand ? subcommand : SKIP_ARG,
                       #define A(i) (args[i] ? args[i] : SKIP_ARG)
                       A(0), A(1), A(2), A(3), A(4), A(5), A(6), 
                       (name_if_not_remote && !file->is_remote) ? file->name : SKIP_ARG,
                       NULL); 
    
    file->file       = stream_from_stream_stdout (input_decompressor);
    file->redirected = true;
    file->effective_codec = streamed_codec; // data received from input_decompressor is in this codec
}

static void file_open_txt_read_bz2 (FileP file)
{
    file->file = file->is_remote  ? BZ2_bzdopen (fileno (url_open_remote_file (NULL, file->name)), READ) // we're abandoning the FILE structure (and leaking it, if libc implementation dynamically allocates it) and working only with the fd
               : file->redirected ? BZ2_bzdopen (STDIN_FILENO, READ)
               :                    BZ2_bzopen (file->name, READ);  // for local files we decompress ourselves   

    ASSERT (file->file, "failed to open BZ2 file %s", file->name);
        
    if (!file->is_remote && !file->redirected) {
        int fd = BZ2_get_fd (file->file);
        
        stream_set_inheritability (fd, false); // Windows: allow file_remove in case of --replace
        #ifdef __linux__
        posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL); // ignore errors
        #endif
    }
}

static void file_open_txt_read_gz (FileP file)
{
    file->file = file->is_remote  ? url_open_remote_file (NULL, file->name)  
               : file->redirected ? fdopen (STDIN_FILENO,  "rb") 
               :                    fopen (file->name, READ);

    ASSERT (file->file, "failed to open %s: %s", file->name, strerror (errno));

    if (!file->is_remote && !file->redirected) {
        stream_set_inheritability (fileno (file->file), false); // Windows: allow file_remove in case of --replace
        #ifdef __linux__
        posix_fadvise (fileno ((FILE *)file->file), 0, 0, POSIX_FADV_SEQUENTIAL); // ignore errors
        #endif
    }

    // case: discovery deferred to the end of segconf when we know segconf.tech
    if (file->data_type == DT_FASTQ) {  // note: even if --no-bgzf: so we can report correct src_codec in stats
        file->effective_codec = file->src_codec = txtfile_is_gzip (file) ? CODEC_GZ : CODEC_NONE; // based on the first 3 bytes
        file->discover_during_segconf = (file->effective_codec == CODEC_GZ);
        txtfile_initialize_igzip (file);
    }

    // run discovery now if not FASTQ. That's because other data types might have header
    // which is read before segconf. luckily FASTQ doesn't.
    else
        txtfile_discover_specific_gz (file); // decide between GZ, BGZF and NONE
}

FileP file_open_txt_read (rom filename)
{
    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype   = TXT_FILE;
    file->redirected  = !filename; // later on, also CRAM, XZ, BCF will be set as redirected
    file->mode        = READ;
    file->is_remote   = filename && url_is_url (filename);
    flag.from_url     = file->is_remote;

    int64_t url_file_size = 0; // will be -1 if the web/ftp site does not provide the file size
    rom error = NULL;
    
    thool is_file_exists = unknown;

    if (file->is_remote) {
        error = url_get_status (filename, &is_file_exists, &url_file_size); // accessing is expensive - get existance and size in one call
        if (!error && url_file_size >= 0) file->disk_size = (uint64_t)url_file_size;
    }
    
    else if (!file->redirected) { // not stdin
        is_file_exists = file_exists (filename);
        error = strerror (errno);

        if (is_file_exists) 
            file->disk_size = file_get_size (filename);
    }

    // return null if genozip input file size is known to be 0, so we can skip it. note: file size of url might be unknown
    if (is_file_exists == yes && !file->disk_size && !url_file_size && 
        !(!file->is_remote && !file->redirected && file_is_fifo (filename)))  // a fifo is allowed is be size 0 (as it always is) 
        goto fail;    

    if (!file->redirected) {
        ASSINP (is_file_exists != no, "Failed to open \"%s\" for reading: %s", filename, error);
    
        file_set_filename (file, filename);
    
        file->basename = filename_base (file->name, false, "", NULL, 0);

        // if user provided the type with --input, we use that, otherwise derive from the file name
        file->type = flag.stdin_type == GNRIC ? file_get_type_of_generic (file->basename)
                   : flag.stdin_type          ? flag.stdin_type 
                   :                            file_get_type (file->basename);
    }

    else {   // stdin 
        file->basename = filename_base (NULL, false, FILENAME_STDIN, NULL, 0);
        file->type = flag.stdin_type; 
    }

    file->data_type = file_get_data_type_of_input_file (file->type);

    // show meaningful error if file is not a supported data type
    if (file_open_txt_read_test_valid_dt (file)) goto fail; // skip this file
    
    // open the file, based on the codec (as guessed by file extension)
    file->src_codec       = file_get_codec_by_txt_ft (file->data_type, file->type, true);
    file->effective_codec = file->src_codec; // initialize: can be changed if streaming or if gz variant

    switch (file->src_codec) { 
        case CODEC_GZ: case CODEC_BGZF: case CODEC_BAM: case CODEC_NONE: gz: 
            file_open_txt_read_gz (file);
            break;
        
        case CODEC_BZ2:
            file_open_txt_read_bz2 (file);
            break;

        case CODEC_CRAM: {
            // note: in CRAM, we read the header in advance in possible, directly (without samtools), so we can handle the case 
            // that the reference file is wrong. In samtools, if we read beyond the header with a wrong ref, samtools will hang.
            if (!file->is_remote && !file->redirected) {
                cram_inspect_file (file); // if file is indeed CRAM, updates file->est_num_lines, file->header_size, and if not, updates file->data_type and file->codec/src_codec
                if (file->src_codec == CODEC_GZ || file->src_codec == CODEC_NONE) goto gz; // actually, this is a GZ file (possibly BAM)
            }
            
            StrTextSuperLong samtools_T_option = cram_get_samtools_option_T();

            file_open_ext_decompessor (file, "samtools", "view", CODEC_BGZF, true, (rom[7]){  
                                       "--bam", "--uncompressed",           // BAM with BGZF blocks in which the payload is not compressed
                                       "--threads=10",                      // in practice, samtools is able to consume ~4 cores
                                       file_samtools_no_PG() ? "--no-PG" : SKIP_ARG, // don't add a PG line to the header 
                                       samtools_T_option.s[0] ? samtools_T_option.s : SKIP_ARG });
            
            if (flag.no_bgzf) {
                file->effective_codec = CODEC_GZ; 
                txtfile_initialize_igzip (file);
            }

            break;
        }

        case CODEC_BCF: {
            file_open_ext_decompessor (file, "bcftools", "view", CODEC_NONE, true, (rom[7]){ 
                                       "--threads=8", "-Ov", "--no-version" }); // BCF: do not append version and command line to the header
            break;
        }

        case CODEC_XZ:
            if (file->redirected) ABORTINP0 ("Compressing piped-in data in xz format is not currently supported");

            file_open_ext_decompessor (file, "xz", NULL, CODEC_NONE, true, (rom[7]){
                                       "--threads=8", "--decompress", "--keep", "--stdout", 
                                       flag.quiet ? "--quiet" : SKIP_ARG }); 
            break;

        case CODEC_ZIP:
            file_open_ext_decompessor (file, "unzip", NULL, CODEC_NONE, true, (rom[7]){
                                       "-p", flag.quiet ? "--quiet" : SKIP_ARG }); 
            break;

        case CODEC_ORA: {
            file_open_ext_decompessor (file, "orad", NULL, CODEC_NONE, false, (rom[7]){  
                                       "--raw", "--quiet", "--stdout",
                                       "--threads", str_int_s (global_max_threads).s, 
                                       (file->is_remote || file->redirected) ? "-" : file->name });    // local file name 
            break;
        }

        default:
            ABORT ("%s: invalid filename extension for %s files: %s", global_cmd, dt_name (file->data_type), file->name);
    }
    
    if (!file->file) goto fail;

    file_initialize_txt_file_fields (file);

    if (file->is_remote) FREE (error); // allocated by url_get_status

    return file;

fail:
    if (file->is_remote) FREE (error); 
    FREE (file->name);
    FREE (file->basename);
    FREE (file);
    return NULL;
}

FileP file_open_txt_write (rom filename, DataType data_type, MgzipLevel bgzf_level)
{
    ASSERT (data_type > DT_NONE && data_type < NUM_DATATYPES ,"invalid data_type=%d", data_type);

    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype  = TXT_FILE;
    file->mode       = WRITE;
    file->data_type  = data_type;
    file->redirected = !filename;

    file->effective_codec = data_type == DT_CRAM       ? CODEC_CRAM 
                          : data_type == DT_BCF        ? CODEC_BCF
                          : bgzf_level != BGZF_NO_BGZF ? CODEC_BGZF // see mgzip_piz_calculate_mgzip_flags
                          : /* BGZF_NO_BGZF */           CODEC_NONE;
    
    if (!file->redirected) { // not stdout
        if (file_exists (filename)   && 
            !file_is_fifo (filename) && // a fifo can be "overwritten" (that's just normal writing to a fifo)
            !flag.force && !flag.test)
            
            file_ask_user_to_confirm_overwrite (filename); // function doesn't return if user responds "no"

        file_set_filename (file, filename);

        file->type = file_get_type_force_dt (filename, data_type);
    }

    else  // stdout
        file->type = txt_out_ft_by_dt[data_type][bgzf_level >= 1]; // plain file or .gz file
    
    file->basename = filename_base (file->name, false, FILENAME_STDOUT, NULL, 0);

    // don't actually open the output file if we're not going to write to it
    if (flag.no_writer) return file;

    // open the file, based on the codec 
    switch (file->effective_codec) { 
        case CODEC_BGZF : 
        case CODEC_NONE : file->file = file->redirected ? fdopen (STDOUT_FILENO, "wb") : fopen (file->name, WRITE); break;

        case CODEC_CRAM : {
            StrTextSuperLong samtools_T_option = cram_get_samtools_option_T();
            file_redirect_output_to_stream (file, "samtools", "view", "-OCRAM", 
                                            file_samtools_no_PG(), 
                                            samtools_T_option.s[0] ? samtools_T_option.s : NULL); 
            break;
        }
        
        case CODEC_BCF  : {
            char comp_level[4] = { '-', 'l', '0' + MIN_(bgzf_level, 9), 0 };

            if (flag_show_bgzf)
                iprintf ("%s: launching external compressor \"bcftools\" with bgzf_level=%d\n", file->basename, bgzf_level);
            
            file_redirect_output_to_stream (file, "bcftools", "view", "-Ob", comp_level, NULL); 
            break;
        }

        default: {} // never reaches here
    }

    ASSINP (file->file, "cannot open file \"%s\": %s", file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    file_initialize_txt_file_fields (file);

    return file;
}

// note: we insert all the z_file buffers into the buffer list in advance and mark them as promiscuous, to avoid this 
// thread satety issue: without this pre-allocation, some of these buffers will be first allocated by the first 
// compute thread to use it, causing buf_alloc to modify evb's buf_list - this is not permitted as the main 
// thread might be doing so concurrently resulting in a corrupted evb.buf_list.

static void file_initialize_z_file_data (FileP file)
{
    init_dict_id_to_did_map (file->d2d_map); 
    profiler_new_z_file();
        
    #define Z_INIT(buf) ({ buf_set_promiscuous (&file->buf, "z_file->" #buf); })

    if (file->mode != READ) { // careful not to use IS_ZIP - which is set when reading aux files
        for (Did did_i=0; did_i < MAX_DICTS; did_i++) {
            ctx_zip_init_promiscuous (&file->contexts[did_i]); // must be done from main thread
            file->contexts[did_i].vb_1_pending_merges = -1;    // uninitialized - will be initialized in ctx_set_vb_1_pending_merges
        }
        __atomic_thread_fence (__ATOMIC_RELEASE); // release all vb_1_pending_merges
        
        // initialize evb "promiscuous" buffers - i.e. buffers that can be allocated by any thread (obviously protected by eg a mutex)
        // promiscuous buffers must be initialized by the main thread, and buffer.c does not verify their integrity.
        Z_INIT (ra_buf);
        Z_INIT (sag_grps);
        Z_INIT (sag_grps_index);
        Z_INIT (sag_alns);
        Z_INIT (sag_qnames);
        Z_INIT (sag_cigars); // union with solo_data
        Z_INIT (sag_seq);
        Z_INIT (sag_qual);
        Z_INIT (deep_index);
        Z_INIT (deep_ents);
        Z_INIT (vb_num_deep_lines);
        Z_INIT (section_list);
        Z_INIT (contexts[CHROM].chrom2ref_map);
        Z_INIT (vb_info[SAM_COMP_MAIN]);
    }
    else { // PIZ
        Z_INIT (sag_qnames);
        Z_INIT (sag_qual);
        Z_INIT (sag_cigars); // union with solo_data
    }

    if (flag_no_biopsy_line) // no need to initialize in --biopsy-line (as destroying it later will error)
        serializer_initialize (file->digest_serializer); 

    clock_gettime (CLOCK_REALTIME, &file->start_time);
}

// get time since creation of z_file object in memory
StrText file_get_z_run_time (FileP file)
{
    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 

    int seconds_so_far = ((tb.tv_sec - file->start_time.tv_sec)*1000 + 
                         ((int64_t)tb.tv_nsec - (int64_t)file->start_time.tv_nsec) / 1000000) / 1000; 

    return str_human_time (seconds_so_far, true);
}

// opens z_file for read or write
FileP file_open_z_read (rom filename)
{
    START_TIMER;

    ASSINP0 (filename, "it is not possible to redirect genozip files from stdin");

    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype = Z_FILE;
    file->mode      = READ;
    file->is_in_tar = (flag.t_offset > 0);

    if (flag.debug_tar) 
        iprintf ("file_open_z_read: t_offset=%"PRIu64" t_size=%"PRIu64" %s\n", flag.t_offset, flag.t_size, filename);

    rom disk_filename = file->is_in_tar ? tar_get_tar_name() : filename;
    ASSINP (file_exists (disk_filename), "cannot open \"%s\" for reading: %s", disk_filename, strerror (errno));

    file->disk_size = file->is_in_tar ? flag.t_size : file_get_size (filename);

    file_set_filename (file, filename);

    file->type = file_get_type (file->name);

    file->basename = filename_base (file->name, false, NULL, NULL, 0);

    // if a FASTA file was given as an argument to --reference or --REFERENCE, get the .ref.genozip file,
    // possobily running --make-reference in a separate process if needed
    if (flag.reading_reference && (file_get_data_type_of_input_file (file_get_type (file->name)) == DT_FASTA) && (file_get_type (file->name) != FASTA_GENOZIP)) 
        disk_filename = ref_fasta_to_ref (file);

    ASSINP (!flag.reading_reference || filename_has_ext (file->name, REF_GENOZIP_), 
            "You specified file \"%s\", however with --reference or --REFERENCE, you must specify a reference file (%s file or FASTA file)\n"
            "Tip: To create a genozip reference file from a FASTA file, use 'genozip --make-reference myfasta.fa'",
            file->name, REF_GENOZIP_);

    if ((!flag.seg_only && !flag.show_bam) || flag_loading_auxiliary) {

        // make sure file is a regular file (not FIFO, directory etc)
        struct stat sb;
        int cause=0, stat_errno=0;
        if (stat (disk_filename, &sb)) {
            cause = 6; // stat failed
            stat_errno = errno;
        }

        if ((sb.st_mode & S_IFMT) != S_IFREG) cause=7; // not regular file

        if (!cause) {
            file->file = fopen (disk_filename, READ);

            stream_set_inheritability (fileno (file->file), false); // Windows: allow file_remove in case of --replace
        }
        
        // verify that this is a genozip file 
        // we read the Magic at the end of the file (as the magic at the beginning may be encrypted)
        uint32_t magic;
        if (cause ||
            (cause = 1 * !file->file) ||
            (cause = 2 * !sb.st_size) ||
            (cause = 3 * !file_seek (file, -(int)sizeof (magic), SEEK_END, READ, SOFT_FAIL)) || 
            (cause = 4 * !fread (&magic, sizeof (magic), 1, file->file)) ||
            (cause = 5 * (BGEN32 (magic) != GENOZIP_MAGIC && !(flag.show_headers && flag.force)))) {
            
            int fail_errno = errno;
            FCLOSE (file->file, disk_filename);

            if (flag.validate == VLD_REPORT_INVALID) 
                flag.validate = VLD_INVALID_FOUND;

            else if (flag.validate == VLD_NO_REPORT)
                exit (EXIT_INVALID_GENOZIP_FILE); // silent exit with error code, if even a single file is not valid

            rom cause_str =   cause==1 ? strerror (fail_errno)
                            : cause==2 ? "file is empty"
                            : cause==3 ? "file_seek failed"
                            : cause==4 ? "fread failed"
                            : cause==5 ? "Not a valid genozip file (bad magic)"
                            : cause==6 ? strerror (stat_errno)
                            : cause==7 ? "Not a regular file"
                            :            "no error";

            if (flag.multiple_files) {

                if (flag.validate == VLD_INVALID_FOUND) // outputs even if --quiet
                    iprintf ("Cannot open %s: %s\n", disk_filename, cause_str);

                else if (flag.validate == VLD_NONE) {    // silenced by --quiet
                    static int once=0;
                    WARN ("Skipping %s: %s%s", file->name, cause_str,
                            !(once++) ? " (--quiet to silence this message)" : "");
                }
                
                file_close (&file);
                
                goto done;
            }
            else { // single file
                if (flag.validate == VLD_REPORT_VALID)
                    exit (EXIT_INVALID_GENOZIP_FILE); // exit quietly - with a return code indicating invalidity
                else
                    ABORTINP ("Cannot open %s: %s %s", disk_filename, cause_str,
                              (cause==3 || cause==4) ? strerror(fail_errno) : "");
            }
        }
        
        // file is valid
        else if (flag.validate == VLD_REPORT_VALID)
            iprintf ("%s\n", file->name); // print just filename, so a script can use this output
    }

    file->data_type = DT_NONE; // we will get the data type from the genozip header, not by the file name
    
    file_initialize_z_file_data (file);

    ASSINP (file->file, "cannot open file \"%s\": %s", file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

done:
    COPY_TIMER_EVB (file_open_z);
    return file;
}

// opens z_file for read or write
FileP file_open_z_write (rom filename, FileMode mode, DataType data_type, Codec src_codec)
{
    START_TIMER;

    ASSINP0 (filename, "it is not possible to redirect genozip files to stdout");

    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype = Z_FILE;
    file->mode      = mode;
    file->is_in_tar = tar_zip_is_tar();
    
    if (file_exists (filename) && 
        !flag.force            && 
        !flag.zip_no_z_file    && // not zip with --seg-only
        !file->is_in_tar)   

        file_ask_user_to_confirm_overwrite (filename); // function doesn't return if user responds "no"

    file_set_filename (file, filename);

    file->type = file_get_type_force_dt (file->name, data_type);
    file->data_type = data_type; 
    file->src_codec = src_codec;

    file->basename = filename_base (file->name, false, NULL, NULL, 0);

    ASSINP (filename_has_ext (file->name, GENOZIP_EXT), 
            "file \"%s\" must have a " GENOZIP_EXT " extension", file_printname (file));
    
    // set file->type according to the data type, overriding the previous setting - i.e. if the user
    // uses the --output option, he is unrestricted in the choice of a file name
    file->type = file_get_z_ft_by_txt_in_ft (file->data_type, txt_file->type); 

    mutex_initialize (file->dicts_mutex);
    mutex_initialize (file->custom_merge_mutex);
    mutex_initialize (file->test_abbrev_mutex);
    mutex_initialize (file->zriter_mutex);
    
    if (!flag.zip_no_z_file) {

        if (flag.force && !file->is_in_tar) 
            unlink (file->name); // delete file if it already exists (needed in weird cases, eg symlink to non-existing file)

        // if we're writing to a tar file, we get the already-openned tar file
        if (file->is_in_tar)
            file->file = tar_open_file (file->name, file->name);
            // note: tar doesn't have a z_reread_file bc --pair and --deep are not yet supported with --tar

        else {
            file->file = fopen (file->name, file->mode);
            
            if (!flag.no_zriter) 
                file->z_reread_file = fopen (file->name, READ);

#ifndef _WIN32
            // set z_file permissions to be the same as the txt_file permissions (if possible)
            if (file->file && txt_file && txt_file->name && !txt_file->is_remote) {
                struct stat st;
                if (stat (txt_file->name, &st))
                    WARN ("FYI: Failed to set permissions of %s because failed to stat(%s): %s", file->name, txt_file->name, strerror(errno));
                
                else 
                    chmod (file->name, st.st_mode); // ignore errors (e.g. this doesn't work on NTFS)
            }
#endif
        }
    }
    
    file->genozip_version   = code_version_major(); // to allow the VER macro to operate consistently across ZIP/PIZ
    file->genozip_minor_ver = code_version_minor();

    file_initialize_z_file_data (file);

    ASSINP (file->file || flag.zip_no_z_file, 
            "cannot open file \"%s\": %s", file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    COPY_TIMER_EVB (file_open_z);
    return file;
}

// PIZ: index file is it is a disk file of a type that can be indexed
static void file_index_txt (ConstFileP file)
{
    ASSERTNOTNULL (file);

    RETURNW (file->name,, "%s: cannot create an index file when output goes to stdout", global_cmd);

    StreamP indexing = NULL;

    switch (file->data_type) {
        case DT_SAM:
        case DT_BAM: 
            RETURNW (file->effective_codec == CODEC_BGZF,, "%s: output file needs to be a .sam.gz or .bam to be indexed", global_cmd); 
            indexing = stream_create (0, 0, 0, 0, 0, 0, 0, "to create an index", "samtools", "index", file->name, NULL); 
            break;
            
        case DT_VCF: 
            RETURNW (file->effective_codec == CODEC_BGZF,, "%s: output file needs to be a .vcf.gz or .bcf to be indexed", global_cmd); 
            RETURNW (vcf_header_get_has_fileformat(),, "%s: file needs to start with ##fileformat=VCF be indexed", global_cmd); 
            indexing = stream_create (0, 0, 0, 0, 0, 0, 0, "to create an index", "bcftools", "index", file->name, NULL); 
            break;

        case DT_FASTQ:
        case DT_FASTA:
            RETURNW (file->effective_codec == CODEC_BGZF || file->effective_codec == CODEC_NONE,, 
                     "%s: To be indexed, the output file cannot be compressed with %s", global_cmd, codec_name (file->effective_codec)); 
            indexing = stream_create (0, 0, 0, 0, 0, 0, 0, "to create an index", "samtools", "faidx", file->name, NULL); 
            break;

        default: break; // we don't know how to create an index for other data types
    }

    if (indexing) {
        progress_new_component (file->basename, "Indexing", false, NULL);

        stream_wait_for_exit (indexing, false);

        progress_finalize_component_time ("Done indexing", DIGEST_NONE);
    }
}

void file_close (FileP *file_p) 
{
    START_TIMER;

    FileP file = *file_p;

    if (z_file && file == z_file && !flag_loading_auxiliary && 
        flag.show_time_comp_i == COMP_ALL && !flag.show_time[0]) // show-time without the optional parameter 
        profiler_add_evb_and_print_report();

    store_relaxed (*file_p, (FileP)NULL); 

    if (!file) return; // nothing to do

    if (file->file && file->supertype == TXT_FILE) {

        if (file->mode == READ && file->effective_codec == CODEC_BZ2)
            BZ2_bzclose((BZFILE *)file->file);
        
        else if (file->mode == READ && is_read_via_ext_decompressor (file)) 
            stream_close (&input_decompressor, STREAM_WAIT_FOR_PROCESS);

        else if (file->mode == WRITE && is_written_via_ext_compressor (file->effective_codec)) 
            stream_close (&output_compressor, STREAM_WAIT_FOR_PROCESS);

        // if its stdout - just flush, don't close - we might need it for the next file
        else if (file->mode == WRITE && flag.to_stdout) 
            fflush ((FILE *)file->file);

        else if (file->is_remote)
            url_close_remote_file_stream ((FILE**)&file->file);

        else
            FCLOSE (file->file, file_printname (file));

        // create an index file using samtools, bcftools etc, if applicable
        if (file->mode == WRITE && flag.index_txt && !flag_loading_auxiliary) 
            file_index_txt (file);
    }

    else if (file->file && file->supertype == Z_FILE) {

        // ZIP note: we need to destory all even if unused, because they were initialized in file_initialize_z_file_data
        if (IS_ZIP)
            for (Did did_i=0; did_i < (IS_ZIP ? MAX_DICTS : file->num_contexts); did_i++) 
                mutex_destroy (file->ctx_mutex[did_i]); 

        if (file->is_in_tar && file->mode != READ)
            tar_close_file (&file->file);
        else {
            FCLOSE (file->file, file_printname (file));
            FCLOSE (file->z_reread_file, file_printname (file));
        }
        serializer_destroy (file->digest_serializer);  
    }

    // free resources if we are NOT near the end of the execution. If we are at the end of the execution
    // it is faster to just let the process die

    if (!flag.let_OS_cleanup_on_exit) {

        if (IS_PIZ && flag.deep && file->supertype == Z_FILE) { // in this case, deep_index and deep_ents are Buffers containing arrays of Buffers
            for_buf (Buffer, buf, file->deep_index) buf_destroy (*buf);
            for_buf (Buffer, buf, file->deep_ents)  buf_destroy (*buf);  
        }

        buflist_destroy_file_bufs (file);

        mutex_destroy (file->zriter_mutex);
        mutex_destroy (file->dicts_mutex);
        mutex_destroy (file->custom_merge_mutex);
        mutex_destroy (file->test_abbrev_mutex);
        mutex_destroy (file->bgzf_discovery_mutex);

        FREE (file->name);
        FREE (file->basename);
        FREE (file);
    }

    COPY_TIMER_EVB (file_close);
}

void file_remove (rom filename, bool fail_quietly)
{
    chmod (filename, S_IRUSR | S_IWUSR); // make sure its +w so we don't get permission denied (ignore errors)

#ifndef _WIN32
    int ret = remove (filename); 
    ASSERTW (!ret || fail_quietly, "Warning: failed to remove %s: %s", filename, strerror (errno));
#else
    ASSERTW (DeleteFile (filename) || fail_quietly, "Warning: failed to remove %s: %s", filename, str_win_error());
#endif
}

bool file_rename (rom old_name, rom new_name, bool fail_quietly)
{
    chmod (old_name, S_IRUSR | S_IWUSR); // make sure its +w so we don't get permission denied (ignore errors)

    int ret = rename (old_name, new_name); 
    ASSERTW (!ret || fail_quietly, "Warning: failed to rename %s to %s: %s", old_name, new_name, strerror (errno));

    return !ret; // true if successful
}

// also updates filename to .gz (but not if .bam)
void file_gzip (char *filename)
{
    unsigned fn_len = strlen (filename);

    char command[fn_len + 50];

    int ret = 1;
    
    snprintf (command, sizeof (command), "bgzip -@%u -f \"%s\" %s", global_max_threads, filename, flag.is_windows ? "" : " > /dev/null 2>&1");
    ret = system (command); // note: runs sh on Unix, and cmd.exe on Windows
    
    if (ret && errno == ENOENT) { // no bgzip - try pigz
        snprintf (command, sizeof (command), "pigz -f \"%s\" %s", filename, flag.is_windows ? "" : " > /dev/null 2>&1");
        ret = system (command);
    }

    if (ret && errno == ENOENT) { // no pigz - try gzip
        snprintf (command, sizeof (command), "gzip -f \"%s\" %s", filename, flag.is_windows ? "" : " > /dev/null 2>&1");
        ret = system (command);
    }

    ASSERTW (!ret, "FYI: \"%s\" returned %d. No harm.", command, ret); 

    if (!ret) {
        // special case: rename .bam.gz -> .bam
        if (fn_len >= 4 && !memcmp (&filename[fn_len-4], ".bam", 4)) {
            char gz_filename[fn_len + 10];
            snprintf (gz_filename, sizeof (gz_filename), "%s.gz", filename);
            file_remove (filename, true);
            file_rename (gz_filename, filename, false);
        }
        else 
            strcpy (&filename[fn_len], ".gz");
    }
}

void file_mkfifo (rom filename)
{
#ifndef _WIN32
    file_remove (filename, true);
    ASSERT (!mkfifo (filename, 0666), "mkfifo failed for %s: %s", filename, strerror (errno));

#else
    ABORT0 ("file_mkfifo not supported on Windows");
#endif
}

bool file_is_fifo (rom filename)
{
    if (flag.is_windows) return false; // we don't support FIFOs in Win32 yet

    struct stat st;
    ASSERT (!stat (filename, &st), "stat failed on %s", filename);

    return S_ISFIFO (st.st_mode);
}

bool file_exists (rom filename)
{

    if (!filename || !filename[0]) return false;
    bool exists = !access (filename, F_OK);

#ifdef _WIN32
    // TO DO: overcome this limitation, see: https://docs.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation
    if (!exists && strlen (filename) > MAX_PATH)
        WARN_ONCE ("Genozip limitation: filenames on Windows are limited to %u characters. Please contact "EMAIL_SUPPORT" for advice: %s", PATH_MAX, filename);
#endif

    return exists;
}

// returns true if successful. depending on soft_fail, a failure will either emit an error 
// (and exit) or a warning (and return).
bool file_seek (FileP file, int64_t offset, 
                int whence, // SEEK_SET, SEEK_CUR or SEEK_END
                rom mode, // READ if seeking before reading, WRITE if before writing 
                FailType fail_type) 
{
    ASSERTNOTNULL (file);
    ASSERTNOTNULL (file->file);
    
    if (file->supertype == Z_FILE) {

        if (whence == SEEK_END && file->is_in_tar && IS_PIZ) { 
            offset += flag.t_offset + flag.t_size;
            whence = SEEK_SET;
            goto test_already_there;
        }

        // in SEEK_SET of a z_file that is being tarred, update the offset to the beginning of the file data in the tar file
        else if (whence == SEEK_SET) {
            offset += (IS_ZIP ? tar_file_offset() : flag.t_offset); // 0 if not using tar

            test_already_there:
            if (ftello64 (GET_FP(file, mode)) == offset) return true; // already at the right offset
        }
    }

    int ret = fseeko64 (GET_FP(file, mode), offset, whence);

    if (fail_type != HARD_FAIL) {
        if (!flag.to_stdout && fail_type==WARNING_FAIL) {
            ASSERTW (!ret, errno == EINVAL ? "Warning: Error while reading file %s (fseeko64 (whence=%d offset=%"PRId64")): it is too small%s" 
                                           : "Warning: fseeko64 failed on file %s (whence=%d offset=%"PRId64"): %s", 
                     file_printname (file), whence, offset, errno == EINVAL ? "" : strerror (errno));
        }
    } 
    else
        ASSERT (!ret, "fseeko64(offset=%"PRId64" whence=%d) failed on file %s (FILE*=%p remote=%s redirected=%s): %s", 
                offset, whence, file_printname (file), file->file, TF(file->is_remote), TF(file->redirected), strerror (errno));

    return !ret;
}

uint64_t file_get_size (rom filename)
{
    struct stat64 st;
    
    int ret = stat64(filename, &st);
    ASSERT (!ret, "stat64 failed on '%s': %s", filename, strerror(errno));
    
    return st.st_size;
}

bool file_is_dir (rom filename)
{
    ASSERTNOTNULL (filename);
    int filename_len = strlen (filename);
    
    // temporarily remove trailing / 
    if (filename[filename_len-1] == '/' || filename[filename_len-1] == '\\')
        filename_len--;

    SAFE_NULT (filename);

    struct stat64 st;
    int ret = stat64 (filename, &st); // 0 if successful

    SAFE_RESTORE;

    return !ret && S_ISDIR (st.st_mode);
}

void file_mkdir (rom dirname)
{
    if (file_is_dir (dirname)) return; // already exists - that's ok

#ifdef _WIN32
    int ret = _mkdir (flag.out_dirname);
#else
    int ret = mkdir (flag.out_dirname, 0777);
#endif
    ASSERT (!ret, "mkdir(%s) failed: %s", flag.out_dirname, strerror (errno));
}

// reads an entire file into a buffer. if filename is "-", reads from stdin
void file_get_file (VBlockP vb, rom filename, BufferP buf, rom buf_name,
                    uint64_t max_size, // 0 to read entire file, or specify for max size
                    FileContentVerificationType ver_type, bool add_string_terminator)
{
    bool is_stdin = !strcmp (filename, "-");
    if (is_stdin && !max_size) max_size = 10000000; // max size for stdin
    
    uint64_t file_size = is_stdin ? 0 : file_get_size (filename);

    uint64_t size = is_stdin   ? max_size
                  : !file_size ? max_size
                  : max_size   ? MIN_(max_size, file_size) 
                  :              file_size;

    buf_alloc (vb, buf, 0, size + add_string_terminator, char, 1, buf_name);

    FILE *file = is_stdin ? stdin : fopen (filename, "rb");
    ASSINP (file, "cannot open \"%s\": %s", filename, strerror (errno));

    buf->len = fread (buf->data, 1, size, file);
    ASSERT (is_stdin || max_size || buf->len == size, "Error reading file %s: %s", filename, strerror (errno));

    ASSINP (ver_type != VERIFY_ASCII || str_is_printable (STRb(*buf)), "Expecting %s to contain text (ASCII)", filename);
    ASSINP (ver_type != VERIFY_UTF8  || str_is_utf8      (STRb(*buf)), "Expecting %s to contain ASCII or UTF-8 text", filename);

    if (add_string_terminator) 
        *BAFTc (*buf) = 0;

    FCLOSE (file, filename);
}

// writes data to a file and flushes it, returns true if successful

static Mutex put_data_mutex = {};
#define MAX_PUT_DATA_FILES_PER_EXECUTION 10 // maximum files deletable at abort
static rom put_data_tmp_filenames[MAX_PUT_DATA_FILES_PER_EXECUTION] = {};
static unsigned num_put_files=0;

bool file_put_data (rom filename, const void *data, uint64_t len, 
                    mode_t chmod_to) // optional - ignored if 0
{
    int fn_len = strlen (filename);

    // remove invalid characters from filename
    if (flag.is_windows)
        for (char *c=(char*)filename ; *c ; c++)
            if (*c == ':' && (c-filename != 1)) *c = '-'; // ':' exist eg in SAM AUX names 

    int tmp_filename_size = fn_len + 5;
    char *tmp_filename = MALLOC (tmp_filename_size);
    // we first write to tmp_filename, and after we complete and flush, we rename to the final name
    // so that if a file exists (in its final name) - then its guaranteed to be fully written
    snprintf (tmp_filename, tmp_filename_size, "%s.tmp", filename);

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

    // write in blocks (Windows hangs if the block is too big, a few GB)
    for (size_t written = 0; written < len; written += 1 MB) {
        uint64_t this_len = MIN_(1 MB, len - written); 
        uint64_t this_written = fwrite ((rom)data + written, 1, this_len, file);

        if (this_written != this_len) {
            WARN ("Failed to write %s: wrote only %"PRIu64" bytes of the expected %"PRIu64" (file removed)", tmp_filename, (uint64_t)written, len);
            put_data_tmp_filenames[my_file_i] = NULL; // no need to lock mutex
            remove (tmp_filename);
            return false;
        }
    }

    SAVE_VALUE (errno);
    fflush (file);    
    FCLOSE (file, tmp_filename); 
    RESTORE_VALUE (errno); // in cases caller wants to print fwrite error

    // we can't enter if file_put_data_abort is active, or it needs to wait for us before deleting tmp files
    mutex_lock (put_data_mutex);

    remove (filename);
    int renamed_failed = rename (tmp_filename, filename);

    put_data_tmp_filenames[my_file_i] = NULL; // remove tmp file name from list 
    
    mutex_unlock (put_data_mutex);

    ASSERT (!renamed_failed, "Failed to rename %s to %s: %s", tmp_filename, filename, strerror (errno));
    FREE (tmp_filename);

    if (chmod_to) 
        ASSERT (!chmod (filename, chmod_to), "Failed to chmod %s: %s", filename, strerror (errno));
    
    return true;
}

// error handling: unlink files currently being written (the actual writing will terminate when the thread terminates)
void file_put_data_abort (void)
{
    if (!put_data_mutex.initialized) return;

    mutex_lock (put_data_mutex);

    for (unsigned i=0; i < num_put_files; i++) 
        if (put_data_tmp_filenames[i]) {
            // TODO: this works on Linux but not Windows - gets "Permission Denied" if file_put_data is in fflush()
            unlink (put_data_tmp_filenames[i]); // ignore errors
            remove (put_data_tmp_filenames[i]); // in case unlinked failed - eg NTFS - ignore errors
        }

    // mutex remains locked - no more files can be put after this point
}

void file_put_data_reset_after_fork (void)
{
    put_data_mutex = (Mutex){};
}

PutLineFn file_put_line (VBlockP vb, STRp(line), rom msg)
{
    PutLineFn fn;
    snprintf (fn.s, sizeof (fn.s), "line.%u.%d.%s%s", vb->vblock_i, vb->line_i, command==ZIP ? "zip" : "piz",
             file_plain_ext_by_dt ((VB_DT(SAM) && z_file->z_flags.txt_is_bin) ? DT_BAM : vb->data_type));
    
    file_put_data (fn.s, STRa(line), 0);

    if (IS_PIZ)
        WARN ("\n%s line=%s line_in_file(1-based)=%"PRId64". Dumped %s (dumping first occurance only)", 
                msg, line_name(vb).s, writer_get_txt_line_i (vb, vb->line_i), fn.s);
    else
        WARN ("\n%s line=%s vb_size=%u MB. Dumped %s", msg, line_name(vb).s, (int)(segconf.vb_size >> 20), fn.s);

    return fn;
}

void file_assert_ext_decompressor (void)
{
    if (!stream_wait_for_exit (input_decompressor, false)) return; // just a normal EOF - all good!

    if (flag.truncate) return; // truncated as requested - all good

    // read error from stderr
    #define INPUT_DECOMPRESSOR_RESPSONSE_LEN 4096
    char error_str[INPUT_DECOMPRESSOR_RESPSONSE_LEN];

    FILE *stderr_pipe = stream_from_stream_stderr (input_decompressor);
    int bytes_read = fread (error_str, 1, INPUT_DECOMPRESSOR_RESPSONSE_LEN-1, stderr_pipe);
    error_str[bytes_read] = 0; // string terminator

    ABORT ("%s: failed to read file: %s\n%s: %s", 
           global_cmd, txt_name, stream_get_exec_name (input_decompressor), error_str);
}

// used when aborting due to an error. avoid the compressors outputting their own errors after our process is gone
void file_kill_external_compressors (void)
{
    stream_close (&input_decompressor, STREAM_KILL_PROCESS);
    stream_close (&output_compressor,  STREAM_KILL_PROCESS);
}

rom ft_name (FileType ft)
{
    return type_name (ft, &file_exts[ft], ARRAY_LEN(file_exts));
}

rom file_plain_ext_by_dt (DataType dt)
{
    FileType plain_ft = txt_in_ft_by_dt[FAF ? DT_FASTA : dt][0].in;

    return file_exts[plain_ft];
}

bool file_buf_locate (FileP file, ConstBufferP buf)
{
    return is_p_in_range (buf, file, sizeof (File));
}
