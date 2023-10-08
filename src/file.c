// ------------------------------------------------------------------
//   file.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif
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
#include "writer.h"
#include "version.h"
#include "filename.h"
#include "huffman.h"
#include "arch.h"

// globals
FileP z_file   = NULL;
FileP txt_file = NULL;

static StreamP input_decompressor = NULL; // bcftools or xz, unzip or samtools - only one at a time
static StreamP output_compressor  = NULL; // samtools (for cram), bcftools

// global pointers - so the can be compared eg "if (mode == READ)"
rom READ      = "rb";  // use binary mode (b) in read and write so Windows doesn't add \r
rom WRITE     = "wb";
rom WRITEREAD = "wb+"; // only supported for z_file and gencomp disk files

rom file_exts[] = FILE_EXTS;

static const struct { FileType in; Codec codec; FileType out; } txt_in_ft_by_dt[NUM_DATATYPES][50] = TXT_IN_FT_BY_DT;
static const FileType txt_out_ft_by_dt[NUM_DATATYPES][20] = TXT_OUT_FT_BY_DT;
static const FileType z_ft_by_dt[NUM_DATATYPES][20] = Z_FT_BY_DT;

// get data type by file type
DataType file_get_data_type (FileType ft, bool is_input)
{
    // note: if make-reference, we scan the array from dt=0 (DT_REF), otherwise we ignore DT_REF
    for (DataType dt=!flag.make_reference; dt < NUM_DATATYPES; dt++) 
        for (unsigned i=0; (is_input ? txt_in_ft_by_dt[dt][i].in : txt_out_ft_by_dt[dt][i]); i++)
            if ((is_input ? txt_in_ft_by_dt[dt][i].in : txt_out_ft_by_dt[dt][i]) == ft) 
                return dt;

    return (is_input && ft == UNKNOWN_FILE_TYPE) ? DT_GENERIC : DT_NONE;
}

rom file_plain_text_ext_of_dt (DataType dt) 
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
    static char s[1000];
        
    int len = 0;
    for (DataType dt=1; dt < NUM_DATATYPES; dt++) { // start from 1, excluding DT_REFERENCE
        
        if (dt == DT_GENERIC || dt == DT_ME23) continue;

        if (plain_only) 
            for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++) {
                Codec codec = txt_in_ft_by_dt[dt][i].codec;
                if (codec != CODEC_BGZF && codec != CODEC_GZ && codec != CODEC_BZ2 && codec != CODEC_XZ && codec != CODEC_ZIP)
                    len += sprintf (&s[len], "%s ", &file_exts[txt_in_ft_by_dt[dt][i].in][1]);
            }
            // sprintf (&s[strlen(s)], "%s ", &file_exts[txt_in_ft_by_dt[dt][0].in][1]);
    
        else {
            sprintf (&s[strlen (s)], "\n%-8s: ", dt_name (dt));

            for (unsigned i=0; txt_in_ft_by_dt[dt][i].in; i++)
                len += sprintf (&s[len], "%s ", &file_exts[txt_in_ft_by_dt[dt][i].in][1]);
        }
    }

    if (plain_only)
        sprintf (&s[strlen(s)], "23andme generic ");
    else
        sprintf (&s[strlen(s)], "\n23andMe : 23andme 23andme.zip"
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
void file_get_raw_name_and_type (rom filename, rom *raw_name, FileType *out_ft)
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
    if (ft != UNKNOWN_FILE_TYPE) 
        (*(char **)raw_name)[len - strlen (file_exts[ft])] = 0;

    if (out_ft) *out_ft = ft;
}

void file_set_input_size (rom size_str)
{   
    unsigned len = strlen (size_str);

    for (unsigned i=0; i < len; i++) {
        ASSINP (IS_DIGIT (size_str[i]), "expecting the file size in bytes to be a positive integer: %s", size_str);

        flag.stdin_size = flag.stdin_size * 10 + (size_str[i] - '0'); 
    }
}

void file_set_input_type (rom type_str)
{
    char ext[strlen (type_str) + 2]; // +1 for . +1 for \0
    sprintf (ext, ".%s", type_str);

    str_tolower (ext, ext); // lower-case to allow case-insensitive --input argument (eg vcf or VCF)

    if (!strcmp (ext, ".23andme")) 
        flag.stdin_type = ME23;

    else if (!strcmp (ext, ".23andme.zip")) 
        flag.stdin_type = ME23_ZIP;

    else if (!strcmp (ext, ".generic")) {
        flag.stdin_type = GNRIC;
        flag.explicitly_generic = true;
    }

    else {
        flag.stdin_type = file_get_type (ext); // we don't enforce 23andMe name format - any .txt or .zip will be considered ME23
        if (flag.stdin_type == GNRIC && !flag.explicitly_generic) 
            flag.stdin_type = UNKNOWN_FILE_TYPE; // This is an unknown input name - not generic
    }

    if (file_get_data_type (flag.stdin_type, true) != DT_NONE) return; // all good 

    // the user's argument is not an accepted input file type - print error message
    ABORT ("%s: --input (or -i) must be ones of these: %s", global_cmd, file_compressible_extensions(false));
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

static void file_redirect_output_to_stream (FileP file, char *exec_name, 
                                            rom stdout_option, rom format_option_1, rom format_option_2, rom format_option_3)
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
static bool file_open_txt_read_test_valid_dt (ConstFileP file)
{ 
    if (file->data_type == DT_NONE) { 

        if (flag.multiple_files || tar_zip_is_tar()) {
            if (filename_has_ext (file->name, ".genozip")) {
            
                // case: --tar - include .genozip files verbatim
                if (tar_zip_is_tar()) {
                    tar_copy_file (file->name);
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

static void file_initialize_txt_file_data (FileP file)
{
    #define TXT_INIT(buf) ({ buf_set_promiscuous (&file->buf, "txt_file->" #buf); })

    if (IS_ZIP) {
        mutex_initialize (file->recon_plan_mutex[0]);
        mutex_initialize (file->recon_plan_mutex[1]);

        // initialize evb "promiscuous" buffers - i.e. buffers that can be allocated by any thread
        // promiscuous buffers must be initialized by the main thread, and buffer.c does not verify their integrity.
        TXT_INIT(line_info[0]);
        TXT_INIT(line_info[1]);
        TXT_INIT(vb_info[0]);
        TXT_INIT(vb_info[1]);
    }
    else {

    }
}

FileP file_open_txt_read (rom filename)
{
    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype   = TXT_FILE;
    file->is_remote   = filename && url_is_url (filename);
    file->redirected  = !filename; // later on, also CRAM, XZ, BCF will be set as redirected
    file->mode        = READ;

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
        ASSINP (is_file_exists != no, "cannot open \"%s\" for reading: %s", filename, error);
    
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

    file->data_type = file_get_data_type (file->type, true);

    // show meaningful error if file is not a supported data type
    if (file_open_txt_read_test_valid_dt (file)) goto fail; // skip this file
    
    // open the file, based on the codec (as guessed by file extension)
    file->codec        = file_get_codec_by_txt_ft (file->data_type, file->type, false);
    file->source_codec = file_get_codec_by_txt_ft (file->data_type, file->type, true);

    switch (file->codec) { 
        case CODEC_CRAM: {
            if (!file->is_remote && !file->redirected)
                file->est_num_lines = cram_estimated_num_alignments (file->name); // for progress indication
                
            rom samtools_T_option = ref_get_cram_ref(gref);
            input_decompressor = stream_create (0, DEFAULT_PIPE_SIZE, DEFAULT_PIPE_SIZE, 0, 0, 
                                                file->is_remote ? file->name : NULL,      // url                                        
                                                file->redirected,
                                                "To decompress a CRAM file", 
                                                "samtools", "view", "-bu", "--threads", "10", "-h", // in practice, samtools is able to consume ~4 cores
                                                file_samtools_no_PG() ? "--no-PG" : "-h", // don't add a PG line to the header (just repeat -h if this is an older samtools without --no-PG - no harm)
                                                file->is_remote ? SKIP_ARG : file->name,  // local file name 
                                                samtools_T_option, NULL);
            file->file = stream_from_stream_stdout (input_decompressor);
            file->redirected = true;
            FREE(samtools_T_option);
            goto fallthrough_from_cram;
        }

        case CODEC_GZ:   // we test the first few bytes of the file to differentiate between NONE, GZ and BGZIP
        case CODEC_BGZF: 
        case CODEC_NONE: {
            file->file = file->is_remote  ? url_open (NULL, file->name)  
                       : file->redirected ? fdopen (STDIN_FILENO,  "rb") 
                       :                    fopen (file->name, READ);
            ASSERT (file->file, "failed to open %s: %s", file->name, strerror (errno));

fallthrough_from_cram: {}
            // read the first potential BGZF block to test if this is GZ or BGZF
            uint8_t block[BGZF_MAX_BLOCK_SIZE]; 
            uint32_t block_size;

            // we can't use libc read buffer unfortunately, because when reading the test bgzf block, if it turns out to be GZIP 
            // we will do gzdopen discarding the internal FILE buffer. We also cannot use setvbuf after the file is already partly read (?).
            // TO DO: this doesn't work over a pipe - so we can't pipe-in non-BGZF GZIP files (bug 243)
            setvbuf (file->file, 0, _IONBF, 0); 

            ASSERTNOTINUSE (evb->scratch);

            int32_t bgzf_uncompressed_size = bgzf_read_block (file, block, &block_size, SOFT_FAIL);

            // case: this is indeed a bgzf - we put the still-compressed data in vb->scratch for later consumption
            // in txtfile_read_block_bgzf
            if (bgzf_uncompressed_size > 0) {
                if (file->source_codec != CODEC_CRAM && file->source_codec != CODEC_BAM && file->source_codec != CODEC_BCF) 
                    file->source_codec = CODEC_BGZF;

                file->codec = CODEC_BGZF;
                
                evb->scratch.count = bgzf_uncompressed_size; // pass uncompressed size in param
                buf_add_more (evb, &evb->scratch, (char*)block, block_size, "scratch");
                
                file->bgzf_flags = bgzf_get_compression_level (file->name ? file->name : FILENAME_STDIN, 
                                                               block, block_size, (uint32_t)bgzf_uncompressed_size);

                // if we're compressing an level-0 BGZF (eg produced by samtools view -u), we don't record its
                // value so that PIZ reconstructs at the default level
                if (file->bgzf_flags.level == 0) 
                    file->bgzf_flags.level = BGZF_COMP_LEVEL_UNKNOWN;
            }

            // for regulars files, we already skipped 0 size files. This can happen in STDIN
            else if (bgzf_uncompressed_size == 0) {

                ASSINP (!flags_pipe_in_process_died(), // only works for Linux
                        "Pipe-in process %s (pid=%u) died without sending any data",
                        flags_pipe_in_process_name(), flags_pipe_in_pid());

                ABORTINP ("No data exists in input file %s", file->name ? file->name : FILENAME_STDIN);
            }

            // Bug 490: if a gzip- pr bz2-compressed file appears to be uncompressed (CODEC_NONE =
            // no .gz or .bz2 extension), the user is unaware that this file is compressed. Since we can't 
            // re-compress CODEC_GZ or CODEC_BZ2 upon decompressed, it will appear as if we're not recreating the original file. 
            
            // case: this is a non-BGZF gzip format - open with zlib and hack back the read bytes 
            // (note: we cannot re-read the bytes from the file as the file might be piped in)
            else if (bgzf_uncompressed_size == BGZF_BLOCK_GZIP_NOT_BGZIP) {
                
                ASSERT0 (!file->redirected, "genozip can't read gzip data from a pipe - piped data must be either plain or in BGZF format - i.e. compressed with bgzip, htslib etc");

                file->codec = file->source_codec = CODEC_GZ;
                
                int fd = dup (fileno((FILE *)file->file));
                ASSERT (fd >= 0, "dup failed: %s", strerror (errno));

                FCLOSE (file->file, file_printname (file));
                
                file->file  = gzdopen (fd, READ); 
                gzinject (file->file, block, block_size); // a hack - adds a 18 bytes of compressed data to the in stream, which will be consumed next, instead of reading from disk
            } 

            // case: this is not GZIP format at all. treat as a plain file, and put the data read in vb->scratch 
            // for later consumption is txtfile_read_block_plain
            else if (bgzf_uncompressed_size == BGZF_BLOCK_IS_NOT_GZIP) {

                #define BZ2_MAGIC "BZh"
                #define XZ_MAGIC  (char[]){ 0xFD, '7', 'z', 'X', 'Z', 0 }
                #define ZIP_MAGIC (char[]){ 0x50, 0x4b, 0x03, 0x04 }

                // we already open the file, so not easy to re-open with BZ2_bzopen as it would require injecting the read data into the BZ2 buffers
                if (str_isprefix_((rom)block, block_size, BZ2_MAGIC, 3)) 
                    ABORTINP0 ("The data seems to be in bz2 format. Please use --input to specify the type (eg: \"genozip --input sam.bz2\")");

                if (str_isprefix_((rom)block, block_size, XZ_MAGIC, 6)) {
                    if (file->redirected) ABORTINP0 ("Compressing piped-in data in xz format is not currently supported");
                    if (file->is_remote) ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\")");
                    ABORTINP0 ("The data seems to be in xz format. Please use --input to specify the type (eg: \"genozip --input sam.xz\")");
                }

                if (str_isprefix_((rom)block, block_size, ZIP_MAGIC, 4)) {
                    if (file->redirected) ABORTINP0 ("Compressing piped-in data in zip format is not currently supported");
                    if (file->is_remote) ABORTINP0 ("The data seems to be in zip format. Please use --input to specify the type (eg: \"genozip --input generic.zip\")");
                    ABORTINP0 ("The data seems to be in zip format. Please use --input to specify the type (eg: \"genozip --input generic.zip\")");
                }

                file->codec = CODEC_NONE;
                buf_add_more (evb, &evb->scratch, (char*)block, block_size, "scratch");
            }

            ASSINP (!file->redirected || file->codec == CODEC_NONE || file->codec == CODEC_BGZF, 
                    "genozip only supports piping in data that is either plain (uncompressed) or compressed in BGZF format (typically with .gz extension) (codec=%s)", 
                    codec_name (file->codec));
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
            file->redirected = true;
            file->codec = CODEC_NONE;
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
            file->redirected = true;
            file->codec = CODEC_NONE;
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
            file->redirected = true;
            file->codec = CODEC_NONE;
            break;
        }

        default:
            ABORT ("%s: invalid filename extension for %s files: %s", global_cmd, dt_name (file->data_type), file->name);
    }

    if (flag.show_codec) 
        iprintf ("%s: source_code=%s\n", file->basename, codec_name (file->source_codec));

    if (!file->file) goto fail;

    file_initialize_txt_file_data (file);

    if (file->is_remote) FREE (error); // allocated by url_get_status
    return file;

fail:
    if (file->is_remote) FREE (error); 
    FREE (file);
    return NULL;
}

FileP file_open_txt_write (rom filename, DataType data_type, int bgzf_level)
{
    ASSERT (data_type > DT_NONE && data_type < NUM_DATATYPES ,"invalid data_type=%u", data_type);

    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype  = TXT_FILE;
    file->mode       = WRITE;
    file->data_type  = data_type;
    file->redirected = !filename;

    file->codec = filename_has_ext(file->name, ".cram") ? CODEC_CRAM
                : filename_has_ext(file->name, ".bcf")  ? CODEC_BCF
                : (bgzf_level != 0)                     ? CODEC_BGZF
                :                                         CODEC_NONE;
    
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
    switch (file->codec) { 
        case CODEC_BGZF : 
        case CODEC_NONE : file->file = file->redirected ? fdopen (STDOUT_FILENO, "wb") : fopen (file->name, WRITE); break;

        case CODEC_CRAM : file_redirect_output_to_stream (file, "samtools", "view", "-OCRAM", file_samtools_no_PG(), NULL); break;
        
        case CODEC_BCF  : {
            char comp_level[4] = { '-', 'l', '0' + MIN_(bgzf_level, 9), 0 };
            file_redirect_output_to_stream (file, "bcftools", "view", "-Ob", "/dev/stdin", comp_level); 
            break;
        }

        default: {} // never reaches here
    }

    ASSINP (file->file, "cannot open file \"%s\": %s", file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    file_initialize_txt_file_data (file);

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

            mutex_initialize (file->wait_for_vb_1_mutex[did_i]);
            mutex_lock (file->wait_for_vb_1_mutex[did_i]); 
        }
        
        // initialize evb "promiscuous" buffers - i.e. buffers that can be allocated by any thread (obviously protected by eg a mutex)
        // promiscuous buffers must be initialized by the main thread, and buffer.c does not verify their integrity.
        Z_INIT (ra_buf);
        Z_INIT (ra_buf_luft);
        Z_INIT (sag_grps);
        Z_INIT (sag_alns);
        Z_INIT (sag_qnames);
        Z_INIT (sag_cigars);
        Z_INIT (sag_seq);
        Z_INIT (sag_qual);
        Z_INIT (deep_index_by[BY_SEQ]);
        Z_INIT (deep_index_by[BY_QNAME]);
        Z_INIT (deep_ents);

        Z_INIT (contexts[CHROM].chrom2ref_map);
        buf_set_shared (&file->contexts[CHROM].chrom2ref_map);
    }
    else {
        Z_INIT (sag_qual);
        Z_INIT (sag_cigars);
    }

    if (flag.biopsy_line.line_i == NO_LINE) // no need to initialize in --biopsy-line (as destroying it later will error)
        serializer_initialize (file->digest_serializer); 

    clock_gettime (CLOCK_REALTIME, &file->start_time);
}

// get time since creation of z_file object in memory
rom file_get_z_run_time (FileP file)
{
    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 

    int seconds_so_far = ((tb.tv_sec - file->start_time.tv_sec)*1000 + 
                         ((int64_t)tb.tv_nsec - (int64_t)file->start_time.tv_nsec) / 1000000) / 1000; 

    static char time_str[16];
    str_human_time (seconds_so_far, true, time_str);

    return time_str;
}

// opens z_file for read or write
FileP file_open_z_read (rom filename)
{
    START_TIMER;

    ASSINP0 (filename, "it is not possible to redirect genozip files from stdin");

    FileP file = (FileP)CALLOC (sizeof(File));

    file->supertype = Z_FILE;
    file->mode      = READ;
    file->is_in_tar = !flag_loading_auxiliary && flag.t_size;
    
    rom disk_filename = file->is_in_tar ? tar_get_tar_name() : filename;
    ASSINP (file_exists (disk_filename), "cannot open \"%s\" for reading: %s", disk_filename, strerror (errno));

    file->disk_size = file->is_in_tar ? flag.t_size : file_get_size (filename);

    file_set_filename (file, filename);

    file->type = file_get_type (file->name);

    file->basename = filename_base (file->name, false, NULL, NULL, 0);

    // if a FASTA file was given as an argument to --reference or --REFERENCE, get the .ref.genozip file,
    // possobily running --make-reference in a separate process if needed
    if (flag.reading_reference && (file_get_data_type (file_get_type (file->name), true) == DT_FASTA) && (file_get_type (file->name) != FASTA_GENOZIP)) 
        disk_filename = ref_fasta_to_ref (file);

    ASSINP (!flag.reading_reference || filename_has_ext (file->name, REF_GENOZIP_), 
            "You specified file \"%s\", however with --reference or --REFERENCE, you must specify a reference file (%s file or FASTA file)\n"
            "Tip: To create a genozip reference file from a FASTA file, use 'genozip --make-reference myfasta.fa'",
            file->name, REF_GENOZIP_);

    ASSINP (!flag.reading_chain || filename_has_ext (file->name, GENOZIP_EXT), 
            "You specified file \"%s\", however with %s, you must specify a genozip chain file (%s extension)\n"
            "Tip: To create a genozip chain file from a chain file, use e.g. 'genozip my-chain-file.chain.gz --reference target-coord-ref.ref.genozip'",
            file->name, (command==ZIP) ? "--chain" : "--show-chain", GENOZIP_EXT);

    if ((!flag.seg_only && !flag.show_bam) || flag_loading_auxiliary) {

        // make sure file is a regular file (not FIFO, directory etc)
        struct stat sb;
        int cause=0, stat_errno=0;
        if (stat (disk_filename, &sb)) {
            cause = 6; // stat failed
            stat_errno = errno;
        }

        if ((sb.st_mode & S_IFMT) != S_IFREG) cause=7; // not regular file

        if (!cause)
            file->file = fopen (disk_filename, READ);

        // verify that this is a genozip file 
        // we read the Magic at the end of the file (as the magic at the beginning may be encrypted)
        uint32_t magic;
        if (cause ||
            (cause = 1 * !file->file) ||
            (cause = 2 * !sb.st_size) ||
            (cause = 3 * !file_seek (file, -(int)sizeof (magic), SEEK_END, SOFT_FAIL)) || 
            (cause = 4 * !fread (&magic, sizeof (magic), 1, file->file)) ||
            (cause = 5 * (BGEN32 (magic) != GENOZIP_MAGIC && !(flag.show_headers && flag.force)))) {
            
            int fail_errno = errno;
            FCLOSE (file->file, disk_filename);

            if (flag.validate == VLD_REPORT_INVALID) flag.validate = VLD_INVALID_FOUND; 

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
FileP file_open_z_write (rom filename, FileMode mode, DataType data_type)
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
    file->data_type = data_type; // note: for READ, data_type is set by file_open_*

    file->basename = filename_base (file->name, false, NULL, NULL, 0);

    ASSINP (filename_has_ext (file->name, GENOZIP_EXT), 
            "file %s must have a " GENOZIP_EXT " extension", file_printname (file));
    
    // set file->type according to the data type, overriding the previous setting - i.e. if the user
    // uses the --output option, he is unrestricted in the choice of a file name
    file->type = file_get_z_ft_by_txt_in_ft (file->data_type, txt_file->type); 

    mutex_initialize (file->dicts_mutex);
    mutex_initialize (file->custom_merge_mutex);

    if (!flag.seg_only && !flag.show_bam) {

        if (flag.force && !file->is_in_tar) 
            unlink (file->name); // delete file if it already exists (needed in weird cases, eg symlink to non-existing file)

        // if we're writing to a tar file, we get the already-openned tar file
        if (file->is_in_tar)
            file->file = tar_open_file (file->name);

        else if (!flag.zip_no_z_file) {
            file->file = fopen (file->name, file->mode);
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
    
    if (chain_is_loaded)
        file->z_flags.has_gencomp = true; // dual-coordinate file

    file->genozip_version = GENOZIP_FILE_FORMAT_VERSION; // to allow the VER macro to operate consistently across ZIP/PIZ
    
    file_initialize_z_file_data (file);

    ASSINP (file->file || flag.zip_no_z_file, 
            "cannot open file \"%s\": %s", file->name, strerror(errno)); // errno will be retrieve even the open() was called through zlib and bzlib 

    COPY_TIMER_EVB (file_open_z);
    return file;
}

// index file is it is a disk file of a type that can be indexed
static void file_index_txt (ConstFileP file)
{
    ASSERTNOTNULL (file);

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

void file_close (FileP *file_p) 
{
    START_TIMER;

    FileP file = *file_p;

    if (z_file && file == z_file && !flag_loading_auxiliary && 
        flag.show_time_comp_i == COMP_ALL && !flag.show_time[0]) // show-time without the optional parameter 
        profiler_add_evb_and_print_report();

    __atomic_store_n (file_p, (FileP)NULL, __ATOMIC_RELAXED); 

    if (!file) return; // nothing to do

    if (file->file && file->supertype == TXT_FILE) {

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

        else if (file->is_remote)
            url_kill_curl ((FILE**)&file->file);

        else
            FCLOSE (file->file, file_printname (file));

        // create an index file using samtools, bcftools etc, if applicable
        if (file->mode == WRITE && flag.index_txt && !flag_loading_auxiliary) 
            file_index_txt (file);
    }

    else if (file->file && file->supertype == Z_FILE) {

        // ZIP note: we need to destory all even if unused, because they were initialized in file_initialize_z_file_data
        if (IS_ZIP)
            for (Did did_i=0; did_i < (IS_ZIP ? MAX_DICTS : file->num_contexts); did_i++) {
                mutex_destroy (file->ctx_mutex[did_i]); 
                mutex_destroy (file->wait_for_vb_1_mutex[did_i]); 
            }

        if (file->is_in_tar && file->mode != READ)
            tar_close_file (&file->file);
        else
            FCLOSE (file->file, file_printname (file));

        serializer_destroy (file->digest_serializer);     
    }

    // free resources if we are NOT near the end of the execution. If we are at the end of the execution
    // it is faster to just let the process die

    if (!flag.let_OS_cleanup_on_exit) {

        if (IS_PIZ && flag.deep && file->supertype == Z_FILE) { // in this case, deep_index and deep_ents are Buffers containing arrays of Buffers
            for_buf (Buffer, buf, file->deep_index) buf_destroy (*buf);
            for_buf (Buffer, buf, file->deep_ents)  buf_destroy (*buf);  
            huffman_destroy (&file->qname_huf);          
        }

        buflist_destroy_file_bufs (file);

        mutex_destroy (file->dicts_mutex);
        mutex_destroy (file->custom_merge_mutex);
        mutex_destroy (file->qname_huf_mutex);
        mutex_destroy (file->recon_plan_mutex[0]);
        
        FREE (file->name);
        FREE (file->basename);
        FREE (file);
    }

    COPY_TIMER_EVB (file_close);
}

void file_write (FileP file, const void *data, unsigned len)
{
    if (!len) return; // nothing to do

    ASSERTNOTNULL (file);
    ASSERTNOTNULL (file->file);
    ASSERTNOTNULL (data);
    
    size_t bytes_written = fwrite (data, 1, len, (FILE *)file->file); // use fwrite - let libc manage write buffers for us

    // if we're streaming our genounzip/genocat/genols output to another process and that process has 
    // ended prematurely then exit quietly. In genozip we display an error because this means the resulting
    // .genozip file will be corrupted
    if (!file->name && command != ZIP && errno == EINVAL) exit (EXIT_DOWNSTREAM_LOST);

    // exit quietly if failed to write to stdout - likely downstream consumer (piped executable or terminal) was closed
    if (bytes_written < len && !file->name) exit (EXIT_DOWNSTREAM_LOST);

    // error if failed to write to file
    ASSERT (bytes_written == len, "wrote only %u of the expected %u bytes to %s: %s", (int)bytes_written, len, file->name, strerror(errno));
}

void file_remove (rom filename, bool fail_quietly)
{
    chmod (filename, S_IRUSR | S_IWUSR); // make sure its +w so we don't get permission denied (ignore errors)

    int ret = remove (filename); 
    ASSERTW (!ret || fail_quietly, "Warning: failed to remove %s: %s", filename, strerror (errno));
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

    if (!flag.is_windows) {
        sprintf (command, "bgzip -@%u -f \"%s\"", global_max_threads, filename);
        ret = system (command);
    }
    
    if ((ret && errno == ENOENT) || flag.is_windows) { // no bgzip - try pigz
        sprintf (command, "pigz -f \"%s\"", filename);
        ret = system (command);
    }

    if (ret && errno == ENOENT) { // no pigz - try gzip
        memcpy (command, "gzip", 4);
        ret = system (command);
    }

    ASSERTW (!ret, "FYI: \"%s\" returned %d. No harm.", command, ret); 

    if (!ret) {
        // special case: rename .bam.gz -> .bam
        if (fn_len >= 4 && !memcmp (&filename[fn_len-4], ".bam", 4)) {
            char gz_filename[fn_len + 10];
            sprintf (gz_filename, "%s.gz", filename);
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
            if (ftello64 ((FILE *)file->file) == offset) return true; // already at the right offset
        }
    }

    int ret = fseeko64 ((FILE *)file->file, offset, whence);

    if (fail_type != HARD_FAIL) {
        if (!flag.to_stdout && fail_type==WARNING_FAIL) {
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

int64_t file_tell_do (FileP file, FailType soft_fail, rom func, unsigned line)
{
    ASSERTNOTNULL (file);
    ASSERTNOTNULL (file->file);

    if (IS_ZIP && file->supertype == TXT_FILE && file->codec == CODEC_GZ)
        return gzconsumed64 ((gzFile)file->file); 
    
    if (IS_ZIP && file->supertype == TXT_FILE && file->codec == CODEC_BZ2)
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

uint64_t file_get_size (rom filename)
{
    struct stat64 st;
    
    int ret = stat64(filename, &st);
    ASSERT (!ret, "stat64 failed on '%s': %s", filename, strerror(errno));
    
    return st.st_size;
}

bool file_is_dir (rom filename)
{
    struct stat64 st;
    int ret = stat64(filename, &st); // 0 if successful
    
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
                    bool verify_textual/*plain ascii*/, bool add_string_terminator)
{
    bool is_stdin = !strcmp (filename, "-");
    if (is_stdin && !max_size) max_size = 10000000; // max size for stdin
    
    uint64_t size = max_size ? max_size : file_get_size (filename);

    buf_alloc (vb, buf, 0, size + add_string_terminator, char, 1, buf_name);

    FILE *file = is_stdin ? stdin : fopen (filename, "rb");
    ASSINP (file, "cannot open \"%s\": %s", filename, strerror (errno));

    buf->len = fread (buf->data, 1, size, file);
    ASSERT (is_stdin || max_size || buf->len == size, "Error reading file %s: %s", filename, strerror (errno));

    ASSINP (!verify_textual || str_is_printable (STRb(*buf)), "Expecting \"%s\" to be a textual file", filename);

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
                    mode_t mode) // optional - ignored if 0
{
    int fn_len = strlen (filename);

    // remove invalid characters from filename
    if (flag.is_windows)
        for (char *c=(char*)filename ; *c ; c++)
            if (*c == ':' && (c-filename != 1)) *c = '-'; // ':' exist eg in SAM AUX names 

    char *tmp_filename = MALLOC (fn_len+5);
    // we first write to tmp_filename, and after we complete and flush, we rename to the final name
    // so that if a file exists (in its final name) - then its guaranteed to be fully written
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

    // write in blocks (Windows hangs if the block is too big, a few GB)
    size_t written = 0;
    const uint64_t block_size = 1 << 24; // 16MB
    for (int i=0; i < (len + block_size - 1) / block_size; i++) // round up
        written += fwrite (&((rom)data)[i * block_size], 1, MIN_(block_size, len - i*block_size), file);
    
    SAVE_VALUE (errno);
    fflush (file);    
    FCLOSE (file, tmp_filename); 
    RESTORE_VALUE (errno); // in cases caller wants to print fwrite error

    if (written != len) {
        WARN ("Failed to write %s: wrote only %"PRIu64" bytes of the expected %"PRIu64, tmp_filename, (uint64_t)written, len);
        put_data_tmp_filenames[my_file_i] = NULL; // no need to lock mutex
        remove (tmp_filename);
        return false;
    }

    // we can't enter if file_put_data_abort is active, or it needs to wait for us before deleting tmp files
    mutex_lock (put_data_mutex);

    remove (filename);
    int renamed_failed = rename (tmp_filename, filename);

    put_data_tmp_filenames[my_file_i] = NULL; // remove tmp file name from list 
    
    mutex_unlock (put_data_mutex);

    ASSERT (!renamed_failed, "Failed to rename %s to %s: %s", tmp_filename, filename, strerror (errno));
    FREE (tmp_filename);

    if (mode) 
        ASSERT (!chmod (filename, mode), "Failed to chmod %s: %s", filename, strerror (errno));
    
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

PutLineFn file_put_line (VBlockP vb, STRp(line), rom msg)
{
    PutLineFn fn;
    sprintf (fn.s, "line.%u.%d.%s%s", vb->vblock_i, vb->line_i, command==ZIP ? "zip" : "piz",
             file_plain_ext_by_dt ((VB_DT(SAM) && z_file->z_flags.txt_is_bin) ? DT_BAM : vb->data_type));
    
    file_put_data (fn.s, STRa(line), 0);

    if (IS_PIZ)
        WARN ("\n%s line=%s line_in_file(1-based)=%"PRIu64". Dumped %s (dumping first occurance only)", 
                msg, line_name(vb).s, writer_get_txt_line_i (vb, vb->line_i), fn.s);
    else
        WARN ("\n%s line=%s vb_size=%u MB. Dumped %s", msg, line_name(vb).s, (int)(segconf.vb_size >> 20), fn.s);

    return fn;
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

rom ft_name (FileType ft)
{
    return type_name (ft, &file_exts[ft], ARRAY_LEN(file_exts));
}

rom file_plain_ext_by_dt (DataType dt)
{
    FileType plain_ft = txt_in_ft_by_dt[dt][0].in;

    return file_exts[plain_ft];
}

bool file_buf_locate (FileP file, ConstBufferP buf)
{
    return is_p_in_range (buf, file, sizeof (File));
}
