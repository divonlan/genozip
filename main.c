// ------------------------------------------------------------------
//   main.c
//   Copyright (C) 2019 Divon Lan <vczip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <stdbool.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <getopt.h>
#include <pthread.h>
#include <errno.h>
#include <inttypes.h>
#if __WIN32__
#include <process.h>
#else
#include <sys/ioctl.h>
#include <sys/sysinfo.h>
#include <termios.h>
#endif
#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#include "vczip.h"

// globals - set it main() and never change
const char *global_cmd = NULL; 
unsigned global_max_threads = DEFAULT_MAX_THREADS;
bool global_little_endian;

// the flags are globals 
int flag_stdout=0, flag_force=0, flag_replace=0, flag_quiet=0, flag_gzip=0, flag_concat_mode=0, 
    flag_show_content=0;

int main_print_license()
{
#include "license.h"

#if __WIN32__
    unsigned line_width = 100;
#else
    // in Linux, we can get the actual terminal width
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    unsigned line_width = w.ws_col;
#endif

    unsigned num_lines = sizeof(license) / sizeof(char*);
    
    for (unsigned i=0; i < num_lines; i++)  {
        char *line = license[i];
        unsigned line_len = strlen (line);
        
        // print line with line wraps
        while (line_len > line_width) {
            int c; for (c=line_width-1; c>=0 && line[c] != ' '; c--); // find 
            printf ("%.*s\n", c, line);
            line += c + 1; // skip space too
            line_len -= c + 1;
        }
        printf ("%s\n\n", line);
    }
    return 0;
}

int main_print_help()
{
    printf ("\n");
    printf ("Usage: vczip  [options]... [files]...\n");
    printf ("       vcpiz  [options]... [files]...\n");
    printf ("       vccat [options]... [files]...\n");
    printf ("\n");
    printf ("Compress or uncompress VCF (Variant Call Format) files\n");
    printf ("\n");
    printf ("Actions - use at most one of these actions:\n");
    printf ("   -z --compress     compress a .vcf or a .vcs.gz file into a .vcz file. The source file is left unchanged.\n");
    printf ("                     this is the default action for vczip\n");
    printf ("   -d --decompress   decompress a .vcz file into a .vcf file. The .vcz file is left unchanged.\n");
    printf ("                     this is the default action for vcpiz\n");
    printf ("   -l --list         list the compression ratios of these .vcz files\n");
    printf ("   -t --test         test vczip. Compress the .vcf file(s), decompress, and then compare the\n");
    printf ("                     result to the original .vcf - all in memory without writing to any file\n");
    printf ("   -h --help         show this help page\n");
    printf ("   -L --license      show the license terms and conditions for this product\n");
    printf ("   -V --version      display version number\n");
    printf ("\n");    
    printf ("Flags:\n");    
    printf ("   -c --stdout       send output to standard output instead of a file. -dc is the default action of dvcat\n");    
    printf ("   -f --force        force overwrite of the output file, or force writing .vcz data to standard output\n");    
    printf ("   -R --replace      replace the source file with the result file, rather than leaving it unchanged\n");    
    printf ("   -o --output       output file name. this option can also be used to concatenate multiple input files\n");
    printf ("                     into a single concatented output file\n");
    printf ("   -g --gzip         in conjunction with -d, generates .vcf.gz files compatabile with bgzip / gzip\n");    
    printf ("   -q --quiet        don't show the progress indicator\n");    
    printf ("   -@ --threads      specify how many threads to use. default is 8. for best performance, this should match\n");    
    printf ("   --show-content    show the information content of VCF files and the compression ratios of each components\n");
    printf ("                     the number of logical CPU cores (which is double the number of physical cores on Intel processors)");

    printf ("\n");
    printf ("One or more file names may be given, or if omitted, standard input/output is used instead\n");
    printf ("\n");
    printf ("vczip is freely available for academic use. Commercial use requires a license. vczip contains pantent-pending inventions\n");
    printf ("\n");
    printf ("For bug reports or license inquires: vczip@blackpawventures.com\n");
    printf ("\n");
    printf ("THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE\n");
    printf ("WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE\n");
    printf ("COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT,\n");
    printf ("TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n");
    printf ("\n");
    return 0;
}

int main_print_version()
{
    printf ("version=%u\n", VCZIP_VERSION);
    return 0;
}

// get basename of a filename 
static char *get_basename(const char *fn)
{
    if (!fn) return NULL;

    // we need to make a copy of the name first, as basename() modifies the string
    unsigned fn_len = strlen(fn);

    char *fn_copy = malloc (fn_len+1);
    ASSERT0 (fn_copy, "Failed to malloc fn_copy");

    strcpy (fn_copy, fn);

    return basename(fn_copy);

    // note: we leak fn_copy and never free it - it's ok, its a small amount
}

static unsigned main_get_num_cores()
{
#if __WIN32__
    char *env = getenv ("NUMBER_OF_PROCESSORS");
    if (!env) return DEFAULT_MAX_THREADS;

    unsigned num_cores;
    int ret = sscanf (env, "%u", &num_cores);
    return ret==1 ? num_cores : DEFAULT_MAX_THREADS; 

#elif defined __APPLE__
    int num_cores;
    size_t len = sizeof(num_cores);
    int mib[2] = { CTL_HW, HW_NCPU };
    if (sysctl(mib, 2, &num_cores, &len, NULL, 0))
        return DEFAULT_MAX_THREADS;

    return num_cores;

#else // Linux etc
    return get_nprocs();
#endif
}

static void main_display_section_stats (const File *vcf_file, const File *z_file)
{
    const SectionType secs[] = {SEC_VCF_HEADER, SEC_VARIANT_DATA, SEC_HAPLOTYPE_DATA, SEC_GENOTYPE_DATA, SEC_PHASE_DATA, SEC_STATS_HT_SEPERATOR};
    const char *sec_names[] = {"VCF header", "Columns 1-9", "Allele values", "Other subfields", "Phasing char","Allele separator char"};
    char vsize[30], zsize[30];
    uint64_t total_vcf=0, total_z=0;

    printf ("\n\nSection stats:\n");
    printf ("Section                 VCF bytes     %%     VCZ bytes     %%  Ratio\n");
    const char *format = "%21s    %8s %5.1f      %8s %5.1f  %5.1f%s\n";
    for (unsigned sec_i=0; sec_i < NUM_SEC_TYPES; sec_i++) {
        printf (format, sec_names[sec_i], 
                buf_human_readable_size(vcf_file->section_bytes[secs[sec_i]], vsize), 100.0 * (double)vcf_file->section_bytes[secs[sec_i]] / (double)vcf_file->vcf_data_size,
                buf_human_readable_size(z_file->section_bytes[secs[sec_i]], zsize), 100.0 * (double)z_file->section_bytes[secs[sec_i]] / (double)z_file->disk_size,
                z_file->section_bytes[secs[sec_i]] ? (double)vcf_file->section_bytes[secs[sec_i]] / (double)z_file->section_bytes[secs[sec_i]] : 0,
                !z_file->section_bytes[secs[sec_i]] ? (vcf_file->section_bytes[secs[sec_i]] ? "\b\b\bInf" : "\b\b\b---") : "");

        total_vcf += vcf_file->section_bytes[sec_i];
        total_z   += z_file->section_bytes[sec_i];
    }
    printf (format, "TOTAL", 
            buf_human_readable_size(total_vcf, vsize), 100.0 * (double)total_vcf / (double)vcf_file->vcf_data_size,
            buf_human_readable_size(total_z, zsize),   100.0 * (double)total_z   / (double)z_file->disk_size,
            (double)total_vcf / (double)total_z, "");

    ASSERTW (total_z == z_file->disk_size, "Hmm... incorrect calculation for VCZ sizes: total section sizes=%"PRIu64" but file size is %"PRIu64" (diff=%"PRId64")", 
             total_z, z_file->disk_size, z_file->disk_size - total_z);

    ASSERTW (total_vcf == vcf_file->vcf_data_size, "Hmm... incorrect calculation for VCF sizes: total section sizes=%"PRIu64" but file size is %"PRIu64" (diff=%"PRId64")", 
             total_vcf, vcf_file->vcf_data_size, vcf_file->vcf_data_size - total_vcf);
}

static void main_vczip (const char *vcf_filename, 
                        char *z_filename,
                        int pipefd_zip_to_unzip,  // send output to pipe (used for testing)
                        unsigned max_threads,
                        bool is_last_file)
{
    File *vcf_file;
    static File *z_file = NULL; // static to support concat mode

    // get input file
    if (vcf_filename) 
        vcf_file = file_open (vcf_filename, READ, VCF);
    else {  // stdin
        vcf_file = file_fdopen (0, READ, STDIN);
        flag_stdout = true; // implicit setting of stdout by using stdin
    }

    // get output FILE
    if (!flag_stdout && pipefd_zip_to_unzip < 0) {

        if (!z_file) { // if we're the second file onwards in concatenation mode - nothing to do
            if (!z_filename) {
                unsigned fn_len = strlen (vcf_filename);
                z_filename = malloc (fn_len + 4);
                ASSERT(z_filename, "Error: Failed to malloc z_filename len=%u", fn_len+4);

                if (vcf_file->type == VCF)
                    sprintf (z_filename, "%.*sz", fn_len-1, vcf_filename); // .vcf -> .vcz
                else // VCF_GZ
                    sprintf (z_filename, "%.*sz", fn_len-4, vcf_filename); // .vcf.gz -> .vcz
            }
            ASSERT (flag_force || access(z_filename, F_OK), "%s: output file %s already exists", global_cmd, z_filename);

            z_file = file_open (z_filename, WRITE, VCZ);
        }
        z_file->disk_at_beginning_of_this_vcf_file = z_file->disk_so_far;
        z_file->vcf_data_so_far                    = 0; // reset these as they relate to the VCF data of the VCF file currently being processed
        z_file->vcf_data_size                      = vcf_file->vcf_data_size; // 0 if vcf is .gz or a pipe or stdin
    }
    else if (flag_stdout || !vcf_filename) { // stdout
#ifdef _WIN32
        ASSERT (isatty(1), "%s: redirecting output is not currently supported on Windows", global_cmd);
#endif
        ASSERT (flag_force || !isatty(1), "%s: you must use --force to output a compressed file to the terminal", global_cmd);

        z_file = file_fdopen (1, WRITE, STDOUT);
    } 
    else if (pipefd_zip_to_unzip >= 0) {
        z_file = file_fdopen (pipefd_zip_to_unzip, WRITE, PIPE);
    }
    else ABORT0 ("Error: No output channel");
    
    zip_dispatcher (flag_quiet ? NULL : get_basename(vcf_filename), vcf_file, z_file, 
                    flag_concat_mode, pipefd_zip_to_unzip >= 0, max_threads);

    if (flag_show_content) main_display_section_stats (vcf_file, z_file);

    bool remove_vcf_file = z_file && flag_replace && vcf_filename;

    file_close (&vcf_file);

    if (is_last_file && z_file) {
        // update the vcf data size in the VCZ vcf_header section
        if (z_file->type == VCZ)
            zfile_update_vcf_header_section_header (z_file, z_file->vcf_concat_data_size, z_file->num_lines);

        file_close (&z_file); 
    }

    if (remove_vcf_file) file_remove (vcf_filename); 
}

static void main_vcpiz (const char *z_filename,
                        char *vcf_filename, 
                        int pipe_from_zip_thread, 
                        int pipe_to_test_thread,
                        unsigned max_threads)
{
    static File *vcf_file = NULL; // static to support concat mode
    File *z_file;

    if (vcf_filename) {
        if (file_has_ext (vcf_filename, ".gz"))
            flag_gzip = true; // output file .gz implicitly turns on flag_gzip

        ASSERT (!flag_gzip || file_has_ext (vcf_filename, ".gz"), "%s: output file %s does not end with gz, which is incompatable with the --gzip / -g option", 
                global_cmd, vcf_filename);
    }

    // get input FILE
    if (z_filename) {
        unsigned fn_len = strlen (z_filename);

        ASSERT (file_has_ext (z_filename, ".vcz"), "%s: file: %s - vczip can only decompress files with a .vcz extension", global_cmd, z_filename);

        if (!vcf_filename && !flag_stdout && pipe_to_test_thread < 0) {
            vcf_filename = malloc(fn_len + 10);
            ASSERT(vcf_filename, "Error: failed to malloc vcf_filename, len=%u", fn_len+10);

            if (flag_gzip)
                sprintf (vcf_filename, "%.*sf.gz", fn_len-1, z_filename);
            else
                sprintf (vcf_filename, "%.*sf", fn_len-1, z_filename);
        }

        z_file = file_open (z_filename, READ, VCZ);    
    }
    else if (pipe_from_zip_thread >= 0) {
        z_file = file_fdopen (pipe_from_zip_thread, READ, VCZ);
    }
    else { // stdin
        z_file = file_fdopen (0, READ, STDIN);
    }

    // get output FILE
    if (vcf_filename) {

        if (!vcf_file) { // in concat mode, for second file onwards, vcf_file is already open
            ASSERT (flag_force || access(vcf_filename, F_OK), "%s: output file %s already exists", global_cmd, vcf_filename);

            vcf_file = file_open (vcf_filename, WRITE, VCF);
        }
    }
    else if (pipe_to_test_thread >= 0) {
        vcf_file = file_fdopen (pipe_to_test_thread, WRITE, VCF);
    }
    else { // stdout
#ifdef _WIN32
        ASSERT (isatty(1), "%s: redirecting output is not currently supported on Windows", global_cmd);
#endif
        vcf_file = file_fdopen (1, WRITE, VCF); // STDOUT
    }

    piz_dispatcher (flag_quiet ? NULL : get_basename (z_filename), z_file, vcf_file, 
                    flag_concat_mode, pipe_from_zip_thread >= 0, max_threads);

    if (!flag_concat_mode) 
        file_close(&vcf_file); // no worries about not closing the concatenated file - it will close with the process exits

    file_close (&z_file);

    if (flag_replace && vcf_filename && z_filename) file_remove (z_filename); 
}

typedef struct {
    const char *vcf_filename;
    int pipe_to_uncompress_thread;
    int max_threads;
} TestToCompressData;

typedef struct {
    int pipe_from_zip_thread;
    int pipe_to_test_thread;
    int max_threads;
} TestToUncompressData;

void *main_test_compress_thread_entry (void *p_)
{
    TestToCompressData *t2c = (TestToCompressData*)p_;

    main_vczip (t2c->vcf_filename, NULL, t2c->pipe_to_uncompress_thread, t2c->max_threads, true);

    return NULL;
}

void *main_test_uncompress_thread_entry (void *p_)
{
    TestToUncompressData *t2u = (TestToUncompressData*)p_;

    main_vcpiz (NULL, NULL, t2u->pipe_from_zip_thread, t2u->pipe_to_test_thread, t2u->max_threads);

    return NULL;
}

static void main_test (const char *vcf_filename)
{
    ASSERT0 (vcf_filename, "vczip: filename missing");

    File *vcf_file = file_open (vcf_filename, READ, VCF);

    int pipefd_zip_to_unzip [2];
    int pipefd_unzip_to_main[2];
#if __WIN32__
    _pipe (pipefd_zip_to_unzip,  25000000, _O_BINARY); // 25MB pipe space to make sure both sides stay busy
    _pipe (pipefd_unzip_to_main, 25000000, _O_BINARY);
#else
    pipe (pipefd_zip_to_unzip);
    pipe (pipefd_unzip_to_main);
#endif

    // open 2 threads - one for compression and for decompression - and the main thread and compares to the original
    // they are connected by two pipes compress thread | decompress thread | main thread

    // thread allocation: 1 thread for this control thread, and the remaining divided between compress and uncompress

    pthread_t test_compress_thread, test_uncompress_thread;
    
    TestToCompressData t2c;
    t2c.vcf_filename              = vcf_filename;
    t2c.pipe_to_uncompress_thread = pipefd_zip_to_unzip[1];
    t2c.max_threads               = (global_max_threads-1) / 2; 

    TestToUncompressData t2u;
    t2u.pipe_from_zip_thread = pipefd_zip_to_unzip[0];
    t2u.pipe_to_test_thread       = pipefd_unzip_to_main[1];
    t2u.max_threads               = (global_max_threads-1) / 2;

    unsigned err = pthread_create(&test_compress_thread, NULL, main_test_compress_thread_entry, &t2c);
    ASSERT (!err, "Error: failed to create test compress thread, err=%u", err);

    err = pthread_create(&test_uncompress_thread, NULL, main_test_uncompress_thread_entry, &t2u);
    ASSERT (!err, "Error: failed to create test decompress thread, err=%u", err);

    FILE *from_pipe = fdopen (pipefd_unzip_to_main[0], "rb");
    ASSERT0 (from_pipe, "Failed to fdopen pipe from decompress");

    vcffile_compare_pipe_to_file (from_pipe, vcf_file);

    fclose (from_pipe);

    file_close(&vcf_file);

    return;
}

static void main_list (const char *z_filename, bool finalize) 
{
    ASSERT (z_filename || finalize, "%s: missing filename", global_cmd);

    static bool first_file = true;
    static unsigned files_listed=0, files_ignored=0;
    static long long total_uncompressed_len=0, total_compressed_len=0;
    char c_str[20], u_str[20];

    if (finalize) {
        if (files_listed > 1) {
            buf_human_readable_size(total_compressed_len, c_str);
            buf_human_readable_size(total_uncompressed_len, u_str);
            unsigned ratio = total_compressed_len ? ((double)total_uncompressed_len / (double)total_compressed_len) : 0;

            printf ("Total:                   %19s %19s  %5uX  \n", c_str, u_str, ratio);
            
            if (files_ignored) printf ("\nIgnored %u files that do not have a .vcz extension or are not readable\n", files_ignored);
        }
        return;
    }

    if (first_file) {
        printf ("Individuals     Variants          Compressed        Uncompressed  Factor  Name\n"); // follow gzip format
        first_file = false;
    }
    
    File *z_file = file_open(z_filename, READ, VCZ_TEST);    
    if (!z_file) {
        files_ignored++;
        return;
    }

    SectionHeaderVCFHeader vcf_header_header;
    bool success = vcf_header_get_vcf_header (z_file, &vcf_header_header);
    if (!success) {
        files_ignored++;
        return;
    }   

    unsigned ratio = z_file->disk_size ? ((double)ENDN64(vcf_header_header.vcf_data_size) / (double)z_file->disk_size) : 0;
    
    buf_human_readable_size(z_file->disk_size, c_str);
    buf_human_readable_size(ENDN64(vcf_header_header.vcf_data_size), u_str);
    printf ("%11u  %11"PRIu64" %19s %19s  %5uX  %s\n", 
            ENDN32(vcf_header_header.num_samples), ENDN64(vcf_header_header.num_lines), 
            c_str, u_str, ratio, z_filename);
            
    total_compressed_len   += z_file->disk_size;
    total_uncompressed_len += ENDN64(vcf_header_header.vcf_data_size);
    
    file_close (&z_file);
}

int main (int argc, char **argv)
{
#   define COMPRESS   'z'
#   define UNCOMPRESS 'd'
#   define TEST       't'
#   define LIST       'l'
#   define LICENSE    'L'
#   define VERSION    'V'
#   define HELP       'h'

    static int command = -1;  // must be static to initialize list_options 
    char *out_filename = NULL;
    char *threads_str = NULL;

    global_cmd = get_basename(argv[0]); // global var

    // verify CPU architecture and compiler is supported
    ASSERT0 (sizeof(char)==1 && sizeof(short)==2 && sizeof (unsigned)==4 && sizeof(long long)==8, 
             "Error: Unsupported C type lengths, check compiler options");
    
    // runtime endianity (byte order) detection - safer than compile-time
    unsigned test_endianity = 0x01020304;
    global_little_endian = *(uint8_t*)&test_endianity==0x04; // in big endian it is 0x01;

    // process command line options
    while (1) {

        static struct option long_options[] =
        {
            // options that set a flag
            {"stdout",     no_argument,       &flag_stdout,  1      },
            {"decompress", no_argument,       &command, UNCOMPRESS  },
            {"force",      no_argument,       &flag_force,   1      },
            {"help",       no_argument,       &command, HELP        },
            {"list",       required_argument, &command, LIST        },
            {"license",    no_argument,       &command, LICENSE     },
            {"quiet",      no_argument,       &flag_quiet, 1        },
            {"replace",    no_argument,       &flag_replace, 1      },
            {"test",       no_argument,       &command, TEST        },
            {"version",    no_argument,       &command, VERSION     },
            {"compress",   no_argument,       &command, COMPRESS    },
            {"threads",    required_argument, 0, '@'                },
            {"gzip",       no_argument,       &flag_gzip, 1         },
            {"output",     required_argument, 0, 'o'                }, 
            {"show-content",no_argument,      &flag_show_content, 1 }, 
            {0, 0, 0, 0                                             },
        };        
        
        int option_index = 0;
        int c = getopt_long (argc, argv, "cdfhlLqRtVz@:go:", long_options, &option_index);

        if (c == -1) break; // no more options

        switch (c) {
            case COMPRESS : case UNCOMPRESS : case LIST : case LICENSE : 
            case TEST     : case VERSION    : case HELP :
                ASSERT(command<0 || command==c, "%s: can't have both -%c and -%c", global_cmd, command, c); 
                command=c; 
                break;

            case 'c' : flag_stdout        = 1 ; break;
            case 'f' : flag_force         = 1 ; break;
            case 'R' : flag_replace       = 1 ; break;
            case 'q' : flag_quiet         = 1 ; break;
            case 'C' : flag_show_content   = 1 ; break;
            case 'g' : flag_gzip          = 1 ; break;
            case '@' : threads_str  = optarg  ; break;
            case 'o' : out_filename = optarg  ; break;

            case 0   : // a long option - already handled; except for 'o' and '@'
                if (long_options[option_index].val == 'o') 
                    out_filename = optarg;

                else if (long_options[option_index].val == '@') 
                    threads_str = optarg;

                else 
                    ASSERT (long_options[option_index].flag != &command || 
                            long_options[option_index].val  == command ||
                            command < 0, 
                            "%s: can't have both --%s and -%c", global_cmd, long_options[option_index].name, command);
                break; 

            case '?' : // unrecognized option - error message already displayed by libc
            default  :
                fprintf(stderr, "Usage: %s [OPTIONS] filename1 filename2...\nTry %s --help for more information.\n", global_cmd, global_cmd);
                abort();  
        }
    }

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        // vczip with no input filename, no output filename, and no output or input redirection - show help
        if (command == -1 && optind == argc && !out_filename && isatty(0) && isatty(1)) command = HELP;
        else if (strstr (argv[0], "vcpiz"))   command = UNCOMPRESS;
        else if (strstr (argv[0], "vccat")) { command = UNCOMPRESS; flag_stdout=1 ; }
        else                                  command = COMPRESS; // default 
    }

    // sanity checks
    ASSERT (!flag_stdout || !out_filename, "%s: option --stdout / -c is incompatable with --output / -o", global_cmd);
    ASSERT (!flag_stdout || !flag_replace, "%s: option --stdout / -c is incompatable with --replace / -R", global_cmd);
    ASSERTW (!flag_stdout       || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring --stdout / -c option", global_cmd);
    ASSERTW (!flag_force        || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring --force / -f option", global_cmd);
    ASSERTW (!flag_replace      || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring --replace / -R option", global_cmd);
    ASSERTW (!flag_quiet        || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring --quiet / -q option", global_cmd);
    ASSERTW (!threads_str       || command == COMPRESS || command == UNCOMPRESS || command == TEST, "%s: ignoring --threads / -@ option", global_cmd);
    ASSERTW (!flag_gzip         ||                        command == UNCOMPRESS, "%s: ignoring --gzip / -g option", global_cmd);
    ASSERTW (!out_filename      || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring --output / -o option", global_cmd);
    ASSERTW (!flag_show_content || command == COMPRESS || command == TEST      , "%s: ignoring --show-content, it only works with -z or -t", global_cmd);
    
    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (threads_str) {
        int ret = sscanf (threads_str, "%u", &global_max_threads);
        ASSERT (ret == 1 && global_max_threads >= 1, "%s: --threads / -@ option requires an integer value of at least 1", global_cmd);
        ASSERT (command != TEST || global_max_threads >= 3, "%s invalid --threads / -@ value: number of threads for --test / -t should be at least 3", global_cmd);
    }
    else    
        global_max_threads = main_get_num_cores();
printf ("num_cores=%u\n", global_max_threads);
    if (command == TEST) {
        flag_stdout = flag_force = flag_replace = flag_quiet = flag_gzip = false;
        out_filename = NULL;
    }

    // take action, depending on the command selected
    if (command == VERSION) return main_print_version();
    if (command == LICENSE) return main_print_license();
    if (command == HELP)    return main_print_help();

    flag_concat_mode = (out_filename != NULL);

    unsigned count=0;
    do {
        char *next_input_file = optind < argc ? argv[optind++] : NULL;  // NULL means stdin

        ASSERTW (next_input_file || !flag_replace, "%s: ignoring --replace / -R option", global_cmd);
        
        ASSERT0 (!count || !flag_show_content, "Error: --showcontent can only work on one file at time");

        switch (command) {
            case COMPRESS   : main_vczip (next_input_file, out_filename, -1, global_max_threads, optind==argc); break;
            case UNCOMPRESS : main_vcpiz (next_input_file, out_filename, -1, -1, global_max_threads); break;
            case TEST       : main_test  (next_input_file); break;
            case LIST       : main_list  (next_input_file, false); break;
            
            default         : ASSERT(false, "%s: unrecognized command %c", global_cmd, command);
        }

        count++;
    } while (optind < argc);
            
    // if this is "list", finalize
    if (command == LIST) main_list (NULL, true);
        
    return 0;
}
