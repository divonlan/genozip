// ------------------------------------------------------------------
//   genozip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#ifndef _WIN32
#include <execinfo.h>
#include <signal.h>
#endif
#ifndef _MSC_VER // Microsoft compiler
#include <getopt.h>
#else
#include "compatibility/visual_c_getopt.h"
#endif

#include "genozip.h"
#include "text_help.h"
#include "version.h" // automatically incremented by the make when we create a new distribution
#include "txtfile.h"
#include "zfile.h"
#include "zip.h"
#include "piz.h"
#include "crypt.h"
#include "file.h"
#include "vblock.h"
#include "endianness.h"
#include "regions.h"
#include "vcf.h"
#include "stream.h"
#include "url.h"
#include "strings.h"
#include "stats.h"
#include "arch.h"
#include "license.h"
#include "vcf.h"
#include "dict_id.h"

// globals - set it main() and never change
const char *global_cmd = NULL; 
ExeType exe_type;
int command = -1;  // must be static or global to initialize list_options 

uint32_t global_max_threads = DEFAULT_MAX_THREADS; 
uint32_t global_max_memory_per_vb = 0; // ZIP only: used for reading text file data

// the flags - representing command line options - available globally
int flag_quiet=0, flag_force=0, flag_concat=0, flag_md5=0, flag_split=0, flag_optimize=0, flag_bgzip=0, flag_bam=0, flag_bcf=0,
    flag_show_alleles=0, flag_show_time=0, flag_show_memory=0, flag_show_dict=0, flag_show_gt_nodes=0, flag_multiple_files=0,
    flag_show_b250=0, flag_show_sections=0, flag_show_headers=0, flag_show_index=0, flag_show_gheader=0, flag_show_threads=0,
    flag_stdout=0, flag_replace=0, flag_test=0, flag_regions=0, flag_samples=0, flag_fast=0,
    flag_drop_genotypes=0, flag_no_header=0, flag_header_only=0, flag_header_one=0, flag_noisy=0,
    flag_show_vblocks=0, flag_gtshark=0, flag_sblock=0, flag_vblock=0, flag_gt_only=0, flag_fasta_sequential=0,
    flag_debug_memory=0, flag_debug_progress=0, flag_show_hash, flag_register=0, flag_debug_no_singletons=0,

    flag_optimize_sort=0, flag_optimize_PL=0, flag_optimize_GL=0, flag_optimize_GP=0, flag_optimize_VQSLOD=0, 
    flag_optimize_QUAL=0, flag_optimize_Vf=0, flag_optimize_ZM=0;


uint64_t flag_stdin_size = 0;
char *flag_grep = NULL;

DictIdType dict_id_show_one_b250 = { 0 },  // argument of --show-b250-one
           dict_id_show_one_dict = { 0 },  // argument of --show-dict-one
           dict_id_dump_one_b250 = { 0 };  // argument of --dump-b250-one

static char *threads_str  = NULL;

void exit_on_error(void) 
{
    buf_test_overflows_all_vbs();

    url_kill_curl();
    file_kill_external_compressors(); 

    // if we're in ZIP - remove failed genozip file (but don't remove partial failed text file in PIZ - it might be still useful to the user)
    if (command == ZIP && z_file && z_file->name) {
        char *save_name = malloc (strlen (z_file->name)+1);
        strcpy (save_name, z_file->name);

        // if we're not the main thread - cancel the main thread before closing z_file, so that the main 
        // thread doesn't attempt to access it (eg. z_file->data_type) and get a segmentation fault.
        if (!arch_am_i_io_thread()) 
            cancel_io_thread(); 

        file_close (&z_file, false); // also frees file->name

        file_remove (save_name, true);
    }

    exit(1);
} 

#ifndef _WIN32
static void main_sigsegv_handler (int sig) 
{
#   define STACK_DEPTH 15
    void *array[STACK_DEPTH];
    size_t size = backtrace(array, STACK_DEPTH);
    
    fprintf (stderr, "\nError: segmentation fault. Call stack:\n");
    backtrace_symbols_fd (array, size, STDERR_FILENO);
    exit(1);
}
#endif

static void main_print_help (bool explicit)
{
    static const char **texts[] = {help_genozip, help_genounzip, help_genols, help_genocat, help_genozip_developer}; // same order as ExeType
    static unsigned sizes[] = {sizeof(help_genozip), sizeof(help_genounzip), sizeof(help_genols), sizeof(help_genocat), sizeof(help_genozip_developer)};
    
    if (flag_force) exe_type = (ExeType)4; // -h -f shows developer help

    str_print_text (texts[exe_type], sizes[exe_type] / sizeof(char*), 
                    flag_force ? "                          " : "                     ",  "\n", 0);
    str_print_text (help_footer, sizeof(help_footer) / sizeof(char*), "", "\n", 0);

// in Windows, we ask the user to click a key - this is so that if the user double clicks on the EXE
// from Windows Explorer - the terminal will open and he will see the help
#ifdef _WIN32
    if (!explicit) {
        printf ("Press any key to continue...\n");
        getc(stdin);
    }
#endif
}

static void main_print_version()
{
    printf ("version=%s\n", GENOZIP_CODE_VERSION);  
}

static void main_list_dir(); // forward declaration

static void main_genols (const char *z_filename, bool finalize, const char *subdir, bool recursive) 
{
    if (!finalize) {
        // no specific filename = show entire directory
        if (!z_filename) {
            main_list_dir("."); 
            return;
        }

        // filename is a directory - show directory contents (but not recursively)
        if (!subdir && file_is_dir (z_filename)) {
            main_list_dir (z_filename);
            return;
        }
    }

    static bool first_file = true;
    static unsigned files_listed=0, files_ignored=0;
    static long long total_uncompressed_len=0, total_compressed_len=0;
    char z_size_str[20], txt_size_str[20], num_lines_str[20], indiv_str[20];

    const unsigned FILENAME_WIDTH = 40;

    const char *head_format   = "\n%5s %11s %10s %10s %6s %s  %*s %s\n";
    const char *foot_format_1 = "\nTotal:            %10s %10s %5uX\n";
    const char *foot_format_2 = "\nTotal:            %10s %10s %5.1fX\n";
    const char *item_format_1 = "%5s %11s %10s %10s %5uX %s  %s%s%*s %s\n";
    const char *item_format_2 = "%5s %11s %10s %10s %5.1fX %s  %s%s%*s %s\n";

    // we accumulate the string in str_buf and print in the end - so it doesn't get mixed up with 
    // warning messages regarding individual files
    static Buffer str_buf = EMPTY_BUFFER; 
    
    if (finalize) {
        if (files_listed > 1) {
            str_size(total_compressed_len, z_size_str);
            str_size(total_uncompressed_len, txt_size_str);
            double ratio = total_compressed_len ? ((double)total_uncompressed_len / (double)total_compressed_len) : 0;

            bufprintf (evb, &str_buf, ratio < 100 ? foot_format_2 : foot_format_1, z_size_str, txt_size_str, ratio);
        }
        
        ASSERTW (!files_ignored, "Ignored %u file%s that %s not have a " GENOZIP_EXT " extension", 
                 files_ignored, files_ignored==1 ? "" : "s", files_ignored==1 ? "does" : "do");
        
        goto finish;
    }

    if (first_file) {
        bufprintf (evb, &str_buf, head_format, "Indiv", "Records", "Compressed", "Original", "Factor", " MD5 of original textual file    ", -(int)FILENAME_WIDTH, "Name", "Creation");
        first_file = false;
    }
    
    if (!file_has_ext (z_filename, GENOZIP_EXT) || access (z_filename, F_OK)!=0) {
        files_ignored++;
        goto finish;
    }

    bool is_subdir = subdir && (subdir[0] != '.' || subdir[1] != '\0');

    z_file = file_open (z_filename, READ, Z_FILE, 0); // open global z_file
    uint64_t txt_data_size, num_lines;
    uint32_t num_samples;
    Md5Hash md5_hash_concat, license_hash;
    char created[FILE_METADATA_LEN];
    bool success = zfile_get_genozip_header (&txt_data_size, &num_samples, &num_lines, 
                                             &md5_hash_concat, created, FILE_METADATA_LEN, &license_hash);
    if (!success) goto finish;

    double ratio = z_file->disk_size ? ((double)txt_data_size / (double)z_file->disk_size) : 0;
    
    str_size (z_file->disk_size, z_size_str);
    str_size (txt_data_size, txt_size_str);
    str_uint_commas (num_lines, num_lines_str);
    str_uint_commas (num_samples, indiv_str);
    
    bufprintf (evb, &str_buf, ratio < 100 ? item_format_2 : item_format_1, indiv_str, num_lines_str, 
               z_size_str, txt_size_str, ratio, 
               md5_display (&md5_hash_concat, true),
               (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
               is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH,
               z_filename, created);
            
    total_compressed_len   += z_file->disk_size;
    total_uncompressed_len += txt_data_size;
    
    files_listed++;

    file_close (&z_file, false);

finish:
    if (!recursive) {
            buf_print (&str_buf, false);
            buf_free (&str_buf);
    }
}
static void main_genounzip (const char *z_filename,
                            const char *txt_filename, 
                            unsigned max_threads,
                            bool is_last_file)
{
    txtfile_header_initialize();
    
    // get input FILE
    ASSERT0 (z_filename, "Error: z_filename is NULL");

    // we cannot work with a remote genozip file because the decompression process requires random access
    ASSERT (!url_is_url (z_filename), 
            "%s: genozip files must be regular files, they cannot be a URL: %s", global_cmd, z_filename);

    ASSERT (!txt_filename || !url_is_url (txt_filename), 
            "%s: output files must be regular files, they cannot be a URL: %s", global_cmd, txt_filename);

    unsigned fn_len = strlen (z_filename);

    // skip this file if its size is 0
    RETURNW (file_get_size (z_filename),, "Cannot decompress file %s because its size is 0 - skipping it", z_filename);

    if (!txt_filename && (!flag_stdout || flag_bgzip || flag_bcf || flag_bam) && !flag_split) {
        txt_filename = (char *)malloc(fn_len + 10);
        ASSERT(txt_filename, "Error: failed to malloc txt_filename, len=%u", fn_len+10);

        // .vcf.genozip -> .vcf or .vcf.gz or .bcf ; .sam.genozip -> .sam or .bam or .sam.gz ; fastq.genozip -> .fastq or .fastq.gz
        sprintf ((char *)txt_filename, "%.*s%s", 
                 (int)(fn_len - strlen(GENOZIP_EXT)), z_filename,
                 flag_bgzip ? ".gz" : flag_bam ? ".bam" : flag_bcf ? ".bcf" : "");    
    }

    z_file = file_open (z_filename, READ, Z_FILE, DT_NONE);    

    // get output FILE 
    if (txt_filename) {
        ASSERT0 (!txt_file || flag_concat, "Error: txt_file is open but not in concat mode");

        if (!txt_file)  // in concat mode, for second file onwards, txt_file is already open
            txt_file = file_open (txt_filename, WRITE, TXT_FILE, z_file->data_type);
    }
    else if (flag_stdout) { // stdout
        txt_file = file_open_redirect (WRITE, TXT_FILE, z_file->data_type); // STDOUT
    }
    else if (flag_split) {
        // do nothing - the component files will be opened by txtfile_genozip_to_txt_header()
    }
    else {
        ABORT0 ("Error: unrecognized configuration for the txt_file");
    }
    
    const char *basename = file_basename (z_filename, false, "(stdin)", NULL, 0);
    
    // a loop for decompressing all components in split mode. in non-split mode, it collapses to one a single iteration.
    bool piz_successful;
    unsigned num_components=0;
    do {
        piz_successful = piz_dispatcher (basename, max_threads, num_components==0, is_last_file);
        if (piz_successful) num_components++;
    } while (flag_split && piz_successful); 

    if (!flag_concat && !flag_stdout && !flag_split) 
        // don't close the concatenated file - it will close with the process exits
        // don't close in split mode - piz_dispatcher() opens and closes each component
        // don't close stdout - in concat mode, we might still need it for the next file
        file_close (&txt_file, false); 

    file_close (&z_file, false);

    FREE ((void *)basename);

    if (flag_replace && txt_filename && z_filename) file_remove (z_filename, true); 
}

// run the test genounzip after genozip - for the most reliable testing that is nearly-perfectly indicative of actually 
// genounzipping, we create a new genounzip process
static void main_test_after_genozip (char *exec_name, char *z_filename)
{
    const char *password = crypt_get_password();

    StreamP test = stream_create (0, 0, 0, 0, 0, 0, 
                                  "To use the --test option",
                                  exec_name, "--decompress", "--test", z_filename,
                                  flag_quiet       ? "--quiet"       : SKIP_ARG,
                                  password         ? "--password"    : SKIP_ARG,
                                  password         ? password        : SKIP_ARG,
                                  flag_show_memory ? "--show-memory" : SKIP_ARG,
                                  flag_show_time   ? "--show-time"   : SKIP_ARG,
                                  threads_str      ? "--threads"     : SKIP_ARG,
                                  threads_str      ? threads_str     : SKIP_ARG,
                                  NULL);

    // wait for child process to finish, so that the shell doesn't print its prompt until the test is done
    int exit_code = stream_wait_for_exit (test);
    ASSERT (!exit_code, "genozip test exited with status %d\n", exit_code);
}

static void main_genozip (const char *txt_filename, 
                          char *z_filename,
                          unsigned max_threads,
                          bool is_first_file, bool is_last_file,
                          char *exec_name)
{
    license_get(); // ask the user to register if she doesn't already have a license (note: only genozip requires registration - unzip,cat,ls do not)

    ASSERT (!z_filename || !url_is_url (z_filename), 
            "%s: output files must be regular files, they cannot be a URL: %s", global_cmd, z_filename);

    // get input file
    if (txt_filename) {
        // open the file
        txt_file = file_open (txt_filename, READ, TXT_FILE, 0);

        // skip this file if its size is 0
        RETURNW (txt_file,, "Cannot compresss file %s because its size is 0 - skipping it", txt_filename);

        if (!txt_file->file) return; // this is the case where multiple files are given in the command line, but this one is not compressible - we skip it
    }
    else {  // stdin
        txt_file = file_open_redirect (READ, TXT_FILE, DT_NONE);
        flag_stdout = (z_filename == NULL); // implicit setting of stdout by using stdin, unless -o was used
    }
 
    ASSERT0 (flag_concat || !z_file, "Error: expecting z_file to be NULL in non-concat mode");

    // get output FILE
    if (!flag_stdout) {

        if (!z_file) { // skip if we're the second file onwards in concatenation mode - nothing to do

            if (!z_filename) {
                bool is_url = url_is_url (txt_filename);
                const char *basename = is_url ? file_basename (txt_filename, false, "", 0,0) : NULL;
                const char *local_txt_filename = basename ? basename : txt_filename;

                unsigned fn_len = strlen (local_txt_filename);
                z_filename = (char *)malloc (fn_len + 30); // add enough the genozip extension e.g. 23andme.genozip
                ASSERT(z_filename, "Error: Failed to malloc z_filename len=%u", fn_len+4);

                // if the file has an extension matching its type, replace it with the genozip extension, if not, just add the genozip extension
                const char *genozip_ext = file_exts[file_get_z_ft_by_txt_in_ft (txt_file->data_type, txt_file->type)];

                if (file_has_ext (local_txt_filename, file_exts[txt_file->type]))
                    sprintf (z_filename, "%.*s%s", (int)(fn_len - strlen (file_exts[txt_file->type])), local_txt_filename,
                             genozip_ext); 
                else 
                    sprintf (z_filename, "%s%s", local_txt_filename, genozip_ext); 

                if (basename) FREE ((char*)basename);
            }

            z_file = file_open (z_filename, WRITE, Z_FILE, txt_file->data_type);
        }
    }
    else if (flag_stdout) { // stdout
#ifdef _WIN32
        // this is because Windows redirection is in text (not binary) mode, meaning Windows edits the output stream...
        ASSERT (isatty(1), "%s: redirecting binary output is not supported on Windows, use --output instead", global_cmd);
#endif
        ASSERT (flag_force || !isatty(1), "%s: you must use --force to output a compressed file to the terminal", global_cmd);

        z_file = file_open_redirect (WRITE, Z_FILE, txt_file->data_type);
    } 
    else ABORT0 ("Error: No output channel");
    
    const char *basename = file_basename (txt_filename, false, "(stdin)", NULL, 0);
    zip_dispatcher (basename, max_threads, is_last_file);

    if (flag_show_sections && is_last_file) stats_show_sections();

    bool remove_txt_file = z_file && flag_replace && txt_filename;

    file_close (&txt_file, !is_last_file);

    if ((is_last_file || !flag_concat) && !flag_stdout && z_file) 
        file_close (&z_file, !is_last_file); 

    if (remove_txt_file) file_remove (txt_filename, true); 

    FREE ((void *)basename);

    // test the compression, if the user requested --test
    if (flag_test && (!flag_concat || is_last_file)) main_test_after_genozip (exec_name, z_filename);
}

static void main_list_dir(const char *dirname)
{
    DIR *dir;
    struct dirent *ent;

    dir = opendir (dirname);
    ASSERT (dir, "Error: failed to open directory: %s", strerror (errno));

    int ret = chdir (dirname);
    ASSERT (!ret, "Error: failed to chdir(%s)", dirname);

    while ((ent = readdir(dir))) 
        if (!file_is_dir (ent->d_name))  // don't go down subdirectories recursively
            main_genols (ent->d_name, false, dirname, true);
    
    closedir(dir);    

    ret = chdir ("..");
    ASSERT0 (!ret, "Error: failed to chdir(..)");
}

void main_warn_if_duplicates (int argc, char **argv, const char *out_filename)
{
    int num_files = argc - optind;
    if (num_files <= 1) return; // nothing to do

    # define BASENAME_LEN 256
    char *basenames = malloc (num_files * BASENAME_LEN);

    for (unsigned i=0; i < num_files; i++)
        file_basename (argv[optind + i], false, "", &basenames[i*BASENAME_LEN], BASENAME_LEN);

    qsort (basenames, num_files, BASENAME_LEN, (int (*)(const void *, const void *))strcmp);

    for (unsigned i=1; i < num_files; i++) 
        ASSERTW (strncmp(&basenames[(i-1) * BASENAME_LEN], &basenames[i * BASENAME_LEN], BASENAME_LEN), 
                 "Warning: two files with the same name '%s' - if you later split with 'genounzip --split %s', these files will overwrite each other", 
                 &basenames[i * BASENAME_LEN], out_filename);

    FREE (basenames);    
}

void genozip_set_global_max_memory_per_vb (const char *mem_size_mb_str)
{
    const char *err_msg = "Error: invalid argument of --vblock: %s. Expecting an integer between 1 and 2048. The file will be read and processed in blocks of this number of megabytes.";

    unsigned len = strlen (mem_size_mb_str);
    ASSERT (len <= 4 || (len==1 && mem_size_mb_str[0]=='0'), err_msg, mem_size_mb_str);

    for (unsigned i=0; i < len; i++) {
        ASSERT (IS_DIGIT (mem_size_mb_str[i]), err_msg, mem_size_mb_str);
    }

    unsigned mem_size_mb = atoi (mem_size_mb_str);
    ASSERT (mem_size_mb <= 2048, err_msg, mem_size_mb_str);

    global_max_memory_per_vb = mem_size_mb * 1024 * 1024;
}

//#include "base64.h"
int main (int argc, char **argv)
{
    arch_initialize();

#ifdef _WIN32
    // lowercase argv[0] to allow case-insensitive comparison in Windows
    str_to_lowercase (argv[0]);
#else
    signal (SIGSEGV, main_sigsegv_handler);   // segmentation fault handler
#endif

    if      (strstr (argv[0], "genols"))    exe_type = EXE_GENOLS;
    else if (strstr (argv[0], "genocat"))   exe_type = EXE_GENOCAT;
    else if (strstr (argv[0], "genounzip")) exe_type = EXE_GENOUNZIP;
    else                                    exe_type = EXE_GENOZIP; // default
    
    buf_initialize();

    char *out_filename = NULL;

    global_cmd = file_basename (argv[0], true, "(executable)", NULL, 0); // global var

    bool is_short[256]; // indexed by character of short option.
    memset (is_short, 0, sizeof(is_short));

    // process command line options
    while (1) {

        #define _i  {"input-type",    required_argument, 0, 'i'                }
        #define _I  {"stdin-size",    required_argument, 0, 'I'                }
        #define _c  {"stdout",        no_argument,       &flag_stdout,       1 }
        #define _d  {"decompress",    no_argument,       &command, UNZIP       }
        #define _f  {"force",         no_argument,       &flag_force,        1 }
        #define _h  {"help",          no_argument,       &command, HELP        }
        #define _l  {"list",          no_argument,       &command, LIST        }
        #define _L1 {"license",       no_argument,       &command, LICENSE     } // US spelling
        #define _L2 {"licence",       no_argument,       &command, LICENSE     } // British spelling
        #define _q  {"quiet",         no_argument,       &flag_quiet,        1 }
        #define _Q  {"noisy",         no_argument,       &flag_noisy,        1 }
        #define _DL {"replace",       no_argument,       &flag_replace,      1 }
        #define _V  {"version",       no_argument,       &command, VERSION     }
        #define _z  {"bgzip",         no_argument,       &flag_bgzip,        1 }
        #define _zb {"bam",           no_argument,       &flag_bam,          1 }
        #define _zc {"bcf",           no_argument,       &flag_bcf,          1 }
        #define _m  {"md5",           no_argument,       &flag_md5,          1 }
        #define _t  {"test",          no_argument,       &flag_test,         1 }
        #define _fa {"fast",          no_argument,       &flag_fast,         1 }
        #define _9  {"optimize",      no_argument,       &flag_optimize,     1 } // US spelling
        #define _99 {"optimise",      no_argument,       &flag_optimize,     1 } // British spelling
        #define _9s {"optimize-sort", no_argument,       &flag_optimize_sort,1 }
        #define _9P {"optimize-PL",   no_argument,       &flag_optimize_PL,  1 }
        #define _9G {"optimize-GL",   no_argument,       &flag_optimize_GL,  1 }
        #define _9g {"optimize-GP",   no_argument,       &flag_optimize_GP,  1 }
        #define _9V {"optimize-VQSLOD", no_argument,     &flag_optimize_VQSLOD, 1 }
        #define _9Q {"optimize-QUAL", no_argument,       &flag_optimize_QUAL,1 } 
        #define _9f {"optimize-Vf",   no_argument,       &flag_optimize_Vf,  1 }
        #define _9Z {"optimize-ZM",   no_argument,       &flag_optimize_ZM,  1 }
        #define _gt {"gtshark",       no_argument,       &flag_gtshark,      1 } 
        #define _th {"threads",       required_argument, 0, '@'                }
        #define _O  {"split",         no_argument,       &flag_split,        1 }
        #define _o  {"output",        required_argument, 0, 'o'                }
        #define _p  {"password",      required_argument, 0, 'p'                }
        #define _B  {"vblock",        required_argument, 0, 'B'                }
        #define _S  {"sblock",        required_argument, 0, 'S'                }
        #define _r  {"regions",       required_argument, 0, 'r'                }
        #define _tg {"targets",       required_argument, 0, 't'                }
        #define _s  {"samples",       required_argument, 0, 's'                }
        #define _g  {"grep",          required_argument, 0, 'g'                }
        #define _G  {"drop-genotypes",no_argument,       &flag_drop_genotypes,1}
        #define _H1 {"no-header",     no_argument,       &flag_no_header,    1 }
        #define _H0 {"header-only",   no_argument,       &flag_header_only,  1 }
        #define _1  {"header-one",    no_argument,       &flag_header_one,   1 }
        #define _GT {"GT-only",       no_argument,       &flag_gt_only,      1 }
        #define _Gt {"gt-only",       no_argument,       &flag_gt_only,      1 }
        #define _fs {"sequential",    no_argument,       &flag_fasta_sequential, 1 }  
        #define _rg {"register",      no_argument,       &flag_register,     1 }
        #define _ss {"show-sections", no_argument,       &flag_show_sections,1 } 
        #define _sd {"show-dict",     no_argument,       &flag_show_dict,    1 } 
        #define _d1 {"show-one-dict", required_argument, 0, '3'                }
        #define _d2 {"show-dict-one", required_argument, 0, '3'                }
        #define _sg {"show-gt-nodes", no_argument,       &flag_show_gt_nodes,1 } 
        #define _s2 {"show-b250",     no_argument,       &flag_show_b250,    1 } 
        #define _s5 {"show-one-b250", required_argument, 0, '2'                }
        #define _s6 {"show-b250-one", required_argument, 0, '2'                }
        #define _s7 {"dump-one-b250", required_argument, 0, '5'                }
        #define _s8 {"dump-b250-one", required_argument, 0, '5'                }
        #define _sa {"show-alleles",  no_argument,       &flag_show_alleles, 1 }
        #define _st {"show-time",     no_argument,       &flag_show_time   , 1 } 
        #define _sm {"show-memory",   no_argument,       &flag_show_memory , 1 } 
        #define _sh {"show-headers",  no_argument,       &flag_show_headers, 1 } 
        #define _si {"show-index",    no_argument,       &flag_show_index  , 1 } 
        #define _sr {"show-gheader",  no_argument,       &flag_show_gheader, 1 }  
        #define _sT {"show-threads",  no_argument,       &flag_show_threads, 1 }  
        #define _sv {"show-vblocks",  no_argument,       &flag_show_vblocks, 1 }  
        #define _dm {"debug-memory",  no_argument,       &flag_debug_memory, 1 }  
        #define _dp {"debug-progress",no_argument,       &flag_debug_progress, 1 }  
        #define _ds {"debug-no-singletons",no_argument,  &flag_debug_no_singletons, 1 }  
        #define _dh {"show-hash",    no_argument,        &flag_show_hash,   1 }  
        #define _00 {0, 0, 0, 0                                                }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _i, _I, _c, _d, _f, _h, _l, _L1, _L2, _q, _Q, _t, _DL, _V,               _m, _th, _O, _o, _p,                                          _ss, _sd, _sT, _d1, _d2, _sg, _s2, _s5, _s6, _s7, _s8, _sa, _st, _sm, _sh, _si, _sr, _sv, _B, _S, _dm, _dp, _dh,_ds, _9, _99, _9s, _9P, _9G, _9g, _9V, _9Q, _9f, _9Z, _gt, _fa,          _rg, _00 };
        static Option genounzip_lo[]  = {         _c,     _f, _h,     _L1, _L2, _q, _Q, _t, _DL, _V, _z, _zb, _zc, _m, _th, _O, _o, _p,                                               _sd, _sT, _d1, _d2,      _s2, _s5, _s6,                _st, _sm, _sh, _si, _sr,              _dm, _dp,                                                                                 _00 };
        static Option genocat_lo[]    = {                 _f, _h,     _L1, _L2, _q, _Q,          _V,                   _th,     _o, _p, _r, _tg, _s, _G, _1, _H0, _H1, _Gt, _GT,      _sd, _sT, _d1, _d2,      _s2, _s5, _s6,                _st, _sm, _sh, _si, _sr,              _dm, _dp,                                                                   _fs, _g,      _00 };
        static Option genols_lo[]     = {                 _f, _h,     _L1, _L2, _q,              _V,                                _p,                                                                                                      _st, _sm,                             _dm,                                                                                      _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static const char *short_options[] = { // same order as ExeType
            "i:I:cdfhlLqQt^Vzm@:Oo:p:B:S:9KWF", // genozip
            "czfhLqQt^V@:Oo:p:m",               // genounzip
            "hLVp:qf",                          // genols
            "hLV@:p:qQ1r:t:s:H1Go:fg:"          // genocat
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case UNZIP : case LIST : case LICENSE : 
            case VERSION  : case HELP :
                ASSERT(command<0 || command==c, "%s: can't have both -%c and -%c", global_cmd, command, c); 
                command=c; 
                break;

            case 'i' : file_set_input_type (optarg); break;
            case 'I' : file_set_input_size (optarg); break;
            case 'c' : flag_stdout        = 1      ; break;
            case 'F' : flag_fast          = 1      ; break;
            case 'z' : flag_bgzip         = 1      ; break;
            case 'f' : flag_force         = 1      ; break;
            case '^' : flag_replace       = 1      ; break;
            case 'q' : flag_quiet         = 1      ; break;
            case 'Q' : flag_noisy         = 1      ; break;
            case '9' : flag_optimize      = 1      ; break;
            case 'W' : flag_show_sections = 1      ; break;
            case 'K' : flag_gtshark       = 1      ; break;
            case 't' : if (exe_type != EXE_GENOCAT) { flag_test = 1 ; break; }
                       // fall through for genocat -r
            case 'r' : flag_regions = true; regions_add (optarg); break;
            case 's' : flag_samples = true; vcf_samples_add  (optarg); break;
            case 'm' : flag_md5           = 1      ; break;
            case 'O' : flag_split         = 1      ; break;
            case 'G' : flag_drop_genotypes= 1      ; break;
            case 'H' : flag_no_header     = 1      ; break;
            case '1' : flag_header_one    = 1      ; break;
            case '@' : threads_str  = optarg       ; break;
            case 'o' : out_filename = optarg       ; break;
            case 'g' : flag_grep    = optarg       ; break;
            case '2' : dict_id_show_one_b250 = dict_id_make (optarg, strlen (optarg)); break;
            case '5' : dict_id_dump_one_b250 = dict_id_make (optarg, strlen (optarg)); break;
            case '3' : dict_id_show_one_dict = dict_id_make (optarg, strlen (optarg)); break;
            case 'B' : genozip_set_global_max_memory_per_vb (optarg); 
                       flag_vblock = true;
                       break;
            case 'S' : vcf_zip_set_global_samples_per_block (optarg); 
                       flag_sblock = true;
                       break;
            case 'p' : crypt_set_password (optarg) ; break;

            case 0   : // a long option - already handled; except for 'o' and '@'

                if (long_options[exe_type][option_index].val == 'o') 
                    out_filename = optarg;

                if (long_options[exe_type][option_index].val == 'p') 
                    crypt_set_password (optarg);

                else if (long_options[exe_type][option_index].val == '@') 
                    threads_str = optarg;

                else 
                    ASSERT (long_options[exe_type][option_index].flag != &command || 
                            long_options[exe_type][option_index].val  == command ||
                            command < 0, 
                            "%s: can't have both --%s and -%c", global_cmd, long_options[exe_type][option_index].name, command);
                break; 

            case '?' : // unrecognized option - error message already displayed by libc
            default  :
                fprintf(stderr, "Usage: %s [OPTIONS] filename1 filename2...\nTry %s --help for more information.\n", global_cmd, global_cmd);
                exit(0);  
        }
    }

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        if (exe_type == EXE_GENOLS) command = LIST; // genols can be run without arguments
        
        // genozip with no input filename, no output filename, and no output or input redirection - show help (unless --register)
        // note: in docker stdin is a pipe even if going to a terminal. so we show the help even if
        // coming from a pipe. the user must use "-" to redirect from stdin
        else if (command == -1 && optind == argc && !out_filename && 
                 (isatty(0) || arch_am_i_in_docker()) && isatty(1)) {
            if (flag_register) 
                license_get();
            else
                main_print_help (false);
            return 0;
        }

        else if (exe_type == EXE_GENOUNZIP) command = UNZIP;
        else if (exe_type == EXE_GENOCAT) { command = UNZIP; flag_stdout = !out_filename ; }
        else command = ZIP; // default 
    }

    // check for incompatabilities between flags
    #define OT(l,s) is_short[(int)s[0]] ? "-"s : "--"l

    ASSERT (!flag_stdout      || !out_filename,                     "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("output", "o"));
    ASSERT (!flag_stdout      || !flag_split,                       "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("split", "O"));
    ASSERT (!flag_split       || !out_filename,                     "%s: option %s is incompatable with %s", global_cmd, OT("split",  "O"), OT("output", "o"));
    ASSERT (!flag_stdout      || !flag_replace,                     "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("replace", "^"));

    // if flag_md5 we need to seek back and update the md5 in the txt header section - this is not possible with flag_stdout
    ASSERT (!flag_stdout      || !flag_md5,                         "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("md5", "m"));
    ASSERT (!flag_test        || !out_filename || command != UNZIP, "%s: option %s is incompatable with %s", global_cmd, OT("output", "o"),  OT("test", "t"));
    ASSERT (!flag_test        || !flag_replace || command != UNZIP, "%s: option %s is incompatable with %s", global_cmd, OT("replace", "^"), OT("test", "t"));
    ASSERT (!flag_test        || !flag_stdout  || command != ZIP,   "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("test", "t"));
    ASSERT (!flag_header_only || !flag_no_header,                   "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), "header-only");
    ASSERT (!flag_no_header   || !flag_header_one,                  "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), OT("header-one", "1"));
    ASSERT (!flag_quiet       || !flag_noisy,                       "%s: option %s is incompatable with %s", global_cmd, OT("quiet", "q"), OT("noisy", "Q"));
    ASSERT (!flag_test        || !flag_optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("test", "t"), OT("optimize", "9"));
    ASSERT (!flag_md5         || !flag_optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("md5", "m"), OT("optimize", "9"));
    ASSERT (!flag_samples     || !flag_drop_genotypes,              "%s: option %s is incompatable with %s", global_cmd, OT("samples", "s"), OT("drop-genotypes", "G"));

    if (flag_gtshark) stream_abort_if_cannot_run ("gtshark", "To use the --gtshark option"); 

    // deal with final setting of flag_quiet that suppresses warnings 
    
    // don't show progress or warning when outputing to the termial
    if (flag_stdout && isatty(1)) flag_quiet=true; 
    
    // don't show progress for flags that output throughout the process. no issue with flags that output only in the end
    if (flag_show_dict || flag_show_gt_nodes || flag_show_b250 || flag_show_headers || flag_show_threads ||
        dict_id_show_one_b250.num || dict_id_show_one_dict.num || dict_id_dump_one_b250.num || 
        flag_show_alleles || flag_show_vblocks || (flag_show_index && command==UNZIP))
        flag_quiet=true; // don't show progress

    // override these ^ if user chose to be --noisy
    if (flag_noisy) flag_quiet=false;

    if (flag_test) flag_md5=true; // test requires md5

    // default values, if not overridden by the user
    if (!flag_vblock) genozip_set_global_max_memory_per_vb (flag_fast ? TXT_DATA_PER_VB_FAST : TXT_DATA_PER_VB_DEFAULT); 
    if (!flag_sblock) vcf_zip_set_global_samples_per_block (VCF_SAMPLES_PER_VBLOCK); 

    // if --optimize was selected, all optimizations are turned on
    if (flag_optimize)
        flag_optimize_sort = flag_optimize_PL = flag_optimize_GL = flag_optimize_GP = flag_optimize_VQSLOD = 
        flag_optimize_QUAL = flag_optimize_Vf = flag_optimize_ZM = true;
    
    // if any optimization flag is on, we turn on flag_optimize
    if (flag_optimize_sort || flag_optimize_PL || flag_optimize_GL || flag_optimize_GP || flag_optimize_VQSLOD ||
        flag_optimize_QUAL || flag_optimize_Vf || flag_optimize_ZM)
        flag_optimize = true;

    // if using the -o option - check that we don't have duplicate filenames (even in different directory) as they
    // will overwrite each other if extracted with --split
    if (command == ZIP && out_filename && !flag_quiet) main_warn_if_duplicates (argc, argv, out_filename);

    unsigned num_files = argc - optind;
    flag_multiple_files = (num_files > 1);
     
    flag_concat = (command == ZIP) && (out_filename != NULL) && (num_files > 1);

    ASSERT (num_files <= 1 || flag_concat || !flag_show_sections, "%s: --show-sections can only work on one file at time", global_cmd);

    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (threads_str) {
        int ret = sscanf (threads_str, "%u", &global_max_threads);
        ASSERT (ret == 1 && global_max_threads >= 1, "%s: %s requires an integer value of at least 1", global_cmd, OT("threads", "@"));
    }
    else global_max_threads = arch_get_num_cores();
    
    // take action, depending on the command selected
    if (command == VERSION) { main_print_version();   return 0; }
    if (command == LICENSE) { license_display();      return 0; }
    if (command == HELP)    { main_print_help (true); return 0; }

    for (unsigned file_i=0; file_i < MAX (num_files, 1); file_i++) {

        char *next_input_file = optind < argc ? argv[optind++] : NULL;  // NULL means stdin
        
        if (next_input_file && !strcmp (next_input_file, "-")) next_input_file = NULL; // "-" is stdin too

        ASSERT (next_input_file || command != UNZIP, 
                "%s: filename(s) required (redirecting from stdin is not possible)", global_cmd);

        ASSERTW (next_input_file || !flag_replace, "%s: ignoring %s option", global_cmd, OT("replace", "^")); 
        
        switch (command) {
            case ZIP   : main_genozip (next_input_file, out_filename, global_max_threads, file_i==0, !next_input_file || file_i==num_files-1, argv[0]); 
                         break;
            
            case UNZIP : main_genounzip (next_input_file, out_filename, global_max_threads, file_i==num_files-1); break;
            
            case LIST  : main_genols  (next_input_file, false, NULL, false); break;
            
            default         : ABORT ("%s: unrecognized command %c", global_cmd, command);
        }
    }
            
    // if this is "list", finalize
    if (command == LIST) main_genols (NULL, true, NULL, false);

    return 0;
}
