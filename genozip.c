// ------------------------------------------------------------------
//   genozip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <sys/types.h>
#ifndef _MSC_VER // Microsoft compiler
#include <getopt.h>
#else
#include "compatability/visual_c_getopt.h"
#endif
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <process.h>
#elif defined __APPLE__
#include <sys/ioctl.h>
#include <sys/sysctl.h>
#include <termios.h>
#else // LINUX
#include <sched.h>
#include <sys/ioctl.h>
#include <sys/sysinfo.h>
#include <termios.h>
#endif

#include "genozip.h"
#include "text_license.h"
#include "text_help.h"

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT } ExeType;

// globals - set it main() and never change
const char *global_cmd = NULL; 

unsigned global_max_threads = DEFAULT_MAX_THREADS;
bool global_little_endian;

// the flags - some are static, some are globals 
static int flag_stdout=0, flag_force=0, flag_replace=0, flag_show_content=0;
int flag_quiet=0, flag_concat_mode=0, flag_md5=0, flag_split=0, flag_show_alleles=0, flag_show_time=0, flag_show_memory=0;

static int main_print (const char **text, unsigned num_lines,
                       const char *wrapped_line_prefix, 
                       const char *newline_separator, 
                       unsigned line_width /* 0=calcuate optimal */)
{                       
    if (!line_width) {
#ifdef _WIN32
        line_width = 120; // default width of cmd window in Windows 10
#else
        // in Linux and Mac, we can get the actual terminal width
        struct winsize w;
        ioctl(0, TIOCGWINSZ, &w);
        line_width = w.ws_col;
#endif
    }

    for (unsigned i=0; i < num_lines; i++)  {
        const char *line = text[i];
        unsigned line_len = strlen (line);
        
        // print line with line wraps
        bool wrapped = false;
        while (line_len + (wrapped ? strlen (wrapped_line_prefix) : 0) > line_width) {
            int c; for (c=line_width-1 - (wrapped ? strlen (wrapped_line_prefix) : 0); c>=0 && line[c] != ' '; c--); // find 
            printf ("%s%.*s\n", wrapped ? wrapped_line_prefix : "", c, line);
            line += c + 1; // skip space too
            line_len -= c + 1;
            wrapped = true;
        }
        printf ("%s%s%s", wrapped ? wrapped_line_prefix : "", line, newline_separator);
    }
    return 0;
}

static int main_print_help (ExeType exe_type, bool explicit)
{
    static const char **texts[] = {help_genozip, help_genounzip, help_genols, help_genocat}; // same order as ExeType
    static unsigned sizes[] = {sizeof(help_genozip), sizeof(help_genounzip), sizeof(help_genols), sizeof(help_genocat)};
    
    main_print (texts[exe_type], sizes[exe_type] / sizeof(char*), "                     ",  "\n", 0);
    main_print (help_footer, sizeof(help_footer) / sizeof(char*), "", "\n", 0);

// in Windows, we ask the user to click a key - this is so that if the user double clicks on the EXE
// from Windows Explorer - the terminal will open and he will see the help
#ifdef _WIN32
    if (!explicit) {
        printf ("Press any key to continue...\n");
        getc(stdin);
    }
#endif

    return 0;
}

int main_print_version()
{
    printf ("version=%s\n", GENOZIP_CODE_VERSION);  
    return 0;
}

static unsigned main_get_num_cores()
{
#ifdef _WIN32
    char *env = getenv ("NUMBER_OF_PROCESSORS");
    if (!env) return DEFAULT_MAX_THREADS;

    unsigned num_cores;
    int ret = sscanf (env, "%u", &num_cores);
    return ret==1 ? num_cores : DEFAULT_MAX_THREADS; 

#elif defined __APPLE__
    int num_cores;
    size_t len = sizeof(num_cores);
    if (sysctlbyname("hw.activecpu", &num_cores, &len, NULL, 0) &&  
        sysctlbyname("hw.ncpu", &num_cores, &len, NULL, 0))
            return DEFAULT_MAX_THREADS; // if both failed

    return (unsigned)num_cores;
 
#else // Linux etc
    // this works correctly with slurm too (get_nprocs doesn't account for slurm core allocation)
    cpu_set_t cpu_set_mask;
    extern int sched_getaffinity (__pid_t __pid, size_t __cpusetsize, cpu_set_t *__cpuset);
    sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);
    unsigned cpu_count = __sched_cpucount (sizeof (cpu_set_t), &cpu_set_mask);
    // TODO - sort out include files so we don't need this extern

    // if failed to get a number - fall back on good ol' get_nprocs
    if (!cpu_count) cpu_count = get_nprocs();

    return cpu_count;
#endif
}

static void main_display_section_stats (const File *vcf_file, const File *z_file)
{
    char vsize[30], zsize[30];
    uint64_t total_vcf=0, total_z=0;

    fprintf (stderr, "\n\n");
    if (vcf_file->name) fprintf (stderr, "File name: %s\n", vcf_file->name);
    fprintf (stderr, 
#ifdef _MSC_VER
             "Individuals: %u   Variants: %I64u   Non-GT subfields: %u\n", 
#else
             "Individuals: %u   Variants: %"PRIu64"   Non-GT subfields: %u\n", 
#endif
             global_num_samples, z_file->num_lines, z_file->num_subfields);

    fprintf (stderr, "Compression stats:\n");
    fprintf (stderr, "                               VCF     %%       GENOZIP     %%  Ratio\n");
    const char *format = "%22s    %8s %5.1f      %8s %5.1f  %5.1f%s\n";

#ifdef DEBUG
    // the order in which we want them displayed
    const SectionType secs[] = {SEC_VCF_HEADER, SEC_VARIANT_DATA, 
                                SEC_HAPLOTYPE_DATA, SEC_PHASE_DATA, SEC_STATS_HT_SEPERATOR, 
                                SEC_GENOTYPE_DATA, SEC_DICTIONARY};

    const char *categories[] = {"VCF header", "Columns 1-9", 
                                "Haplotype data", "Phasing char", "HT separator char", 
                                "Other subfields", "Other subfields: dicts"};

    uint64_t vbytes[NUM_SEC_TYPES], zbytes[NUM_SEC_TYPES];
    for (unsigned sec_i=0; sec_i < NUM_SEC_TYPES; sec_i++) {
        vbytes[sec_i] = vcf_file->section_bytes[secs[sec_i]];
        zbytes[sec_i] = z_file->section_bytes[secs[sec_i]];
    }

#else // show a summary only
    const char *categories[] = {"Haplotype data", "Other sample data", "Header and columns 1-9"};
    
    uint64_t vbytes[] = {vcf_file->section_bytes[SEC_HAPLOTYPE_DATA],
                         vcf_file->section_bytes[SEC_PHASE_DATA] + vcf_file->section_bytes[SEC_GENOTYPE_DATA] + vcf_file->section_bytes[SEC_STATS_HT_SEPERATOR],
                         vcf_file->section_bytes[SEC_VCF_HEADER] + vcf_file->section_bytes[SEC_VARIANT_DATA]
                        };

    uint64_t zbytes[] = {z_file->section_bytes[SEC_HAPLOTYPE_DATA],
                         z_file->section_bytes[SEC_PHASE_DATA] + z_file->section_bytes[SEC_GENOTYPE_DATA] + z_file->section_bytes[SEC_DICTIONARY],
                         z_file->section_bytes[SEC_VCF_HEADER] + z_file->section_bytes[SEC_VARIANT_DATA]
                        };
#endif

    for (unsigned i=0; i < sizeof(vbytes)/sizeof(vbytes[0]); i++) {
        fprintf (stderr, format, categories[i], 
                 buf_human_readable_size(vbytes[i], vsize), 100.0 * (double)vbytes[i] / (double)vcf_file->vcf_data_size,
                 buf_human_readable_size(zbytes[i], zsize), 100.0 * (double)zbytes[i] / (double)z_file->disk_size,
                 zbytes[i] ? (double)vbytes[i] / (double)zbytes[i] : 0,
                 !zbytes[i] ? (vbytes[i] ? "\b\b\bInf" : "\b\b\b---") : "");

        total_vcf += vbytes[i];
        total_z   += zbytes[i];
    }

    fprintf (stderr, format, "TOTAL", 
             buf_human_readable_size(total_vcf, vsize), 100.0 * (double)total_vcf / (double)vcf_file->vcf_data_size,
             buf_human_readable_size(total_z, zsize),   100.0 * (double)total_z   / (double)z_file->disk_size,
             (double)total_vcf / (double)total_z, "");

    ASSERTW (total_z == z_file->disk_size, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%"PRIu64" but file size is %"PRIu64" (diff=%"PRId64")", 
             total_z, z_file->disk_size, z_file->disk_size - total_z);

    ASSERTW (total_vcf == vcf_file->vcf_data_size, "Hmm... incorrect calculation for VCF sizes: total section sizes=%"PRIu64" but file size is %"PRIu64" (diff=%"PRId64")", 
             total_vcf, vcf_file->vcf_data_size, vcf_file->vcf_data_size - total_vcf);
}

static void main_genozip (const char *vcf_filename, 
                          char *z_filename,
                          int pipefd_zip_to_unzip,  // send output to pipe (used for testing)
                          unsigned max_threads,
                          bool is_first_file, bool is_last_file)
{
    File *vcf_file;
    static File *z_file = NULL; // static to support concat mode

    // get input file
    if (vcf_filename) 
        vcf_file = file_open (vcf_filename, READ, VCF);
    else {  // stdin
        vcf_file = file_fdopen (0, READ, STDIN, false);
        flag_stdout = (z_filename == NULL); // implicit setting of stdout by using stdin, unless -o was used
    }
 
    ASSERTW (!flag_stdout || !flag_md5, "%s: ignoring --md5 / -m option - it cannot be used when output is redirected", global_cmd);
  
    ASSERT0 (flag_concat_mode || !z_file, "Error: expecting z_file to be NULL in non-concat mode");

    // get output FILE
    if (!flag_stdout && pipefd_zip_to_unzip < 0) {

        if (!z_file) { // skip if we're the second file onwards in concatenation mode - nothing to do
            if (!z_filename) {
                unsigned fn_len = strlen (vcf_filename);
                z_filename = (char *)malloc (fn_len + strlen (GENOZIP_EXT) + 1);
                ASSERT(z_filename, "Error: Failed to malloc z_filename len=%u", fn_len+4);

                if (vcf_file->type == VCF)
                    sprintf (z_filename, "%s" GENOZIP_EXT, vcf_filename); // .vcf -> .vcf.genozip
                else if (vcf_file->type == VCF_GZ)
                    sprintf (z_filename, "%.*s" GENOZIP_EXT, fn_len-3, vcf_filename); // .vcf.gz -> .vcf.genozip
                else if (vcf_file->type == VCF_BZ2)
                    sprintf (z_filename, "%.*s" GENOZIP_EXT, fn_len-4, vcf_filename); // .vcf.bz2 -> .vcf.genozip
                else
                    ABORT ("Invalid file type: %u", vcf_file->type);
            }
            ASSERT (flag_force || access(z_filename, F_OK), 
                    "%s: output file %s already exists", global_cmd, z_filename);

            z_file = file_open (z_filename, WRITE, GENOZIP);
        }
        z_file->disk_at_beginning_of_this_vcf_file = z_file->disk_so_far;
        z_file->vcf_data_so_far                    = 0; // reset these as they relate to the VCF data of the VCF file currently being processed
        z_file->vcf_data_size                      = vcf_file->vcf_data_size; // 0 if vcf is .gz or a pipe or stdin
    }
    else if (flag_stdout) { // stdout
#ifdef _WIN32
        // this is because Windows redirection is in text (not binary) mode, meaning Windows edits the output stream...
        ASSERT (isatty(1), "%s: redirecting binary output is not supported on Windows", global_cmd);
#endif
        ASSERT (flag_force || !isatty(1), "%s: you must use --force to output a compressed file to the terminal", global_cmd);

        z_file = file_fdopen (1, WRITE, STDOUT, false);
    } 
    else if (pipefd_zip_to_unzip >= 0) {
        z_file = file_fdopen (pipefd_zip_to_unzip, WRITE, PIPE, true);
    }
    else ABORT0 ("Error: No output channel");
    
    const char *basename = file_basename (vcf_filename, false, "(stdin)", NULL, 0);
    zip_dispatcher (basename, vcf_file, z_file, pipefd_zip_to_unzip >= 0, max_threads, is_last_file);

    if (flag_show_content) main_display_section_stats (vcf_file, z_file);

    bool remove_vcf_file = z_file && flag_replace && vcf_filename;

    file_close (&vcf_file);

    if ((is_last_file || !flag_concat_mode) && !flag_stdout && z_file) 
        file_close (&z_file); 

    if (remove_vcf_file) file_remove (vcf_filename); 

    free ((void *)basename);
}

static void main_genounzip (const char *z_filename,
                            char *vcf_filename, 
                            int pipe_from_zip_thread, 
                            int pipe_to_test_thread,
                            unsigned max_threads)
{
    static File *vcf_file = NULL; 
    File *z_file;
    
    // get input FILE
    if (z_filename) {
        unsigned fn_len = strlen (z_filename);

        ASSERT (file_has_ext (z_filename, ".vcf" GENOZIP_EXT), "%s: file: %s - expecting a file with a .vcf" GENOZIP_EXT " extension", global_cmd, z_filename);

        if (!vcf_filename && !flag_stdout && pipe_to_test_thread < 0) {
            vcf_filename = (char *)malloc(fn_len + 10);
            ASSERT(vcf_filename, "Error: failed to malloc vcf_filename, len=%u", fn_len+10);

            sprintf (vcf_filename, "%.*s", (int)(fn_len - strlen(GENOZIP_EXT)), z_filename);    // .vcf.genozip -> .vcf
        }

        z_file = file_open (z_filename, READ, GENOZIP);    
    }
    else if (pipe_from_zip_thread >= 0) {
        z_file = file_fdopen (pipe_from_zip_thread, READ, GENOZIP, false);
    }
    else { // stdin
#ifdef _WIN32
        // this is because Windows redirection is in text (not binary) mode, meaning Windows edits the input stream...
        ABORT ("%s: redirecting binary input is not supported on Windows", global_cmd);
#endif
        z_file = file_fdopen (0, READ, STDIN, false);
    }

    // get output FILE
    if (vcf_filename) {

        ASSERT0 (!vcf_file || flag_concat_mode, "Error: vcf_file is open but not in concat mode");

        if (!vcf_file) { // in concat mode, for second file onwards, vcf_file is already open
            ASSERT (flag_force || access(vcf_filename, F_OK), "%s: output file %s already exists", global_cmd, vcf_filename);

            vcf_file = file_open (vcf_filename, WRITE, VCF);
        }
    }
    else if (pipe_to_test_thread >= 0) {
        vcf_file = file_fdopen (pipe_to_test_thread, WRITE, VCF, false);
    }
    else if (flag_stdout) { // stdout
        vcf_file = file_fdopen (1, WRITE, VCF, false); // STDOUT
    }
    else if (flag_split) {
        // do nothing - the files will be opened by zip_dispatcher
    }
    else {
        ABORT0 ("Error: unrecognized configuration for the vcf_file");
    }
    
    const char *basename = file_basename (z_filename, false, "(stdin)", NULL, 0);
    piz_dispatcher (basename, z_file, vcf_file, pipe_from_zip_thread >= 0, max_threads);

    if (!flag_concat_mode && !flag_stdout) 
        file_close (&vcf_file); // no worries about not closing the concatenated file - it will close with the process exits

    file_close (&z_file);

    free ((void *)basename);

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

    main_genozip (t2c->vcf_filename, NULL, t2c->pipe_to_uncompress_thread, t2c->max_threads, true, true);

    return NULL;
}

void *main_test_uncompress_thread_entry (void *p_)
{
    TestToUncompressData *t2u = (TestToUncompressData*)p_;

    main_genounzip (NULL, NULL, t2u->pipe_from_zip_thread, t2u->pipe_to_test_thread, t2u->max_threads);

    return NULL;
}

static void main_test (const char *vcf_filename)
{
    ASSERT (vcf_filename, "%s: filename missing", global_cmd);

    File *vcf_file = file_open (vcf_filename, READ, VCF);

    int pipefd_zip_to_unzip [2];
    int pipefd_unzip_to_main[2];
#ifdef _WIN32
    _pipe (pipefd_zip_to_unzip,  25000000, _O_BINARY); // 25MB pipe space to make sure both sides stay busy
    _pipe (pipefd_unzip_to_main, 25000000, _O_BINARY);
#else
    ASSERT (pipe (pipefd_zip_to_unzip)  >= 0, "Error: failed to create pipe 1: %s", strerror(errno));
    ASSERT (pipe (pipefd_unzip_to_main) >= 0, "Error: failed to create pipe 2: %s", strerror(errno));
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
    t2u.pipe_from_zip_thread      = pipefd_zip_to_unzip[0];
    t2u.pipe_to_test_thread       = pipefd_unzip_to_main[1];
    t2u.max_threads               = (global_max_threads-1) / 2;

    unsigned err = pthread_create(&test_compress_thread, NULL, main_test_compress_thread_entry, &t2c);
    ASSERT (!err, "Error: failed to create test compress thread, err=%u", err);

    err = pthread_create(&test_uncompress_thread, NULL, main_test_uncompress_thread_entry, &t2u);
    ASSERT (!err, "Error: failed to create test unzip thread, err=%u", err);

    FILE *from_pipe = fdopen (pipefd_unzip_to_main[0], "rb");
    ASSERT0 (from_pipe, "Failed to fdopen pipe from unzip");

    vcffile_compare_pipe_to_file (from_pipe, vcf_file);

    fclose (from_pipe);

    file_close(&vcf_file);

    return;
}

static void main_list_dir(); // forward declaration

static void main_list (const char *z_filename, bool finalize, const char *subdir) 
{
    if (!finalize) {
        // no specific filename = show entire directory
        if (!z_filename) {
            main_list_dir("."); 
            return;
        }

        // filename is a directory - show directory contents (but not recursively)
        if (!subdir) {
            struct stat st;
            int ret = stat (z_filename, &st);
            ASSERT (!ret, "Error: failed to stat(%s): %s", z_filename, strerror (errno));

            if (S_ISDIR (st.st_mode)) {
                main_list_dir (z_filename);
                return;
            }
        }
    }

    static bool first_file = true;
    static unsigned files_listed=0, files_ignored=0;
    static long long total_uncompressed_len=0, total_compressed_len=0;
    char c_str[20], u_str[20];

    const unsigned FILENAME_WIDTH = 50;

    const char *head_format = "%5s %8s %10s %10s %6s%s %*s %s\n";
    const char *foot_format = "\nTotal:         %10s %10s %5uX\n";
    const char *encr_format = "<encrypted>      %s                           %s%s%*s\n";
    const char *item_format = "%5u %s %10s %10s %5uX%s %s%s%*s %s\n";

    if (finalize) {
        if (files_listed > 1) {
            buf_human_readable_size(total_compressed_len, c_str);
            buf_human_readable_size(total_uncompressed_len, u_str);
            unsigned ratio = total_compressed_len ? ((double)total_uncompressed_len / (double)total_compressed_len) : 0;

            printf (foot_format, c_str, u_str, ratio);
        }
        
        if (files_ignored) printf ("\nIgnored %u files that do not have a .vcf" GENOZIP_EXT " extension\n\n", files_ignored);
        
        return;
    }

    if (first_file) {
        printf (head_format, "Indiv", "Sites", "Compressed", "Original", "Factor", flag_md5 ? " MD5                             " : "", -FILENAME_WIDTH, "Name", "Creation");
        first_file = false;
    }
    
    File *z_file = file_open(z_filename, READ, GENOZIP_TEST);    
    if (!z_file) {
        files_ignored++;
        return;
    }

    static SectionHeaderVCFHeader *vcf_header_header = NULL;
    if (!vcf_header_header) vcf_header_header = malloc (crypt_padded_len (sizeof (SectionHeaderVCFHeader)));
    ASSERT0 (vcf_header_header, "failed to malloc vcf_header_header");

    bool is_subdir = subdir && (subdir[0] != '.' || subdir[1] != '\0');

    bool encrypted;
    bool success = vcf_header_get_vcf_header (z_file, vcf_header_header, &encrypted);
    if (!success) {
        if (encrypted)
            printf (encr_format, 
                    flag_md5 ? "                                 " : "",
                    (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
                    is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH,
                    z_filename);
        else
            files_ignored++;
        return;
    }   

    uint64_t vcf_data_size = BGEN64 (vcf_header_header->vcf_data_size);
    uint32_t num_samples   = BGEN32 (vcf_header_header->num_samples);
    uint64_t num_lines     = BGEN64 (vcf_header_header->num_lines);

    char num_lines_str[50];
    if (num_lines != NUM_LINES_UNKNOWN)
#ifdef _MSC_VER        
        sprintf (num_lines_str, "%8I64u", num_lines);
#else
        sprintf (num_lines_str, "%8"PRIu64, num_lines);
#endif
    else
        sprintf (num_lines_str, "%-8s", "N/A");

    unsigned ratio = z_file->disk_size ? ((double)vcf_data_size / (double)z_file->disk_size) : 0;
    
    buf_human_readable_size(z_file->disk_size, c_str);
    buf_human_readable_size(vcf_data_size, u_str);
    printf (item_format, num_samples, num_lines_str, 
            c_str, u_str, ratio, 
            flag_md5 ? md5_display (vcf_header_header->md5_hash, true) : "",
            (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
            is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH,
            z_filename, vcf_header_header->created);
            
    total_compressed_len   += z_file->disk_size;
    total_uncompressed_len += vcf_data_size;
    
    files_listed++;

    file_close (&z_file);
}

static void main_list_dir(const char *dirname)
{
    DIR *dir;
    struct dirent *ent;

    dir = opendir (dirname);
    ASSERT (dir, "Error: failed to open directory: %s", strerror (errno));

    int ret = chdir (dirname);
    ASSERT (!ret, "Error: failed to chdir(%s)", dirname);

    while ((ent = readdir(dir))) {
        
        struct stat st;
        ret = stat (ent->d_name, &st);
        ASSERT (!ret, "Error: failed to stat(%s): %s", ent->d_name, strerror (errno));

        if (!S_ISDIR (st.st_mode))  // don't go down subdirectories recursively
            main_list (ent->d_name, false, dirname);
    
    }
    closedir(dir);    

    ret = chdir ("..");
    ASSERT0 (!ret, "Error: failed to chdir(..)");
}

int main (int argc, char **argv)
{
    #define COMPRESS   'z'
    #define UNCOMPRESS 'd'
    #define TEST       't'
    #define LIST       'l'
    #define LICENSE    'L'
    #define VERSION    'V'
    #define HELP       'h'

    ExeType exe_type;
    if      (strstr (argv[0], "genols"))    exe_type = EXE_GENOLS;
    else if (strstr (argv[0], "genocat"))   exe_type = EXE_GENOCAT;
    else if (strstr (argv[0], "genounzip")) exe_type = EXE_GENOUNZIP;
    else                                    exe_type = EXE_GENOZIP; // default
    
    buf_initialize();

    static int command = -1;  // must be static to initialize list_options 
    char *out_filename = NULL;
    char *threads_str = NULL;

    global_cmd = file_basename (argv[0], true, "(executable)", NULL, 0); // global var
    
    // verify CPU architecture and compiler is supported
    ASSERT0 (sizeof(char)==1 && sizeof(short)==2 && sizeof (unsigned)==4 && sizeof(long long)==8, 
             "Error: Unsupported C type lengths, check compiler options");
    
    // runtime endianity (byte order) detection - safer than compile-time
    unsigned test_endianity = 0x01020304;
    global_little_endian = *(uint8_t*)&test_endianity==0x04; // in big endian it is 0x01;

    bool is_short[256]; // indexed by character of short option.

    // process command line options
    while (1) {

        #define _c  {"stdout",     no_argument,       &flag_stdout,  1      }
        #define _d  {"decompress", no_argument,       &command, UNCOMPRESS  }
        #define _f  {"force",      no_argument,       &flag_force,   1      }
        #define _h  {"help",       no_argument,       &command, HELP        }
        #define _l  {"list",       required_argument, &command, LIST        }
        #define _L1 {"license",    no_argument,       &command, LICENSE     } // US spelling
        #define _L2 {"licence",    no_argument,       &command, LICENSE     } // British spelling
        #define _q  {"quiet",      no_argument,       &flag_quiet, 1        }
        #define _DL {"replace",    no_argument,       &flag_replace, 1      }
        #define _t  {"test",       no_argument,       &command, TEST        }
        #define _V  {"version",    no_argument,       &command, VERSION     }
        #define _z  {"compress",   no_argument,       &command, COMPRESS    }
        #define _m  {"md5",        no_argument,       &flag_md5, 1          }
        #define _th {"threads",    required_argument, 0, '@'                }
        #define _O  {"split",      no_argument,       &flag_split, 1        }
        #define _o  {"output",     required_argument, 0, 'o'                }
        #define _p  {"password",   required_argument, 0, 'p'                }
        #define _sc {"show-content",no_argument,      &flag_show_content, 1 } 
        #define _sa {"show-alleles",no_argument,      &flag_show_alleles, 1 }
        #define _st {"show-time"   ,no_argument,      &flag_show_time   , 1 } 
        #define _sm {"show-memory" ,no_argument,      &flag_show_memory , 1 } 
        #define _00 {0, 0, 0, 0                                             }

        typedef const struct option Option;
        static Option genozip_lo[]   = { _c, _d, _f, _h, _l, _L1, _L2, _q, _DL, _t, _V, _z, _m, _th,     _o, _p, _sc, _sa, _st, _sm, _00 };
        static Option genounzip_lo[] = { _c,     _f, _h,     _L1, _L2, _q, _DL,     _V,         _th, _O, _o, _p,           _st, _sm, _00 };
        static Option genols_lo[]    = {             _h,     _L1, _L2,              _V,     _m,              _p,                     _00 };
        static Option genocat_lo[]   = {             _h,     _L1, _L2,              _V,         _th,         _p,                     _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        static const char *short_options[] = { // same order as ExeType
            "cdfhlLq^tVzm@:o:p:", // genozip
            "cfhLq^V@:Oo:p:",     // genounzip
            "hLVmp:",             // genols
            "hLV@:p:"             // genocat
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case COMPRESS : case UNCOMPRESS : case LIST : case LICENSE : 
            case TEST     : case VERSION    : case HELP :
                ASSERT(command<0 || command==c, "%s: can't have both -%c and -%c", global_cmd, command, c); 
                command=c; 
                break;

            case 'c' : flag_stdout        = 1      ; break;
            case 'f' : flag_force         = 1      ; break;
            case '^' : flag_replace       = 1      ; break;
            case 'q' : flag_quiet         = 1      ; break;
            case 'm' : flag_md5           = 1      ; break;
            case 'O' : flag_split         = 1      ; break;
            case '@' : threads_str  = optarg       ; break;
            case 'o' : out_filename = optarg       ; break;
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
                abort();  
        }
    }

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        // genozip with no input filename, no output filename, and no output or input redirection - show help
        if (exe_type == EXE_GENOLS) command = LIST; // genols can be run without arguments
        
        else if (command == -1 && optind == argc && !out_filename && isatty(0) && isatty(1)) 
            return main_print_help (exe_type, false);

        else if (exe_type == EXE_GENOUNZIP) command = UNCOMPRESS;
        else if (exe_type == EXE_GENOCAT) { command = UNCOMPRESS; flag_stdout=1 ; }
        else command = COMPRESS; // default 
    }

    // sanity checks
    #define OT(l,s) is_short[(int)s[0]] ? "-"s : "--"l
    
    ASSERT (!flag_stdout || !out_filename,      "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("output", "o"));
    ASSERT (!flag_stdout || !flag_split,        "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("split", "O"));
    ASSERT (!flag_split  || !out_filename,      "%s: option %s is incompatable with %s", global_cmd, OT("split",  "O"), OT("output", "o"));
    ASSERT (!flag_stdout || !flag_replace,      "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("replace", "^"));
    ASSERT (!flag_stdout || !flag_show_content, "%s: option %s is incompatable with --show-content", global_cmd, OT("stdout", "c"));
    ASSERT (!flag_stdout || !flag_show_alleles, "%s: option %s is incompatable with --show-alleles", global_cmd, OT("stdout", "c"));
    ASSERTW (!flag_stdout       || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring %s option", global_cmd, OT("stdout", "c"));
    ASSERTW (!flag_force        || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring %s option", global_cmd, OT("force", "f"));
    ASSERTW (!flag_replace      || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring %s option", global_cmd, OT("replace", "^"));
    ASSERTW (!flag_quiet        || command == COMPRESS || command == UNCOMPRESS || command == TEST, "%s: ignoring %s option", global_cmd, OT("quiet", "q"));
    ASSERTW (!flag_md5          || command == COMPRESS || command == LIST      , "%s: ignoring %s option %s", global_cmd, OT("md5","m"),
             command==UNCOMPRESS ? "- decompress always verifies MD5 if the file was compressed with --md5" : "");
    ASSERTW (!threads_str       || command == COMPRESS || command == UNCOMPRESS || command == TEST, "%s: ignoring %s option", global_cmd, OT("threads", "@"));
    ASSERTW (!out_filename      || command == COMPRESS || command == UNCOMPRESS, "%s: ignoring %s option", global_cmd, OT("output", "o"));
    ASSERTW (!flag_show_content || command == COMPRESS || command == TEST      , "%s: ignoring --show-content, it only works with --compress or --test", global_cmd);
    ASSERTW (!flag_show_alleles || command == COMPRESS || command == TEST      , "%s: ignoring --show-alleles, it only works with --compress or --test", global_cmd);
    
    if (command != COMPRESS && command != LIST) flag_md5=false;
    
    if (command == UNCOMPRESS && flag_stdout) flag_quiet=true; // don't show progress when outputing to stdout

    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (threads_str) {
        int ret = sscanf (threads_str, "%u", &global_max_threads);
        ASSERT (ret == 1 && global_max_threads >= 1, "%s: %s requires an integer value of at least 1", global_cmd, OT("threads", "@"));
        ASSERT (command != TEST || (global_max_threads >= 3 && global_max_threads != 4), "%s invalid %s value: number of threads for %s should be at least 3 (but not 4)", global_cmd, OT("threads", "@"), OT("test", "t"));
    }
    else {
        global_max_threads = main_get_num_cores();
        
        if (global_max_threads < 3 && command == TEST) // with -t, we allow 3 threads even if we have only 1 or 2 cores
            global_max_threads = 3; 
    }

    if (command == TEST) {
        flag_stdout = flag_force = flag_replace = false;
        out_filename = NULL;
    }

    // take action, depending on the command selected
    if (command == VERSION) 
        return main_print_version();

    if (command == LICENSE) 
        return main_print (license, sizeof(license) / sizeof(char*), "", "\n\n", 
                           flag_force ? 60 : 0); // --license --force means output license in Windows installer format (used by Makefile) - 60 is width of InstallForge license text field

    if (command == HELP) 
        return main_print_help (exe_type, true);

    flag_concat_mode = (out_filename != NULL);

    unsigned count=0;
    do {
        char *next_input_file = optind < argc ? argv[optind++] : NULL;  // NULL means stdin
        
        if (next_input_file && !strcmp (next_input_file, "-")) next_input_file = NULL; // "-" is stdin too

        ASSERTW (next_input_file || !flag_replace, "%s: ignoring %s option", global_cmd, OT("replace", "^")); 
        
        ASSERT0 (!count || !flag_show_content, "Error: --show-content can only work on one file at time");

        switch (command) {
            case COMPRESS   : main_genozip (next_input_file, out_filename, -1, global_max_threads, !count, optind==argc); break;
            case UNCOMPRESS : main_genounzip (next_input_file, out_filename, -1, -1, global_max_threads); break;
            case TEST       : main_test  (next_input_file); break; // returns if successful, displays error and exits if not
            case LIST       : main_list  (next_input_file, false, NULL); break;
            
            default         : ASSERT(false, "%s: unrecognized command %c", global_cmd, command);
        }

        count++;
    } while (optind < argc);
            
    // if this is "list", finalize
    if (command == LIST) main_list (NULL, true, NULL);

    return 0;
}

