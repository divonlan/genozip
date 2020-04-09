// ------------------------------------------------------------------
//   genozip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef _MSC_VER // Microsoft compiler
#include <getopt.h>
#else
#include "compatibility/visual_c_getopt.h"
#endif
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#ifndef WIN32
#include <sys/ioctl.h>
#include <termios.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#else // LINUX
#include <sched.h>
#include <sys/sysinfo.h>
#endif
#endif

#include "genozip.h"
#include "text_license.h"
#include "text_help.h"
#include "version.h" // automatically incremented by the make when we create a new distribution
#include "vcffile.h"
#include "zfile.h"
#include "zip.h"
#include "piz.h"
#include "crypt.h"
#include "file.h"
#include "dict_id.h"
#include "vcf_vb.h"
#include "endianness.h"
#include "regions.h"
#include "samples.h"
#include "gtshark.h"
#include "stream.h"
#include "url.h"

// globals - set it main() and never change
const char *global_cmd = NULL; 
ExeType exe_type;
int command = -1;  // must be static or global to initialize list_options 

uint32_t global_max_threads = DEFAULT_MAX_THREADS; 

// the flags - representing command line options - available globally
int flag_quiet=0, flag_force=0, flag_concat=0, flag_md5=0, flag_split=0, flag_optimize=0, flag_bgzip=0,
    flag_show_alleles=0, flag_show_time=0, flag_show_memory=0, flag_show_dict=0, flag_show_gt_nodes=0,
    flag_show_b250=0, flag_show_sections=0, flag_show_headers=0, flag_show_index=0, flag_show_gheader=0, flag_show_threads=0,
    flag_stdout=0, flag_replace=0, flag_show_content=0, flag_test=0, flag_regions=0, flag_samples=0,
    flag_drop_genotypes=0, flag_no_header=0, flag_header_only=0, flag_header_one=0, flag_noisy=0, flag_strip=0,
    flag_debug_memory=0, flag_show_vblocks=0, flag_gtshark=0, flag_sblock=0, flag_vblock=0, flag_gt_only=0;

DictIdType dict_id_show_one_b250 = { 0 },  // argument of --show-b250-one
           dict_id_show_one_dict = { 0 },  // argument of --show-dict-one
           dict_id_dump_one_b250 = { 0 };  // argument of --dump-b250-one

static char *threads_str  = NULL;

void my_exit(void) 
{
    url_kill_curl();
    file_kill_external_compressors(); 

    // if we're in ZIP - remove failed genozip file (but don't remove partial failed VCF file in PIZ - it might be still useful to the user)
    if (command == ZIP && z_file && z_file->name) {
        char *save_name = malloc (strlen (z_file->name)+1);
        strcpy (save_name, z_file->name);

        file_close (&z_file, false); // also frees file->name

        file_remove (save_name, true);
    }

    exit(1);
} 

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
        line_width = MAX (40, w.ws_col); // our wrapper cannot work with to small line widths 
#endif
    }

    for (unsigned i=0; i < num_lines; i++)  {
        const char *line = text[i];
        unsigned line_len = strlen (line);
        
        // print line with line wraps
        bool wrapped = false;
        while (line_len + (wrapped ? strlen (wrapped_line_prefix) : 0) > line_width) {
            int c; for (c=line_width-1 - (wrapped ? strlen (wrapped_line_prefix) : 0); 
                        c>=0 && ((line[c] >= 'A' && line[c] <= 'Z') || 
                                 (line[c] >= 'a' && line[c] <= 'z') ||
                                 (line[c] >= '0' && line[c] <= '9')); // wrap lines at - and | too, so we can break very long regex strings like in genocat
                        c--); // find 
            printf ("%s%.*s\n", wrapped ? wrapped_line_prefix : "", c, line);
            line += c + (line[c]==' '); // skip space too
            line_len -= c + (line[c]==' ');
            wrapped = true;
        }
        printf ("%s%s%s", wrapped ? wrapped_line_prefix : "", line, newline_separator);
    }
    return 0;
}

static int main_print_help (bool explicit)
{
    static const char **texts[] = {help_genozip, help_genounzip, help_genols, help_genocat, help_genozip_developer}; // same order as ExeType
    static unsigned sizes[] = {sizeof(help_genozip), sizeof(help_genounzip), sizeof(help_genols), sizeof(help_genocat), sizeof(help_genozip_developer)};
    
    if (flag_force) exe_type = (ExeType)4; // -h -f shows developer help

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

static bool main_am_i_in_docker (void)
{
    FILE *fp = fopen ("/.dockerenv", "r");
    if (!fp) return false;

    FCLOSE (fp, "/.dockerenv");
    return true;
}

static void main_show_file_metadata (void)
{
    fprintf (stderr, "\n\n");
    if (txt_file->name) fprintf (stderr, "File name: %s\n", txt_file->name);
    fprintf (stderr, 
#ifdef _MSC_VER
             "Individuals: %u   Variants: %I64u   Non-GT subfields: %u\n", 
#else
             "Individuals: %u   Variants: %"PRIu64"   INFO & FORMAT tags: %u\n", 
#endif
             global_num_samples, z_file->num_lines, z_file->num_dict_ids-8);
}

static void main_show_sections (void)
{
    main_show_file_metadata();

    char vsize[30], zsize[30], zentries_str[30];

    fprintf (stderr, "Sections stats:\n");
    fprintf (stderr, "                           #Sec   #Entries                VCF     %%        GENOZIP     %%   Ratio\n");
    const char *format = "%22s    %6u  %16s  %9s %5.1f      %9s %5.1f  %6.1f%s\n";

    // the order in which we want them displayed
    const SectionType secs[] = {
        SEC_GENOZIP_HEADER, SEC_RANDOM_ACCESS,
        SEC_VCF_HEADER, SEC_VB_HEADER,
        SEC_CHROM_B250, SEC_CHROM_DICT, SEC_POS_B250, SEC_POS_DICT, 
        SEC_ID_B250, SEC_ID_DICT, SEC_REFALT_B250, SEC_REFALT_DICT, SEC_QUAL_B250, SEC_QUAL_DICT,
        SEC_FILTER_B250, SEC_FILTER_DICT, SEC_INFO_B250, SEC_INFO_DICT,
        SEC_INFO_SUBFIELD_B250, SEC_INFO_SUBFIELD_DICT, SEC_FORMAT_B250, SEC_FORMAT_DICT,
        SEC_GENOTYPE_DATA, SEC_FRMT_SUBFIELD_DICT,
        SEC_HAPLOTYPE_DATA, SEC_STATS_HT_SEPERATOR, SEC_PHASE_DATA
    };

    static const char *categories[] = {
        "Genozip header", "Random access index", "VCF header", "Variant block metadata", 
        "CHROM b250", "CHROM dict", "POS b250", "POS dict", "ID b250", "ID dict", "REF+ALT b250", "REF+ALT dict", 
        "QUAL b250", "QUAL dict", "FILTER b250", "FILTER dict",
        "INFO names b250", "INFO names dict", "INFO values b250", "INFO values dict", 
        "FORMAT b250", "FORMAT dict", "Tags b250", "Tags dict",
        "Haplotype data", "HT separator char", "Phasing char"
    };

    unsigned num_secs = sizeof(secs)/sizeof(secs[0]);
    ASSERT0 (sizeof(categories)/sizeof(categories[0]) == num_secs, "Error: categories and secs are not the same length");

    int64_t total_vcf=0, total_z=0, total_entries=0;
    uint32_t total_sections=0;

    for (unsigned sec_i=0; sec_i < num_secs; sec_i++) {
        int64_t vbytes    = txt_file->section_bytes[secs[sec_i]];
        int64_t zbytes    = z_file->section_bytes[secs[sec_i]];
        int64_t zentries  = z_file->section_entries[secs[sec_i]];
        int32_t zsections = z_file->num_sections[secs[sec_i]];

        char *vcf_size_str = (vbytes || section_type_is_dictionary (sec_i)) ? buf_display_size(vbytes, vsize) : "       ";
        
        fprintf (stderr, format, categories[sec_i], zsections, buf_display_uint (zentries, zentries_str),
                 vcf_size_str, 100.0 * (double)vbytes / (double)txt_file->txt_data_size_single,
                 buf_display_size(zbytes, zsize), 100.0 * (double)zbytes / (double)z_file->disk_size,
                 zbytes ? (double)vbytes / (double)zbytes : 0,
                 !zbytes ? (vbytes ? "\b\b\bInf" : "\b\b\b---") : "");

        total_sections += zsections;
        total_entries  += zentries;
        total_vcf      += vbytes;
        total_z        += zbytes;
    }

    fprintf (stderr, format, "TOTAL", total_sections, buf_display_uint (total_entries, zentries_str),
             buf_display_size(total_vcf, vsize), 100.0 * (double)total_vcf / (double)txt_file->txt_data_size_single,
             buf_display_size(total_z, zsize),   100.0 * (double)total_z   / (double)z_file->disk_size,
             (double)total_vcf / (double)total_z, "");

    fprintf (stderr, "\nDictionaries:\n");
    fprintf (stderr, "Name     Type         #Words        #Uniq         Hash    uncomp_dict_size\n");
    for (uint32_t i=0; i < z_file->num_dict_ids; i++) { // don't show CHROM-FORMAT as they are already showed above
        const MtfContext *ctx = &z_file->mtf_ctx[i];
    
        char s1[20], s2[20], s3[20];
        fprintf (stderr, "%*.*s %-6.6s %12s %12s %12s %9s\n", -DICT_ID_LEN, DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, 
                 dict_id_display_type (ctx->dict_id), buf_display_uint (ctx->mtf_i.len, s1), buf_display_uint (ctx->mtf.len, s2), 
                 buf_display_uint (ctx->global_hash_prime, s3), buf_display_size(ctx->dict.len, vsize));
    }

    char s1[20], s2[20];
    ASSERTW (total_z == z_file->disk_size, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
             buf_display_uint (total_z, s1), buf_display_uint (z_file->disk_size, s2), (int32_t)(z_file->disk_size - total_z));

    ASSERTW (total_vcf == txt_file->txt_data_size_single, "Hmm... incorrect calculation for VCF sizes: total section sizes=%s but file size is %s (diff=%d)", 
             buf_display_uint (total_vcf, s1), buf_display_uint (txt_file->txt_data_size_single, s2), (int32_t)(txt_file->txt_data_size_single - total_vcf));

}

static void main_show_content (void)
{
    main_show_file_metadata();

    char vsize[30], zsize[30];
    int64_t total_vcf=0, total_z=0;

    fprintf (stderr, "Compression stats:\n");
    fprintf (stderr, "                              VCF     %%       GENOZIP     %%  Ratio\n");
    const char *format = "%22s   %8s %5.1f      %8s %5.1f  %5.1f%s\n";

    const char *categories[] = {"Haplotype data", "Other sample data", "Header and columns 1-9", "Index"};

#define NUM_CATEGORIES 4

    int sections_per_category[NUM_CATEGORIES][30] = { 
        { SEC_HAPLOTYPE_DATA, NIL },
        { SEC_PHASE_DATA, SEC_GENOTYPE_DATA, SEC_FRMT_SUBFIELD_DICT, SEC_STATS_HT_SEPERATOR, NIL},
        { SEC_VCF_HEADER, SEC_VB_HEADER, SEC_CHROM_B250, SEC_POS_B250, SEC_ID_B250, SEC_REFALT_B250, 
          SEC_QUAL_B250, SEC_FILTER_B250, SEC_INFO_B250, SEC_FORMAT_B250, SEC_INFO_SUBFIELD_B250, 
          SEC_CHROM_DICT, SEC_POS_DICT, SEC_ID_DICT, SEC_REFALT_DICT, SEC_QUAL_DICT,
          SEC_FILTER_DICT, SEC_INFO_DICT, SEC_INFO_SUBFIELD_DICT, SEC_FORMAT_DICT, NIL },
        { SEC_RANDOM_ACCESS, SEC_GENOZIP_HEADER, NIL }
    };

    for (unsigned i=0; i < NUM_CATEGORIES; i++) {

        int64_t vbytes=0, zbytes=0;
        for (int *sec_i = sections_per_category[i]; *sec_i != NIL; sec_i++) {
            vbytes += txt_file->section_bytes[*sec_i];
            zbytes += z_file->section_bytes[*sec_i];
        }

        fprintf (stderr, format, categories[i], 
                 buf_display_size(vbytes, vsize), 100.0 * (double)vbytes / (double)txt_file->txt_data_size_single,
                 buf_display_size(zbytes, zsize), 100.0 * (double)zbytes / (double)z_file->disk_size,
                 zbytes ? (double)vbytes / (double)zbytes : 0,
                 !zbytes ? (vbytes ? "\b\b\bInf" : "\b\b\b---") : "");

        total_vcf      += vbytes;
        total_z        += zbytes;
    }

    fprintf (stderr, format, "TOTAL", 
             buf_display_size(total_vcf, vsize), 100.0 * (double)total_vcf / (double)txt_file->txt_data_size_single,
             buf_display_size(total_z, zsize),   100.0 * (double)total_z   / (double)z_file->disk_size,
             (double)total_vcf / (double)total_z, "");

    ASSERTW (total_z == z_file->disk_size, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%"PRId64" but file size is %"PRId64" (diff=%"PRId64")", 
             total_z, z_file->disk_size, z_file->disk_size - total_z);

    ASSERTW (total_vcf == txt_file->txt_data_size_single, "Hmm... incorrect calculation for VCF sizes: total section sizes=%"PRId64" but file size is %"PRId64" (diff=%"PRId64")", 
             total_vcf, txt_file->txt_data_size_single, txt_file->txt_data_size_single - total_vcf);

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
    char c_str[20], u_str[20], s_str[20];

    const unsigned FILENAME_WIDTH = 50;

    const char *head_format = "\n%5s %10s %10s %10s %6s %s  %*s %s\n";
    const char *foot_format = "\nTotal:           %10s %10s %5uX\n";
    const char *item_format = "%5u %10s %10s %10s %5uX %s  %s%s%*s %s\n";

    // we accumulate the string in str_buf and print in the end - so it doesn't get mixed up with 
    // warning messages regarding individual files
    static Buffer str_buf = EMPTY_BUFFER; 
    
    if (finalize) {
        if (files_listed > 1) {
            buf_display_size(total_compressed_len, c_str);
            buf_display_size(total_uncompressed_len, u_str);
            unsigned ratio = total_compressed_len ? ((double)total_uncompressed_len / (double)total_compressed_len) : 0;

            bufprintf (evb, &str_buf, foot_format, c_str, u_str, ratio);
        }
        
        ASSERTW (!files_ignored, "\nIgnored %u file%s that %s not have a .vcf" GENOZIP_EXT " extension\n\n", 
                 files_ignored, files_ignored==1 ? "" : "s", files_ignored==1 ? "does" : "do");
        
        goto finish;
    }

    if (first_file) {
        bufprintf (evb, &str_buf, head_format, "Indiv", "Sites", "Compressed", "Original", "Factor", " MD5 (of original VCF)           ", -(int)FILENAME_WIDTH, "Name", "Creation");
        first_file = false;
    }
    
    if (!file_has_ext (z_filename, VCF_GENOZIP_) || access (z_filename, F_OK)!=0) {
        files_ignored++;
        goto finish;
    }

    bool is_subdir = subdir && (subdir[0] != '.' || subdir[1] != '\0');

    z_file = file_open (z_filename, READ, Z_FILE, 0); // open global z_file
    uint64_t txt_data_size, num_lines;
    uint32_t num_samples;
    Md5Hash md5_hash_concat;
    char created[FILE_METADATA_LEN];
    bool success = zfile_get_genozip_header (&txt_data_size, &num_samples, &num_lines, 
                                             &md5_hash_concat, created, FILE_METADATA_LEN);
    if (!success) goto finish;

    unsigned ratio = z_file->disk_size ? ((double)txt_data_size / (double)z_file->disk_size) : 0;
    
    buf_display_size (z_file->disk_size, c_str);
    buf_display_size (txt_data_size, u_str);
    buf_display_uint (num_lines, s_str);

    bufprintf (evb, &str_buf, item_format, num_samples, s_str, 
               c_str, u_str, ratio, 
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
    vcf_header_initialize();

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

    if (!txt_filename && (!flag_stdout || flag_bgzip) && !flag_split) {
        txt_filename = (char *)malloc(fn_len + 10);
        ASSERT(txt_filename, "Error: failed to malloc txt_filename, len=%u", fn_len+10);

        // .vcf.genozip -> .vcf or .vcf.gz
        sprintf ((char *)txt_filename, "%.*s%s", (int)(fn_len - strlen(GENOZIP_EXT)), z_filename,
                 flag_bgzip ? ".gz" : "");    
    }

    z_file = file_open (z_filename, READ, Z_FILE, 0);    

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
        // do nothing - the vcf component files will be opened by vcf_header_genozip_to_vcf()
    }
    else {
        ABORT0 ("Error: unrecognized configuration for the txt_file");
    }
    
    const char *basename = file_basename (z_filename, false, "(stdin)", NULL, 0);
    
    // a loop for decompressing all vcf components in split mode. in non-split mode, it collapses to one a single iteration.
    bool piz_successful;
    unsigned num_vcf_components=0;
    do {
        piz_successful = piz_dispatcher (basename, max_threads, num_vcf_components==0, is_last_file);
        if (piz_successful) num_vcf_components++;
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
    ASSERTW (!exit_code, "genozip test exited with status %d\n", exit_code);
}

static void main_genozip (const char *txt_filename, 
                          char *z_filename,
                          unsigned max_threads,
                          bool is_first_file, bool is_last_file,
                          char *exec_name)
{
    ASSERT (!z_filename || !url_is_url (z_filename), 
            "%s: output files must be regular files, they cannot be a URL: %s", global_cmd, z_filename);

    // get input file
    if (txt_filename) {
        // open the file
        txt_file = file_open (txt_filename, READ, TXT_FILE, 0);

        // skip this file if its size is 0
        RETURNW (txt_file,, "Cannot compresss file %s because its size is 0 - skipping it", txt_filename);
    }
    else {  // stdin
        txt_file = file_open_redirect (READ, TXT_FILE, UNKNOWN_FILE_TYPE);
        flag_stdout = (z_filename == NULL); // implicit setting of stdout by using stdin, unless -o was used
    }
 
    ASSERT0 (flag_concat || !z_file, "Error: expecting z_file to be NULL in non-concat mode");

    // get output FILE
    if (!flag_stdout) {

        if (!z_file) { // skip if we're the second file onwards in concatenation mode - nothing to do

            if (!z_filename) {
                const char *basename = url_is_url (txt_filename) ? file_basename (txt_filename, false, "", 0,0) : NULL;
                const char *local_vcf_filename = basename ? basename : txt_filename;

                unsigned fn_len = strlen (local_vcf_filename);
                z_filename = (char *)malloc (fn_len + strlen (GENOZIP_EXT) + 1);
                ASSERT(z_filename, "Error: Failed to malloc z_filename len=%u", fn_len+4);

                // get name, e.g. xx.bcf.gz -> xx.vcf.genozip
                sprintf (z_filename, "%.*s" VCF_GENOZIP_, (int)(fn_len - strlen (file_exts[txt_file->type])), local_vcf_filename); 

                if (basename) FREE ((char*)basename);
            }
            z_file = file_open (z_filename, WRITE, Z_FILE, txt_file->data_type);
        }
    }
    else if (flag_stdout) { // stdout
#ifdef _WIN32
        // this is because Windows redirection is in text (not binary) mode, meaning Windows edits the output stream...
        ASSERT (isatty(1), "%s: redirecting binary output is not supported on Windows", global_cmd);
#endif
        ASSERT (flag_force || !isatty(1), "%s: you must use --force to output a compressed file to the terminal", global_cmd);

        z_file = file_open_redirect (WRITE, Z_FILE, txt_file->data_type);
    } 
    else ABORT0 ("Error: No output channel");
    
    const char *basename = file_basename (txt_filename, false, "(stdin)", NULL, 0);
    zip_dispatcher (basename, max_threads, is_last_file);

    if (flag_show_sections && is_last_file) main_show_sections();
    if (flag_show_content && is_last_file)  main_show_content();

    bool remove_vcf_file = z_file && flag_replace && txt_filename;

    file_close (&txt_file, !is_last_file);

    if ((is_last_file || !flag_concat) && !flag_stdout && z_file) 
        file_close (&z_file, !is_last_file); 

    if (remove_vcf_file) file_remove (txt_filename, true); 

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

void verify_architecture()
{
    // verify CPU architecture and compiler is supported
    ASSERT0 (sizeof(char)==1 && sizeof(short)==2 && sizeof (unsigned)==4 && sizeof(long long)==8, 
             "Error: Unsupported C type lengths, check compiler options");
    
    // verify endianity is as expected
    uint16_t test_endianity = 0x0102;
#if defined __LITTLE_ENDIAN__
    ASSERT0 (*(uint8_t*)&test_endianity==0x02, "Error: expected CPU to be Little Endian but it is not");
#elif defined __BIG_ENDIAN__
    ASSERT0 (*(uint8_t*)&test_endianity==0x01, "Error: expected CPU to be Big Endian but it is not");
#else
#error  "Neither __BIG_ENDIAN__ nor __LITTLE_ENDIAN__ is defined - is endianness.h included?"
#endif

// Verify that this Windows is 64 bit
#ifdef _WIN32
#ifndef _WIN64
#error Compilation error - on Windows, genozip must be compiled as a 64 bit application
#endif
#endif
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
    const char *err_msg = "Error: invalid argument of --vblock: %s. Expecting an integer between 1 and 2048. The VCF file will be read and processed in blocks of this number of megabytes.";

    unsigned len = strlen (mem_size_mb_str);
    ASSERT (len <= 4 || (len==1 && mem_size_mb_str[0]=='0'), err_msg, mem_size_mb_str);

    for (unsigned i=0; i < len; i++) {
        ASSERT (mem_size_mb_str[i] >= '0' && mem_size_mb_str[i] <= '9', err_msg, mem_size_mb_str);
    }

    unsigned mem_size_mb = atoi (mem_size_mb_str);
    ASSERT (mem_size_mb <= 2048, err_msg, mem_size_mb_str);

    global_max_memory_per_vb = mem_size_mb * 1024 * 1024;
}

int main (int argc, char **argv)
{
#ifdef _WIN32
    // lowercase argv to allow case-insensitive comparison in Windows
    for (char *c=argv[0]; *c; c++) 
        if (*c >= 'A' && *c <= 'Z') 
            *c += 'a' - 'A';
#endif

    if      (strstr (argv[0], "genols"))    exe_type = EXE_GENOLS;
    else if (strstr (argv[0], "genocat"))   exe_type = EXE_GENOCAT;
    else if (strstr (argv[0], "genounzip")) exe_type = EXE_GENOUNZIP;
    else                                    exe_type = EXE_GENOZIP; // default
    
    verify_architecture();

    buf_initialize();

    char *out_filename = NULL;

    global_cmd = file_basename (argv[0], true, "(executable)", NULL, 0); // global var

    bool is_short[256]; // indexed by character of short option.

    // process command line options
    while (1) {

        #define _i  {"stdin",         required_argument, 0, 'i'                }
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
        #define _m  {"md5",           no_argument,       &flag_md5,          1 }
        #define _t  {"test",          no_argument,       &flag_test,         1 }
        #define _9  {"optimize",      no_argument,       &flag_optimize,     1 } // US spelling
        #define _9a {"optimise",      no_argument,       &flag_optimize,     1 } // British spelling
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
        #define _G  {"drop-genotypes",no_argument,       &flag_drop_genotypes,1}
        #define _H1 {"no-header",     no_argument,       &flag_no_header,    1 }
        #define _H0 {"header-only",   no_argument,       &flag_header_only,  1 }
        #define _1  {"header-one",    no_argument,       &flag_header_one,   1 }
        #define _sp {"strip",         no_argument,       &flag_strip,        1 }
        #define _GT {"GT-only",       no_argument,       &flag_gt_only,      1 }
        #define _Gt {"gt-only",       no_argument,       &flag_gt_only,      1 }
        #define _sc {"show-content",  no_argument,       &flag_show_content, 1 } 
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
        #define _00 {0, 0, 0, 0                                                }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _i, _c, _d, _f, _h, _l, _L1, _L2, _q, _Q, _t, _DL, _V,     _m, _th, _O, _o, _p,                                               _sc, _ss, _sd, _sT, _d1, _d2, _sg, _s2, _s5, _s6, _s7, _s8, _sa, _st, _sm, _sh, _si, _sr, _sv, _B, _S, _dm, _9, _9a, _gt, _00 };
        static Option genounzip_lo[]  = {     _c,     _f, _h,     _L1, _L2, _q, _Q, _t, _DL, _V, _z, _m, _th, _O, _o, _p,                                                         _sd, _sT, _d1, _d2,      _s2, _s5, _s6,      _st, _sm, _sh, _si,                             _dm,               _00 };
        static Option genocat_lo[]    = {             _f, _h,     _L1, _L2, _q, _Q,          _V,         _th,     _o, _p, _r, _tg, _s, _G, _1, _H0, _H1, _sp, _Gt, _GT,                _sT,                                    _st, _sm,                                       _dm,               _00 };
        static Option genols_lo[]     = {             _f, _h,     _L1, _L2, _q,              _V,                      _p,                                                                                                      _st, _sm,                                       _dm,               _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static const char *short_options[] = { // same order as ExeType
            "icdfhlLqQt^Vzm@:Oo:p:B:S:9K", // genozip
            "czfhLqQt^V@:Oo:p:m",          // genounzip
            "hLVp:qf",                     // genols
            "hLV@:p:qQ1r:t:s:H1Go:f"       // genocat
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

            case 'i' : file_set_stdin_type (optarg); break;
            case 'c' : flag_stdout        = 1      ; break;
            case 'z' : flag_bgzip         = 1      ; break;
            case 'f' : flag_force         = 1      ; break;
            case '^' : flag_replace       = 1      ; break;
            case 'q' : flag_quiet         = 1      ; break;
            case 'Q' : flag_noisy         = 1      ; break;
            case '9' : flag_optimize      = 1      ; break;
            case 'K' : flag_gtshark       = 1      ; break;
            case 't' : if (exe_type != EXE_GENOCAT) { flag_test = 1 ; break; }
                       // fall through for genocat -r
            case 'r' : flag_regions = true; regions_add (optarg); break;
            case 's' : flag_samples = true; samples_add (optarg); break;
            case 'm' : flag_md5           = 1      ; break;
            case 'O' : flag_split         = 1      ; break;
            case 'G' : flag_drop_genotypes= 1      ; break;
            case 'H' : flag_no_header     = 1      ; break;
            case '1' : flag_header_one    = 1      ; break;
            case '@' : threads_str  = optarg       ; break;
            case 'o' : out_filename = optarg       ; break;
            case '2' : dict_id_show_one_b250 = dict_id_make (optarg, strlen (optarg)); break;
            case '5' : dict_id_dump_one_b250 = dict_id_make (optarg, strlen (optarg)); break;
            case '3' : dict_id_show_one_dict = dict_id_make (optarg, strlen (optarg)); break;
            case 'B' : genozip_set_global_max_memory_per_vb (optarg); 
                       flag_vblock = true;
                       break;
            case 'S' : zip_set_global_samples_per_block (optarg); 
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
        
        // genozip with no input filename, no output filename, and no output or input redirection - show help
        // note: in docker stdin is a pipe even if going to a terminal. so we show the help even if
        // coming from a pipe. the user must use "-" to redirect from stdin
        else if (command == -1 && optind == argc && !out_filename && 
                 (isatty(0) || main_am_i_in_docker()) && isatty(1)) 
            return main_print_help (false);

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
    ASSERT (!flag_stdout      || !flag_md5,                         "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("md5", "m"));
    ASSERT (!flag_test        || !out_filename || command != UNZIP, "%s: option %s is incompatable with %s", global_cmd, OT("output", "o"),  OT("test", "t"));
    ASSERT (!flag_test        || !flag_replace || command != UNZIP, "%s: option %s is incompatable with %s", global_cmd, OT("replace", "^"), OT("test", "t"));
    ASSERT (!flag_test        || !flag_stdout  || command != ZIP,   "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("test", "t"));
    ASSERT (!flag_header_only || !flag_no_header,              "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), "header-only");
    ASSERT (!flag_quiet       || !flag_noisy,                       "%s: option %s is incompatable with %s", global_cmd, OT("quiet", "q"), OT("noisy", "Q"));
    ASSERT (!flag_test        || !flag_optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("test", "t"), OT("optimize", "9"));
    ASSERT (!flag_md5         || !flag_optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("md5", "m"), OT("optimize", "9"));
    ASSERT (!flag_samples     || !flag_drop_genotypes,              "%s: option %s is incompatable with %s", global_cmd, OT("samples", "s"), OT("drop-genotypes", "G"));
    ASSERT (!flag_no_header   || !flag_header_one,                  "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), OT("header-one", "1"));

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
    if (!flag_vblock) genozip_set_global_max_memory_per_vb (VCF_DATA_PER_VB); 
    if (!flag_sblock) zip_set_global_samples_per_block     (SAMPLES_PER_BLOCK); 

    // if using the -o option - check that we don't have duplicate filenames (even in different directory) as they
    // will overwrite each other if extracted with --split
    if (command == ZIP && out_filename && !flag_quiet) main_warn_if_duplicates (argc, argv, out_filename);

    unsigned num_files = argc - optind;

    flag_concat = (command == ZIP) && (out_filename != NULL) && (num_files > 1);

    ASSERT (num_files <= 1 || flag_concat || !flag_show_sections, "%s: --show-sections can only work on one file at time", global_cmd);
    ASSERT (num_files <= 1 || flag_concat || !flag_show_content,  "%s: --show-content can only work on one file at time", global_cmd);

    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (threads_str) {
        int ret = sscanf (threads_str, "%u", &global_max_threads);
        ASSERT (ret == 1 && global_max_threads >= 1, "%s: %s requires an integer value of at least 1", global_cmd, OT("threads", "@"));
    }
    else global_max_threads = main_get_num_cores();
    
    // take action, depending on the command selected
    if (command == VERSION) 
        return main_print_version();

    if (command == LICENSE) 
        return main_print (license, sizeof(license) / sizeof(char*), "", "\n\n", 
                           flag_force ? 60 : 0); // --license --force means output license in Windows installer format (used by Makefile) - 60 is width of InstallForge license text field

    if (command == HELP) 
        return main_print_help (true);

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

