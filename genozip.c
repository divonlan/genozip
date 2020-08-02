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
#include "reference.h"
#include "context.h"

// globals - set it main() and never change
const char *global_cmd = NULL; 
ExeType exe_type;

// primary_command vs command: primary_command is what the user typed on the command line. command is what is 
// running now - for example, when ZIP is unzipping a reference, primary_command=ZIP and command=PIZ
CommandType command = NO_COMMAND, primary_command = NO_COMMAND; 

uint32_t global_max_threads = DEFAULT_MAX_THREADS; 
uint32_t global_max_memory_per_vb = 0; // ZIP only: used for reading text file data

// the flags - representing command line options - available globally
int flag_quiet=0, flag_force=0, flag_bind=0, flag_md5=0, flag_unbind=0, flag_optimize=0, flag_bgzip=0, flag_bam=0, flag_bcf=0,
    flag_show_alleles=0, flag_show_time=0, flag_show_memory=0, flag_show_dict=0, flag_show_gt_nodes=0, flag_multiple_files=0,
    flag_show_b250=0, flag_show_stats=0, flag_show_headers=0, flag_show_index=0, flag_show_gheader=0, flag_show_threads=0,
    flag_stdout=0, flag_replace=0, flag_test=0, flag_regions=0, flag_samples=0, flag_fast=0, flag_list_chroms=0,
    flag_drop_genotypes=0, flag_no_header=0, flag_header_only=0, flag_header_one=0, flag_noisy=0, flag_show_aliases=0,
    flag_show_vblocks=0, flag_gtshark=0, flag_sblock=0, flag_vblock=0, flag_gt_only=0, flag_fasta_sequential=0,
    flag_debug_memory=0, flag_debug_progress=0, flag_show_hash, flag_register=0, flag_debug_no_singletons=0, flag_genocat_info_only=0,
    flag_reading_reference=0, flag_make_reference=0, flag_show_reference=0, flag_show_ref_index=0, flag_show_ref_hash=0, flag_ref_whole_genome=0,
    flag_optimize_sort=0, flag_optimize_PL=0, flag_optimize_GL=0, flag_optimize_GP=0, flag_optimize_VQSLOD=0, 
    flag_optimize_QUAL=0, flag_optimize_Vf=0, flag_optimize_ZM=0, flag_optimize_DESC=0, flag_optimize_SEQ=0,
    flag_pair=NOT_PAIRED_END;

ReferenceType flag_reference = REF_NONE;
uint64_t flag_stdin_size = 0;
char *flag_grep = NULL;
char *flag_show_is_set = NULL;

DictId dict_id_show_one_b250  = DICT_ID_NONE,  // argument of --show-b250-one
       dict_id_show_one_dict  = DICT_ID_NONE,  // argument of --show-dict-one
       dump_one_b250_dict_id  = DICT_ID_NONE,  // argument of --dump-b250-one
       dump_one_local_dict_id = DICT_ID_NONE;  // argument of --dump-local-one
static char *threads_str  = NULL;
static char *out_filename = NULL; // the filename specified in --output

static void print_call_stack (void) 
{
#ifndef _WIN32
#   define STACK_DEPTH 15
    void *array[STACK_DEPTH];
    size_t size = backtrace(array, STACK_DEPTH);
    
    fprintf (stderr, "Call stack:\n");
    backtrace_symbols_fd (array, size, STDERR_FILENO);
#endif
}

static bool im_in_main_exit = false, exit_completed = false;
void main_exit (bool show_stack, bool is_error) 
{
    im_in_main_exit = true;

    if (false /* show_stack */) print_call_stack(); //this is useless - doesn't print function names

    buf_test_overflows_all_vbs("exit_on_error");

    url_kill_curl();
    file_kill_external_compressors(); 

    // if we're in ZIP - remove failed genozip file (but don't remove partial failed text file in PIZ - it might be still useful to the user)
    if (primary_command == ZIP && z_file && z_file->name && !flag_reading_reference) {
        char save_name[strlen (z_file->name)+1];
        strcpy (save_name, z_file->name);

        // if we're not the main thread - cancel the main thread before closing z_file, so that the main 
        // thread doesn't attempt to access it (eg. z_file->data_type) and get a segmentation fault.
        if (!arch_am_i_io_thread()) 
            cancel_io_thread(); 

        file_close (&z_file, false); // also frees file->name

        // note: logic to avoid a race condition causing the file not to be removed - if another thread seg-faults
        // because it can't access a z_file filed after z_file is freed, and main_sigsegv_handler aborts
        file_remove (save_name, true);
    }

    exit_completed = true;

    if (is_error)
        abort();
    else
        exit (0);
} 

#ifndef _WIN32
static void main_sigsegv_handler (int sig) 
{
    // note: during exit_on_error, we close z_file which might cause compute threads access z_file fields to 
    // throw a segmentation fault. we don't show it in this case, as we have already displayed the error itself
    if (!im_in_main_exit) 
        fprintf (stderr, "\nError: Segmentation fault\n");

    // busy-wait for exit_on_error to complete before exiting cleanly
    else {
        //fprintf (stderr, "Segmentation fault might appear now - it is not an error - you can safely ignore it\n");
        while (!exit_completed) 
            usleep (10000); // 10 millisec
        exit (1);
    }

    //print_call_stack(); //this is useless - doesn't print function names
    abort();
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
    static int64_t total_uncompressed_len=0, total_compressed_len=0;
    char z_size_str[20], txt_size_str[20], num_lines_str[20], indiv_str[20];

    const unsigned FILENAME_WIDTH = 40;

    const char *head_format   = "\n%5s %11s %10s %10s %6s %s  %*s %s\n";
    const char *foot_format_1 = "\nTotal:            %10s %10s %5uX\n";
    const char *foot_format_2 = "\nTotal:            %10s %10s %5.1fX\n";
    const char *item_format_1 = "%5s %11s %10s %10s %5uX  %s  %s%s%*s %s\n";
    const char *item_format_2 = "%5s %11s %10s %10s %5.1fX  %s  %s%s%*s %s\n";

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
    Md5Hash md5_hash_bound, license_hash, ref_file_md5;
    char created[FILE_METADATA_LEN];
    char ref_file_name[REF_FILENAME_LEN];
    bool success = zfile_get_genozip_header (z_file, NULL, &txt_data_size, &num_samples, &num_lines, 
                                             &md5_hash_bound, created, FILE_METADATA_LEN, &license_hash,
                                             ref_file_name, REF_FILENAME_LEN, &ref_file_md5);
    if (!success) goto finish;

    double ratio = z_file->disk_size ? ((double)txt_data_size / (double)z_file->disk_size) : 0;
    
    str_size (z_file->disk_size, z_size_str);
    str_size (txt_data_size, txt_size_str);
    str_uint_commas (num_lines, num_lines_str);
    str_uint_commas (num_samples, indiv_str);
    
    // TODO: have an option to print ref_file_name and ref_file_md5
    
    bufprintf (evb, &str_buf, ratio < 100 ? item_format_2 : item_format_1, indiv_str, num_lines_str, 
               z_size_str, txt_size_str, ratio, 
               md5_display (md5_hash_bound),
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
                            bool is_last_file)
{
    txtfile_header_initialize();
    
    // get input FILE
    ASSINP0 (z_filename, "Error: z_filename is NULL");

    // we cannot work with a remote genozip file because the decompression process requires random access
    ASSINP (!url_is_url (z_filename), 
            "%s: genozip files must be regular files, they cannot be a URL: %s", global_cmd, z_filename);

    ASSINP (!txt_filename || !url_is_url (txt_filename), 
            "%s: output files must be regular files, they cannot be a URL: %s", global_cmd, txt_filename);

    unsigned fn_len = strlen (z_filename);

    // skip this file if its size is 0
    RETURNW (file_get_size (z_filename),, "Cannot decompress file %s because its size is 0 - skipping it", z_filename);

    if (!txt_filename && (!flag_stdout || flag_bgzip || flag_bcf || flag_bam) && !flag_unbind) {
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
        ASSERT0 (!txt_file || flag_bind, "Error: txt_file is unexpectedly already open"); // note: in bound mode, we expect it to be open for 2nd+ file

        if (!txt_file)  // in bound mode, for second file onwards, txt_file is already open
            txt_file = file_open (txt_filename, WRITE, TXT_FILE, z_file->data_type);
    }
    else if (flag_stdout) { // stdout
        txt_file = file_open_redirect (WRITE, TXT_FILE, z_file->data_type); // STDOUT
    }
    else if (flag_unbind) {
        // do nothing - the component files will be opened by txtfile_genozip_to_txt_header()
    }
    else {
        ABORT0 ("Error: unrecognized configuration for the txt_file");
    }
    
    z_file->basename = file_basename (z_filename, false, FILENAME_STDIN, NULL, 0); // memory freed in file_close
    
    // save flag as it might be modified by PIZ (for example, REF_STORED overriding REF_EXTERNAL)
    SAVE_FLAG (flag_reference);

    // a loop for decompressing all components in unbind mode. in non-unbind mode, it collapses to one a single iteration.
    bool piz_successful;
    unsigned num_components=0;
    do {
        piz_successful = piz_dispatcher (num_components==0, is_last_file);
        if (piz_successful) num_components++;
    } while (flag_unbind && piz_successful); 

    if (!flag_bind && !flag_stdout && !flag_unbind) 
        // don't close the bound file - it will close with the process exits
        // don't close in unbind mode - piz_dispatcher() opens and closes each component
        // don't close stdout - in bound mode, we might still need it for the next file
        file_close (&txt_file, false); 

    file_close (&z_file, false);

    RESTORE_FLAG (flag_reference);

    if (flag_replace && txt_filename && z_filename) file_remove (z_filename, true); 
}

// run the test genounzip after genozip - for the most reliable testing that is nearly-perfectly indicative of actually 
// genounzipping, we create a new genounzip process
static void main_test_after_genozip (char *exec_name, char *z_filename, bool is_last_file)
{
    const char *password = crypt_get_password();

    // is we have a loaded reference and it is no longer needed, unload it now, to free the memory for the testing process
    if (is_last_file) ref_unload_reference (true);

    StreamP test = stream_create (0, 0, 0, 0, 0, 0, 
                                  "To use the --test option",
                                  exec_name, "--decompress", "--test", z_filename,
                                  flag_quiet       ? "--quiet"       : SKIP_ARG,
                                  password         ? "--password"    : SKIP_ARG,
                                  password         ? password        : SKIP_ARG,
                                  flag_show_memory ? "--show-memory" : SKIP_ARG,
                                  flag_show_time   ? "--show-time"   : SKIP_ARG,
                                  threads_str      ? "--threads"     : SKIP_ARG,
                                  flag_reference == REF_EXTERNAL ? "--reference" : SKIP_ARG,
                                  flag_reference == REF_EXTERNAL ? ref_filename  : SKIP_ARG,
                                  threads_str      ? threads_str     : SKIP_ARG,
                                  NULL);

    // wait for child process to finish, so that the shell doesn't print its prompt until the test is done
    int exit_code = stream_wait_for_exit (test);

    primary_command = TEST_AFTER_ZIP; // make exit_on_error NOT delete the genozip file in this case, so its available for debugging
    ASSINP (!exit_code, "genozip test exited with status %d\n", exit_code);
    primary_command = ZIP; // recover in case of more non-concatenated files
}

static void main_genozip (const char *txt_filename, 
                          char *z_filename,
                          bool is_first_file, bool is_last_file,
                          char *exec_name)
{
    license_get(); // ask the user to register if she doesn't already have a license (note: only genozip requires registration - unzip,cat,ls do not)

    ASSINP (!z_filename || !url_is_url (z_filename), 
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

    stats_add_txt_name (txt_name);

    ASSERT0 (flag_bind || !z_file, "Error: expecting z_file to be NULL in non-bound mode");

    // get output FILE
    if (!flag_stdout) {

        if (!z_file) { // skip if we're the second file onwards in bind mode - nothing to do

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

                FREE (basename);
            }

            z_file = file_open (z_filename, flag_pair ? WRITEREAD : WRITE, Z_FILE, txt_file->data_type);
        }
    }
    else if (flag_stdout) { // stdout
#ifdef _WIN32
        // this is because Windows redirection is in text (not binary) mode, meaning Windows edits the output stream...
        ASSINP (isatty(1), "%s: redirecting binary output is not supported on Windows, use --output instead", global_cmd);
#endif
        ASSINP (flag_force || !isatty(1), "%s: you must use --force to output a compressed file to the terminal", global_cmd);

        z_file = file_open_redirect (WRITE, Z_FILE, txt_file->data_type);
    } 
    else ABORT0 ("Error: No output channel");
    
    txt_file->basename = file_basename (txt_filename, false, FILENAME_STDIN, NULL, 0);
    zip_dispatcher (txt_file->basename, is_last_file);

    if (flag_show_stats && is_last_file) stats_show_stats();

    bool remove_txt_file = z_file && flag_replace && txt_filename;

    file_close (&txt_file, !is_last_file);

    if ((is_last_file || !flag_bind) && !flag_stdout && z_file) 
        file_close (&z_file, !is_last_file); 

    if (remove_txt_file) file_remove (txt_filename, true); 

    // test the compression, if the user requested --test
    if (flag_test && (!flag_bind || is_last_file)) main_test_after_genozip (exec_name, z_filename, is_last_file);
}

static void main_list_dir(const char *dirname)
{
    DIR *dir;
    struct dirent *ent;

    dir = opendir (dirname);
    ASSINP (dir, "Error: failed to open directory: %s", strerror (errno));

    int ret = chdir (dirname);
    ASSINP (!ret, "Error: failed to chdir(%s)", dirname);

    while ((ent = readdir(dir))) 
        if (!file_is_dir (ent->d_name))  // don't go down subdirectories recursively
            main_genols (ent->d_name, false, dirname, true);
    
    closedir(dir);    

    ret = chdir ("..");
    ASSINP0 (!ret, "Error: failed to chdir(..)");
}

static void main_warn_if_duplicates (int num_files, char **filenames)
{
    if (num_files <= 1) return; // nothing to do

    # define BASENAME_LEN 256
    char basenames[num_files * BASENAME_LEN];

    for (unsigned i=0; i < num_files; i++)
        file_basename (filenames[i], false, "", &basenames[i*BASENAME_LEN], BASENAME_LEN);

    qsort (basenames, num_files, BASENAME_LEN, (int (*)(const void *, const void *))strcmp);

    for (unsigned i=1; i < num_files; i++) 
        ASSERTW (strncmp(&basenames[(i-1) * BASENAME_LEN], &basenames[i * BASENAME_LEN], BASENAME_LEN), 
                 "Warning: two files with the same name '%s' - if you later unbind with 'genounzip --unbind %s', these files will overwrite each other", 
                 &basenames[i * BASENAME_LEN], out_filename);
}

static DataType main_get_file_dt (const char *filename)
{   
    DataType dt = DT_NONE;

    if (command == ZIP) {
        FileType ft = file_get_stdin_type(); // check for --input option

        if (ft == UNKNOWN_FILE_TYPE) // no --input - get file type from filename
            ft = file_get_type (filename, false);

        dt = file_get_data_type (ft, true);
    }

    else if (command == PIZ) {
        FileType ft = file_get_type (filename, false);
        dt = file_get_dt_by_z_ft (ft);

        // case: we don't know yet what file type this is - we need to read the genozip header to determine
        if (dt == DT_NONE && filename) {
            File *genozip_file = file_open (filename, READ, Z_FILE, DT_NONE);
            zfile_get_genozip_header (genozip_file, &dt, 0,0,0,0,0,0,0,0,0,0);
            file_close (&genozip_file, false);
        }
    }

    return dt;
}

static int main_sort_input_filenames (const void *fn1, const void *fn2)
{
    DataType dt1 = main_get_file_dt (*(char **)fn1);
    DataType dt2 = main_get_file_dt (*(char **)fn2);
    
    bool use_refhash1 = (dt1 == DT_FASTQ);
    bool use_refhash2 = (dt2 == DT_FASTQ);

    // refhash users - at the end of the list
    if (use_refhash1 != use_refhash2) return (int)use_refhash1 - (int)use_refhash2;

    // within refhash users, and within refhash non-users, sort by data type
    if (dt1 != dt2) return (int)dt1 - (int)dt2;

    // within files of the same data type, keep original order
    return fn1 - fn2;
}

// if two filenames differ by one character only, which is '1' and '2', creates a combined filename with "1+2"
static char *main_get_fastq_pair_filename (const char *fn1, const char *fn2)
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

    char *pair_fn = malloc (len+20);
    sprintf (pair_fn, "%.*s1+2%s" FASTQ_GENOZIP_, df, rn1, &rn1[df+1]);
    
    return pair_fn;
}

static void main_load_reference (const char *filename, bool is_first_file)
{
    int old_flag_ref_whole_genome = flag_ref_whole_genome;
    flag_ref_whole_genome = flag_pair ||  // this flag is also used when PIZ reads a stored reference
                            (main_get_file_dt (filename) == DT_FASTQ || main_get_file_dt (filename) == DT_FASTA); 

    // no need to load the reference if genocat just wants to see some sections 
    if (flag_genocat_info_only) return;

    if (flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE) {

        if (is_first_file)
            ref_load_external_reference (false); 

        // Reload the reference in some cases. TODO: eliminate reloads (bug 157)
        else if (!is_first_file && 
            ((old_flag_ref_whole_genome == false && flag_ref_whole_genome == true) ||  // case 1: fastq follows a non-fastq
            (flag_reference == REF_EXT_STORE && !flag_bind))) {             // case 2: --REFERENCE for non-binding files
            ref_unload_reference (true); 
            ref_load_external_reference (false); 
        }
    }
}

static void main_set_flags_from_command_line (int argc, char **argv, bool *is_short)
{
    // process command line options
    while (1) {

        #define _i  {"input-type",    required_argument, 0, 'i'                    }
        #define _I  {"stdin-size",    required_argument, 0, 'I'                    }
        #define _c  {"stdout",        no_argument,       &flag_stdout,           1 }
        #define _d  {"decompress",    no_argument,       &command, PIZ           }
        #define _f  {"force",         no_argument,       &flag_force,            1 }
        #define _h  {"help",          no_argument,       &command, HELP            }
        #define _l  {"list",          no_argument,       &command, LIST            }
        #define _L1 {"license",       no_argument,       &command, LICENSE         } // US spelling
        #define _L2 {"licence",       no_argument,       &command, LICENSE         } // British spelling
        #define _q  {"quiet",         no_argument,       &flag_quiet,            1 }
        #define _Q  {"noisy",         no_argument,       &flag_noisy,            1 }
        #define _DL {"replace",       no_argument,       &flag_replace,          1 }
        #define _V  {"version",       no_argument,       &command, VERSION         }
        #define _z  {"bgzip",         no_argument,       &flag_bgzip,            1 }
        #define _zb {"bam",           no_argument,       &flag_bam,              1 }
        #define _zc {"bcf",           no_argument,       &flag_bcf,              1 }
        #define _m  {"md5",           no_argument,       &flag_md5,              1 }
        #define _t  {"test",          no_argument,       &flag_test,             1 }
        #define _fa {"fast",          no_argument,       &flag_fast,             1 }
        #define _9  {"optimize",      no_argument,       &flag_optimize,         1 } // US spelling
        #define _99 {"optimise",      no_argument,       &flag_optimize,         1 } // British spelling
        #define _9s {"optimize-sort", no_argument,       &flag_optimize_sort,    1 }
        #define _9P {"optimize-PL",   no_argument,       &flag_optimize_PL,      1 }
        #define _9G {"optimize-GL",   no_argument,       &flag_optimize_GL,      1 }
        #define _9g {"optimize-GP",   no_argument,       &flag_optimize_GP,      1 }
        #define _9V {"optimize-VQSLOD", no_argument,     &flag_optimize_VQSLOD,  1 }
        #define _9Q {"optimize-QUAL", no_argument,       &flag_optimize_QUAL,    1 } 
        #define _9f {"optimize-Vf",   no_argument,       &flag_optimize_Vf,      1 }
        #define _9Z {"optimize-ZM",   no_argument,       &flag_optimize_ZM,      1 }
        #define _9D {"optimize-DESC", no_argument,       &flag_optimize_DESC,    1 }
        #define _9S {"optimize-SEQ",  no_argument,       &flag_optimize_SEQ,     1 }
        #define _gt {"gtshark",       no_argument,       &flag_gtshark,          1 } 
        #define _pe {"pair",          no_argument,       &flag_pair,   PAIR_READ_1 } 
        #define _th {"threads",       required_argument, 0, '@'                    }
        #define _u  {"unbind",        no_argument,       &flag_unbind,           1 }
        #define _o  {"output",        required_argument, 0, 'o'                    }
        #define _p  {"password",      required_argument, 0, 'p'                    }
        #define _B  {"vblock",        required_argument, 0, 'B'                    }
        #define _S  {"sblock",        required_argument, 0, 'S'                    }
        #define _r  {"regions",       required_argument, 0, 'r'                    }
        #define _s  {"samples",       required_argument, 0, 's'                    }
        #define _e  {"reference",     required_argument, 0, 'e'                    }
        #define _E  {"REFERENCE",     required_argument, 0, 'E'                    }
        #define _me {"make-reference",no_argument,       &flag_make_reference,   1 }
        #define _g  {"grep",          required_argument, 0, 'g'                    }
        #define _G  {"drop-genotypes",no_argument,       &flag_drop_genotypes,   1 }
        #define _H1 {"no-header",     no_argument,       &flag_no_header,        1 }
        #define _H0 {"header-only",   no_argument,       &flag_header_only,      1 }
        #define _1  {"header-one",    no_argument,       &flag_header_one,       1 }
        #define _GT {"GT-only",       no_argument,       &flag_gt_only,          1 }
        #define _Gt {"gt-only",       no_argument,       &flag_gt_only,          1 }
        #define _fs {"sequential",    no_argument,       &flag_fasta_sequential, 1 }  
        #define _rg {"register",      no_argument,       &flag_register,         1 }
        #define _ss {"show-stats",    no_argument,       &flag_show_stats,       1 } 
        #define _SS {"SHOW-STATS",    no_argument,       &flag_show_stats,       2 } 
        #define _sd {"show-dict",     no_argument,       &flag_show_dict,        1 } 
        #define _d1 {"show-one-dict", required_argument, 0, '\3'                   }
        #define _d2 {"show-dict-one", required_argument, 0, '\3'                   }
        #define _lc {"list-chroms",   no_argument,       &flag_list_chroms,      1 } // identical to --show-one-dict for the CHROM dict of the data type
        #define _sg {"show-gt-nodes", no_argument,       &flag_show_gt_nodes,    1 } 
        #define _s2 {"show-b250",     no_argument,       &flag_show_b250,        1 } 
        #define _s5 {"show-one-b250", required_argument, 0, '\2'                   }
        #define _s6 {"show-b250-one", required_argument, 0, '\2'                   }
        #define _s7 {"dump-one-b250", required_argument, 0, '\5'                   }
        #define _s8 {"dump-b250-one", required_argument, 0, '\5'                   }
        #define _S7 {"dump-one-local",required_argument, 0, '\6'                   }
        #define _S8 {"dump-local-one",required_argument, 0, '\6'                   }
        #define _sa {"show-alleles",  no_argument,       &flag_show_alleles,     1 }
        #define _st {"show-time",     no_argument,       &flag_show_time   ,     1 } 
        #define _sm {"show-memory",   no_argument,       &flag_show_memory ,     1 } 
        #define _sh {"show-headers",  no_argument,       &flag_show_headers,     1 } 
        #define _si {"show-index",    no_argument,       &flag_show_index  ,     1 } 
        #define _Si {"show-ref-index",no_argument,       &flag_show_ref_index,   1 } 
        #define _Sh {"show-ref-hash" ,no_argument,       &flag_show_ref_hash,    1 } 
        #define _sr {"show-gheader",  no_argument,       &flag_show_gheader,     1 }  
        #define _sT {"show-threads",  no_argument,       &flag_show_threads,     1 }  
        #define _sv {"show-vblocks",  no_argument,       &flag_show_vblocks,     1 }  
        #define _sR {"show-reference",no_argument,       &flag_show_reference,   1 }  
        #define _sI {"show-is-set",   required_argument, 0, '~',                   }  
        #define _sA {"show-aliases",  no_argument,       &flag_show_aliases,     1 }  
        #define _dm {"debug-memory",  no_argument,       &flag_debug_memory,     1 }  
        #define _dp {"debug-progress",no_argument,       &flag_debug_progress,   1 }  
        #define _ds {"debug-no-singletons",no_argument,  &flag_debug_no_singletons, 1 }  
        #define _dh {"show-hash",     no_argument,       &flag_show_hash,        1 }  
        #define _00 {0, 0, 0, 0                                                    }

        typedef const struct option Option;
        static Option genozip_lo[]    = { _i, _I, _c, _d, _f, _h, _l, _L1, _L2, _q, _Q, _t, _DL, _V, _z, _zb, _zc, _m, _th, _u, _o, _p, _e, _E,                                    _ss, _SS, _sd, _sT, _d1, _d2, _lc, _sg, _s2, _s5, _s6, _s7, _s8, _S7, _S8, _sa, _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv, _B, _S, _dm, _dp, _dh,_ds, _9, _99, _9s, _9P, _9G, _9g, _9V, _9Q, _9f, _9Z, _9D, _9S, _gt, _pe, _fa,          _rg, _sR, _me, _sA, _sI, _00 };
        static Option genounzip_lo[]  = {         _c,     _f, _h,     _L1, _L2, _q, _Q, _t, _DL, _V, _z, _zb, _zc, _m, _th, _u, _o, _p, _e,                                                  _sd, _sT, _d1, _d2, _lc,      _s2, _s5, _s6, _s7, _s8, _S7, _S8,      _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv,         _dm, _dp,                                                                                                  _sR,      _sA, _sI, _00 };
        static Option genocat_lo[]    = {                 _f, _h,     _L1, _L2, _q, _Q,          _V,                   _th,     _o, _p,         _r, _s, _G, _1, _H0, _H1, _Gt, _GT,          _sd, _sT, _d1, _d2, _lc,      _s2, _s5, _s6, _s7, _s8, _S7, _S8,      _st, _sm, _sh, _si, _Si, _Sh, _sr, _sv,         _dm, _dp,                                                                                    _fs, _g,      _sR,      _sA, _sI, _00 };
        static Option genols_lo[]     = {                 _f, _h,     _L1, _L2, _q,              _V,                                _p, _e,                                                                                                                        _st, _sm,                                       _dm,                                                                                                                           _00 };
        static Option *long_options[] = { genozip_lo, genounzip_lo, genols_lo, genocat_lo }; // same order as ExeType

        // include the option letter here for the short version (eg "-t") to work. ':' indicates an argument.
        static const char *short_options[] = { // same order as ExeType
            "i:I:cdfhlLqQt^Vzm@:o:p:B:S:9KwWFe:E:2uz",  // genozip
            "czfhLqQt^V@:uo:p:me:",                     // genounzip
            "hLVp:qf",                                  // genols
            "hLV@:p:qQ1r:s:H1Go:fg:e:E:"              // genocat
        };

        int option_index = -1;
        int c = getopt_long (argc, argv, short_options[exe_type], long_options[exe_type], &option_index);

        if (c == -1) break; // no more options

        is_short[c] = (option_index == -1); // true if user provided a short option - eg -c rather than --stdout

        switch (c) {
            case PIZ : case LIST : case LICENSE : 
            case VERSION  : case HELP :
                ASSINP (command<0 || command==c, "%s: can't have both -%c and -%c", global_cmd, command, c); 
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
            case 'w' : flag_show_stats    = 1      ; break;
            case 'W' : flag_show_stats    = 2      ; break;
            case 'K' : flag_gtshark       = 1      ; break;
            case '2' : flag_pair    = PAIR_READ_1  ; break;
            case 't' : flag_test          = 1      ; break; 
                       // fall through for genocat -r
            case 'r' : flag_regions       = 1      ; regions_add     (optarg); break;
            case 's' : flag_samples       = 1      ; vcf_samples_add (optarg); break;
            case 'e' : flag_reference     = REF_EXTERNAL  ; ref_set_reference (optarg); break;
            case 'E' : flag_reference     = REF_EXT_STORE ; ref_set_reference (optarg); break;
            case 'm' : flag_md5           = 1      ; break;
            case 'u' : flag_unbind        = 1      ; break;
            case 'G' : flag_drop_genotypes= 1      ; break;
            case 'H' : flag_no_header     = 1      ; break;
            case '1' : flag_header_one    = 1      ; break;
            case '@' : threads_str      = optarg   ; break;
            case 'o' : out_filename     = optarg   ; break;
            case 'g' : flag_grep        = optarg   ; break;
            case '~' : flag_show_is_set = optarg   ; break;
            case '\2' : dict_id_show_one_b250  = dict_id_make (optarg, strlen (optarg)); break;
            case '\5' : mtf_initialize_binary_dump (optarg, &dump_one_b250_dict_id,  "b250");  break;
            case '\6' : mtf_initialize_binary_dump (optarg, &dump_one_local_dict_id, "local"); break;
            case '\3' : dict_id_show_one_dict  = dict_id_make (optarg, strlen (optarg)); break;
            case 'B' : vb_set_global_max_memory_per_vb (optarg); 
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
                    ASSINP (long_options[exe_type][option_index].flag != &command || 
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
}

static void main_process_flags (unsigned num_files, char **filenames, const bool *is_short)
{
    // check for incompatabilities between flags
    #define OT(l,s) is_short[(int)s[0]] ? "-"s : "--"l

    ASSINP (!flag_stdout      || !out_filename,                     "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("output", "o"));
    ASSINP (!flag_stdout      || !flag_unbind,                      "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("unbind", "u"));
    ASSINP (!flag_unbind      || !out_filename,                     "%s: option %s is incompatable with %s", global_cmd, OT("unbind",  "u"), OT("output", "o"));
    ASSINP (!flag_stdout      || !flag_replace,                     "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("replace", "^"));

    // if flag_md5 we need to seek back and update the md5 in the txt header section - this is not possible with flag_stdout
    ASSINP (!flag_stdout      || !flag_md5,                         "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("md5", "m"));
    ASSINP (!flag_test        || !out_filename || command != PIZ,   "%s: option %s is incompatable with %s", global_cmd, OT("output", "o"),  OT("test", "t"));
    ASSINP (!flag_test        || !flag_replace || command != PIZ,   "%s: option %s is incompatable with %s", global_cmd, OT("replace", "^"), OT("test", "t"));
    ASSINP (!flag_test        || !flag_stdout  || command != ZIP,   "%s: option %s is incompatable with %s", global_cmd, OT("stdout", "c"), OT("test", "t"));
    ASSINP (!flag_header_only || !flag_no_header,                   "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), "header-only");
    ASSINP (!flag_no_header   || !flag_header_one,                  "%s: option %s is incompatable with %s", global_cmd, OT("no-header", "H"), OT("header-one", "1"));
    ASSINP (!flag_quiet       || !flag_noisy,                       "%s: option %s is incompatable with %s", global_cmd, OT("quiet", "q"), OT("noisy", "Q"));
    ASSINP (!flag_test        || !flag_optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("test", "t"), OT("optimize", "9"));
    ASSINP (!flag_md5         || !flag_optimize,                    "%s: option %s is incompatable with %s", global_cmd, OT("md5", "m"), OT("optimize", "9"));
    ASSINP (!flag_samples     || !flag_drop_genotypes,              "%s: option %s is incompatable with %s", global_cmd, OT("samples", "s"), OT("drop-genotypes", "G"));
    ASSINP (!flag_test        || !flag_make_reference,              "%s: option %s is incompatable with --make-reference", global_cmd, OT("test", "t"));
    ASSINP (flag_reference != REF_EXTERNAL  || !flag_make_reference,"%s: option %s is incompatable with --make-reference", global_cmd, OT("reference", "e"));
    ASSINP (flag_reference != REF_EXT_STORE || !flag_make_reference,"%s: option %s is incompatable with --make-reference", global_cmd, OT("REFERENCE", "E"));
    ASSINP (flag_reference != REF_EXT_STORE || exe_type != EXE_GENOCAT, "%s: option %s supported only for viewing the reference file itself", global_cmd, OT("REFERENCE", "E"));
    ASSINP (!dump_one_b250_dict_id.num || !dump_one_local_dict_id.num, "%s: option --dump-one-b250 is incompatable with --dump-one-local", global_cmd);

    // some genozip flags are allowed only in combination with --decompress 
    if (exe_type == EXE_GENOZIP && command == ZIP) {
        ASSINP (!flag_bgzip,        "%s: option %s can only be used if --decompress is used too", global_cmd, OT("bgzip", "z"));
        ASSINP (!flag_bam,          "%s: option --flag_bam can only be used if --decompress is used too", global_cmd);
        ASSINP (!flag_bcf,          "%s: option --flag_bcf can only be used if --decompress is used too", global_cmd);
        ASSINP (!flag_unbind,       "%s: option %s can only be used if --decompress is used too", global_cmd, OT("unbind", "u"));
        ASSINP (!flag_show_aliases, "%s: option --flag_show_aliases can only be used if --decompress is used too", global_cmd);
        ASSINP (!flag_show_is_set,  "%s: option --flag_show_is_set can only be used if --decompress is used too", global_cmd);
    }

    // --paired_end: verify an even number of fastq files, --output, and --reference/--REFERENCE
    if (flag_pair) {

        for (unsigned i=0; i < num_files; i++)
            ASSERT (main_get_file_dt (filenames[i]) == DT_FASTQ, "%s: when using %s, all input files are expected to be FASTQ files, but %s is not", global_cmd, OT("pair", "2"), filenames[i]);

        // in case of a flag_pair with 2 files, in which --output is missing, we attempt to figure it out if possible
        if (!out_filename && num_files == 2) 
            out_filename = main_get_fastq_pair_filename (filenames[0], filenames[1]);

        ASSINP (out_filename,       "%s: --output must be specified when using %s", global_cmd, OT("pair", "2"));
        ASSINP (flag_reference,     "%s: either --reference or --REFERENCE must be specified when using %s", global_cmd, OT("pair", "2"));
        ASSINP (num_files % 2 == 0, "%s: when using %s, expecting an even number of FASTQ input files, each consecutive two being a pair", global_cmd, OT("pair", "2"));
    }

    if (flag_gtshark) stream_abort_if_cannot_run ("gtshark", "To use the --gtshark option"); 

    // deal with final setting of flag_quiet that suppresses warnings 
    
    // don't show progress or warning when outputing to stdout (note: we are "quiet" even if output doesn't go to the terminal
    // because often it will be piped and ultimately go the terminal)
    if (flag_stdout) flag_quiet=true; 
    
    // don't show progress for flags that output throughout the process. no issue with flags that output only in the end
    if (flag_show_dict || flag_show_gt_nodes || flag_show_b250 || flag_show_headers || flag_show_threads ||
        dict_id_show_one_b250.num || dict_id_show_one_dict.num || flag_show_reference ||
        flag_show_alleles || flag_show_vblocks || (flag_show_index && command==PIZ))
        flag_quiet=true; // don't show progress

    // override these ^ if user chose to be --noisy
    if (flag_noisy) flag_quiet=false;

    if (flag_test) flag_md5=true; // test requires md5

    // default values, if not overridden by the user
    if (!flag_vblock) vb_set_global_max_memory_per_vb (flag_fast ? TXT_DATA_PER_VB_FAST : TXT_DATA_PER_VB_DEFAULT); 
    if (!flag_sblock) vcf_zip_set_global_samples_per_block (VCF_SAMPLES_PER_VBLOCK); 

    // --make-reference implies --md5 --B1 (unless --vblock says otherwise), and not encrypted. 
    // in addition, txtfile_read_vblock() limits each VB to have exactly one contig.
    if (flag_make_reference) {
        ASSINP (!crypt_have_password(), "%s: option --make-reference is incompatable with %s", global_cmd, OT("password", "p"));

        flag_md5 = true;
        if (!flag_vblock) vb_set_global_max_memory_per_vb("1");

        ASSINP (num_files <= 1, "%s: you can specify only one FASTA file when using --make-reference.\n"
                "To create a reference from multiple FASTAs use something like this:\n"
                "cat *.fa | %s --make-reference --input fasta --output myref.ref.genozip -", global_cmd, global_cmd);
    }

    // if --optimize was selected, all optimizations are turned on
    if (flag_optimize)
        flag_optimize_sort = flag_optimize_PL = flag_optimize_GL = flag_optimize_GP   = flag_optimize_VQSLOD = 
        flag_optimize_QUAL = flag_optimize_Vf = flag_optimize_ZM = flag_optimize_DESC = flag_optimize_SEQ = true;
    
    // if any optimization flag is on, we turn on flag_optimize
    if (flag_optimize_sort || flag_optimize_PL || flag_optimize_GL || flag_optimize_GP   || flag_optimize_VQSLOD ||
        flag_optimize_QUAL || flag_optimize_Vf || flag_optimize_ZM || flag_optimize_DESC || flag_optimize_SEQ)
        flag_optimize = true;

    if (flag_fast) flag_optimize_SEQ = false; // if --fast, SEQ is compressed with ACGT (if we left flag_optimize_SEQ==true, it would be set to LZMA and fast will change it to BZ2 - very lousy compression for SEQ data)

    // if using the -o option - check that we don't have duplicate filenames (even in different directory) as they
    // will overwrite each other if extracted with --unbind
    if (command == ZIP && out_filename && !flag_quiet) main_warn_if_duplicates (num_files, filenames);

    flag_multiple_files = (num_files > 1);

    flag_bind = (command == ZIP) && (out_filename != NULL) && (num_files > 1);

    // cases where genocat is used to view some information, but not the file contents
    flag_genocat_info_only = exe_type == EXE_GENOCAT &&
                             (flag_show_dict || flag_show_b250 || flag_list_chroms || dict_id_show_one_dict.num ||
                              flag_show_index || dump_one_local_dict_id.num || dump_one_b250_dict_id.num);

    ASSINP (num_files <= 1 || flag_bind || !flag_show_stats, "%s: --show-stats can only work on one file at time", global_cmd);
}

void TEST()
{
}

int main (int argc, char **argv)
{
    //TEST();exit(0);
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

    global_cmd = file_basename (argv[0], true, "(executable)", NULL, 0); // global var

    bool is_short[256] = { 0 }; // indexed by character of short option.
    main_set_flags_from_command_line (argc, argv, is_short);

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        if (exe_type == EXE_GENOLS) command = LIST; // genols can be run without arguments
        
        // genozip with no input filename, no output filename, and no output or input redirection 
        // note: in docker stdin is a pipe even if going to a terminal. so we show the help even if
        // coming from a pipe. the user must use "-" to redirect from stdin
        else if (command == -1 && optind == argc && !out_filename && 
                 (isatty(0) || arch_am_i_in_docker()) && isatty(1)) {
            // case: --register
            if (flag_register) 
                license_get();

            // case: requesting to display the reference: genocat --reference <ref-file> and optionally --regions
            if (exe_type == EXE_GENOCAT && flag_reference) 
                ref_load_external_reference (true);

            // otherwise: show help
            else
                main_print_help (false);

            return 0;
        }

        else if (exe_type == EXE_GENOUNZIP) command = PIZ;
        else if (exe_type == EXE_GENOCAT) { command = PIZ; flag_stdout = !out_filename ; }
        else command = ZIP; // default 
    }

    unsigned num_files = argc - optind;

    main_process_flags (num_files, &argv[optind], is_short);

    // sort files by data type to improve VB re-using, and refhash-using files in the end to improve reference re-using
    qsort (&argv[optind], num_files, sizeof (argv[0]), main_sort_input_filenames);
    
    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (threads_str) {
        int ret = sscanf (threads_str, "%u", &global_max_threads);
        ASSINP (ret == 1 && global_max_threads >= 1, "%s: %s requires an integer value of at least 1", global_cmd, OT("threads", "@"));
    }
    else global_max_threads = arch_get_num_cores();
    
    // handle call commands except for ZIP, PIZ or LIST
    if (command == VERSION) { main_print_version();   return 0; }
    if (command == LICENSE) { license_display();      return 0; }
    if (command == HELP)    { main_print_help (true); return 0; }

    primary_command = command; 
    
    for (unsigned file_i=0; file_i < MAX (num_files, 1); file_i++) {

        char *next_input_file = optind < argc ? argv[optind++] : NULL;  // NULL means stdin
        
        if (next_input_file && !strcmp (next_input_file, "-")) next_input_file = NULL; // "-" is stdin too

        ASSINP (next_input_file || command != PIZ, 
                "%s: filename(s) required (redirecting from stdin is not possible)", global_cmd);

        ASSERTW (next_input_file || !flag_replace, "%s: ignoring %s option", global_cmd, OT("replace", "^")); 

        main_load_reference (next_input_file, !file_i);
        
        switch (command) {
            case ZIP  : main_genozip (next_input_file, out_filename, file_i==0, !next_input_file || file_i==num_files-1, argv[0]); break;
            case PIZ  : main_genounzip (next_input_file, out_filename, file_i==num_files-1); break;           
            case LIST : main_genols (next_input_file, false, NULL, false); break;
            default   : ABORT ("%s: unrecognized command %c", global_cmd, command);
        }

        if (flag_pair) flag_pair = 3 - flag_pair; // alternate between PAIR_READ_1 and PAIR_READ_2
    }

    // if this is "list", finalize
    if (command == LIST) main_genols (NULL, true, NULL, false);

    return 0;
}
