// ------------------------------------------------------------------
//   genozip.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <getopt.h>
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>

#include "genozip.h"
#include "text_help.h"
#include "version.h" // automatically incremented by the make when we create a new distribution
#include "zfile.h"
#include "zip.h"
#include "piz.h"
#include "crypt.h"
#include "file.h"
#include "filename.h"
#include "vblock.h"
#include "url.h"
#include "strings.h"
#include "stats.h"
#include "arch.h"
#include "license.h"
#include "tip.h"
#include "reference.h"
#include "refhash.h"
#include "random_access.h"
#include "codec.h"
#include "threads.h"
#include "recon_plan_io.h"
#include "bases_filter.h"
#include "genols.h"
#include "tar.h"
#include "gencomp.h"

// globals - set it main() and never change
rom global_cmd = NULL; 
ExeType exe_type;

// primary_command vs command: primary_command is what the user typed on the command line. command is what is 
// running now - for example, when ZIP is unzipping a reference, primary_command=ZIP and command=PIZ
CommandType command = NO_COMMAND, primary_command = NO_COMMAND; 

uint32_t global_max_threads = DEFAULT_MAX_THREADS; 
static bool tip_printed = false;

#define MAIN(format, ...) ({ if (flag.debug_top) { progress_newline(); fprintf (stderr, "%s[%u]: ", command_name(), getpid()); fprintf (stderr, (format), __VA_ARGS__); fprintf (stderr, "\n"); } })
#define MAIN0(string)     ({ if (flag.debug_top) { progress_newline(); fprintf (stderr, "%s[%u]: %s\n", command_name(), getpid(), string); } })

rom command_name (void) // see CommandType
{
    switch (command) {
        case ZIP        : return "ZIP";
        case PIZ        : return "PIZ";
        case LIST       : return "LIST";
        case VERSION    : return "VERSION";
        case HELP       : return "HELP";
        case LICENSE    : return "LICENSE";
        case NO_COMMAND : return "NO_COMMAND";
        default         : return "Invalid command";
    }
}

void main_exit (bool show_stack, bool is_error) 
{
    buf_set_cleanup_on_exit();  

    // finish dumping reference and/or refhash to cache
    if (!is_error) {
        ref_create_cache_join (gref, false);
        ref_create_cache_join (prim_ref, false);
        refhash_create_cache_join(false);
        
        version_print_notice_if_has_newer();

        if (flag.show_time && !flag.show_time[0]) { // show-time without the optional parameter 
            profiler_add (evb);
            profiler_print_report();
        }

        else if (IS_ZIP || flag.check_latest/*PIZ - test after compress*/)
            tip_print();
    }

    if (is_error && flag.debug_threads)
        threads_write_log (true);
        
    if (show_stack) 
        threads_print_call_stack(); // this works ok on mac, but does not print function names on Linux (even when compiled with -g)

    buf_test_overflows_all_vbs (is_error ? "exit_on_error" : "on_exit");

    if (flag.log_filename) {
        iprintf ("%s - execution ended %s\n", str_time().s, is_error ? "with an error" : "normally");
        fclose (info_stream);
    } 

    if (is_error) {
        close (1);   // prevent other threads from outputting to terminal (including buffered output), obscuring our error message
        close (2);
        url_kill_curl();  /* <--- BREAKPOINT BRK */
        file_kill_external_compressors(); 
    
        // cancel all other threads before closing z_file, so other threads don't attempt to access it 
        // (eg. z_file->data_type) and get a segmentation fault.
        threads_cancel_other_threads();
}

    // if we're in ZIP - rename failed genozip file (but not in PIZ - as its not expected to change the file)
    if (primary_command == ZIP && z_file && z_file->name && !flag_loading_auxiliary) {
        STRli (save_name, strlen (z_file->name)+1);

        if (z_file && z_file->name) 
            strcpy (save_name, z_file->name);

        file_close (&z_file, false, false); // also frees file->name

        // note: logic to avoid a race condition causing the file not to be removed - if another thread seg-faults
        // because it can't access a z_file filed after z_file is freed, and threads_sigsegv_handler aborts
        if (is_error && !flag.debug_or_test && file_exists (save_name)) 
            file_remove (save_name, true);
    }

    if (is_error) // call after canceling the writing threads 
        file_put_data_abort();

    fflush (stdout);
    fflush (stderr);
    
    if (!is_error) MAIN0 ("Exiting : Success"); 

    exit (is_error ? EXIT_GENERAL_ERROR : EXIT_OK);
} 

static void main_print_help (bool explicit)
{
    static rom *texts[NUM_EXE_TYPES] = { help_genozip, help_genounzip, help_genols, help_genocat }; // same order as ExeType
    static unsigned sizes[NUM_EXE_TYPES] = {sizeof(help_genozip), sizeof(help_genounzip), sizeof(help_genols), sizeof(help_genocat) };

    if (flag.help && !strcmp (flag.help, "genozip")) 
        str_print_text (help_genozip,      ARRAY_LEN(help_genozip),      "                          ",  "\n", NULL, 0);

    else if (flag.help && !strcmp (flag.help, "genounzip")) 
        str_print_text (help_genounzip,    ARRAY_LEN(help_genounzip),    "                          ",  "\n", NULL, 0);

    else if (flag.help && !strcmp (flag.help, "genols")) 
        str_print_text (help_genols,       ARRAY_LEN(help_genols),       "                          ",  "\n", NULL, 0);

    else if (flag.help && !strcmp (flag.help, "genocat")) 
        str_print_text (help_genocat,      ARRAY_LEN(help_genocat),      "                          ",  "\n", NULL, 0);

    else if (flag.help && !strcmp (flag.help, "attributions")) 
        str_print_text (help_attributions, ARRAY_LEN(help_attributions), "                          ",  "\n", NULL, 0);

    else if (flag.help && !strcmp (flag.help, "input")) 
        iprintf ("Supported file types for --input:\n%s\n", file_compressible_extensions (false));
    
    else 
        str_print_text (texts[exe_type], sizes[exe_type] / sizeof(char*),  "                     ",  "\n", NULL, 0);
    
    str_print_text (help_footer, ARRAY_LEN(help_footer), "", "\n", NULL, 0);

    // in Windows, we ask the user to click a key - this is so that if the user double clicks on the EXE
    // from Windows Explorer - the terminal will open and he will see the help
    if (flag.is_windows && !explicit && isatty(0) && isatty (1)) {
        printf ("Press Enter to continue...\n");
        getc(stdin);
    }
}

static void main_print_version()
{
    iprintf ("version=%s distribution=%s\n", GENOZIP_CODE_VERSION, arch_get_distribution());  
}

static void main_genounzip (rom z_filename, rom txt_filename, int z_file_i, unsigned n_z_files, bool is_last_z_file)
{
    if (n_z_files > 1) MAIN ("main_genounzip (%u/%u): %s", z_file_i+1, n_z_files, z_filename);
    else               MAIN ("main_genounzip: %s", z_filename);

    // save flag as it might be modified - so that next file has the same flags
    SAVE_FLAGS;

    static bool is_first_z_file = true;
    
    // get input FILE
    ASSERTNOTNULL (z_filename);

    // we cannot work with a remote genozip file because the decompression process requires random access
    ASSINP (!url_is_url (z_filename), 
            "genozip files must be regular files, they cannot be a URL: %s", z_filename);

    ASSINP (!txt_filename || !url_is_url (txt_filename), 
            "output files must be regular files, they cannot be a URL: %s", txt_filename);

    z_file = file_open (z_filename, READ, Z_FILE, DT_NONE);    

    if (flag.validate) {
        file_close (&z_file, false, false);
        return; // we're only validating that z_file is valid
    }

    // read the genozip header:
    // 1) verify the data type deduced from the file name, or set the data type if it wasn't deduced
    // 2) if an external reference is not specified, check if the file needs one, and if it does - set it from the header
    // 3) identify skip cases (DT_NONE returned) - empty file, unzip of a reference
    // 4) reset flag.unbind if file contains only one component
    TEMP_FLAG (genocat_no_reconstruct, true); // so zfile_read_genozip_header doesn't fail on a reference file before flags_update_piz_one_file is run
    if (!z_file->file || !zfile_read_genozip_header (0)) goto done; 
    RESTORE_FLAG(genocat_no_reconstruct);

    // case: reference not loaded yet bc --reference wasn't specified, and we got the ref name from zfile_read_genozip_header()   
    if (IS_REF_EXTERNAL && !ref_is_external_loaded(gref)) {
        ASSINP0 (!IS_REF_EXTERNAL || !flag.show_ref_seq, "--show-ref-seq cannot be used on a file that requires a reference file: use genocat --show-ref-seq on the reference file itself instead");

        if (!flag.genocat_no_reconstruct || 
            (flag.collect_coverage && Z_DT(FASTQ))) { // in collect_coverage with FASTQ we read the non-data sections of the reference
            SAVE_VALUE (z_file); // actually, read the reference first
            ref_load_external_reference (gref, NULL);
            RESTORE_VALUE (z_file);
        }
    }

    flags_update_piz_one_file (z_file_i);
   
    // set txt_filename from genozip file name (inc. extensions if translating or --bgzf)
    if (!txt_filename && !flag.to_stdout && !flag.unbind) 
        txt_filename = txtfile_piz_get_filename (z_filename, "", true);

    // open output txt file (except if unbinding - see comment below)
    // note: we need the txt_file object for recontructing, even if not writing (the disk file won't be opened if flag.no_writer)
    if ((txt_filename && !flag.unbind) || (flag.to_stdout && !flag.genocat_no_reconstruct)) {
        ASSERT0 (!txt_file, "txt_file is unexpectedly already open"); // note: in bound mode, we expect it to be open for 2nd+ file
        txt_file = file_open (txt_filename, WRITE, TXT_FILE, flag.out_dt);
    }
    else if (flag.unbind || flag.no_writer || (txt_filename && flag.one_component >= 2)) {
        // do nothing - in unbind, the component files will be opened by txtheader_piz_read_and_reconstruct()
    }
    else 
        ABORT0 ("Error: unrecognized configuration for the txt_file");

    Dispatcher dispatcher = piz_z_file_initialize();  

    // generate txt_file(s)
    if (dispatcher) {
        if (!flag.unbind) {    // single txt_file
            piz_one_txt_file (dispatcher, is_first_z_file, is_last_z_file, COMP_NONE, COMP_NONE);
            file_close (&txt_file, flag.index_txt, !is_last_z_file); 
        }

        else if (Z_DT(FASTQ))  // paired FASTQ
            for (CompIType comp_i=FQ_COMP_R1; comp_i <= FQ_COMP_R2; comp_i++) {
                piz_one_txt_file (dispatcher, is_first_z_file, is_last_z_file, comp_i, comp_i);
                file_close (&txt_file, flag.index_txt, comp_i==FQ_COMP_R1 || !is_last_z_file); 
            }

        else if (Z_DT(SAM) && flag.deep)    // Deep: one SAM/BAM and one or more FASTQ
            for (CompIType comp_i=SAM_COMP_MAIN; comp_i < z_file->num_txt_files + 2; comp_i++) {
                if (comp_i==SAM_COMP_MAIN) {
                    piz_one_txt_file (dispatcher, is_first_z_file, is_last_z_file, SAM_COMP_MAIN, SAM_COMP_DEPN);
                    comp_i += 2;
                }
                else {
                    flag.out_dt = DT_FASTQ;
                    piz_one_txt_file (dispatcher, is_first_z_file, is_last_z_file, comp_i, comp_i);
                }

                file_close (&txt_file, flag.index_txt, comp_i < z_file->num_txt_files + 1 || !is_last_z_file); // skip cleanup in final component, except if there are more z_files
            }

        else
            ABORT0 ("Invalid unbinding mode");
    }

    tip_dt_encountered (z_file->data_type);

    file_close (&z_file, false, false);
    is_first_z_file = false;

    if (digest_piz_has_it_failed()) exit_on_error (false); // error message already displayed in piz_verify_digest_one_txt_file
    
    // case --replace: now that the file was reconstructed, we can remove the genozip file
    if (flag.replace && (txt_filename || flag.unbind) && z_filename) file_remove (z_filename, true); 

done:
    RESTORE_FLAGS;
}

// run the test genounzip after genozip - for the most reliable testing that is nearly-perfectly indicative of actually 
// genounzipping, we create a new genounzip process
static void main_test_after_genozip (rom z_filename, DataType z_dt, bool is_last_txt_file, bool is_chain)
{
    rom password = crypt_get_password();

    // finish dumping reference and/or refhash to cache before destroying it...
    ref_create_cache_join (gref, true);
    ref_create_cache_join (prim_ref, true);
    refhash_create_cache_join (false);

    // On Windows and Mac that usually have limited memory, if ZIP consumed more than 2GB, free memory before PIZ. 
    // Note: on Windows, freeing memory takes considerable time.
    if ((flag.is_windows || flag.is_mac || flag.is_wsl) && buf_get_memory_usage () > (1ULL<<31)) {
        ref_destroy_reference (gref, true); // on Windows I observed a race condition: if we unmap mapped memory here, and remap it in the test process, and the system is very slow due to low memory, then "MapViewOfFile" in the test processs will get "Access is Denied". That's why destroy_only_if_not_mmap=true.
        ref_destroy_reference (prim_ref, true); 
        refhash_destroy (true);
        kraken_destroy();
        chain_destroy();
        vb_destroy_pool_vbs (POOL_MAIN);
    }

    rom exec_path = arch_get_executable();

    if (!is_last_txt_file || flag.is_windows) {

        StreamP test = stream_create (NULL, 0, 0, 0, 0, 0, 0,
                                      "To use the --test option",
                                      exec_path, 
                                      "--decompress", 
                                      "--test", 
                                      z_filename,
                                      flag.quiet         ? "--quiet"         : SKIP_ARG,
                                      password           ? "--password"      : SKIP_ARG,
                                      password           ? password          : SKIP_ARG,
                                      flag.show_digest   ? "--show-digest"   : SKIP_ARG,
                                      flag.log_digest    ? "--log-digest"    : SKIP_ARG,
                                      flag.show_memory   ? "--show-memory"   : SKIP_ARG,
                                      flag.show_time     ? "--show-time"     : SKIP_ARG,
                                      flag.threads_str   ? "--threads"       : SKIP_ARG,
                                      flag.threads_str   ? flag.threads_str  : SKIP_ARG,
                                      flag.xthreads      ? "--xthreads"      : SKIP_ARG,
                                      flag.show_alleles  ? "--show-alleles"  : SKIP_ARG,
                                      flag.show_aligner  ? "--show-aligner"  : SKIP_ARG,
                                      flag.show_threads  ? "--show-threads"  : SKIP_ARG,
                                      flag.debug_threads ? "--debug-threads" : SKIP_ARG,
                                      flag.echo          ? "--echo"          : SKIP_ARG,
                                      flag.verify_codec  ? "--verify-codec"  : SKIP_ARG,
                                      flag.debug_lines   ? "--debug-lines"   : SKIP_ARG,
                                      flag.show_buddy    ? "--show-buddy"    : SKIP_ARG,
                                      flag.no_tip        ? "--no-tip"        : SKIP_ARG,
                                      flag.debug_latest  ? "--debug-latest"  : SKIP_ARG,
                                      flag.license_filename           ? "--licfile"            : SKIP_ARG,
                                      flag.license_filename           ? flag.license_filename  : SKIP_ARG,
                                      is_last_txt_file && !flag.debug ? "--check-latest"       : SKIP_ARG,
                                      IS_REF_EXTERNAL && !is_chain    ? "--reference"          : SKIP_ARG, // normal pizzing of a chain file doesn't require a reference
                                      IS_REF_EXTERNAL && !is_chain    ? ref_get_filename(gref) : SKIP_ARG, 
                                      NULL);
                                      // ↓↓↓ Don't forget to add below too ↓↓↓

        tip_printed = is_last_txt_file; // --check-latest sent to PIZ causes it to print the tip
        
        // wait for child process to finish, so that the shell doesn't print its prompt until the test is done
        int exit_code = stream_wait_for_exit (test);

        TEMP_VALUE (primary_command, TEST_AFTER_ZIP); // make exit_on_error NOT delete the genozip file in this case, so its available for debugging
        ASSERT (!exit_code, "%s: test exited with status: \"%s\"\n", global_cmd, exit_code_name (exit_code)); // exit with error status 
        RESTORE_VALUE (primary_command); // recover in case of more non-concatenated files
    
        FREE (exec_path);
    }

    // case: nothing more to do on the compression side - run test replacing the current process, without forking
    // note: in Windows, we can't use execv because it creates a new process with a new pid, while the current process exits
    // and returns an exit code to the shell.
    else {
      rom argv[32]; 
        int argc = 0;
                                argv[argc++] = exec_path;
                                argv[argc++] = "--decompress";
                                argv[argc++] = "--test";
                                argv[argc++] = "--check-latest";
                                argv[argc++] = z_filename;
        if (flag.quiet)         argv[argc++] = "--quiet";
        if (password)         { argv[argc++] = "--password"; 
                                argv[argc++] = password; }
        if (flag.show_digest)   argv[argc++] = "--show-digest";
        if (flag.show_memory)   argv[argc++] = "--show-memory";
        if (flag.show_time)     argv[argc++] = "--show-time";
        if (flag.threads_str) { argv[argc++] = "--threads"; argv[argc++] = flag.threads_str; }
        if (flag.xthreads)      argv[argc++] = "--xthreads";
        if (flag.show_alleles)  argv[argc++] = "--show-alleles";
        if (flag.show_aligner)  argv[argc++] = "--show-aligner";
        if (flag.show_threads)  argv[argc++] = "--show-threads";
        if (flag.debug_threads) argv[argc++] = "--debug-threads";
        if (flag.echo)          argv[argc++] = "--echo";
        if (flag.verify_codec)  argv[argc++] = "--verify-codec";
        if (flag.debug_lines)   argv[argc++] = "--debug-lines";
        if (flag.show_buddy)    argv[argc++] = "--show-buddy";
        if (flag.no_tip)        argv[argc++] = "--no-tip";
        if (flag.debug_latest)  argv[argc++] = "--debug-latest";
        if (flag.license_filename) { argv[argc++] = "--licfile"; argv[argc++] = flag.license_filename; }
        if (is_last_txt_file && !flag.debug) 
                                argv[argc++] = "--check-latest";
        if (IS_REF_EXTERNAL && !is_chain) 
                              { argv[argc++] = "--reference"; 
                                argv[argc++] = ref_get_filename(gref); }
        argv[argc] = NULL;
                                // ↑↑↑ Don't forget to add above too ↑↑↑

        if (flag.show_time && !flag.show_time[0]) { // show-time without the optional parameter 
            profiler_add (evb);
            profiler_print_report();
        }

        execv (exec_path, (char **)argv);
        ABORT ("Failed to run the executable \"%s\" for testing the compression: %s", exec_path, strerror (errno));
    }
}

static void main_genozip (rom txt_filename, 
                          rom next_txt_filename,      // ignored unless we are of pair_1 in a --pair
                          unsigned txt_file_i,        // 0-based
                          unsigned n_txt_files,
                          bool is_last_user_txt_file) // very last file in this execution 
{
    if (n_txt_files > 1) MAIN ("main_genozip (%u/%u): %s", txt_file_i+1, n_txt_files, txt_filename ? txt_filename : "stdin");
    else                 MAIN ("main_genozip: %s", txt_filename ? txt_filename : "stdin");
    
    SAVE_FLAGS;

    ASSINP (!flag.out_filename || !url_is_url (flag.out_filename), 
            "output files must be regular files, they cannot be a URL: %s", flag.out_filename);

    // get input file
    if (!txt_file)  // open the file - possibly already open from main_load_reference
        txt_file = file_open (txt_filename, READ, TXT_FILE, DT_NONE); 

    // skip this file if its size is 0
    if (!txt_file) {
        if (tar_is_tar()) {
            tar_copy_file (txt_filename);
            WARN ("Copied %s to the tar file", txt_filename);
        }    
        else 
            WARN ("Cannot compress file %s because its size is 0 - skipping it", txt_filename);
        return;
    }

    if (!txt_file->file) {
        file_close (&txt_file, false, false);
        goto done; // this is the case where multiple files are given in the command line, but this one is not compressible (eg its a .genozip file) - we skip it
    }

    bool has_zfile = !!z_file;

    // get output FILE
    if (!z_file) { // skip if we're the second file onwards in bind mode, or pair_2 in unbound list of pairs - nothing to do        
        
        rom z_filename = flag.out_filename                 ? filename_z_by_flag() // given with --output
                       : (flag.pair && !flag.out_filename) ? filename_z_pair (txt_filename, next_txt_filename, false) // first file in a FASTQ pair
                       : (flag.deep && !flag.out_filename) ? filename_z_deep (txt_file->name) // SAM/BAM file in --deep
                       :                                     filename_z_normal (txt_file->name, txt_file->data_type, txt_file->type);

        z_file = file_open (z_filename, flag.pair ? WRITEREAD : WRITE, Z_FILE, txt_file->data_type);
    }

    if (flag.zip_comp_i == COMP_MAIN) {
        stats_add_txt_name (txt_name); // add txt_name (inluding stdin) to stats data stored in z_file
        flags_update_zip_one_file();
    }
    
    ASSERT0 (flag.bind || !has_zfile, "Unexpectedly, z_file was already open despite BIND_NONE");

    z_file->z_closes_after_me = (flag.bind == BIND_FQ_PAIR) ? (flag.zip_comp_i == FQ_COMP_R2)
                              : (flag.bind == BIND_DEEP)    ? is_last_user_txt_file
                              :                               true;

    zip_one_file (txt_file->basename, is_last_user_txt_file);

    // increment zip_comp_i for next file 
    flag.zip_comp_i = (flag.bind == BIND_FQ_PAIR) ? (flag.zip_comp_i + 1) % 2
                    : (flag.bind == BIND_DEEP)    ? ((SamComponentType[]){ SAM_COMP_MAIN, SAM_COMP_FQ00, SAM_COMP_FQ01, COMP_MAIN })[txt_file_i + 1]
                    :                               COMP_MAIN;

    tip_dt_encountered (z_file->data_type);

    if (flag.show_stats && z_file->z_closes_after_me)
        stats_display();

    bool remove_txt_file = z_file && flag.replace && txt_filename;

    file_close (&txt_file, false, !(is_last_user_txt_file && z_file->z_closes_after_me));  // no need to waste time closing the last file, the process termination will do that

    bool is_chain = (Z_DT(CHAIN));
    
    // close the file if its an open disk file AND we need to close it
    if (!flag.to_stdout && z_file && z_file->z_closes_after_me) {

        // take over the name 
        rom z_filename = z_file->name; 
        z_file->name = NULL;
        DataType z_dt = z_file->data_type; 

        file_close (&z_file, false, !is_last_user_txt_file); 

        license_one_file_compressed(z_dt);

        // test the compression, if the user requested --test
        if (flag.test) 
            main_test_after_genozip (z_filename, z_dt, is_last_user_txt_file, is_chain);

        FREE (z_filename); // allocated in filename_z_*
    }

    // remove after test (if any) was successful
    if (remove_txt_file) {
        
        // add file to remove_list, don't actually remove it yet
        static Buffer remove_list = { .name = "remove_list" };
        buf_append_string (evb, &remove_list, txt_filename);
        remove_list.len++;   // include the \0 separator added by buf_append_string
        remove_list.count++; // count files

        // case: z_file has closed and tested if needed - remove all files included in this z_file 
        if (!z_file) {
            str_split (remove_list.data, remove_list.len-1, remove_list.count, '\0', rm_file, true); // -1 to remove last \0

            for (unsigned i=0; i < n_rm_files; i++)
                file_remove (rm_files[i], true); 
            
            buf_free (remove_list);
        }
    }

done: 
    // if next file is a fresh file, restore flags
    if (flag.zip_comp_i == COMP_MAIN) { // next file
        SAVE_FLAG (data_modified); // propagate up
        SAVE_FLAG (aligner_available);
        SAVE_FLAG (no_tip);
        SAVE_FLAG (bind);
        
        RESTORE_FLAGS; // overwrite all flags with saved flag struct

        if (flag.bind) RESTORE_FLAG (data_modified); 
        RESTORE_FLAG (aligner_available);
        RESTORE_FLAG (no_tip);
        RESTORE_FLAG (bind);
        RESTORE_FLAG (bind);
        flag.zip_comp_i = COMP_MAIN; // re-assign in case RESTORE_FLAG modified
    }
} 

static inline DataType main_get_file_dt (rom filename)
{   
    switch (command) {
        case ZIP: return txtfile_get_file_dt (filename);
        case PIZ: return zfile_get_file_dt (filename);
        default : return DT_NONE;
    }
}

static int main_sort_input_filenames (const void *fn1, const void *fn2)
{
    DataType dt1 = main_get_file_dt (*(char **)fn1);
    DataType dt2 = main_get_file_dt (*(char **)fn2);

    int sizeof_vb1 = ((dt1 != DT_NONE && dt_props[dt1].sizeof_vb) ? dt_props[dt1].sizeof_vb : def_vb_size)(dt1);
    int sizeof_vb2 = ((dt2 != DT_NONE && dt_props[dt2].sizeof_vb) ? dt_props[dt2].sizeof_vb : def_vb_size)(dt2);

    // sort by VB type (assuming their size is unique)
    if (sizeof_vb1 != sizeof_vb2) return sizeof_vb1 - sizeof_vb2;

    // within files of the same sizeof_vb, keep original order (note: we assume the order of filename pointers is consistent
    // with their order on the command line. However, this isn't guaranteed. TODO: improve this by having a struct with the explict order in it)
    return *(char **)fn1 - *(char **)fn2;
}

static void main_get_filename_list (unsigned num_files, char **filenames,  // in 
                                    unsigned *num_z_files, BufferP fn_buf) // out
{
    if (!num_files && !flag.files_from) {
        flags_update (0, NULL);
        return; // no files
    }

    // add names from command line
    buf_append (evb, *fn_buf, char *, filenames, num_files, NULL);

    // set argv to file names
    if (flag.files_from) {
        file_split_lines (flag.files_from, "files-from", true);
        
        // handle case of "tar xvf myfile.tar |& genounzip --files-from -" when using bsdtar:
        // the output has an added "x " prefix eg "x junk.sam.genozip"
        if (n_lines && lines[0][0]=='x' && lines[0][1]==' ')
            for (int i=0; i < n_lines; i++)
                lines[i] += 2;

        buf_append (evb, *fn_buf, char *, lines, n_lines, NULL);
    }

    // expand directories if --subdirs 
    for (unsigned i=0; i < fn_buf->len; i++) { // len increases during the loop as we add more files which may themselves be directories

        char *fn = *B(char *, *fn_buf, i);
        unsigned fn_len = strlen (fn);
        if (fn[fn_len-1] == '/') 
            fn[--fn_len] = 0; // change "mydir/" to "mydir"

        if (!file_is_dir (fn)) continue;

        if (flag.subdirs) {
            DIR *dir = opendir (fn);

            ASSERT (dir, "failed to open directory %s: %s", fn, strerror (errno));

            struct dirent *ent;
            while ((ent = readdir(dir))) {
                if (ent->d_name[0] =='.') continue; // skip . ..  and hidden files

                char *ent_fn = MALLOC (strlen(ent->d_name) + fn_len + 2);
                sprintf (ent_fn, "%.*s/%s", fn_len, fn, ent->d_name);
                buf_append (evb, *fn_buf, char *, &ent_fn, 1, NULL);
            }
            closedir(dir);    
        } 
        else
            WARN ("Skipping directory %s. Use --subdirs to compress directories.", fn);


        // remove the directory from the list of files to compress
        buf_remove (*fn_buf, char *, i, 1);
        i--;
    }

    ASSINP0 (fn_buf->len, "No work for me :( all files listed are directories.");

    flags_update (fn_buf->len, B1ST (rom , *fn_buf));

    *num_z_files = flag.pair ? (fn_buf->len / 2) : fn_buf->len; // note: flag.pair is only available in ZIP

    // sort files by data type to improve VB re-using, and refhash-using files in the end to improve reference re-using
    if (!flag.deep) {
        SET_FLAG (quiet); // suppress warnings
        qsort (B1ST (char *, *fn_buf), fn_buf->len, sizeof (char *), main_sort_input_filenames);
        RESTORE_FLAG (quiet);
    }
}

static void main_load_reference (rom filename, bool is_first_file, bool is_last_z_file)
{
    if (!IS_REF_EXTERNAL && !IS_REF_EXT_STORE) return;

    int old_aligner_available = flag.aligner_available;
    DataType dt = main_get_file_dt (filename);
    flag.aligner_available = primary_command == ZIP && 
                            (old_aligner_available || dt == DT_FASTQ || dt == DT_FASTA ||
                             ((dt==DT_SAM || dt==DT_BAM) && flag.best)); // load refhash only in --best

    // no need to load the reference if not needed (unless its genocat of the refernece file itself)
    if (flag.genocat_no_ref_file && dt != DT_REF) return;

    // no need to load the reference if just collecting coverage except FASTQ for which we need the contigs
    if (flag.collect_coverage && dt != DT_FASTQ) return;

    RESET_VALUE (txt_file); // save and reset - for use by reference loader

    if (!ref_is_external_loaded (gref)) {
        MAIN ("Loading external reference to gref: %s", ref_get_filename(gref));
        ref_load_external_reference (gref, NULL); // also loads refhash if needed
    }

    if (dt == DT_CHAIN && !ref_is_external_loaded (prim_ref)) {

        ref_set_reference (prim_ref, NULL, REF_EXTERNAL, false); // set the prim_ref from $GENOZIP_REFERENCE if applicable

        ASSINP0 (ref_get_filename (prim_ref), 
                "When compressing a chain file, you must also specify two --reference arguments: the first is the reference file in Primary coordinates "
                "(i.e. those of the VCF files to be lifted), and the second is the reference file in Luft coordinates (i.e. the coordinates aftering lifting). See "WEBSITE_DVCF);

        MAIN ("Loading external reference to prim_ref: %s", ref_get_filename (prim_ref));
        ref_load_external_reference (prim_ref, NULL);
    }

    // Read the refhash and calculate the reverse compliment genome for the aligner algorithm - it was not used before and now it is
    if (!old_aligner_available && flag.aligner_available) {
        MAIN0 ("Loading refhash");
        refhash_load_standalone();
    }

    RESTORE_VALUE (txt_file);
}

static void set_exe_type (rom argv0)
{
    rom bn = filename_base (argv0, false, NULL, NULL, 0);
    
    if      (strstr (bn, "genols"))    exe_type = EXE_GENOLS;
    else if (strstr (bn, "genocat"))   exe_type = EXE_GENOCAT;
    else if (strstr (bn, "genounzip")) exe_type = EXE_GENOUNZIP;
    else                               exe_type = EXE_GENOZIP; // default

    FREE (bn);
}

int main (int argc, char **argv)
{    
    MAIN0 ("Starting main");

    // When debugging with Visual Studio Code, its debugger wrapper adds 3 args: "2>CON", "1>CON", "<CON". We remove them.
    if (argc >= 4 && !strcmp (argv[argc-1], "<CON")) argc -= 3;

    global_cmd = filename_base (argv[0], true, "(executable)", NULL, 0); // global var

    set_exe_type(argv[0]);
    info_stream = stdout; // may be changed during intialization
    profiler_initialize();
    buf_initialize(); 
    arch_initialize (argv[0]);
    evb = vb_initialize_nonpool_vb(VB_ID_EVB, DT_NONE, "main_thread");
    threads_initialize(); // requires evb
    random_access_initialize();
    codec_initialize();

    flags_init_from_command_line (argc, argv);
    if (is_genozip && IS_PIZ) exe_type = EXE_GENOUNZIP; // treat "genozip -d" as genounzip

    // --make-reference might be called by genocat or genounzip from ref_fasta_to_ref - we treat it as genozip
    if (flag.make_reference) exe_type = EXE_GENOZIP;

    flags_store_command_line (argc, argv); // can only be called after --password is processed

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        if (is_genols) command = LIST; // genols can be run without arguments
        
        else if (is_genocat && flag.show_headers) command = SHOW_HEADERS;

        // genozip with no input filename, no output filename, and no input redirection 
        // note: in docker stdin is a pipe even if going to a terminal. so we show the help even if
        // coming from a pipe. the user must use "-" to redirect from stdin
        else if (command == -1 && optind == argc && !flag.out_filename && !flag.files_from &&
                 (isatty(0) || !strcmp (arch_get_distribution(), "Docker"))) {
            // case: --register
            if (flag.do_register) {
                license_register();
                exit (EXIT_OK);
            }

            // case: requesting to display the reference: genocat --reference <ref-file> and optionally --regions
            if (is_genocat && (IS_REF_EXTERNAL || IS_REF_EXT_STORE)) {
                command = PIZ;
                flags_update (0, NULL);

                if (flag.show_ref_diff) 
                    ref_diff_ref();                    
                else
                    ref_display_ref (gref);
            }

            // genozip with no parameters and not registered yet - register now
            else if (is_genozip && argc == 1 && isatty(0) && !license_is_registered())
                license_register();
                
            // otherwise: show help
            else
                main_print_help (false);

            return 0;
        }

        else if (is_genounzip) command = PIZ;
        else if (is_genocat) { command = PIZ; flag.to_stdout = !flag.out_filename ; }
        else command = ZIP; // default 
    }

    unsigned num_z_files=0;
    static Buffer input_files_buf = { .name = "input_files" };
    main_get_filename_list (argc - optind, &argv[optind], &num_z_files, &input_files_buf);
    ARRAY (rom , input_files, input_files_buf);

    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (flag.threads_str) 
        ASSINP (str_get_int_range32 (flag.threads_str, 0, 1, MAX_GLOBAL_MAX_THREADS, (int32_t *)&global_max_threads),
                "invalid argument of --threads: \"%s\". Expecting an integer between 1 and %u.", flag.threads_str, MAX_GLOBAL_MAX_THREADS);

    else global_max_threads = MIN_(MAX_GLOBAL_MAX_THREADS,
                                   ((flag.is_windows || flag.is_mac || flag.is_wsl) 
                                        ? ((float)arch_get_num_cores() * 0.75)    // under-subscribe on Windows / Mac to maintain UI interactivity
                                        : ((float)arch_get_num_cores() * 1.1 ))); // over-subscribe to keep all cores busy even when some threads are waiting on mutex or join

    // handle all commands except for ZIP, PIZ or LIST
    if (command == VERSION) { main_print_version();   return 0; }
    if (command == LICENSE) { license_display();      return 0; }
    if (command == HELP)    { main_print_help (true); return 0; }

    ASSINP (input_files_len || !isatty(0) || command != ZIP, "missing input file. Example: %s myfile.bam", global_cmd);
    ASSINP (input_files_len || !isatty(0) || command != PIZ, "missing input file. Example: %s myfile.bam.genozip", global_cmd);

    primary_command = command; 

#ifdef DISTUNZIPMSG
    if (is_genounzip && strlen (DISTUNZIPMSG))
        iprintf ("%s\n\n", DISTUNZIPMSG);
#endif

    // we test for a newer version if its a single file compression (if --test is used, we test after PIZ - check_for_newer is set)
    if (!flag.quiet && ((IS_ZIP && input_files_len == 1 && !flag.test) || (IS_PIZ && flag.check_latest/*PIZ - test after compress*/)))
        version_background_test_for_newer();

    // IF YOU'RE CONSIDERING EDITING THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
    // and might put you personally as well as your organization at legal and financial risk - see "Severly unauthorized use of Genozip"
    // section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
    if (IS_ZIP) license_load(); 

    if (flag.reading_chain) {
        vcf_tags_cmdline_rename_option(); 
        vcf_tags_cmdline_drop_option();

        MAIN0 ("Loading chain file");
        chain_load();

        if (is_genocat) exit(0);
    }

    if (flag.reading_kraken) {
        MAIN0 ("Loading kraken file");
        kraken_load();
    }

    // if we're genozipping with tar, initialize tar file
    if (tar_is_tar()) tar_initialize();

    unsigned n_files = MAX_(input_files_len, 1);
    for (unsigned file_i=0, z_file_i=0; file_i < n_files; file_i++) {

        // get file name
        rom next_input_file = input_files_len ? input_files[file_i] : NULL;  // NULL means stdin

        if (flag.show_filename) iprintf ("%s:\n", next_input_file ? next_input_file : "(standard input)");

        if (next_input_file && !strcmp (next_input_file, "-")) next_input_file = NULL; // "-" is stdin too

        // // (bug 339) On Windows, we redirect files from stdin - For binary files \n is converted to \r\n. For
        // // textual files, if a \r\n appears in the first characters read directly in file_open_txt_read, then we wierdly
        // // receive only \n, and the byte immediately following the last character, is never read - not by
        // // the fread here, and not by the subsequent read - losing it. Likely a bug in the mingw libc.
        // ASSINP0 (!flag.is_windows || next_input_file, "genozip on Windows does not support piping in data");

        ASSINP0 (next_input_file || command != PIZ, "filename(s) required (redirecting from stdin is not possible)");

        ASSERTW (next_input_file || !flag.replace, "%s: ignoring %s option", global_cmd, OT("replace", "^")); 

        bool is_last_txt_file = (file_i==input_files_len-1);
        bool is_last_z_file = (z_file_i==num_z_files-1);

        main_load_reference (next_input_file, !file_i, is_last_z_file);

        switch (command) {
            case ZIP  : main_genozip (next_input_file, 
                                      (next_input_file && file_i < input_files_len-1) ? input_files[file_i+1] : NULL, // file name of next file, if there is one
                                      file_i, n_files, !next_input_file || is_last_txt_file); 
                        break;

            case PIZ  : main_genounzip (next_input_file, flag.out_filename, file_i, n_files, is_last_z_file); 
                        break;           

            case SHOW_HEADERS : genocat_show_headers(next_input_file); break;
                        
            case LIST : genols (next_input_file, false, NULL, false); break;

            default   : ABORTINP ("unrecognized command %c", command);
        }

        if (flag.pair) flag.pair = 3 - flag.pair; // alternate between PAIR_READ_1 and PAIR_READ_2

        if (!z_file) z_file_i++; // z_file was closed, meaning we're done with it
    }

    // if this is "list", finalize
    if (command == LIST) genols (NULL, true, NULL, false);

    // if we're genozipping with tar, finalize tar file
    if (tar_is_tar()) tar_finalize();

    if (flag.multiple_files && flag.validate==VLD_REPORT_INVALID /* reporting invalid files, and none found */) 
        WARN0 ("All files are valid genozip files");

    ASSINP0 (!flag.biopsy, "biopsy VB not found, try a lower number");
    ASSINP0 (flag.biopsy_line.line_i == NO_LINE, "biopsy line not found, try a lower number");

    exit_ok();
    return 0;
}
