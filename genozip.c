// ------------------------------------------------------------------
//   genozip.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <getopt.h>

#include "genozip.h"
#include "text_help.h"
#include "version.h" // automatically incremented by the make when we create a new distribution
#include "zfile.h"
#include "zip.h"
#include "piz.h"
#include "crypt.h"
#include "file.h"
#include "vblock.h"
#include "url.h"
#include "strings.h"
#include "stats.h"
#include "arch.h"
#include "license.h"
#include "reference.h"
#include "refhash.h"
#include "random_access.h"
#include "codec.h"
#include "threads.h"
#include "linesorter.h"
#include "bases_filter.h"
#include "genols.h"

// globals - set it main() and never change
const char *global_cmd = NULL; 
ExeType exe_type;

// primary_command vs command: primary_command is what the user typed on the command line. command is what is 
// running now - for example, when ZIP is unzipping a reference, primary_command=ZIP and command=PIZ
CommandType command = NO_COMMAND, primary_command = NO_COMMAND; 

uint32_t global_max_threads = DEFAULT_MAX_THREADS; 

void main_exit (bool show_stack, bool is_error) 
{
    if (is_error && flag.debug_threads)
        threads_write_log (true);
        
    if (show_stack) threads_print_call_stack(); // this works ok on mac, but seems to not print function names on Linux

    buf_test_overflows_all_vbs("exit_on_error");

    if (flag.log_filename) {
        iprintf ("%s - execution ended %s\n", str_time().s, is_error ? "with an error" : "normally");
        fclose (info_stream);
    } 

    if (is_error) {
        close (1);   // prevent other threads from outputting to terminal (including buffered output), obscuring our error message
        close (2);
        url_kill_curl();  /* <--- BREAKPOINT */
        file_kill_external_compressors(); 
    }

    // if we're in ZIP - remove failed genozip file (but don't remove partial failed text file in PIZ - it might be still useful to the user)
    if (primary_command == ZIP && z_file && z_file->name && !flag_loading_auxiliary) {
        char *save_name = NULL;
        if (z_file && z_file->name) {
            save_name = malloc (strlen (z_file->name)+1);
            strcpy (save_name, z_file->name);
        }

        // cancel all other threads before closing z_file, so other threads don't attempt to access it 
        // (eg. z_file->data_type) and get a segmentation fault.
        threads_cancel_other_threads();

        file_close (&z_file, false, false); // also frees file->name

        // note: logic to avoid a race condition causing the file not to be removed - if another thread seg-faults
        // because it can't access a z_file filed after z_file is freed, and threads_sigsegv_handler aborts
        if (save_name) file_remove (save_name, true);
    }

    exit (is_error ? EXIT_GENERAL_ERROR : EXIT_OK);
} 

static void main_print_help (bool explicit)
{
    static const char **texts[NUM_EXE_TYPES] = {help_genozip, help_genounzip, help_genols, help_genocat }; // same order as ExeType
    static unsigned sizes[] = {sizeof(help_genozip), sizeof(help_genounzip), sizeof(help_genols), sizeof(help_genocat), sizeof(help_genozip_developer)};
    
    if (flag.help && !strcmp (flag.help, "dev")) 
        str_print_text (help_genozip_developer, sizeof(help_genozip_developer) / sizeof(char*), 
                        "                          ",  "\n", 0);
    else if (flag.help && !strcmp (flag.help, "genozip")) 
        str_print_text (help_genozip, sizeof(help_genozip) / sizeof(char*), 
                        "                          ",  "\n", 0);

    else if (flag.help && !strcmp (flag.help, "genounzip")) 
        str_print_text (help_genounzip, sizeof(help_genounzip) / sizeof(char*), 
                        "                          ",  "\n", 0);

    else if (flag.help && !strcmp (flag.help, "genols")) 
        str_print_text (help_genols, sizeof(help_genols) / sizeof(char*), 
                        "                          ",  "\n", 0);

    else if (flag.help && !strcmp (flag.help, "genocat")) 
        str_print_text (help_genocat, sizeof(help_genocat) / sizeof(char*), 
                        "                          ",  "\n", 0);

    else if (flag.help && !strcmp (flag.help, "input")) 
        iprintf ("Supported file types for --input:\n%s\n", file_compressible_extensions (false));
    
    else 
        str_print_text (texts[exe_type], sizes[exe_type] / sizeof(char*), 
                        "                     ",  "\n", 0);
    
    str_print_text (help_footer, sizeof(help_footer) / sizeof(char*), "", "\n", 0);

    // in Windows, we ask the user to click a key - this is so that if the user double clicks on the EXE
    // from Windows Explorer - the terminal will open and he will see the help
    if (flag.is_windows && !explicit && isatty(0) && isatty (1)) {
        printf ("Press any key to continue...\n");
        getc(stdin);
    }
}

static void main_print_version()
{
    iprintf ("version=%s\n", GENOZIP_CODE_VERSION);  
}

static void main_genounzip (const char *z_filename, const char *txt_filename, int z_file_i, bool is_last_z_file)
{
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

    if (z_dual_coords)
        z_file->max_conc_writing_vbs = linesorter_get_max_conc_writing_vbs(); // get data by reading all SEC_RECON_PLAN headers

    if (z_file->max_conc_writing_vbs < 3) z_file->max_conc_writing_vbs = 3;

    // case: reference not loaded yet bc --reference wasn't specified, and we got the ref name from zfile_read_genozip_header()   
    if (flag.reference == REF_EXTERNAL && !ref_is_external_loaded(gref)) {
        ASSINP0 (flag.reference != REF_EXTERNAL || !flag.show_ref_seq, "--show-ref-seq cannot be used on a file that requires a reference file: use genocat --show-ref-seq on the reference file itself instead");

        if (!flag.genocat_no_reconstruct || 
            (flag.collect_coverage && z_file->data_type == DT_FASTQ)) { // in collect_coverage with FASTQ we read the non-data sections of the reference
            SAVE_VALUE (z_file); // actually, read the reference first
            ref_load_external_reference (gref, NULL);
            RESTORE_VALUE (z_file);
        }
    }

    flags_update_piz_one_file (z_file_i);
    
    // set txt_filename from genozip file name (inc. extensions if translating or --bgzf)
    if (!txt_filename && !flag.to_stdout && !flag.unbind) 
        txt_filename = txtfile_piz_get_filename (z_filename, "", true);

    // open output txt file (except if unbinding)
    if (txt_filename || (flag.to_stdout && !flag.no_writer)) {
        ASSERT0 (!txt_file, "txt_file is unexpectedly already open"); // note: in bound mode, we expect it to be open for 2nd+ file
        txt_file = file_open (txt_filename, WRITE, TXT_FILE, flag.out_dt);
    }
    else if (flag.unbind || flag.no_writer) {
        // do nothing - in unbind, the component files will be opened by txtheader_piz_read_and_reconstruct()
    }
    else {
        ABORT0 ("Error: unrecognized configuration for the txt_file");
    }

    Dispatcher dispatcher = piz_z_file_initialize (is_last_z_file);  

    // a loop for decompressing all txt_files (1 if concatenating, possibly more if unbinding)
    if (dispatcher)
        while (z_file->num_txt_components_so_far < z_file->num_components) {
            piz_one_txt_file (dispatcher, is_first_z_file);
            file_close (&txt_file, flag.index_txt, flag.unbind || !is_last_z_file); 
        }

    file_close (&z_file, false, false);
    is_first_z_file = false;

    if (flag.replace && txt_filename && z_filename) file_remove (z_filename, true); 

done:
    RESTORE_FLAGS;
}

// run the test genounzip after genozip - for the most reliable testing that is nearly-perfectly indicative of actually 
// genounzipping, we create a new genounzip process
static void main_test_after_genozip (const char *exec_name, const char *z_filename, bool is_last_txt_file, bool is_chain)
{
    const char *password = crypt_get_password();

    // finish dumping reference and/or refhash to cache before destroying it...
    ref_create_cache_join (gref);
    ref_create_cache_join (prim_ref);

    refhash_create_cache_join();
    
    // On Windows and Mac that usually have limited memory, if ZIP consumed more than 2GB, free memory before PIZ. 
    // Note: on Windows, freeing memory takes considerable time.
#if defined _WIN32 || defined __APPLE__
    if (buf_get_memory_usage () > (1ULL<<31)) {
        ref_destroy_reference (gref, true); // on Windows I observed a race condition: if we unmap mapped memory here, and remap it in the test process, and the system is very slow due to low memory, then "MapViewOfFile" in the test processs will get "Access is Denied". That's why destroy_only_if_not_mmap=true.
        ref_destroy_reference (prim_ref, true); 
        kraken_destroy();
        chain_destroy();
        vb_destroy_pool_vbs();
    }
#endif

    StreamP test = stream_create (0, 0, 0, 0, 0, 0, 0,
                                  "To use the --test option",
                                  exec_name, "--decompress", "--test", z_filename,
                                  flag.quiet         ? "--quiet"         : SKIP_ARG,
                                  password           ? "--password"      : SKIP_ARG,
                                  password           ? password          : SKIP_ARG,
                                  flag.show_digest   ? "--show-digest"   : SKIP_ARG,
                                  flag.show_memory   ? "--show-memory"   : SKIP_ARG,
                                  flag.show_time     ? "--show-time"     : SKIP_ARG,
                                  flag.threads_str   ? "--threads"       : SKIP_ARG,
                                  flag.threads_str   ? flag.threads_str  : SKIP_ARG,
                                  flag.xthreads      ? "--xthreads"      : SKIP_ARG,
                                  flag.show_alleles  ? "--show-alleles"  : SKIP_ARG,
                                  flag.debug_threads ? "--debug-threads" : SKIP_ARG,
                                  flag.echo          ? "--echo"          : SKIP_ARG,
                                  flag.reference == REF_EXTERNAL && !is_chain ? "--reference" : SKIP_ARG, // normal pizzing of a chain file doesn't require a reference
                                  flag.reference == REF_EXTERNAL && !is_chain ? ref_get_filename(gref) : SKIP_ARG, 
                                  NULL);

    // wait for child process to finish, so that the shell doesn't print its prompt until the test is done
    int exit_code = stream_wait_for_exit (test);

    TEMP_VALUE (primary_command, TEST_AFTER_ZIP); // make exit_on_error NOT delete the genozip file in this case, so its available for debugging
    ASSINP (!exit_code, "test exited with status %d\n", exit_code);
    RESTORE_VALUE (primary_command); // recover in case of more non-concatenated files
}

static void main_genozip_open_z_file_write (const char **z_filename)
{
    DataType z_data_type = txt_file->data_type;

    if (! *z_filename) {

        ASSINP0 (txt_file->name, "use --output to specify the output filename (with a .genozip extension)");

        bool is_url = url_is_url (txt_file->name);
        const char *basename = is_url ? file_basename (txt_file->name, false, "", 0,0) : NULL;
        const char *local_txt_filename = basename ? basename : txt_file->name; 

        unsigned fn_len = strlen (local_txt_filename);
        *z_filename = (char *)MALLOC (fn_len + 30); // add enough the genozip extension e.g. 23andme.genozip

        // if the file has an extension matching its type, replace it with the genozip extension, if not, just add the genozip extension
        const char *genozip_ext = file_exts[file_get_z_ft_by_txt_in_ft (txt_file->data_type, txt_file->type)];

        if (chain_is_loaded && txt_file->data_type == DT_VCF)
            genozip_ext = DVCF_GENOZIP_;

        if (file_has_ext (local_txt_filename, file_exts[txt_file->type]))
            sprintf ((char *)*z_filename, "%.*s%s", (int)(fn_len - strlen (file_exts[txt_file->type])), local_txt_filename, genozip_ext); 
        else 
            sprintf ((char *)*z_filename, "%s%s", local_txt_filename, genozip_ext); 

        FREE (basename);
    }

    z_file = file_open (*z_filename, flag.pair ? WRITEREAD : WRITE, Z_FILE, z_data_type);
}

static void main_genozip (const char *txt_filename, 
                          const char *next_txt_filename, // ignored unless we are of pair_1 in a --pair
                          const char *z_filename,
                          unsigned txt_file_i, bool is_last_txt_file,
                          Coords txt_file_coords,
                          char *exec_name)
{
    SAVE_FLAGS;

    ASSINP (!z_filename || !url_is_url (z_filename), 
            "output files must be regular files, they cannot be a URL: %s", z_filename);

    ASSINP0 (!z_file || !z_dual_coords || flag.rejects_coord, 
             "it is not possible to concatenate dual-coordinate files");

    // get input file
    if (!txt_file) { // open the file - possibly already open from main_load_reference
        txt_file = file_open (txt_filename, READ, TXT_FILE, DT_NONE); 
        if (txt_file) txt_file->coords = txt_file_coords; // maybe overridden by ##dual_coordinates in the header.
    }

    // skip this file if its size is 0
    RETURNW (txt_file,, "Cannot compress file %s because its size is 0 - skipping it", txt_filename);

    if (!txt_file->file) goto done; // this is the case where multiple files are given in the command line, but this one is not compressible - we skip it

    ASSERT0 (flag.bind || !z_file, "expecting z_file to be NULL in non-bound mode");

    // get output FILE
    if (!z_file) { // skip if we're the second file onwards in bind mode, or pair_2 in unbound list of pairs - nothing to do
        
        if (flag.pair && !z_filename)  // case: if we're the first file in a pair
            z_filename = file_get_fastq_pair_filename (txt_filename, next_txt_filename, false);
        
        main_genozip_open_z_file_write (&z_filename);
    }

    if (!flag.rejects_coord)
        stats_add_txt_name (txt_name); // add txt file name to stats data stored in z_file

    TEMP_FLAG (quiet, flag.rejects_coord ? true : flag.quiet); // no warnings when (re) processing rejects

    flags_update_zip_one_file();

    bool z_closes_after_me = is_last_txt_file || flag.bind==BIND_NONE || (flag.bind==BIND_PAIRS && txt_file_i%2);

    zip_one_file (txt_file->basename, &is_last_txt_file, z_closes_after_me);

    if (flag.show_stats && z_closes_after_me && 
        (!z_dual_coords || flag.rejects_coord)) {
        stats_display();
    }

    bool remove_txt_file = z_file && flag.replace && txt_filename;

    file_close (&txt_file, false, !is_last_txt_file);  // no need to waste time closing the last file, the process termination will do that

    bool is_chain = (z_file->data_type == DT_CHAIN);

    // close the file if its an open disk file AND we need to close it
    if (!flag.to_stdout && z_file && z_closes_after_me) {

        // if this is a dual coords file, recursively call to add the rejects file, and close z_file there.
        if (z_dual_coords && !flag.rejects_coord) {
            flag.sort = false;
            z_closes_after_me = false;                
            
            for (flag.rejects_coord = DC_PRIMARY ; flag.rejects_coord <= DC_LUFT ; flag.rejects_coord++)
                main_genozip (z_file->rejects_file_name[flag.rejects_coord-1], 0, z_file->name, txt_file_i, flag.rejects_coord==DC_LUFT, 
                              flag.rejects_coord, exec_name);
        }
        
        else if (!z_dual_coords || flag.rejects_coord == DC_LUFT) {
            if (!z_filename) { z_filename = z_file->name ; z_file->name = 0; } // take over the name if we don't have it (eg 2nd file in a pair)
            file_close (&z_file, false, !is_last_txt_file); 
        }
    }

    RESTORE_FLAG (quiet); // don't pass on quiet to test, just because we turned it on for rejects 

    // test the compression, if the user requested --test
    if (flag.test && z_closes_after_me) main_test_after_genozip (exec_name, z_filename, is_last_txt_file, is_chain);

    // remove after test (if any) was successful - NOT GOOD ENOUGH - bug 342
    if (remove_txt_file) file_remove (txt_filename, true); 

done: {
    // propogate up some flags
    SAVE_FLAG(vblock_memory);  
    SAVE_FLAG(data_modified);

    RESTORE_FLAGS;
    
    if (flag.pair) RESTORE_FLAG (vblock_memory); // FASTQ R2 must have the same vblock_memory as R1, so it should not re-calculate in zip_dynamically_set_max_memory - otherwise txtfile_read_vblock will fail
    if (flag.bind) RESTORE_FLAG (data_modified); // When binding, if any of the compressed files is a Luft, we don't store digests
    
}}

static inline DataType main_get_file_dt (const char *filename)
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
    
    bool use_refhash1 = (dt1 == DT_FASTQ);
    bool use_refhash2 = (dt2 == DT_FASTQ);

    // refhash users - at the end of the list
    if (use_refhash1 != use_refhash2) return (int)use_refhash1 - (int)use_refhash2;

    // within refhash users, and within refhash non-users, sort by data type
    if (dt1 != dt2) return (int)dt1 - (int)dt2;

    // within files of the same data type, keep original order
    return fn1 - fn2;
}

static void main_load_reference (const char *filename, bool is_first_file, bool is_last_z_file)
{
    if (flag.reference != REF_EXTERNAL && flag.reference != REF_EXT_STORE) return;

    int old_ref_use_aligner = flag.ref_use_aligner;
    DataType dt = main_get_file_dt (filename);
    flag.ref_use_aligner = (old_ref_use_aligner || dt == DT_FASTQ || dt == DT_FASTA) && primary_command == ZIP;

    // no need to load the reference if not needed (unless its genocat of the refernece file itself)
    if (flag.genocat_no_ref_file && dt != DT_REF) return;

    // no need to load the reference if just collecting coverage except FASTQ for which we need the contigs
    if (flag.collect_coverage && dt != DT_FASTQ) return;

    // we also need the aligner if this is an unaligned SAM 
    if (!flag.ref_use_aligner && dt==DT_SAM && primary_command==ZIP) {

        // open here instead of in main_genozip
        txt_file = file_open (filename, READ, TXT_FILE, 0);

        // use the aligner if over 5 of the 100 first lines of the file are unaligned
        flag.ref_use_aligner = txt_file && txt_file->file && txtfile_test_data ('@', 100, 0.05, sam_zip_is_unaligned_line); 
    }

    RESET_VALUE (txt_file); // save and reset - for use by reference loader

    if (!ref_is_external_loaded (gref))
        ref_load_external_reference (gref, NULL); // also loads refhash if needed

    if (dt == DT_CHAIN && !ref_is_external_loaded (prim_ref)) {

        ref_set_reference (prim_ref, NULL, REF_EXTERNAL, false); // set the prim_ref from $GENOZIP_REFERENCE if applicable

        ASSINP0 (ref_get_filename (prim_ref), 
                "When compressing a chain file, you must also specify two --reference arguments: the first is the reference file in Primary coordinates "
                "(i.e. those of the VCF files to be lifted), and the second is the reference file in Luft coordinates (i.e. the coordinates aftering lifting). See "WEBSITE_DVCF);

        ref_load_external_reference (prim_ref, NULL);
    }

    // Read the refhash and calculate the reverse compliment genome for the aligner algorithm - it was not used before and now it is
    if (!old_ref_use_aligner && flag.ref_use_aligner) 
        refhash_load_standalone();

    RESTORE_VALUE (txt_file);
}

void TEST(char *str) {
}


int main (int argc, char **argv)
{    
    //TEST ("3.2E-003");

    info_stream = stdout; // may be changed during intialization
    profiler_initialize();
    buf_initialize(); 
    arch_initialize (argv[0]);
    evb = vb_initialize_nonpool_vb(EVB);
    threads_initialize(); // requires evb
    random_access_initialize();
    codec_initialize();

    if      (strstr (argv[0], "genols"))    exe_type = EXE_GENOLS;
    else if (strstr (argv[0], "genocat"))   exe_type = EXE_GENOCAT;
    else if (strstr (argv[0], "genounzip")) exe_type = EXE_GENOUNZIP;
    else                                    exe_type = EXE_GENOZIP; // default
    
    global_cmd = file_basename (argv[0], true, "(executable)", NULL, 0); // global var

    flags_init_from_command_line (argc, argv);
    flags_store_command_line (argc, argv); // can only be called after --password is processed

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        if (exe_type == EXE_GENOLS) command = LIST; // genols can be run without arguments
        
        // genozip with no input filename, no output filename, and no input redirection 
        // note: in docker stdin is a pipe even if going to a terminal. so we show the help even if
        // coming from a pipe. the user must use "-" to redirect from stdin
        else if (command == -1 && optind == argc && !flag.out_filename && 
                 (isatty(0) || arch_am_i_in_docker())) {
            // case: --register
            if (flag.do_register) {
                license_get();
                exit (EXIT_OK);
            }

            // case: requesting to display the reference: genocat --reference <ref-file> and optionally --regions
            if (exe_type == EXE_GENOCAT && (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE)) 
                ref_display_ref (gref);

            // otherwise: show help
            else
                main_print_help (false);

            return 0;
        }

        else if (exe_type == EXE_GENOUNZIP) command = PIZ;
        else if (exe_type == EXE_GENOCAT) { command = PIZ; flag.to_stdout = !flag.out_filename ; }
        else command = ZIP; // default 
    }

    unsigned num_files = argc - optind;

    flags_update (num_files, (const char **)&argv[optind]);

    unsigned num_z_files = command != ZIP          ? num_files     :
                           flag.bind == BIND_ALL   ? 1             :  // flag.bind was set by flags_update
                           flag.bind == BIND_PAIRS ? num_files / 2 :
                                                     num_files;

    // sort files by data type to improve VB re-using, and refhash-using files in the end to improve reference re-using
    SET_FLAG (quiet); // suppress warnings
    qsort (&argv[optind], num_files, sizeof (argv[0]), main_sort_input_filenames);
    RESTORE_FLAG (quiet);
    
    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (flag.threads_str) 
        ASSINP (str_get_int_range32 (flag.threads_str, 0, 1, 10000, (int32_t *)&global_max_threads),
                "invalid argument of --threads: \"%s\". Expecting an integer between 1 and 10000.", flag.threads_str);

    else global_max_threads = 
#if defined _WIN32 || defined __APPLE__
        (double)arch_get_num_cores() * 0.75; // under-subscribe on Windows / Mac to maintain UI interactivity
#else
        (double)arch_get_num_cores() * 1.1;  // over-subscribe to keep all cores busy even when some threads are waiting on mutex or join
#endif

    // handle all commands except for ZIP, PIZ or LIST
    if (command == VERSION) { main_print_version();   return 0; }
    if (command == LICENSE) { license_display();      return 0; }
    if (command == HELP)    { main_print_help (true); return 0; }

    primary_command = command; 

    // IF YOU'RE CONSIDERING EDITING THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
    // and might put you personally as well as your organization at legal and financial risk - Genozip copyright holders
    // view themselves as entitled to any or all of the financial gains you and your organization make while using an 
    // unlicensed copy of Genozip. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
    if (command == ZIP) license_get(); 

    if (flag.reading_chain) {
        chain_load();
        if (exe_type == EXE_GENOCAT) exit(0);
    }

    if (flag.reading_kraken) kraken_load();

    for (unsigned file_i=0, z_file_i=0; file_i < MAX (num_files, 1); file_i++) {

        // get file name
        char *next_input_file = optind < argc ? argv[optind++] : NULL;  // NULL means stdin

        if (flag.show_filename) iprintf ("%s:\n", next_input_file ? next_input_file : "(standard input)");

        if (next_input_file && !strcmp (next_input_file, "-")) next_input_file = NULL; // "-" is stdin too

        ASSINP0 (next_input_file || command != PIZ, "filename(s) required (redirecting from stdin is not possible)");

        ASSERTW (next_input_file || !flag.replace, "%s: ignoring %s option", global_cmd, OT("replace", "^")); 

        bool is_last_txt_file = (file_i==num_files-1);
        bool is_last_z_file = (z_file_i==num_z_files-1);

        main_load_reference (next_input_file, !file_i, is_last_z_file);
        
        switch (command) {
            case ZIP  : main_genozip (next_input_file, 
                                      optind < argc ? argv[optind] : NULL, // file name of next file, if there is one
                                      flag.out_filename, file_i, !next_input_file || is_last_txt_file, DC_NONE, argv[0]); 
                        break;

            case PIZ  : main_genounzip (next_input_file, flag.out_filename, file_i, is_last_z_file); break;           

            case LIST : genols (next_input_file, false, NULL, false); break;

            default   : ABORTINP ("unrecognized command %c", command);
        }

        if (flag.pair) flag.pair = 3 - flag.pair; // alternate between PAIR_READ_1 and PAIR_READ_2

        if (!z_file) z_file_i++; // z_file was closed, meaning we're done with it
    }

    // if this is "list", finalize
    if (command == LIST) genols (NULL, true, NULL, false);

    if (flag.multiple_files && flag.validate==VLD_REPORT_INVALID /* reporting invalid files, and none found */) 
        WARN0 ("All files are valid genozip files");

    // finish dumping reference and/or refhash to cache
    ref_create_cache_join (gref);
    ref_create_cache_join (prim_ref);
    refhash_create_cache_join();

    exit_ok;
    return 0;
}
