// ------------------------------------------------------------------
//   dispatcher.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <time.h>
#ifndef _MSC_VER // Microsoft compiler
#include <pthread.h>
#else
#include "compatibility/visual_c_pthread.h"
#include "compatibility/visual_c_gettime.h"
#endif
#if defined __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif

#include "genozip.h"
#include "dispatcher.h"
#include "vblock.h"
#include "file.h"
#include "profiler.h"

typedef struct {
    pthread_t thread_id;
    VBlock *vb;
    void (*func)(VBlock *);
} Thread;

typedef struct {
    unsigned max_vb_id_so_far; 
    Buffer compute_threads_buf;
    Thread *compute_threads;
    VBlock *next_vb; // next vb to be dispatched
    VBlock *processed_vb; // last vb for which caller got the processing results

    bool input_exhausted;

    unsigned next_thread_to_dispatched;
    unsigned next_thread_to_be_joined;

    unsigned num_running_compute_threads;
    unsigned next_vb_i;
    unsigned max_threads;
    bool test_mode;
    bool is_last_file;

    // progress indicator stuff
    TimeSpecType start_time; 
    bool show_progress;
    double last_percent;
    unsigned last_len;
    unsigned last_seconds_so_far;
    const char *filename;
} DispatcherData;

// variables that persist across multiple dispatchers run sequentially
static TimeSpecType profiler_timer; // wallclock
static bool ever_time_initialized = false;
static TimeSpecType ever_time;

void dispatcher_show_time (const char *stage, int32_t thread_index, uint32_t vb_i)
{
    static bool initialized = false;
    static const char *prev_stage;
    static int32_t prev_thread_index;
    static uint32_t prev_vb_i;
    static TimeSpecType prev_timer;

    TimeSpecType timer; 
    clock_gettime(CLOCK_REALTIME, &timer); 

    int diff_micro = 0;
    if (initialized) {
        diff_micro = 1000000 *(timer.tv_sec - prev_timer.tv_sec) + (int)((int64_t)timer.tv_nsec - (int64_t)prev_timer.tv_nsec) / 1000;
        fprintf (stderr, "TH=%-2d VB=%-3u Stage='%s' Microsec_in_this_stage=%u\n", prev_thread_index, prev_vb_i, prev_stage, diff_micro);
    }

    initialized = true;
    prev_stage        = stage;
    prev_timer        = timer;
    prev_thread_index = thread_index;
    prev_vb_i         = vb_i;
}

static void dispatcher_human_time (unsigned secs, char *str /* out */)
{
    unsigned hours = secs / 3600;
    unsigned mins  = (secs % 3600) / 60;
             secs  = secs % 60;

    if (hours) 
        sprintf (str, "%u %s %u %s", hours, hours==1 ? "hour" : "hours", mins, mins==1 ? "minute" : "minutes");
    else if (mins)
        sprintf (str, "%u %s %u %s", mins, mins==1 ? "minute" : "minutes", secs, secs==1 ? "second" : "seconds");
    else 
        sprintf (str, "%u %s", secs, secs==1 ? "second" : "seconds");
}

const char *dispatcher_ellapsed_time (Dispatcher dispatcher, bool ever)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    TimeSpecType start = ever ? ever_time : dd->start_time;

    int seconds_so_far = ((tb.tv_sec - start.tv_sec)*1000 + (tb.tv_nsec - start.tv_nsec) / 1000000) / 1000; 

    static char time_str[70];
    dispatcher_human_time (seconds_so_far, time_str);

    return time_str;
}

static void dispatcher_show_start (DispatcherData *dd)
{
    if (!dd->show_progress) return; 

    const char *progress = (command == ZIP && txt_file->redirected) ? "Compressing...\b\b\b\b\b\b\b\b\b\b\b\b\b\b" : "0\%"; // we can't show % when compressing from stdin as we don't know the file size
    
    if (dd->test_mode) 
        fprintf (stderr, "testing: %s%s --test %s : %s", global_cmd, strstr (global_cmd, "genozip") ? " --decompress" : "", 
                 dd->filename, progress); 
    else
        fprintf (stderr, "%s %s : %s", global_cmd, dd->filename, progress); 
             
    dd->last_len = strlen (progress); // so dispatcher_show_progress knows how many characters to erase
}

static void dispatcher_show_progress (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    
    if (!dd->show_progress && !flag_debug_progress) return; 

    TimeSpecType tb; 
    clock_gettime(CLOCK_REALTIME, &tb); 
    
    int seconds_so_far = ((tb.tv_sec-dd->start_time.tv_sec)*1000 + (tb.tv_nsec-dd->start_time.tv_nsec) / 1000000) / 1000; 

    uint64_t total=0, sofar=0;
    
    // case: genozip of plain txt files (including if decompressed by an external compressor) 
    // - we go by the amount of txt content processed 
    if (command == ZIP && txt_file->disk_size && file_is_plain_or_ext_decompressor (txt_file)) { 
        total = txt_file->txt_data_size_single; // if its a physical plain VCF file - this is the file size. if not - its an estimate done after the first VB
        sofar = z_file->txt_data_so_far_single;
    } 
    
    // case: UNZIP: always ; ZIP: locally decompressed files - eg .vcf.gz or .fq.bz2 - 
    // we go by the physical disk size and how much has been consumed from it so far
    else if (command == UNZIP || 
             (command == ZIP && file_is_read_via_int_decompressor (txt_file))) {
        File *input_file  = (command == ZIP ? txt_file : z_file);
        total = input_file->disk_size; 
        sofar = input_file->disk_so_far; 
    } 
        
    // case: we have no idea what is the disk size of the VCF file - for example, because its STDIN, or because
    // its coming from a URL that doesn't provide the size
    else if (command == ZIP && !txt_file->disk_size)
        return; // we can't show anything if we don't know the file size

    else ABORT ("Error in dispatcher_show_progress: unsupported case: command=%u txt_file->type=%s", command, ft_name (txt_file->type));

    double percent;
    if (total > 10000000) // gentle handling of really big numbers to avoid integer overflow
        percent = MIN (((double)(sofar/100000ULL)*100) / (double)(total/100000ULL), 100.0); // divide by 100000 to avoid number overflows
    else
        percent = MIN (((double)sofar*100) / (double)total, 100.0); // divide by 100000 to avoid number overflows
    
    // need to update progress indicator, max once a second or if 100% is reached
    const char *eraser = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";

    // in split mode - dispatcher is not done if there's another component after this one
    bool done = (dispatcher_is_done (dispatcher)) && (!flag_split || !buf_is_allocated (&z_file->v1_next_vcf_header));

    // case: we've reached 99% prematurely... we under-estimated the time
    if (!done && percent > 99 && (dd->last_seconds_so_far < seconds_so_far)) {
        const char *progress = "Finalizing...";

        // note we have spaces at the end to make sure we erase the previous string, if it is longer than the current one
        fprintf (stderr, "%.*s%s            %.12s", dd->last_len, eraser, progress, eraser);

        dd->last_len = strlen (progress);
    }

    // case: we're making progress... show % and time remaining
    else if (!done && percent && (dd->last_seconds_so_far < seconds_so_far)) { 

        if (!dispatcher_is_done (dispatcher) ||
            (flag_split && buf_is_allocated (&z_file->v1_next_vcf_header))) { 

            // time remaining
            char time_str[70], progress[100];
            unsigned secs = (100.0 - percent) * ((double)seconds_so_far / (double)percent);
            dispatcher_human_time (secs, time_str);
            sprintf (progress, "%u%% (%s)", (unsigned)percent, time_str);

            // note we have spaces at the end to make sure we erase the previous string, if it is longer than the current one
            if (!flag_debug_progress)
                fprintf (stderr, "%.*s%s            %.12s", dd->last_len, eraser, progress, eraser);
            else
                fprintf (stderr, "%u%% (%s) sofar=%"PRIu64" total=%"PRIu64" seconds_so_far=%d\n", (unsigned)percent, time_str, sofar, total, seconds_so_far);

            dd->last_len = strlen (progress);
        }
    }

    // case: we're done - caller will print the "Done" message after finalizing the genozip header etc
    else if (done) fprintf (stderr, "%.*s", dd->last_len, eraser); 

    dd->last_percent = percent;
    dd->last_seconds_so_far = seconds_so_far;
}

Dispatcher dispatcher_init (unsigned max_threads, unsigned previous_vb_i,
                            bool test_mode, bool is_last_file, const char *filename)
{
    clock_gettime(CLOCK_REALTIME, &profiler_timer);

    DispatcherData *dd = (DispatcherData *)calloc (1, sizeof(DispatcherData));
    ASSERT0 (dd, "failed to calloc DispatcherData");

    clock_gettime(CLOCK_REALTIME, &dd->start_time); 

    if (!ever_time_initialized) {
        ever_time = dd->start_time;
        ever_time_initialized = true;
    }

    dd->next_vb_i     = previous_vb_i;  // used if we're concatenating files - the vblock_i will continue from one file to the next
    dd->max_threads   = max_threads;
    dd->test_mode     = test_mode;
    dd->is_last_file  = is_last_file;
    dd->show_progress = !flag_quiet && !!isatty(2);
    dd->filename      = filename;

    vb_create_pool (MAX (2,max_threads+1 /* one for evb */));

    buf_alloc (evb, &dd->compute_threads_buf, sizeof(Thread) * MAX (1, max_threads), 1, "compute_threads_buf", 0);
    dd->compute_threads = (Thread *)dd->compute_threads_buf.data;

    if (!flag_split) dispatcher_show_start (dd); // note: for flag_split, we print this in dispatcher_resume() 

    return dd;
}

void dispatcher_pause (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->next_vb_i--;
}

// reinit dispatcher, used when splitting a genozip to its vcf components, using a single dispatcher object
void dispatcher_resume (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    clock_gettime(CLOCK_REALTIME, &dd->start_time); 

    dd->input_exhausted = false;
    dd->last_len        = 2;
    dd->filename        = txt_file->name;
    
    dispatcher_show_start (dd);    
}

void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i)
{
    DispatcherData *dd = (DispatcherData *)*dispatcher;

    COPY_TIMER (evb->profile.wallclock);

    if (flag_show_time) 
        profiler_print_report (&evb->profile, 
                               dd->max_threads, dd->max_vb_id_so_far+1,
                               dd->filename, dd->next_vb_i + (command != ZIP)); // in ZIP, the last VB is empty

    // must be before vb_cleanup_memory() 
    if (flag_show_memory) buf_display_memory_usage (false, dd->max_threads, dd->max_vb_id_so_far);    

    buf_destroy (&dd->compute_threads_buf); // we need to destroy (not marely free) because we are about to free dd

    // note: we can only test evb when no compute thread is running as compute threads might modify evb buffers
    // mid-way through test causing a buffer to have an inconsiset state and for buf_test_overflows to therefore report an error
    buf_test_overflows(evb, "dispatcher_finish"); 

    // free memory allocations that assume subsequent files will have the same number of samples.
    // (we assume this if the files are being concatenated). don't bother freeing (=same time) if this is the last file
    if (!flag_concat && !dd->is_last_file) {
        vb_cleanup_memory(); 
        vb_release_vb (evb);
    }
    
    if (last_vb_i) *last_vb_i = dd->next_vb_i; // for continuing vblock_i count between subsequent concatented files

    FREE (dd);

    *dispatcher = NULL;
}

static void *dispatcher_thread_entry (void *thread_)
{
    Thread *th = (Thread *)thread_;

    th->func (th->vb);
    
    return NULL;
}

VBlock *dispatcher_generate_next_vb (Dispatcher dispatcher, uint32_t vb_i)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->next_vb_i = vb_i ? vb_i : dd->next_vb_i+1;

    if (flag_show_threads) dispatcher_show_time ("Generate vb", -1, dd->next_vb_i);

    dd->next_vb = vb_get_vb (dd->next_vb_i);
    dd->max_vb_id_so_far = MAX (dd->max_vb_id_so_far, dd->next_vb->id);

    return dd->next_vb;
}

void dispatcher_compute (Dispatcher dispatcher, void (*func)(VBlockP))
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    Thread *th = &dd->compute_threads[dd->next_thread_to_dispatched];

    th->vb = dd->next_vb;
    th->func = func;

    if (flag_show_threads) dispatcher_show_time ("Start compute", dd->next_thread_to_dispatched, th->vb->vblock_i);

    if (dd->max_threads > 1) {
        unsigned err = pthread_create(&th->thread_id, NULL, dispatcher_thread_entry, th);
        ASSERT (!err, "Error: failed to create thread for next_vb_i=%u, err=%u", dd->next_vb->vblock_i, err);

        dd->next_thread_to_dispatched = (dd->next_thread_to_dispatched + 1) % dd->max_threads;
    }
    else func(dd->next_vb); // single thread
                    
    dd->next_vb = NULL;
    dd->num_running_compute_threads++;
}

bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final) 
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (!dd->num_running_compute_threads) return false; // no running compute threads 

    Thread *th = &dd->compute_threads[dd->next_thread_to_be_joined];

    bool my_is_final = dd->input_exhausted && !dd->next_vb && dd->num_running_compute_threads == 1; // this is the last vb to be processed

    if (is_final) *is_final = my_is_final;

    return my_is_final || (th->vb && th->vb->is_processed);
}

VBlock *dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (dd->max_threads > 1 && !dd->num_running_compute_threads) return NULL; // no running compute threads 

    Thread *th = &dd->compute_threads[dd->next_thread_to_be_joined];

    if (flag_show_threads) dispatcher_show_time ("Wait for thread", dd->next_thread_to_be_joined, th->vb->vblock_i);

    if (dd->max_threads > 1) 
        // wait for thread to complete (possibly it completed already)
        pthread_join(th->thread_id, NULL);

    if (flag_show_threads) dispatcher_show_time ("Join (end compute)", dd->next_thread_to_be_joined, th->vb->vblock_i);

    dd->processed_vb = th->vb;
    
    memset (th, 0, sizeof(Thread));
    dd->num_running_compute_threads--;
    dd->next_thread_to_be_joined = (dd->next_thread_to_be_joined + 1) % MAX (1, dd->max_threads);

    return dd->processed_vb;
}

bool dispatcher_has_free_thread (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    return dd->num_running_compute_threads < MAX(1, dd->max_threads);
}

VBlock *dispatcher_get_next_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    return dd->next_vb;
}

void dispatcher_abandon_next_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (!dd->next_vb) return;

    buf_test_overflows(dd->next_vb, "dispatcher_abandon_next_vb"); 

    if (flag_show_time) profiler_add (&evb->profile, &dd->next_vb->profile);

    vb_release_vb (dd->next_vb); 
    dd->next_vb = NULL;
}

void dispatcher_finalize_one_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (dd->processed_vb) {

        buf_test_overflows(dd->processed_vb, "dispatcher_finalize_one_vb"); // just to be safe, this isn't very expensive

        if (flag_show_time) profiler_add (&evb->profile, &dd->processed_vb->profile);

        vb_release_vb (dd->processed_vb); // cleanup vb and get it ready for another usage (without freeing memory)
        dd->processed_vb = NULL;
    }

    dispatcher_show_progress (dispatcher);
}                           

void dispatcher_input_exhausted (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    dd->input_exhausted = true;

    vb_release_vb (dd->next_vb);
    dd->next_vb = NULL;

    dd->next_vb_i--; // we didn't use this vb_i
}    

bool dispatcher_is_done (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    return dd->input_exhausted && !dd->next_vb && !dd->processed_vb && !dd->num_running_compute_threads;
}

bool dispatcher_is_input_exhausted (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    return dd->input_exhausted;
}

void dispatcher_fan_out_task (const char *task_name, bool test_mode, 
                              DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output)
{
    Dispatcher dispatcher = dispatcher_init (global_max_threads, 0, test_mode, true, task_name);
    do {
        VBlock *next_vb = dispatcher_get_next_vb (dispatcher);
        bool has_vb_ready_to_compute = next_vb && next_vb->ready_to_dispatch;
        bool has_free_thread = dispatcher_has_free_thread (dispatcher);

        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
        if (has_vb_ready_to_compute && has_free_thread)
            dispatcher_compute (dispatcher, compute);
        
        // PRIORITY 2: output completed vbs, so they can be released and re-used
        else if (dispatcher_has_processed_vb (dispatcher, NULL) ||  // case 1: there is a VB who's compute processing is completed
                 (has_vb_ready_to_compute && !has_free_thread)) {   // case 2: a VB ready to dispatch but all compute threads are occupied. wait here for one to complete
           
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); // this will block until one is available
            if (!processed_vb) continue; // no running compute threads 

            if (output) output (processed_vb);
            
            dispatcher_finalize_one_vb (dispatcher);
        }        
        
        // PRIORITY 3: If there is no variant block available to compute or to output, but input is not exhausted yet - get one range
        else if (!next_vb && !dispatcher_is_input_exhausted (dispatcher)) {

            next_vb = dispatcher_generate_next_vb (dispatcher, 0);

            prepare (next_vb);

            if (!next_vb->ready_to_dispatch) {
                dispatcher_input_exhausted (dispatcher);
                dispatcher_finalize_one_vb (dispatcher); 
            }
        }

        // if no condition was met, we're either done, or we still have some threads processing that are not done yet.
        // we wait a bit to avoid a tight busy loop
        else usleep (5000); 

    } while (!dispatcher_is_done (dispatcher));
}
