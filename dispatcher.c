// ------------------------------------------------------------------
//   dispatcher.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt


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
#include "progress.h"

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
    bool is_last_file;
    ProgressType prog;
    const char *filename;
} DispatcherData;

// variables that persist across multiple dispatchers run sequentially
static TimeSpecType profiler_timer; // wallclock

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

static void dispatcher_show_progress (Dispatcher dispatcher)
{
    uint64_t total=0, sofar=0;
    
    // case: ZIP of plain txt files (including if decompressed by an external compressor) 
    // - we go by the amount of txt content processed 
    if (command == ZIP && txt_file->disk_size && file_is_plain_or_ext_decompressor (txt_file)) { 
        total = txt_file->txt_data_size_single; // if its a physical plain VCF file - this is the file size. if not - its an estimate done after the first VB
        sofar = z_file->txt_data_so_far_single;
    } 
    
    // case: PIZ: always ; ZIP: locally decompressed files - eg .vcf.gz or .fq.bz2 - 
    // we go by the physical disk size and how much has been consumed from it so far
    else if (command == PIZ || 
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

    // in unbind mode - dispatcher is not done if there's another component after this one
    bool done = dispatcher_is_done (dispatcher);

    progress_update (sofar, total, done);
}

Dispatcher dispatcher_init (unsigned max_threads, unsigned previous_vb_i,
                            bool test_mode, bool is_last_file, 
                            const char *filename, // filename, or NULL if filename is unchanged
                            ProgressType prog, const char *prog_msg /* used if prog=PROGRESS_MESSAGE */)   
{
    clock_gettime (CLOCK_REALTIME, &profiler_timer);

    DispatcherData *dd = (DispatcherData *)CALLOC (sizeof(DispatcherData));
    dd->next_vb_i    = previous_vb_i;  // used if we're binding files - the vblock_i will continue from one file to the next
    dd->max_threads  = max_threads;
    dd->is_last_file = is_last_file;
    dd->prog         = prog;

    if (filename)
        dd->filename = filename;

    ASSERT (max_threads <= global_max_threads, "Error in dispatcher_init: expecting max_threads=%u <= global_max_threads=%u", max_threads, global_max_threads);
    
    // always create the pool based on global_max_threads, not max_threads, because it is the same pool throughout the execution
    vb_create_pool (MAX (2,global_max_threads+1 /* one for evb */));

    buf_alloc (evb, &dd->compute_threads_buf, sizeof(Thread) * MAX (1, max_threads), 1, "compute_threads_buf");
    dd->compute_threads = (Thread *)dd->compute_threads_buf.data;

    if (!flag.unbind && filename) // note: for flag.unbind (in main file), we print this in dispatcher_resume() 
        progress_new_component (filename, 
                                command == ZIP && txt_file->redirected ? "Compressing..." : "0\%", // we can't show % when compressing from stdin as we don't know the file size
                                test_mode); 

    if (prog == PROGRESS_MESSAGE)
        progress_update_status (prog_msg);

    return dd;
}

void dispatcher_pause (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->next_vb_i--;
}

// PIZ: reinit dispatcher, used when splitting a genozip to its components, using a single dispatcher object
void dispatcher_resume (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->input_exhausted = false;
    dd->filename        = txt_file->name;
    
    progress_new_component (dd->filename, "0\%", -1);    
}

void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i)
{
    DispatcherData *dd = (DispatcherData *)*dispatcher;

    COPY_TIMER_VB (evb, wallclock);

    if (flag.show_time && !flag.show_time[0]) // show-time without the optional parameter 
        profiler_print_report (&evb->profile, 
                               dd->max_threads, dd->max_vb_id_so_far+1,
                               dd->filename, dd->next_vb_i + (command != ZIP)); // in ZIP, the last VB is empty

    // must be before vb_cleanup_memory() 
    if (flag.show_memory) buf_display_memory_usage (false, dd->max_threads, dd->max_vb_id_so_far);    

    buf_destroy (&dd->compute_threads_buf); // we need to destroy (not marely free) because we are about to free dd

    // note: we can only test evb when no compute thread is running as compute threads might modify evb buffers
    // mid-way through test causing a buffer to have an inconsiset state and for buf_test_overflows to therefore report an error
    buf_test_overflows(evb, "dispatcher_finish"); 

    // free memory allocations that assume subsequent files will have the same number of samples.
    // (we assume this if the files are being bound). don't bother freeing (=same time) if this is the last file, unless we're going to test and need the memory
    if (!flag.bind && (!dd->is_last_file || flag.test)) {
        vb_cleanup_memory(); 
        vb_release_vb (evb);
    }
    
    if (last_vb_i) *last_vb_i = dd->next_vb_i; // for continuing vblock_i count between subsequent bound files

    FREE (*dispatcher);
}

static void *dispatcher_thread_entry (void *thread_)
{
    Thread *th = (Thread *)thread_;

    ASSERT0 (th->vb->vblock_i, "Error in dispatcher_thread_entry: vb_i=0");

    th->func (th->vb);
    
    return NULL;
}

VBlock *dispatcher_generate_next_vb (Dispatcher dispatcher, uint32_t vb_i)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->next_vb_i = vb_i ? vb_i : dd->next_vb_i+1;

    if (flag.show_threads) dispatcher_show_time ("Generate vb", -1, dd->next_vb_i);

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

    if (flag.show_threads) dispatcher_show_time ("Start compute", dd->next_thread_to_dispatched, th->vb->vblock_i);

    ASSERT0 (dd->next_vb->vblock_i, "Error in dispatcher_compute: vb_i=0");

    if (dd->max_threads > 1) {
        unsigned err = pthread_create (&th->thread_id, NULL, dispatcher_thread_entry, th);
        ASSERT (!err, "Error: failed to create thread for next_vb_i=%u, err=%u", dd->next_vb->vblock_i, err);

        dd->next_thread_to_dispatched = (dd->next_thread_to_dispatched + 1) % dd->max_threads;
    }
    else  
        func(dd->next_vb); // single thread

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

    if (flag.show_threads) dispatcher_show_time ("Wait for thread", dd->next_thread_to_be_joined, th->vb->vblock_i);

    if (dd->max_threads > 1) 
        // wait for thread to complete (possibly it completed already)
        pthread_join(th->thread_id, NULL);

    if (flag.show_threads) dispatcher_show_time ("Join (end compute)", dd->next_thread_to_be_joined, th->vb->vblock_i);

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

    if (flag.show_time) profiler_add (&evb->profile, &dd->next_vb->profile);

    vb_release_vb (dd->next_vb); 
    dd->next_vb = NULL;
}

void dispatcher_finalize_one_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (dd->processed_vb) {

        buf_test_overflows(dd->processed_vb, "dispatcher_finalize_one_vb"); // just to be safe, this isn't very expensive

        if (flag.show_time) profiler_add (&evb->profile, &dd->processed_vb->profile);

        vb_release_vb (dd->processed_vb); // cleanup vb and get it ready for another usage (without freeing memory)
        dd->processed_vb = NULL;
    }

    if (dd->prog == PROGRESS_PERCENT)
        dispatcher_show_progress (dispatcher);
}                           

void dispatcher_set_input_exhausted (Dispatcher dispatcher, bool exhausted)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    dd->input_exhausted = exhausted;

    if (exhausted) {
        vb_release_vb (dd->next_vb);
        dd->next_vb = NULL;

        dd->next_vb_i--; // we didn't use this vb_i
    }
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

// returns the number of VBs successfully outputted
uint32_t dispatcher_fan_out_task (const char *filename,  // NULL to continue with previous filename
                                  ProgressType prog,
                                  const char *prog_msg,  // used if prog=PROGRESS_MESSAGE 
                                  bool test_mode,
                                  bool force_single_thread, 
                                  DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output)
{
    Dispatcher dispatcher = dispatcher_init (force_single_thread ? 1 : global_max_threads, 0, test_mode, true, filename, prog, prog_msg);
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
                dispatcher_set_input_exhausted (dispatcher, true);
                dispatcher_finalize_one_vb (dispatcher); 
            }
        }

        // if no condition was met, we're either done, or we still have some threads processing that are not done yet.
        // we wait a bit to avoid a tight busy loop
        else usleep (5000); 

    } while (!dispatcher_is_done (dispatcher));

    return ((DispatcherData *)dispatcher)->next_vb_i;
}
