// ------------------------------------------------------------------
//   dispatcher.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#if defined __APPLE__ 
#include "compatibility/mac_gettime.h"
#endif

#include "genozip.h"
#include "dispatcher.h"
#include "vblock.h"
#include "file.h"
#include "profiler.h"
#include "progress.h"
#include "threads.h"
#include "segconf.h"

#define MAX_COMPUTED_VBS 4096
typedef struct {
    const char *task_name;
    unsigned max_vb_id_so_far; 
    VBlock *vbs[MAX_COMPUTED_VBS]; // VBs currently in the pipeline (with a compute thread or just before or after)
    VBlock *processed_vb; // processed VB returned to caller (VB moved from the "vbs" field to this field)

    bool input_exhausted;
    bool paused;

    unsigned next_dispatched; // index into vbs
    unsigned next_joined;     // index into vbs

    unsigned num_running_compute_threads;
    unsigned next_vb_i;
    unsigned max_threads;
    ProgressType prog;
    const char *filename;
    char *progress_prefix;
} DispatcherData;

// variables that persist across multiple dispatchers run sequentially
static TimeSpecType profiler_timer; // wallclock

static void dispatcher_show_progress (Dispatcher dispatcher)
{
    uint64_t total=0, sofar=0;
    
    // case: ZIP of plain txt files (including if decompressed by an external compressor) 
    // - we go by the amount of txt content processed 
    if (command == ZIP && txt_file->disk_size && file_is_plain_or_ext_decompressor (txt_file)) { 
        sofar = z_file->txt_data_so_far_single; 
        total = MAX (sofar, txtfile_get_seggable_size()); 
    } 
    
    // case ZIP: locally decompressed files - eg .vcf.gz or .fq.bz2 - 
    // we go by the physical disk size and how much has been consumed from it so far
    else if (command == ZIP && file_is_read_via_int_decompressor (txt_file)) {
        total = txt_file->disk_size; 
        sofar = txt_file->disk_so_far; 
    } 
        
    // case PIZ: physical z_file size, minus any bytes skipped 
    else if (command == PIZ || 
             (command == ZIP && file_is_read_via_int_decompressor (txt_file))) {
        total = z_file->disk_size_minus_skips; 
        sofar = z_file->disk_so_far; 
    } 
        
    // case: we have no idea what is the disk size of the VCF file - for example, because its STDIN, or because
    // its coming from a URL that doesn't provide the size
    else if (command == ZIP && !txt_file->disk_size)
        return; // we can't show anything if we don't know the file size
    
    else ABORT ("Error in dispatcher_show_progress: unsupported case: command=%u txt_file->type=%s", command, ft_name (txt_file->type));

    // in unbind mode - dispatcher is not done if there's another component after this one
    bool done = dispatcher_is_done (dispatcher);

    progress_update (&((DispatcherData *)dispatcher)->progress_prefix, sofar, total, done);
}

void dispatcher_start_wallclock (void)
{
    clock_gettime (CLOCK_REALTIME, &profiler_timer);
}

Dispatcher dispatcher_init (const char *task_name, unsigned max_threads, unsigned previous_vb_i,
                            bool test_mode, 
                            const char *filename, // filename, or NULL if filename is unchanged
                            ProgressType prog, const char *prog_msg /* used if prog=PROGRESS_MESSAGE */)   
{
    DispatcherData *dd   = (DispatcherData *)CALLOC (sizeof(DispatcherData));
    dd->task_name        = task_name;
    dd->next_vb_i        = previous_vb_i;  // used if we're binding files - the vblock_i will continue from one file to the next
    dd->max_threads      = MIN_(max_threads, MAX_COMPUTED_VBS);
    dd->prog             = prog;

    if (filename)
        dd->filename = filename;

    ASSERT (max_threads <= global_max_threads, "expecting max_threads=%u <= global_max_threads=%u", max_threads, global_max_threads);
    
    // always create the pool based on global_max_threads, not max_threads, because it is the same pool for all fan-outs throughout the execution
    vb_create_pool (MAX_(1, global_max_threads)    // compute thread VBs
                  + (command == PIZ)               // txt header VB (for PIZ) or 
                  + z_file->max_conc_writing_vbs + // writer thread VBs 
                  + 2);                            // background cache creation of (gref + primref) or (gref + gref refhash)
    if (!flag.unbind && filename) // note: for flag.unbind (in main file), we print this in dispatcher_resume() 
        dd->progress_prefix = progress_new_component (filename, prog_msg, test_mode); 

    return dd;
}

void dispatcher_pause (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    dd->paused = true;
    dd->input_exhausted = false;
    dd->next_vb_i--;
}

// PIZ: reinit dispatcher, used when splitting a genozip to its components, using a single dispatcher object
void dispatcher_resume (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->input_exhausted = false;
    dd->paused          = false;
    dd->filename        = txt_file->basename;
    
    dd->progress_prefix = progress_new_component (dd->filename, "0\%", -1);    
}

uint32_t dispatcher_get_next_vb_i (Dispatcher dispatcher)
{
    return ((DispatcherData *)dispatcher)->next_vb_i;
}

void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i, bool cleanup_after_me)
{
    if (! *dispatcher) return; // nothing to do

    DispatcherData *dd = (DispatcherData *)*dispatcher;

    // must be before vb_cleanup_memory() 
    if (flag.show_memory) buf_show_memory (false, dd->max_threads, dd->max_vb_id_so_far);    

    // free memory allocations between files, when compressing multiple non-bound files or 
    // decompressing multiple files. 
    // don't bother freeing (=save time) if this is the last file, unless we're going to test and need the memory
    if (cleanup_after_me) {
        vb_cleanup_memory(); 
        vb_release_vb (&evb); // reset memory 
    }
    
    // only relevant to the ZIP dispatcher
    if (last_vb_i) *last_vb_i = dd->next_vb_i; // for continuing vblock_i count between subsequent bound files

    FREE (dd->progress_prefix);
    FREE (*dispatcher);
}

VBlock *dispatcher_generate_next_vb (Dispatcher dispatcher, uint32_t vb_i)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    dd->next_vb_i = vb_i ? vb_i : dd->next_vb_i+1;

    dd->vbs[dd->next_dispatched] = vb_get_vb (dd->task_name, dd->next_vb_i);
    dd->max_vb_id_so_far = MAX_(dd->max_vb_id_so_far, dd->vbs[dd->next_dispatched]->id);

    return dd->vbs[dd->next_dispatched];
}

void dispatcher_compute (Dispatcher dispatcher, void (*func)(VBlockP))
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    VBlockP vb = dd->vbs[dd->next_dispatched];
    ASSERTNOTNULL (vb);
    ASSERT0 (vb->vblock_i, "dispatcher_compute: cannot compute a VB because vb->vblock_i=0");

    if (dd->max_threads > 1) {
        threads_create (func, vb);
        dd->next_dispatched = (dd->next_dispatched + 1) % dd->max_threads;
    }
    else  
        func(vb); // single thread

    dd->num_running_compute_threads++;
}

void dispatcher_abandon_next_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    VBlockP vb = dd->vbs[dd->next_dispatched];
    ASSERTNOTNULL (vb);
    ASSERT0 (vb->vblock_i, "dispatcher_compute: cannot dont_compute a VB because vb->vblock_i=0");

    dd->processed_vb = vb;
    dd->vbs[dd->next_dispatched] = NULL;
}

bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final) 
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (!dd->num_running_compute_threads) return false; // no running compute threads 

    VBlockP vb = dd->vbs[dd->next_joined];

    bool my_is_final = dd->input_exhausted && 
                       !dd->vbs[dd->next_dispatched] && 
                       dd->num_running_compute_threads == 1; // this is the last vb to be processed

    if (is_final) *is_final = my_is_final;

    return my_is_final || (vb && vb->is_processed);
}

// returns the next processed VB, or NULL if non-blocking the VBlock is not ready yet, or no running compute threads
VBlock *dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final, bool blocking)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (!dd->num_running_compute_threads) return NULL; // no running compute threads 

    if (dd->max_threads > 1) 
        // wait for thread to complete (possibly it completed already)
        if (!threads_join (&dd->vbs[dd->next_joined]->compute_thread_id, blocking ? NULL : dd->vbs[dd->next_joined]))
            return NULL; // only happens if non-blocking

    // move VB from "vbs" array to processed_vb
    dd->processed_vb = dd->vbs[dd->next_joined];
    dd->vbs[dd->next_joined] = NULL;
    dd->num_running_compute_threads--;

    dd->next_joined = (dd->next_joined + 1) % MAX_(1, dd->max_threads);

    return dd->processed_vb; 
}

bool dispatcher_has_free_thread (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    return dd->num_running_compute_threads < MAX_(1, dd->max_threads);
}

bool dispatcher_has_active_threads (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    return dd->num_running_compute_threads > 0;
}

VBlock *dispatcher_get_next_vb (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    return dd->vbs[dd->next_dispatched];
}

void dispatcher_recycle_vbs (Dispatcher dispatcher, bool release_vb)
{
    START_TIMER;

    DispatcherData *dd = (DispatcherData *)dispatcher;

    if (dd->processed_vb) {

        if (release_vb) { 
            // WORKAROUND to bug 343: there is a race condition of unknown cause is flag.no_writer_thread=true (eg --coverage, --count) crashes
            if (flag.no_writer_thread && !flag.test && !strcmp (dd->task_name, "piz")) usleep (1000); 
            vb_release_vb (&dd->processed_vb); // cleanup vb and get it ready for another usage (without freeing memory)
        }

        // case: VB dispatched to the writer thread, and released there
        else
            dd->processed_vb = NULL;
    }

    if (dd->prog == PROGRESS_PERCENT)
        dispatcher_show_progress (dispatcher);

    COPY_TIMER_VB (evb, dispatcher_recycle_vbs);
}                           

void dispatcher_set_input_exhausted (Dispatcher dispatcher, bool exhausted)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;
    dd->input_exhausted = exhausted;

    if (exhausted) {
        vb_release_vb (&dd->vbs[dd->next_dispatched]);
        dd->next_vb_i--; // we didn't use this vb_i
    }
}    

bool dispatcher_is_done (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    return dd->input_exhausted 
        && !dd->vbs[dd->next_dispatched]
        && !dd->processed_vb  
        && !dd->num_running_compute_threads;
}

bool dispatcher_is_input_exhausted (Dispatcher dispatcher)
{
    DispatcherData *dd = (DispatcherData *)dispatcher;

    return dd->input_exhausted;
}

// returns the number of VBs successfully outputted
Dispatcher dispatcher_fan_out_task (const char *task_name,
                                    const char *filename,   // NULL to continue with previous filename
                                    ProgressType prog,
                                    const char *prog_msg,   // used if prog=PROGRESS_MESSAGE 
                                    bool test_mode,
                                    bool force_single_thread, 
                                    uint32_t previous_vb_i, // used if binding file
                                    uint32_t idle_sleep_microsec,
                                    DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output)
{
    Dispatcher dispatcher = dispatcher_init (task_name, force_single_thread ? 1 : global_max_threads, 
                                             previous_vb_i, test_mode, filename, prog, prog_msg ? prog_msg : "0%");
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

            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL, true); // this will block until one is available
            if (!processed_vb) continue; // no running compute threads 

            if (output) output (processed_vb);
            
            dispatcher_recycle_vbs (dispatcher, true);
        }        
        
        // PRIORITY 3: If there is no variant block available to compute or to output, but input is not exhausted yet - get one range
        else if (!next_vb && !dispatcher_is_input_exhausted (dispatcher)) {

            next_vb = dispatcher_generate_next_vb (dispatcher, 0);

            prepare (next_vb);

            if (!next_vb->ready_to_dispatch) 
                dispatcher_set_input_exhausted (dispatcher, true);
        }

        // if no condition was met, we're either done, or we still have some threads processing that are not done yet.
        // we wait a bit to avoid a tight busy loop
        else usleep (idle_sleep_microsec); 

    } while (!dispatcher_is_done (dispatcher));

    return dispatcher;
}
