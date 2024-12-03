// ------------------------------------------------------------------
//   dispatcher.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "dispatcher.h"
#include "vblock.h"
#include "file.h"
#include "progress.h"
#include "threads.h"
#include "segconf.h"
#include "arch.h"
#include "zip.h"
#include "txtheader.h"

#define RR(x) ((x) % d->max_threads)

#define MAX_COMPUTED_VBS 4096
typedef struct DispatcherData {
    rom task_name;
    rom preproc_task_name;
    VBlockPoolType pool_type;
    uint32_t pool_in_use_at_init;  // Used to verify that all VBs allocated by dispatcher are released by finish time
    uint32_t max_vb_id_so_far; 
    VBlockP vbs[MAX_COMPUTED_VBS]; // VBs currently in the pipeline (with a compute thread or just before or after)
    VBlockP processed_vb;     // processed VB returned to caller (VB moved from the "vbs" field to this field)

    bool input_exhausted;
    bool paused;
    bool out_of_order;        // completed VBs may be retrieved out of order

    int next_dispatched;      // index into vbs of VB generated but not yet dispatched to a thread. -1 otherwise.
    int last_dispatched;      // index into vbs of last VB dispatched to a thread
    int last_joined;          // index into vbs - last VB joined - -1 if none was joined yet 

    uint32_t num_running_compute_threads;
    uint32_t next_vb_i;
    uint32_t max_threads;
    enum { PROGRESS_PERCENT, PROGRESS_MESSAGE, PROGRESS_NONE } progress_type;
    StrTextLong filename;
    
    uint64_t progress;        // progress towards target_progress
    uint64_t target_progress; // progress reaches this, it is at 100%

    Timestamp init_timestamp;
    uint32_t total_compute_time; // in msec
} DispatcherData;

// variables that persist across multiple dispatchers run sequentially
static TimeSpecType start_time; // wallclock
static bool start_time_initialized = false;

static Dispatcher main_dispatcher = 0; // dispatcher that updates percentage progress

// called from main and also compute threads

void dispatcher_increment_progress (rom where, int64_t increment)
{
    DispatcherData *d = (DispatcherData *)main_dispatcher;
    if (!d) return; // not PROGRESS_PERCENT

    // updated from the main thread and writer thread
    if (increment > 0)      __atomic_add_fetch (&d->progress, (uint64_t)increment,  __ATOMIC_RELAXED);
    else if (increment < 0) __atomic_sub_fetch (&d->progress, (uint64_t)-increment, __ATOMIC_RELAXED);

    if (!threads_am_i_main_thread()) return; // only main thread updates the output message

    // update target
    if (IS_ZIP && !txt_file->est_num_lines)
        d->target_progress = zip_get_target_progress();

    // in unbind mode - dispatcher is not done if there's another component after this one
    bool done = dispatcher_is_done (main_dispatcher);

    progress_update (d->task_name, load_relaxed (d->progress), d->target_progress, done);
}

void dispatcher_start_wallclock (void)
{
    clock_gettime (CLOCK_REALTIME, &start_time);
    start_time_initialized = true;
}

Dispatcher dispatcher_init (rom task_name, 
                            rom preproc_task_name, // if different
                            VBlockPoolType pool_type, uint32_t max_threads, uint32_t previous_vb_i,
                            bool out_of_order, // completed VBs may be retrieved out of order
                            bool test_mode, 
                            rom filename, // filename, or NULL if filename is unchanged
                            uint64_t target_progress, // implies PROGRESS_PERCENT 
                            rom prog_msg) // implies progress_type=PROGRESS_MESSAGE
{
    DispatcherData *d    = (DispatcherData *)CALLOC (sizeof(DispatcherData));
    d->task_name         = task_name;
    d->preproc_task_name = preproc_task_name ? preproc_task_name : task_name; 
    d->next_vb_i         = previous_vb_i;  // used if we're binding files - the vblock_i will continue from one file to the next
    d->max_threads       = MIN_(max_threads, MAX_COMPUTED_VBS);
    d->progress_type     = prog_msg?PROGRESS_MESSAGE : target_progress?PROGRESS_PERCENT : PROGRESS_NONE;
    d->target_progress   = target_progress;
    d->pool_type         = pool_type;
    d->out_of_order      = out_of_order && (d->max_threads > 1);  // max_threads=1 implies processing in order
    d->last_joined       = (d->max_threads > 1) ? -1 /* none joined yet */ : 0;
    d->last_dispatched   = -1; // none dispatched yet
    d->next_dispatched   = -1; // none generated yet
    d->init_timestamp    = arch_timestamp();

    if (d->progress_type == PROGRESS_PERCENT)
        main_dispatcher = d;

    if (filename)
        strncpy (d->filename.s, filename, sizeof (StrTextLong)-1);

    ASSERT (max_threads <= global_max_threads, "expecting max_threads=%u <= global_max_threads=%u", max_threads, global_max_threads);

    // always create the pool based on global_max_threads, not max_threads, because it is the same pool for all fan-outs throughout the execution
    vb_create_pool (pool_type, pool_type == POOL_MAIN ? "MAIN" : "BGZF");

    if (filename) // note: when unbinding, we print this in dispatcher_resume() 
        progress_new_component (filename, prog_msg, test_mode, start_time_initialized ? &start_time : NULL); 

    d->pool_in_use_at_init = vb_pool_get_num_in_use (pool_type, NULL); // we need to return the pool after dispatcher in the condition we received it...

    return d;
}

// allow out of order joining of VBs (normally set in dispatcher_init)
void dispatcher_allow_out_of_order (Dispatcher d)
{
    store_relaxed (d->out_of_order, (bool)(d->max_threads > 1)); // max_threads=1 implies processing in order
}

void dispatcher_pause (Dispatcher d)
{
    ASSERT0 (!progress_has_component(), "expecting progress to be finalized");

    d->paused = true;
    d->input_exhausted = false;
    d->next_vb_i--;

    start_time_initialized = false;
}

// PIZ: reinit dispatcher, used when splitting a genozip file to its components, using a single dispatcher object
void dispatcher_resume (Dispatcher d, uint32_t target_progress, CompIType comp_i)
{
    d->input_exhausted = false;
    d->filename        = txtheader_get_txt_filename_from_section (comp_i);
    d->progress        = 0;
    d->target_progress = target_progress;

    if (d->paused) 
        progress_new_component (d->filename.s, "0\%", flag.test, start_time_initialized ? &start_time : NULL);    

    d->paused          = false;
}

uint32_t dispatcher_get_next_vb_i (Dispatcher d)
{
    return d->next_vb_i;
}

void dispatcher_set_task_name (Dispatcher d, rom task_name)
{
    d->task_name = task_name;
}

void dispatcher_calc_avg_compute_vbs (Dispatcher d)
{
    uint32_t dispatcher_lifetime = arch_time_lap (d->init_timestamp);
    
    profiler_set_avg_compute_vbs (dispatcher_lifetime ? ((double)d->total_compute_time / (double)dispatcher_lifetime) : 0);
}

void dispatcher_finish (Dispatcher *dd_p, uint32_t *last_vb_i, bool cleanup_after_me,
                        bool show_memory)
{
    if (! *dd_p) return; // nothing to do

    DispatcherData *d = *dd_p;

    // must be before memory release (in ZIP - show in final component of final file)
    if (show_memory) 
        buflist_show_memory (false, d->max_threads, d->max_vb_id_so_far);    

    // free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files. 
    // don't bother freeing (=save time) if this is the last file, unless we're going to test and need the memory
    if (cleanup_after_me) {
        if ((IS_ZIP && (IS_REF_INTERNAL || IS_REF_EXT_STORE)) ||
            (IS_PIZ && IS_REF_STORED_PIZ))
            ref_unload_reference();

        vb_release_vb (&evb, d->task_name); // reset memory 
    }
    
    // only relevant to the ZIP dispatcher
    if (last_vb_i) *last_vb_i = d->next_vb_i; // for continuing vblock_i count between subsequent bound files

    if (main_dispatcher == *dd_p) main_dispatcher = 0;

    int32_t id_in_use;
    uint32_t pool_in_use_at_finish = vb_pool_get_num_in_use (d->pool_type, &id_in_use);
    // note for piz: we have 1 at start (wvb) and 0 at finish
    ASSERT (pool_in_use_at_finish <= d->pool_in_use_at_init, "Dispatcher \"%s\" leaked VBs: pool_in_use_at_init=%u pool_in_use_at_finish=%u (one VB in use is vb_id=%d)",
            d->task_name, d->pool_in_use_at_init, pool_in_use_at_finish, id_in_use);

    FREE (*dd_p);

    start_time_initialized = false;
}

VBlockP dispatcher_generate_next_vb (Dispatcher d, VBIType vb_i, CompIType comp_i)
{
    ASSERT (d->next_dispatched == -1, "next_dispatched=%d - VB already generated but not dispatched yet (task=%s)", d->next_dispatched, d->task_name);

    // update next_dispatched to first next entry (round robin). note: if not out_of_order, this will always be the next entry
    for (uint32_t i=0, nd=RR(d->last_dispatched + 1); i < d->max_threads; i++, nd = RR(nd + 1)) 
        if (!d->vbs[nd]) {
            d->next_dispatched = nd;
            break;
        }

    ASSERT (d->next_dispatched >= 0, "VB array is full (task=%s)", d->task_name);

    d->next_vb_i = vb_i ? vb_i : d->next_vb_i+1;

    d->vbs[d->next_dispatched] = vb_get_vb (d->pool_type, d->task_name, d->next_vb_i, comp_i);
    d->max_vb_id_so_far = MAX_(d->max_vb_id_so_far, d->vbs[d->next_dispatched]->id);

    return d->vbs[d->next_dispatched];
}

void dispatcher_compute (Dispatcher d, void (*func)(VBlockP))
{
    VBlockP vb = d->vbs[d->next_dispatched];
    ASSERTNOTNULL (vb);
    ASSERT0 (vb->vblock_i, "dispatcher_compute: cannot compute a VB because vb->vblock_i=0");

    vb->start_compute_timestamp = arch_timestamp();

    if (d->max_threads > 1) 
        threads_create (func, vb);
    else  
        func (vb); // single thread

    d->last_dispatched = d->next_dispatched;
    d->next_dispatched = -1;
    d->num_running_compute_threads++;
}

void dispatcher_abandon_next_vb (Dispatcher d)
{
    VBlockP vb = d->vbs[d->next_dispatched];
    ASSERTNOTNULL (vb);
    ASSERT0 (vb->vblock_i, "dispatcher_compute: cannot abandon a VB because vb->vblock_i=0");

    d->processed_vb = vb;
    d->vbs[d->next_dispatched] = NULL;
    d->next_dispatched = -1;
}

// returns index into vbs of a VB ready for joining, or -1 if there is none. 
static int32_t get_joinable (DispatcherData *d, bool next_if_none)
{
    // search for a thread to join - start after last_joined
    for (uint32_t i=0, nj=RR(d->last_joined + 1); 
         i < (load_relaxed (d->out_of_order) ? d->max_threads : 1); 
         i++, nj = RR(nj + 1)) 

        if (d->vbs[nj] && vb_is_processed (d->vbs[nj])) 
            return nj;

    // note: if next_if_none and no VBs are processed, we return the next one
    return next_if_none ? RR(d->last_joined + 1) : -1; 
}

bool dispatcher_has_processed_vb (Dispatcher d, bool *is_final) 
{
    if (!d->num_running_compute_threads) return false; // no running compute threads 

    bool my_is_final = d->input_exhausted && 
                       d->next_dispatched == -1 && 
                       d->num_running_compute_threads == 1; // this is the last vb to be processed

    if (is_final) *is_final = my_is_final;

    return my_is_final || get_joinable (d, false) >= 0;
}

// returns the next processed VB, or NULL if non-blocking the VBlock is not ready yet, or no running compute threads
VBlockP dispatcher_get_processed_vb (Dispatcher d, bool *is_final, bool blocking)
{
    ASSERT0 (!blocking || !load_relaxed (d->out_of_order), "out_of_order requires non blocking");

    if (!d->num_running_compute_threads) {
        if (is_final) *is_final = d->input_exhausted && d->next_dispatched == -1;
        return NULL; // no running compute threads 
    }

    if (d->max_threads > 1) {
        int32_t nj = get_joinable (d, blocking);

        if (nj == -1) return NULL; // blocking=false and no processed VBs 

        threads_join2 (&d->vbs[nj]->compute_thread_id, d->preproc_task_name, d->task_name); // blocking if no processed VB yet and blocking=true
        d->last_joined = nj;
    }

    if (is_final)
        *is_final = d->input_exhausted &&   
                    d->next_dispatched == -1 &&          // no VB generated but not dispatched
                    d->num_running_compute_threads == 1; // this is the last vb to be processed

    // move VB from "vbs" array to processed_vb
    d->processed_vb = d->vbs[d->last_joined];
    d->vbs[d->last_joined] = NULL;
    d->num_running_compute_threads--;
    d->total_compute_time += arch_time_lap (d->processed_vb->start_compute_timestamp);

    return d->processed_vb; 
}

bool dispatcher_has_free_thread (Dispatcher d)
{
    return d->num_running_compute_threads < MAX_(1, d->max_threads);
}

static bool dispatcher_has_free_vb_slot (Dispatcher d)
{
    return d->num_running_compute_threads + (d->next_dispatched >= 0) < d->max_threads;
}

uint32_t dispatcher_get_num_running_compute_threads (Dispatcher d)
{
    return d->num_running_compute_threads;
}

void dispatcher_recycle_vbs (Dispatcher d, bool release_vb)
{
    START_TIMER;

    if (d->processed_vb) {

        if (release_vb) { 
            // WORKAROUND to bug 343: there is a race condition of unknown cause if flag.no_writer_thread=true (eg --coverage, --count) crashes
            if (flag.no_writer_thread && !flag.test && !strcmp (d->task_name, PIZ_TASK_NAME)) usleep (1000); 
            vb_release_vb (&d->processed_vb, d->task_name); // cleanup vb and get it ready for another usage (without freeing memory)
        }

        // case: VB dispatched to the writer thread, and released there
        else
            d->processed_vb = NULL;
    }

    if (d->progress_type == PROGRESS_PERCENT)
        dispatcher_increment_progress (0, 0);

    COPY_TIMER_EVB (dispatcher_recycle_vbs);
}                           

void dispatcher_set_no_data_available (Dispatcher d, bool abandon_next_vb, DispatchStatus dispatch_status)
{
    if (dispatch_status == DATA_EXHAUSTED)
        d->input_exhausted = true;

    if (abandon_next_vb) {
        ASSERT (d->next_dispatched >= 0, "there is no next VB (task=%s)", d->task_name);
        vb_release_vb (&d->vbs[d->next_dispatched], d->task_name);
        d->next_vb_i--; // we didn't use this vb_i
        d->next_dispatched = -1;
    }
}    

bool dispatcher_is_done (Dispatcher d)
{
    return d->input_exhausted 
        && d->next_dispatched == -1
        && !d->processed_vb  
        && !d->num_running_compute_threads;
}

bool dispatcher_is_input_exhausted (Dispatcher d)
{
    ASSERTNOTNULL (d);
    
    return d->input_exhausted;
}

Dispatcher dispatcher_fan_out_task (rom task_name,
                                    rom filename,   // NULL to continue with previous filename
                                    uint64_t target_progress, // used if progress_type=PROGRESS_PERCENT 
                                    rom prog_msg,   // implies progress_type=PROGRESS_MESSAGE 
                                    bool out_of_order,
                                    bool test_mode,
                                    bool force_single_thread, 
                                    uint32_t previous_vb_i, // used if binding file
                                    uint32_t idle_sleep_microsec,
                                    bool free_when_done,
                                    DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output)
{
    Dispatcher d = dispatcher_init (task_name, NULL, POOL_MAIN, force_single_thread ? 1 : global_max_threads, 
                                    previous_vb_i, out_of_order, test_mode, filename, target_progress, 
                                    (prog_msg || target_progress) ? prog_msg : "0%");

    do {
        VBlockP next_vb = (d->next_dispatched >= 0) ? d->vbs[d->next_dispatched] : NULL;
        bool has_vb_ready_to_compute = next_vb && (next_vb->dispatch == READY_TO_COMPUTE);
        bool has_free_thread = dispatcher_has_free_thread (d);
        bool can_generate = !next_vb && !dispatcher_is_input_exhausted (d);
        bool has_free_vb_slot = dispatcher_has_free_vb_slot (d);

        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
        if (has_vb_ready_to_compute && has_free_thread) 
            dispatcher_compute (d, compute);

        // PRIORITY 2: output completed vbs, so they can be released and re-used
        else if (dispatcher_has_processed_vb (d, NULL)         ||  // case 1: there is a VB who's compute processing is completed
                 (has_vb_ready_to_compute && !has_free_thread) ||  // case 2: a VB ready to dispatch but all compute threads are occupied. wait here for one to complete
                 (can_generate && !has_free_vb_slot)) {            // case 3: we can read another VB but we have no slot for it

            VBlockP processed_vb = dispatcher_get_processed_vb (d, NULL, !load_relaxed (d->out_of_order)); // this will block until one is available
            if (!processed_vb) continue; // no running compute threads 

            if (output) output (processed_vb);
            
            dispatcher_recycle_vbs (d, true);
        }        
        
        // PRIORITY 3: If there is no VBs available to compute or to output, but input is not exhausted yet - get one range
        else if (can_generate) {

            next_vb = dispatcher_generate_next_vb (d, 0, COMP_NONE);

            prepare (next_vb);
            __atomic_thread_fence (__ATOMIC_RELEASE);

            if (next_vb->dispatch != READY_TO_COMPUTE) 
                dispatcher_set_no_data_available (d, true, next_vb->dispatch);
        }

        // if no condition was met, we're either done, or we still have some threads processing that are not done yet.
        // we wait a bit to avoid a tight busy loop
        else {
            if (d == main_dispatcher) dispatcher_increment_progress (0, 0); // ZIP: update progress

            START_TIMER;
            usleep (idle_sleep_microsec); 
            if (!strcmp (d->task_name, ZIP_TASK_NAME)) COPY_TIMER_EVB (zip_main_loop_idle);
        }

    } while (!dispatcher_is_done (d));

    if (free_when_done) FREE(d);
    
    // make sure memory writes by compute threads are visible to the main thread (not sure if this is needed or does pthread_join already do this)
    __atomic_thread_fence (__ATOMIC_ACQUIRE); 
 
    return d;
}
