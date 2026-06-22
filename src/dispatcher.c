// ------------------------------------------------------------------
//   dispatcher.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
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

rom _task_names[NUM_TASKS] = TASK_NAMES;
#define ធ(x) (task == TASK_##x)
#define TASK_IS_ZIP (d->task == TASK_ZIP)
#define TASK_IS_PIZ (d->task == TASK_PIZ)

#define RR(x) ((x) % d->max_threads)

#define MAX_COMPUTED_VBS 4096
typedef struct DispatcherData {
    Task task;
    VBlockPoolType pool_type;
    uint32_t pool_in_use_at_init;  // Used to verify that all VBs allocated by dispatcher are released by finish time
    uint32_t max_vb_id_so_far; 
    VBlockP vbs[MAX_COMPUTED_VBS]; // VBs currently in the pipeline (with a compute thread or just before or after)
    VBlockP processed_vb;     // processed VB returned to caller (VB moved from the "vbs" field to this field)

    bool input_exhausted;
    bool paused;
    DispatcherJoinMode join_mode; // whether completed VBs may be joined out of order

    int next_dispatched;      // index into vbs of VB generated but not yet dispatched to a thread. -1 otherwise.
    int last_dispatched;      // index into vbs of last VB dispatched to a thread
    int last_joined;          // index into vbs - last VB joined - -1 if none was joined yet 

    int task_pc_of_total;     // progress: how much of the total work on the file is account for by this dispatcher
    uint64_t progress;        // progress towards target_progress
    uint64_t target_progress; // progress reaches this, it is at 100% of its task_pc_of_total 

    uint32_t num_running_compute_threads;
    uint32_t next_vb_i;
    uint32_t max_threads;
    StrText1K filename;
    
    Timestamp init_timestamp;
    uint32_t total_compute_time; // in msec
} DispatcherData;

// variables that persist across multiple dispatchers run sequentially
static TimeSpecType start_time; // wallclock
static bool start_time_initialized = false;

// progress counter: zip: zeroed for new z_file, piz: zeroed for new txt_file
static enum { PROG_UNSTARTED, PROG_PREPROC_MSG, PROG_PERCENT } progress_stage;
static int pc_of_total_completed; // what % of total was completed by pervious dispatchers on the same file (integer arithmetic to avoid rounding errors when comparing to 100%)
static DispatcherData *curr;  // POOL_MAIN dispatcher currently running (only one concurrent dispatcher allowed per pool) 

static double task_time[NUM_TASKS] = {}; // time in seconds of each task

Task task_by_name (rom name)
{
    for (Task task=1; task < NUM_TASKS; task++)
        if (!strcmp (_task_names[task], name)) 
            return task;

    ABORTINP ("Unrecognized task name: %s", name);
}

void dispatcher_new_progress_bar (void)
{
    curr = NULL;
    pc_of_total_completed = 0;
    progress_stage = PROG_UNSTARTED;
}

// called from main and also compute threads

void dispatcher_increment_progress (rom where, int64_t increment)
{
    DispatcherData *d = curr;

    if (!d || !d->task_pc_of_total) return; // not contributing to progresss

    // updated from the main thread and writer thread
    if (increment > 0)      add_relaxed (d->progress, increment);
    else if (increment < 0) sub_relaxed (d->progress, -increment);

    if (!threads_am_i_main_thread()) return; // only main thread updates the output message

    // update target
    if (TASK_IS_ZIP && !txt_file->est_num_lines)
        d->target_progress = zip_get_target_progress();

    double portion_of_task  = d->target_progress ? MIN_(1.0, ((double)load_relaxed (d->progress) / (double)d->target_progress)) : 0;
    double portion_of_total = ((double)pc_of_total_completed / 100.0) + portion_of_task * ((double)d->task_pc_of_total / 100.0);

    if (progress_stage == PROG_PERCENT)
        progress_update (d->task, portion_of_total, portion_of_task, false);
}

void dispatcher_start_wallclock (void)
{
    clock_gettime (CLOCK_REALTIME, &start_time);
    start_time_initialized = true;
}

// messages to display before progress % begins
static rom dispatcher_pre_progress_msg (Task task)
{
    switch (task) {
        case TASK_ZIP           : return flag.skip_segconf ? "Compressing (skipped segconf)..."/*can't calc target*/ : NULL;
        case TASK_LOAD_EXT_REF  : 
        case TASK_LOAD_REFHASH  : return ref_cache_is_populating() ? "Caching reference file" : "Reading reference file";
        case TASK_BAMASS_READ   :
        case TASK_BAMASS_LINK   : return "Inspecting BAM...";
        case TASK_SCAN_FOR_DEPN : return "Preprocessing...";
        default                 : return NULL;
    }       
}

// estimated % of time spend on the ZIP task: progress bar will allocate this % as the "100%" of each task
static int dispatcher_task_percent_of_total (Task task)
{
    if (IS_PIZ && ធ(PIZ)) // note: reading internal reference is about ~0.2% in a WGS BAM
        return 100;

    else if (IS_PIZ) // all PIZ other tasks
        return 0;

    else if (flag.seg_only) 
        return ធ(ZIP) ? 100 : 0;
    
    else if (IN_RANGX (flag.make_reference, MAKE_REF_TINY, MAKE_REF_LARGE))
        return ធ(ZIP)           ? 31 // very roughly divide equally to 4 tasks - varies a lot
             : ធ(MAKE_REF_COMP) ? 23
             : ធ(MRH_OCCUPY)    ? 23
             : ធ(MRH_COMPRESS)  ? 23
             :                   0;

    else if (flag.make_reference == MAKE_REF_MINIMAL) 
        return ធ(ZIP)           ? 72 // very roughly divide equally to 4 tasks - varies a lot
             : ធ(MAKE_REF_COMP) ? 28
             :                   0;

    else if (z_file && IS_REF_INTERNAL) // note: not IS_REF_EXT_STORE as sections are copied as-is from external reference, so very fast
        return ធ(ZIP)           ? 97
             : ធ(COMPRESS_REF)  ? 3
             :                    0;

    else 
        return ធ(ZIP) ? 100 : 0;
}

// tasks for which we don't free the dispatcher after ending
static bool dispatcher_keep_after_ending (Task task)
{
    return ធ(ZIP) || ធ(PIZ) || 
           (ធ(COMPRESS_REF) && (Z_DT(BAM) || Z_DT(SAM))); // see ref_compress_ref
}

static VBlockPoolType dispatcher_task_pool_type (Task task)
{
    switch (task) {
        case TASK_BGZF              : return POOL_BGZF;
        case TASK_READ_TXT_HEADER   :
        case TASK_UNCOMP_RECON_PLAN : return POOL_MISC; // small tasks that are run while PIZ/ZIP are also running (in POOL_MAIN)
        default                     : return POOL_MAIN;
    }
}

Dispatcher dispatcher_init (Task task, 
                            uint32_t max_threads, uint32_t previous_vb_i,
                            DispatcherJoinMode join_mode, // completed VBs may be joined out of order
                            rom filename, // filename, or NULL if filename is unchanged
                            uint64_t target_progress) // of this work of THIS dispatcher. implies PROGRESS_PERCENT 
{
    DispatcherData *d   = (DispatcherData *)CALLOC (sizeof(DispatcherData));
    d->task             = task;
    d->next_vb_i        = previous_vb_i;  // used if we're binding files - the vblock_i will continue from one file to the next
    d->max_threads      = MIN_(max_threads, MAX_COMPUTED_VBS);
    d->target_progress  = target_progress;
    d->task_pc_of_total = dispatcher_task_percent_of_total (task);
    d->pool_type        = dispatcher_task_pool_type (task);
    d->join_mode      = (d->max_threads > 1) ? join_mode : JOIN_IN_ORDER;  // max_threads=1 implies processing in order
    d->last_joined      = (d->max_threads > 1) ? -1 /* none joined yet */ : 0;
    d->last_dispatched  = -1; // none dispatched yet
    d->next_dispatched  = -1; // none generated yet
    d->init_timestamp   = arch_timestamp();

    if (d->pool_type == POOL_MAIN) {
        ASSERTISNULL (curr); // only one POOL_MAIN dispatcher at a time allowed
        curr = d;
    }

    task_time[task] = arch_timestamp() / 1000000000.0; // seconds

    if (filename)
        strncpy (d->filename.s, filename, sizeof (StrText1K)-1);

    ASSERT (max_threads <= global_max_threads, "expecting max_threads=%u <= global_max_threads=%u", max_threads, global_max_threads);

    // always create the pool based on global_max_threads, not max_threads, because it is the same pool for all fan-outs throughout the execution
    vb_create_pool (d->pool_type);

    rom preproc_msg = dispatcher_pre_progress_msg (task);

    if (flag.show_tasks) 
        iprintf ("Start task %s (task_%%_of_total=%u)\n", task_name (d->task), d->task_pc_of_total);

    if (d->pool_type == POOL_MAIN) {
        if (progress_stage == PROG_UNSTARTED && filename) { // note: there might be short tasks with no filename before starting
            progress_new_component (filename, preproc_msg ? preproc_msg : "0\%", start_time_initialized ? &start_time : NULL); 
            progress_stage = preproc_msg ? PROG_PREPROC_MSG : PROG_PERCENT;
        }

        else if (progress_stage == PROG_PREPROC_MSG && !preproc_msg)
            progress_stage = PROG_PERCENT;
    }

    d->pool_in_use_at_init = vb_pool_get_num_in_use (d->pool_type, NULL); // we need to return the pool after dispatcher in the condition we received it...

    return d;
}

// allow out of order joining of VBs (normally set in dispatcher_init)
void dispatcher_allow_out_of_order (Dispatcher d)
{
    ASSERTNOTNULL (d);
    store_relaxed (d->join_mode, (d->max_threads > 1) ? JOIN_OUT_OF_ORDER : JOIN_IN_ORDER); // max_threads=1 implies processing in order
}

void dispatcher_pause (Dispatcher d)
{
    ASSERTNOTNULL (d);
    ASSERT0 (!progress_has_component(), "expecting progress to be finalized");

    d->paused = true;
    d->input_exhausted = false;
    d->next_vb_i--;

    start_time_initialized = false;
    pc_of_total_completed  = 0; // we will start a new progress bar upon resume
}

// PIZ: reinit dispatcher, used when splitting a genozip file to its components, using a single dispatcher object
void dispatcher_resume (Dispatcher d, uint32_t target_progress, CompIType comp_i)
{
    ASSERTNOTNULL (d);
    d->input_exhausted = false;
    d->filename        = txtheader_get_txt_filename_from_section (comp_i, false);
    d->progress        = 0;
    d->target_progress = target_progress;
    curr = d;

    if (d->paused) {
        progress_new_component (d->filename.s, "0\%", start_time_initialized ? &start_time : NULL);    
        d->paused = false;
    }
}

uint32_t dispatcher_get_next_vb_i (Dispatcher d)
{
    ASSERTNOTNULL (d);
    return d->next_vb_i;
}

void dispatcher_set_task (Dispatcher d, Task task)
{
    ASSERTNOTNULL (d);
    d->task = task;
}

void dispatcher_calc_avg_compute_vbs (Dispatcher d)
{
    ASSERTNOTNULL (d);
    uint32_t dispatcher_lifetime = arch_time_lap (d->init_timestamp);
    
    profiler_set_avg_compute_vbs (dispatcher_lifetime ? ((double)d->total_compute_time / (double)dispatcher_lifetime) : 0);
}

static void dispatcher_show_task_times (void)
{
    double total=0;
    for (Task task=0; task < NUM_TASKS; task++)
        if (task != TASK_BGZF) // BGZF runs in parallel
            total += task_time[task];

    iprint0 ("\nBreakdown of task times (in seconds):\n");
    for (Task task=0; task < NUM_TASKS; task++)
        if (task_time[task])
            iprintf ("%-20s  %5.2f\"  %1.1f%%\n", 
                     task_name(task), task_time[task], (task != TASK_BGZF ? percent (task_time[task], total) : 0));    
    
    iprintf ("%-20s  %5.2f\"  %1.1f%%\n", "TOTAL", total, 100.0);
}

VBlockP dispatcher_generate_next_vb (Dispatcher d, VBIType vb_i, CompIType comp_i)
{
    ASSERT (d->next_dispatched == -1, "next_dispatched=%d - VB already generated but not dispatched yet (task=%s)", 
            d->next_dispatched, task_name (d->task));

    // update next_dispatched to first next entry (round robin). note: if not out_of_order, this will always be the next entry
    for (uint32_t i=0, nd=RR(d->last_dispatched + 1); i < d->max_threads; i++, nd = RR(nd + 1)) 
        if (!d->vbs[nd]) {
            d->next_dispatched = nd;
            break;
        }

    ASSERT (d->next_dispatched >= 0, "VB array is full (task=%s)", task_name (d->task));

    d->next_vb_i = vb_i ? vb_i : d->next_vb_i+1;

    d->vbs[d->next_dispatched] = vb_get_vb (d->pool_type, d->task, d->next_vb_i, comp_i);
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
         i < (load_relaxed (d->join_mode) == JOIN_OUT_OF_ORDER ? d->max_threads : 1); 
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
    ASSERT0 (!blocking || load_relaxed (d->join_mode) != JOIN_OUT_OF_ORDER, "out_of_order requires non blocking");

    if (!d->num_running_compute_threads) {
        if (is_final) *is_final = d->input_exhausted && d->next_dispatched == -1;
        return NULL; // no running compute threads 
    }

    if (d->max_threads > 1) {
        int32_t nj = get_joinable (d, blocking);

        if (nj == -1) return NULL; // blocking=false and no processed VBs 

        threads_join (&d->vbs[nj]->compute_thread_id, d->task); // blocking if no processed VB yet and blocking=true
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
            if (flag.no_writer_thread && !flag.test && TASK_IS_PIZ) usleep (1000); 
            vb_release_vb (&d->processed_vb); // cleanup vb and get it ready for another usage (without freeing memory)
        }

        // case: VB dispatched to the writer thread, and released there
        else
            d->processed_vb = NULL;
    }

    if (progress_stage == PROG_PERCENT)
        dispatcher_increment_progress (0, 0);

    COPY_TIMER_EVB (dispatcher_recycle_vbs);
}                           

void dispatcher_set_no_data_available (Dispatcher d, bool abandon_next_vb, DispatchStatus dispatch_status)
{
    if (dispatch_status == DATA_EXHAUSTED)
        d->input_exhausted = true;

    if (abandon_next_vb) {
        ASSERT (d->next_dispatched >= 0, "there is no next VB (task=%s)", task_name (d->task));
        vb_release_vb (&d->vbs[d->next_dispatched]);
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

void dispatcher_finish (Dispatcher *dd_p, uint32_t *last_vb_i, bool cleanup_after_me,
                        bool show_memory)
{
    if (! *dd_p) return; // nothing to do

    DispatcherData *d = *dd_p;
    ASSERT0 (TASK_IS_ZIP || TASK_IS_PIZ, "this should only be called for the main task");

    // must be before memory release (in ZIP - show in final component of final file)
    if (show_memory) 
        buflist_show_memory (false, d->max_threads, d->max_vb_id_so_far);    

    // free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files. 
    // don't bother freeing (=save time) if this is the last file, unless we're going to test and need the memory
    if (cleanup_after_me) {
        if ((IS_ZIP && (IS_REF_INTERNAL || IS_REF_EXT_STORE)) ||
            (IS_PIZ && IS_REF_STORED_PIZ))
            ref_unload_reference();

        vb_release_vb (&evb); // reset memory 
    }
    
    // only relevant to the ZIP dispatcher
    if (last_vb_i) *last_vb_i = d->next_vb_i; // for continuing vblock_i count between subsequent bound files

    if (flag.show_tasks) {
        iprintf ("Finish component %s\n", d->filename.s);
        dispatcher_show_task_times();
    }
    
    start_time_initialized = false;
    curr = NULL;

    FREE (*dd_p);
}

void dispatcher_end_task (Dispatcher d)
{
    if (d->pool_type == POOL_MAIN) {
        pc_of_total_completed += d->task_pc_of_total;
        if (d->task_pc_of_total)
            progress_update (d->task, (double)pc_of_total_completed / 100.0, 1, pc_of_total_completed == 100);
        curr = 0;
    }

    task_time[d->task] = (arch_timestamp() / 1000000000.0) - task_time[d->task]; // time on this task (in seconds)

    int32_t id_in_use;
    uint32_t pool_in_use_at_finish = vb_pool_get_num_in_use (d->pool_type, &id_in_use);
    // note for piz: we have 1 at start (wvb) and 0 at finish
    ASSERT (pool_in_use_at_finish <= d->pool_in_use_at_init, "Dispatcher \"%s\" leaked VBs: pool_in_use_at_init=%u pool_in_use_at_finish=%u (one VB in use is vb_id=%d)",
            task_name (d->task), d->pool_in_use_at_init, pool_in_use_at_finish, id_in_use);

    if (flag.show_tasks) {
        progress_newline();
        iprintf ("End task %s 🕑=%1.2f seconds\n", task_name (d->task), task_time[d->task]);
    }

    if (!dispatcher_keep_after_ending (d->task)) // free dispatcher, unless caller will take care of it
        FREE (d);
}

Dispatcher dispatcher_fan_out_task (Task task,
                                    rom filename,             // NULL to continue with previous filename
                                    uint64_t target_progress, 
                                    DispatcherJoinMode join_mode,
                                    bool force_single_thread, 
                                    uint32_t previous_vb_i,   // used if binding file
                                    uint32_t idle_sleep_microsec,
                                    DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output)
{
    Dispatcher d = dispatcher_init (task, force_single_thread ? 1 : global_max_threads, 
                                    previous_vb_i, join_mode, filename, target_progress);

    do {
        VBlockP next_vb = (d->next_dispatched >= 0) ? d->vbs[d->next_dispatched] : NULL;
        bool has_vb_ready_to_compute = next_vb && (next_vb->dispatch == READY_TO_COMPUTE);
        bool has_free_thread         = dispatcher_has_free_thread (d);
        bool can_generate            = !next_vb && !dispatcher_is_input_exhausted (d);
        bool has_free_vb_slot        = dispatcher_has_free_vb_slot (d);

        // PRIORITY 1: is there a block available and a compute thread available? in that case dispatch it
        if (has_vb_ready_to_compute && has_free_thread) 
            dispatcher_compute (d, compute);

        // PRIORITY 2: output completed vbs, so they can be released and re-used
        else if (dispatcher_has_processed_vb (d, NULL)         ||  // case 1: there is a VB who's compute processing is completed
                 (has_vb_ready_to_compute && !has_free_thread) ||  // case 2: a VB ready to dispatch but all compute threads are occupied. wait here for one to complete
                 (can_generate && !has_free_vb_slot)) {            // case 3: we can read another VB but we have no slot for it

            VBlockP processed_vb = dispatcher_get_processed_vb (d, NULL, load_relaxed(d->join_mode) == JOIN_IN_ORDER); // this will block until one is available
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
            if (d->pool_type == POOL_MAIN) 
                dispatcher_increment_progress (0, 0); 

            START_TIMER;
            usleep (idle_sleep_microsec); 
            if (TASK_IS_ZIP) COPY_TIMER_EVB (zip_main_loop_idle);
        }

    } while (!dispatcher_is_done (d));

    bool is_misc = (d->pool_type == POOL_MISC); // save before destroying d
    Dispatcher ret = dispatcher_keep_after_ending (d->task) ? d : NULL;
    
    dispatcher_end_task (d); // destroys d for most tasks

    if (is_misc)
        vb_destroy_pool (POOL_MISC, true); // this was short-term pool for a specific task

    // make sure memory writes by compute threads are visible to the main thread (not sure if this is needed or does pthread_join already do this)
    __atomic_thread_fence (__ATOMIC_ACQUIRE); 
 
    return ret;
}

