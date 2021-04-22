// ------------------------------------------------------------------
//   threads.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#define __USE_GNU // needed for pthread_tryjoin_np on Linux
#include <pthread.h>
#undef __USE_GNU
#include <string.h>
#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <execinfo.h>
#include <signal.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#else // LINUX
#include <sched.h>
#include <sys/sysinfo.h>
#include <pthread.h>
#endif
#endif
#include "genozip.h"
#include "buffer.h"
#include "mutex.h"
#include "flags.h"
#include "strings.h"
#include "threads.h"
#include "vblock.h"

static Buffer threads = EMPTY_BUFFER;
static Mutex threads_mutex = {};
static pthread_t main_thread;

static Buffer log = EMPTY_BUFFER; // for debugging thread issues, activated with --debug-threads
static Mutex log_mutex = {};

typedef struct {
    bool in_use;
    pthread_t pthread;
    const char *task_name;
    uint32_t vb_i, vb_id;
} ThreadEnt;

void threads_print_call_stack (void) 
{
#if ! defined _WIN32 && defined DEBUG
#   define STACK_DEPTH 15
    void *array[STACK_DEPTH];
    size_t size = backtrace(array, STACK_DEPTH);
    
    fprintf (stderr, "Call stack:\n");
    backtrace_symbols_fd (array, size, STDERR_FILENO);
#endif
}

#ifndef _WIN32

static void threads_sigsegv_handler (void) 
{
    threads_print_call_stack(); // this works ok on mac, but seems to not print function names on Linux
    threads_cancel_other_threads();
    abort(); // dumps core
}

static void threads_sighup_handler (void) 
{
    threads_cancel_other_threads();
    exit(1);
}

// explicitly catch signals, that we otherwise blocked
static void *threads_signal_handler (void *sigset)
{    
    while (1) {
        int sig, err;
        ASSERT (!(err = sigwait ((sigset_t *)sigset, &sig)), "sigwait failed: %s", strerror (err));

        switch (sig) {
            case SIGSEGV : threads_sigsegv_handler();            break;
            case SIGHUP  : threads_sighup_handler();             break;
            case SIGUSR1 : buf_display_memory_usage_handler();   break;
            default      : ABORT ("Unexpected signal %s", strsignal (sig));
        }
    }
    return 0; // never reaches here`
}

#endif

void threads_initialize (void)
{
    mutex_initialize (threads_mutex);

    main_thread = pthread_self(); // we can't put it in buffer yet, because evb is not yet initialized

    // note: when this is called, flags are not set yet, so we just allocate regardless of flag.debug_threads
    mutex_initialize (log_mutex);

#ifndef _WIN32

    // block some signals - we will wait for them explicitly in threads_signal_handler(). this is inherited by all threads
    sigset_t sigset;
    sigemptyset (&sigset);
    sigaddset (&sigset, SIGSEGV);
    sigaddset (&sigset, SIGHUP);
    sigaddset (&sigset, SIGUSR1);
    int err = pthread_sigmask(SIG_BLOCK, &sigset, NULL);
    ASSERT (!err, "pthread_sigmask failed: %s", strerror (err));

    pthread_t signal_handler_thread;
    err = pthread_create (&signal_handler_thread, NULL, threads_signal_handler, &sigset); // can't use threads_create as evb isn't initialized yet
    ASSERT (!err, "failed to create signal_handler_thread: %s", strerror(err));
#endif
}

bool threads_am_i_main_thread (void)
{
    return pthread_self() == main_thread;
}

void threads_display_log (void)
{
    iprintf ("\nThread log (activated by --debug-threads):\n%.*s", (int)log.len, log.data);
}

// called by MAIN thread only
static void threads_log_by_thread_id (ThreadId thread_id, const ThreadEnt *ent, const char *event)
{
    bool has_vb = ent->vb_i != (uint32_t)-1;

    if (flag.show_threads)  {
        if (has_vb) iprintf ("%s: vb_i=%u vb_id=%u %s thread_id=%d pthread=%u\n", ent->task_name, ent->vb_i, ent->vb_id, event, thread_id, (unsigned)ent->pthread);
        else        iprintf ("%s: %s: thread_id=%d pthread=%u\n", ent->task_name, event, thread_id, (unsigned)ent->pthread);
    }
    
    if (flag.debug_threads) {
        mutex_lock (log_mutex);
        if (has_vb) bufprintf (evb, &log, "%s: vb_i=%u vb_id=%u %s thread_id=%d pthread=%u\n", ent->task_name, ent->vb_i, ent->vb_id, event, thread_id, (unsigned)ent->pthread);
        else        bufprintf (evb, &log, "%s: %s thread_id=%d pthread=%u\n", ent->task_name, event, thread_id, (unsigned)ent->pthread);
        mutex_unlock (log_mutex);
    }
}

// called by any thread
void threads_log_by_vb (ConstVBlockP vb, const char *task_name, const char *event, 
                        int time_usec /* optional */)
{
    if (flag.show_threads) {
        if (time_usec)
            iprintf ("%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64" compute_thread_time=%s usec\n", 
                     task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, pthread_self(), str_uint_commas (time_usec).s);
        else
            iprintf ("%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64"\n", 
                     task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, pthread_self());
    }

    if (flag.debug_threads) {
        mutex_lock (log_mutex);
        
        if (threads_am_i_main_thread()) // counting on the first call to be from main thread
            buf_alloc (evb, &log, 0, 1000000, char, 2, "log");
        else
            // if thread other than main allocates evb it could cause corruption. we allocating
            ASSERT (log.size - log.len > 100, "Thread log is out of space, log.size=%u log.len=%u", (unsigned)log.size, (unsigned)log.len);
        
        if (time_usec)
            bufprintf (evb, &log, "%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64" compute_thread_time=%s usec\n", 
                       task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, pthread_self(), str_uint_commas (time_usec).s);
        else
            bufprintf (evb, &log, "%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64"\n", 
                       task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, pthread_self());

        mutex_unlock (log_mutex);
    }
}

// thread entry point
static void *thread_entry_caller (void *vb_)
{
    VBlockP vb = (VBlockP)vb_;

    // wait for threads_create to complete updating VB
    mutex_wait (vb->vb_ready_for_compute_thread); 

    threads_log_by_vb (vb, vb->compute_task, "STARTED", 0);

    TimeSpecType start_time, end_time; 
    clock_gettime(CLOCK_REALTIME, &start_time); 

    const char *task_name = vb->compute_task; // save, as vb might be released by compute_func

    // call entry point
    vb->compute_func (vb); 

    clock_gettime(CLOCK_REALTIME, &end_time); 
    int time_usec = 1000000 *(end_time.tv_sec - start_time.tv_sec) + (int)((int64_t)end_time.tv_nsec - (int64_t)start_time.tv_nsec) / 1000;

    threads_log_by_vb (vb, task_name, "COMPLETED", time_usec);
    
    return NULL;
}

// called from main thread only
ThreadId threads_create (void (*func)(VBlockP), VBlockP vb)
{
    ASSERTMAINTHREAD;

    ASSERT (evb, "evb is NULL. task=%s", vb->compute_task);
    ASSERT0 (threads_am_i_main_thread(), "threads_create can only be called from the main thread");

    if (vb) threads_log_by_vb (vb, vb->compute_task, "ABOUT TO CREATE", 0);

    mutex_lock (threads_mutex);

    int thread_id;
    for (thread_id=0; thread_id < threads.len; thread_id++)
        if (!ENT (ThreadEnt, threads, thread_id)->in_use) break;

    buf_alloc (evb, &threads, 1, global_max_threads + 3, ThreadEnt, 2, "threads");
    threads.len = MAX (threads.len, thread_id+1);

    mutex_initialize (vb->vb_ready_for_compute_thread);
    mutex_lock (vb->vb_ready_for_compute_thread); // thread_entry_caller will wait on this until VB is initialized

    pthread_t pthread;
    unsigned err = pthread_create (&pthread, NULL, thread_entry_caller, vb);
    ASSERT (!err, "failed to create thread task=\"%s\" vb_i=%d: %s", vb->compute_task, vb->vblock_i, strerror(err));

    ThreadEnt ent = (ThreadEnt){ 
        .in_use    = true, 
        .task_name = vb->compute_task, 
        .vb_i      = vb->vblock_i,
        .vb_id     = vb->id,
        .pthread   = pthread 
    };
    *ENT (ThreadEnt, threads, thread_id) = ent;

    vb->compute_thread_id = thread_id; // assigned while thread_entry_caller is waiting on mutex
    vb->compute_func      = func;
    mutex_unlock (vb->vb_ready_for_compute_thread); // alloc thread_entry_caller to call entry point

    mutex_unlock (threads_mutex);

    threads_log_by_thread_id (thread_id, &ent, "CREATED");

    return thread_id;
}

// returns success if joined (which is always the case if blocking)
bool threads_join (ThreadId *thread_id, bool blocking)
{
    ASSERTMAINTHREAD;

    ASSERT0 (*thread_id != THREAD_ID_NONE, "thread not created or already joined");

    mutex_lock (threads_mutex);
    const ThreadEnt ent = *ENT (ThreadEnt, threads, *thread_id); // make a copy as array be realloced
    mutex_unlock (threads_mutex);

    static ThreadId last_joining = THREAD_ID_NONE;
    if (*thread_id != last_joining) { // show only once in a busy wait
        threads_log_by_thread_id (*thread_id, &ent, "JOINING: Wait for");
        last_joining = *thread_id;
    }

    // case: wait for thread to complete (possibly it completed already)
    if (blocking)
        pthread_join (ent.pthread, NULL);
    else {
#ifdef _WIN32
        int err = _pthread_tryjoin (ent.pthread, NULL);
#else
        int err = pthread_tryjoin_np (ent.pthread, NULL);
#endif
        if (err == EBUSY) return false;
        ASSERT (!err, "Error in pthread_tryjoin_np: pthread=%u thread_id=%u task=%s vb_i=%u: %s", 
                (unsigned)ent.pthread, *thread_id, ent.task_name, ent.vb_i, strerror (err));
    }
    
    threads_log_by_thread_id (*thread_id, &ent, "JOINED");

    mutex_lock (threads_mutex);
    ENT (ThreadEnt, threads, *thread_id)->in_use = false;
    mutex_unlock (threads_mutex);

    *thread_id = THREAD_ID_NONE;
    return true;
}

// kills all other threads
void threads_cancel_other_threads (void)
{
    mutex_lock (threads_mutex);

    ARRAY (ThreadEnt, th, threads);
    for (unsigned i=0; i < th_len; i++)
        if (th->in_use && th->pthread != pthread_self()) {
            pthread_cancel (th->pthread);    
            th->in_use = false;
        }

    mutex_unlock (threads_mutex);

    usleep (500000); // wait 0.5 second for the all threads to die. pthread_join here hangs on Windows (not tested on others)
}
