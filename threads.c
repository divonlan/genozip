// ------------------------------------------------------------------
//   threads.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include <pthread.h>
#include <string.h>
#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <execinfo.h>
#include <signal.h>
#ifdef __APPLE__
#include "compatibility/mac_gettime.h"
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
static Mutex threads_mutex = { .name = "threads_mutex-not-initialized" };
static pthread_t main_thread, writer_thread;
static bool writer_thread_is_set = false;

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

#ifdef __linux__

static void threads_sigsegv_handler (void) 
{
    iprintf ("Segmentation fault, exiting. Time: %s\n", str_time().s);
    threads_print_call_stack(); // this works ok on mac, but seems to not print function names on Linux
    threads_cancel_other_threads();
    exit (EXIT_SIGSEGV); // dumps core
}

static void threads_sighup_handler (void) 
{
    iprintf ("Process %u received SIGHUP, ignoring. Use 'kill -9' to kill the process. Time: %s\n", getpid(), str_time().s);
}

// explicitly catch signals, that we otherwise blocked
static void *threads_signal_handler (void *sigset)
{    
    while (1) {
        int sig=0, err=0;
        ASSERT (!(err = sigwait ((sigset_t *)sigset, &sig)), "sigwait failed: %s", strerror (err));

        switch (sig) {
            case SIGSEGV : threads_sigsegv_handler(); break;
            case SIGHUP  : threads_sighup_handler();  break;
            case SIGUSR1 : fprintf (stderr, "Caught signal SIGUSR1, showing memory\n");
                           buf_show_memory_handler(); break;
            case SIGUSR2 : fprintf (stderr, "Caught signal SIGUSR2, writing threads log\n");
                           threads_write_log (false); break;
            default      : ABORT ("Unexpected signal %s", strsignal (sig));
        }
    }
    return 0; // never reaches here
}

#endif

void threads_initialize (void)
{
    mutex_initialize (threads_mutex);
    mutex_initialize (log_mutex);
    main_thread = pthread_self(); // we can't put it in buffer yet, because evb is not yet initialized

    buf_alloc (evb, &log, 500, 1000000, char, 2, "log");

// only Linux, as in Mac sigwait sometimes returns "invalid argument" - TODO solve this
#ifdef __linux__

    // block some signals - we will wait for them explicitly in threads_signal_handler(). this is inherited by all threads
    sigset_t sigset;
    sigemptyset (&sigset);
    sigaddset (&sigset, SIGSEGV);
    sigaddset (&sigset, SIGHUP);
    sigaddset (&sigset, SIGUSR1);
    sigaddset (&sigset, SIGUSR2);
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

void threads_set_writer_thread (void)
{
    writer_thread_is_set = true;
    writer_thread = pthread_self();
}

void threads_unset_writer_thread (void)
{
    writer_thread_is_set = false;
}

bool threads_am_i_writer_thread (void)
{
    return writer_thread_is_set && pthread_self() == writer_thread;    
}

void threads_write_log (bool to_info_stream)
{
    if (!log.len)
        iprint0 ("\nNo thread log - activate with --debug-threads\n");

    else if (to_info_stream)
        iprintf ("\nThread log (activated by --debug-threads):\n%.*s", (int)log.len, log.data);
    
    else {
        char filename[100];
        sprintf (filename, "genozip.threads-log.%u", getpid());
        buf_dump_to_file (filename, &log, 1, false, false, true, false);
    }
}

// called by MAIN thread only
static void threads_log_by_thread_id (ThreadId thread_id, const ThreadEnt *ent, const char *event)
{
    bool has_vb = ent->vb_i != (uint32_t)-1;

    if (flag.show_threads)  {
        if (has_vb) iprintf ("%s: vb_i=%u vb_id=%u %s thread_id=%d pthread=%"PRIu64"\n", ent->task_name, ent->vb_i, ent->vb_id, event, thread_id, (uint64_t)ent->pthread);
        else        iprintf ("%s: %s: thread_id=%d pthread=%"PRIu64"\n", ent->task_name, event, thread_id, (uint64_t)ent->pthread);
    }
    
    if (flag.debug_threads) {
        mutex_lock (log_mutex);
        buf_alloc (evb, &log, 10000, 1000000, char, 2, "log");
        
        if (has_vb) bufprintf (evb, &log, "%s: vb_i=%u vb_id=%u %s thread_id=%d pthread=%"PRIu64"\n", ent->task_name, ent->vb_i, ent->vb_id, event, thread_id, (uint64_t)ent->pthread);
        else        bufprintf (evb, &log, "%s: %s thread_id=%d pthread=%"PRIu64"\n", ent->task_name, event, thread_id, (uint64_t)ent->pthread);
        mutex_unlock (log_mutex);
    }
}

// called by any thread
void threads_log_by_vb (ConstVBlockP vb, const char *task_name, const char *event, 
                        int time_usec /* optional */)
{
    if (flag.show_threads) {
        unsigned pthread = (unsigned)(((uint64_t)pthread_self()) % 100000); // 5 digits
        if (time_usec) {
            if (vb->compute_thread_id >= 0)
                iprintf ("%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%u compute_thread_time=%s usec\n", 
                        task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, pthread, str_uint_commas (time_usec).s);
            else
                iprintf ("%s: vb_i=%d vb_id=%d %s compute_thread_time=%s usec\n", 
                        task_name, vb->vblock_i, vb->id, event, str_uint_commas (time_usec).s);
        } else {
            if (vb->compute_thread_id >= 0)
                iprintf ("%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%u\n", 
                        task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, pthread);
            else
                iprintf ("%s: vb_i=%d vb_id=%d %s\n", 
                        task_name, vb->vblock_i, vb->id, event);
        }
    }

    if (flag.debug_threads) {
        mutex_lock (log_mutex);

        // if thread other than main allocates evb it could cause corruption. we allocating
        if (log.size - log.len < 500) {
            WARN ("FYI: Thread log is out of space, log.size=%u log.len=%u", (unsigned)log.size, (unsigned)log.len);
            threads_write_log (false);
            return;
        }

        if (time_usec)
            bufprintf (evb, &log, "%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64" compute_thread_time=%s usec\n", 
                        task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, (uint64_t)pthread_self(), str_uint_commas (time_usec).s);
        else
            bufprintf (evb, &log, "%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64"\n", 
                        task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, (uint64_t)pthread_self());

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
    threads.len = MAX_(threads.len, thread_id+1);

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
bool threads_join_do (ThreadId *thread_id, 
                      VBlockP vb, // optional: if given, will return false if VB is not ready, if NULL, join will block
                      const char *func)
{
    ASSERTMAINTHREAD;

    ASSERT (*thread_id != THREAD_ID_NONE, "called from %s: thread not created or already joined", func);

    mutex_lock (threads_mutex);
    const ThreadEnt ent = *ENT (ThreadEnt, threads, *thread_id); // make a copy as array be realloced
    mutex_unlock (threads_mutex);

    static ThreadId last_joining = THREAD_ID_NONE;
    if (*thread_id != last_joining) { // show only once in a busy wait
        threads_log_by_thread_id (*thread_id, &ent, "JOINING: Wait for");
        last_joining = *thread_id;
    }

    // non-blocking - return false if the VB is not ready yet
    if (vb && !vb->is_processed) return false;

    // wait for thread to complete (possibly it completed already)
    pthread_join (ent.pthread, NULL);
    
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
