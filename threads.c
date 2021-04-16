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
static unsigned next_thread_number=0;
static pthread_t main_thread;

static Buffer log = EMPTY_BUFFER; // for debugging thread issues, activated with --debug-threads
static Mutex log_mutex = {};

typedef struct {
    bool in_use;
    pthread_t pthread;
    unsigned number; // unique number of thread in execution history
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

    if (flag.debug_threads) {
        buf_alloc (evb, &log, 0, 10000000, char, 1, "log");
        mutex_initialize (log_mutex);
    }

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

// called by MAIN thread only
static void threads_log_by_ent (const ThreadEnt *ent, const char *event)
{
    bool has_vb = ent->vb_i != (uint32_t)-1;

    if (flag.show_threads)  {
        if (has_vb) iprintf ("%s: vb_i=%u vb_id=%u %s thread=%u pthread=%u\n", ent->task_name, ent->vb_i, ent->vb_id, event, ent->number, (unsigned)ent->pthread);
        else        iprintf ("%s: %s: thread=%u pthread=%u\n", ent->task_name, event, ent->number, (unsigned)ent->pthread);
    }
    
    if (flag.debug_threads) {
        mutex_lock (log_mutex);
        if (has_vb) bufprintf (evb, &log, "%s: vb_i=%u vb_id=%u %s thread=%u pthread=%u\n", ent->task_name, ent->vb_i, ent->vb_id, event, ent->number, (unsigned)ent->pthread);
        else        bufprintf (evb, &log, "%s: %s thread %u pthread=%u\n", ent->task_name, event, ent->number, (unsigned)ent->pthread);
        mutex_unlock (log_mutex);
    }
}

// called by any thread
void threads_log_by_vb (ConstVBlockP vb, const char *task_name, const char *event)
{
    if (flag.show_threads)  
        iprintf ("%s: vb_i=%u vb_id=%u %s\n", task_name, vb->vblock_i, vb->id, event);
    
    if (flag.debug_threads) {
        mutex_lock (log_mutex);
        
        // if thread other than main allocates evb it could cause corruption. we allocating
        ASSERT0 (log.size - log.len < 100 && !threads_am_i_main_thread(), "Thread log is out of space");
        
        bufprintf (evb, &log, "%s: vb_i=%u vb_id=%u %s\n", task_name, vb->vblock_i, vb->id, event);
        mutex_unlock (log_mutex);
    }
}

static void (*thread_entry)(void *arg); // global - counting on threads_create to only be called from the main thread

static void *thread_entry_caller (void *arg)
{
    thread_entry (arg);

    return NULL;
}

// called from main thread only
int threads_create (void (*func)(void *), void *arg, const char *task_name, ConstVBlockP vb)
{
    ASSERT (evb, "evb is NULL. task=%s", task_name);
    ASSERT0 (threads_am_i_main_thread(), "threads_create can only be called from the main thread");

    pthread_t thread;
    thread_entry = func;
    unsigned err = pthread_create (&thread, NULL, thread_entry_caller, arg);
    ASSERT (!err, "failed to create thread task=\"%s\" vb_i=%d: %s", task_name, vb ? (int)vb->vblock_i : -1, strerror(err));

    mutex_lock (threads_mutex);

    int thread_id;
    for (thread_id=0; thread_id < threads.len; thread_id++)
        if (!ENT (ThreadEnt, threads, thread_id)->in_use) break;

    buf_alloc (evb, &threads, 1, global_max_threads + 3, ThreadEnt, 2, "threads");
    threads.len = MAX (threads.len, thread_id+1);

    ThreadEnt ent = (ThreadEnt){ 
        .in_use    = true, 
        .number    = next_thread_number++,
        .task_name = task_name, 
        .vb_i      = vb ? vb->vblock_i : (uint32_t)-1,
        .vb_id     = vb ? vb->id       : (uint32_t)-1,
        .pthread   = thread 
    };
    *ENT (ThreadEnt, threads, thread_id) = ent;

    mutex_unlock (threads_mutex);

    threads_log_by_ent (&ent, "CREATED");

    return thread_id;
}

// returns success if joined (which is always the case if blocking)
bool threads_join (int thread_id, bool blocking)
{
    mutex_lock (threads_mutex);
    const ThreadEnt ent = *ENT (ThreadEnt, threads, thread_id); // make a copy as array be realloced
    mutex_unlock (threads_mutex);

    threads_log_by_ent (&ent, "JOINING: Wait for");

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
        ASSERT (!err, "Error in pthread_tryjoin_np: thread=%u pthread=%u thread_id=%u task=%s vb_i=%u: %s", 
                ent.number, (unsigned)ent.pthread, thread_id, ent.task_name, ent.vb_i, strerror (err));
    }
    
    threads_log_by_ent (&ent, "JOINED:");

    mutex_lock (threads_mutex);
    ENT (ThreadEnt, threads, thread_id)->in_use = false;
    mutex_unlock (threads_mutex);

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
