// ------------------------------------------------------------------
//   threads.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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
#include <sys/sysctl.h>
#else // LINUX
#include <sched.h>
#include <sys/sysinfo.h>
#endif
#endif
#include "genozip.h"
#include "buffer.h"
#include "mutex.h"
#include "flags.h"
#include "strings.h"
#include "threads.h"

static Buffer threads = EMPTY_BUFFER;
static Mutex threads_mutex = {};
static unsigned next_thread_number=0;
static pthread_t main_thread;

typedef struct {
    bool in_use;
    pthread_t thread;
    unsigned number; // unique number of thread in execution history
    const char *task_name;
    uint32_t vb_i;
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
    abort();
}

static void threads_sighup_handler (void) 
{
    threads_cancel_other_threads();
    abort();
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

#ifndef _WIN32

    // block some signals - we will wait for them explicitly in threads_signal_handler(). this is inherited by all threads
    sigset_t sigset;
    sigemptyset (&sigset);
    sigaddset (&sigset, SIGSEGV);
    sigaddset (&sigset, SIGHUP);
    sigaddset (&sigset, SIGUSR1);
    int err = pthread_sigmask(SIG_BLOCK, &sigset, NULL);
    ASSERT (!err, "pthread_sigmask failed: %s", strerror (err));

    threads_create (threads_signal_handler, 0, "signal_handler", THREADS_NO_VB);
#endif
}

bool threads_am_i_main_thread (void)
{
    return pthread_self() == main_thread;
}

int threads_create (void *(*func)(void *), void *arg, const char *task_name, uint32_t vb_i)
{
    if (flag.show_threads) 
        iprintf ("%s: Creating thread %u (vb_i=%s)\n", task_name, next_thread_number, vb_i != THREADS_NO_VB ? str_int_s (vb_i).s : "NONE");

    pthread_t thread;
    unsigned err = pthread_create (&thread, NULL, func, arg);
    ASSERT (!err, "failed to create thread task=\"%s\" vb_i=%d: error=%u", task_name, (int)vb_i, err);

    mutex_lock (threads_mutex);

    int thread_id;
    for (thread_id=0; thread_id < threads.len; thread_id++)
        if (!ENT (ThreadEnt, threads, thread_id)->in_use) break;

    buf_alloc (evb, &threads, 0, thread_id+1, ThreadEnt, 2, "threads");
    threads.len = MAX (threads.len, thread_id+1);

    *ENT (ThreadEnt, threads, thread_id) = (ThreadEnt){ 
        .in_use    = true, 
        .number    = next_thread_number++,
        .task_name = task_name, 
        .vb_i      = vb_i, 
        .thread    = thread 
    };
    
    mutex_unlock (threads_mutex);

    return thread_id;
}

// returns success if joined (which is always the case if blocking)
bool threads_join (int thread_id, bool blocking)
{
    ThreadEnt *ent = ENT (ThreadEnt, threads, thread_id);

    if (flag.show_threads) 
        iprintf ("%s: Wait for thread %u (vb_i=%s)\n", ent->task_name, ent->number, ent->vb_i != THREADS_NO_VB ? str_int_s (ent->vb_i).s : "NONE");

    // case: wait for thread to complete (possibly it completed already)
    if (blocking)
        pthread_join (ent->thread, NULL);
    else {
#ifdef _WIN32
        int err = _pthread_tryjoin (ent->thread, NULL);
#else
        int err = pthread_tryjoin_np (ent->thread, NULL);
#endif
        if (err == EBUSY) return false;
        ASSERT (!err, "Error in pthread_tryjoin_np: %s", strerror (err));
    }
    
    if (flag.show_threads) 
        iprintf ("%s: Joined thread %u (vb_i=%s)\n", ent->task_name, ent->number, ent->vb_i != THREADS_NO_VB ? str_int_s (ent->vb_i).s : "NONE");

    ent->in_use = 0; // recycle entry
    return true;
}

// kills all other threads
void threads_cancel_other_threads (void)
{
    mutex_lock (threads_mutex);

    ARRAY (ThreadEnt, th, threads);
    for (unsigned i=0; i < th_len; i++)
        if (th->in_use && th->thread != pthread_self()) {
            pthread_cancel (th->thread);    
            th->in_use = false;
        }

    mutex_unlock (threads_mutex);

    usleep (500000); // wait 0.5 second for the all threads to die. pthread_join here hangs on Windows (not tested on others)
}
