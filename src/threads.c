// ------------------------------------------------------------------
//   threads.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include <pthread.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#pragma pack(push, imagehlp, 8)
#include <imagehlp.h>
// #ifdef DEBUG  
// #include <backtrace.h> // pacman -S mingw-w64-x86_64-libbacktrace
// #endif
#pragma pack(pop, imagehlp)
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
#include "piz.h"
#include "sections.h"
#include "arch.h"
#include "license.h"

static Buffer threads = {};
static Mutex threads_mutex = { .name = "threads_mutex-not-initialized" };
static Mutex print_call_stack_mutex = { };
pthread_t main_thread = 0; 
static pthread_t writer_thread, zriter_thread;
static bool writer_thread_is_set = false, zriter_thread_is_set = false;

static Buffer log = {}; // for debugging thread issues, activated with --debug-threads
static Mutex log_mutex = {};

typedef struct {
    bool in_use;
    bool canceled;
    pthread_t pthread;
    rom task_name;
    VBIType vb_i;
    VBID vb_id;
} ThreadEnt;

static rom __attribute__((unused)) threads_get_task_name (void)
{
    pthread_t pthread = pthread_self();

    if (pthread == main_thread) return "main";
    
    for_buf (ThreadEnt, ent, threads)
        if (ent->pthread == pthread) return ent->task_name; // note: this will also detect the writer thread

    return "unregistered";
}

#ifdef _WIN32

// #ifdef DEBUG // requires debug info 

// static void backtrace_error_cb (void *data, const char *msg, int errnum)
// {
//     WARN ("%s error: %s\n", (rom)data, msg);
// }

// static int frame;
// static int backtrace_stackwalk_cb (void *data, uintptr_t pc, rom filename, int lineno, rom function)
// {
//     fprintf (stderr, "xxx #%u pc=%p %s at %s:%u\n", frame++, pc, function, filename, lineno);

//     return 0;
// }   

// #endif

// doesn't work properly yet
static void threads_print_call_stack_do (CONTEXT thread_ctx)
{
// bug 985 
// #ifdef DEBUG // requires debug info    
//     static struct backtrace_state *state = NULL;
//     static StrTextSuperLong exe_name; // must be static

//     if (!state) {
//         exe_name = arch_get_executable(); 
//         state = backtrace_create_state (exe_name.s, true, backtrace_error_cb, "backtrace_create_state");
//     }

//     frame = 0;
//     backtrace_full (state, 0, backtrace_stackwalk_cb, backtrace_error_cb, "backtrace_full");
//     return;
// #endif

    // loosly following: https://theorangeduck.com/page/printing-stack-trace-mingw
    // issue: shows only Windows DLL symbols, not ours

    fprintf (stderr, "\nCall stack (%s thread):\n",  // note: SUBMIT is not a separate thread in Windows
               threads_am_i_main_thread()   ? "MAIN" 
             : threads_am_i_writer_thread() ? "WRITER"
             : threads_am_i_zriter_thread() ? "ZRITER"
             :                                threads_get_task_name());

//     HANDLE process = GetCurrentProcess();
//     HANDLE thread  = GetCurrentThread();

//     ASSGOTO (SymInitialize (process, NULL, true), "SymInitialize failed: %s", str_win_error());

//     SymSetOptions (SymGetOptions() | SYMOPT_LOAD_LINES);

// #ifdef _M_X64 // most modern Intel and AMD processors
//     DWORD image_file = IMAGE_FILE_MACHINE_AMD64;
//     STACKFRAME64 sf = { .AddrPC     = { .Mode = AddrModeFlat, .Offset = thread_ctx.Rip },
//                         .AddrStack  = { .Mode = AddrModeFlat, .Offset = thread_ctx.Rsp },
//                         .AddrFrame  = { .Mode = AddrModeFlat, .Offset = thread_ctx.Rbp/* Rsp*/ } };

// #elif defined _M_IA64 // Intel Itanium processors
//     DWORD image_file = IMAGE_FILE_MACHINE_IA64;
//     STACKFRAME64 sf = { .AddrPC     = { .Mode = AddrModeFlat, .Offset = thread_ctx.StIIP },
//                         .AddrStack  = { .Mode = AddrModeFlat, .Offset = thread_ctx.IntSp },
//                         .AddrBStore = { .Mode = AddrModeFlat, .Offset = thread_ctx.RsBSP },
//                         .AddrFrame  = { .Mode = AddrModeFlat, .Offset = thread_ctx.IntSp } };
// #else
//     fprintf (stderr, "Stack tracing not available for this CPU architecture"); // eg _M_IX86
// #endif
    
//     int frame=0; 
//     while (StackWalk64 (image_file, process, thread, &sf, &thread_ctx, NULL, SymFunctionTableAccess64, SymGetModuleBase64, NULL) &&
//            sf.AddrPC.Offset) {

//         union symbol {
//             IMAGEHLP_SYMBOL64;
//             char data[sizeof(IMAGEHLP_SYMBOL64) + MAX_SYM_NAME-1];
//         } symbol = { { .SizeOfStruct = sizeof (union symbol), .MaxNameLength = MAX_SYM_NAME } };

//         DWORD64 displacement = 0;
//         // ASSGOTO (
//         SymGetSymFromAddr64 (process, sf.AddrPC.Offset, &displacement, (PIMAGEHLP_SYMBOL64)&symbol);
//             // , "SymGetSymFromAddr64 failed: %s", str_win_error());

//         IMAGEHLP_LINE64 line = { .SizeOfStruct = sizeof (IMAGEHLP_LINE64) };
//         DWORD displacement2 = 0;
//         // ASSGOTO (
//         SymGetLineFromAddr64 (process, sf.AddrPC.Offset, &displacement2, &line);
//             // , "SymGetLineFromAddr64 failed: %s", str_win_error());

//         fprintf (stderr, "#%u %s at %s:%u\n", frame++, symbol.Name, line.FileName, (uint32_t)line.LineNumber);
//     }

// error: 
//     CloseHandle (thread);
//     SymCleanup (process);
}

void threads_print_call_stack (void) 
{
    fflush (stdout); fflush (stderr);

#if defined _M_X64 || defined _M_IA64
    CONTEXT current_thread_ctx = { .ContextFlags = CONTEXT_FULL };
    RtlCaptureContext (&current_thread_ctx); // note: GetThreadContext() doesn't work for the current thread

    threads_print_call_stack_do (current_thread_ctx);
#else
    fprintf (stderr, "Stack tracing not available for this CPU architecture"); // RtlCaptureContext not supported on _M_IX86
#endif

}

static LONG WINAPI windows_exception_handler (EXCEPTION_POINTERS *ep)
{
    threads_print_call_stack_do (*ep->ContextRecord);

    return EXCEPTION_EXECUTE_HANDLER;
}

static void threads_init_signal_handlers (void)
{
    SetUnhandledExceptionFilter (windows_exception_handler);
}

#else

void threads_print_call_stack (void) 
{
    DO_ONCE {
        fflush (stdout); fflush (stderr);
        
#       define STACK_DEPTH 100
        void *array[STACK_DEPTH];
        size_t count = backtrace (array, STACK_DEPTH);
        
        // note: need to pass -rdynamic to the linker in order to see function names
        fprintf (stderr, "\nCall stack "); // separate printing, in case next fprintf hangs if stack is bad
        fflush (stderr);

        fprintf (stderr, "(%s thread):\n", 
                am_i_submit()                ? "SUBMIT"
              : threads_am_i_main_thread()   ? "MAIN" 
              : threads_am_i_writer_thread() ? "WRITER"
              : threads_am_i_zriter_thread() ? "ZRITER" 
              :                                 threads_get_task_name());
        fflush (stderr);

        backtrace_symbols_fd (array, count, STDERR_FILENO);
        fflush (stderr);
    }
}

// signal handler of SIGINT (CTRL-C) - debug only 
static void noreturn threads_sigint_handler (int signum) 
{
    fprintf (stderr, "\n%s process %u Received SIGINT (usually caused by Ctrl-C):", global_cmd, getpid()); 

    threads_print_call_stack(); // this works ok on mac, but seems to not print function names on Linux

    // look for deadlocks
    if (!flag.show_time) fprintf (stderr, "Tip: use --show-time to see locked mutexes with Ctrl-C\n");

    __atomic_thread_fence (__ATOMIC_ACQUIRE);
    __atomic_signal_fence (__ATOMIC_ACQUIRE);

    mutex_who_is_locked(); // works only if --show_time
    buflist_who_is_locked();
    
    // flush remaining output (eg in genocat)
    fflush (stdout);
    fflush (stderr);

    exit (128 + SIGINT); 
}

// signal handler of SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGSYS
static void noreturn threads_bug_signal_handler (int signum) 
{
    iprintf ("\n\nSignal \"%s\" received, exiting. Time: %s%s%s", 
             strsignal (signum), str_time().s, cond_str(flags_command_line(), "\nCommand line: ", flags_command_line()), report_support_if_unexpected());
    
    threads_print_call_stack(); // this works ok on mac (in debug only), but seems to not print function names on Linux

    buflist_test_overflows_all_vbs ("threads_bug_signal_handler", false);

    exit (128 + signum); // convention for exit code due to signal - both Linux and MacOS
}

// signal handler of SIGUSR1 
static void threads_sigusr1_handler (int signum) 
{
    __atomic_thread_fence (__ATOMIC_ACQUIRE);
    __atomic_signal_fence (__ATOMIC_ACQUIRE);
    
    iprint0 ("\n\n");
    buflist_show_memory (false, 0, 0);
    
    ctx_show_zctx_big_consumers (info_stream);
}

// signal handler of SIGUSR2 
static void threads_sigusr2_handler (int signum) 
{
    __atomic_thread_fence (__ATOMIC_ACQUIRE);
    __atomic_signal_fence (__ATOMIC_ACQUIRE);

    iprint0 ("\n\n");
    threads_write_log (true);
}

static void threads_init_signal_handlers (void)
{
    signal (SIGSEGV, threads_bug_signal_handler);
    signal (SIGBUS,  threads_bug_signal_handler);
    signal (SIGFPE,  threads_bug_signal_handler);
    signal (SIGILL,  threads_bug_signal_handler);
    signal (SIGSYS,  threads_bug_signal_handler); // used on Mac, ignored on Linux

    signal (SIGUSR1, threads_sigusr1_handler);
    signal (SIGUSR2, threads_sigusr2_handler);

#ifndef santize_thread // looks like this signal handler gets called all the time when compiled with --sanitize-thread
    if (flag.debug_or_test) 
        signal (SIGINT, threads_sigint_handler);
#endif
}

#endif

// called by main thread at system initialization
void threads_initialize (void)
{
    mutex_initialize (threads_mutex);
    mutex_initialize (log_mutex);
    mutex_initialize (print_call_stack_mutex);

    main_thread = pthread_self(); 

    buf_set_promiscuous (&log, "log");
    buf_alloc (evb, &log, 500, 1000000, char, 2, NULL);
    
    buf_set_promiscuous (&threads, "threads");
    buf_alloc (evb, &threads, 0, global_max_threads + 3, ThreadEnt, 2, NULL);

    threads_init_signal_handlers();
}

void threads_finalize (void)
{
    // normally, we intentionally leak the signal handler thread, letting it cancel implicitly when
    // the process exits. Under valgrind, we cancel it explicitly so that we only display unintentional leaks, not this intentional one 
    if (arch_is_valgrind()) {
        buf_destroy (log);
        buf_destroy (threads);
    }
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

void threads_set_zriter_thread (void)
{
    zriter_thread_is_set = true;
    zriter_thread = pthread_self();
}

void threads_unset_zriter_thread (void)
{
    zriter_thread_is_set = false;
}

bool threads_am_i_zriter_thread (void)
{
    return zriter_thread_is_set && pthread_self() == zriter_thread;    
}

void threads_write_log (bool to_info_stream)
{
    if (!log.len)
        iprint0 ("\nNo thread log - activate with --debug-threads\n");

    else if (to_info_stream)
        iprintf ("\nThread log (activated by --debug-threads):\n%.*s", (int)log.len, log.data);
    
    else {
        char filename[100];
        snprintf (filename, sizeof (filename), "genozip.threads-log.%u", getpid());
        buf_dump_to_file (filename, &log, 1, false, false, true, false);
    }
}

// called by main and writer threads
static void threads_log_by_thread_id (ThreadId thread_id, const ThreadEnt *ent, rom event)
{
    bool has_vb = ent->vb_i != (uint32_t)-1;

    if (flag_show_threads)  {
        if (has_vb) iprintf ("%s: vb_i=%u vb_id=%d %s thread_id=%d pthread=%"PRIu64"\n", ent->task_name, ent->vb_i, ent->vb_id, event, thread_id, (uint64_t)ent->pthread);
        else        iprintf ("%s: %s: thread_id=%d pthread=%"PRIu64"\n", ent->task_name, event, thread_id, (uint64_t)ent->pthread);
    }
    
    if (flag.debug_threads) {
        mutex_lock (log_mutex);
        buf_alloc (NULL, &log, 10000, 1000000, char, 2, NULL);
        
        if (has_vb) bufprintf (NULL, &log, "%s: vb_i=%u vb_id=%d %s thread_id=%d pthread=%"PRIu64"\n", ent->task_name, ent->vb_i, ent->vb_id, event, thread_id, (uint64_t)ent->pthread);
        else        bufprintf (NULL, &log, "%s: %s thread_id=%d pthread=%"PRIu64"\n", ent->task_name, event, thread_id, (uint64_t)ent->pthread);
        mutex_unlock (log_mutex);
    }
}

// called by any thread
void threads_log_by_vb (ConstVBlockP vb, rom task_name, rom event, 
                        int time_usec /* optional */)
{
    if (flag_show_threads) {
        unsigned pthread = (unsigned)(((uint64_t)pthread_self()) % 100000); // 5 digits

        #define COMP (vb->comp_i != COMP_NONE ? " comp=" : ""), \
                     (vb->comp_i != COMP_NONE ? comp_name (vb->comp_i) : "")

        if (time_usec) {
            if (vb->compute_thread_id >= 0)
                iprintf ("%s: vb_i=%d%s%s vb_id=%d %s vb->compute_thread_id=%d pthread=%u compute_thread_time=%s usec\n", 
                        task_name, vb->vblock_i, COMP, vb->id, event, vb->compute_thread_id, pthread, str_int_commas (time_usec).s);
            else
                iprintf ("%s: vb_i=%d%s%s vb_id=%d %s compute_thread_time=%s usec\n", 
                        task_name, vb->vblock_i, COMP, vb->id, event, str_int_commas (time_usec).s);
        } else {
            if (vb->compute_thread_id >= 0)
                iprintf ("%s: vb_i=%d%s%s vb_id=%d %s vb->compute_thread_id=%d pthread=%u\n", 
                        task_name, vb->vblock_i, COMP, vb->id, event, vb->compute_thread_id, pthread);
            else
                iprintf ("%s: vb_i=%d%s%s vb_id=%d %s\n", 
                        task_name, vb->vblock_i, COMP, vb->id, event);
        }
    }

    if (flag.debug_threads) {
        mutex_lock (log_mutex);

        // if thread other than main allocates evb it could cause corruption. we allocating
        if (log.size - log.len < 500) {
            WARN ("FYI: Thread log is out of space, log.size=%u log.len=%u", (unsigned)log.size, log.len32);
            threads_write_log (false);
            return;
        }

        if (time_usec)
            bufprintf (NULL, &log, "%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64" compute_thread_time=%s usec\n", 
                       task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, (uint64_t)pthread_self(), str_int_commas (time_usec).s);
        else
            bufprintf (NULL, &log, "%s: vb_i=%d vb_id=%d %s vb->compute_thread_id=%d pthread=%"PRIu64"\n", 
                       task_name, vb->vblock_i, vb->id, event, vb->compute_thread_id, (uint64_t)pthread_self());

        mutex_unlock (log_mutex);
    }
}

// thread entry point
static void *thread_entry_caller (void *vb_)
{
    VBlockP vb = (VBlockP)vb_;
    pthread_setcanceltype (PTHREAD_CANCEL_ASYNCHRONOUS, NULL); // thread can be canceled at any time

    // wait for threads_create to complete updating VB
    mutex_wait (vb->ready_for_compute, true); 

    threads_log_by_vb (vb, vb->compute_task, "STARTED", 0);

    // wait for VB initialzation data to be visible to this thread
    __atomic_thread_fence (__ATOMIC_ACQUIRE); 

    TimeSpecType start_time, end_time; 
    clock_gettime(CLOCK_REALTIME, &start_time); 

    rom task_name = vb->compute_task; // save, as vb might be released by compute_func

    // call entry point
    vb->compute_func (vb); 

    clock_gettime(CLOCK_REALTIME, &end_time); 
    int time_usec = 1000000 *(end_time.tv_sec - start_time.tv_sec) + (int)((int64_t)end_time.tv_nsec - (int64_t)start_time.tv_nsec) / 1000;

    threads_log_by_vb (vb, task_name, "COMPLETED", time_usec);
    
    // wait for data written by this thread to be visible to other threads
    __atomic_thread_fence (__ATOMIC_RELEASE); 

    return NULL;
}

// called from main thread or writer threads only
ThreadId threads_create (void (*func)(VBlockP), VBlockP vb)
{
    if (vb) threads_log_by_vb (vb, vb->compute_task, "ABOUT TO CREATE", 0);

    mutex_lock (threads_mutex);

    int thread_id;
    for (thread_id=0; thread_id < threads.len; thread_id++)
        if (!B(ThreadEnt, threads, thread_id)->in_use) break;

    buf_alloc (NULL, &threads, 1, global_max_threads + 3, ThreadEnt, 2, "threads");
    threads.len = MAX_(threads.len, thread_id+1);

    // release all data (inc. VB initialization data) to be visible to the new thread
    __atomic_thread_fence (__ATOMIC_RELEASE); 

    // set thread stack size: Mac fails without this (in reconstruct and also in longr zip)
    pthread_attr_t tattr;
    int err;
    ASSERT (!(err = pthread_attr_init(&tattr)), "pthread_attr_init: %s", strerror (err)); // initialized with default attributes
    ASSERT (!(err = pthread_attr_setstacksize (&tattr, 4 MB)), "pthread_attr_setstacksize: %s", strerror (err)); 

    mutex_initialize (vb->ready_for_compute);
    mutex_lock (vb->ready_for_compute); // thread_entry_caller will wait on this until VB is initialized

    pthread_t pthread;
    err = pthread_create (&pthread, &tattr, thread_entry_caller, vb);
    ASSERT (!err, "pthread_create: task=\"%s\" vb_i=%d: %s", vb->compute_task, vb->vblock_i, strerror(err));

    *B(ThreadEnt, threads, thread_id) = (ThreadEnt){ 
        .in_use    = true, 
        .task_name = vb->compute_task, 
        .vb_i      = vb->vblock_i,
        .vb_id     = vb->id,
        .pthread   = pthread 
    };

    vb->compute_thread_id = thread_id; // assigned while thread_entry_caller is waiting on mutex
    vb->compute_func      = func;
    mutex_unlock (vb->ready_for_compute); // alloc thread_entry_caller to call entry point

    mutex_unlock (threads_mutex);

    threads_log_by_thread_id (thread_id, B(ThreadEnt, threads, thread_id), "CREATED");

    return thread_id;
}

// returns success if joined (which is always the case if blocking)
void threads_join_do (ThreadId *thread_id, rom expected_task, rom expected_task2, rom func)
{
    ASSERT (*thread_id != THREAD_ID_NONE, "called from %s: thread not created or already joined", func);

    mutex_lock (threads_mutex);
    const ThreadEnt ent = *B(ThreadEnt, threads, *thread_id); // make a copy as array be realloced
    mutex_unlock (threads_mutex);

    ASSERT (*thread_id < threads.len32, "called from %s: thread_id=%u ∉ [0,%u]", func, *thread_id, threads.len32);
    ASSERT (ent.task_name, "called from %s: entry for thread_id=%u has task_name=NULL", func, *thread_id);
    ASSERT (!strcmp (ent.task_name, expected_task) || !strcmp (ent.task_name, expected_task2), 
            "called from %s: Expected thread_id=%u to have task=\"%s\" or task=\"%s\", but it has \"%s\"", func, *thread_id, expected_task, expected_task2, ent.task_name);

    static ThreadId last_joining = THREAD_ID_NONE;
    if (*thread_id != last_joining) { // show only once in a busy wait
        threads_log_by_thread_id (*thread_id, &ent, "JOINING: Wait for");
        last_joining = *thread_id;
    }

    // wait for thread to complete (no wait if it completed already)
    PTHREAD_JOIN (ent.pthread, ent.task_name);
    
    // wait for data from this thread to arrive
    __atomic_thread_fence (__ATOMIC_ACQUIRE); 

    threads_log_by_thread_id (*thread_id, &ent, "JOINED");

    mutex_lock (threads_mutex);
    B(ThreadEnt, threads, *thread_id)->in_use = false;
    mutex_unlock (threads_mutex);

    *thread_id = THREAD_ID_NONE;
}

// kills all other threads
void threads_cancel_other_threads (void)
{
#ifndef _WIN32 // segfaults on Windows, see bug 708
    if (!threads_am_i_main_thread())
        pthread_cancel (main_thread);

    mutex_lock (threads_mutex);

    ARRAY (ThreadEnt, th, threads);

    // send cancellation request to all threads
    for (unsigned i=0; i < th_len; i++) 
        if (th[i].in_use && th[i].pthread != pthread_self() && !pthread_cancel (th[i].pthread))
            th[i].canceled = true;

    // give time for all threads to terminate. note: we don't use pthread_join() here because sometimes it hangs
    usleep (500000);

    // // join all threads, which happens after cancellation is processed 
    // for (unsigned i=0; i < th_len; i++) 
    //     if (th[i].canceled) {
    //         pthread_join (th[i].pthread, NULL);
    //         th[i].in_use = th[i].canceled = false;
    //     }
        
    mutex_unlock (threads_mutex);
#endif
}
