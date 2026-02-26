// ------------------------------------------------------------------
//   error.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#include <tlhelp32.h>
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
#endif
#endif
#include "genozip.h"
#include "threads.h"
#include "flags.h"
#include "buffer.h"
#include "strings.h"
#include "mutex.h"
#include "context.h"
#include "arch.h"
#include "chrom.h"
#include "url.h"
#include "progress.h"
#include "tip.h"

// mechanism to print a specific message upon exception
static _Thread_local rom catch_msg=0, catch_func=0;
static _Thread_local uint32_t catch_line=0;

static int process_argc = 0;
static char **process_argv = 0; 

#ifdef _WIN32

// #ifdef DEBUG // requires debug info 

// static void backtrace_error_cb (void *data, const char *msg, int errnum)
// {
//     WARN ("%s error: %s\n", (rom)data, msg);
// }

// static int frame;
// static int backtrace_stackwalk_cb (void *data, uintptr_t pc, rom filename, int lineno, rom function)
// {
//     fprintf (stderr, "#%u pc=%p %s at %s:%u\n", frame++, pc, function, filename, lineno);

//     return 0;
// }   

// #endif

// doesn't work properly yet
static void error_print_call_stack_do (CONTEXT thread_ctx)
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

static void error_print_call_stack (void) 
{
    fflush (stdout); fflush (stderr);

#if defined _M_X64 || defined _M_IA64
    CONTEXT current_thread_ctx = { .ContextFlags = CONTEXT_FULL };
    RtlCaptureContext (&current_thread_ctx); // note: GetThreadContext() doesn't work for the current thread

    error_print_call_stack_do (current_thread_ctx);
#else
    fprintf (stderr, "Stack tracing not available for this CPU architecture"); // RtlCaptureContext not supported on _M_IX86
#endif

}

static LONG WINAPI windows_exception_handler (EXCEPTION_POINTERS *ep)
{
    if (catch_msg) {
        progress_newline(); 
        fprintf (stderr, "Error in %s:%u: %s\n", catch_func, catch_line, catch_msg);
    }

    else
        error_print_call_stack_do (*ep->ContextRecord);

    return EXCEPTION_EXECUTE_HANDLER;
}

static void error_init_signal_handlers (void)
{
    SetUnhandledExceptionFilter (windows_exception_handler);
}

#else

extern bool am_i_submit (void);

static void error_print_call_stack (void) 
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
              :                                threads_get_task_name());
        fflush (stderr);

        backtrace_symbols_fd (array, count, STDERR_FILENO);
        fflush (stderr);
    }
}

// signal handler of SIGINT (CTRL-C) - debug only 
static void noreturn signal_handler_INT (int signum) 
{
    fprintf (stderr, "\n%s process %u Received SIGINT (usually caused by Ctrl-C):", global_cmd, getpid()); 

    error_print_call_stack(); // this works ok on mac, but seems to not print function names on Linux

    // look for deadlocks
    if (!flag.show_time) fprintf (stderr, _TIP "use --show-time to see locked mutexes with Ctrl-C\n");

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
static void noreturn signal_handler_bug (int signum) 
{
    if (catch_msg) {
        progress_newline(); 
        fprintf (stderr, "Error in %s:%u: %s\n", catch_func, catch_line, catch_msg);
    }

    else {
        iprintf ("\n\nSignal \"%s\" received, exiting. Time: %s%s%s", 
                 strsignal (signum), str_time().s, cond_str(flags_command_line(), "\nCommand line: ", flags_command_line()), report_support_if_unexpected());
        
        error_print_call_stack(); // this works ok on mac (in debug only), but seems to not print function names on Linux

        buflist_test_overflows_all_vbs ("signal_handler_bug", false);
    }

    exit (128 + signum); // convention for exit code due to signal - both Linux and MacOS
}

// signal handler of SIGUSR1 
static void signal_handler_USR1 (int signum) 
{
    __atomic_thread_fence (__ATOMIC_ACQUIRE);
    __atomic_signal_fence (__ATOMIC_ACQUIRE);
    
    iprint0 ("\n\n");
    buflist_show_memory (false, 0, 0);
    
    ctx_show_zctx_big_consumers (info_stream);
}

// signal handler of SIGUSR2 
static void signal_handler_USR2 (int signum) 
{
    __atomic_thread_fence (__ATOMIC_ACQUIRE);
    __atomic_signal_fence (__ATOMIC_ACQUIRE);

    iprint0 ("\n\n");
    threads_write_log (true);
}

static void error_init_signal_handlers (void)
{
    signal (SIGSEGV, signal_handler_bug);
    signal (SIGBUS,  signal_handler_bug);
    signal (SIGFPE,  signal_handler_bug);
    signal (SIGILL,  signal_handler_bug);
    signal (SIGSYS,  signal_handler_bug); // used on Mac, ignored on Linux

    signal (SIGUSR1, signal_handler_USR1);
    signal (SIGUSR2, signal_handler_USR2);

#ifndef santize_thread // looks like this signal handler gets called all the time when compiled with --sanitize-thread
    if (flag.debug_or_test) 
        signal (SIGINT, signal_handler_INT);
#endif
}

#endif

// set the error message to be displayed in case of an exception (segfault etc)
void catch_exception_do (rom msg, FUNCLINE)
{
    catch_func = func;  // note: catch_* are thread-local variables
    catch_line = code_line; 
    catch_msg  = msg; 
}

void uncatch_exception (void)
{
    catch_msg  = NULL; 
}

noreturn void stall (void) 
{ 
    while (1) 
        usleep (100000); // cancelation point for pthread_cancel per POSIX requirement (tested upon entry)
}

static void error_delete_zfile (bool is_error)
{
    STRli (save_name, strlen (z_file->name)+1);

    if (z_file && z_file->name) 
        strcpy (save_name, z_file->name);

    file_close (&z_file); // also frees file->name

    // note: logic to avoid a race condition causing the file not to be removed - if another thread seg-faults
    // because it can't access a z_file filed after z_file is freed, and signal_handler_bug aborts
    if (is_error && !flag.debug_or_test && file_exists (save_name)) 
        file_remove (save_name, true);
}

static void error_free_all (void)
{
    flags_finalize();
    chrom_finalize();
    ref_finalize (true);
    vb_destroy_pool_vbs (POOL_MAIN, true);
    vb_destroy_pool_vbs (POOL_BGZF, true);
    vb_destroy_vb (&evb);
}

extern rom command_name (void);

void noreturn error_exit (bool show_stack, bool is_error) 
{
    DO_ONCE { // prevent recursive entry due to a failed ASSERT in the cleanup process

        // case: skip freeing buffers since we're exiting any, UNLESS run under valgrind - in which case we free, so we can detect true memory leaks
        if (!arch_is_valgrind())
            flag.let_OS_cleanup_on_exit = true;  

        if (!is_error && (IS_ZIP || (IS_PIZ && flag.check_latest)/*non-returning test after compress*/)) {
            version_print_notice_if_has_newer();
            tip_print();
        }

        if (is_error && flag.debug_threads)
            threads_write_log (true);
            
        if (show_stack) 
            error_print_call_stack(); // this works ok on mac, but does not print function names on Linux (even when compiled with -g)

        // this is normally blocked, it runs only with certain flags (see code)
        buflist_test_overflows_all_vbs (is_error ? "exit_on_error" : "on_exit", true);

        if (flag.log_filename) {
            iprintf ("%s - execution ended %s\n", str_time().s, is_error ? "with an error" : "normally");
            fclose (info_stream);
        } 

        if (is_error) {
            close (1);   // prevent other threads from outputting to terminal (including buffered output), obscuring our error message
            close (2);
            url_close_remote_file_stream (NULL);  // <--- BREAKPOINT BRK
            file_kill_external_compressors(); 
        
            // cancel all other compute threads before closing z_file, so other threads don't attempt to access it 
            // (eg. z_file->data_type) and get a segmentation fault.
            threads_cancel_other_threads();
        }

        // if we're in ZIP - delete z_file
        if (primary_command == ZIP && z_file && z_file->name && !flag_loading_auxiliary) 
            error_delete_zfile (is_error);

        if (is_error) // call after canceling the writing threads 
            file_put_data_abort();

        fflush (stdout);
        fflush (stderr);
        
        if (!is_error && !flag.explicit_quiet && (flag.echo || flag.test_i)) { 
            progress_newline(); 
            fprintf (stderr, "%s[%u]: %s\n", command_name(), getpid(), "Exiting : Success"); 
        }

        threads_finalize();

        // we normally intentionally leak the VBs at process termination, but if we're running valgrind we free so to not display intentional leaks
        if (arch_is_valgrind()) 
            error_free_all();
    }

    else
        sleep (10000); // wait to be killed by first thread that is executing in DO_ONCE

    exit (is_error ? EXIT_GENERAL_ERROR : EXIT_OK);
} 

void noreturn error_restart (rom add_cmd_option)
{
    if (flag.restarted) exit_on_error (true); // if this is already a restarted process - we don't restart again

    DO_ONCE { // prevent recursive entry due to a failed ASSERT in the cleanup process
        fprintf (stderr, "\nRecovering by restarting with %s\n\n", add_cmd_option); 

        // case: skip freeing buffers since we're exiting any, UNLESS run under valgrind - in which case we free, so we can detect true memory leaks
        if (!arch_is_valgrind())
            flag.let_OS_cleanup_on_exit = true;  

        if (flag.debug_threads)
            threads_write_log (true);

        // this is normally blocked, it runs only with certain flags (see code)
        buflist_test_overflows_all_vbs ("error_restart", true);

        if (flag.log_filename) {
            iprintf ("%s - execution errored, restarting with %s\n", str_time().s, add_cmd_option);
            fclose (info_stream);
        } 

        // prevent other threads from outputting to terminal (including buffered output), obscuring our error message
        int save_stdout = dup (1);
        int save_stderr = dup (2);
        close (1);   
        close (2);

        url_close_remote_file_stream (NULL);  
        file_kill_external_compressors(); 

        // cancel all other compute threads before closing z_file, so other threads don't attempt to access it 
        // (eg. z_file->data_type) and get a segmentation fault.
        threads_cancel_other_threads();

        // if we're in ZIP - delete z_file
        if (primary_command == ZIP && z_file && z_file->name && !flag_loading_auxiliary)
            error_delete_zfile (true);

        file_put_data_abort();

        dup2 (save_stdout, 1); // restore - inherited through execvp
        dup2 (save_stderr, 2);
        
        fflush (stdout);
        fflush (stderr);

#ifndef _WIN32
        rom argv[process_argc + 3];
        memcpy (argv, process_argv, process_argc * sizeof (rom));
        argv[process_argc]     = add_cmd_option; // extra command line option(s)
        argv[process_argc + 1] = "--restarted";  // prevent recursive restarting
        argv[process_argc + 2] = NULL;  // list terminator

        execvp (argv[0], (char **)argv); // doesn't normally return
        ABORT ("Failed to restart process: execvp(%s) failed: %s", argv[0], strerror (errno));

#else
        // free most memory, as this process will still be alive while the child process runs
        error_free_all(); 
        return_freed_memory_to_kernel();
        
        STARTUPINFO si = { .cb = sizeof (STARTUPINFO) };
        PROCESS_INFORMATION pi = {};
        rom exec = arch_get_genozip_executable().s;
        
        STRl(cmd, 32 KB - 1)=0; // max command line length per Windows API
        cmd_len += snprintf (&cmd[cmd_len], sizeof(cmd)-1 - cmd_len, "\"%s\" ", exec);
        cmd_len += snprintf (&cmd[cmd_len], sizeof(cmd)-1 - cmd_len, "--restarted %s ", add_cmd_option);

        for (int i=1; i < process_argc; i++)
            cmd_len += snprintf (&cmd[cmd_len], sizeof(cmd)-1 - cmd_len, "\"%s\" ", process_argv[i]);
        
        cmd[cmd_len] = 0;
        ASSERTINRANGX (cmd_len, 0, sizeof(cmd)-1); // not ideal, testing after memory potentially already overwritten...
        
        ASSERT (CreateProcessA (exec, cmd, NULL, NULL, true, NORMAL_PRIORITY_CLASS | CREATE_NEW_PROCESS_GROUP, NULL, NULL, &si, &pi),
                "Failed to restart process: \"%s\": %s", cmd, str_win_error());
                
        // wait for "background" process so we can return its exit code
        WaitForSingleObject (pi.hProcess, INFINITE);

        DWORD exit_code = 0;
        exit (GetExitCodeProcess (pi.hProcess, &exit_code) ? exit_code : 1);
#endif
    }

    else
        stall(); // wait to be killed by first thread that is executing in DO_ONCE
}

void error_initialize (int argc, char *argv[])
{
    info_stream  = stdout; // may be changed during initialization
    process_argv = argv;
    process_argc = argc;

    error_init_signal_handlers();
}
