// ------------------------------------------------------------------
//   stream.c
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include <errno.h>
#include <sys/types.h>
#include <stdarg.h>

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#include <fcntl.h> // for _O_BINARY
#else // non-Windows
#include <sys/ioctl.h>
#include <sys/wait.h>
#endif

#ifdef __APPLE__
#include <signal.h> // for kill()
#endif

#ifdef __linux__ 
extern int fcntl (int __fd, int __cmd, ...); // defined in fcntl.h but not linux/fcntl.h
#include <linux/fcntl.h> // for F_SETPIPE_SZ - if needed for porting, you can delete this include (the code with F_SETPIPE_SZ is #ifdefed)
#endif

#include "stream.h"
#include "file.h"
#include "url.h"
#include "buffer.h"

typedef struct stream_ {
    rom exec_name;
    StreamP substream; // when an input_decompressor itself has a curl input stream
    FILE *from_stream_stdout;
    FILE *from_stream_stderr;
    FILE *to_stream_stdin;
#ifdef _WIN32
    HANDLE pid;
    DWORD exit_status;
#else
    pid_t pid;
    int exit_status; // if the process already exited and we got its code, the code will be here and pid==0;
#endif
} Stream;

#define EXIT_OK                   0
#define EXIT_GENERAL_ERROR        1
#define EXIT_INVALID_GENOZIP_FILE 2
#define EXIT_DOWNSTREAM_LOST      3
#define EXIT_STREAM               4
#define EXIT_SIGHUP               5
#define EXIT_SIGSEGV              6

#define EXIT_ABNORMAL             7

rom exit_code_name (int exit_code)
{
    static rom names[NUM_EXIT_CODES] = 
        { "EXIT_OK", "EXIT_GENERAL_ERROR", "EXIT_INVALID_GENOZIP_FILE", "EXIT_DOWNSTREAM_LOST", "EXIT_STREAM",
          "EXIT_SIGHUP", "EXIT_SIGSEGV", "EXIT_ABNORMAL" };

#ifndef _WIN32 
    if (WIFSIGNALED(exit_code))
        return strsignal (WTERMSIG(exit_code)); // name of signal
#endif

    return (exit_code >= 0 && exit_code < NUM_EXIT_CODES) ? names[exit_code] : "Invalid exit code";
}

#ifdef _WIN32
static rom stream_windows_error (void)
{
    char *msg;
    FormatMessage (FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                   NULL, GetLastError(),
                   MAKELANGID (LANG_NEUTRAL, SUBLANG_DEFAULT),
                   (LPTSTR)&msg, 0, NULL);
    return msg;
}
#endif

static void stream_pipe (int *fds, uint32_t pipe_size, bool is_stream_to_genozip)
{
#ifdef _WIN32
    ASSERT (!_pipe (fds,  pipe_size, _O_BINARY), "failed to create pipe: %s", strerror (errno));

#else // Not Windows
    ASSERT (!pipe (fds), "failed to create pipe: %s (errno=%u)", strerror (errno), errno);

#ifdef F_SETPIPE_SZ // this is Linux-specific, and available since kernel 2.6.35
    
    // lower pipe size requested to the maximum allowed system size (fcntl will return EPERM if size is larger than limit)
    FILE *pipe_max_size_file = fopen ("/proc/sys/fs/pipe-max-size", "rb");
    if (pipe_max_size_file) {
        char max_pipe_size_str[30] = {};
        size_t bytes = fread (max_pipe_size_str, 1, 29, pipe_max_size_file);
        if (bytes) pipe_size = MIN_(pipe_size, atoi (max_pipe_size_str));
        FCLOSE (pipe_max_size_file, "pipe_max_size_file");
    }

    // set the pipe size
    fcntl(fds[0], F_SETPIPE_SZ, pipe_size); // ignore errors
#endif // F_SETPIPE_SZ

#endif // not Windows
}

static void stream_abort_cannot_exec (rom exec_name, rom reason)
{
    url_kill_curl(NULL);

    fprintf (stderr, "\n%s: %s, %s needs to be in the execution path.\n", global_cmd, reason, exec_name);  

    if (!strstr (exec_name, "genozip")) // this is NOT genozip run from main_test_after_genozip
        fprintf (stderr, "Note that %s is a separate software package that is not affiliated with Genozip in any way.\n", exec_name);  

    exit (EXIT_STREAM); // special code - if we are a child existing, the parent will catch this error code
}

#ifdef _WIN32
void stream_set_inheritability (int fd, bool is_inheritable)
{
    ASSERT (SetHandleInformation ((HANDLE)_get_osfhandle (fd), HANDLE_FLAG_INHERIT, is_inheritable), 
            "SetHandleInformation failed: %s", stream_windows_error());
}

static HANDLE stream_exec_child (int *stream_stdout_to_genozip, int *stream_stderr_to_genozip, int *genozip_to_stream_stdin,
                                 FILE *redirect_stdout_file, FILE *redirect_stdin_pipe, 
                                 unsigned argc, char * const *argv, rom reason)
{
    int cmd_line_len = 0;
    for (int i=0; i < argc; i++)
        if (argv[i]) cmd_line_len += strlen (argv[i]) + 3; // +1 for the space (or \0) after +2 for the quotes

    char cmd_line[cmd_line_len];
    char *fmt_tmpl = "%s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s %s%s%s"; // # of %s = MAX_ARGC
    char fmt[MAX_ARGC * 7];
    strcpy (fmt, fmt_tmpl);
    fmt[argc*7 - 1] = '\0'; // cut short the fmt string

    // surround each argument with quotes just in case it contains spaces
#   define Q(i) argv[i] ? "\"" : "", argv[i], argv[i] ? "\"" : ""
    sprintf (cmd_line, fmt, Q(0), Q(1), Q(2), Q(3), Q(4), Q(5), Q(6), Q(7), Q(8), Q(9), Q(10), Q(11), Q(12), Q(13), Q(14), 
             Q(15), Q(16), Q(17), Q(18), Q(19), Q(20), Q(21), Q(22), Q(23), Q(24), Q(25), Q(26), Q(27), Q(28), Q(29));
#   undef Q

    STARTUPINFO startup_info = {};
    startup_info.cb         = sizeof startup_info;
    startup_info.hStdError  = stream_stderr_to_genozip ? (HANDLE)_get_osfhandle (stream_stderr_to_genozip[1]) : GetStdHandle (STD_ERROR_HANDLE);
    startup_info.dwFlags    = STARTF_USESTDHANDLES;

    // determine where the child's stdout will go 
    if (stream_stdout_to_genozip)  startup_info.hStdOutput = (HANDLE)_get_osfhandle (stream_stdout_to_genozip[1]);
    else if (redirect_stdout_file) startup_info.hStdOutput = (HANDLE)_get_osfhandle (fileno (redirect_stdout_file));
    else                           startup_info.hStdOutput = GetStdHandle (STD_OUTPUT_HANDLE);

    // determine where the child's stdin comes from 
    if (genozip_to_stream_stdin)   startup_info.hStdInput = (HANDLE)_get_osfhandle (genozip_to_stream_stdin[0]);
    else if (redirect_stdin_pipe)  startup_info.hStdInput = (HANDLE)_get_osfhandle (fileno (redirect_stdin_pipe));
    else                           startup_info.hStdInput = INVALID_HANDLE_VALUE;

    // ensure that "our" side of the pipe isn't inherited by the child process
    if (stream_stderr_to_genozip) stream_set_inheritability (stream_stderr_to_genozip[0], false);
    if (stream_stdout_to_genozip) stream_set_inheritability (stream_stdout_to_genozip[0], false);
    if (genozip_to_stream_stdin)  stream_set_inheritability (genozip_to_stream_stdin[1],  false);

    PROCESS_INFORMATION proc_info = {}; // the output of CreateProcess goes in here

    bool success = CreateProcess (NULL, cmd_line, NULL, NULL, TRUE, NORMAL_PRIORITY_CLASS, 
                                  NULL, NULL, &startup_info, &proc_info);
    if (!success) stream_abort_cannot_exec (argv[0], reason);
    
    CloseHandle (proc_info.hThread); // child process's main thread

    // recover inheritability in case we need it in future
    if (stream_stderr_to_genozip) stream_set_inheritability (stream_stderr_to_genozip[0], true);
    if (stream_stdout_to_genozip) stream_set_inheritability (stream_stdout_to_genozip[0], true);
    if (genozip_to_stream_stdin)  stream_set_inheritability (genozip_to_stream_stdin[1],  true);

    return proc_info.hProcess;
}

#else // not Windows

void stream_set_inheritability (int fd, bool is_inheritable) {} // does nothing

static pid_t stream_exec_child (int *stream_stdout_to_genozip, int *stream_stderr_to_genozip, int *genozip_to_stream_stdin,
                                FILE *redirect_stdout_file, FILE *redirect_stdin_pipe, 
                                unsigned argc, char *const *argv, rom reason)
{
    pid_t child_pid = fork();
 
    if (child_pid == -1 && errno == ENOMEM)
        buflist_show_memory (true, 0, 0);

    ASSERT (child_pid != -1, "fork() failed when attempting to run %s: %s (errno=%d)", argv[0], strerror (errno), errno);

    if (child_pid) return child_pid; // parent returns

    // determine what happens to the child's stdin
    
    if (redirect_stdin_pipe) // redirect child stdin from a pipe to another process
        dup2 (fileno (redirect_stdin_pipe), STDIN_FILENO);
 
    else if (genozip_to_stream_stdin) { // redirect child stdin from a pipe to genozip
        dup2 (genozip_to_stream_stdin[0], STDIN_FILENO);
        CLOSE (genozip_to_stream_stdin[1], "genozip_to_stream_stdin[1]", false); // child closes writing side of the pipe  
    }

    else CLOSE (STDIN_FILENO, "STDIN_FILENO", true); // no stdin for child. quiet bc stdin may already be closed - if we just finished reading txt input from redirected stdin

    // determine what happens to the child's stdout
    if (redirect_stdout_file)  // redirect child stdout to a file
        dup2 (fileno (redirect_stdout_file), STDOUT_FILENO);
 
    else if (stream_stdout_to_genozip) {   // redirect child stdout to a pipe to genozip
        dup2 (stream_stdout_to_genozip[1], STDOUT_FILENO);
        close (stream_stdout_to_genozip[0]); // ignore errors
    }
    else {
        // stdout continues to go to be shared with genozip (e.g. the terminal)
    }

    // determine what happens to the child's stderr
    
    if (stream_stderr_to_genozip) { // redirect child stdout to a pipe to genozip
        dup2 (STDERR_FILENO, 3); // keep the original stderr in fd=3 in case we can't execute and need to report an error
        dup2 (stream_stderr_to_genozip[1], STDERR_FILENO);
        close (stream_stderr_to_genozip[0]); // ignore errors
    } 
    else {
        // stderr continues to go to be shared with genozip (e.g. the terminal)
    }

    execvp (argv[0], (char * const *)argv); // if successful, doesn't return

    if (stream_stderr_to_genozip) dup2 (3, STDERR_FILENO); // get back the original stderr

    stream_abort_cannot_exec (argv[0], reason);
    
    return 0; // squash compiler warning - never reach here
}

#endif

StreamP stream_create (StreamP parent_stream, uint32_t from_stream_stdout, uint32_t from_stream_stderr, uint32_t to_stream_stdin,
                       FILE *redirect_stdout_file, 
                       rom input_url_name, // input to exec is coming from a URL
                       bool input_stdin,   // input to exec is coming from stdin 
                       rom reason,
                       rom exec_name, ...)
{
    ASSERT0 (!from_stream_stdout || !redirect_stdout_file, "cannot redirect child output to both genozip and a file");
    ASSERT0 (!to_stream_stdin    || !input_url_name,       "cannot redirect child input from both genozip and a url");

    int stream_stdout_to_genozip[2], stream_stderr_to_genozip[2], genozip_to_stream_stdin[2];
    
    if (from_stream_stdout) stream_pipe (stream_stdout_to_genozip, from_stream_stdout, true);
    if (from_stream_stderr) stream_pipe (stream_stderr_to_genozip, from_stream_stderr, true);
    if (to_stream_stdin)    stream_pipe (genozip_to_stream_stdin,  to_stream_stdin, false);

    Stream *stream = CALLOC (sizeof (Stream));
    stream->exec_name = exec_name;
    if (parent_stream) parent_stream->substream = stream;

    // copy our function arguments to argv
    va_list argp;
    va_start (argp, exec_name);

    rom argv[MAX_ARGC] = {};
    unsigned argc = 0;

    argv[argc++] = exec_name;

    rom arg;
    while (argc < MAX_ARGC) {
        arg = va_arg (argp, rom );
        if (arg == SKIP_ARG) continue;
        if (arg == NULL) break;
        argv[argc++] = arg;
    }

    ASSERT (!arg, "too many arguments when executing %s - max is %u", exec_name, MAX_ARGC-1); // MAX_ARGC-1 real args and last must be NULL

    FILE *input_pipe = NULL;
    if (input_url_name) 
        input_pipe = url_open (stream, input_url_name);
    else if (input_stdin)
        input_pipe = fdopen (STDIN_FILENO,  "rb");

    // execute the child in a separate process
    stream->pid = stream_exec_child (from_stream_stdout ? stream_stdout_to_genozip : NULL, 
                                     from_stream_stderr ? stream_stderr_to_genozip : NULL, 
                                     to_stream_stdin    ? genozip_to_stream_stdin  : NULL, 
                                     redirect_stdout_file,
                                     input_pipe,
                                     argc, (char * const *)argv,
                                     reason); 

    if (redirect_stdout_file) 
        FCLOSE (redirect_stdout_file, "redirect_stdout_file"); // the child has this file open, we don't need it

    if (input_pipe && input_url_name) 
        url_disconnect_from_curl (&input_pipe); // curl now belongs to the child process - we can close the stream to it 

    if (input_pipe && !input_url_name) 
        FCLOSE (input_pipe, "input_pipe"); // the child has this pipe open, we don't need it

    // store our side of the pipe, and close the side used by the child
    if (from_stream_stdout) {
        stream->from_stream_stdout = fdopen (stream_stdout_to_genozip[0], READ);        
        CLOSE (stream_stdout_to_genozip[1], "stream_stdout_to_genozip[1]", false);
    }

    if (from_stream_stderr) {
        stream->from_stream_stderr = fdopen (stream_stderr_to_genozip[0], READ);
        CLOSE (stream_stderr_to_genozip[1], "stream_stderr_to_genozip[1]", false);
    }

    if (to_stream_stdin) {
        stream->to_stream_stdin    = fdopen (genozip_to_stream_stdin[1], WRITE);
        CLOSE (genozip_to_stream_stdin[0], "genozip_to_stream_stdin[0]", false);
    }

    return stream;
}

void stream_close_pipes (Stream *stream)
{
    // fclose and fail silently (pipes might be gone already after we killed the counter part process?)
    if (stream->from_stream_stdout) FCLOSE (stream->from_stream_stdout, "stream->from_stream_stdout");
    if (stream->from_stream_stderr) FCLOSE (stream->from_stream_stderr, "stream->from_stream_stderr");
    if (stream->to_stream_stdin)    FCLOSE (stream->to_stream_stdin,    "stream->to_stream_stdin");
}

int stream_close (Stream **stream, StreamCloseMode close_mode)
{
    if (! *stream) return 0; // nothing to do

    if ((*stream)->substream) 
        stream_close (&(*stream)->substream, close_mode);

    stream_close_pipes (*stream);

    if ((*stream)->pid && close_mode == STREAM_KILL_PROCESS) 
#ifdef WIN32
        TerminateProcess ((*stream)->pid, 9); // ignore errors
#else
        kill ((*stream)->pid, 9); // ignore errors
#endif
    
    int exit_code = 0; 
    if ((*stream)->pid && close_mode != STREAM_DONT_WAIT_FOR_PROCESS)
        // TerminateProcess is asynchronous - we need to wait to make sure the process is terminated (not sure about kill)
        exit_code = stream_wait_for_exit (*stream);

    url_reset_if_curl (*stream);
    FREE (*stream);

    return exit_code;
}

// returns exit status of child
int stream_wait_for_exit (Stream *stream)
{
    if (!stream->pid) return stream->exit_status; // we've already waited for this process and already have the exit status;

#ifdef _WIN32
    // wait for child, so that the terminal doesn't print the prompt until the child is done
    WaitForSingleObject (stream->pid, INFINITE);
    
    ASSERTW (GetExitCodeProcess (stream->pid, &stream->exit_status), "Warning: GetExitCodeProcess() failed: %s", stream_windows_error());
    CloseHandle (stream->pid);

#else
    START_TIMER_ALWAYS;

    int exit_status;
    waitpid (stream->pid, &exit_status, 0); 
    
    stream->exit_status = WIFEXITED (exit_status) ? WEXITSTATUS (exit_status) : exit_status;

    ASSERTW (WIFEXITED (exit_status), "Child process pid=%d exited abnormally (i.e. not via exit()). I waited for it %u milliseconds",          
             stream->pid, (unsigned)(CHECK_TIMER / 1000000ULL));

    // in Windows, the main process fails to CreateProcess it exits. In Unix, it is the child process that 
    // fails to execv, and exits and code EXIT_STREAM. The main process catches it here, and exits silently.
    if (stream->exit_status == EXIT_STREAM) 
        exit(EXIT_GENERAL_ERROR); // child process failed to exec and displayed error message, we can exit silently

#endif
    stream->pid = 0;

    return (int)stream->exit_status;
}

FILE *stream_from_stream_stdout (Stream *stream) { return stream->from_stream_stdout; }
FILE *stream_from_stream_stderr (Stream *stream) { return stream->from_stream_stderr; }
FILE *stream_to_stream_stdin    (Stream *stream) { return stream->to_stream_stdin;    }

void stream_abort_if_cannot_run (rom exec_name, rom reason)
{
    StreamP stream = stream_create (0, 1024, 1024, 0, 0, 0, 0, reason, exec_name, NULL); // will abort if cannot run

// in Windows, the main process fails to CreateProcess and exits. In Unix, it is the child process that 
// fails to execv, and exits and code EXIT_STREAM. The main process catches it in stream_wait_for_exit, and exits.
#ifndef _WIN32
    stream_wait_for_exit (stream); // exits if child process exited with code EXIT_STREAM
#endif

    // if we reach here, everything's good - the exec can run.
    stream_close (&stream, STREAM_KILL_PROCESS);
}

// check if executable is in the path. Note: in Windows it actually runs the executable, so 
// only suitable for executables that would terminate immediately.
bool stream_is_exec_in_path (rom exec)
{
#ifdef _WIN32
    // our own instance of curl - to not conflict with url_open
    StreamP curl = stream_create (0, 0, 0, 0, 0, 0, 0,
                                  "where.exe", // reason in case of failure to execute curl
                                  "where.exe", "/Q", exec, NULL); 

    return stream_close (&curl, STREAM_WAIT_FOR_PROCESS) == 0;

#else
    char run[32 + strlen(exec)];
    sprintf (run, "which %s > /dev/null 2>&1", exec);
    return !system ("which wget > /dev/null 2>&1") && file_exists ("/dev/stdout");
#endif
}
