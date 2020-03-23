// ------------------------------------------------------------------
//   stream.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <errno.h>
#include <sys/types.h>
#include <stdarg.h>

#ifdef _WIN32
#include <windows.h>
#include <process.h>
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

#ifdef _WIN32
static const char *stream_windows_error (void)
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
    HANDLE read_handle, write_handle;
    ASSERT (CreatePipe (&read_handle, &write_handle, NULL, pipe_size), 
            "Error: CreatePipe failed: %s", stream_windows_error());
    
    // ensure that "our" side of the pipe isn't inherited by the child process
    ASSERT (SetHandleInformation (is_stream_to_genozip ? read_handle : write_handle, HANDLE_FLAG_INHERIT, 0),
            "Error: SetHandleInformation failed: %s", stream_windows_error());

    fds[0] = _open_osfhandle (PtrToUlong (read_handle), 0); // note: closing the fd with close() closes the HANDLE too: https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/close?view=vs-2019
    ASSERT0 (fds[0] >= 0, "Error: _open_osfhandle failed of fds[0]");

    fds[1] = _open_osfhandle (PtrToUlong (write_handle), 0);
    ASSERT0 (fds[1] >= 0, "Error: _open_osfhandle failed of fds[1]");

#else // Not Windows
    ASSERT (!pipe (fds), "Error: failed to create pipe: %s", strerror (errno));

#ifdef F_SETPIPE_SZ // this is Linux-specific, and available since kernel 2.6.35
    
    // lower pipe size requested to the maximum allowed system size (fcntl will return EPERM if size is larger than limit)
    FILE *pipe_max_size_file = fopen ("/proc/sys/fs/pipe-max-size", "rb");
    if (pipe_max_size_file) {
        char max_pipe_size_str[30];
        memset (max_pipe_size_str, 0, 30);
        size_t bytes = fread (max_pipe_size_str, 1, 29, pipe_max_size_file);
        if (bytes) pipe_size = MIN (pipe_size, atoi (max_pipe_size_str));
        FCLOSE (pipe_max_size_file, "pipe_max_size_file");
    }

    // set the pipe size
    fcntl(fds[0], F_SETPIPE_SZ, pipe_size); // ignore errors
#endif // F_SETPIPE_SZ

#endif // not Windows
}

static void stream_abort_cannot_exec (const char *exec_name)
{
    fprintf (stderr, "\n%s failed because it could not execute %s.\n"
                     "%s needs to be installed and in the execution path.\n"
                     "%s is a separate software package that is not affiliated with genozip in any way.\n"
                     "Please see 'genozip --help' for more details\n\n",
             global_cmd, exec_name, exec_name, exec_name);  

    exit (99);
}

#ifdef _WIN32
static HANDLE stream_exec_child (int *stream_stdout_to_genozip, int *stream_stderr_to_genozip, int *genozip_to_stream_stdin,
                                 FILE *redirect_stdout_file, unsigned argc, char * const *argv)
{
    int cmd_line_len = 0;
    for (int i=0; i < argc; i++)
        if (argv[i]) cmd_line_len += strlen (argv[i]) + 1; // +1 for the space (or \0) after

    char *cmd_line = malloc (cmd_line_len);
    char *fmt_tmpl = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"; // # of %s = MAX_ARGC
    char fmt[MAX_ARGC * 3];
    strcpy (fmt, fmt_tmpl);
    fmt[argc*3 - 1] = '\0'; // cut short the fmt string
    sprintf (cmd_line, fmt, argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12], argv[13], argv[14], argv[15], argv[16], argv[17], argv[18], argv[19], argv[20], argv[21], argv[22], argv[23], argv[24], argv[25], argv[26], argv[27], argv[28], argv[29]);

    STARTUPINFO startup_info;
    memset (&startup_info, 0, sizeof startup_info);
    startup_info.cb         = sizeof startup_info;
    startup_info.hStdError  = stream_stderr_to_genozip ? (HANDLE)_get_osfhandle (stream_stderr_to_genozip[1]) : GetStdHandle (STD_ERROR_HANDLE);
    startup_info.hStdInput  = genozip_to_stream_stdin  ? (HANDLE)_get_osfhandle (genozip_to_stream_stdin[0])  : INVALID_HANDLE_VALUE;
    startup_info.dwFlags    = STARTF_USESTDHANDLES;

    // determine where the child's stdout will go 
    if (stream_stdout_to_genozip)  startup_info.hStdOutput = (HANDLE)_get_osfhandle (stream_stdout_to_genozip[1]);
    else if (redirect_stdout_file) startup_info.hStdOutput = (HANDLE)_get_osfhandle (fileno (redirect_stdout_file));
    else                           startup_info.hStdOutput = GetStdHandle (STD_OUTPUT_HANDLE);

    PROCESS_INFORMATION proc_info;
    memset (&proc_info, 0, sizeof proc_info);

    bool success = CreateProcess (NULL, cmd_line, NULL, NULL, TRUE, NORMAL_PRIORITY_CLASS, 
                                  NULL, NULL, &startup_info, &proc_info);
    if (!success) stream_abort_cannot_exec (argv[0]);
    
    // close handles we (the parent process) don't need
    if (stream_stdout_to_genozip) CLOSE (stream_stdout_to_genozip[1], "stream_stdout_to_genozip[1]"); // this also closes the Windows HANDLE
    if (stream_stderr_to_genozip) CLOSE (stream_stderr_to_genozip[1], "stream_stderr_to_genozip[1]");
    if (genozip_to_stream_stdin)  CLOSE (genozip_to_stream_stdin[0],  "genozip_to_stream_stdin[0]");
    CloseHandle (proc_info.hThread); // child process's main thread

    return proc_info.hProcess;
}

#else // not Windows
static pid_t stream_exec_child (int *stream_stdout_to_genozip, int *stream_stderr_to_genozip, int *genozip_to_stream_stdin,
                                FILE *redirect_stdout_file, unsigned argc, char * const *argv)
{
    pid_t child_pid = fork();
    if (child_pid) return child_pid; // parent returns

    // redirect child stdin to our pipe if needed
    if (genozip_to_stream_stdin) {         
        dup2 (genozip_to_stream_stdin[0], STDIN_FILENO);
        CLOSE (genozip_to_stream_stdin[1], "genozip_to_stream_stdin[1]"); // child closes writing side of the pipe  
    }

    // redirect child stdout to our pipe if needed
    if (redirect_stdout_file)          
        dup2 (fileno (redirect_stdout_file), STDOUT_FILENO);
 
    // or - redirect child stdout to a file if requested
    else if (stream_stdout_to_genozip) {         
        dup2 (stream_stdout_to_genozip[1], STDOUT_FILENO);
        CLOSE (stream_stdout_to_genozip[0], "stream_stdout_to_genozip[0]"); // child closes reading side of the pipe  
    }

    // redirect child stderr to our pipe if needed
    if (stream_stderr_to_genozip) {
        dup2 (STDERR_FILENO, 3); // keep the original stderr in fd=3 in case we can't execute and need to report an error
        dup2 (stream_stderr_to_genozip[1], STDERR_FILENO);
        CLOSE (stream_stderr_to_genozip[0], "stream_stderr_to_genozip[0]"); // child closes reading side of the pipe
    }

    execvp (argv[0], (char * const *)argv); // if successful, doesn't return

    if (stream_stderr_to_genozip) dup2 (3, STDERR_FILENO); // get back the original stderr

    stream_abort_cannot_exec (argv[0]);
    
    return 0; // squash compiler warning - never reach here
}

#endif

Stream stream_create (uint32_t from_stream_stdout, uint32_t from_stream_stderr, uint32_t to_stream_stdin,
                      FILE *redirect_stdout_file, const char *exec_name, ...)
{
    ASSERT0 (!from_stream_stdout || !redirect_stdout_file, "Error in stream_create: cannot redirect child output to both genozip and a file");

    int stream_stdout_to_genozip[2], stream_stderr_to_genozip[2], genozip_to_stream_stdin[2];
    
    if (from_stream_stdout) stream_pipe (stream_stdout_to_genozip, from_stream_stdout, true);
    if (from_stream_stderr) stream_pipe (stream_stderr_to_genozip, from_stream_stderr, true);
    if (to_stream_stdin)    stream_pipe (genozip_to_stream_stdin,  to_stream_stdin, false);

    Stream stream;
    memset (&stream, 0, sizeof (Stream));

    // copy our function arguments to argv
    va_list argp;
    va_start (argp, exec_name);

    const char *argv[MAX_ARGC];
    memset (argv, 0, sizeof(argv));
    unsigned argc = 0;

    argv[argc++] = exec_name;

    const char *arg;
    while (argc < MAX_ARGC) {
        arg = va_arg (argp, const char *);
        if (arg == SKIP_ARG) continue;
        if (arg == NULL) break;
        argv[argc++] = arg;
    }

    ASSERT (!arg, "Error: too many arguments - max is %u", MAX_ARGC-1); // MAX_ARGC-1 real args and last must be NULL

    // execute the child in a separate process
    stream.pid = stream_exec_child (from_stream_stdout ? stream_stdout_to_genozip : NULL, 
                                    from_stream_stderr ? stream_stderr_to_genozip : NULL, 
                                    to_stream_stdin    ? genozip_to_stream_stdin  : NULL, 
                                    redirect_stdout_file,
                                    argc, (char * const *)argv); 

    if (redirect_stdout_file) 
        FCLOSE (redirect_stdout_file, "redirect_stdout_file"); // the child has this file open, we don't need it

    // store our side of the pipe, and close the side used by the child
    if (from_stream_stdout) {
        stream.from_stream_stdout = fdopen (stream_stdout_to_genozip[0], READ);
        CLOSE (stream_stdout_to_genozip[1], "stream_stdout_to_genozip[1]");
    }

    if (from_stream_stderr) {
        stream.from_stream_stderr = fdopen (stream_stderr_to_genozip[0], READ);
        CLOSE (stream_stderr_to_genozip[1], "stream_stderr_to_genozip[1]");
    }

    if (to_stream_stdin) {
        stream.to_stream_stdin = fdopen (genozip_to_stream_stdin[1], WRITE);
        CLOSE (genozip_to_stream_stdin[0], "genozip_to_stream_stdin[0]");
    }

    return stream;
}

void stream_close (Stream *stream)
{
#ifdef WIN32
    TerminateProcess (stream->pid, 9); // ignore errors
#else
    kill (stream->pid, 9); // ignore errors
#endif

    if (stream->from_stream_stdout) FCLOSE (stream->from_stream_stdout, "stream->from_stream_stdout");
    if (stream->from_stream_stderr) FCLOSE (stream->from_stream_stderr, "stream->from_stream_stderr");
    if (stream->to_stream_stdin)    FCLOSE (stream->to_stream_stdin,    "stream->to_stream_stdin");

    memset (stream, 0, sizeof (Stream));
}

// returns exit status of child
int stream_wait_for_exit (Stream stream)
{
#ifdef _WIN32
    // wait for child, so that the terminal doesn't print the prompt until the child is done
    WaitForSingleObject (stream.pid, INFINITE);
    
    DWORD exit_status = 0;
    ASSERTW (GetExitCodeProcess (stream.pid, &exit_status), "Error: GetExitCodeProcess() failed: %s", stream_windows_error());
    CloseHandle (stream.pid);

#else
    int exit_status = 0;
    waitpid (stream.pid, &exit_status, 0); 

    // in Windows, the main process fails to CreateProcess it exits. In Unix, it is the child process that 
    // fails to execv, and exits and code 99. The main process catches it here, and exits silently.
    if (WEXITSTATUS (exit_status) == 99) exit(1); // child process failed to exec and displayed error message, we can exit silently

#endif

    return (int)exit_status;
}

void stream_abort_if_cannot_run (const char *exec_name)
{
    Stream stream = stream_create (1024, 1024, 0, 0, exec_name, NULL); // will abort if cannot run

// in Windows, the main process fails to CreateProcess and exits. In Unix, it is the child process that 
// fails to execv, and exits and code 99. The main process catches it in stream_wait_for_exit, and exits.
#ifndef _WIN32
    stream_wait_for_exit (stream); // exits if child process exited with code 99
#endif

    // if we reach here, everything's good - the exec can run.
    stream_close (&stream);
}
