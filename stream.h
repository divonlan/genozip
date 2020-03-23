// ------------------------------------------------------------------
//   stream.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef STREAM_INCLUDED
#define STREAM_INCLUDED

#include <stdio.h>
#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "genozip.h"

typedef struct {
    FILE *from_stream_stdout;
    FILE *from_stream_stderr;
    FILE *to_stream_stdin;
#ifdef _WIN32
    HANDLE pid;
#else
    pid_t pid;
#endif
} Stream;

#define MAX_ARGC 30 // maximum number of arguments including exec_name

#define DEFAULT_PIPE_SIZE 65536
#define SKIP_ARG ((const char *)-1)
extern Stream stream_create (uint32_t from_stream, uint32_t from_stream_stderr, uint32_t to_stream, 
                             FILE *redirect_stdout_file, const char *exec_name, ...);

extern void stream_close (Stream *stream);

extern int stream_wait_for_exit (Stream stream);

extern void stream_abort_if_cannot_run (const char *exec_name);

#endif