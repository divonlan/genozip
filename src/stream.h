// ------------------------------------------------------------------
//   stream.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include <stdio.h>
#include "genozip.h"

typedef struct stream_ *StreamP;

#define MAX_ARGC 64 // maximum number of arguments including exec_name

#define DEFAULT_PIPE_SIZE 65536
#define SKIP_ARG ((rom )-1)
extern StreamP stream_create (StreamP parent_stream, uint32_t from_stream_stdout, uint32_t from_stream_stderr, uint32_t to_stream, 
                              FILE *redirect_stdout_file, rom input_url_name, bool input_stdin, rom reason,
                              rom exec_name, ...);

extern void stream_set_inheritability (int fd, bool is_inheritable);

extern void stream_close_pipes (StreamP stream);

typedef enum { STREAM_KILL_PROCESS, STREAM_WAIT_FOR_PROCESS, STREAM_DONT_WAIT_FOR_PROCESS } StreamCloseMode;
extern int stream_close (StreamP *stream, StreamCloseMode close_mode);

extern int stream_wait_for_exit (StreamP stream, bool killed);

extern void stream_abort_if_cannot_run (rom exec_name, rom reason);

extern FILE *stream_from_stream_stdout (StreamP stream);
extern FILE *stream_from_stream_stderr (StreamP stream);
extern FILE *stream_to_stream_stdin    (StreamP stream);

extern rom stream_get_exec_name (StreamP stream);

extern rom exit_code_name (int exit_code);

extern bool stream_is_exec_in_path (rom exec);
