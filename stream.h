// ------------------------------------------------------------------
//   stream.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef STREAM_INCLUDED
#define STREAM_INCLUDED

#include <stdio.h>
#include "genozip.h"

typedef struct stream_ *StreamP;

#define MAX_ARGC 30 // maximum number of arguments including exec_name

#define DEFAULT_PIPE_SIZE 65536
#define SKIP_ARG ((const char *)-1)
extern StreamP stream_create (StreamP parent_stream, uint32_t from_stream, uint32_t from_stream_stderr, uint32_t to_stream, 
                              FILE *redirect_stdout_file, const char *input_url_name,
                              const char *exec_name, ...);

extern void stream_close_pipes (StreamP stream);

extern void stream_close (StreamP *stream);

extern int stream_wait_for_exit (StreamP stream);

extern void stream_abort_if_cannot_run (const char *exec_name);

extern FILE *stream_from_stream_stdout (StreamP stream);
extern FILE *stream_from_stream_stderr (StreamP stream);
extern FILE *stream_to_stream_stdin    (StreamP stream);

#endif