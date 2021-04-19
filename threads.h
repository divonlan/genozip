// ------------------------------------------------------------------
//   threads.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef THREADS_INCLUDED
#define THREADS_INCLUDED

extern void threads_initialize (void);

extern ThreadId threads_create (void (*func)(VBlockP), VBlockP vb);
extern bool threads_join (ThreadId *thread_id, bool blocking);
extern void threads_cancel_other_threads (void);
extern bool threads_am_i_main_thread (void);
extern void threads_print_call_stack (void);

// for debugging thread issues, activated with --debug-threads or --show-threads
void threads_log_by_vb (ConstVBlockP vb, const char *task_name, const char *event, int time_usec);
void threads_display_log (void);

#define ASSERTMAINTHREAD ASSERT0 (threads_am_i_main_thread(), "expected to be running in main thread")
#endif
