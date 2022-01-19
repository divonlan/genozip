// ------------------------------------------------------------------
//   threads.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

extern void threads_initialize (void);

extern ThreadId threads_create (void (*func)(VBlockP), VBlockP vb);
extern bool threads_join_do (ThreadId *thread_id, VBlockP vb, const char *func);
#define threads_join(threads_id,vb) threads_join_do((threads_id),(vb),__FUNCTION__)
extern void threads_cancel_other_threads (void);
extern bool threads_am_i_main_thread (void);
extern void threads_print_call_stack (void);

// writer thread stuff
extern void threads_set_writer_thread (void);
extern void threads_unset_writer_thread (void);
extern bool threads_am_i_writer_thread (void);

// for debugging thread issues, activated with --debug-threads or --show-threads
void threads_log_by_vb (ConstVBlockP vb, const char *task_name, const char *event, int time_usec);
void threads_write_log (bool to_info_stream);

#define ASSERTMAINTHREAD ASSERT0 (threads_am_i_main_thread(), "expected to be running in main thread")

