// ------------------------------------------------------------------
//   threads.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

extern void threads_initialize (void);
extern void threads_finalize (void);

extern ThreadId threads_create (void (*func)(VBlockP), VBlockP vb);
extern void threads_join_do (ThreadId *thread_id, rom expected_task, rom expected_task2, rom func);
#define threads_join(threads_id, expected_task) threads_join_do((threads_id), (expected_task), (expected_task), __FUNCTION__)
#define threads_join2(threads_id, expected_task, expected_task2) threads_join_do((threads_id), (expected_task), (expected_task2), __FUNCTION__)
extern void threads_cancel_other_threads (void);
extern void threads_print_call_stack (void);

extern pthread_t main_thread;
static inline bool threads_am_i_main_thread (void) { return pthread_self() == main_thread; }

// writer / zriter thread stuff
extern void threads_set_writer_thread (void);
extern void threads_set_zriter_thread (void);
extern void threads_unset_writer_thread (void);
extern void threads_unset_zriter_thread (void);
extern bool threads_am_i_writer_thread (void);
extern bool threads_am_i_zriter_thread (void);

// for debugging thread issues, activated with --debug-threads or --show-threads
void threads_log_by_vb (ConstVBlockP vb, rom task_name, rom event, int time_usec);
void threads_write_log (bool to_info_stream);

#define ASSERTMAINTHREAD ASSERT (threads_am_i_main_thread(), "%s can only be called in main thread", __FUNCTION__)

