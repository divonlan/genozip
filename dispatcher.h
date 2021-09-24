// ------------------------------------------------------------------
//   dispatcher.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef DISPATCHER_INCLUDED
#define DISPATCHER_INCLUDED

#include "genozip.h"
#include "buffer.h"

typedef void *Dispatcher;

typedef enum { PROGRESS_PERCENT, PROGRESS_MESSAGE, PROGRESS_NONE } ProgressType;
extern Dispatcher dispatcher_init (const char *task_name, unsigned max_threads, unsigned previous_vb_i,
                                   bool test_mode, bool is_last_file, bool cleanup_after_me, const char *filename, ProgressType prog, const char *prog_msg);
extern void dispatcher_start_wallclock (void);
extern void dispatcher_pause (Dispatcher dispatcher);
extern void dispatcher_resume (Dispatcher dispatcher);
extern void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i);

typedef void (*DispatcherFunc)(VBlockP);
extern void dispatcher_compute (Dispatcher dispatcher, DispatcherFunc func);
extern VBlockP dispatcher_generate_next_vb (Dispatcher dispatcher, uint32_t vb_i);       
extern bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final);                                  
extern VBlockP dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final, bool blocking);
extern bool dispatcher_has_free_thread (Dispatcher dispatcher);
extern bool dispatcher_has_active_threads (Dispatcher dispatcher);
extern VBlockP dispatcher_get_next_vb (Dispatcher dispatcher);
extern uint32_t dispatcher_get_next_vb_i (Dispatcher dispatcher);
extern void dispatcher_recycle_vbs (Dispatcher dispatcher, bool release_vb);
extern void dispatcher_abandon_next_vb (Dispatcher dispatcher);
extern void dispatcher_set_input_exhausted (Dispatcher dispatcher, bool exhausted);
extern bool dispatcher_is_input_exhausted (Dispatcher dispatcher);
extern bool dispatcher_is_done (Dispatcher dispatcher);
extern Dispatcher dispatcher_fan_out_task (const char *task_name, const char *filename, ProgressType prog, const char *prog_msg, bool test_mode, bool is_last_file, bool cleanup_after_me, bool force_single_thread, uint32_t previous_vb_i, uint32_t idle_sleep_microsec, DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output);

#endif
