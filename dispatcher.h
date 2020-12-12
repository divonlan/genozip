// ------------------------------------------------------------------
//   dispatcher.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DISPATCHER_INCLUDED
#define DISPATCHER_INCLUDED

#include "genozip.h"
#include "buffer.h"

typedef void *Dispatcher;

typedef enum { PROGRESS_PERCENT, PROGRESS_MESSAGE, PROGRESS_NONE } ProgressType;
extern Dispatcher dispatcher_init (const char *task_name, unsigned max_threads, unsigned previous_vb_i,
                                   bool test_mode, bool is_last_file, bool cleanup_after_me, const char *filename, ProgressType prog, const char *prog_msg);
extern void dispatcher_pause (Dispatcher dispatcher);
extern void dispatcher_resume (Dispatcher dispatcher);
extern void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i);

typedef void (*DispatcherFunc)(VBlockP);
extern void dispatcher_compute (Dispatcher dispatcher, DispatcherFunc func);
extern VBlockP dispatcher_generate_next_vb (Dispatcher dispatcher, uint32_t vb_i);       
extern bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final);                                  
extern VBlockP dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final);
extern bool dispatcher_has_free_thread (Dispatcher dispatcher);
extern VBlockP dispatcher_get_next_vb (Dispatcher dispatcher);
extern void dispatcher_recycle_vbs (Dispatcher dispatcher);
extern void dispatcher_abandon_next_vb (Dispatcher dispatcher);
extern void dispatcher_set_input_exhausted (Dispatcher dispatcher, bool exhausted);
extern bool dispatcher_is_input_exhausted (Dispatcher dispatcher);
extern bool dispatcher_is_done (Dispatcher dispatcher);
extern void dispatcher_show_time (const char *stage, int32_t thread_index, uint32_t vb_i);
extern uint32_t dispatcher_fan_out_task_do (const char *task_name, const char *filename, ProgressType prog, const char *prog_msg, bool test_mode, bool force_single_thread, DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output);
#define dispatcher_fan_out_task(filename, prog, prog_msg, test_mode, force_single_thread, prepare, compute, output) \
    dispatcher_fan_out_task_do (__FUNCTION__, (filename), (prog), (prog_msg), (test_mode), (force_single_thread), (prepare), (compute), (output))

#endif
