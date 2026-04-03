// ------------------------------------------------------------------
//   dispatcher.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"
#include "flags.h"

typedef void (*DispatcherFunc)(VBlockP);
typedef enum { JOIN_IN_ORDER, JOIN_OUT_OF_ORDER } DispatcherJoinMode;

extern Dispatcher dispatcher_init (Task task, unsigned max_threads, unsigned previous_vb_i, DispatcherJoinMode join_mode, rom filename, uint64_t target_progress);
extern Dispatcher dispatcher_fan_out_task (Task task, rom filename, uint64_t target_progress, DispatcherJoinMode join_mode, bool force_single_thread, uint32_t previous_vb_i, uint32_t idle_sleep_microsec, DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output);
extern void dispatcher_new_progress_bar (void);
extern void dispatcher_start_wallclock (void);
extern void dispatcher_allow_out_of_order (Dispatcher dispatcher);
extern void dispatcher_pause (Dispatcher dispatcher);
extern void dispatcher_resume (Dispatcher dispatcher, uint32_t target_progress, CompIType comp_i);
extern void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i, bool cleanup_after_me, bool show_memory);
extern void dispatcher_end_task (Dispatcher dispatcher);
extern void dispatcher_compute (Dispatcher dispatcher, DispatcherFunc func);
extern VBlockP dispatcher_generate_next_vb (Dispatcher dispatcher, VBIType vb_i, CompIType comp_i);       
extern bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final);                                  
extern VBlockP dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final, bool blocking);
extern bool dispatcher_has_free_thread (Dispatcher dispatcher);
extern uint32_t dispatcher_get_num_running_compute_threads (Dispatcher dispatcher);
extern uint32_t dispatcher_get_next_vb_i (Dispatcher dispatcher);
extern void dispatcher_recycle_vbs (Dispatcher dispatcher, bool release_vb);
extern void dispatcher_abandon_next_vb (Dispatcher dispatcher);
extern void dispatcher_set_no_data_available (Dispatcher dispatcher, bool abandon_next_vb, DispatchStatus dispatch_status);
extern bool dispatcher_is_input_exhausted (Dispatcher dispatcher);
extern bool dispatcher_is_done (Dispatcher dispatcher);
extern void dispatcher_set_task (Dispatcher dispatcher, Task task);
#define PROGRESS_UNIT (txt_file->est_num_lines ? vb->lines.len : vb->txt_size) // ZIP
extern void dispatcher_increment_progress (rom where, int64_t increment);
extern void dispatcher_calc_avg_compute_vbs (Dispatcher d);

extern rom _task_names[NUM_TASKS];
static inline rom task_name (Task task) { return (task >= 0 && task < NUM_TASKS) ? _task_names[task] : "Invalid_Task"; }
extern Task task_by_name (rom name);

