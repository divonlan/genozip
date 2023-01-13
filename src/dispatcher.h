// ------------------------------------------------------------------
//   dispatcher.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"

extern Dispatcher dispatcher_init (rom task_name, VBlockPoolType pool_type, unsigned max_threads, unsigned previous_vb_i,
                                   bool out_of_order, bool test_mode, rom filename, uint64_t target_progress, rom prog_msg);
extern void dispatcher_start_wallclock (void);
extern void dispatcher_pause (Dispatcher dispatcher);
extern void dispatcher_resume (Dispatcher dispatcher);
extern void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i, bool cleanup_after_me, bool show_memory);

typedef void (*DispatcherFunc)(VBlockP);
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
extern void dispatcher_set_task_name (Dispatcher dispatcher, rom task_name);
extern Dispatcher dispatcher_fan_out_task (rom task_name, rom filename, uint32_t target_progress, rom prog_msg, bool out_of_order, bool test_mode, bool force_single_thread, uint32_t previous_vb_i, uint32_t idle_sleep_microsec, DispatcherFunc prepare, DispatcherFunc compute, DispatcherFunc output);
extern void dispatcher_increment_progress (rom where, int64_t increment);
