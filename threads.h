// ------------------------------------------------------------------
//   threads.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef THREADS_INCLUDED
#define THREADS_INCLUDED

extern void threads_initialize (void);

#define THREADS_NO_VB ((uint32_t)-1)
extern int threads_create (void *(*func)(void *), void *arg, const char *task_name, uint32_t vb_i);
extern void threads_join (int thread_id);
extern void threads_cancel_other_threads (void);
extern bool threads_am_i_main_thread (void);
extern void threads_print_call_stack (void);

#endif
