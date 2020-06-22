// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ARCH_INCLUDED
#define ARCH_INCLUDED

#include <pthread.h>

extern void arch_initialize (void);

extern unsigned arch_get_num_cores (void);
extern const char *arch_get_ip_addr (void);
extern const char *arch_get_os (void);
extern bool arch_am_i_in_docker (void);
extern bool arch_am_i_io_thread (void);
extern void cancel_io_thread (void);

// -----------
// mutex stuff
// -----------
extern void mutex_initialize (pthread_mutex_t *mutex, bool *initialized);

#define mutex_lock(m)   { int ret = pthread_mutex_lock (&m); \
                          ASSERT (!ret, "Error in %s: pthread_mutex_lock failed: %s", __FUNCTION__, strerror (ret)); }

#define mutex_unlock(m) { int ret = pthread_mutex_unlock (&m); \
                          ASSERT (!ret, "Error in %s: pthread_mutex_lock failed: %s", __FUNCTION__, strerror (ret)); }

#endif
