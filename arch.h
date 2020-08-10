// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ARCH_INCLUDED
#define ARCH_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <pthread.h>
#else
#include "compatibility/visual_c_pthread.h"
#endif
#ifdef __APPLE__
#include <os/lock.h>
#endif

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
#define MUTEX(name) \
    static pthread_mutex_t name; \
    static bool name##_initialized = false;

extern void mutex_initialize_do (const char *name, pthread_mutex_t *mutex, bool *initialized);
#define mutex_initialize(name) mutex_initialize_do (#name, &name, &name##_initialized)
#define mutex_destroy(name) if (name##_initialized) { pthread_mutex_destroy (&name); name##_initialized = false; }

#define mutex_lock(m)   { int ret = pthread_mutex_lock (&m); \
                          ASSERT (!ret, "Error in %s: pthread_mutex_lock failed: %s", __FUNCTION__, strerror (ret)); }

#define mutex_unlock(m) { int ret = pthread_mutex_unlock (&m); \
                          ASSERT (!ret, "Error in %s: pthread_mutex_lock failed: %s", __FUNCTION__, strerror (ret)); }

// -----------
// spinlock stuff
// -----------

#ifdef __APPLE__
#define pthread_spinlock_t os_unfair_lock
#define pthread_spin_destroy(x) // do nothing
#endif // __APPLE__

#define SPINLOCK(name) \
    static pthread_spinlock_t name; \
    static bool name##_initialized = false;

extern void spin_initialize_do (const char *name, pthread_spinlock_t *spin, bool *initialized);
#define spin_initialize(name) spin_initialize_do (#name, &name, &name##_initialized)
#define spin_destroy(name) if (name##_initialized) { pthread_spin_destroy (&name); name##_initialized = false; }

#ifdef __APPLE__
#define spin_lock(m)   os_unfair_lock_lock(&m)
#define spin_unlock(m) os_unfair_lock_unlock (&m)
#else
#define spin_lock(m)   { int ret = pthread_spin_lock (&m); \
                         ASSERT (!ret, "Error in %s: pthread_spin_lock failed: %s", __FUNCTION__, strerror (ret)); }

#define spin_unlock(m) { int ret = pthread_spin_unlock (&m); \
                         ASSERT (!ret, "Error in %s: pthread_spin_lock failed: %s", __FUNCTION__, strerror (ret)); }
#endif
#endif
