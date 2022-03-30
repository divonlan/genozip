// ------------------------------------------------------------------
//   mutex.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include <pthread.h>
#ifdef __APPLE__
#include <availability.h>
#ifdef __MAC_10_12
#include <os/lock.h>
#else // spinlock was available prior to 10.12
#include <libkern/OSAtomic.h>
#endif
#endif
#include "genozip.h"

// -----------
// mutex stuff
// -----------

typedef struct Mutex {
    rom name, initialized, lock_func;
    pthread_mutex_t mutex;
    VBIType vb_i_last; // used by mutex_lock_by_vb_order
} Mutex;

extern void mutex_initialize_do (MutexP mutex, rom name, rom func);
#define mutex_initialize(mutex) mutex_initialize_do (&(mutex), #mutex, __FUNCTION__)

extern void mutex_destroy_do (MutexP mutex, rom func);
#define mutex_destroy(mutex) mutex_destroy_do (&(mutex), __FUNCTION__)

extern bool mutex_lock_do (MutexP mutex, bool blocking, rom func);
#define mutex_lock(mutex) mutex_lock_do (&(mutex), true, __FUNCTION__)
#define mutex_trylock(mutex) mutex_lock_do (&(mutex), false, __FUNCTION__)

extern void mutex_unlock_do (MutexP mutex, FUNCLINE);
#define mutex_unlock(mutex) mutex_unlock_do (&(mutex), __FUNCLINE)

extern void mutex_wait_do (MutexP mutex, FUNCLINE);
#define mutex_wait(mutex) mutex_wait_do (&(mutex), __FUNCLINE)

extern void mutex_lock_by_vb_order_do (VBIType vb_i, MutexP mutex, FUNCLINE);
#define mutex_lock_by_vb_order(vb_i, mutex) mutex_lock_by_vb_order_do (vb_i, &(mutex), __FUNCLINE)

#define mutex_is_show(name) (flag.show_mutex && (flag.show_mutex==(char*)1 || !strncmp ((name), flag.show_mutex, 8))) // only 8 chars so we can catch all genome_muteces[%u]

// --------------
// spinlock stuff
// --------------

#ifdef __APPLE__
#ifdef __MAC_10_12
#define pthread_spinlock_t os_unfair_lock
#define spin_lock(m)   os_unfair_lock_lock(&m)
#define spin_unlock(m) os_unfair_lock_unlock(&m)
#define spin_initialize(name) { if (! name##_initialized) { name = OS_UNFAIR_LOCK_INIT; name##_initialized = true; } }
#define spin_destroy(name) name##_initialized = false; 

#else // spinlock was available prior to 10.12
#define pthread_spinlock_t OSSpinLock
#define spin_initialize(name) { if (! name##_initialized) { name = 0; name##_initialized = true; } }
#define spin_destroy(name) name##_initialized = false; 
#define spin_lock(m)   OSSpinLockLock(&m)
#define spin_unlock(m) OSSpinLockUnlock(&m)
#endif // __MAC_10_12

#else // not mac
#define spin_initialize(name) do { if (! name##_initialized) { \
                                   int ret; ASSERT (!(ret = pthread_spin_init (&name, PTHREAD_PROCESS_PRIVATE)), "failed for %s: %s", #name, strerror (ret)); \
                                   name##_initialized = true;  \
                                 } } while (0)
#define spin_destroy(name) if (name##_initialized) { pthread_spin_destroy (&name); name##_initialized = false; }

#endif // __APPLE__

#define SPINLOCK(name) \
    pthread_spinlock_t name; \
    bool name##_initialized

#ifdef __APPLE__
#ifdef __MAC_10_12
#else // spinlock was available prior to 10.12
#endif // __MAC_10_12

#else
#define spin_lock(m)   do { int ret = pthread_spin_lock (&m); \
                            ASSERT (!ret, "pthread_spin_lock failed: %s", strerror (ret)); } while(0)

#define spin_unlock(m) do { int ret = pthread_spin_unlock (&m); \
                            ASSERT (!ret, "pthread_spin_lock failed: %s", strerror (ret)); } while(0)
#endif

