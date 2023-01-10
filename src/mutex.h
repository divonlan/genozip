// ------------------------------------------------------------------
//   mutex.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

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
#include "buffer.h"

// -----------
// mutex stuff
// -----------

typedef struct Mutex {
    rom name, initialized, lock_func;
    pthread_mutex_t mutex;
} Mutex;

extern void mutex_initialize_do (MutexP mutex, rom name, rom func);
#define mutex_initialize(mutex)    mutex_initialize_do (&(mutex), #mutex, __FUNCTION__)

extern void mutex_destroy_do (MutexP mutex, rom func);
#define mutex_destroy(mutex) mutex_destroy_do (&(mutex), __FUNCTION__)

extern bool mutex_lock_do (MutexP mutex, bool blocking, FUNCLINE);
#define mutex_lock(mutex) mutex_lock_do (&(mutex), true, __FUNCLINE)
#define mutex_trylock(mutex) mutex_lock_do (&(mutex), false, __FUNCLINE)

extern void mutex_unlock_do (MutexP mutex, FUNCLINE);
#define mutex_unlock(mutex)    mutex_unlock_do (&(mutex), __FUNCLINE)

extern bool mutex_wait_do (MutexP mutex, bool blocking, FUNCLINE);
#define mutex_wait(mutex, blocking) mutex_wait_do (&(mutex), (blocking), __FUNCLINE)

extern void mutex_show_bottleneck_analsyis (void);

#define mutex_is_show(name) (flag.show_mutex && (flag.show_mutex==(char*)1 || !strncmp ((name), flag.show_mutex, 8))) // only 8 chars so we can catch all genome_muteces[%u]

// -------------------
// VB serializer stuff
// -------------------

typedef struct Serializer {
    Mutex;             // note: gcc/clang flag -fms-extensions is needed for this type of anonymous struct use
    VBIType vb_i_last; // used by serializer_lock
    Buffer skips;      // a bytemap of vb_i's to skip
} Serializer;

extern void serializer_initialize_do (SerializerP ser, rom name, rom func);
#define serializer_initialize(ser) serializer_initialize_do (&(ser), #ser, __FUNCTION__)

extern void serializer_destroy_do (SerializerP ser, rom func);
#define serializer_destroy(ser) serializer_destroy_do (&(ser), __FUNCTION__)

extern void serializer_lock_do (SerializerP ser, VBIType vb_i, FUNCLINE);
#define serializer_lock(ser, vb_i) serializer_lock_do (&(ser), vb_i, __FUNCLINE)

#define serializer_unlock(ser) mutex_unlock_do ((MutexP)&(ser), __FUNCLINE)

extern void serializer_skip_do (SerializerP ser, VBIType vb_i, FUNCLINE);
#define serializer_skip(ser, vb_i) serializer_skip_do (&(ser), vb_i, __FUNCLINE)

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
#define spin_initialize(name) ({ if (! name##_initialized) { \
                                   int ret; ASSERT (!(ret = pthread_spin_init (&name, PTHREAD_PROCESS_PRIVATE)), "failed for %s: %s", #name, strerror (ret)); \
                                   name##_initialized = true;\
                                 } })
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
#define spin_lock(m)   ({ int ret = pthread_spin_lock (&m); \
                          ASSERT (!ret, "pthread_spin_lock failed: %s", strerror (ret)); })

#define spin_unlock(m) ({ int ret = pthread_spin_unlock (&m); \
                          ASSERT (!ret, "pthread_spin_lock failed: %s", strerror (ret)); })
#endif

