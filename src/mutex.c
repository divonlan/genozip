// ------------------------------------------------------------------
//   mutex.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include "mutex.h"
#include "flags.h"
#include "buffer.h"
#include "profiler.h"
#include "sorter.h"

typedef struct { 
    int64_t accumulator;
    rom mutex_name;
    Caller caller; 
    uint32_t code_line; 
    uint32_t lock_count; // multiple muteces (e.g. same mutex in different contexts of VBs) can be locked concurrently
} LockPoint;

#define MAX_CODE_LINE 4095
static LockPoint lp[MAX_CODE_LINE+1]; // note: a static array, becauses its hard to use a Buffer, because it uses muteces...

void mutex_initialize_do (Mutex *mutex, rom name, rom func)
{ 
    if (!mutex->initialized) {
        int ret = pthread_mutex_init (&mutex->mutex, 0); 
        ASSERT (!ret || errno == EBUSY,  // EBUSY is not an error - the failure is bc a race condition and the mutex is already initialized - all good
                "pthread_mutex_init failed for %s from %s: %s", name, func, strerror (ret)); 
    }

    mutex->name = (uintptr_t)name;
    mutex->initialized = true;
}

void mutex_destroy_do (Mutex *mutex, rom func) 
{
    if (!mutex->initialized) return;
        
    pthread_mutex_destroy (&mutex->mutex); 
    memset (mutex, 0, sizeof (Mutex));
}

bool mutex_lock_do (Mutex *mutex, bool blocking, Caller caller)   
{ 
    ASSERT (mutex->initialized, "called from %s:%u: mutex not initialized", CALLERf);

    bool show = mutex_is_show ((rom)(uintptr_t)mutex->name);

    if (show) iprintf ("LOCKING : Mutex %s by thread %"PRIu64" %s:%u\n", (rom)(uintptr_t)mutex->name, (uint64_t)pthread_self(), CALLERf);

    int ret;
    if (blocking) {
        𝓅𝓇ℴ𝒻𝒾𝓁ℯ (ProfilerTime start_time = (flag.show_time_comp_i == COMP_ALL) ? get_timer_start() : NULL_TIMER;)

        ret = pthread_mutex_lock (&mutex->mutex);

        if (__builtin_expect (flag.show_time_comp_i != COMP_NONE, false)) { // test same condition as START_TIMER 
            if (!lp[caller.code_line].mutex_name) { // first lock at this lockpoint
                ASSERT (caller.code_line <= MAX_CODE_LINE, "mutex_lock at %s:%u: cannot lock a mutex in a code_line > %u", CALLERf, MAX_CODE_LINE);
                lp[caller.code_line] = (LockPoint){ .mutex_name = (rom)(uintptr_t)mutex->name, .caller = caller };
            }

            else 
                if (lp[caller.code_line].caller.funcר != caller.funcר) 
                    WARN_ONCE (_FYI "Two calls to mutex_lock exist on the same code_line: %s @ %s:%u and %s @ %s:%u - --show-time will show their combined time. To solve, add an empty line to shift the code line number of one of them",
                               lp[caller.code_line].mutex_name, CALLERff(lp[caller.code_line].caller), (rom)(uintptr_t)mutex->name, CALLERf);
                
            if (flag.show_time_comp_i == COMP_ALL)
                𝓅𝓇ℴ𝒻𝒾𝓁ℯ (lp[caller.code_line].accumulator += get_timer_delta (start_time)); // luckily, we're protected by the mutex...
        }
    }
    
    else {
        ret = pthread_mutex_trylock (&mutex->mutex);
        if (ret == EBUSY) return false;
    }

    increment_relaxed (lp[caller.code_line].lock_count);

    ASSERT (!ret, "called from %s:%u by thread=%"PRIu64": pthread_mutex_lock failed on mutex->name=%s: %s", 
            CALLERf, (uint64_t)pthread_self(), mutex && mutex->name ? (rom)(uintptr_t)mutex->name : "(null)", strerror (ret)); 

    mutex->locked = true; // mutex->locked is protected by the mutex

    if (show) iprintf ("LOCKED  : Mutex %s by thread %"PRIu64"\n", (rom)(uintptr_t)mutex->name, (uint64_t)pthread_self());

    return true;
}

void mutex_unlock_do (Mutex *mutex, Caller caller) 
{ 
    ASSERT (mutex->initialized, "called from %s:%u mutex not initialized", CALLERf);
    ASSERT (mutex->locked, "called from %s:%u by thread=%"PRIu64": mutex %s is not locked", 
            CALLERf, (uint64_t)pthread_self(), (rom)(uintptr_t)mutex->name);

    mutex->locked = false; // mutex->locked is protected by the mutex

    decrement_relaxed (lp[caller.code_line].lock_count);

    int ret = pthread_mutex_unlock (&mutex->mutex); 
    ASSERT (!ret, "called from %s:%u: pthread_mutex_unlock failed for %s: %s", CALLERf, mutex->name, strerror (ret)); 

    if (mutex_is_show ((rom)(uintptr_t)mutex->name))
        iprintf ("UNLOCKED: Mutex %s by thread %"PRIu64" %s:%u\n", (rom)(uintptr_t)mutex->name, (uint64_t)pthread_self(), CALLERf);
}

bool mutex_wait_do (Mutex *mutex, bool blocking, Caller caller)   
{
    if (mutex_lock_do (mutex, blocking, caller)) {
        mutex_unlock_do (mutex, caller);
        return true;
    }
    
    else
        return false; // didn't lock (can only happen if non-blocking)
}

void serializer_initialize_do (SerializerP ser, rom name, rom func)
{
    ASSERT (!ser->mutex.initialized, "called from %s: serializer already initialized", func);
    mutex_initialize_do (&ser->mutex, name, func);
}

void serializer_destroy_do (SerializerP ser, rom func)
{
    if (!ser->mutex.initialized) return; // nothing to do

    mutex_destroy_do (&ser->mutex, func);
}

void serializer_lock_do (SerializerP ser, VBIType vb_i, Caller caller)
{
    #define WAIT_TIME_USEC 5000
    #define TIMEOUT (30*60) // 30 min

    for (unsigned i=0; ; i++) {
        mutex_lock_do (&ser->mutex, true, caller);

        ASSERT (ser->vb_i_last < vb_i, "called from %s:%u: Expecting vb_i_last=%u < vb->vblock_i=%u. serializer=%s", 
                CALLERf, ser->vb_i_last, vb_i, ser->mutex.name);
        
        if (ser->vb_i_last == vb_i - 1) { // its our turn now
            ser->vb_i_last++; // next please
            return;           // return with mutex locked
        }
        
        // not our turn, wait 5ms and try again
        mutex_unlock_do (&ser->mutex, caller);
        usleep (WAIT_TIME_USEC);

        // timeout after approx 30 minutes
        ASSERT (i < TIMEOUT * (1000000 / WAIT_TIME_USEC), "called from %s:%u: Timeout (%u sec) while waiting for serializer %s in vb=%u. vb_i_last=%u", 
                CALLERf, TIMEOUT, ser->mutex.name, vb_i, ser->vb_i_last);
    }
}

static DESCENDING_SORTER (mutex_sort_by_accumulator, LockPoint, accumulator)

void mutex_bottleneck_analysis_init (void)
{
    memset (lp, 0, sizeof(lp));
}

void mutex_show_bottleneck_analsyis (void)
{
    qsort (lp, MAX_CODE_LINE+1, sizeof(LockPoint), mutex_sort_by_accumulator);

    iprint0 ("Bottleneck analysis - Time waiting on locks:\n"
             "Millisec  Mutex / Join            LockPoint\n");

    for (int i=0; i <= MAX_CODE_LINE; i++) {
        if (!lp[i].accumulator) break; // done, since its sorted

        iprintf ("%-9s %-23s %s:%u\n", str_int_commas (lp[i].accumulator / 1000000).s, 
                 lp[i].mutex_name, CALLERff(lp[i].caller));
    }
}

// this is called from Ctrl-C. Works only if --show-time is used as well.
void mutex_who_is_locked (void)
{
    for (int i=0; i <= MAX_CODE_LINE; i++) {
        LockPoint my_lp = lp[i]; // make a copy for a bit of thread safety
        if (my_lp.mutex_name && my_lp.lock_count)
            printf ("Mutex locked: %s locked in %s:%u. %s\n", 
                    my_lp.mutex_name, CALLERff(my_lp.caller), 
                    cond_int (my_lp.lock_count > 1, "num_locks_from_different_objects=", my_lp.lock_count));
    }
}

// call so join time will reported in by profiler (in --show-time)
void thread_join_lock_point (rom thread_name, ProfilerTime start_time, Caller caller)
{
    if (!lp[caller.code_line].mutex_name) { // first lock at this lockpoint
        ASSERT (caller.code_line <= MAX_CODE_LINE, "pthreads_join at %s:%u: cannot lock a mutex in a code_line > %u", CALLERf, MAX_CODE_LINE);
        lp[caller.code_line] = (LockPoint){ .mutex_name = thread_name, .caller = caller };
    }

    else 
        if (lp[caller.code_line].caller.funcר != caller.funcר) 
            WARN_ONCE (_FYI "Two calls to mutex_lock/pthreads_join exist on the same code_line: %s @ %s:%u and %s @ %s:%u - --show-time will show their combined time. To solve, add an empty line to shift the code line number of one of them",
                        lp[caller.code_line].mutex_name, CALLERff(lp[caller.code_line].caller), thread_name, CALLERf);

    add_relaxed (lp[caller.code_line].accumulator, get_timer_delta (start_time)); 
}
