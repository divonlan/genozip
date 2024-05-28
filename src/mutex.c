// ------------------------------------------------------------------
//   mutex.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <errno.h>
#include "mutex.h"
#include "flags.h"
#include "buffer.h"
#include "profiler.h"

typedef struct { 
    rom mutex_name;  
    rom func; 
    int64_t accumulator;
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

    mutex->name        = name;
    mutex->initialized = func;
}

void mutex_destroy_do (Mutex *mutex, rom func) 
{
    if (!mutex->initialized) return;
        
    pthread_mutex_destroy (&mutex->mutex); 
    memset (mutex, 0, sizeof (Mutex));
}

bool mutex_lock_do (Mutex *mutex, bool blocking, rom func, uint32_t code_line)   
{ 
    ASSERT (mutex->initialized, "called from %s: mutex not initialized", func);

    bool show = mutex_is_show (mutex->name);

    if (show) iprintf ("LOCKING : Mutex %s by thread %"PRIu64" %s\n", mutex->name, (uint64_t)pthread_self(), func);

    int ret;
    if (blocking) {
        START_TIMER;
        ret = pthread_mutex_lock (&mutex->mutex);

        if (flag.show_time_comp_i != COMP_NONE) { // test same condition as START_TIMER 
            if (!lp[code_line].mutex_name) { // first lock at this lockpoint
                ASSERT (code_line <= MAX_CODE_LINE, "mutex_lock at %s:%u: cannot lock a mutex in a code_line > %u", func, code_line, MAX_CODE_LINE);
                lp[code_line] = (LockPoint){ .mutex_name = mutex->name, .func = func, .code_line = code_line };
            }

            else 
                if (lp[code_line].func != func) 
                    WARN_ONCE ("FYI: Two calls to mutex_lock exist on the same code_line: %s @ %s:%u and %s @ %s:%u - --show-time will show their combined time. To solve, add an empty line to shift the code line number of one of them",
                               lp[code_line].mutex_name, lp[code_line].func, lp[code_line].code_line, mutex->name, func, code_line);
                
            lp[code_line].accumulator += CHECK_TIMER; // luckily, we're protected by the mutex...
        }
    }
    
    else {
        ret = pthread_mutex_trylock (&mutex->mutex);
        if (ret == EBUSY) return false;
    }

    __atomic_add_fetch (&lp[code_line].lock_count, (uint32_t)1, __ATOMIC_RELAXED);

    ASSERT (!ret, "called from %s by thread=%"PRIu64": pthread_mutex_lock failed on mutex->name=%s: %s", 
            func, (uint64_t)pthread_self(), mutex && mutex->name ? mutex->name : "(null)", strerror (ret)); 

    mutex->lock_func = func; // mutex->lock_func is protected by the mutex

    if (show) iprintf ("LOCKED  : Mutex %s by thread %"PRIu64"\n", mutex->name, (uint64_t)pthread_self());

    return true;
}

void mutex_unlock_do (Mutex *mutex, FUNCLINE) 
{ 
    ASSERT (mutex->initialized, "called from %s:%u mutex not initialized", func, code_line);
    ASSERT (mutex->lock_func, "called from %s:%u by thread=%"PRIu64": mutex %s is not locked", 
            func, code_line, (uint64_t)pthread_self(), mutex->name);

    mutex->lock_func = NULL; // mutex->lock_func is protected by the mutex

    __atomic_sub_fetch (&lp[code_line].lock_count, (uint32_t)1, __ATOMIC_RELAXED);

    int ret = pthread_mutex_unlock (&mutex->mutex); 
    ASSERT (!ret, "called from %s:%u: pthread_mutex_unlock failed for %s: %s", func, code_line, mutex->name, strerror (ret)); 

    if (mutex_is_show (mutex->name))
        iprintf ("UNLOCKED: Mutex %s by thread %"PRIu64" %s\n", mutex->name, (uint64_t)pthread_self(), func);
}

bool mutex_wait_do (Mutex *mutex, bool blocking, FUNCLINE)   
{
    if (mutex_lock_do (mutex, blocking, func, code_line)) {
        mutex_unlock_do (mutex, func, code_line);
        return true;
    }
    
    else
        return false; // didn't lock (can only happen if non-blocking)
}

void serializer_initialize_do (SerializerP ser, rom name, rom func)
{
    ASSERT (!ser->initialized, "called from %s: serializer already initialized", func);
    mutex_initialize_do ((MutexP)ser, name, func);
}

void serializer_destroy_do (SerializerP ser, rom func)
{
    if (!ser->initialized) return; // nothing to do

    mutex_destroy_do ((MutexP)ser, func);
}

void serializer_lock_do (SerializerP ser, VBIType vb_i, FUNCLINE)
{
    #define WAIT_TIME_USEC 5000
    #define TIMEOUT (30*60) // 30 min

    for (unsigned i=0; ; i++) {
        mutex_lock_do ((MutexP)ser, true, func, code_line);

        ASSERT (ser->vb_i_last < vb_i, "called from %s:%u: Expecting vb_i_last=%u < vb->vblock_i=%u. serializer=%s", 
                func, code_line, ser->vb_i_last, vb_i, ser->name);
        
        if (ser->vb_i_last == vb_i - 1) { // its our turn now
            ser->vb_i_last++; // next please
            return;           // return with mutex locked
        }
        
        // not our turn, wait 5ms and try again
        mutex_unlock_do ((MutexP)ser, func, code_line);
        usleep (WAIT_TIME_USEC);

        // timeout after approx 30 minutes
        ASSERT (i < TIMEOUT * (1000000 / WAIT_TIME_USEC), "called from %s:%u: Timeout (%u sec) while waiting for serializer %s in vb=%u. vb_i_last=%u", 
                func, code_line, TIMEOUT, ser->name, vb_i, ser->vb_i_last);
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
             "Millisec Mutex                   LockPoint\n");

    for (int i=0; i <= MAX_CODE_LINE; i++) {
        if (!lp[i].accumulator) break; // done, since its sorted

        iprintf ("%-8s %-23s %s:%u\n", str_int_commas (lp[i].accumulator / 1000000).s, lp[i].mutex_name, lp[i].func, lp[i].code_line);
    }
}

// this is called from Ctrl-C. Works only if --show-time is used as well.
void mutex_who_is_locked (void)
{
    for (int i=0; i <= MAX_CODE_LINE; i++) {
        LockPoint my_lp = lp[i]; // make a copy for a bit of thread safety
        if (my_lp.mutex_name && my_lp.lock_count)
            printf ("Mutex locked: %s locked in %s:%u. %s\n", 
                    my_lp.mutex_name, (my_lp.func ? my_lp.func : ""), my_lp.code_line, 
                    cond_int (my_lp.lock_count > 1, "num_locks_from_different_objects=", my_lp.lock_count));
    }
}