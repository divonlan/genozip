// ------------------------------------------------------------------
//   mutex.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include "genozip.h"
#include "mutex.h"
#include "flags.h"

void mutex_initialize_do (Mutex *mutex, const char *name, const char *func)
{ 
    if (mutex->initialized) return;

    int ret = pthread_mutex_init (&mutex->mutex, 0); 
    ASSERT (!ret || errno == EBUSY,  // EBUSY is not an error - the failure is bc a race condition and the mutex is already initialized - all good
            "pthread_mutex_init failed for %s from %s: %s", name, func, strerror (ret)); 

    mutex->name = name;
    mutex->initialized = func;
    mutex->vb_i_last = 0; // reset
}

void mutex_destroy_do (Mutex *mutex, const char *func) 
{
    if (!mutex->initialized) return;
        
    pthread_mutex_destroy (&mutex->mutex); 
    memset (mutex, 0, sizeof (Mutex));
}

bool mutex_lock_do (Mutex *mutex, bool blocking, const char *func)   
{ 
    ASSERT (mutex->initialized, "called from %s: mutex not initialized", func);

    bool show = mutex_is_show (mutex->name);

    if (show) iprintf ("LOCKING : Mutex %s by thread %"PRIu64" %s\n", mutex->name, (uint64_t)pthread_self(), func);

    int ret = blocking ? pthread_mutex_lock (&mutex->mutex)
                       : pthread_mutex_trylock (&mutex->mutex);

    if (!blocking && ret == EBUSY) return false;

    ASSERT (!ret, "called from %s by thread=%"PRIu64": pthread_mutex_lock failed on mutex->name=%s: %s", 
            func, (uint64_t)pthread_self(), mutex && mutex->name ? mutex->name : "(null)", strerror (ret)); 

    mutex->lock_func = func; // mutex->lock_func is protected by the mutex

    if (show) iprintf ("LOCKED  : Mutex %s by thread %"PRIu64"\n", mutex->name, (uint64_t)pthread_self());

    return true;
}

void mutex_unlock_do (Mutex *mutex, const char *func, uint32_t line) 
{ 
    ASSERT (mutex->initialized, "called from %s:%u mutex not initialized", func, line);
    ASSERT (mutex->lock_func, "called from %s:%u by thread=%"PRIu64": mutex %s is not locked", 
            func, line, (uint64_t)pthread_self(), mutex->name);

    mutex->lock_func = NULL; // mutex->lock_func is protected by the mutex

    int ret = pthread_mutex_unlock (&mutex->mutex); 
    ASSERT (!ret, "called from %s:%u: pthread_mutex_unlock failed for %s: %s", func, line, mutex->name, strerror (ret)); 

    if (mutex_is_show (mutex->name))
        iprintf ("UNLOCKED: Mutex %s by thread %"PRIu64" %s\n", mutex->name, (uint64_t)pthread_self(), func);
}

void mutex_wait_do (Mutex *mutex, const char *func, uint32_t line)   
{
    mutex_lock_do (mutex, true, func);
    mutex_unlock_do (mutex, func, line);
}

void serializer_lock_do (MutexP ser, uint32_t vb_i, const char *func, uint32_t code_line)
{
    #define WAIT_TIME_USEC 5000
    #define TIMEOUT (30*60) // 30 min

    for (unsigned i=0; ; i++) {
        mutex_lock_do ((MutexP)ser, true, func);

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
        ASSERT (i < TIMEOUT*(1000000/WAIT_TIME_USEC), "called from %s:%u: Timeout (%u sec) while waiting for serializer %s in vb=%u. vb_i_last=%u", 
                func, code_line, TIMEOUT, ser->name, vb_i, ser->vb_i_last);
    }
}
