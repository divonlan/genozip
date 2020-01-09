// ------------------------------------------------------------------
//   win32_pthread.c
//   Copyright (C) 2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <windows.h>
#include "win32_pthread.h"
#include "../genozip.h"

#ifdef _WIN32

typedef struct {
    void *(*start_routine) (void *);
    void *arg;
} ThreadPayload;

// we need an intermediate start routine, because in Windows it returns DWORD (32 bit) and in POSIX
// it returns void* (might be 64 depending on the architecture)
static DWORD pthread_entry (void *payload_)
{
    ThreadPayload *payload = (ThreadPayload*)payload_;
    payload->start_routine (payload->arg); // we ignore the return value - we don't use it in genozip
    return 0;
}
int pthread_create (pthread_t *newthread, const void *unused, void *(*start_routine) (void *), void *arg)
{
    // allocate the payload on the heap to avoid weird problems (and leak the memory... its tiny an infrequent)
    ThreadPayload *payload = malloc (sizeof (ThreadPayload));
    payload->start_routine = start_routine;
    payload->arg = arg;

    *newthread = CreateThread (NULL, 0, pthread_entry, payload, 0, NULL);
    ASSERT (*newthread, "Error: CreateThread() failed: GetLastError()=%lu", GetLastError());
    return 0;
}

int pthread_join (pthread_t th, void **unused)
{
    DWORD ret = WaitForSingleObject (th, INFINITE);  
    ASSERT (ret != WAIT_FAILED, "Error: WaitForSingleObject() failed in pthread_join(): GetLastError()=%lu", GetLastError());
    return 0;
}

int pthread_mutex_init (pthread_mutex_t *mutex, void *unused)
{
    *mutex = CreateMutex (NULL, false, NULL);
    ASSERT (*mutex, "Error: CreateMutex() failed: GetLastError()=%lu", GetLastError());
    return 0;
}

int pthread_mutex_destroy (pthread_mutex_t *mutex)
{
    CloseHandle (*mutex);
    *mutex = NULL;
    return 0;
}

int pthread_mutex_lock (pthread_mutex_t *mutex)
{
    DWORD ret = WaitForSingleObject (*mutex, INFINITE);  
    ASSERT (ret != WAIT_FAILED, "Error: WaitForSingleObject() failed in pthread_mutex_lock(): GetLastError()=%lu", GetLastError());
    return 0;
}

int pthread_mutex_unlock (pthread_mutex_t *mutex)
{
    ASSERT (ReleaseMutex (*mutex), "Error: ReleaseMutex() failed: GetLastError()=%lu", GetLastError());
    return 0;
}
#endif