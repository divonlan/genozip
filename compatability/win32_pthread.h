// ------------------------------------------------------------------
//   win32_pthread.h
//   Copyright (C) 2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// this is an implementation of a subset of pthread using Win32

#ifdef _WIN32

#ifndef WIN32_PTHREAD_INCLUDED
#define WIN32_PTHREAD_INCLUDED

// thread functions
typedef void *pthread_t; // same as Windows HANDLE
extern int pthread_create (pthread_t *newthread, const void *unused, 
                           void *(*start_routine) (void *), void *arg);

extern int pthread_join (pthread_t th, void **unused);


// mutex functions
typedef void *pthread_mutex_t; // same as Windows HANDLE
extern int pthread_mutex_init (pthread_mutex_t *mutex, void *unused);
extern int pthread_mutex_destroy (pthread_mutex_t *mutex);
extern int pthread_mutex_lock (pthread_mutex_t *mutex);
extern int pthread_mutex_unlock (pthread_mutex_t *mutex);

#endif

#endif