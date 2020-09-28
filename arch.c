// ------------------------------------------------------------------
//   arch.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/utsname.h>
#include <termios.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#else // LINUX
#include <sched.h>
#include <sys/sysinfo.h>
#endif
#endif

#include "genozip.h"
#include "endianness.h"
#include "url.h"
#include "arch.h"
#include "sections.h"

static pthread_t io_thread_id = 0; // thread ID of I/O thread (=main thread) - despite common wisdom, it is NOT always 0 (on Windows it is 1)

void arch_initialize(void)
{
    // verify CPU architecture and compiler is supported
    ASSERT0 (sizeof(char)==1 && sizeof(short)==2 && sizeof (unsigned)==4 && sizeof(long long)==8, 
             "Error: Unsupported C type lengths, check compiler options");
    
    // verify endianity is as expected
    uint16_t test_endianity = 0x0102;
#if defined __LITTLE_ENDIAN__
    ASSERT0 (*(uint8_t*)&test_endianity==0x02, "Error: expected CPU to be Little Endian but it is not");
#elif defined __BIG_ENDIAN__
    ASSERT0 (*(uint8_t*)&test_endianity==0x01, "Error: expected CPU to be Big Endian but it is not");
#else
#error  "Neither __BIG_ENDIAN__ nor __LITTLE_ENDIAN__ is defined - is endianness.h included?"
#endif

// Verify that this Windows is 64 bit
#ifdef _WIN32
#ifndef _WIN64
#error Compilation error - on Windows, genozip must be compiled as a 64 bit application
#endif
#endif

    // verify that type sizes are as required (types that appear in section headers written to the genozip format)
    ASSERT0 (sizeof (SectionType)    ==1, "Error: expecting sizeof (SectionType)==1");
    ASSERT0 (sizeof (Codec) ==1, "Error: expecting sizeof (Codec)==1");
    ASSERT0 (sizeof (LocalType)      ==1, "Error: expecting sizeof (LocalType)==1");

    io_thread_id = pthread_self();
}

bool arch_am_i_io_thread (void)
{
    return pthread_self() == io_thread_id;
}

void cancel_io_thread (void)
{
    pthread_cancel (io_thread_id);     // cancel the main thread
    usleep (200000);                   // wait 200ms for the main thread to die. pthread_join here hangs on Windows (not tested on others)
    //pthread_join (io_thread_id, NULL); // wait for the thread cancelation to complete
}

unsigned arch_get_num_cores (void)
{
#ifdef _WIN32
    char *env = getenv ("NUMBER_OF_PROCESSORS");
    if (!env) return DEFAULT_MAX_THREADS;

    unsigned num_cores;
    int ret = sscanf (env, "%u", &num_cores);
    return ret==1 ? num_cores : DEFAULT_MAX_THREADS; 

#elif defined __APPLE__
    int num_cores;
    size_t len = sizeof(num_cores);
    if (sysctlbyname("hw.activecpu", &num_cores, &len, NULL, 0) &&  
        sysctlbyname("hw.ncpu", &num_cores, &len, NULL, 0))
            return DEFAULT_MAX_THREADS; // if both failed

    return (unsigned)num_cores;
 
#else // Linux etc
    // this works correctly with slurm too (get_nprocs doesn't account for slurm core allocation)
    cpu_set_t cpu_set_mask;
    extern int sched_getaffinity (__pid_t __pid, size_t __cpusetsize, cpu_set_t *__cpuset);
    sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);
    unsigned cpu_count = __sched_cpucount (sizeof (cpu_set_t), &cpu_set_mask);
    // TODO - sort out include files so we don't need this extern

    // if failed to get a number - fall back on good ol' get_nprocs
    if (!cpu_count) cpu_count = get_nprocs();

    return cpu_count;
#endif
}

const char *arch_get_os (void)
{
    static char os[256];

#ifdef _WIN32
    uint32_t windows_version = GetVersion();

    sprintf (os, "Windows %u.%u.%u", LOBYTE(LOWORD(windows_version)), HIBYTE(LOWORD(windows_version)), HIWORD(windows_version));
#else

    struct utsname uts;
    ASSERT (!uname (&uts), "Error: uname failed: %s", strerror (errno));

    sprintf (os, "%s %s", uts.sysname, uts.release);

#endif

    return os;
}

const char *arch_get_ip_addr (void)
{
    static char ip_str[1000];

    url_read_string ("https://api.ipify.org", ip_str, sizeof(ip_str));

    return ip_str;
}



bool arch_am_i_in_docker (void)
{
    FILE *fp = fopen ("/.dockerenv", "r");
    if (!fp) return false;

    fclose (fp);
    return true;
}

// initialize mutex, if its not initialized already
void mutex_initialize_do (const char *name, pthread_mutex_t *mutex, bool *initialized)
{
    if (*initialized) return;

    unsigned ret = pthread_mutex_init (mutex, NULL);
    ASSERT (!ret, "Error: pthread_mutex_init failed for %s: %s", name, strerror (ret));

    *initialized = true;
}
