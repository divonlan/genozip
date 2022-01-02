// ------------------------------------------------------------------
//   arch.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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
#include <sys/utsname.h>
#endif
#endif

#include "genozip.h"
#include "endianness.h"
#include "url.h"
#include "arch.h"
#include "sections.h"
#include "flags.h"

#ifdef _WIN32
// add the genozip path to the user's Path environment variable, if its not already there. 
// Note: We do all string operations in Unicode, so as to preserve any Unicode characters in the existing Path
static void arch_add_to_windows_path (const char *argv0)
{
    const char *backslash = strrchr (argv0, '\\');
    if (!backslash) return; // no directory

    unsigned genozip_path_len = backslash - argv0;
    
    WCHAR genozip_path[genozip_path_len+1];
    MultiByteToWideChar(CP_OEMCP, 0, argv0, -1, genozip_path, genozip_path_len);
    genozip_path[genozip_path_len] = 0; 
    
    HKEY key;
    if (RegOpenKeyEx (HKEY_CURRENT_USER, "Environment", 0, KEY_READ | KEY_WRITE, &key))
        return; // fail silently
    
    DWORD value_type;
    WCHAR value[16384];
    DWORD value_size = sizeof (value) - (genozip_path_len + 1)*2; // size in bytes, not unicode characters
    LSTATUS ret = RegQueryValueExW (key, L"Path", 0, &value_type, (BYTE *)value, &value_size);

    if (ret == ERROR_FILE_NOT_FOUND)
        { value[0]=0; value_size=2; }

    else if (ret != ERROR_SUCCESS || value_type != REG_EXPAND_SZ)
        return; // fail silently

    if (wcsstr (value, genozip_path))
        return; // path already exists
    
    wsprintfW (&value[value_size/2-1], L";%s", genozip_path);

    ret = RegSetValueExW (key, L"Path", 0, REG_EXPAND_SZ, (BYTE *)value, value_size + (genozip_path_len + 1)*2); // ignore errors
    if (ret == ERROR_SUCCESS)
        WARN ("%ls has been add to your Path. It will take effect after the next Windows restart.", genozip_path);
} 
#endif

void arch_initialize (const char *argv0)
{
    // verify CPU architecture and compiler is supported
    ASSERT0 (sizeof(char)==1 && sizeof(short)==2 && sizeof (unsigned)==4 && sizeof(long long)==8, 
             "Unsupported C type lengths, check compiler options");
    
    // verify endianity is as expected
    arch_get_endianity();

// Verify that this Windows is 64 bit
#ifdef _WIN32
#ifndef _WIN64
#error Compilation error - on Windows, genozip must be compiled as a 64 bit application
#endif
#endif

    // verify that type sizes are as required (types that appear in section headers written to the genozip format)
    ASSERT0 (sizeof (SectionType) == 1,  "expecting sizeof (SectionType)==1");
    ASSERT0 (sizeof (Codec)       == 1,  "expecting sizeof (Codec)==1");
    ASSERT0 (sizeof (LocalType)   == 1,  "expecting sizeof (LocalType)==1");
    ASSERT0 (sizeof (uint128_t)   == 16, "expecting sizeof (uint128_t)==16");

    // verify that order of bit fields in a structure is as expected (this is compiler-implementation dependent, and we go by gcc)
    // it might be endianity-dependent, and we haven't implemented big-endian yet, see: http://mjfrazer.org/mjfrazer/bitfields/
    union {
        uint8_t byte;
        struct __attribute__ ((__packed__)) { uint8_t a : 1; uint8_t b : 1; } bit_1;
        struct __attribute__ ((__packed__)) { uint8_t a : 3; } bit_3;
    } bittest = { .bit_1 = { .a = 1 } }; // we expect this to set the LSb of .byte and of .bit_3.a
    ASSERT0 (bittest.byte == 1, "unsupported bit order in a struct, please use gcc to compile (1)");
    ASSERT0 (bittest.bit_3.a == 1, "unsupported bit order in a struct, please use gcc to compile (2)");

#ifdef _WIN32
    arch_add_to_windows_path (argv0);
#endif
}

const char *arch_get_endianity (void)
{
    // verify endianity is as expected
    uint16_t test_endianity = 0x0102;
#if defined __LITTLE_ENDIAN__
    ASSERT0 (*(uint8_t*)&test_endianity==0x02, "expected CPU to be Little Endian but it is not");
    return "little";
#elif defined __BIG_ENDIAN__
    ASSERT0 (*(uint8_t*)&test_endianity==0x01, "expected CPU to be Big Endian but it is not");
    return "big";
#else
#error  "Neither __BIG_ENDIAN__ nor __LITTLE_ENDIAN__ is defined - is endianness.h included?"
#endif    
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

    sprintf (os, "Windows_%u.%u.%u", LOBYTE(LOWORD(windows_version)), HIBYTE(LOWORD(windows_version)), HIWORD(windows_version));
#else

    struct utsname uts;
    ASSERT (!uname (&uts), "uname failed: %s", strerror (errno));

    sprintf (os, "%s_%s", uts.sysname, uts.release);

#endif

    return os;
}

const char *arch_get_ip_addr (const char *reason) // optional text in case curl execution fails
{
    static char ip_str[ARCH_IP_LEN] = "0.0.0.0"; // default in case of failure

    url_read_string ("https://api.ipify.org", ip_str, sizeof(ip_str)); // ignore failure

    return ip_str;
}

const char *arch_get_user_host (void)
{
    static char user_host[200];

#ifdef _WIN32
    const char *host = getenv ("HOSTNAME");
    if (!host) host = "";
#else
    char host[100] = "";
    gethostname (host, sizeof (host)-1);
#endif

    const char *user = getenv (flag.is_windows ? "USERNAME" : "USER");

    sprintf (user_host, "%.99s@%.99s", user ? user : "", host);
    return user_host;
}

bool arch_is_wsl (void)
{
#ifdef __linux__    
    struct utsname uts = {};
    if (uname(&uts)) return false; // uname doesn't work

    return (strstr (uts.release, "Microsoft"/*WSL1*/) || strstr (uts.release, "microsoft-standard"/*WSL2*/));

#else
    return false;
#endif
}

#ifndef DISTRIBUTION
    #define DISTRIBUTION "unknown" // this occurs if the code is built not using the genozip Makefile
#endif    
const char *arch_get_distribution (void)
{
    return DISTRIBUTION[0] ? DISTRIBUTION : "github"; // DISTRIBUTION is "" if genozip is built with "make" without defining DISTRIBUTION - re-write as "github"
}
 