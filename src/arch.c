// ------------------------------------------------------------------
//   arch.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <locale.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#include <fcntl.h>
#include <psapi.h>
#else // Mac and Linux
#include <sys/utsname.h>
#include <termios.h>
#include <sys/resource.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach-o/dyld.h>
#include <sys/param.h>
#include <sys/mount.h>
#else // LINUX
#include <sys/sysinfo.h>
#include <sys/vfs.h>
#include <gnu/libc-version.h>
#endif
#endif

#include "genozip.h"
#include "endianness.h"
#include "url.h"
#include "arch.h"
#include "sections.h"
#include "flags.h"
#include "strings.h"
#include "file.h"

static rom argv0 = NULL, base_argv0 = NULL;
Timestamp arch_start_time = 0;

#ifdef _WIN32

#ifdef _M_IX86
#error Genozip does not work on 32-bit Windows
#endif

// add the genozip path to the user's Path environment variable, if its not already there. 
// Note: We do all string operations in Unicode, so as to preserve any Unicode characters in the existing Path
static void arch_add_to_windows_path (void)
{
    rom backslash = strrchr (argv0, '\\');
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

void arch_set_locale (void)
{
#ifdef _WIN32 // see bug 679
    // ASSERTWD (SetThreadLocale (LOCALE_USER_DEFAULT),
    ASSERTWD (SetThreadLocale (LOCALE_INVARIANT), // same as En_US
              "Warning: failed SetThreadLocale: %s", str_win_error());
#else 
    ASSERTWD0 (setlocale (LC_CTYPE,   flag.is_windows ? ".UTF-8" : "en_US.UTF-8"), "Warning: failed to setlocale of LC_CTYPE");   // accept and print UTF-8 text (to do: doesn't work for Windows)
#endif
    ASSERTWD0 (setlocale (LC_NUMERIC, flag.is_windows ? "english" : "en_US.UTF-8"), "Warning: failed to setlocale of LC_NUMERIC"); // force printf's %f to use '.' as the decimal separator (not ',') (required by the Genozip file format)
}

static bool arch_is_wsl (void)
{
#ifdef __linux__    
    struct utsname uts = {};
    if (uname(&uts)) return false; // uname doesn't work

    return (strstr (uts.release, "Microsoft"/*WSL1*/) || strstr (uts.release, "microsoft-standard"/*WSL2*/));

#else
    return false;
#endif
}

void arch_initialize (rom my_argv0)
{
    argv0 = my_argv0;

    rom slash = strrchr (argv0, '/');
    if (!slash && flag.is_windows) slash = strrchr (argv0, '\\');

    base_argv0 = slash ? slash + 1 : argv0;

    // verify CPU architecture and compiler is supported
    ASSERT0 (sizeof(char)==1 && sizeof(short)==2 && sizeof (unsigned)==4 && sizeof(long long)==8, 
             "Unsupported C type lengths, check compiler options");
    
    // verify endianity is as expected
    ASSERT0 (!strcmp (arch_get_endianity(), "little"), "Genozip is currently not supported on big endian architectures");

// Verify that this Windows is 64 bit
#ifdef _WIN32
#ifndef _WIN64
#error Compilation error - on Windows, genozip must be compiled as a 64 bit application
#endif
#endif

    // verify that type sizes are as required (types that appear in section headers written to the genozip format)
    ASSERT0 (sizeof (SectionType)   == 1,  "expecting sizeof (SectionType)==1");
    ASSERT0 (sizeof (Codec)         == 1,  "expecting sizeof (Codec)==1");
    ASSERT0 (sizeof (LocalType)     == 1,  "expecting sizeof (LocalType)==1");
    ASSERT0 (sizeof (uint128_t)     == 16, "expecting sizeof (uint128_t)==16");
    ASSERT0 (sizeof (ReconPlanItem) == 12, "expecting sizeof (ReconPlanItem)==12");
    ASSERT0 (sizeof (void *)        <= 8,  "expecting sizeof (void *)<=8"); // important bc void* is a member of ValueType, and also counting on it in huffman_uncompress
    ASSERT0 (sizeof (ValueType)     == 8,  "expecting sizeof (ValueType)==8");

    // Note: __builtin_clzl is inconsistent between Windows and Linux, even on the same host, so we don't use it
    ASSERT0 (__builtin_clz(5)   == 29, "expecting __builtin_clz to be 32 bit");
    ASSERT0 (__builtin_clzll(5) == 61, "expecting __builtin_clzll to be 64 bit");

    // verify that order of bit fields in a structure is as expected (this is compiler-implementation dependent, and we go by gcc)
    // it might be endianity-dependent, and we haven't implemented big-endian yet, see: http://mjfrazer.org/mjfrazer/bitfields/
    union {
        uint8_t byte;
        struct __attribute__ ((packed)) { uint8_t a : 1; uint8_t b : 1; } bit_1;
        struct __attribute__ ((packed)) { uint8_t a : 3; } bit_3;
    } bittest = { .bit_1 = { .a = 1 } }; // we expect this to set the LSb of .byte and of .bit_3.a
    ASSERT0 (bittest.byte == 1, "unsupported bit order in a struct, please use gcc to compile (1)");
    ASSERT0 (bittest.bit_3.a == 1, "unsupported bit order in a struct, please use gcc to compile (2)");

    arch_set_locale();

#ifdef _WIN32
    _setmode(_fileno(stdin),  _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);

    arch_add_to_windows_path();
#endif
    
    flag.is_wsl = arch_is_wsl();
    arch_start_time = arch_timestamp();
}

rom arch_get_endianity (void)
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

// physical RAM size in GB
double arch_get_physical_mem_size (void)
{
    static double mem_size = 0;
    if (mem_size) return mem_size;

#ifdef __linux__    
    FILE *fp = fopen ("/proc/meminfo", "rb");
    if (!fp) return 0;

    char meminfo[100] = {}; // note: can't use a Buffer because called from a signal handler - we don't know which is the running thread
    fread (meminfo, 1, sizeof(meminfo)-1, fp); // -1 to guarantee that (at least) last character is \0 

    int num_start = strcspn (meminfo, "0123456789");
    mem_size = (double)atoll(&meminfo[num_start]) / (1024.0*1024.0); // convert KB to GB
    
    fclose (fp);

#elif defined _WIN32
    ULONGLONG kb = 0;
    GetPhysicallyInstalledSystemMemory (&kb);
    mem_size = (double)kb / (1024.0*1024.0);

#elif defined __APPLE__
    int64_t bytes = 0;
    size_t length = sizeof (int64_t);
    sysctl((int[]){ CTL_HW, HW_MEMSIZE }, 2, &bytes, &length, NULL, 0);
    mem_size = (double)bytes / (1024.0*1024.0*1024.0);

#else
    return 0;

#endif

    return mem_size;
}


StrText arch_get_filesystem_type (FileP file)
{
    StrText s = { "unknown" };
    int save_errno = errno; // save errno, as this function is often used in ASSERT.

    if (file && file->is_remote) {
        strcpy (s.s, "remote");
        goto done;
    }

    if (file && file->redirected && !file->name) {
        strcpy (s.s, "pipe");
        goto done;
    }

    if (!file || !file->file || !file->name) 
        goto done;

#ifdef __linux__    
    struct statfs fs;
    if (statfs (file->name, &fs)) goto done; 

    rom name = NULL;
    #define NAME(magic, name_s) case magic: name = name_s; break
    switch (fs.f_type) {
        NAME (0xff534d42, "CIFS");
        NAME (0xef53,     "ext2/3/4");
        NAME (0x6969,     "nfs");
        NAME (0x5346544e, "NTFS");
        NAME (0x858458f6, "ramfs");
        NAME (0x58465342, "xfs");
        NAME (0x01021997, "v9fs");     // Used by WSL2
        NAME (0x0bd00bd0, "Lustre");   // HPC filesystem: https://www.lustre.org/
        NAME (0x65735546, "FUSE");     // Filesystem in user space: https://www.kernel.org/doc/html/next/filesystems/fuse.html
        NAME (0xaad7aaea, "PanFS");    // Clustered filesystem: https://www.panasas.com/products/panfs/
        NAME (0xc36400,   "CephFS");   // Distributed filesystem: https://docs.ceph.com/en/latest/cephfs/
        NAME (0x47504653, "GPFS");     // IBM Spectrum Scale GPFS: https://www.ibm.com/docs/en/storage-scale/4.2.0?topic=scale-overview-gpfs
        NAME (0xfe534d42, "SMB2");     // Windows file sharing
        NAME (0x2fc12fc1, "ZFS");      // Oracle ZFS (originally in Solaris) https://docs.oracle.com/cd/E19253-01/819-5461/zfsover-2/
        NAME (0x19830326, "FhGFS");    // https://www.beegfs.io/docs/SC13_FHGFS_Presentation.pdf
        NAME (0x53464846, "wslfs");    // WSL1: https://github.com/MicrosoftDocs/WSL/issues/465
        NAME (0x1021994,  "tmpfs");    // Heap Backing Filesystem
        NAME (0x2011bab0, "exFAT");    // Filesystem for flash memory: https://en.wikipedia.org/wiki/ExFAT
        NAME (0x9123683e, "brtfs");    // Copy-on-write filesystem for Linux: https://docs.kernel.org/filesystems/btrfs.html
        NAME (0x794C7630, "OverlayFS");// A union-mount filesystem: https://en.wikipedia.org/wiki/OverlayFS
        NAME (0xf15f,     "eCryptfs"); // A cryptographic filesystem for Linux: https://www.ecryptfs.org/
        default: snprintf (s.s, sizeof (s.s), "0x%lx", fs.f_type); 
    }

    if (name) strcpy (s.s, name);

#elif defined __APPLE__
    struct statfs fs;
    if (statfs (file->name, &fs)) goto done; 

    memcpy (s.s, fs.f_fstypename, MIN_(sizeof(fs.f_fstypename), sizeof(s)-1));

#elif defined _WIN32
    WCHAR ws[100];
    if (!GetVolumeInformationByHandleW ((HANDLE)_get_osfhandle(fileno (file->file)), 0, 0, 0, 0, 0, ws, ARRAY_LEN(ws))) goto done;

    if (wcstombs (s.s, ws, sizeof(s.s)-1) == (size_t)-1)
        strcpy (s.s, "failed-wcstombs"); // can happen if locale is set to non-english
#endif    

done:
    errno = save_errno;
    return s;
} 

StrText arch_get_txt_filesystem (void)
{
    return arch_get_filesystem_type (txt_file);
}

StrText arch_get_z_filesystem (void)
{
    return arch_get_filesystem_type (z_file);
}

// returns value in bytes
uint64_t arch_get_max_resident_set (void)
{
#ifndef _WIN32
    // Linux and MacOS - get maximal RSS this process ever had (TO DO: get current resident set)
    struct rusage usage;
    if (getrusage (RUSAGE_SELF, &usage) != 0) return 0; // failed
    return usage.ru_maxrss KB;

#else  
    // Windows - get *current* working set
    HANDLE process = OpenProcess (PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
    if (!process) return 0;

    PROCESS_MEMORY_COUNTERS_EX mem_counters = {};
    GetProcessMemoryInfo (process, (PROCESS_MEMORY_COUNTERS*)&mem_counters, sizeof (mem_counters));
    CloseHandle (process);

    return mem_counters.WorkingSetSize; // still 0 if error
#endif
}

rom arch_get_os (void)
{
    static char os[64];

#ifdef _WIN32
    uint32_t windows_version = GetVersion();

    snprintf (os, sizeof (os), "Windows_%u.%u.%u", LOBYTE(LOWORD(windows_version)), HIBYTE(LOWORD(windows_version)), HIWORD(windows_version));
#else

    struct utsname uts;
    ASSERT (!uname (&uts), "uname failed: %s", strerror (errno));

    snprintf (os, sizeof (os), "%.30s_%.30s", uts.sysname, uts.release);

#endif

    return os;
}

rom arch_get_scheduler (void)
{
    if (getenv ("SLURM_JOB_ID"))            return "slurm";
    if (getenv ("KUBERNETES_SERVICE_HOST")) return "kubernetes";
    if (getenv ("LSB_JOBID"))               return "LSF";
    if (arch_is_docker())                   return "docker";

    return NULL; // not scheduler identified
}

rom arch_get_glibc (void)
{
#ifdef __linux__
    return gnu_get_libc_version();
#else
    return "not_glibc";
#endif
}

// good summary here: https://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe/1024937#1024937
// returns nul-terminated executable path
StrTextSuperLong arch_get_executable (void) 
{
    StrTextSuperLong path = {};

#ifdef __linux__    
    ssize_t path_len = readlink ("/proc/self/exe", path.s, sizeof(path.s) - 1); // doesn't nul-terminate
    ASSGOTO (path_len > 0, "readlink() failed to get executable path from /proc/self/exe: %s", strerror(errno));

#elif defined _WIN32
    DWORD path_len = GetModuleFileNameA (NULL, path.s, sizeof(path.s) - 1); 
    if (GetLastError() == ERROR_INSUFFICIENT_BUFFER) path_len = 0; // this is also an error

    ASSGOTO (path_len, "GetModuleFileNameA() failed: %s", str_win_error());

#elif defined __APPLE__
    // see: https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/dyld.3.html
    uint32_t path_len = 0;
    _NSGetExecutablePath (NULL, &path_len); // get path len - possibly more than MAXPATHLEN if has symlinks
    path_len = MIN_(path_len, sizeof(path.s) - 1);

    ASSGOTO0 (!_NSGetExecutablePath (path.s, &path_len), "_NSGetExecutablePath() failed"); 

#else // another OS
    goto error;

#endif
    path.s[sizeof(path.s)-1] = 0;
    return path;

error:
    ASSERT (strlen (argv0) < sizeof (path.s), "executable name by argv[0] is longer than Genozip's maximum of %u characters: \"%s\"", (int)sizeof(path.s)-1, argv0);
    strcpy (path.s, argv0);
    return path;
}

StrTextSuperLong arch_get_genozip_executable (void)
{
    StrTextSuperLong fn = arch_get_executable();

    if (!is_genozip) {
        rom bn       = is_genounzip?"genounzip" : is_genocat?"genocat" : "genols";
        int bn_len   = is_genounzip?9           : is_genocat?7         : 6;

        // replace filename if possible
        char *loc = strstr (fn.s, bn);
        if (loc) {  
            memmove (loc + 7, loc + bn_len, strlen (loc+bn_len) + 1/*\0*/);
            memcpy (loc, "genozip", 7);
        }
        
        // note: do nothing is is_genounzip - this is likely genozip --decompress
        else if (!is_genounzip)
            ABORT ("Cannot find substring %s in %s", bn, fn.s);
    }

    return fn;
}



rom arch_get_argv0 (void)
{
    return base_argv0;
}
 
// true if running under valgrind
bool arch_is_valgrind (void)
{
    static thool is_valgrind = unknown;

    if (is_valgrind == unknown) {
        rom p = getenv ("LD_PRELOAD");
        is_valgrind = flag.debug_valgrind || (p && (strstr (p, "/valgrind/") || strstr (p, "/vgpreload")));
    }

    return is_valgrind;
}

bool arch_is_docker (void)
{
    static thool is_docker = unknown; 
    
    return (is_docker != unknown) ? is_docker : (is_docker = file_exists ("/.dockerinit"));
}

Timestamp inline arch_timestamp (void) 
{
    struct timespec tb;
    clock_gettime (CLOCK_REALTIME, &tb);
    return (uint128_t)tb.tv_sec * 1000000000 + (uint128_t)tb.tv_nsec;
}

bool arch_is_process_alive (uint32_t pid)
{
#ifndef _WIN32
    bool is_alive = (getpgid(pid) >= 0); // test its process group id which is always possible even for processes belong to other users
#else
    HANDLE process = OpenProcess (PROCESS_QUERY_LIMITED_INFORMATION, false, pid);
    
    DWORD exit_code;
    bool is_alive = process && GetExitCodeProcess (process, &exit_code) && (exit_code == STILL_ACTIVE);

    CloseHandle (process);
#endif
    return is_alive;
}

// check if executable is in the path. Note: in Windows it actually runs the executable, so 
// only suitable for executables that would terminate immediately.
static bool arch_is_exec_in_path (rom exec)
{
#ifdef _WIN32
    StreamP where = stream_create (0, 0, 0, 0, 0, 0, 0, "where.exe", "where.exe", "/Q", exec, NULL); 

    return stream_close (&where, STREAM_WAIT_FOR_PROCESS) == 0;

#else
    char run[32 + strlen(exec)];
    snprintf (run, sizeof (run), "which %s > /dev/null 2>&1", exec);
    return !system (run) && file_exists ("/dev/stdout");
#endif
}

bool wget_available (void)
{   
    static thool installed = unknown;
    if (installed == unknown)
        // note: wget not used on Windows, bc I can't get it to output to stdout, and also earlier wget versions may be adding \r ... : https://stackoverflow.com/questions/8522983/wget-of-binary-file-piped-into-other-commands-on-windows-breaks-the-binary
        installed = !flag.is_windows && arch_is_exec_in_path ("wget");

    return installed;
}

bool curl_available (void)
{
    static thool installed = unknown;
    if (installed == unknown)
        installed = arch_is_exec_in_path ("curl");

    return installed;
}

#ifdef santize_thread
void *__gxx_personality_v0; // overcome "undefined reference to '__gxx_personality_v0'" when linking with --sanitize=thread
#endif
