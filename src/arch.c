// ------------------------------------------------------------------
//   arch.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
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
        WARN (_FYI "%ls has been add to your Path. It will take effect after the next Windows restart.", genozip_path);
} 
#endif

void arch_set_locale (void)
{
#ifdef _WIN32 
    // set Windows Console to UTF-8 input/output
    SetConsoleCP(65001);       
    SetConsoleOutputCP(65001);
#else 
    // accept and print UTF-8 text (doesn't work for Windows - although printing UTF-8 works fine anyway, and we don't allow non-ASCII in --activate)
    ASSERTWD (setlocale (LC_CTYPE, "C.UTF-8"), _WRN "failed to setlocale of LC_CTYPE", NULL);   
#endif
    // force printf's %f to use '.' as the decimal separator (not ',') (required by the Genozip file format)
    ASSERTWD (setlocale (LC_NUMERIC, 𝓌𝒾𝓃("english") X𝓌𝒾𝓃("C.UTF-8")), _WRN "failed to setlocale of LC_NUMERIC", NULL); 
}

static bool arch_is_wsl (void)
{
#ifdef __linux__    
    // case: Linux executable is running in WSL
    struct utsname uts = {};
    if (uname(&uts)) return false; // uname doesn't work
    return (strstr (uts.release, "Microsoft"/*WSL1*/) || strstr (uts.release, "microsoft-standard"/*WSL2*/));

#elif defined _WIN32
    // case: Windows executable is called from WSL
    rom env = getenv("WSLENV");
    return env && env[0]; // existing and not ""

#else
    return false;

#endif
}

void arch_initialize (rom my_argv0)
{
    argv0 = my_argv0;

    // We require AVX2 support in x86
    χ64 (ASSINP0 (__builtin_cpu_supports("avx2"), "Genozip, running on Intel/AMD CPUs, requires AVX2 support");)
    χ64 (ASSINP0 (__builtin_cpu_supports("bmi2"), "Genozip, running on Intel/AMD CPUs, requires BMI2 support");)

    rom slash = strrchr (argv0, '/');
    𝓌𝒾𝓃 (if (!slash) slash = strrchr (argv0, '\\');)

    base_argv0 = slash ? slash + 1 : argv0;

    // verify CPU architecture and compiler is supported: compile-time tests
    ASSERT_SIZEOF (char,          1);
    ASSERT_SIZEOF (short,         2);
    ASSERT_SIZEOF (unsigned,      4);
    ASSERT_SIZEOF (long long,     8);
    ASSERT_SIZEOF (uint128_t,     16);
    ASSERT_SIZEOF (size_t,        8);
    ASSERT_SIZEOF (time_t,        8);
    ASSERT_SIZEOF (SectionType,   1);
    ASSERT_SIZEOF (Codec,         1);
    ASSERT_SIZEOF (LocalType,     1);
    ASSERT_SIZEOF (StoreType,     1);
    ASSERT_SIZEOF (void *,        8);
    ASSERT_SIZEOF (ValueType,     8);        
    ASSERT_SIZEOF (Buffer,        64);
    ASSERT_SIZEOF (SnipIterator,  8);

    // file-format structs
    ASSERT_SIZEOF (struct FlagsCtx,           1);
    ASSERT_SIZEOF (struct FlagsDict,          1);
    ASSERT_SIZEOF(SectionFlags,               1);
    ASSERT_SIZEOF(SectionHeader,              28);
    ASSERT_SIZEOF(SectionHeaderGenozipHeader, 720);
    ASSERT_SIZEOF(SectionFooterGenozipHeader, 12);
    ASSERT_SIZEOF(QnameFlavorProp,            2);
    ASSERT_SIZEOF(SectionHeaderTxtHeader,     400);
    ASSERT_SIZEOF(SectionHeaderVbHeader,      84);
    ASSERT_SIZEOF(VbPlanItem,                 2);
    ASSERT_SIZEOF(SectionHeaderDictionary,    40);
    ASSERT_SIZEOF(SectionHeaderCounts,        44);
    ASSERT_SIZEOF(SectionHeaderSubDicts,      44);
    ASSERT_SIZEOF(SectionHeaderHuffman,       36);
    ASSERT_SIZEOF(SectionHeaderGzDigests,     37);
    ASSERT_SIZEOF(SectionHeaderCtx,           40);
    ASSERT_SIZEOF(SectionHeaderReference,     52);
    ASSERT_SIZEOF(SectionHeaderRefHash,       36);
    ASSERT_SIZEOF(SectionHeaderReconPlan,     36);
    ASSERT_SIZEOF(ReconPlanItem,              12);
    ASSERT_SIZEOF(GencompSecItem,             12);
    ASSERT_SIZEOF(SectionEntFileFormat,       19);
    ASSERT_SIZEOF(SectionEntFileFormatV14,    24);
    ASSERT_SIZEOF(RAEntry,                    24);
    ASSERT_SIZEOF(Iupac,                      9);    
    
    // verify endianity is as expected
    ASSINP0 (!strcmp (arch_get_endianity(), "little"), "Genozip is currently not supported on big endian architectures");

// Verify that this Windows is 64 bit
#ifdef _WIN32
#ifndef _WIN64
#error Compilation error - on Windows, genozip must be compiled as a 64 bit application
#endif
#endif

    // Note: __builtin_clzl is inconsistent between Windows and Linux, even on the same host, so we don't use it
    ASSINP0 (__builtin_clz(5)   == 29, "expecting __builtin_clz to be 32 bit");
    ASSINP0 (__builtin_clzll(5) == 61, "expecting __builtin_clzll to be 64 bit");

    // verify that order of bit fields in a structure is as expected (this is compiler-implementation dependent, and we go by gcc)
    // it might be endianity-dependent, and we haven't implemented big-endian yet, see: http://mjfrazer.org/mjfrazer/bitfields/
    union {
        uint8_t byte;
        struct __attribute__ ((packed)) { uint8_t a : 1; uint8_t b : 1; } bit_1;
        struct __attribute__ ((packed)) { uint8_t a : 3; } bit_3;
    } bittest = { .bit_1 = { .a = 1 } }; // we expect this to set the LSb of .byte and of .bit_3.a
    ASSINP0 (bittest.byte == 1, "unsupported bit order in a struct, please use gcc to compile (1)");
    ASSINP0 (bittest.bit_3.a == 1, "unsupported bit order in a struct, please use gcc to compile (2)");

    // verify gcc / SYS-V bit packing, not Microsoft
    // gcc packs to the smallest possible size (32bit in this case) but Windows ABI (inc. MingW) will start a new packing upon type change.
    // in gcc sizeof() is 4, whereas in Windows it is 8 ("a" is padded). gcc option -mno-ms-bitfields fixes this for MingW.
    struct { uint8_t a : 2; uint32_t b : 4; } ms_bitfields_test; 
    ASSINP0 (sizeof (ms_bitfields_test) == 4, "expecting gcc-style bit packing");

    // verify that malloced memory is always on 16B-boundary (assumed for 64B padding in buf_alloc_do)
    void *ptrs[32];
    for (int i=0; i < ARRAY_LEN(ptrs); i++) {
        ptrs[i] = malloc (1999/*prime*/ / (i+1)); // don't free after allocation, to prevent reuse of same block
        ASSERT ((uintptr_t)ptrs[i] % 16 == 0, "expecting malloced memory to be 16B-aligned, but ptr=%p", ptrs[i]);
    }
    for (int i=0; i < ARRAY_LEN(ptrs); i++) free (ptrs[i]);

    arch_set_locale();

#ifdef _WIN32
    _setmode(_fileno(stdin),  O_BINARY);
    _setmode(_fileno(stdout), O_BINARY);

    arch_add_to_windows_path();
#endif
    
    flag.is_wsl = arch_is_wsl();
    arch_start_time = arch_timestamp();

    // initialize the random number generator
    struct timespec ts;
    srand (!clock_gettime (CLOCK_MONOTONIC, &ts) ? ts.tv_nsec : time (NULL));

    // test for valgrind
    rom p = getenv ("LD_PRELOAD");
    flag.is_valgrind |= (p && (strstr (p, "/valgrind/") || strstr (p, "/vgpreload")));
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
    int num_cores = 0;
    if (num_cores) return num_cores;

#ifdef _WIN32
    char *env = getenv ("NUMBER_OF_PROCESSORS");
    if (!env) return DEFAULT_MAX_THREADS;

    int ret = sscanf (env, "%d", &num_cores);
    if (ret != 1) num_cores = DEFAULT_MAX_THREADS;

#elif defined __APPLE__
    size_t len = sizeof(num_cores);
    if (sysctlbyname("hw.activecpu", &num_cores, &len, NULL, 0) &&  
        sysctlbyname("hw.ncpu", &num_cores, &len, NULL, 0))
            num_cores = DEFAULT_MAX_THREADS; // if both failed
 
#else // Linux etc
    // this works correctly with slurm too (get_nprocs doesn't account for slurm core allocation)
    cpu_set_t cpu_set_mask;
    extern int sched_getaffinity (__pid_t __pid, size_t __cpusetsize, cpu_set_t *__cpuset);
    sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);
    num_cores = __sched_cpucount (sizeof (cpu_set_t), &cpu_set_mask);
    // TODO - sort out include files so we don't need this extern

    // if failed to get a number - fall back on good ol' get_nprocs
    if (!num_cores) num_cores = get_nprocs();
#endif

    if (num_cores <= 0) num_cores = DEFAULT_MAX_THREADS; // safety
    
    return (unsigned)num_cores; 
}

// physical RAM size in GB
double arch_get_physical_mem_size (void)
{
    static double mem_size = 0;
    if (mem_size) return mem_size;

#ifdef __linux__    
    FILE *fp = fopen ("/proc/meminfo", READ);
    if (!fp) return 0;

    char meminfo[100] = {}; // note: can't use a Buffer because called from a signal handler - we don't know which is the running thread
    if (fread (meminfo, 1, sizeof(meminfo)-1, fp) > 0) { // -1 to guarantee that (at least) last character is \0 

        int num_start = strcspn (meminfo, "0123456789");
        mem_size = (double)atoll(&meminfo[num_start]) / (1024.0*1024.0); // convert KB to GB
    }
    
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

uint64_t arch_get_shmmax(void) 
{
    uint64_t shmmax = 0;

#if defined __linux__
    ASSERTNOTINUSE (evb->scratch);
    file_get_file (evb, "/proc/sys/kernel/shmmax", &evb->scratch, "scratch", 1 KB, VERIFY_ASCII, true);
    shmmax = atoll (B1STc(evb->scratch));
    buf_free (evb->scratch);

#elif defined __APPLE__
    size_t len = sizeof (shmmax);
    ASSERT (sysctlbyname ("kern.sysv.shmmax", &shmmax, &len, NULL, 0) == 0,
                          "sysctlbyname(kern.sysv.shmmax) failed: %s", strerror(errno));

#elif defined _WIN32
    // Windows has no fixed SHMMAX constraint; it is bound by the system commit limit.
    PERFORMANCE_INFORMATION pi = { .cb = sizeof(PERFORMANCE_INFORMATION) } ;
    
    if (GetPerformanceInfo(&pi, sizeof(pi))) 
        shmmax = (uint64_t)pi.CommitLimit * (uint64_t)pi.PageSize; // the absolute ceiling of virtual memory pages the system can allocate across RAM + Pagefile without running dry. 
    else {    
        MEMORYSTATUSEX statex = { .dwLength = sizeof (statex) };
        if (GlobalMemoryStatusEx (&statex)) 
            shmmax = statex.ullTotalPhys;
    }
    
#endif

    return shmmax;
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
        NAME (0x9123683e, "Btrfs‍");    // A copy-on-write B-tree filesystem for Linux: https://docs.kernel.org/filesystems/btrfs.html
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
    if (!GetVolumeInformationByHandleW ((HANDLE)_get_osfhandle(fileno (file->os_file)), 0, 0, 0, 0, 0, ws, ARRAY_LEN(ws))) goto done;

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

    union {
        PROCESS_MEMORY_COUNTERS base;
        PROCESS_MEMORY_COUNTERS_EX ex;
    } mem_counters;

    GetProcessMemoryInfo (process, &mem_counters.base, sizeof (mem_counters));
    CloseHandle (process);

    return mem_counters.ex.WorkingSetSize; // still 0 if error
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

// true if an environment variable exists with a name that starts with prefix
#ifdef __linux__
static bool has_env_prefix (rom prefix)
{
    extern char **environ;
    if (!environ) return false;

    size_t prefix_len = strlen (prefix);
    for (char **env=environ; *env; env++) 
        if (!strncmp (*env, prefix, prefix_len)) 
            return true;
    return false;
}
#endif    

#define HAS(e) ({ rom val=getenv(e); val && val[0]; })

rom arch_get_scheduler (void)
{
    static rom sched = NULL;

#ifdef __linux__
    DO_ONCE sched = 
           // Traditional HPC Schedulers
          (HAS("CONDOR_JOB_ID") || HAS("_CONDOR_JOB_PIDS"))    ? "htcondor"
         : HAS("SLURM_JOB_ID")            ? "slurm"
         : HAS("PBS_JOBID")               ? "pbs"
         : HAS("LSB_JOBID")               ? "lsf"
         : HAS("FLUX_JOB_ID")             ? "flux"
         : HAS("COBALT_JOBID")            ? "cobalt"
         : HAS("PJM_JOBID")               ? "pjm"
         : has_env_prefix ("SGE_")        ? "sge"

           // Orchestrators & Cloud Batch
         :(HAS("BATCH_TASK_INDEX") || HAS("BATCH_TASK_COUNT")) ? "gcp-batch"
         : HAS("KUBERNETES_SERVICE_HOST") ? "kubernetes"
         : HAS("NOMAD_ALLOC_ID")          ? "nomad"
         : HAS("AWS_BATCH_JOB_ID")        ? "aws-batch"
         : HAS("AZ_BATCH_JOB_ID")         ? "azure-batch"
         : HAS("RAY_JOB_ID")              ? "ray"

           // Workflow Engines
         : HAS("NXF_TASK_WORKDIR")        ? "nextflow"
         : has_env_prefix ("SNAKEMAKE_")  ? "snakemake"  
         : NULL;
#endif
        return sched;
}

// IaaS - Infrastructure as a Service
rom arch_get_cloud (void)
{
    static rom cloud = NULL;

#ifdef __linux__
    DO_ONCE cloud = 
          (HAS("GCE_METADATA_HOST") || (HAS("GOOGLE_CLOUD_PROJECT") && HAS("K_SERVICE"))) ? "gcp"
         :(HAS("AWS_EXECUTION_ENV") || HAS("AWS_REGION")) ? "aws"
         : HAS("OCI_RESOURCE_PRINCIPAL_VERSION")          ? "oci"
         : HAS("ALIBABA_CLOUD_ACCOUNT_ID")                ? "alibaba"
         : HAS("IBM_CLOUD_REGION")                        ? "ibm"
         : HAS("AZURE_HTTP_USER_AGENT ")                  ? "azure"
         : NULL;
#endif
        return cloud;
}

// Platform as a Service
rom arch_get_PaaS (void)
{
    static rom paas = NULL;

#ifdef __linux__
    DO_ONCE paas = 
         // Platform / Serverless / Edge Providers (PaaS) - run on some IaaS platform (AWS etc)
          (HAS("RAILWAY_STATIC_URL") || HAS("RAILWAY_ENVIRONMENT ")) ? "railway"
         : HAS("VERCEL")              ? "vercel"
         : HAS("NETLIFY")             ? "netlify"
         : HAS("RENDER")              ? "render"
         : HAS("FLY_APP_NAME")        ? "fly.io"
         : has_env_prefix ("HEROKU_") ? "heroku"
         : HAS("DIGITALOCEAN")        ? "digitalocean"
         : NULL;
#endif
        return paas;
}

rom arch_get_glibc (void)
{
    return ℓ𝒾𝓃𝓊𝓍(gnu_get_libc_version()) 
           Xℓ𝒾𝓃𝓊𝓍("not_glibc");
}

// good summary here: https://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe/1024937#1024937
// returns nul-terminated executable path
StrText4K arch_get_executable (void) 
{
    StrText4K path = {};

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

    ASSGOTO (!_NSGetExecutablePath (path.s, &path_len), "_NSGetExecutablePath() failed", NULL); 

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

StrText4K arch_get_genozip_executable (void)
{
    StrText4K fn = arch_get_executable();

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
 
Timestamp inline arch_timestamp (void) 
{
    struct timespec tb;
    clock_gettime (CLOCK_REALTIME, &tb);
    return (uint128_t)tb.tv_sec * 1000000000 + (uint128_t)tb.tv_nsec;
}

// seconds since Unix epoch, same as time()
uint64_t arch_time (void)
{
    // in WSL2 there is a issue of a stale clock after sleep / hybernation
    if (!arch_is_wsl()) return time (NULL);

    // bypass WSL clock and get the Windows clock instead
    uint64_t win_time = 0;
    FILE *pipe = popen ("powershell.exe -Command \"[DateTimeOffset]::UtcNow.ToUnixTimeSeconds()\"", "r");
    
    if (pipe) {
        char str[32];
        if (fgets (str, sizeof (str), pipe)) 
            win_time = (uint64_t)atoll (str);
        
        pclose (pipe);
    }
    
    return win_time ? win_time : time (NULL); // fallback to time() if failed to get from Windows
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
    
    // note: wget not used on Windows, bc I can't get it to output to stdout, and also earlier wget versions may be adding \r ... : https://stackoverflow.com/questions/8522983/wget-of-binary-file-piped-into-other-commands-on-windows-breaks-the-binary
    X𝓌𝒾𝓃 (if (installed == unknown) installed = arch_is_exec_in_path ("wget");)
        
    return installed;
}

bool curl_available (void)
{
    static thool installed = unknown;
    if (installed == unknown)
        installed = arch_is_exec_in_path ("curl");

    return installed;
}

#ifdef sanitize_thread
void *__gxx_personality_v0; // overcome "undefined reference to '__gxx_personality_v0'" when linking with --sanitize=thread
#endif

rom arch_str_error (void)
{
    return 𝓌𝒾𝓃(str_win_error()) 
           X𝓌𝒾𝓃(strerror (errno));    
}

static size_t arch_get_l3_cache_size (void)
{
#ifdef __linux__
    return sysconf (_SC_LEVEL3_CACHE_SIZE);

#elif defined __APPLE__
    uint64_t cache_size = 0; 
    size_t size = sizeof(cache_size); 
    
    ASSERT (!sysctlbyname("hw.l3cachesize", &cache_size, &size, NULL, 0),
            "sysctlbyname failed: %s", strerror(errno));
    return cache_size; 

#elif defined _WIN32
    DWORD len = 0;

    GetLogicalProcessorInformationEx (RelationCache, NULL, &len);

    SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *info = MALLOC (len);

    ASSERT (GetLogicalProcessorInformationEx (RelationCache, info, &len),
            "GetLogicalProcessorInformationEx failed: %s", str_win_error());

    size_t largest_l3 = 0;

    for (DWORD offset=0; offset < len; ) {
        SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *p = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *)((uint8_t *)info + offset);

        if (p->Relationship == RelationCache &&
            p->Cache.Level == 3 &&
            p->Cache.CacheSize > largest_l3)
            largest_l3 = p->Cache.CacheSize;

        offset += p->Size;
    }

    FREE (info);

    return largest_l3;
#endif
}

// flush CPU cache (most importantly, L3 cache) - ahead of timing algorithms for consistent timing
void arch_flush_cpu_cache (void)
{
    flag.flush_cpu_cache = true;
    
    size_t l3_bytes = arch_get_l3_cache_size();
    
    // volatile prevents the compiler from eliminating the loop.
    volatile uint8_t *flush_buffer = (volatile uint8_t *)MALLOC(l3_bytes);

    // touch every cache line   
    for (size_t i=0; i < l3_bytes; i += 64)
        flush_buffer[i]++;

    __sync_synchronize();
    FREE (flush_buffer);

    iprint0 ("CPU cache flush complete.\n");
}