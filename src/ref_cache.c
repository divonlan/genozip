// ------------------------------------------------------------------
//   ref_cache.c
//   Copyright (C) 2023-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <fcntl.h>
#include <errno.h>
#include <pwd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>
#endif
#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/shared_region.h>
#endif
#include "ref_private.h"
#include "file.h"
#include "refhash_friend.h"
#include "arch.h"
#include "filename.h"

bool ref_cache_is_cached     (void) { return gref.cache_state == CACHE_READY;      }
bool ref_cache_is_populating (void) { return gref.cache_state == CACHE_POPULATING; }

#define ALOAD(cache_field) load_acquire (gref.cache->cache_field)
#define ASTORE(cache_field, value) store_release (gref.cache->cache_field, (value))

#define NO_SHM ((void *)-1)

#ifdef _WIN32
#define CACHE_ITERATOR_CB(func) bool func (uint32_t ref_i, int pid)

// messages genozip -> holder process (Windows)
#define MSG_LIST      "List" 
#define MSG_TERMINATE "Terminate"

#else
#define CACHE_ITERATOR_CB(func) bool func (uint32_t ref_i, int shmid, RefCache *cache)
#endif
typedef CACHE_ITERATOR_CB ((*RefCacheIteratorCallback));


static RefCacheState ref_cache_set_ready (void)
{
    void *old_attachment = gref.cache;

    // re-attach as read-only first attach new, then detach old, to prevent cache from being deleted if marked for removal
#ifndef _WIN32  
    gref.cache = shmat (gref.cache_shm, NULL, SHM_RDONLY); // sometimes fails in Mac, bug 1095
    if (gref.cache != NO_SHM) 
        ASSERT (!shmdt (old_attachment), "shmdt failed: %s", strerror (errno));
    else
        WARN (_WRN "shmat (read-only) failed: %s. shm remains RW. No harm.", strerror (errno)); 

#else
    gref.cache = MapViewOfFile (gref.cache_shm, FILE_MAP_READ, 0, 0, 0);
    ASSERT (gref.cache, "MapViewOfFile (read-only) failed: %s", str_win_error()); 

    ASSERT (UnmapViewOfFile (old_attachment), "UnmapViewOfFile failed: %s", str_win_error());
#endif

    if (flag.show_cache) iprintf ("show-cache: cache of %s: attached read-only + detached read-write. READY.\n", gref.filename);

    // for data integrity: the only place we set the state to CACHE_READY is here, after attaching as read-only
    return (gref.cache_state = CACHE_READY); 
}

static RefCacheState ref_cache_handle_existing (uint64_t data_size, 
                                                uint64_t original_creation_ts, uint64_t now)
{    
    // note: this is the pid of the original creator, or a third process that grabbed the pid in the mean time
    uint32_t creator_pid = ALOAD (creator_pid);

    // case: original creator died without completing populating that cache 
    if (!arch_is_process_alive (creator_pid) && !ALOAD (is_populated)) {    // the creator did not finish its work

        attempt_grab:
        // case: we successfully beat any other process to setting creaton_ts, so we can now populate the cache.
        if (__atomic_compare_exchange_n (&gref.cache->creation_ts, &original_creation_ts, now, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {

            // reset for safety (note: original shm allocation is reset by OS)
            memset (gref.cache->genome_data, 0, data_size);
            return CACHE_POPULATING;
        }

        // a third party beat us to changing the ts. that process is now populating so we shall wait for it
        else {}
    }

    // either the original creator process or a third process that beat us to grabbing the 
    // timestamp are alive. We shall wait up to 1 minute for the cache to be populated
    for (int i=0; !ALOAD (is_populated) && (arch_timestamp() / 1000000000) - gref.cache->creation_ts < 60;   i++) {

        // process died without completing
        if (!arch_is_process_alive (creator_pid) && !ALOAD (is_populated)) 
            goto attempt_grab; 

        // message after 5 seconds
        if (i==5*4) iprint0 ("Waiting for another Genozip process to populate the reference cache. Going to wait for a minute before giving up.\n"); 
        usleep (250000); // 1/4 sec
    }

    if (ALOAD (is_populated)) 
        return ref_cache_set_ready();

    ABORTINP ("Another process pid=%u is populating the cache for this reference file but it appears\n"
              "that that process has stalled. To resume using the cache, remove the abandoned cache from shared memory with:\n"
              "genocat%s --no-cache -e %s\n", creator_pid, flag.is_windows ? ".exe" : "", gref.filename);
    
    return 0; // bogus return
}

// maps cache to reference if shm exists, or creates a new shm if it doesn't, or creates genome in memory if no shm
// returns true if ref_cache can be used
bool ref_cache_initialize_genome (void)
{
    static rom tip = _TIP "run: 'genozip --no-cache' to clear cache";

    uint64_t refhash_size = refhash_get_refhash_size();  

    uint64_t genome_size  = roundup_bits2bytes64 (gref.genome_nbases * 2);

    if (gref.cache_state == CACHE_READY) goto cache_ok; // cache already loaded

    #define FAIL_MSG _FYI "Reference file cache is not used. Technical details: "
    ASSERT (gref.cache_state == CACHE_INITITAL, "unexpected cache_state=%d", gref.cache_state);
    ASSERT0 (flag.reading_reference && z_file && z_file->file, "not reading reference");

    uint64_t shm_size = sizeof (RefCache) + genome_size + refhash_size;
    uint32_t holder_pid = 0;
    
#ifndef _WIN32
    struct stat st;
    ASSGOTO (stat64 (gref.filename, &st) >= 0, FAIL_MSG "stat (%s) failed: %s", gref.filename, strerror(errno));

    // since 15.0.82: key by inode (up to 15.0.81: ftok (ref_basename, 20010802), ref_basename as in command line)
    // note: inode is consistent across docker containers mounting the same filesystem, 
    // and in case of a symlink, it is the inode of the actual target file
    key_t key = fibonacci (st.st_ino + 20010802/*salt*/, 31); // 31 to avoid negative keys, which cause mis-indentation in /proc/sysvipc/shm that we parse
    int permissions = 0600 | (st.st_mode & 066); // RW permissions to "groups" and "other" copied from the reference file

    // note: a new shm segment is initialized by the OS to 0.
    gref.cache_shm = shmget (key, shm_size, (flag.removing_cache ? 0 : IPC_CREAT) | permissions); 

    if (flag.is_mac && gref.cache_shm == -1 && (errno == EINVAL/*shmmax issue*/ || errno == ENOMEM/*shmall issue*/)) {
        WARN_ONCE (_FYI "Failed to cache the reference file in shared memory, because shared memory limits are too small. To increase limits temporarily until next reboot:\n====\n"
                   "genozip --no-cache\n"
                   "sudo sysctl -w kern.sysv.shmmax=%"PRIu64"\n"
                   "sudo sysctl -w kern.sysv.shmall=%"PRIu64"\n====\n"
                   "To increase limits permanently, modify /etc/sysctl.conf\n",
                   ROUNDUP1M (shm_size),                                  // max size of single shm (in bytes)
                   ROUNDUP1M (shm_size) / getpagesize() + 1/*arch shm*/); // max size of all shm (in pages).
        goto error;
    }

    bool cache_did_not_exist = (errno == ENOENT);

    ASSGOTO (gref.cache_shm >= 0 || (flag.removing_cache && cache_did_not_exist), 
             FAIL_MSG "shmget (%s key=0x%08x size=%"PRIu64") failed: %s.%s", 
             gref.filename, key, shm_size, strerror(errno), tip);

    if (flag.show_cache) iprintf ("show-cache: shmget of shm id %u\n", gref.cache_shm);

    if (gref.cache_shm != CACHE_SHM_NONE) {
        gref.cache = shmat (gref.cache_shm, NULL, 0);
        ASSGOTO (gref.cache != (void*)-1, FAIL_MSG "shmat (%s) failed: %s.%s", gref.filename, strerror(errno), tip);
    }

#else // Windows
    bool cache_did_not_exist;
    if (!flag.removing_cache) {
        gref.cache_shm = CreateFileMappingA (INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, shm_size >> 32, shm_size & 0xffffffff, gref.filename);
        cache_did_not_exist = (GetLastError() == ERROR_SUCCESS);
    
        ASSGOTO (gref.cache_shm, FAIL_MSG "CreateFileMapping (%s, size=%"PRIu64") failed: %s.%s", 
                 gref.filename, shm_size, str_win_error(), tip);
        if (flag.show_cache) iprintf ("show-cache: CreateFileMapping %s\n", gref.filename);
    }
    else { 
        gref.cache_shm = OpenFileMappingA (FILE_MAP_WRITE, false, gref.filename);
        cache_did_not_exist = !gref.cache_shm && GetLastError() == ERROR_FILE_NOT_FOUND;

        ASSERT (gref.cache_shm || cache_did_not_exist, FAIL_MSG "OpenFileMapping (%s) failed: %s.%s", 
                gref.filename, str_win_error(), tip);
        if (flag.show_cache) iprintf ("show-cache: OpenFileMapping %s: %s\n", gref.filename, str_win_error());
    }

    // if we just created new shm, also created the holder process
    if (cache_did_not_exist && !flag.removing_cache) {
        STARTUPINFO si = { .cb = sizeof (STARTUPINFO) };
        PROCESS_INFORMATION pi = {};
        rom exec = arch_get_genozip_executable().s;
        char cmd[strlen(exec) + 64];
        snprintf (cmd, sizeof (cmd), "%s --hold-cache=%"PRIu64, exec, (uint64_t)gref.cache_shm);
        
        ASSERT (SetHandleInformation (gref.cache_shm, HANDLE_FLAG_INHERIT, HANDLE_FLAG_INHERIT), "SetHandleInformation failed: %s", str_win_error());

        ASSERT (CreateProcessA (exec, cmd, NULL, NULL, true, 
                                BELOW_NORMAL_PRIORITY_CLASS |  // also used by ref_cache_remove_all to idetify the holder processes
                                CREATE_NEW_PROCESS_GROUP |     // without this, some shells will kill ALL cache holder processes created in this shell, if a genozip process is stopped with Ctrl-C
                                DETACHED_PROCESS | CREATE_BREAKAWAY_FROM_JOB, // Running a Windows exe on WSL2: this prevents Linux from killing the Windows cache holder process when the main genozip.exe process finishes
                                NULL, NULL, &si, &pi),
                "Failed to create cache holder process for %s: \"%s\": %s", gref.filename, cmd, str_win_error());
        
        CloseHandle (pi.hThread); // child process's main thread
        holder_pid = pi.dwProcessId;
    }

    if (gref.cache_shm != CACHE_SHM_NONE) {
        gref.cache = MapViewOfFile (gref.cache_shm, FILE_MAP_WRITE, 0, 0, 0);
        ASSGOTO (gref.cache, FAIL_MSG "MapViewOfFile (read-write) (%s) failed: %s.%s", gref.filename, str_win_error(), tip);
    }

#endif
    
    if (flag.removing_cache) {
        ref_cache_remove_do (!cache_did_not_exist, true); 
        exit_ok;
    }

    if (flag.show_cache) iprint0 ("show-cache: shmat read-write\n");

    uint64_t now = arch_timestamp() / 1000000000; // convert nanosec to seconds
    uint64_t expected = 0;

    // save the ts - this may change if creating process crashes and a third process grabs the cache
    uint64_t original_creation_ts = ALOAD (creation_ts); 

    // grab (=set) the creation_ts field - that would indicate that we are the creators
    if (!original_creation_ts &&
        __atomic_compare_exchange_n (&gref.cache->creation_ts, &expected, now, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE))

        gref.cache_state = CACHE_POPULATING;

    else
        // returns CACHE_READY, CACHE_POPULATING or aborts
        gref.cache_state = ref_cache_handle_existing (genome_size + refhash_size, original_creation_ts, now);

    // our process needs to populate this shm
    if (gref.cache_state == CACHE_POPULATING) {
        // set our other data *after* the creation_ts
        gref.cache->creator_pid = getpid();
        gref.cache->magic = GENOZIP_MAGIC;
        gref.cache->genozip_version = code_version().major;
        gref.cache->RefCache_ver = 0;
        gref.cache->shm_size = shm_size;
        gref.cache->holder_pid = holder_pid; // 15.0.82
        gref.cache->unused = 0;
        filename_base (gref.filename, false, "<unknown>", gref.cache->ref_basename, sizeof (gref.cache->ref_basename));

        __atomic_thread_fence (__ATOMIC_RELEASE); 

        if (flag.show_cache) iprint0 ("show-cache: POPULATING\n");
    }

    gref.genome = (BitsP)&gref.genome_buf;

cache_ok:
    buf_attach_bits_to_shm (evb, &gref.genome_buf, gref.cache->genome_data, gref.genome_nbases * 2, "genome_buf");
    if (flag.show_cache) iprintf ("show-cache: attached genome_buf (%"PRIu64" bases) to %s shm\n", gref.genome_nbases, gref.cache_state == CACHE_READY ? "READONLY" : "READWRITE");

    // attach the cache data (read-write if CACHE_POPULATING - we will switch to read-only in ref_cache_done_populating; read-only if CACHE_READY)
    if (refhash_exists()) {
        buf_attach_to_shm (evb, &refhash_buf, gref.cache->genome_data + genome_size, refhash_size, "refhash_buf");
        refhash_buf.len = refhash_size;
        
        if (flag.show_cache) iprintf ("show-cache: attached refhash_buf (len=%"PRIu64") to %s shm\n", refhash_buf.len, gref.cache_state == CACHE_READY ? "READONLY" : "READWRITE");
    }
    
    return true;

error:
    gref.cache_shm = 0;
    gref.cache_state = CACHE_NONE;
    flag.no_cache = true;

    return false;
}

void ref_cache_done_populating (void)
{
    ASTORE (is_populated, (bool)true); // immutable once set
    ASTORE (creator_pid, (uint32_t)0); 
    if (flag.show_cache) iprint0 ("show-cache: done populating shm\n");

    ref_cache_set_ready();
}

// mark cache shm for removal. it will be removed when n_attached nattch decrements to 0.
void ref_cache_remove (void)
{
    flag.no_cache = false; // so that the cache loads
    flag.removing_cache = true;
    
    static Context chrom_ctx = {}; // static because the contained buffers cannot be on the stack as they are added to buf_list
    ref_load_external_reference (&chrom_ctx); // this will call ref_cache_remove_do
}

unsigned ref_cache_iterator (RefCacheIteratorCallback (callback), bool dormant_only)
{
    unsigned count = 0;

#ifdef _WIN32
    // find and terminate all holder process
    DWORD pids[16384], size; //  a large number
    ASSERT (EnumProcesses (pids, sizeof (pids), &size), "EnumProcesses failed: %s", str_win_error());

    for (int pid_i=0; pid_i < size / sizeof(DWORD); pid_i++) {
        HANDLE process = OpenProcess (PROCESS_QUERY_INFORMATION | PROCESS_TERMINATE, false, pids[pid_i]);
        if (!process) continue;

        // we use the priority to differentiate holder processes from normal genozip processes
        if (GetPriorityClass(process) != BELOW_NORMAL_PRIORITY_CLASS) {
            CloseHandle (process);
            continue;
        }

        char filename[MAX_PATH+1];
        uint32_t filename_len = GetProcessImageFileNameA (process, filename, MAX_PATH);
        filename[filename_len] = 0;
        CloseHandle (process);

        if (strstr (filename, "genozip") || strstr (filename, "genounzip") || strstr (filename, "genocat")) 
            count += callback (count, pids[pid_i]);
    }

#else

#ifdef __linux__
    ASSERTNOTINUSE (evb->scratch);
    file_get_file (evb, "/proc/sysvipc/shm", &evb->scratch, "scratch", 1 MB, VERIFY_ASCII, false);

    str_split_by_lines (evb->scratch.data, evb->scratch.len, 1000);

    for (int i=1; i < n_lines; i++) { // note: skipping first line - its a header
        if (line_lens[i] < 20) continue; 

        char *shmid_str = (char *)&lines[i][11];
        uint32_t shmid_str_len = 10;
        shmid_str[10] = 0;
        str_trim (qSTRa(shmid_str));

        int shmid = atoi (shmid_str);

#else // Mac and Unixes
    // brute-force scan for valid shm_id's... surely there is a better way
    for (int shmid = 0; shmid <= 1 MB; shmid++) { 
#endif // mac & unix

        // filter out non-dormant if needed. note: must test for dormancy before attaching!
        if (dormant_only) { 
            struct shmid_ds ds = {};
            if (shmctl (shmid, IPC_STAT, &ds) == -1) continue; // fail silently

            time_t ago = time(0) - MAX_(ds.shm_atime, ds.shm_dtime);

            if (ago < 60 * 60 * 24) // not dormant (attached or detached in the past 24 hours)
                continue;
        }

        RefCache *cache = shmat (shmid, NULL, SHM_RDONLY);
        if ((cache != NO_SHM) && (cache->magic == GENOZIP_MAGIC)) 
            count += callback (count, shmid, cache);
        
        shmdt (cache);
    }
#endif // not windows

    return count;
}

#ifdef _WIN32
StrText ref_cache_win_pipe_name (int pid)
{
    if (!pid) // server side
        pid = GetProcessId (GetCurrentProcess());

    StrText s;
    snprintf (s.s, sizeof(s.s), "\\\\.\\pipe\\genozip_%d", pid);
    return s;
}

static bool ref_cache_msg_holder_process (int pid, rom request, StrText1K *response)
{

    // 2. Open the pipe as a file
    HANDLE named_pipe = CreateFileA (
        ref_cache_win_pipe_name (pid).s,
        GENERIC_READ | GENERIC_WRITE,
        0,              // No sharing
        NULL,           // Default security
        OPEN_EXISTING,  // Open existing pipe
        0,              // Default attributes
        NULL);          // No template file

    ASSRET(named_pipe != INVALID_HANDLE_VALUE, false, 
           "opening pipe to cache holder failed: %s", str_win_error());

    DWORD bytes;
    ASSRET (WriteFile(named_pipe, request, strlen (request)+1, &bytes, NULL), false,
           "sending request to cache holder failed: %s", str_win_error());
           
    ASSRET (ReadFile(named_pipe, response->s, sizeof(response->s), &bytes, NULL), false,
           "reading response from cache holder failed: %s", str_win_error());

    CloseHandle (named_pipe); // pipe leaked in case of failure
    return true;
}
#endif


static CACHE_ITERATOR_CB(do_remove)
{
    rom name;

#ifdef _WIN32
    StrText1K response;
    bool success = ref_cache_msg_holder_process (pid, MSG_TERMINATE, &response);
    name = response.s;
#else
    bool success = (shmctl (shmid, IPC_RMID, NULL) != -1); 
    ASSERTW (success, _WRN "shmctl failed: %s", strerror (errno));

    name = cache->ref_basename;
#endif

    if (IS_RM_CACHE || flag.show_cache)
        WARN (_FYI "Unloaded reference cache \"%s\"", name);

    return success;
}

void ref_cache_remove_do (bool cache_exists, bool verbose)
{
    ASSINP (cache_exists, "There is currently no cache for reference file %s", gref.filename);
    if (!gref.cache) return; // fail silently (could happen with --no-cache)

#ifndef _WIN32    
    do_remove (0, gref.cache_shm, gref.cache);
#else // Windows  
    do_remove (0, gref.cache->holder_pid);    
#endif

    if (verbose) 
        iprintf ("Removing cache of reference file %s from shared memory.\n"
                 "If there are running genozip processes using this reference, it will be removed when the last process exits.\n",
                 gref.filename);
}

// remove all Genozip cache shm segments 
void ref_cache_remove_all (RefCacheRemoveType rm_type)
{
    if (rm_type == REF_CACHE_REMOVE_DORMANT && flag.is_windows) return; // not supported for Windows
    
    unsigned n_removed = ref_cache_iterator (do_remove, rm_type == REF_CACHE_REMOVE_DORMANT);
    
    WARN_IF (!n_removed && rm_type == REF_CACHE_REMOVE_ALL, "No in-memory cached reference files found", NULL);
}

static CACHE_ITERATOR_CB (do_list)
{
#ifdef _WIN32
    if (ref_i == 0) // first reference
        iprint0 ("pid    size         loaded               name\n");

    StrText1K response;
    if (!ref_cache_msg_holder_process (pid, MSG_LIST, &response))
        return false;

    iprintf ("%s", response.s);
    return true;
    
#else
    if (ref_i == 0) // first reference
        iprint0 ("shmid  owner     perm size         loaded               name\n");

    struct shmid_ds ds = {};
    ASSERT (shmctl (shmid, IPC_STAT, &ds) != -1, "shmctl failed: %s", strerror (errno));

    StrText time_str = {};
    struct tm *time_info = localtime ((time_t *)&cache->creation_ts); // assumes time_t is 64 bit
    strftime (time_str.s, sizeof(time_str)-1, "%Y-%m-%d %H:%M:%S", time_info);
    iprintf ("%-5u  %-8s  %03o  %-11"PRIu64"  %19s  %s  %s\n", 
             shmid, getpwuid (ds.shm_perm.uid)->pw_name, ds.shm_perm.mode & 0777, cache->shm_size,
             time_str.s, cache->ref_basename, (cache->is_populated ? "" : " NOT READY"));
    return true;
#endif
}

void ref_cache_ls (void)
{
    if (!ref_cache_iterator (do_list, false))
        WARN (_FYI "No in-memory cached reference files found", NULL);
}

// note: buffers attached to this cache become invalid after shm is detached. best to free/destroy them before.
void ref_cache_detach (void)
{
    if (gref.cache_state != CACHE_READY ||
        gref.genome_buf.type == BUF_SHM || refhash_buf.type == BUF_SHM) // actually detach only when both genome_buf and refhash are freed
        return;

#ifndef _WIN32 
    shmdt (gref.cache);
#else
    UnmapViewOfFile (gref.cache);
    CloseHandle (gref.cache_shm);
#endif

    gref.cache = NULL;
    gref.cache_shm = 0;
    gref.cache_state = CACHE_INITITAL;

    if (flag.show_cache) iprint0 ("show-cache: detached shm\n");
}

rom cache_state_name (RefCacheState cs)
{
    return IN_RANGE (cs, 0, NUM_CACHE_STATES) ? (rom[])CACHE_STATE_NAMES[cs] : "InvalidRefCacheState";
}

// Windows-only, run by cache-hold process: since it inherited the file mapping handle, the file mapping
// will remain in memory as long as this process is alive. It dies when it is asked to.
void noreturn ref_cache_hold (rom handle_str)
{
#ifdef _WIN32
    HANDLE h = (HANDLE)atoll (handle_str); // inherited from cache creator process

    RefCache *cache = MapViewOfFile (h, FILE_MAP_READ, 0, 0, sizeof (RefCache));

    // wait indefinitely for requests from user
    while (1) {
        #define REQUEST_SIZE 64/*bytes*/
        HANDLE pipe_handle = CreateNamedPipeA ( // a new instance of the pipe
            ref_cache_win_pipe_name(0).s, PIPE_ACCESS_DUPLEX,
            PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_WAIT,
            PIPE_UNLIMITED_INSTANCES, sizeof (StrText1K), REQUEST_SIZE, 0, NULL);
        
        if (pipe_handle == INVALID_HANDLE_VALUE) stall(); // can't receive messages, but is held

        // wait for a client to connect to this instance
        if (ConnectNamedPipe (pipe_handle, NULL) || GetLastError() == ERROR_PIPE_CONNECTED) {
            char req[REQUEST_SIZE] = {};
            DWORD bytes;
            bool terminate = false;

            if (ReadFile (pipe_handle, req, REQUEST_SIZE-1, &bytes, NULL) && bytes) {
                StrText1K response = {};
                
                if (!strcmp (req, MSG_TERMINATE)) {
                    terminate = true;
                    strcpy (response.s, cache->ref_basename);
                }

                else if (!strcmp (req, MSG_LIST)) {
                    StrText time_str = {};
                    struct tm *time_info = localtime ((time_t *)&cache->creation_ts); // assumes time_t is 64 bit
                    strftime (time_str.s, sizeof(time_str)-1, "%Y-%m-%d %H:%M:%S", time_info);

                    snprintf (response.s, sizeof(response)-1, "%-5u  %-11"PRIu64"  %19s  %s  %s\n", 
                              cache->holder_pid, cache->shm_size, time_str.s, cache->ref_basename, (cache->is_populated ? "" : " NOT READY"));
                }

                else
                    strcpy (response.s, "Invalid request to cache holder process");

                WriteFile (pipe_handle, response.s, strlen (response.s)+1/*NUL*/, &bytes, NULL);

                FlushFileBuffers (pipe_handle);
                DisconnectNamedPipe (pipe_handle);
                CloseHandle (pipe_handle);

                if (terminate) break;
            }
       }    
    }
#endif
    exit (0); // now that the holder process has exited, the shared memory will be destroyed by Windows when last user process terminates
}
