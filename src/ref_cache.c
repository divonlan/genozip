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
#include "refhash.h"
#include "arch.h"
#include "filename.h"

bool ref_cache_is_cached     (void) { return gref.cache_state == CACHE_READY;      }
bool ref_cache_is_populating (void) { return gref.cache_state == CACHE_POPULATING; }

#define ALOAD(cache_field) load_acquire (gref.cache->cache_field)
#define ASTORE(cache_field, value) store_release (gref.cache->cache_field, (value))

#define NO_SHM ((void *)-1)

void ref_cache_remove_do (bool cache_exists, bool verbose)
{
    ASSINP (cache_exists, "There is currently no cache for reference file %s", gref.filename);

#ifndef _WIN32    
    ASSERT (shmctl (gref.cache_shm, IPC_RMID, NULL) != -1,
            "Failed to remove cache of reference file %s: %s (try: ipcrm -m %u)", 
            gref.filename, strerror(errno), gref.cache_shm);

#else // Windows
    gref.cache->terminate_holder = true; // Windows memory model ensures this will be immediately visible to the holder process

#endif

    if (flag.show_cache) iprint0 ("show-cache: shmctl removed shm.\n");

    if (verbose) 
        iprintf ("Removing cache of reference file %s from shared memory.\n"
                 "If there are running genozip processes using this reference, it will be removed when the last process exits.\n",
                 gref.filename);

}

static RefCacheState ref_cache_set_ready (void)
{
    void *old_attachment = gref.cache;

    // re-attach as read-only first attach new, then detach old, to prevent cache from being deleted if marked for removal
#ifndef _WIN32  
    gref.cache = shmat (gref.cache_shm, NULL, SHM_RDONLY); // sometimes fails in Mac, bug 1095
    if (gref.cache != NO_SHM) 
        ASSERT (!shmdt (old_attachment), "shmdt failed: %s", strerror (errno));
    else
        WARN ("shmat (read-only) failed: %s. shm remains RW. No harm.", strerror (errno)); 

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

// Windows-only, run by cache-hold process: since it inherited the file mapping handle, the file mapping
// will remain in memory as long as this process is alive. It dies when it is asked to.
void noreturn ref_cache_hold (rom handle_str)
{
#ifdef _WIN32
    HANDLE h = (HANDLE)atoll (handle_str); // inherited from cache creator process

    const RefCache *ref_cache = MapViewOfFile (h, FILE_MAP_READ, 0, 0, sizeof (RefCache));
    
    while (!ref_cache->terminate_holder) 
        sleep(2); // it will take us up to 2 seconds get the termination order. that's fine - user is not waiting on this.
#endif
    exit (0);
}

// maps cache to reference if shm exists, or creates a new shm if it doesn't, or creates genome in memory if no shm
// returns true if ref_cache can be used
bool ref_cache_initialize_genome (void)
{
    uint64_t refhash_size = refhash_get_refhash_size();    
    uint64_t genome_size  = roundup_bits2bytes64 (gref.genome_nbases * 2);

    if (gref.cache_state == CACHE_READY) goto cache_ok; // cache already loaded

    #define FAIL_MSG "FYI: Reference file cache is not used. Technical details: "
    ASSERT (gref.cache_state == CACHE_INITITAL, "unexpected cache_state=%d", gref.cache_state);
    ASSERT0 (flag.reading_reference && z_file && z_file->file, "not reading reference");

    uint64_t shm_size = sizeof (RefCache) + genome_size + refhash_size;
    
#ifndef _WIN32
    key_t key = ftok (gref.filename, 20010802);

    struct stat st;
    ASSGOTO (stat64 (gref.filename, &st) >= 0, FAIL_MSG "stat (%s) failed: %s", gref.filename, strerror(errno));
    
    int permissions = 0600 | (st.st_mode & 066); // RW permissions to "groups" and "other" copied from the reference file

    // note: a new shm segment is initialized by the OS to 0.
    gref.cache_shm = shmget (key, shm_size, (flag.removing_cache ? 0 : IPC_CREAT) | permissions); 

    if (flag.is_mac && gref.cache_shm == -1 && (errno == EINVAL/*shmmax issue*/ || errno == ENOMEM/*shmall issue*/)) {
        WARN_ONCE ("FYI: Failed to cache the reference file in shared memory, because shared memory limits are too small. To increase limits temporarily until next reboot:\n====\n"
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
                FAIL_MSG "shmget (%s) failed: %s", gref.filename, strerror(errno));

    if (flag.show_cache) iprintf ("show-cache: shmget of shm id %u\n", gref.cache_shm);

    if (gref.cache_shm != CACHE_SHM_NONE) {
        gref.cache = shmat (gref.cache_shm, NULL, 0);
        ASSGOTO (gref.cache != (void*)-1, FAIL_MSG "shmat (%s) failed: %s", gref.filename, strerror(errno));
    }

#else // Windows
    bool cache_did_not_exist;
    if (!flag.removing_cache) {
        gref.cache_shm = CreateFileMappingA (INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, shm_size >> 32, shm_size & 0xffffffff, gref.filename);
        cache_did_not_exist = (GetLastError() == ERROR_SUCCESS);
    
        ASSGOTO (gref.cache_shm, FAIL_MSG "CreateFileMapping (%s) failed: %s", gref.filename, str_win_error());
        if (flag.show_cache) iprintf ("show-cache: CreateFileMapping %s\n", gref.filename);
    }
    else { 
        gref.cache_shm = OpenFileMappingA (FILE_MAP_WRITE, false, gref.filename);
        cache_did_not_exist = !gref.cache_shm && GetLastError() == ERROR_FILE_NOT_FOUND;

        ASSERT (gref.cache_shm || cache_did_not_exist, FAIL_MSG "OpenFileMapping (%s) failed: %s", gref.filename, str_win_error());
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

        // Limitation: when running a *Windows* genozip.exe executable from a WSL session, the holder process
        // will die when the parent process dies.
        ASSERT (CreateProcessA (exec, cmd, NULL, NULL, true, 
                                BELOW_NORMAL_PRIORITY_CLASS |  // also used by ref_cache_remove_all to idetify the holder processes
                                CREATE_NEW_PROCESS_GROUP,      // without this, some shells will kill ALL cache holder processes created in this shell, if a genozip process is stopped with Ctrl-C
                                NULL, NULL, &si, &pi),
                "Failed to create cache holder process for %s: \"%s\": %s", gref.filename, cmd, str_win_error());
        
        CloseHandle (pi.hThread); // child process's main thread
    }

    if (gref.cache_shm != CACHE_SHM_NONE) {
        gref.cache = MapViewOfFile (gref.cache_shm, FILE_MAP_WRITE, 0, 0, 0);
        ASSGOTO (gref.cache, FAIL_MSG "MapViewOfFile (read-write) (%s) failed: %s", gref.filename, str_win_error());
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
        gref.cache->genozip_version = code_version_major();
        gref.cache->shm_size = shm_size;
        filename_base (gref.filename, false, "<unknown>", gref.cache->ref_basename, sizeof (gref.cache->ref_basename));

        __atomic_thread_fence (__ATOMIC_RELEASE); 

        if (flag.show_cache) iprint0 ("show-cache: POPULATING\n");
    }

    gref.genome = (BitsP)&gref.genome_buf;

cache_ok:
    buf_attach_bits_to_shm (evb, &gref.genome_buf, gref.cache->genome_data, gref.genome_nbases * 2, "genome_buf");
    if (flag.show_cache) iprintf ("show-cache: attached genome_buf (%"PRIu64" bases) to %s shm\n", gref.genome_nbases, gref.cache_state == CACHE_READY ? "READONLY" : "READWRITE");

    // attach the cache data (read-write if CACHE_POPULATING - we will switch to read-only in ref_cache_done_populating; read-only if CACHE_READY)
    buf_attach_to_shm (evb, &refhash_buf, gref.cache->genome_data + genome_size, refhash_size, "refhash_buf");
    refhash_buf.len = refhash_size;
    
    if (flag.show_cache) iprintf ("show-cache: attached refhash_buf (len=%"PRIu64") to %s shm\n", refhash_buf.len, gref.cache_state == CACHE_READY ? "READONLY" : "READWRITE");
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

unsigned ref_cache_iterator (bool (*callback)(int shmid, RefCache *cache))
{
    unsigned count = 0;

#ifdef __linux__
    ASSERTNOTINUSE (evb->scratch);
    file_get_file (evb, "/proc/sysvipc/shm", &evb->scratch, "scratch", 1 MB, VERIFY_ASCII, false);

    str_split_by_lines (evb->scratch.data, evb->scratch.len, 1000);

    for (int i=1; i < n_lines; i++) { // note: skipping first line - its a header
        if (line_lens[i] < 20) continue; 

        char *shmid_str = (char *)&lines[i][11];
        uint32_t shmid_str_len = 10;
        shmid_str[10] = 0;
        str_trim (qSTRa(shmid_str));;

        int shmid = atoi(shmid_str);

#else // Mac and Unixes
    // brute-force scan for valid shm_id's... surely there is a better way
    for (int shmid = 0; shmid <= 1 MB; shmid++) { 
#endif

#ifndef _WIN32
        RefCache *cache = shmat (shmid, NULL, SHM_RDONLY);
        if ((cache != NO_SHM) && (cache->magic == GENOZIP_MAGIC)) 
            count += callback (shmid, cache);
        
        shmdt (cache);
#endif
    }

    if (!count) WARN0 ("No in-memory cached reference files found");

    return count;
}

#ifndef _WIN32
static bool do_remove (int shmid, RefCache *cache)
{
    WARN ("Unloading reference cache \"%s\"", cache->ref_basename);
    bool success = (shmctl (shmid, IPC_RMID, NULL) != -1); 

    return success;
}
#endif

// remove all Genozip cache shm segments 
void ref_cache_remove_all (void)
{
#ifndef _WIN32
    ref_cache_iterator (do_remove);

#elif defined _WIN32
    // find and terminate all holder process
    DWORD pids[16384], size; //  a large number
    ASSERT (EnumProcesses (pids, sizeof (pids), &size), "EnumProcesses failed: %s", str_win_error());

    int count = 0;
    for (int pid_i=0; pid_i < size / sizeof(DWORD); pid_i++) {
        HANDLE process = OpenProcess (PROCESS_QUERY_INFORMATION | PROCESS_TERMINATE, false, pids[pid_i]);
        if (!process) continue;

        // we use the priority to differentiate holder processes from normal genozip processes
        if (GetPriorityClass(process) != BELOW_NORMAL_PRIORITY_CLASS) continue;

        char filename[MAX_PATH+1];
        uint32_t filename_len = GetProcessImageFileNameA (process, filename, MAX_PATH);
        filename[filename_len] = 0;

        if (strstr (filename, "genozip") || strstr (filename, "genounzip") || strstr (filename, "genocat")) {
            bool success = TerminateProcess (process, 0); 
            ASSERTW (success, "TerminateProcess of pid=%u exec=%s failed: %s", (uint32_t)pids[pid_i], filename, str_win_error());
            if (success) count++;
        }
    }

    if (count) WARN ("Unloading reference cache of %d reference file%s", STRfN(count));
#endif 
}

static bool do_ls (int shmid, RefCache *cache)
{
    iprintf ("%s (shmid=%u size=%"PRIu64")%s\n", cache->ref_basename, shmid,cache->shm_size, (cache->is_populated ? "" : " NOT READY"));
    
    return true;
}

void ref_cache_ls (void)
{
    ASSINP0 (!flag.is_windows, "genols --cache option is not available on Windows");

    ref_cache_iterator (do_ls);
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
