// ------------------------------------------------------------------
//   ref_cache.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#ifndef _WIN32
#include <fcntl.h>
#include <errno.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>
#endif
#include "genozip.h"
#include "mutex.h"
#include "ref_private.h"
#include "flags.h"
#include "file.h"
#include "refhash.h"
#include "arch.h"

bool ref_cache_is_cached  (Reference ref) { return ref->cache_state == CACHE_OK;      }
bool ref_cache_is_loading (Reference ref) { return ref->cache_state == CACHE_LOADING; }

void ref_cache_get_refhash (Reference ref, void **cache_refhash_data, uint64_t *cache_refhash_size, uint32_t *base_layer_bits)
{
    ASSERT (ref->cache_state == CACHE_OK, "Expecting ref=%s to be CACHE_OK", ref->filename);

    *cache_refhash_data = ref->cache_refhash_data;
    *cache_refhash_size = ref->cache->refhash_size;
    *base_layer_bits    = ref->cache->base_layer_bits;
}

// maps cache to reference if shm exists, or creates a new shm if it doesn't, or creates genome in memory if no shm
// returns true if ref_cache can be used
bool ref_cache_initialize_genome (Reference ref)
{
    if (flag.is_windows) goto error;

#ifndef _WIN32

    #define FAIL_MSG "FYI: Reference file cache is not used. Technical details: "
    ASSERT0 (ref->cache_state == CACHE_INITITAL, "ref_cache already initialized");
    ASSERT0 (flag.reading_reference && z_file && z_file->file, "not reading reference");

    key_t key = ftok (ref->filename, 20010802);

    struct stat st;
    ASSERTGOTO (stat64 (ref->filename, &st) >= 0, FAIL_MSG "stat (%s) failed: %s", ref->filename, strerror(errno));

    uint64_t genome_size  = roundup_bits2bytes64 (ref->genome_nbases * 2);
    
    uint8_t base_layer_bits;
    uint64_t refhash_size = refhash_load_init (&base_layer_bits);    

    int permissions = 0600 | (st.st_mode & 066); // RW permissions to "groups" and "other" copied from the reference file

    // note: a new shm segment is initialized by the OS to 0.
    ref->cache_shm = shmget (key, sizeof (struct ref_cache) + genome_size + refhash_size, IPC_CREAT | permissions); 
    ASSERTGOTO (ref->cache_shm >= 0, FAIL_MSG "shmget (%s) failed: %s", ref->filename, strerror(errno));

    ASSERTGOTO ((ref->cache = shmat (ref->cache_shm, NULL, 0)) != (void*)-1, FAIL_MSG "shmat (%s) failed: %s", ref->filename, strerror(errno));

    uint64_t now = arch_timestamp() / 1000000000; // convert nanosec to seconds
    uint64_t expected = 0;
    bool we_load = __atomic_compare_exchange_n (&ref->cache->creation_ts, &expected, now, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE);

    // if someone else has loaded or is loading now, we busy-wait for is_populated, up to 1 minute.
    if (!we_load) {
        while (!__atomic_load_n (&ref->cache->is_populated, __ATOMIC_ACQUIRE) &&
               (arch_timestamp() / 1000000000) - ref->cache->creation_ts < 60) 
            usleep (250000); // 1/4 sec

        ASSERTGOTO (ref->cache->is_populated, FAIL_MSG " another genozip process is populating the cache for the reference file %s "
            "but it appears that that process has stalled. To resume using the cache, manually remove the old shared memory with: "
            "\"ipcrm -m %u\"", ref->filename, ref->cache_shm);

        ref->cache_state = CACHE_OK;
    }

    // our process needs to populate this shm
    else {
        // we can't populate prim_ref, because we only have one refhash data structure, which represents gref data
        if (ref != gref) {
            __atomic_store_n (&ref->cache->creation_ts, (uint64_t)0, __ATOMIC_RELEASE);
            shmdt (ref->cache);
            ref->cache = NULL;
            goto error;
        }

        // note: mutex remains locked until loading is done
        ref->cache_state            = CACHE_LOADING;

        // data set in shm for all processes to consume. immutable once set.
        ref->cache->genome_size     = genome_size;
        ref->cache->refhash_size    = refhash_size;
        ref->cache->base_layer_bits = base_layer_bits;
        ref->cache->ref_genozip_ver = z_file->genozip_version;
    }

    // set data private to our process
    buf_attach_to_shm (evb, &ref->genome_buf, ref->cache->genome_data, genome_size, "genome_buf");
    ref->genome         = (BitsP)&ref->genome_buf;
    ref->genome->nbits  = ref->genome_nbases * 2;   
    ref->genome->nwords = roundup_bits2words64 (ref->genome->nbits);
    ref->cache_refhash_data = ref->cache->genome_data + genome_size;
#endif

    return true;

error:
    ref->cache_shm = 0;
    ref->cache_state = CACHE_NONE;
    return false;
}

void ref_cache_done_loading (Reference ref)
{
    __atomic_store_n (&ref->cache->is_populated, (bool)true, __ATOMIC_RELEASE); // immutable once set
    ref->cache_state = CACHE_OK;
}

// mark cache shm for removal. it will be removed when n_attached nattch decrements to 0.
void ref_cache_remove (Reference ref)
{
#ifndef _WIN32    
    flag.no_cache = false; // so that the cache loads
    
    static Context chrom_ctx = {}; // static because the contained buffers cannot be on the stack as they are added to buf_list
    ref_load_external_reference (ref, &chrom_ctx);

    if (ref->cache_state == CACHE_OK) {
        ASSERT (shmctl (ref->cache_shm, IPC_RMID, NULL) != -1,
                "Failed to remove cache of reference file %s: %s (try: ipcrm -m %u)", 
                ref->filename, strerror(errno), ref->cache_shm);

        iprintf ("Removing cache of reference file %s from shared memory.\n"
                 "If there are running genozip processes using this reference, it will be removed when the last process exits.\n",
                 ref->filename);

        exit_ok(); 
    }
#else
    ABORTINP0 ("reference file caching is currently not supported on Windows");
#endif
}

// note: buffers attached to this cache become invalid after shm is detached. best to free/destroy them before.
void ref_cache_detach (Reference ref)
{
#ifndef _WIN32 
    if (!ref || ref->cache_state != CACHE_OK) return;

    shmdt (ref->cache);
    ref->cache = NULL;
    ref->cache_state = CACHE_INITITAL;
#endif
}