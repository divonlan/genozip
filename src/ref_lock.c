// ------------------------------------------------------------------
//   ref_lock.c
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "mutex.h"
#include "ref_private.h"
#include "buffer.h"
#include "flags.h"

#define GENOME_BASES_PER_MUTEX (64 KB) // note: number of mutexes (per process) on Windows is limited by the total maximum HANDLEs - 16M, on Linux to 2G.

void ref_lock_initialize (void)
{
    buf_alloc_exact_zero (evb, gref.genome_muteces, (gref.genome_nbases + GENOME_BASES_PER_MUTEX-1) / GENOME_BASES_PER_MUTEX/*round up*/, Mutex, "genome_nbases");
    
    bool show = mutex_is_show ("genome_muteces"); // use --show-mutex=genome_muteces to see these
    if (show) buf_alloc (evb, &gref.genome_mutex_names, 32 * gref.genome_muteces.len, 0, char, 0, "genome_mutex_names");

    for_buf2 (Mutex, mutex_p, i, gref.genome_muteces) {
        char *name = BAFTc (gref.genome_mutex_names);
        
        if (show && name) 
            gref.genome_mutex_names.len += snprintf (name, gref.genome_mutex_names.size - gref.genome_mutex_names.len, "genome_muteces[%u]", i) + 1;

        mutex_initialize_do (mutex_p, show ? name : "genome_muteces[]", __FUNCTION__);
    }
}

void ref_lock_free (void)
{
    for_buf (Mutex, mutex_p, gref.genome_muteces) 
        mutex_destroy (*mutex_p);

    buf_free (gref.genome_muteces);
    buf_free (gref.genome_mutex_names);
}

// lock a region that includes the region given and the flanking regions 
RefLock ref_lock (PosType64 gpos_start, uint32_t seq_len)
{
    PosType64 last_gpos = MIN_(gref.genome_nbases, gpos_start + (PosType64)seq_len - 1);
    gpos_start = MAX_(0, gpos_start);

    RefLock lock = { .first_mutex = gpos_start / GENOME_BASES_PER_MUTEX,
                     .last_mutex  = last_gpos  / GENOME_BASES_PER_MUTEX };

    ASSERT (lock.first_mutex >= 0 && lock.first_mutex <= gref.genome_muteces.len32, "lock.first_mutex=%u ∉ [0,%u]", lock.first_mutex, gref.genome_muteces.len32);
    ASSERT (lock.last_mutex  >= 0 && lock.last_mutex  <= gref.genome_muteces.len32, "lock.last_mutex=%u ∉ [0,%u]",  lock.last_mutex,  gref.genome_muteces.len32);
    
    // this happens when locking a region of over 4 Mb e.g. in SAM with cigar = 1M100000000N1M
    // when building with sanitize_thread (for debugging), we won't lock in this case. 
    if (flag.is_sanitize_thread && (lock.last_mutex - lock.first_mutex + 1) > 64) {
        fprintf (stderr, "\nFYI: ref_lock: not locking range because I am asked to lock %u muteces, but sanitize_thread only supports 64 concurrent muteces. In rare cases, this might cause a corrupt file.\n",
                 lock.last_mutex - lock.first_mutex + 1);
        return REFLOCK_NONE;
    }

    // lock muteces in order
    for (int i=lock.first_mutex; i <= lock.last_mutex; i++) 
        mutex_lock (*B(Mutex, gref.genome_muteces, i));

    return lock;
}

void ref_unlock (RefLock *lock)
{
    // unlock muteces in reverse order
    if (lock->first_mutex >= 0) { // first_mutex==-1 if REFLOCK_NONE
        for (int i=lock->last_mutex; i >= lock->first_mutex; i--)
            mutex_unlock (*B(Mutex, gref.genome_muteces, i));
    }

    *lock = REFLOCK_NONE;
}
