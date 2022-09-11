// ------------------------------------------------------------------
//   ref_lock.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "mutex.h"
#include "ref_private.h"
#include "buffer.h"
#include "flags.h"

#define GENOME_BASES_PER_MUTEX (1 << 16) // 2^16 = 64K - note: number of mutexes (per process) on Windows is limited by the total maximum HANDLEs - 16M, on Linux to 2G.

void ref_lock_initialize (Reference ref)
{
    buf_alloc_exact_zero (evb, ref->genome_muteces, (ref->genome_nbases + GENOME_BASES_PER_MUTEX-1) / GENOME_BASES_PER_MUTEX/*round up*/, Mutex, "genome_nbases");
    
    bool show = mutex_is_show ("genome_muteces"); // use --show-mutex=genome_muteces to see these
    if (show) buf_alloc (evb, &ref->genome_mutex_names, 24 * ref->genome_muteces.len, 0, char, 0, "genome_mutex_names");

    for_buf2 (Mutex, mutex_p, i, ref->genome_muteces) {
        char *name = BAFTc (ref->genome_mutex_names);
        if (show) 
            ref->genome_mutex_names.len += sprintf (name, "genome_muteces[%u]", i) + 1;
        
        mutex_initialize_do (mutex_p, show ? name : "genome_muteces[]", __FUNCTION__);
    }
}

void ref_lock_free (Reference ref)
{
    for_buf (Mutex, mutex_p, ref->genome_muteces) 
        mutex_destroy (*mutex_p);

    buf_free (ref->genome_muteces);
    buf_free (ref->genome_mutex_names);
}

// lock a region that includes the region given and the flanking regions 
RefLock ref_lock (Reference ref, PosType gpos_start, uint32_t seq_len)
{
    PosType last_gpos = MIN_(ref->genome_nbases, gpos_start + (PosType)seq_len - 1);
    gpos_start = MAX_(0, gpos_start);

    RefLock lock = { .first_mutex = gpos_start / GENOME_BASES_PER_MUTEX,
                     .last_mutex  = last_gpos  / GENOME_BASES_PER_MUTEX };

    ASSERT (lock.first_mutex >= 0 && lock.first_mutex <= ref->genome_muteces.len32, "lock.first_mutex=%u out of range: [0,%u]", lock.first_mutex, ref->genome_muteces.len32);
    ASSERT (lock.last_mutex  >= 0 && lock.last_mutex  <= ref->genome_muteces.len32, "lock.last_mutex=%u out of range: [0,%u]",  lock.last_mutex,  ref->genome_muteces.len32);
    
    // lock muteces in order
    for (int i=lock.first_mutex; i <= lock.last_mutex; i++)
        mutex_lock (*B(Mutex, ref->genome_muteces, i));

    return lock;
}

void ref_unlock (Reference ref, RefLock *lock)
{
    // unlock muteces in reverse order
    if (lock->first_mutex >= 0) { // first_mutex==-1 if REFLOCK_NONE
        for (int i=lock->last_mutex; i >= lock->first_mutex; i--)
            mutex_unlock (*B(Mutex, ref->genome_muteces, i));
    }

    *lock = REFLOCK_NONE;
}
