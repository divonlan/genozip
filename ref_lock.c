// ------------------------------------------------------------------
//   ref_lock.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "mutex.h"
#include "ref_private.h"
#include "buffer.h"
#include "flags.h"

#define GENOME_BASES_PER_MUTEX (1 << 16) // 2^16 = 64K

static void ref_lock_initialize_do (Reference ref, uint32_t num_muteces)
{
    #define GM_NAME "genome_muteces[%u]"
    #define GM_NAME_LEN 24
    ref->genome_num_muteces = num_muteces;
    ref->genome_muteces = CALLOC (num_muteces * sizeof (Mutex));
    ref->genome_mutex_names = MALLOC (GM_NAME_LEN * num_muteces);

    bool create_names = mutex_is_show (GM_NAME);

    for (unsigned i=0; i < num_muteces; i++) {
        if (create_names) sprintf (&ref->genome_mutex_names[i*GM_NAME_LEN], GM_NAME, i);
        mutex_initialize_do (&ref->genome_muteces[i], create_names ? &ref->genome_mutex_names[i*GM_NAME_LEN] : GM_NAME, __FUNCTION__);
    }
}

void ref_lock_initialize_loaded_genome (Reference ref) 
{ 
    ref_lock_initialize_do (ref, (ref->genome_nbases + GENOME_BASES_PER_MUTEX-1) / GENOME_BASES_PER_MUTEX); // round up
}

void ref_lock_initialize_denovo_genome (Reference ref) 
{ 
    ref_lock_initialize_do (ref, REF_NUM_DENOVO_RANGES); 
}

void ref_lock_free (Reference ref)
{
    if (ref->genome_muteces) {
        for (unsigned i=0; i < ref->genome_num_muteces; i++)
            mutex_destroy (ref->genome_muteces[i]);

        FREE (ref->genome_muteces);
        ref->genome_num_muteces = 0;

        FREE (ref->genome_mutex_names);
    }
}

// lock a region that includes the region given and the flanking regions 
RefLock ref_lock (Reference ref, PosType gpos_start, uint32_t seq_len)
{
    PosType last_pos = gpos_start + seq_len - 1;

    // round to 64 before and after the request region
    RefLock lock = { .first_mutex = gpos_start / GENOME_BASES_PER_MUTEX,
                     .last_mutex  = last_pos   / GENOME_BASES_PER_MUTEX };

    ASSERT (lock.first_mutex >= 0 && lock.first_mutex <= ref->genome_num_muteces, "lock.first_mutex=%u out of range: [0,%u]", lock.first_mutex, ref->genome_num_muteces);
    ASSERT (lock.last_mutex  >= 0 && lock.last_mutex  <= ref->genome_num_muteces, "lock.last_mutex=%u out of range: [0,%u]", lock.last_mutex, ref->genome_num_muteces);
    
    // lock muteces in order
    for (int i=lock.first_mutex; i <= lock.last_mutex; i++)
        mutex_lock (ref->genome_muteces[i]);

    return lock;
}

RefLock ref_unlock (Reference ref, RefLock lock)
{
    if (lock.first_mutex >= 0) {
        // unlock muteces in reverse order
        for (int i=lock.last_mutex; i >= lock.first_mutex; i--)
            mutex_unlock (ref->genome_muteces[i]);
    }

    return REFLOCK_NONE;
}

// used for RT_DENOVO - single mutex er range
RefLock ref_lock_range (Reference ref, int32_t range_id)
{
    RefLock lock = { .first_mutex = range_id, .last_mutex = range_id };
    mutex_lock (ref->genome_muteces[range_id]);

    return lock;
}
