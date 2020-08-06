// ------------------------------------------------------------------
//   ref_lock.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "arch.h"
#include "ref_private.h"

#define GENOME_BASES_PER_MUTEX (1 << 16) // 2^16 = 64K

static pthread_mutex_t *genome_muteces = NULL; // one spinlock per 16K bases - protects genome->is_set
static uint32_t genome_num_muteces=0;

static void ref_lock_initialize_do (uint32_t num_muteces)
{
    genome_num_muteces = num_muteces;
    genome_muteces = calloc (num_muteces, sizeof (pthread_mutex_t));
    for (unsigned i=0; i < num_muteces; i++)
        pthread_mutex_init (&genome_muteces[i], NULL);
}

void ref_lock_initialize_loaded_genome (void) 
{ 
    ref_lock_initialize_do ((genome_size + GENOME_BASES_PER_MUTEX-1) / GENOME_BASES_PER_MUTEX); // round up
}

void ref_lock_initialize_denovo_genome (void) 
{ 
    ref_lock_initialize_do (REF_NUM_DENOVO_RANGES); 
}

void ref_lock_free (void)
{
    if (genome_muteces) {
        for (unsigned i=0; i < genome_num_muteces; i++)
            pthread_mutex_destroy (&genome_muteces[i]);

        FREE (genome_muteces);
        genome_num_muteces = 0;
    }
}

// lock a region that includes the region given and the flanking regions 
RefLock ref_lock (PosType gpos_start, uint32_t seq_len)
{
    PosType last_pos = gpos_start + seq_len - 1;

    // round to 64 before and after the request region
    RefLock lock = { .first_mutex = gpos_start / GENOME_BASES_PER_MUTEX,
                     .last_mutex  = last_pos   / GENOME_BASES_PER_MUTEX };

    ASSERT (lock.first_mutex >= 0 && lock.first_mutex <= genome_num_muteces, "Error in ref_lock: lock.first_mutex=%u out of range: [0,%u]", lock.first_mutex, genome_num_muteces);
    ASSERT (lock.last_mutex  >= 0 && lock.last_mutex  <= genome_num_muteces, "Error in ref_lock: lock.last_mutex=%u out of range: [0,%u]", lock.last_mutex, genome_num_muteces);
    
    // lock muteces in order
    for (int i=lock.first_mutex; i <= lock.last_mutex; i++)
        mutex_lock (genome_muteces[i]);

    return lock;
}

RefLock ref_unlock (RefLock lock)
{
    if (lock.first_mutex >= 0) {
        // unlock muteces in reverse order
        for (int i=lock.last_mutex; i >= lock.first_mutex; i--)
            mutex_unlock (genome_muteces[i]);
    }

    return REFLOCK_NONE;
}

// used for RT_DENOVO - single mutex er range
RefLock ref_lock_range (int32_t range_id)
{
    RefLock lock = { .first_mutex = range_id, .last_mutex = range_id };
    mutex_lock (genome_muteces[range_id]);

    return lock;
}
