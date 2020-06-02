
// ------------------------------------------------------------------
//   sam_ref.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <pthread.h>
#include "sam_private.h"
#include "buffer.h"

#define NUM_SITES_PER_RANGE (1 << 10) // 1 MB

typedef struct {
    const char *chrom_name;
    unsigned chrom_name_len;
    uint64_t pos;
    Buffer *ref; // note that if we realloc ranges, the ref buffer itself isn't realloced (since this is a pointer)
    pthread_mutex_t mutex;
} Range;

// reference sequences - one per chrom
static Buffer ranges = EMPTY_BUFFER; // a buffer containing a array of Range
static pthread_mutex_t ranges_mutex;
static ranges_initialized = false;

// Allocated and initializes the ref and mutex buffers for the given chrom/range
// Note: this is the only function that accesses the ranges buffer directly during the compute threads stage
static void sam_ref_get_range (const char *chrom_name, unsigned chrom_name_len, uint64_t pos, 
                               Buffer **ref_buf, pthread_mutex_t *range_mutex) // out
{
    // one buffer per chrom - indexed by the 5 LSb of the last 3 characters of the chrom name.
    // each buffer contains an array of uint32 indeces into ranges (index for each range)
    // note that in very rare case that two chroms map to the same chrom_ranges element - then both will share the same ref data
    // for this range - causing ineffecient compression for this range, but no incorrectnesss
    static Buffer chrom_ranges[32][32][32]; 
    
    if (!ranges_initialized) { // this happens during vb_i=1 seg, which is run under a vb=1 lock in zip_compress_one_vb()
        pthread_mutex_init (&ranges_mutex, NULL);
        memset (chrom_ranges, 0, sizeof (chrom_ranges));
        ranges_initialized = true;
    }

    pthread_mutex_lock (&range_mutex);

    Range *range = get_range();

    if (!range->ref) {
        buf_alloc (evb, range->ref, NUM_SITES_PER_RANGE, 1, "range->ref", 0);
        pthread_mutex_init (&range->mutex, NULL);
    }

    *ref_buf     = range->ref;
    *range_mutex = range->mutex;

    pthread_mutex_unlock (&range_mutex);
}

void sam_ref_refy_one_seq (char *seq, uint32_t seq_len, uint64_t pos,
                           const char *chrom_name, unsigned chrom_name_len)
{
    // we split the seq into up to two ranges for simpliticy - if there is ever a need for more (super long seqs) we can fix this
    ASSERT (seq_len < NUM_SITES_PER_RANGE, "Error in sam_ref_refy_one_seq: seq_len=%u exceeds maximum supported of %u",
            seq_len, NUM_SITES_PER_RANGE);

    // if seq spans two ranges, call recursively for each range
    if (pos / NUM_SITES_PER_RANGE != (pos + (uint64_t)seq_len -1) / NUM_SITES_PER_RANGE) { 
        uint32_t sub1_len = NUM_SITES_PER_RANGE - (seq_len % NUM_SITES_PER_RANGE);
        sam_ref_refy_one_seq (seq, sub1_len, pos, chrom_name, chrom_name_len);
        sam_ref_refy_one_seq (seq + sub1_len, seq_len - sub1_len, pos + sub1_len, chrom_name, chrom_name_len);
        return;
    }        

    Buffer *ref_buf=NULL;
    pthread_mutex_t range_mutex;
    bool range_mutex_is_locked = false;

    sam_ref_get_range (chrom_name, chrom_name_len, pos, &ref_buf, &range_mutex);
    
    ARRAY (char, ref, *ref_buf);    

    uint64_t final_pos   = (pos + (uint64_t)seq_len - 1);
    uint32_t first_range = (uint32_t)(pos / NUM_SITES_PER_RANGE);
    uint32_t final_range = (uint32_t)(final_pos / NUM_SITES_PER_RANGE);

    for (uint32_t i=0; i < seq_len; i++) {

        if (ref[pos+i] == seq[i])
            seq[i] = '-'; // ref
        
        else if (!ref[pos+i]) { // no ref yet for this site - we will be the ref
            if (range_mutex_is_locked) {
                pthread_mutex_lock (&range_mutex); // lock on our first change to this ref
                range_mutex_is_locked = true;
            }

            ref[pos+1] = seq[i];
        }

        else {} // if ref exists but is different than SEQ - do nothing. SEQ remains
    }

    if (range_mutex_is_locked)
        pthread_mutex_unlock (&range_mutex);
}
