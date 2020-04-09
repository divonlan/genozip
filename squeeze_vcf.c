// ------------------------------------------------------------------
//   squeeze.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// A module for squeezing an index array. An index array of length len, in an array in which each entry containing a unique number 0->(len-1). 
// Squeezing the array means using the minimal number of bits for each entry.

#include <math.h>
#include "genozip.h"
#include "profiler.h"
#include "squeeze_vcf.h"
#include "vblock.h"

unsigned squeeze_bits_per_entry (unsigned len)
{
    return ceil (log2(len));
}

// return the number of bytes needed for a squeezed array
unsigned squeeze_len(unsigned int len)
{
    return ceil(((float)(squeeze_bits_per_entry(len) * len)) / 8.0); 
}

void squeeze (VBlockVCF *vb,
              uint8_t *squeezed, // memory should be pre-allocated by caller
              uint16_t *squeezed_checksum,
              const unsigned *normal, 
              unsigned normal_len)
{
    START_TIMER;

    unsigned squeezed_len = squeeze_len (normal_len);
    memset (squeezed, 0, squeezed_len);

    unsigned num_bits_per_entry = squeeze_bits_per_entry(normal_len);

    unsigned next_squeezed_byte_index=0;
    unsigned next_squeezed_bit_index=0;

    for (unsigned i=0; i < normal_len; i++) {

        unsigned val = normal[i];
        unsigned num_remaining_bits_this_entry = num_bits_per_entry;

        ASSERT (val < normal_len, "Error: invalid value in permutation index. val=%u normal_len=%u", val, normal_len);

        while (num_remaining_bits_this_entry) {
            unsigned num_bits_this_iteration = MIN (num_remaining_bits_this_entry, 8 - next_squeezed_bit_index);
            unsigned mask = ((1<<num_bits_this_iteration)-1);
            uint8_t bits_to_add = (uint8_t)(val & mask);

            squeezed[next_squeezed_byte_index] |= bits_to_add << next_squeezed_bit_index;

            next_squeezed_bit_index += num_bits_this_iteration;

            if (next_squeezed_bit_index == 8) {
                next_squeezed_bit_index = 0;
                next_squeezed_byte_index++;
            }

            val >>= num_bits_this_iteration;
            num_remaining_bits_this_entry -= num_bits_this_iteration;
        }
    }

    *squeezed_checksum = 0;
    for (unsigned i=0; i < squeezed_len; i++)
        (*squeezed_checksum) += (uint16_t)squeezed[i];

    COPY_TIMER (vb->profile.squeeze)
}

static unsigned unsqueeze_get_number(const uint8_t *squeezed, unsigned start_bit, unsigned end_bit)
{
    unsigned number = 0;

    unsigned start_byte = start_bit / 8;
    unsigned end_byte   = end_bit   / 8;

    // if number is contained entirely within a single byte
    if (start_byte == end_byte) {
        unsigned num_high_bits_to_mask_out = ((end_bit+1) % 8);
        uint8_t mask = num_high_bits_to_mask_out ? (1 << num_high_bits_to_mask_out) - 1 : 0xff;
        return (mask & squeezed[start_byte]) >> (start_bit % 8);
    }

    // number spans multiple bytes
    unsigned bits_so_far=0;
    for (unsigned byte_i=start_byte; byte_i <= end_byte ; byte_i++) {

        // in start_byte, we might only use the MSb's
        if (byte_i == start_byte) {
            number = squeezed[byte_i] >> (start_bit % 8);
            bits_so_far = 8 - (start_bit % 8);
        }

        // in end_byte, we might only use the LSb's
        else if (byte_i == end_byte) {
            unsigned num_bits_this_iteration = (end_bit+1) % 8 ? (end_bit+1) % 8 : 8;
            unsigned mask = ((1<<num_bits_this_iteration)-1);
            number += (mask & squeezed[byte_i]) << bits_so_far;
            //number += (((1 << ((end_bit+1) % 8)) -1) & squeezed[byte_i]) << bits_so_far;
        }
        // middle byte - use whole byte
        else {
            number += squeezed[byte_i]  << bits_so_far;
            bits_so_far += 8;
        }
    }

    return number;
}

void unsqueeze (VBlockVCF *vb,
                unsigned *normal, // memory should be pre-allocated by caller
                const uint8_t *squeezed, 
                uint16_t squeezed_checksum,
                unsigned normal_len)
{
    START_TIMER;
    
    // handle the trivial case where there is only one haplotype per line, and therefore no need for an index.
    // we construct an index for consistency anyway
    if (vb->num_haplotypes_per_line == 1) {
        normal[0] = 0;
        return;
    }    
    
    ASSERT (squeezed, "Error: squeezed is NULL in vblock_i=%u", vb->vblock_i);

    unsigned squeezed_len = squeeze_len (normal_len);

    uint16_t my_squeezed_checksum = 0;
    for (unsigned i=0; i < squeezed_len; i++)
        my_squeezed_checksum += (uint16_t)squeezed[i];

    ASSERT (squeezed_checksum == my_squeezed_checksum, "Error: failed squeezed checksum. Expected: %u seeing: %u", squeezed_checksum, my_squeezed_checksum);

    unsigned num_bits_per_entry = squeeze_bits_per_entry(normal_len);

    for (unsigned i=0; i < normal_len; i++) {

        unsigned start_bit = i * num_bits_per_entry;
        unsigned end_bit   = (i+1) * num_bits_per_entry - 1;

        normal[i] = unsqueeze_get_number (squeezed, start_bit, end_bit);
        ASSERT (normal[i] < normal_len, "Error: corruption while unsqueezing i=%u normal[i]=%u normal_len=%u", i, normal[i], normal_len);
    }

    // check that the index contains all numbers (if adding and multiplying all is the expected value)
    unsigned expected_sum  = 0, normal_sum  = 0;
    unsigned expected_prod = 1, normal_prod = 1; // product may overflow - that's ok for our purposes
    for (unsigned i=0; i < normal_len; i++) {
        expected_sum  += i       ; normal_sum  += normal[i];
        expected_prod *= (i | 1) ; normal_prod *= (normal[i] | 1);
    }
    
    ASSERT0 (expected_sum == normal_sum && expected_prod == normal_prod,
             "Error: unsequeeze failed - index is not a permutation");

   COPY_TIMER (vb->profile.squeeze)
}

void squeeze_unit_test()
{
    unsigned normal[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};
    uint8_t squeezed[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned unsqueezed[32];

    uint16_t checksum;

    VBlockVCF vb;
    vb.num_haplotypes_per_line = sizeof (unsqueezed) / sizeof (unsqueezed[0]);
    squeeze(&vb, squeezed, &checksum, normal, 32);


    fprintf (stderr, "original: %d %d %d %d %d\n", normal[0], normal[1], normal[2], normal[3], normal[4]);

    fprintf (stderr, "squeezed: %x %x %x %x %x\n", squeezed[0], squeezed[1], squeezed[2], squeezed[3], squeezed[4]);

    unsqueeze (&vb, unsqueezed, squeezed, checksum, 32);

    fprintf (stderr, "unsqueezed: ");
    for (unsigned i=0; i<32; i++) fprintf (stderr, "%d ", unsqueezed[i]);
    fprintf (stderr, "\n");
}