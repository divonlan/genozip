// ------------------------------------------------------------------
//   random_access.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef RANDOM_ACCESS_INCLUDED
#define RANDOM_ACCESS_INCLUDED

#include <stdint.h>

// we pack this, because it gets written to disk in SEC_RANDOM_ACCESS. On disk, it is stored in Big Endian.

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

typedef struct {
    uint32_t chrom     : 31;              // before merge: node index into chrom context mtf, after merge - word index in CHROM dictionary
    uint32_t is_sorted : 1;               // is this range sorted in a non-decreasing order?
    uint32_t first_pos;                   // POS field value of first position
    int32_t last_pos;                     // in VB, this is in an absolute value. On disk, it is the detla vs first_pos
    uint32_t variant_block_i;             // the vb_i in which this range appears
    uint32_t start_vb_line, num_vb_lines; // corresponds to the line with this VB
} RAEntry; 

#pragma pack(pop)

extern void random_access_merge_in_vb (VariantBlockP vb);

#endif