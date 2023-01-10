// ------------------------------------------------------------------
//   hash.h
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern uint32_t hash_get_estimated_entries (VBlockP merging_vb, ContextP zctx, ConstContextP first_merging_vb_ctx);

extern void hash_alloc_global (ContextP zctx, uint32_t estimated_entries);
extern uint32_t hash_next_size_up (uint64_t size, bool allow_huge);

typedef enum { HASH_NEW_OK_SINGLETON_IN_VB, HASH_NEW_OK_NOT_SINGLETON, HASH_READ_ONLY } HashGlobalGetEntryMode; 
extern WordIndex hash_global_get_entry (ContextP zctx, STRp(snip), HashGlobalGetEntryMode mode,
                                        CtxNodeP *old_node);

extern WordIndex hash_get_entry_for_seg (VBlockP segging_vb, ContextP vctx, STRp(snip), 
                                         WordIndex new_node_i_if_no_old_one,
                                         CtxNodeP *node);

// tested hash table sizes up to 5M. turns out smaller tables (up to a point) are faster, despite having longer
// average linked lists. probably bc the CPU can store the entire hash and nodes arrays in L1 or L2
// memory cache during segmentation
#define NO_NEXT 0xffffffff
static inline uint32_t hash_do (uint32_t hash_len, STRp(snip))
{
    if (!hash_len) return NO_NEXT; // hash table does not exist

    // spread the snip throughout the 64bit word before taking a mod - to ensure about-even distribution 
    // across the hash table
    uint64_t result=0;
    for (unsigned i=0; i < snip_len; i++) 
        result = ((result << 23) | (result >> 41)) ^ (uint64_t)((uint8_t)snip[i]);

    return (uint32_t)(result % hash_len);
}
