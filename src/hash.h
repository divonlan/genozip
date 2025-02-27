// ------------------------------------------------------------------
//   hash.h
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern uint32_t hash_get_estimated_entries (VBlockP merging_vb, ContextP zctx, ConstContextP first_merging_vb_ctx);

extern void hash_alloc_global (ContextP zctx, uint32_t estimated_entries);
extern uint32_t hash_next_size_up (uint64_t size, bool allow_huge);

extern WordIndex hash_global_get_entry (ContextP zctx, STRp(snip), bool allow_singleton, bool snip_is_definitely_new, CtxNodeP *please_update_index);

extern WordIndex hash_find_snip_in_ol_nodes (VBlockP vb, ContextP vctx, STRp(snip), rom *snip_in_dict_out);

extern WordIndex hash_get_entry_for_seg (VBlockP segging_vb, ContextP vctx, STRp(snip), rom *snip_in_dict_out);

// tested hash table sizes up to 5M. turns out smaller tables (up to a point) are faster, despite having longer
// average linked lists. probably bc the CPU can store the entire hash and nodes arrays in L1 or L2
// memory cache during segmentation
#define NO_NEXT 0xffffffff
static inline uint32_t hash_do (uint32_t hash_len, STRp(snip))
{
    uint32_t hash;
    if (!hash_len) 
        hash = NO_NEXT; // hash table does not exist

    // spread the snip throughout the 64bit word before taking a mod - to ensure about-even distribution 
    // across the hash table
    else {
        uint64_t result=0;
        for (unsigned i=0; i < snip_len; i++) 
            result = ((result << 23) | (result >> 41)) ^ (uint64_t)((uint8_t)snip[i]);

        hash = result % hash_len;
    }

    return hash;
}
