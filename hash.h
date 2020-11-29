// ------------------------------------------------------------------
//   hash.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HASH_INCLUDED
#define HASH_INCLUDED

#include "genozip.h"

extern uint32_t hash_get_estimated_entries (VBlockP merging_vb, ContextP zf_ctx, ConstContextP first_merging_vb_ctx);

extern void hash_alloc_global (ContextP zf_ctx, uint32_t estimated_entries);

typedef enum { HASH_NEW_OK_SINGLETON_IN_VB, HASH_NEW_OK_NOT_SINGLETON, HASH_READ_ONLY } HashGlobalGetEntryMode; 
extern WordIndex hash_global_get_entry (ContextP zf_ctx, const char *snip, unsigned snip_len, HashGlobalGetEntryMode mode,
                                        MtfNodeP *old_node);

extern WordIndex hash_get_entry_for_seg (VBlockP segging_vb, ContextP vb_ctx,
                                         const char *snip, unsigned snip_len, 
                                         WordIndex new_node_i_if_no_old_one,
                                         MtfNodeP *node);

#endif
