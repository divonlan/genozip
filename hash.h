// ------------------------------------------------------------------
//   hash.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HASH_INCLUDED
#define HASH_INCLUDED

#include "genozip.h"

extern void hash_alloc_global (VBlock *merging_vb, Context *zf_ctx, const Context *first_merging_vb_ctx);

typedef enum { HASH_NEW_OK_SINGLETON_IN_VB, HASH_NEW_OK_NOT_SINGLETON, HASH_MUST_EXIST } HashGlobalGetEntryMode; 
extern int32_t hash_global_get_entry (Context *zf_ctx, const char *snip, unsigned snip_len, HashGlobalGetEntryMode mode,
                                      MtfNode **old_node);

extern int32_t hash_get_entry_for_seg (VBlock *segging_vb, Context *vb_ctx,
                                       const char *snip, unsigned snip_len, 
                                       int32_t new_mtf_i_if_no_old_one,
                                       MtfNode **node);

#endif
