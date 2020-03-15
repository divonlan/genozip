// ------------------------------------------------------------------
//   hash.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef HASH_INCLUDED
#define HASH_INCLUDED

#include "genozip.h"

extern void hash_alloc_global (VariantBlock *merging_vb, MtfContext *zf_ctx, const MtfContext *first_merging_vb_ctx);

extern int32_t hash_get_entry_for_merge (MtfContext *zf_ctx, const char *snip, unsigned snip_len, 
                                         int32_t new_mtf_i_if_no_old_one,
                                         MtfNode **node);

extern int32_t hash_get_entry_for_seg (VariantBlock *segging_vb, MtfContext *vb_ctx,
                                       const char *snip, unsigned snip_len, 
                                       int32_t new_mtf_i_if_no_old_one,
                                       MtfNode **node);

#endif
