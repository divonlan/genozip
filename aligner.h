// ------------------------------------------------------------------
//   aligner.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ALIGNER_INCLUDED
#define ALIGNER_INCLUDED

#include "genozip.h"

extern void aligner_seg_seq (VBlockP vb, ContextP bitmap_ctx, const char *seq, uint32_t seq_len);
extern void aligner_reconstruct_seq (VBlockP vb, ContextP bitmap_ctx, uint32_t seq_len, bool is_pair_2);

extern BitArray aligner_seq_to_bitmap (const char *seq, word_t seq_len, word_t *bitmap_words, bool *seq_is_all_actg);

#endif
