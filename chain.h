// ------------------------------------------------------------------
//   chain.h
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CHAIN_INCLUDED
#define CHAIN_INCLUDED

#include "genozip.h"
#include "sections.h"

// txtfile
extern int32_t chain_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */);

// seg
extern void chain_seg_initialize (VBlockP vb);
extern void chain_seg_finalize (VBlockP vb);
extern const char *chain_seg_txt_line (VBlockP vb, const char *line, uint32_t remaining_txt_len, bool *has_13);
extern bool chain_seg_is_small (ConstVBlockP vb, DictId dict_id);

// importing and applying chain 
extern void chain_load_chain (void);

// piz
extern CONTAINER_FILTER_FUNC (chain_piz_filter);

#define CHAIN_SPECIAL { chain_piz_special_BACKSPACE, chain_piz_special_QEND, chain_piz_special_SIZE }
SPECIAL (CHAIN, 0, BACKSPACE, chain_piz_special_BACKSPACE);
SPECIAL (CHAIN, 1, QEND,      chain_piz_special_QEND);
SPECIAL (CHAIN, 2, SIZE,      chain_piz_special_SIZE);
#define NUM_CHAIN_SPECIAL 3

#endif
