// ------------------------------------------------------------------
//   linesorter.h
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// --------
// ZIP side
// --------
extern void linesorter_merge_vb (VBlockP vb);
extern void linesorter_compress_recon_plan (void);
extern void linesorter_show_recon_plan (FileP file, bool is_luft, uint32_t conc_writing_vbs, uint32_t vblock_mb);

typedef struct {
    WordIndex chrom_wi; 
    uint32_t tie_breaker;
    PosType start_pos, end_pos;
} LineCmpInfo; 
extern int linesorter_cmp (LineCmpInfo a, LineCmpInfo b);

// --------
// PIZ side
// --------
extern uint32_t linesorter_get_max_conc_writing_vbs (void);
