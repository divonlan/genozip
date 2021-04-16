// ------------------------------------------------------------------
//   linesorter.h
//   Copyright (C) 2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef LINESORTER_INCLUDED
#define LINESORTER_INCLUDED

#include "genozip.h"

// a recon plan entry - one per a group of consecutive lines of Primary or Luft
typedef struct {
    uint32_t vblock_i;
    uint32_t start_line, num_lines;
    WordIndex chrom_wi;
    PosType start_pos, end_pos;
} LineInfo;

// --------
// ZIP side
// --------
extern void linesorter_merge_vb (VBlockP vb);
extern void linesorter_compress_recon_plan (void);
extern void linesorter_show_recon_plan (FileP file, bool is_luft, uint32_t conc_writing_vbs, uint32_t vblock_mb);

// --------
// PIZ side
// --------
extern uint32_t linesorter_get_max_conc_writing_vbs (void);

#endif

