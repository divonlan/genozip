// ------------------------------------------------------------------
//   bed.h
//   Copyright (C) 2023-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#pragma GENDICT_PREFIX BED

#pragma GENDICT BED_CHROM=DTYPE_FIELD=CHROM             // must be 1st as this is BED's CHROM. 
#pragma GENDICT BED_START=DTYPE_FIELD=START
#pragma GENDICT BED_END=DTYPE_FIELD=END
#pragma GENDICT BED_NAME=DTYPE_FIELD=NAME
#pragma GENDICT BED_SCORE=DTYPE_FIELD=SCORE
#pragma GENDICT BED_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT BED_TSTART=DTYPE_FIELD=TSTART
#pragma GENDICT BED_TEND=DTYPE_FIELD=TEND
#pragma GENDICT BED_RGB=DTYPE_FIELD=RGB
#pragma GENDICT BED_BCOUNT=DTYPE_FIELD=BCOUNT
#pragma GENDICT BED_BSIZES=DTYPE_FIELD=BSIZES
#pragma GENDICT BED_BSTARTS=DTYPE_FIELD=BSTARTS
#pragma GENDICT BED_EOL=DTYPE_FIELD=EOL
#pragma GENDICT BED_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT BED_DEBUG_LINES=DTYPE_FIELD=DBGLINES    // used by --debug-lines

// ZIP Stuff
extern void bed_zip_initialize (void);
extern bool is_bed (STRp(header), bool *need_more);
extern int32_t bed_is_header_done (bool is_eof);
extern rom bed_seg_txt_line (VBlockP vb_, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern void bed_seg_initialize (VBlockP vb_);
extern void bed_seg_finalize (VBlockP vb);
extern bool bed_seg_is_small (ConstVBlockP vb, DictId dict_id);

#define BED_SPECIAL {  bed_piz_special_BCOUNT }
SPECIAL (BED, 0, BCOUNT, bed_piz_special_BCOUNT);
#define NUM_BED_SPECIAL 1
