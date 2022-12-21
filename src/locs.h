// ------------------------------------------------------------------
//   locs.h
//   Copyright (C) 2021-2022 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once

#include "genozip.h"

#pragma GENDICT_PREFIX LOCS
#pragma GENDICT LOCS_X=DTYPE_FIELD=X
#pragma GENDICT LOCS_Y=DTYPE_FIELD=Y
#pragma GENDICT LOCS_TOPLEVEL=DTYPE_FIELD=TOPLEVEL
#pragma GENDICT LOCS_DEBUG_LINES=DTYPE_FIELD=DBGLINES      // used by --debug-lines

extern int32_t locs_is_header_done (bool is_eof);
extern int32_t locs_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern bool locs_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern void locs_seg_initialize (VBlockP vb);
extern void locs_seg_finalize (VBlockP vb);
extern rom locs_seg_txt_line (VBlockP vb, rom field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);

#define LOCS_SPECIAL { locs_piz_special_DELTA_FLOAT }
SPECIAL (LOCS, 0, DELTA_FLOAT, locs_piz_special_DELTA_FLOAT);
#define NUM_LOCS_SPECIAL 1

// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (they are included in the TOPLEVEL container)
// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (LOCS, LOCS, 1, LTEN_F32, container_translate_LTEN_F32) 
#define NUM_LOCS_TRANS   2 // including "none"
#define LOCS_TRANSLATORS { NULL /* none */, container_translate_LTEN_F32 }
