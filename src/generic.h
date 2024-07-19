// ------------------------------------------------------------------
//   generic.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

#pragma GENDICT_PREFIX GNRIC
#pragma GENDICT GNRIC_DATA=DTYPE_FIELD=DATA
#pragma GENDICT GNRIC_TOPLEVEL=DTYPE_FIELD=TOPLEVEL

extern int32_t generic_unconsumed (VBlockP vb, uint32_t first_i);
extern int32_t generic_is_header_done (bool is_eof);
extern void generic_seg_initialize (VBlockP vb);
extern rom generic_seg_txt_line (VBlockP vb, rom next_line, uint32_t remaining_txt_len, bool *has_13);
extern rom generic_assseg_line (VBlockP vb);
extern void generic_seg_finalize (VBlockP vb);
extern bool generic_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern StrTextLong generic_get_magic (void);
extern rom generic_get_ext (void);
extern rom fallback_to_generic (VBlockP vb);

// SPECIALS
SPECIAL (GNRIC, 0, TOPLEVEL, generic_piz_TOPLEVEL);
