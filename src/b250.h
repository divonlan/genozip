// ------------------------------------------------------------------
//   b250.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "sections.h"

// ZIP
extern WordIndex b250_seg_get_wi (const uint8_t *msb_p);
#define b250_seg_get_last(ctx) b250_seg_get_wi (BLST8 ((ctx)->b250))

extern void b250_seg_append (VBlockP vb, ContextP vctx, WordIndex node_index);
extern void b250_seg_remove_last (VBlockP vb, ContextP ctx, WordIndex node_index);
extern bool b250_zip_generate (VBlockP vb, ContextP ctx);


// PIZ
extern WordIndex b250_piz_decode (bytes *b, bool advance, B250Size b250_size, rom ctx_name);
