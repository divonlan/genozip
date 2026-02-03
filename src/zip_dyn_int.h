// ------------------------------------------------------------------
//   dyn_int.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "local_type.h"

extern void dyn_int_init_ctx (VBlockP vb, ContextP ctx, int64_t value);
extern void dyn_int_append (VBlockP vb, ContextP ctx, int64_t value, unsigned add_bytes);
extern void dyn_int_append_nothing_char (VBlockP vb, ContextP ctx, unsigned add_bytes);
extern LocalType dyn_int_get_ltype (ContextP ctx);
extern rom dyn_int_lt_order_name (uint8_t dyn_lt_order);
extern void dyn_int_transpose (VBlockP vb, ContextP ctx);
