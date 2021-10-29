// ------------------------------------------------------------------
//   lookback.h
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "sections.h"

extern void lookback_init (VBlockP vb, ContextP ctx, StoreType store_type);

#define TAKE_LAST_VALUE -10
extern void lookback_insert (VBlockP vb, DidIType did_i, int64_t value, bool is_word_index);

extern const void *lookback_get_do (VBlockP vb, ContextP ctx, uint32_t lookback, bool is_word_index);
#define lookback_get_index(vb, ctx, lookback) (*(WordIndex *)lookback_get_do ((VBlockP)vb, (ctx), lookback, true))
#define lookback_get_value(vb, ctx, lookback) (*(int64_t   *)lookback_get_do ((VBlockP)vb, (ctx), lookback, false))

extern uint32_t lookback_get_next (VBlockP vb, ContextP ctx, WordIndex search_for, int64_t *iterator);

extern void lookback_flush (VBlockP vb, ContextP ctx);
