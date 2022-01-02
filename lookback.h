// ------------------------------------------------------------------
//   lookback.h
//   Copyright (C) 2021-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "sections.h"

extern void lookback_init (VBlockP vb, ContextP lb_ctx, ContextP ctx, StoreType store_type);

extern void lookback_insert (VBlockP vb, DidIType lb_did_i, DidIType did_i, bool copy_last_value, ValueType value, bool is_word_index);

extern void lookback_insert_txt (VBlockP vb, DidIType lb_did_i, DidIType did_i, STRp(txt));

extern const void *lookback_get_do (VBlockP vb, ContextP lb_ctx, ContextP ctx, uint32_t lookback, bool is_word_index);
#define lookback_get_index(vb, lb_ctx, ctx, lookback) (*(WordIndex *)lookback_get_do ((VBlockP)vb, (lb_ctx), (ctx), lookback, true))
#define lookback_get_value(vb, lb_ctx, ctx, lookback) (*(ValueType *)lookback_get_do ((VBlockP)vb, (lb_ctx), (ctx), lookback, false))

extern bool lookback_is_same_txt (VBlockP vb, DidIType lb_did_i, ContextP ctx, uint32_t lookback, STRp(str));

extern uint32_t lookback_get_next (VBlockP vb, ContextP lb_ctx, ContextP ctx, WordIndex search_for, int64_t *iterator);

extern void lookback_flush (VBlockP vb, ContextP ctx);

extern uint8_t lookback_size_to_local_param (uint32_t size);
