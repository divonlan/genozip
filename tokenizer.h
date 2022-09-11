// ------------------------------------------------------------------
//   tokenizer.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

extern const char sep_with_space[], sep_without_space[];

extern void tokenizer_seg (VBlockP vb, ContextP field_ctx, STRp(field), rom is_sep, unsigned add_additional_bytes);
