// ------------------------------------------------------------------
//   tokenizer.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef TOKENIZER_INCLUDED
#define TOKENIZER_INCLUDED

#include "genozip.h"

extern const char sep_with_space[], sep_without_space[];

extern void tokenizer_seg (VBlockP vb, ContextP field_ctx, STRp(field), const char *is_sep, unsigned add_additional_bytes);

#endif
