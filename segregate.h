// ------------------------------------------------------------------
//   segregate.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SEGREGATE_INCLUDED
#define SEGREGATE_INCLUDED

#include "genozip.h"
#include "dict_id.h"

extern void seg_all_data_lines (VariantBlockP vb);

extern DictIdType seg_get_format_subfield (const char **data, uint32_t *len, unsigned line_i);

#endif