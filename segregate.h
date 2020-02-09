// ------------------------------------------------------------------
//   segregate.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SEGREGATE_INCLUDED
#define SEGREGATE_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "subfield.h"

extern void seg_all_data_lines (VariantBlockP vb, Buffer *lines_orig /* for testing */);
extern SubfieldIdType seg_get_subfield (const char **data, unsigned len, unsigned line_i);

#endif