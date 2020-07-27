// ------------------------------------------------------------------
//   domqual.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DOMQUAL_INCLUDED
#define DOMQUAL_INCLUDED

#include "genozip.h"
#include "data_types.h"

extern bool domqual_convert_qual_to_domqual (VBlockP vb, LocalGetLineCallback get_line, int qual_field);
extern void domqual_reconstruct (VBlockP vb, ContextP qual_ctx);

#endif
