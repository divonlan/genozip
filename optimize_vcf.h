// ------------------------------------------------------------------
//   optimize_vcf.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef OPTIMIZE_INCLUDED
#define OPTIMIZE_INCLUDED

#define OPTIMIZE_MAX_SNIP_LEN 300

#include "dict_id.h"
#include <stdbool.h>

extern bool optimize_format (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);
extern bool optimize_info (DictIdType dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

#endif
