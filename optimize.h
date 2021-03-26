// ------------------------------------------------------------------
//   optimize.h
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef OPTIMIZE_INCLUDED
#define OPTIMIZE_INCLUDED

#define OPTIMIZE_MAX_SNIP_LEN 300

#include "genozip.h"

extern bool optimize_float_2_sig_dig (const char *snip, unsigned len, double cap_value_at /* 0 if no cap */,
                                      char *optimized_snip, unsigned *optimized_snip_len);

extern bool optimize_vector_2_sig_dig (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

extern bool optimize_vcf_format (DictId dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

extern void optimize_phred_quality_string (char *str, unsigned len);

extern bool optimize_vcf_pl (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

#endif
