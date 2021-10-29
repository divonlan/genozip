// ------------------------------------------------------------------
//   optimize.h
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

extern bool optimize_float_2_sig_dig (const char *snip, unsigned len, float cap_value_at /* 0 if no cap */,
                                      char *optimized_snip, unsigned *optimized_snip_len);

extern bool optimize_vector_2_sig_dig (const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

extern bool optimize_vcf_format (DictId dict_id, const char *snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

extern void optimize_phred_quality_string (char *str, unsigned len);
