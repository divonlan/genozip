// ------------------------------------------------------------------
//   optimize.h
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern bool optimize_float_2_sig_dig (rom snip, unsigned len, float cap_value_at /* 0 if no cap */,
                                      char *optimized_snip, unsigned *optimized_snip_len);

extern bool optimize_vector_2_sig_dig (rom snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

extern bool optimize_vcf_format (DictId dict_id, rom snip, unsigned len, char *optimized_snip, unsigned *optimized_snip_len);

extern void optimize_phred_quality_string (char *str, unsigned len);
