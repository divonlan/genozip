// ------------------------------------------------------------------
//   optimize.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern bool optimize_float_2_sig_dig (STRp(snip), float cap_value_at, STRe(optimized_snip));

extern bool optimize_vector_2_sig_dig (STRp(snip), STRe(optimized_snip));

extern void optimize_phred_quality_string (STRc(qual));
