// ------------------------------------------------------------------
//   qname_filter.h
//   Copyright (C) 2024-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

extern void qname_filter_initialize_from_file (rom filename);
extern void qname_filter_initialize_from_opt (rom opt);
extern bool qname_filter_does_line_survive (STRp(qname));
