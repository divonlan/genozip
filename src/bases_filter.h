// ------------------------------------------------------------------
//   bases_filter.h
//   Copyright (C) 2021-2026 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#pragma once 

extern void iupac_set (rom optarg);
extern void iupac_show (void);
extern bool iupac_is_included_ascii (STRp(seq));
extern bool iupac_is_included_bam (STRp(seq));
