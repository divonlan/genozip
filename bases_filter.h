// ------------------------------------------------------------------
//   bases_filter.h
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once 

extern void iupac_set (rom optarg);
extern void iupac_show (void);
extern bool iupac_is_included_ascii (STRp(seq));
extern bool iupac_is_included_bam (STRp(seq));
