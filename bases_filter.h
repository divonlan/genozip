// ------------------------------------------------------------------
//   bases_filter.h
//   Copyright (C) 2021-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once 

extern void iupac_set (const char *optarg);
extern void iupac_show (void);
extern bool iupac_is_included_ascii (const char *seq, unsigned seq_len);
extern bool iupac_is_included_bam (const char *seq, unsigned seq_len);
