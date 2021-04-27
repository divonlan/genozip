// ------------------------------------------------------------------
//   iupac.h
//   Copyright (C) 2021-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef IUPAC_INCLUDED
#define IUPAC_INCLUDED

extern void iupac_set (const char *optarg);
extern void iupac_show (void);
extern bool iupac_is_included_ascii (const char *seq, unsigned seq_len);
extern bool iupac_is_included_bam (const char *seq, unsigned seq_len);

#endif
