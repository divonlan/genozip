// ------------------------------------------------------------------
//   stats.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef STATS_INCLUDED
#define STATS_INCLUDED

extern void stats_compress (void); // ZIP
extern void stats_display (void);  // PIZ and ZIP
extern void stats_read_and_display (void);     // PIZ
extern void stats_add_txt_name (const char *fn);
extern void stats_set_consolidation (VBlockP vb, DidIType parent, unsigned num_deps, ...);

#endif
