// ------------------------------------------------------------------
//   stats.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

extern void stats_generate (void); // ZIP
extern void stats_display (void);  // PIZ and ZIP
extern void stats_read_and_display (void);     // PIZ
extern void stats_add_txt_name (rom fn);
