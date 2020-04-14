// ------------------------------------------------------------------
//   piz.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef PIZ_INCLUDED
#define PIZ_INCLUDED

#include "genozip.h"

extern bool piz_dispatcher (const char *z_basename, unsigned max_threads, bool is_first_vcf_component, bool is_last_file);

extern int32_t piz_decode_pos (int32_t last_pos, const char *delta_snip, unsigned delta_snip_len, char *pos_str, unsigned *pos_len);

extern void piz_vcf_uncompress_one_vb (VBlockP vb);
extern void piz_sam_uncompress_one_vb (VBlockP vb);

#endif
