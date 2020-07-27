// ------------------------------------------------------------------
//   me23.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt


#ifndef ME23_INCLUDED
#define ME23_INCLUDED

#include "genozip.h"
#include "sections.h"

extern const char *me23_seg_txt_line (VBlockP vb_, const char *field_start_line, bool *has_special_eol);
extern void me23_piz_reconstruct_vb (VBlockP vb);
extern void me23_seg_initialize (VBlockP vb);
extern bool me23_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);

#define ME23_DICT_ID_ALIASES

#endif
