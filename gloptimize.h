// ------------------------------------------------------------------
//   gloptimize.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GLOPTIMIZE_INCLUDED
#define GLOPTIMIZE_INCLUDED

#include "genozip.h"

extern const char *gl_optimize_dictionary (VariantBlockP vb, BufferP dict, MtfNodeP nodes, unsigned dict_start_char, unsigned num_words);
extern void gl_deoptimize_dictionary (char *data, int len);

extern bool gl_optimize (const char *snip, unsigned len, char *updated_snip, unsigned *updated_len);

#endif