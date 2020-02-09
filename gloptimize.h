// ------------------------------------------------------------------
//   gloptimize.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GLOPTIMIZE_INCLUDED
#define GLOPTIMIZE_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "move_to_front.h"

extern const char *gl_optimize_dictionary (VariantBlockP vb, Buffer *dict, MtfNode *nodes, unsigned dict_start_char, unsigned num_words);
extern void gl_deoptimize_dictionary (char *data, int len);

#endif