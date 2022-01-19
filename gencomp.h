// ------------------------------------------------------------------
//   gencomp.h - "generated component"
//   Copyright (C) 2022-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// same values to Coords and SamComponentType
#define MAX_GENCOMP_NUM 2
typedef int GenCompNum; // 0=NONE, 1,2=Generated Components

extern void gencomp_initialize_file (GenCompNum file_gc, const char *name);
extern void gencomp_append_file (VBlockP vb, GenCompNum file_gc);
