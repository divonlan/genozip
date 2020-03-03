// ------------------------------------------------------------------
//   regions.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "regions.h"
#include "buffer.h"

typedef struct {
    const char *reg_str;
} Region;

static Buffer regions = EMPTY_BUFFER;

void regions_add (const char *reg_str)
{
    while (1) {
        char *one_reg_str = strtok (reg_str, ",");
        if (!one_reg_str) break;

        buf_alloc ()
        char *before_colon = strtok (one_reg_str, ":");


    }
    unsigned len = strlen (reg_str);

    // look for first colon
    for (unsigned i=0; i<len; i++)
}