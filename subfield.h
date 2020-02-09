// ------------------------------------------------------------------
//   subfield.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SUBFIELD_INCLUDED
#define SUBFIELD_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#else
#include "compatability/visual_c_stdint.h"
#endif

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

#define SUBFIELD_ID_LEN    ((int)sizeof(uint64_t))    // VCF spec doesn't limit the ID length, we limit it to 8 chars. zero-padded.
typedef union {
    char id[SUBFIELD_ID_LEN];   // \0-padded IDs 
    uint64_t num;               // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
} SubfieldIdType;

#define EMPTY_SUBFIELD_ID { {0,0,0,0,0,0,0,0} }

#pragma pack(pop)

#endif
