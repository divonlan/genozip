// ------------------------------------------------------------------
//   dict_id.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DICT_ID_INCLUDED
#define DICT_ID_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#else
#include "compatability/visual_c_stdint.h"
#endif

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF spec doesn't limit the ID length, we limit it to 8 chars. zero-padded.
typedef union {
    char id[DICT_ID_LEN];   // \0-padded IDs 
    uint64_t num;           // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
} DictIdType;

#define EMPTY_DICT_ID { {0,0,0,0,0,0,0,0} }

#pragma pack(pop)

// 2 MSb of first byte determine dictionary type
#define dict_id_is_gt_subfield(dict_id)   ((dict_id.id[0] >> 6) == 1)
#define dict_id_is_vardata_field(dict_id) ((dict_id.id[0] >> 6) == 0)
#define dict_id_is_info_subfield(dict_id) ((dict_id.id[0] >> 6) == 3)

#endif
