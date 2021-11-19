// ------------------------------------------------------------------
//   dict_id.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include <inttypes.h>
#include "genozip.h"
#include "endianness.h"

#define DTYPE_PLAIN DTYPE_2
extern DictId dict_id_make (STRp(str), DictIdType dict_id_type);

#define dict_id_is(dict_id, str) (dict_id_make (str, strlen(str)).num == dict_id_typeless (dict_id).num)
static inline DictIdType dict_id_type (DictId dict_id) { return ((dict_id.id[0] >> 6) == 0) ? DTYPE_FIELD   
                                                              : ((dict_id.id[0] >> 6) == 1) ? DTYPE_2 
                                                              :                               DTYPE_1; } 
static inline bool dict_id_is_field (DictId dict_id) { return dict_id_type(dict_id) == DTYPE_FIELD; } 
static inline bool dict_id_is_type_1(DictId dict_id) { return dict_id_type(dict_id) == DTYPE_1;     }
static inline bool dict_id_is_type_2(DictId dict_id) { return dict_id_type(dict_id) == DTYPE_2;     }

static inline DictId dict_id_typeless(DictId dict_id) { dict_id.id[0] = (dict_id.id[0] & 0x7f) | 0x40; return dict_id; } // set 2 Msb to 01

#define DICT_ID_NONE ((DictId)(uint64_t)0)

typedef struct { DictId alias, dst; } DictIdAlias;
extern const DictIdAlias *dict_id_aliases;
extern uint32_t dict_id_num_aliases;

extern BufferP dict_id_create_aliases_buf (void);
extern void dict_id_read_aliases (void) ;

// template can be 0 - anything OR a type - must 2 MSb of id[0] are used OR a specific dict_id
// candidate is a specific dict_id that we test for its matching of the template
extern bool dict_id_is_match (DictId template, DictId candidate);

extern const char *dict_id_display_type (DataType dt, DictId dict_id);

typedef struct { char s[20]; } DisplayPrintId;
extern DisplayPrintId dis_dict_id (DictId dict_id);
#define dtype_name_z(dict_id) DTPZ(dtype_names)[dict_id_type(dict_id)]

// constant expressions (if s is a string literal) - these generate the same dict_id.num as dict_id_make()

// Type 1 - first character is 192->255
#define DICT_ID_MAKE1_1(s) ((((uint64_t)(s[0] | 0xc0))))
#define DICT_ID_MAKE1_2(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8)))
#define DICT_ID_MAKE1_3(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16)))
#define DICT_ID_MAKE1_4(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24)))
#define DICT_ID_MAKE1_5(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32)))
#define DICT_ID_MAKE1_6(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32) | ((uint64_t)s[5] << 40)))
#define DICT_ID_MAKE1_7(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32) | ((uint64_t)s[5] << 40) | ((uint64_t)s[6] << 48)))
#define DICT_ID_MAKE1_L(s) ((((uint64_t)(s[0] | 0xc0)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[sizeof s-5] << 32) | ((uint64_t)s[sizeof s-4] << 40) | ((uint64_t)s[sizeof s-3] << 48) | ((uint64_t)s[sizeof s-2] << 56)))

// Type 2 - first character is unchanged - 64-127 (lower or uppercase, but not numbers)
#define DICT_ID_MAKE2_1(s) (((uint64_t)s[0]))
#define DICT_ID_MAKE2_2(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8)))
#define DICT_ID_MAKE2_3(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16)))
#define DICT_ID_MAKE2_4(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24)))
#define DICT_ID_MAKE2_5(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32)))
#define DICT_ID_MAKE2_6(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32) | ((uint64_t)s[5] << 40)))
#define DICT_ID_MAKE2_7(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32) | ((uint64_t)s[5] << 40) | ((uint64_t)s[6] << 48)))
#define DICT_ID_MAKE2_L(s) (((uint64_t)s[0] | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[sizeof s-5] << 32) | ((uint64_t)s[sizeof s-4] << 40) | ((uint64_t)s[sizeof s-3] << 48) | ((uint64_t)s[sizeof s-2] << 56)))

// Field - first character is 0->63
#define DICT_ID_MAKEF_1(s) ((((uint64_t)(s[0] & 0x3f))))
#define DICT_ID_MAKEF_2(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8)))
#define DICT_ID_MAKEF_3(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16)))
#define DICT_ID_MAKEF_4(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24)))
#define DICT_ID_MAKEF_5(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32)))
#define DICT_ID_MAKEF_6(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32) | ((uint64_t)s[5] << 40)))
#define DICT_ID_MAKEF_7(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[4] << 32) | ((uint64_t)s[5] << 40) | ((uint64_t)s[6] << 48)))
#define DICT_ID_MAKEF_L(s) ((((uint64_t)(s[0] & 0x3f)) | ((uint64_t)s[1] << 8) | ((uint64_t)s[2] << 16) | ((uint64_t)s[3] << 24) | ((uint64_t)s[sizeof s-5] << 32) | ((uint64_t)s[sizeof s-4] << 40) | ((uint64_t)s[sizeof s-3] << 48) | ((uint64_t)s[sizeof s-2] << 56)))
