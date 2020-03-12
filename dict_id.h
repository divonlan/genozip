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

#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF spec doesn't limit the ID length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union {
    uint64_t num;            // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN]; // \0-padded IDs 
} DictIdType;

#pragma pack(pop)

static inline DictIdType dict_id_make(const char *str, unsigned str_len) { DictIdType dict_id = {0}; memcpy (dict_id.id, str, MIN (str_len, DICT_ID_LEN)); return dict_id;}

#define dict_id_is(dict_id, str) (dict_id_make(str, strlen(str)).num == dict_id_printable (dict_id).num)

// 2 MSb of first byte determine dictionary type

#define dict_id_is_format_subfield(dict_id) ((dict_id.id[0] >> 6) == 1)
#define dict_id_is_info_subfield(dict_id)   ((dict_id.id[0] >> 6) == 3)
#define dict_id_is_vardata_field(dict_id)   ((dict_id.id[0] >> 6) == 0)

static inline DictIdType dict_id_format_subfield(DictIdType dict_id) {                                       return dict_id; } // no change - keep Msb 01
static inline DictIdType dict_id_info_subfield(DictIdType dict_id)   { dict_id.id[0] = dict_id.id[0] | 0xc0; return dict_id; } // set 2 Msb to 11
static inline DictIdType dict_id_vardata_field(DictIdType dict_id)   { dict_id.id[0] = dict_id.id[0] & 0x3f; return dict_id; } // set 2 Msb to 00

static inline DictIdType dict_id_printable(DictIdType dict_id) { dict_id.id[0] = (dict_id.id[0] & 0x7f) | 0x40; return dict_id; } // set 2 Msb to 01

extern DictIdType dict_id_show_one_b250, dict_id_show_one_dict; // arguments of --show-b250-one and --show-dict-one (defined in genozip.c)
extern DictIdType dict_id_dump_one_b250;                        // arguments of --dump-b250-one (defined in genozip.c)

extern uint64_t dict_id_PL, dict_id_GL; // declared in segregate.c

static inline void dict_id_initialize(void) { dict_id_PL = dict_id_make ("PL", 2).num;\
                                              dict_id_GL = dict_id_make ("GL", 2).num;} // note: these values will be different in big and little endian machines

#endif
