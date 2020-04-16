// ------------------------------------------------------------------
//   dict_id.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DICT_ID_INCLUDED
#define DICT_ID_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#else
#include "compatibility/visual_c_stdint.h"
#endif
#include "genozip.h"
#include "header.h"

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF spec doesn't limit the ID length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union {
    uint64_t num;            // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN]; // \0-padded IDs 
} DictIdType;

#pragma pack(pop)

static inline DictIdType dict_id_make(const char *str, unsigned str_len) { DictIdType dict_id = {0}; memcpy (dict_id.id, str, MIN (str_len, DICT_ID_LEN)); return dict_id;}
static inline DictIdType dict_id_field(    DictIdType dict_id) { dict_id.id[0] = dict_id.id[0] & 0x3f; return dict_id; } // set 2 Msb to 00

#define dict_id_is(dict_id, str) (dict_id_make(str, strlen(str)).num == dict_id_printable (dict_id).num)
#define dict_id_is_field(dict_id)         ((dict_id.id[0] >> 6) == 0) // 2 MSb of first byte determine dictionary type

// VCF field types
#define dict_id_is_vcf_info_sf(dict_id)   ((dict_id.id[0] >> 6) == 3)
#define dict_id_is_vcf_format_sf(dict_id) ((dict_id.id[0] >> 6) == 1)

static inline DictIdType dict_id_vcf_info_sf(  DictIdType dict_id) { dict_id.id[0] = dict_id.id[0] | 0xc0; return dict_id; } // set 2 Msb to 11
static inline DictIdType dict_id_vcf_format_sf(DictIdType dict_id) {                                       return dict_id; } // no change - keep Msb 01

// SAM field types - overload the VCF dict id types
#define dict_id_is_sam_qname_sf dict_id_is_vcf_info_sf
#define dict_id_is_sam_optnl_sf dict_id_is_vcf_format_sf

#define dict_id_sam_qname_sf dict_id_vcf_info_sf
#define dict_id_sam_optnl_sf dict_id_vcf_format_sf

static inline DictIdType dict_id_printable(DictIdType dict_id) { dict_id.id[0] = (dict_id.id[0] & 0x7f) | 0x40; return dict_id; } // set 2 Msb to 01

extern int dict_id_get_field (DictIdType dict_id);

extern DictIdType dict_id_show_one_b250, dict_id_show_one_dict; // arguments of --show-b250-one and --show-dict-one (defined in genozip.c)
extern DictIdType dict_id_dump_one_b250;                        // arguments of --dump-b250-one (defined in genozip.c)

extern uint64_t dict_id_vcf_fields[], dict_id_sam_fields[], 
                dict_id_FORMAT_PL, dict_id_FORMAT_GL, dict_id_FORMAT_GP, // some VCF FORMAT subfields
                dict_id_INFO_AC, dict_id_INFO_AF, dict_id_INFO_AN, dict_id_INFO_DP, dict_id_INFO_VQSLOD, // some VCF INFO subfields
                dict_id_INFO_13,

                // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
                dict_id_OPTION_AM, dict_id_OPTION_AS, dict_id_OPTION_CM, dict_id_OPTION_E2, dict_id_OPTION_LB, dict_id_OPTION_FI, 
                dict_id_OPTION_H0, dict_id_OPTION_H1, dict_id_OPTION_H2, dict_id_OPTION_MQ, dict_id_OPTION_NH, dict_id_OPTION_NM, 
                dict_id_OPTION_OA, dict_id_OPTION_OC, dict_id_OPTION_PG, dict_id_OPTION_PQ, dict_id_OPTION_PU, dict_id_OPTION_RG, 
                dict_id_OPTION_SA, dict_id_OPTION_SM, dict_id_OPTION_TC, dict_id_OPTION_U2, dict_id_OPTION_UQ, dict_id_OPTION_CC, 
                dict_id_OPTION_CG, dict_id_OPTION_MC,
    
                // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
                dict_id_OPTION_X0, dict_id_OPTION_X1, dict_id_OPTION_XA, dict_id_OPTION_XN, dict_id_OPTION_XM, dict_id_OPTION_XO,
                dict_id_OPTION_XG, dict_id_OPTION_XS, dict_id_OPTION_XE,

                dict_id_OPTION_STRAND;

extern void dict_id_initialize (DataType data_type);

extern const char *dict_id_display_type (DictIdType dict_id);

#endif
