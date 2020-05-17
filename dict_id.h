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

#define DICT_ID_LEN    ((int)sizeof(uint64_t))    // VCF/SAM spec don't limit the ID length, we limit it to 8 chars. zero-padded. (note: if two fields have the same 8-char prefix - they will just share the same dictionary)
typedef union DictIdType {
    uint64_t num;            // num is just for easy comparisons - it doesn't have a numeric value and endianity should not be changed
    uint8_t id[DICT_ID_LEN]; // \0-padded IDs 
    uint16_t map_key;        // we use the first two bytes as they key into vb/z_file->dict_id_mapper
} DictIdType;

#pragma pack(pop)

extern DictIdType dict_id_make (const char *str, unsigned str_len);
static inline DictIdType dict_id_field (DictIdType dict_id) { dict_id.id[0] = dict_id.id[0] & 0x3f; return dict_id; } // set 2 Msb to 00

#define dict_id_is(dict_id, str) (dict_id_make (str, strlen(str)).num == dict_id_printable (dict_id).num)
static inline bool dict_id_is_field (DictIdType dict_id) { return ((dict_id.id[0] >> 6) == 0); } // 2 MSb of first byte determine dictionary type
static inline bool dict_id_is_type_1(DictIdType dict_id) { return ((dict_id.id[0] >> 6) == 3); }
static inline bool dict_id_is_type_2(DictIdType dict_id) { return ((dict_id.id[0] >> 6) == 1); }

static inline DictIdType dict_id_type_1(DictIdType dict_id) { dict_id.id[0] = dict_id.id[0] | 0xc0; return dict_id; } // set 2 Msb to 11
static inline DictIdType dict_id_type_2(DictIdType dict_id) { return dict_id; } // no change - keep Msb 01

// VCF field types -
#define dict_id_is_vcf_info_sf   dict_id_is_type_1
#define dict_id_is_vcf_format_sf dict_id_is_type_2

#define dict_id_vcf_info_sf      dict_id_type_1
#define dict_id_vcf_format_sf    dict_id_type_2

// SAM field types 
#define dict_id_is_sam_qname_sf  dict_id_is_type_1
#define dict_id_is_sam_optnl_sf  dict_id_is_type_2

#define dict_id_sam_qname_sf     dict_id_type_1
#define dict_id_sam_optnl_sf     dict_id_type_2

// FASTQ/FASTA field types 
#define dict_id_is_fast_desc_sf dict_id_is_type_2
#define dict_id_fast_desc_sf dict_id_type_2

// GFF3 field types
#define dict_id_is_gff3_attr_sf dict_id_is_type_1
#define dict_id_gff3_attr_sf dict_id_type_1

static inline DictIdType dict_id_printable(DictIdType dict_id) { dict_id.id[0] = (dict_id.id[0] & 0x7f) | 0x40; return dict_id; } // set 2 Msb to 01
#define DICT_ID_NONE ((DictIdType)0ULL)

extern DictIdType dict_id_show_one_b250, dict_id_show_one_dict; // arguments of --show-b250-one and --show-dict-one (defined in genozip.c)
extern DictIdType dict_id_dump_one_b250;                        // arguments of --dump-b250-one (defined in genozip.c)

extern uint64_t dict_id_fields[MAX_NUM_FIELDS_PER_DATA_TYPE],
                
                dict_id_FORMAT_PL, dict_id_FORMAT_GL, dict_id_FORMAT_GP, dict_id_FORMAT_DP, dict_id_FORMAT_MIN_DP, // some VCF FORMAT subfields
                dict_id_INFO_AC,  dict_id_INFO_AF, dict_id_INFO_AN, dict_id_INFO_DP, dict_id_INFO_VQSLOD, // some VCF INFO subfields
                dict_id_INFO_END, dict_id_WindowsEOL,

                // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
                dict_id_OPTION_AM, dict_id_OPTION_AS, dict_id_OPTION_CM, dict_id_OPTION_E2, dict_id_OPTION_LB, dict_id_OPTION_FI, 
                dict_id_OPTION_H0, dict_id_OPTION_H1, dict_id_OPTION_H2, dict_id_OPTION_MQ, dict_id_OPTION_NH, dict_id_OPTION_NM, 
                dict_id_OPTION_OA, dict_id_OPTION_OC, dict_id_OPTION_PG, dict_id_OPTION_PQ, dict_id_OPTION_PU, dict_id_OPTION_RG, 
                dict_id_OPTION_SA, dict_id_OPTION_SM, dict_id_OPTION_TC, dict_id_OPTION_U2, dict_id_OPTION_UQ, dict_id_OPTION_CC, 
                dict_id_OPTION_MC, dict_id_OPTION_MD,
    
                // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
                dict_id_OPTION_X0, dict_id_OPTION_X1, dict_id_OPTION_XA, dict_id_OPTION_XN, dict_id_OPTION_XM, dict_id_OPTION_XO,
                dict_id_OPTION_XG, dict_id_OPTION_XS, dict_id_OPTION_XE,

                // biobambam tags 
                dict_id_OPTION_mc, dict_id_OPTION_ms,

                // GATK tags
                dict_id_OPTION_BD, dict_id_OPTION_BI,
                
                // our own
                dict_id_OPTION_STRAND,

                // GVF attributes - standard
                dict_id_ATTR_ID, dict_id_ATTR_Variant_seq, dict_id_ATTR_Reference_seq, dict_id_ATTR_Variant_freq,

                // GVF attributes - from GRCh37/38 etc
                dict_id_ATTR_Dbxref, // example: "dbSNP_151:rs1282280967"
                dict_id_ATTR_ancestral_allele,
                dict_id_ATTR_Variant_effect, // example: "Variant_effect=non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
                dict_id_ATTR_sift_prediction,
                dict_id_ATTR_polyphen_prediction,
                dict_id_ATTR_variant_peptide,

                dict_id_ATTR_SEQ, dict_id_ENSTid;  // private genozip dict

extern void dict_id_initialize (DataType data_type);

// template can be 0 - anything OR a type - must 2 MSb of id[0] are used OR a specific dict_id
// candidate is a specific dict_id that we test for its matching of the template
extern bool dict_id_is_match (DictIdType template, DictIdType candidate);

extern const char *dict_id_display_type (DataType dt, DictIdType dict_id);

// print the dict_id - NOT thread safe, for use in execution-termination messages
extern const char *err_dict_id (DictIdType dict_id);

#endif
