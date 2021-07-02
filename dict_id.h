// ------------------------------------------------------------------
//   dict_id.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DICT_ID_INCLUDED
#define DICT_ID_INCLUDED

#include <inttypes.h>
#include "genozip.h"
#include "data_types.h"

typedef enum { DTYPE_FIELD, DTYPE_1, DTYPE_2 } DictIdType;
#define DTYPE_PLAIN DTYPE_2
extern DictId dict_id_make (const char *str, unsigned str_len, DictIdType dict_id_type);

#define dict_id_is(dict_id, str) (dict_id_make (str, strlen(str)).num == dict_id_typeless (dict_id).num)
static inline bool dict_id_is_field (DictId dict_id) { return ((dict_id.id[0] >> 6) == 0); } // 2 MSb of first byte determine dictionary type
static inline bool dict_id_is_type_1(DictId dict_id) { return ((dict_id.id[0] >> 6) == 3); }
static inline bool dict_id_is_type_2(DictId dict_id) { return ((dict_id.id[0] >> 6) == 1); }

static inline DictId dict_id_typeless(DictId dict_id) { dict_id.id[0] = (dict_id.id[0] & 0x7f) | 0x40; return dict_id; } // set 2 Msb to 01

#define DICT_ID_NONE ((DictId)(uint64_t)0)

typedef struct { DictId alias, dst; } DictIdAlias;
extern const DictIdAlias *dict_id_aliases;
uint32_t dict_id_num_aliases;

extern BufferP dict_id_create_aliases_buf (void);
extern void dict_id_read_aliases (void) ;

extern uint64_t dict_id_fields[MAX_NUM_FIELDS_PER_DATA_TYPE],
                
                dict_id_FORMAT_PL, dict_id_FORMAT_GL, dict_id_FORMAT_GP, dict_id_FORMAT_DP, dict_id_FORMAT_AF, // some VCF FORMAT subfields
                dict_id_FORMAT_PS, dict_id_FORMAT_GT, dict_id_FORMAT_PRI, dict_id_FORMAT_PP,
                dict_id_FORMAT_GT_HT, dict_id_FORMAT_GT_HT_INDEX,
                dict_id_PBWT_RUNS, dict_id_PBWT_FGRC,
                dict_id_FORMAT_AD, dict_id_FORMAT_ADF, dict_id_FORMAT_ADR, dict_id_FORMAT_ADALL, 
                dict_id_FORMAT_GQ, dict_id_FORMAT_DS, dict_id_FORMAT_SAC, dict_id_FORMAT_SB, dict_id_FORMAT_MB, 
                dict_id_INFO_AC,  dict_id_INFO_AF, dict_id_INFO_AN, dict_id_INFO_DP, dict_id_INFO_AA, dict_id_INFO_VQSLOD, // some VCF INFO subfields
                dict_id_INFO_DP4, dict_id_INFO_SF, dict_id_INFO_SVLEN, dict_id_WindowsEOL,
                dict_id_INFO_LUFT, dict_id_INFO_PRIM, dict_id_INFO_LREJ, dict_id_INFO_PREJ, 
                dict_id_INFO_BaseCounts, dict_id_INFO_MAX_AF,

                // see: https://support.illumina.com/help/BS_App_DRAGEN_Enrichment_OLH_1000000095374/Content/Source/Informatics/Apps/VCFAnnotations_swBS_appDNAA_appDRNA_appDRAGE_appDRAGGP.htm
                dict_id_FORMAT_F1R2, dict_id_FORMAT_F2R1,

                // tags from HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
                dict_id_INFO_END, dict_id_INFO_MLEAC, dict_id_INFO_MLEAF, dict_id_INFO_MQ0, dict_id_INFO_LDAF,
                dict_id_FORMAT_MIN_DP, 
                
                // tags from VEP (Varient Effect Predictor) and similar tools
                dict_id_INFO_CSQ, dict_id_INFO_vep, dict_id_INFO_DP_HIST, dict_id_INFO_GQ_HIST, 
                dict_id_INFO_AGE_HISTOGRAM_HET, dict_id_INFO_AGE_HISTOGRAM_HOM,

                // standard tags, see here: https://samtools.github.io/hts-specs/SAMtags.pdf
                dict_id_OPTION_AM, dict_id_OPTION_AS, dict_id_OPTION_CM, dict_id_OPTION_E2, dict_id_OPTION_LB, dict_id_OPTION_FI, 
                dict_id_OPTION_H0, dict_id_OPTION_H1, dict_id_OPTION_H2, dict_id_OPTION_MQ, dict_id_OPTION_NH, dict_id_OPTION_NM, 
                dict_id_OPTION_OA, dict_id_OPTION_OC, dict_id_OPTION_PG, dict_id_OPTION_PQ, dict_id_OPTION_PU, dict_id_OPTION_RG, 
                dict_id_OPTION_SA, dict_id_OPTION_SM, dict_id_OPTION_TC, dict_id_OPTION_U2, dict_id_OPTION_UQ, dict_id_OPTION_CC, 
                dict_id_OPTION_MC, dict_id_OPTION_MD,
    
                // bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
                dict_id_OPTION_X0, dict_id_OPTION_X1, dict_id_OPTION_XA, dict_id_OPTION_XA_RNAME, dict_id_OPTION_XN, dict_id_OPTION_XM, dict_id_OPTION_XO,
                dict_id_OPTION_XG, dict_id_OPTION_XS, dict_id_OPTION_XE,

                // 
                dict_id_OPTION_ZM,

                // biobambam tags 
                dict_id_OPTION_mc, dict_id_OPTION_ms,

                // GATK tags
                dict_id_OPTION_BD, dict_id_OPTION_BI, dict_id_OPTION_BD_BI,
                
                // our own
                dict_id_OPTION_STRAND, dict_id_OPTION_RNAME, dict_id_OPTION_POS, dict_id_OPTION_CIGAR, dict_id_OPTION_MAPQ,
                dict_id_OPTION_TX,
                
                // GVF attributes - standard
                dict_id_ATTR_ID, dict_id_ATTR_Variant_seq, dict_id_ATTR_Reference_seq, dict_id_ATTR_Variant_freq,

                // GVF attributes - from GRCh37/38 etc
                dict_id_ATTR_Dbxref, // example: "dbSNP_151:rs1282280967"
                dict_id_ATTR_ancestral_allele,
                dict_id_ATTR_Variant_effect, // example: "Variant_effect=non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
                dict_id_ATTR_sift_prediction,
                dict_id_ATTR_polyphen_prediction,
                dict_id_ATTR_variant_peptide,

                dict_id_ENSTid;  // private genozip dict

extern void dict_id_initialize (DataType data_type);

// template can be 0 - anything OR a type - must 2 MSb of id[0] are used OR a specific dict_id
// candidate is a specific dict_id that we test for its matching of the template
extern bool dict_id_is_match (DictId template, DictId candidate);

extern const char *dict_id_display_type (DataType dt, DictId dict_id);

typedef struct { char s[20]; } DisplayPrintId;
extern DisplayPrintId dis_dict_id_ex (DictId dict_id, bool with_type_if_vcf);
#define dis_dict_id(dict_id) dis_dict_id_ex ((dict_id), false)
#define dis_dict_id_name(dict_id) dis_dict_id_ex ((dict_id), true) // display with FORMAT/ or INFO/ if VCF/BCF
#endif
