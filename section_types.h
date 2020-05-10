// ------------------------------------------------------------------
//   section_types.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// note - separated from sections.h bc of circular include dependencies

#ifndef SECTION_TYPES_INCLUDED
#define SECTION_TYPES_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <stdbool.h>
#else
#include "compatibility/visual_c_stdint.h"
#include "compatibility/visual_c_stdbool.h"
#endif

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef enum {
    SEC_EOF                = -1, // doesn't appear in the file - just a value to indicate there are no more sections

    // data sections - statring in v1
    SEC_TXT_HEADER         = 0,                                // used by VCF, SAM, ME23
    SEC_VB_HEADER          = 1,                                // used by all data types

    SEC_VCF_FRMT_SF_DICT   = 2,  
    SEC_VCF_GT_DATA        = 3,  SEC_VCF_PHASE_DATA     = 4,
    SEC_HT_DATA            = 5,                                // used by VCF, ME23

    // data sections added in v2
    SEC_GENOZIP_HEADER     = 6,   SEC_RANDOM_ACCESS      = 7,

    SEC_CHROM_DICT         = 8,   SEC_CHROM_B250         = 9,  // used by VCF, ME23
    SEC_POS_DICT           = 10,  SEC_POS_B250           = 11, // used by VCF, ME23
    SEC_ID_DICT            = 12,  SEC_ID_B250            = 13, // used by VCF until v4 (VCF's ID moved to SEC_NUMERIC_ID_DATA in v5)
    SEC_VCF_REFALT_DICT    = 14,  SEC_VCF_REFALT_B250    = 15,  
    SEC_VCF_QUAL_DICT      = 16,  SEC_VCF_QUAL_B250      = 17, 
    SEC_VCF_FILTER_DICT    = 18,  SEC_VCF_FILTER_B250    = 19, 
    SEC_VCF_INFO_DICT      = 20,  SEC_VCF_INFO_B250      = 21, 
    SEC_VCF_FORMAT_DICT    = 22,  SEC_VCF_FORMAT_B250    = 23,
    SEC_VCF_INFO_SF_DICT   = 24,  SEC_VCF_INFO_SF_B250   = 25,
    
    // added gtshark support in file version 3
    SEC_HT_GTSHARK_DB_DB   = 26,  SEC_HT_GTSHARK_DB_GT   = 27,
    SEC_HT_GTSHARK_X_LINE  = 28,  SEC_HT_GTSHARK_X_HTI   = 29,
    SEC_HT_GTSHARK_X_ALLELE= 30,

    // added SAM, FASTQ, FASTA, GFF3, 23andMe support in file version 5
    SEC_RANDOM_POS_DATA    = 31,  SEC_SAM_MD_DATA        = 32, 
    SEC_SAM_BD_DATA        = 33,  SEC_SAM_BI_DATA        = 34,
    SEC_SEQ_DATA           = 35,                               // used by SAM, FASTA, FASTQ 
    SEC_QUAL_DATA          = 36,                               // used by SAM, FASTQ
    SEC_NUMERIC_ID_DATA    = 37,                               // used by VCF (starting v5) and ME23
    SEC_SAM_QNAME_SF_DICT  = 38,  SEC_SAM_QNAME_SF_B250  = 39,
    SEC_SAM_OPTNL_SF_DICT  = 40,  SEC_SAM_OPTNL_SF_B250  = 41,
    SEC_SAM_QNAME_DICT     = 42,  SEC_SAM_QNAME_B250     = 43,
    SEC_SAM_FLAG_DICT      = 44,  SEC_SAM_FLAG_B250      = 45,
    SEC_SAM_RNAME_DICT     = 46,  SEC_SAM_RNAME_B250     = 47, 
    SEC_SAM_POS_DICT       = 48,  SEC_SAM_POS_B250       = 49, 
    SEC_SAM_MAPQ_DICT      = 50,  SEC_SAM_MAPQ_B250      = 51, 
    SEC_SAM_CIGAR_DICT     = 52,  SEC_SAM_CIGAR_B250     = 53, 
    SEC_SAM_PNEXT_DICT     = 54,  SEC_SAM_PNEXT_B250     = 55, 
    SEC_SAM_TLEN_DICT      = 56,  SEC_SAM_TLEN_B250      = 57, 
    SEC_SAM_OPTIONAL_DICT  = 58,  SEC_SAM_OPTIONAL_B250  = 59, 

    SEC_FAST_DESC_SF_DICT  = 60,  SEC_FAST_DESC_SF_B250  = 61, // used by FASTQ & FASTA
    SEC_FAST_DESC_DICT     = 62,  SEC_FAST_DESC_B250     = 63, // used by FASTQ & FASTA
    SEC_FAST_LINEMETA_DICT = 64,  SEC_FAST_LINEMETA_B250 = 65, // used by FASTQ & FASTA
    SEC_FASTA_COMMENT_DATA = 66,

    SEC_GFF3_SEQID_DICT    = 67,  SEC_GFF3_SEQID_B250    = 68, 
    SEC_GFF3_SOURCE_DICT   = 69,  SEC_GFF3_SOURCE_B250   = 70,
    SEC_GFF3_TYPE_DICT     = 71,  SEC_GFF3_TYPE_B250     = 72, 
    SEC_GFF3_START_DICT    = 73,  SEC_GFF3_START_B250    = 74,
    SEC_GFF3_END_DICT      = 75,  SEC_GFF3_END_B250      = 76, 
    SEC_GFF3_SCORE_DICT    = 77,  SEC_GFF3_SCORE_B250    = 78,
    SEC_GFF3_STRAND_DICT   = 79,  SEC_GFF3_STRAND_B250   = 80, 
    SEC_GFF3_PHASE_DICT    = 81,  SEC_GFF3_PHASE_B250    = 82,
    SEC_GFF3_ATTRS_DICT    = 83,  SEC_GFF3_ATTRS_B250    = 84, 
    SEC_GFF3_ATTRS_SF_DICT = 85,  SEC_GFF3_ATTRS_SF_B250 = 86, 

    SEC_ENST_DATA          = 87,
    // This sections is not a real section - it doesn't appear in the genozip file. It can be changed if needed.
    SEC_STATS_HT_SEPERATOR
} SectionType;

// this data must be perfectly aligned with SectionType. it contains: 1. name 2. 1 if the section is stripped out in --strip 
#define SECTIONTYPE_ABOUT { \
    {"SEC_TXT_HEADER",          0},  {"SEC_VB_HEADER",          0},\
    \
    {"SEC_VCF_FRMT_SF_DICT",    1},  {"SEC_VCF_GT_DATA",        1},\
    {"SEC_VCF_PHASE_DATA",      0},  {"SEC_HT_DATA",            0},\
    \
    {"SEC_GENOZIP_HEADER",      0},  {"SEC_RANDOM_ACCESS",      0},\
    \
    {"SEC_CHROM_DICT",          0},  {"SEC_CHROM_B250",         0},\
    {"SEC_POS_DICT",            0},  {"SEC_POS_B250",           0},\
    {"SEC_ID_DICT",             1},  {"SEC_ID_B250",            1},\
    {"SEC_VCF_REFALT_DICT",     0},  {"SEC_VCF_REFALT_B250",    0},\
    {"SEC_VCF_QUAL_DICT",       1},  {"SEC_VCF_QUAL_B250",      1},\
    {"SEC_VCF_FILTER_DICT",     1},  {"SEC_VCF_FILTER_B250",    1},\
    {"SEC_VCF_INFO_DICT",       1},  {"SEC_VCF_INFO_B250",      1},\
    {"SEC_VCF_FORMAT_DICT",     1},  {"SEC_VCF_FORMAT_B250",    1},\
    {"SEC_VCF_INFO_SF_DICT",    1},  {"SEC_VCF_INFO_SF_B250",   1},\
    \
    {"SEC_HT_GTSHARK_DB_DB",    0},  {"SEC_HT_GTSHARK_DB_GT",   0},\
    {"SEC_HT_GTSHARK_X_LINE",   0},  {"SEC_HT_GTSHARK_X_HTI",   0},\
    {"SEC_HT_GTSHARK_X_ALLELE", 0},\
    \
    {"SEC_RANDOM_POS_DATA",     0},  {"SEC_SAM_MD_DATA",        1},\
    {"SEC_SAM_BD_DATA",         1},  {"SEC_SAM_BI_DATA",        1},\
    {"SEC_SEQ_DATA",            0},\
    {"SEC_QUAL_DATA",           1},\
    {"SEC_NUMERIC_ID_DATA",     1},\
    {"SEC_SAM_QNAME_SF_DICT",   1},  {"SEC_SAM_QNAME_SF_B250",  1},\
    {"SEC_SAM_OPTNL_SF_DICT",   0},  {"SEC_SAM_OPTNL_SF_B250",  0},\
    {"SEC_SAM_QNAME_DICT",      1},  {"SEC_SAM_QNAME_B250",     1},\
    {"SEC_SAM_FLAG_DICT",       1},  {"SEC_SAM_FLAG_B250",      1},\
    {"SEC_SAM_RNAME_DICT",      0},  {"SEC_SAM_RNAME_B250",     0},\
    {"SEC_SAM_POS_DICT",        0},  {"SEC_SAM_POS_B250",       0},\
    {"SEC_SAM_MAPQ_DICT",       1},  {"SEC_SAM_MAPQ_B250",      1},\
    {"SEC_SAM_CIGAR_DICT",      0},  {"SEC_SAM_CIGAR_B250",     0},\
    {"SEC_SAM_PNEXT_DICT",      1},  {"SEC_SAM_PNEXT_B250",     1},\
    {"SEC_SAM_TLEN_DICT",       1},  {"SEC_SAM_TLEN_B250",      1},\
    {"SEC_SAM_OPTIONAL_DICT",   0},  {"SEC_SAM_OPTIONAL_B250",  0},\
    \
    {"SEC_FAST_DESC_SF_DICT",   1},  {"SEC_FAST_DESC_SF_B250",  1},\
    {"SEC_FAST_DESC_DICT",      1},  {"SEC_FAST_DESC_B250",     1},\
    {"SEC_FAST_LINEMETA_DICT",  0},  {"SEC_FAST_LINEMETA_B250", 0},\
    {"SEC_FASTA_COMMENT_DATA",  1},\
    \
    {"SEC_GFF3_SEQID_DICT",     0}, {"SEC_GFF3_SEQID_B250",     0},\
    {"SEC_GFF3_SOURCE_DICT",    0}, {"SEC_GFF3_SOURCE_B250",    0},\
    {"SEC_GFF3_TYPE_DICT",      0}, {"SEC_GFF3_TYPE_B250",      0},\
    {"SEC_GFF3_START_DICT",     0}, {"SEC_GFF3_START_B250",     0},\
    {"SEC_GFF3_END_DICT",       0}, {"SEC_GFF3_END_B250",       0},\
    {"SEC_GFF3_SCORE_DICT",     0}, {"SEC_GFF3_SCORE_B250",     0},\
    {"SEC_GFF3_STRAND_DICT",    0}, {"SEC_GFF3_STRAND_B250",    0},\
    {"SEC_GFF3_PHASE_DICT",     0}, {"SEC_GFF3_PHASE_B250",     0},\
    {"SEC_GFF3_ATTRS_DICT",     0}, {"SEC_GFF3_ATTRS_B250",     0},\
    {"SEC_GFF3_ATTRS_SF_DICT",  0}, {"SEC_GFF3_ATTRS_SF_B250",  0},\
    \
    {"SEC_ENST_DATA",           0},\
    \
    {"SEC_STATS_HT_SEPERATOR",  0} \
}

#define NUM_SEC_TYPES (SEC_STATS_HT_SEPERATOR+1) 

#define section_type_is_dictionary(s) (((s) >= SEC_CHROM_DICT        && (s) <= SEC_VCF_INFO_SF_DICT   && ((s)%2) == (SEC_CHROM_DICT % 2)) ||       \
                                        (s) == SEC_VCF_FRMT_SF_DICT || \
                                       ((s) >= SEC_SAM_QNAME_SF_DICT && (s) <= SEC_FAST_LINEMETA_DICT && ((s)%2) == (SEC_SAM_QNAME_SF_DICT % 2)) ||\
                                       ((s) >= SEC_GFF3_SEQID_DICT   && (s) <= SEC_GFF3_ATTRS_SF_DICT && ((s)%2) == (SEC_GFF3_SEQID_DICT % 2)))

#define section_type_is_b250(s)       (section_type_is_dictionary((s)-1) && (s) != SEC_VCF_GT_DATA /* one exception */)

#define FIELD_TO_DICT_SECTION(dt,f) (dt_fields[dt].first_dict_sec + (f)*2)
#define FIELD_TO_B250_SECTION(dt,f) (FIELD_TO_DICT_SECTION((dt), (f)) + 1)

#endif
