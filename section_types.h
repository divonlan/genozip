// ------------------------------------------------------------------
//   section_types.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// note - separated from sections.h bc of circular include dependencies

#ifndef SECTION_TYPES_INCLUDED
#define SECTION_TYPES_INCLUDED

#include "genozip.h"

// note: the numbering of the sections cannot be modified, for backward compatibility
typedef enum {
    SEC_NONE           = -1, // doesn't appear in the file 

    SEC_RANDOM_ACCESS   = 0,
    SEC_DICT_ID_ALIASES = 1,
    SEC_REFERENCE       = 2,
    SEC_TXT_HEADER      = 3, 
    SEC_VB_HEADER       = 4,
    SEC_DICT            = 5, 
    SEC_GENOZIP_HEADER  = 6, // SEC_GENOZIP_HEADER remains 6 as in v2-v5 to be able to read old versions' genozip header
    SEC_B250            = 7, 
    SEC_LOCAL           = 8, 

    // vcf specific    
    SEC_VCF_GT_DATA     = 9,  
    SEC_VCF_PHASE_DATA  = 10,
    SEC_VCF_HT_DATA     = 11,                               
    SEC_VCF_HT_GTSHARK  = 12,

    NUM_SEC_TYPES // fake section for counting
} SectionType;

// this data must be perfectly aligned with SectionType.
#define SECTIONTYPE_ABOUT {  \
    {"SEC_RANDOM_ACCESS"},   \
    {"SEC_DICT_ID_ALIASES"}, \
    {"SEC_REFERENCE"},       \
    {"SEC_TXT_HEADER"},      \
    {"SEC_VB_HEADER"},       \
    {"SEC_DICT"},            \
    {"SEC_GENOZIP_HEADER"},  \
    {"SEC_B250"},            \
    {"SEC_LOCAL"},           \
    {"SEC_VCF_GT_DATA"},     \
    {"SEC_VCF_PHASE_DATA"},  \
    {"SEC_VCF_HT_GTSHARK"},  \
}

#endif
