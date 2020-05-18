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
    SEC_NONE               = -1, // doesn't appear in the file 

    // data sections - statring in v1
    SEC_TXT_HEADER         = 0,                                // used by VCF, SAM, ME23
    SEC_VB_HEADER          = 1,                                // used by all data types

    SEC_VCF_FRMT_SF_DICT_legacy   = 2,  
    
    SEC_VCF_GT_DATA        = 3,  SEC_VCF_PHASE_DATA     = 4,
    SEC_HT_DATA            = 5,                                // used by VCF, ME23

    // data sections added in v2
    SEC_GENOZIP_HEADER     = 6,   SEC_RANDOM_ACCESS      = 7,

    // these (2 + 8-25) were used in v1-v4, no longer used for files created starting v5 
    // we keep them here to be able to decompress files compressed with older versions
    SEC_CHROM_DICT_legacy         = 8,   SEC_CHROM_B250_legacy         = 9,  
    SEC_POS_DICT_legacy           = 10,  SEC_POS_B250_legacy           = 11, 
    SEC_ID_DICT_legacy            = 12,  SEC_ID_B250_legacy            = 13, 
    SEC_VCF_REFALT_DICT_legacy    = 14,  SEC_VCF_REFALT_B250_legacy    = 15,  
    SEC_VCF_QUAL_DICT_legacy      = 16,  SEC_VCF_QUAL_B250_legacy      = 17, 
    SEC_VCF_FILTER_DICT_legacy    = 18,  SEC_VCF_FILTER_B250_legacy    = 19, 
    SEC_VCF_INFO_DICT_legacy      = 20,  SEC_VCF_INFO_B250_legacy      = 21, 
    SEC_VCF_FORMAT_DICT_legacy    = 22,  SEC_VCF_FORMAT_B250_legacy    = 23,
    SEC_VCF_INFO_SF_DICT_legacy   = 24,  SEC_VCF_INFO_SF_B250_legacy   = 25,
    
    // added gtshark support in file version 3
    SEC_HT_GTSHARK_DB_DB   = 26,  SEC_HT_GTSHARK_DB_GT   = 27,
    SEC_HT_GTSHARK_X_LINE  = 28,  SEC_HT_GTSHARK_X_HTI   = 29,
    SEC_HT_GTSHARK_X_ALLELE= 30,

    SEC_DICT = 31, SEC_B250 = 32, SEC_LOCAL = 33,
    NUM_SEC_TYPES // fake section for counting
} SectionType;

// this data must be perfectly aligned with SectionType.
#define SECTIONTYPE_ABOUT { \
    {"SEC_TXT_HEADER",          },  {"SEC_VB_HEADER",          },\
    \
    {"SEC_VCF_FRMT_SF_DICT_legacy",    },  {"SEC_VCF_GT_DATA",        },\
    {"SEC_VCF_PHASE_DATA",      },  {"SEC_HT_DATA",            },\
    \
    {"SEC_GENOZIP_HEADER",      },  {"SEC_RANDOM_ACCESS",      },\
    \
    {"SEC_CHROM_DICT_legacy",          },  {"SEC_CHROM_B250_legacy",         },\
    {"SEC_POS_DICT_legacy",            },  {"SEC_POS_B250_legacy",           },\
    {"SEC_ID_DICT_legacy",             },  {"SEC_ID_B250_legacy",            },\
    {"SEC_VCF_REFALT_DICT_legacy",     },  {"SEC_VCF_REFALT_B250_legacy",    },\
    {"SEC_VCF_QUAL_DICT_legacy",       },  {"SEC_VCF_QUAL_B250_legacy",      },\
    {"SEC_VCF_FILTER_DICT_legacy",     },  {"SEC_VCF_FILTER_B250_legacy",    },\
    {"SEC_VCF_INFO_DICT_legacy",       },  {"SEC_VCF_INFO_B250_legacy",      },\
    {"SEC_VCF_FORMAT_DICT_legacy",     },  {"SEC_VCF_FORMAT_B250_legacy",    },\
    {"SEC_VCF_INFO_SF_DICT_legacy",    },  {"SEC_VCF_INFO_SF_B250_legacy",   },\
    \
    {"SEC_HT_GTSHARK_DB_DB",    },  {"SEC_HT_GTSHARK_DB_GT",   },\
    {"SEC_HT_GTSHARK_X_LINE",   },  {"SEC_HT_GTSHARK_X_HTI",   },\
    {"SEC_HT_GTSHARK_X_ALLELE", },\
    \
    {"SEC_DICT",                }, {"SEC_B250",                }, {"SEC_LOCAL",               },\
}

#define section_type_is_dictionary(s) ((s) == SEC_DICT          || (s) == SEC_VCF_FRMT_SF_DICT_legacy || (s) == SEC_CHROM_DICT_legacy    || (s) == SEC_POS_DICT_legacy || \
                                       (s) == SEC_ID_DICT_legacy       || (s) == SEC_VCF_REFALT_DICT_legacy  || (s) == SEC_VCF_QUAL_DICT_legacy || (s) == SEC_VCF_FILTER_DICT_legacy || \
                                       (s) == SEC_VCF_INFO_DICT_legacy || (s) == SEC_VCF_FORMAT_DICT_legacy  || (s) == SEC_VCF_INFO_SF_DICT_legacy)

#define section_type_is_b250(s)  ((s)==SEC_B250 || (section_type_is_dictionary((s)-1) && (s) != SEC_VCF_GT_DATA /* one exception */))

#endif
