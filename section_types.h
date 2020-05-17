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

    SEC_VCF_FRMT_SF_DICT   = 2,  
    SEC_VCF_GT_DATA        = 3,  SEC_VCF_PHASE_DATA     = 4,
    SEC_HT_DATA            = 5,                                // used by VCF, ME23

    // data sections added in v2
    SEC_GENOZIP_HEADER     = 6,   SEC_RANDOM_ACCESS      = 7,

    // these were used in v1-v4, no longer used for files created starting v5
    SEC_CHROM_DICT         = 8,   SEC_CHROM_B250         = 9,  
    SEC_POS_DICT           = 10,  SEC_POS_B250           = 11, 
    SEC_ID_DICT            = 12,  SEC_ID_B250            = 13, 
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
    SEC_NUMERIC_ID_DATA    = 37,                               // used by VCF (starting v5), GFF3 and ME23
    SEC_FASTA_COMMENT_DATA = 38,

    SEC_ENST_DATA          = 39,

    SEC_DICT = 40, SEC_B250 = 41, SEC_LOCAL = 42,
    NUM_SEC_TYPES // fake section for counting
} SectionType;

// this data must be perfectly aligned with SectionType.
#define SECTIONTYPE_ABOUT { \
    {"SEC_TXT_HEADER",          },  {"SEC_VB_HEADER",          },\
    \
    {"SEC_VCF_FRMT_SF_DICT",    },  {"SEC_VCF_GT_DATA",        },\
    {"SEC_VCF_PHASE_DATA",      },  {"SEC_HT_DATA",            },\
    \
    {"SEC_GENOZIP_HEADER",      },  {"SEC_RANDOM_ACCESS",      },\
    \
    {"SEC_CHROM_DICT",          },  {"SEC_CHROM_B250",         },\
    {"SEC_POS_DICT",            },  {"SEC_POS_B250",           },\
    {"SEC_ID_DICT",             },  {"SEC_ID_B250",            },\
    {"SEC_VCF_REFALT_DICT",     },  {"SEC_VCF_REFALT_B250",    },\
    {"SEC_VCF_QUAL_DICT",       },  {"SEC_VCF_QUAL_B250",      },\
    {"SEC_VCF_FILTER_DICT",     },  {"SEC_VCF_FILTER_B250",    },\
    {"SEC_VCF_INFO_DICT",       },  {"SEC_VCF_INFO_B250",      },\
    {"SEC_VCF_FORMAT_DICT",     },  {"SEC_VCF_FORMAT_B250",    },\
    {"SEC_VCF_INFO_SF_DICT",    },  {"SEC_VCF_INFO_SF_B250",   },\
    \
    {"SEC_HT_GTSHARK_DB_DB",    },  {"SEC_HT_GTSHARK_DB_GT",   },\
    {"SEC_HT_GTSHARK_X_LINE",   },  {"SEC_HT_GTSHARK_X_HTI",   },\
    {"SEC_HT_GTSHARK_X_ALLELE", },\
    \
    {"SEC_RANDOM_POS_DATA",     },  {"SEC_SAM_MD_DATA",        },\
    {"SEC_SAM_BD_DATA",         },  {"SEC_SAM_BI_DATA",        },\
    {"SEC_SEQ_DATA",            },\
    {"SEC_QUAL_DATA",           },\
    {"SEC_NUMERIC_ID_DATA",     },\
    {"SEC_FASTA_COMMENT_DATA",  },\
      \
    {"SEC_ENST_DATA",           },\
    {"SEC_DICT",                }, {"SEC_B250",                }, {"SEC_LOCAL",               },\
}

#define section_type_is_dictionary(s) ((s) == SEC_DICT          || (s) == SEC_VCF_FRMT_SF_DICT || (s) == SEC_CHROM_DICT    || (s) == SEC_POS_DICT || \
                                       (s) == SEC_ID_DICT       || (s) == SEC_VCF_REFALT_DICT  || (s) == SEC_VCF_QUAL_DICT || (s) == SEC_VCF_FILTER_DICT || \
                                       (s) == SEC_VCF_INFO_DICT || (s) == SEC_VCF_FORMAT_DICT  || (s) == SEC_VCF_INFO_SF_DICT)

#define section_type_is_b250(s)  ((s)==SEC_B250 || (section_type_is_dictionary((s)-1) && (s) != SEC_VCF_GT_DATA /* one exception */))

#endif
