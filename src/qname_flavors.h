// ------------------------------------------------------------------
//   qname_flavors.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once 

#include "container.h"
#include "dict_id_gen.h"
#include "segconf.h"
#include "aliases.h"
#include "qname.h" // for editor

#define PX_MATE_FIXED_0_PAD (char[]){ CI0_SKIP } // add to prefix for mate item IFF preceeding item has CI0_FIXED_0_PAD
#define I_AM_MATE .separator = { CI0_SKIP }

//-----------------------------------------------------------------------------------------
// Illumina-7 format: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  
// Example: HWI-D00360:5:H814YADXX:1:1106:10370:52569
// Example FASTQ: A00488:61:HMLGNDSXX:4:1101:1940:1000 2:N:0:CTGAAGCT+ATAGAGGC
// @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
// "sample number" can sometimes be the sequence of the index read instead of a number
// See here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
// interpretation of the first field here: https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py#L12-L45
//-----------------------------------------------------------------------------------------

static SmallContainer con_illumina_7 = {
    .repeats   = 1,
    .nitems_lo = 8,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: G10321:222:ASCASDFDA:1:2299:15331:1995#CTGGGAAG
static SmallContainer con_illumina_7i = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = "#"          },
                   { .dict_id = { _SAM_Q7NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: A00488:61:HMLGNDSXX:4:1101:4345:1000;umi=ACCTTCCAA
static SmallContainer con_illumina_7umi = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ";"          },
                   { .dict_id = { _SAM_Q7NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

#define PX_illumina_7umi { "", "", "", "", "", "", "", "umi=" }

// Example: A00180:28:HC3F5DRXX:2:2110:27453:21981_1:N:0:ATTACTCGATCT+GGCTCTGA
// Observed as QNAME2 in an SRR FASTQ QNAME and in SAM
static SmallContainer con_illumina_X_2bc = {
    .repeats   = 1,
    .nitems_lo = 13,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = "_"          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q8NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q9NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_QANAME }, .separator = "+"          },
                   { .dict_id = { _SAM_QBNAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: ST-E00314:354:H7J2YCCXY:1:1101:7080:1450_1:N:0:NAGGCG
// Observed as QNAME2 in an SRR FASTQ QNAME and in SAM
static SmallContainer con_illumina_X_1bc = {
    .repeats   = 1,
    .nitems_lo = 12,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = "_"          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q8NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q9NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_QANAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: A00488:61:HMLGNDSXX:4:1101:4345:1000_2:N:0
static SmallContainer con_illumina_X_0bc = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = "_"          },
                   { .dict_id = { _SAM_Q7NAME },                           }, // entire part after underscore is segged together 
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// Example: A01214:45:HF2GTDSX2:1:1405:25997:27383 1:N:0:GAACGCAATA+ACAGTAAGAT
static SmallContainer con_illumina_S_2bc = {
    .repeats   = 1,
    .nitems_lo = 13,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = " "          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q8NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q9NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_QANAME }, .separator = "+"          },
                   { .dict_id = { _SAM_QBNAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: ST-E00314:354:H7J2YCCXY:1:1101:7080:1450 1:N:0:NAGGCG
static SmallContainer con_illumina_S_1bc = {
    .repeats   = 1,
    .nitems_lo = 12,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = " "          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q8NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q9NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_QANAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: A00488:61:HMLGNDSXX:4:1101:4345:1000_2:N:0
static SmallContainer con_illumina_S_0bc = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = " "          },
                   { .dict_id = { _SAM_Q7NAME },                           }, // entire part after underscore is segged together 
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// Example: SDF-02:GFH-0166::1:13435:2311:1233:GTAGCCAATCA
static SmallContainer con_illumina_7bc = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q7NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: A00488:61:HMLGNDSXX:4:1101:4345:1000:TGCTGGG+ACTTTTA
static SmallContainer con_illumina_7bc2 = {
    .repeats   = 1,
    .nitems_lo = 10,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = "+"          },
                   { .dict_id = { _SAM_Q8NAME },                           }, 
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: ATATA-ATGCATAG|ab|A00488:61:HMLGNDSXX:4:1101:4345:1000|1 (BAM only)
static SmallContainer con_illumina_7gs = {
    .repeats   = 1,
    .nitems_lo = 10,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = "|"          }, // the two parts of the barcode are correletated and hence segged together
                   { .dict_id = { _SAM_Q1NAME }, .separator = "|"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q8NAME }, .separator = "|"          },
                   { .dict_id = { _SAM_Q9NAME },                           } } // number of reads merged together during the consensus analysis step 
};

// Example: ATATA-ATGCATAG|ab|A00488:61:HMLGNDSXX:4:1101:4345:1000 (FASTQ only)
static SmallContainer con_illumina_7gsFQ = {
    .repeats   = 1,
    .nitems_lo = 10,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = "|"          },
                   { .dict_id = { _SAM_Q1NAME }, .separator = "|"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q8NAME }                            }, 
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } } 
};

// Example: 2:N:0:CTGAAGCT+ATAGAGGC
static SmallContainer con_illumina_2bc = {
    .repeats   = 1,
    .nitems_lo = 5,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q3NAME }, .separator = "+"                    },  // first barcode
                   { .dict_id = { _SAM_Q4NAME },                                     } } // second barcode
};

#define PX_illumina_2bc { "", "", ":", "", "" } 

// Example: "1:N:0:0" 
static SmallContainer con_illumina_0bc = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q3NAME }                                      } }  // barcode
};

#define PX_illumina_0bc { "", "", ":", "" } 

// Example: 2:N:0:CTGAAGCT
static SmallContainer con_illumina_1bc = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q3NAME }                                      } }  // barcode
};

#define PX_illumina_1bc { "", "", ":", "" } 

//--------------------------------------------------------------------------------------------------------------
// BGI_E format: 27 characters: E, FlowCellSerialNumber[various], L, Lane[1], C, Column[3], R, Row[3] Tile[6-8]
// Example: E100020409L1C001R0030000234
// See: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
//
// Examples:
// R6 @8A_V100004684L3C001R029011637/1   https://db.cngb.org/search/experiment/CNX0094951/
// R6 @V300014293BL2C001R027005967/1     https://db.cngb.org/search/experiment/CNX0045105/  https://db.cngb.org/search/experiment/CNX0048121/ (observed A or B before L)
// R7 @V300017009_8AL2C001R0030001805/1  https://db.cngb.org/search/experiment/CNX0094808/
// R8 @V300046476L1C001R00100001719/1    https://db.cngb.org/search/experiment/CNX0142208/
// R7 @V300022116L2C001R0010002968/1     https://db.cngb.org/search/experiment/CNX0084815/
// R6 @V300003413L4C001R016000000/1      https://db.cngb.org/search/experiment/CNX0094811/
// R7 @V300014296L2C001R0010000027/1     https://db.cngb.org/search/experiment/CNX0094807/
// R7 @E100001117L1C001R0030000000/1
// R7 @E1000536L1C002R0020000005/1       https://db.cngb.org/search/experiment/CNX0145973/
//--------------------------------------------------------------------------------------------------------------
#define CON_BGI_R(n) \
static SmallContainer con_bgi_R##n = {  \
    .repeats             = 1,           \
    .nitems_lo           = 6,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_BGI_R(6);
CON_BGI_R(7);
CON_BGI_R(8);

#define PX_bgi_R { "", "", "C", "R", "", PX_MATE_FIXED_0_PAD }

#define CON_BGI_Rgs(n) /* BAM only*/ \
static SmallContainer con_bgi_Rgs##n = {  \
    .repeats             = 1,           \
    .nitems_lo           = 8,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "|"                    }, /* the two parts of the barcode are correletated and hence segged together */                \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "|"                    },                 \
                             { .dict_id = { _SAM_Q2NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_Q7NAME }                                      } }/* number of reads merged together during the consensus analysis step */ \
}
CON_BGI_Rgs(8);

#define PX_bgi_Rgs { "", "", "", "", "C", "R", "", "|" }

// Example: CGGTCT-AACCT|ab|E200003777L1C001R00100888074
#define CON_BGI_RgsFQ(n) /* FASTQ only */ \
static SmallContainer con_bgi_RgsFQ##n = {  \
    .repeats             = 1,           \
    .nitems_lo           = 8,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "|"                    }, /* the two parts of the barcode are correletated and hence segged together */                \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "|"                    },                 \
                             { .dict_id = { _SAM_Q2NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_BGI_RgsFQ(8);

#define PX_bgi_RgsFQ { "", "", "", "", "C", "R", "", PX_MATE_FIXED_0_PAD }

// variant of BGI flavor where Q4NAME is a variable-length integer rather than a fixed-length zero-padded numeric
static SmallContainer con_bgi_varlen = {  \
    .repeats             = 1,           \
    .nitems_lo           = 6,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME },                                     }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
};

#define PX_bgi_varlen { "", "", "C", "R" }

//--------------------------------------------------------------------------------------------------------------
// Same as CON_BGI_R, but the first separator is two L
//  DP8400010271TLL1C005R0511863479

#define CON_BGI_LL(n) \
static SmallContainer con_bgi_LL##n = {  \
    .repeats             = 1,            \
    .nitems_lo           = 6,            \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { 'L', 'L'           } }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_BGI_LL(7); // only encountered with 7 so far

#define PX_bgi_LL PX_bgi_R

//--------------------------------------------------------------------------------------------------------------------
// BGI_CL format: CL, FlowCellSerialNumber[9], L, Lane[1], C, Column[3], R, Row[3], _, variable-length number
//  CL100025298L1C002R050_244547 reported by https://en.wikipedia.org/wiki/File:BGI_seq_platform_read_name_description.png
// @CL100072652L2C001R001_12/1 for FASTQ reported by https://www.biostars.org/p/395715/
//--------------------------------------------------------------------------------------------------------------------
static SmallContainer con_bgi_CL = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "L" },                     /* Flow cell */
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  /* Lane      */
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  /* Column    */
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  /* Row       */
                             { .dict_id = { _SAM_Q4NAME },                                     },  /* Tile      */
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } /* Mate      */
};

#define PX_bgi_CL { "CL", "", "C", "R", "_" } // "CL" is a prefix, and the 2nd L is seprator 

// ULTIMA_1 format
// Example: 004733_1-X0003-1345904491
static SmallContainer con_ultima_a = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } } } 
};

#define PX_ULTIMA_A { "", "_", "", "-" }  

// Example: 004733_1-X0003-0072646116_TCGTCACTCGAAAACT
static SmallContainer con_ultima_a_bc = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_Q4NAME }                                       } } 
};

#define PX_ULTIMA_A_BC { "", "_", "", "-", "_" }  

// Example: V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:
static SmallContainer con_ultima_c = {
    .repeats   = 1,
    .nitems_lo = 14,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q8NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q9NAME }, .separator = ":" },
                   { .dict_id = { _SAM_QANAME }, .separator = ":" },
                   { .dict_id = { _SAM_QBNAME }, .separator = ":" },
                   { .dict_id = { _SAM_QCNAME }, .separator = ":" },
                   { .dict_id = { _SAM_QDNAME }, .separator = ":" } }
};

// Example: V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:_GCTGCTGACA
static SmallContainer con_ultima_c_bc = {
    .repeats   = 1,
    .nitems_lo = 15,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"  },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q7NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q8NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_Q9NAME }, .separator = ":"  },
                   { .dict_id = { _SAM_QANAME }, .separator = ":"  },
                   { .dict_id = { _SAM_QBNAME }, .separator = ":"  },
                   { .dict_id = { _SAM_QCNAME }, .separator = ":"  },
                   { .dict_id = { _SAM_QDNAME }, .separator = ":_" },
                   { .dict_id = { _SAM_QENAME }                    } }
};

// QF_ULTIMA_b format
// Example: 014214_2-UFOx1-143-9871314132
static SmallContainer con_ultima_b = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 10 } } } 
};

#define PX_ULTIMA_B { "", "_", "", "-", "" }  

// QF_ULTIMA_b_bc format
// Example: 014214_2-UFOx1-143-9871314132_AGTAC
static SmallContainer con_ultima_b_bc = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_Q5NAME }                                       } } 
};

#define PX_ULTIMA_B_BC { "", "_", "", "-", "", "_" }  

// QF_ULTIMA_d format
// Example: 014214-UFOx1-143-9871314132
static SmallContainer con_ultima_d = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } } } 
};

#define PX_ULTIMA_D { "", "-", "-", "" }  

// QF_ULTIMA_d_bc format
// Example: 014214-UFOx1-143-9871314132_AGTAC
static SmallContainer con_ultima_d_bc = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_Q4NAME }                                       } } 
};

#define PX_ULTIMA_D_BC { "", "-", "-", "", "_" }  

// Singular Genomics
// example: B05:000:FC2:4:1:272670:483
static SmallContainer con_singular = {
    .repeats   = 1,
    .nitems_lo = 8,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q6NAME },                                     },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } 
};

#define PX_SINGULAR { "", "", ":" }  

// Element Biosciences
// examples: 
// PLT-03:BBS-0174:2140948523:1:10102:0293:0058
// PLT-16:APP-0316:UNKNOWN_FLOWCELL:1:10102:0582:0027
static SmallContainer con_element = {
    .repeats   = 1,
    .nitems_lo = 8,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    }, 
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"                    }, 
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"                    }, 
                   { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } 
};

#define PX_ELEMENT { "", "", "", "", "", "", ":", PX_MATE_FIXED_0_PAD }  

// example: PLT-04:APP-0289::1:10102:0562:0050
static SmallContainer con_element_6 = {
    .repeats   = 1,
    .nitems_lo = 7,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = {':',':'}              },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    }, 
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"                    }, 
                   { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } 
};

#define PX_ELEMENT_6 { "", "", "", "", "", ":", PX_MATE_FIXED_0_PAD }  

//--------------------------------------------------------------------------------------------------------------
// ION_TORRENT_3 format: 17 characters: RunID[5] : Row[5] : Column[5]
// Example: ZEWTM:10130:07001
// See: https://www.biostars.org/p/330644/
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_ion_torrent_3 = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":" },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":" },
                   { .dict_id = { _SAM_Q2NAME }                   },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE        } }
};
// static SmallContainer con_ion_torrent_3 = {
//     .repeats   = 1,
//     .nitems_lo = 4,
//     .items     = { { .dict_id = { _SAM_Q0NAME },                                     },  
//                    { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
//                    { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
//                    { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }
// };

// #define PX_ion_torrent_3 { "", ":", ":", "" } 

//--------------------------------------------------------------------------------------------------------------
// Illumina-5 (old, including Solexa) format: <machine_id>:<lane>:<tile>:<x_coord>:<y_coord> 
// Example: SOLEXA-1GA-1_4_FC20ENL:7:258:737:870 (all fields are variable-length)
// Example: HWI-ST550_0201:3:1101:1626:2216#ACAGTG
// See: https://www.biostars.org/p/330644/
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_illumina_5 = {
    .repeats   = 1,
    .nitems_lo = 6,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// HWI-ST550_0201:3:1101:1626:2216#ACAGTG
static SmallContainer con_illumina_5i = {
    .repeats   = 1,
    .nitems_lo = 7,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = "#"          },
                   { .dict_id = { _SAM_Q5NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// NOVID_3053_FC625AGAAXX:6:1:1069:11483:0,84
static SmallContainer con_illumina_5rng = {
    .repeats   = 1,
    .nitems_lo = 8,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ","          }, // 5,6 appear to be a range [start_base,end_base)
                   { .dict_id = { _SAM_Q6NAME },                           }, // we speculatively set this as seq_len as normally Q5NAME is 0
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

//--------------------------------------------------------------------------------------------------------------------
// ROCHE_454 format: 
// Example: 000050_1712_0767
// See: https://www.sequencing.uio.no/illumina-services/4.%20Data%20delivery/results_454.html
//--------------------------------------------------------------------------------------------------------------------
static SmallContainer con_roche_454 = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }
};

#define PX_roche_454 { "", "_", "_", PX_MATE_FIXED_0_PAD }

// See: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats
// example: VHE-242383071011-15-1-0-2
static SmallContainer con_helicos = {
    .repeats             = 1,
    .nitems_lo           = 7,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"          },
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-"          },
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-"          },
                             { .dict_id = { _SAM_Q4NAME }, .separator = "-"          },
                             { .dict_id = { _SAM_Q5NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

static SmallContainer con_pacbio_3 = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 8 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "_"                    },
                             { .dict_id = { _SAM_Q2NAME }                                      },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }
};

#define PX_pacbio_3 { "", "_" }


// {movieName}/{holeNumber}/{qStart}_{qEnd} 
// example: m130802_221257_00127_c100560082550000001823094812221334_s1_p0/128361/872_4288
// See: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
static SmallContainer con_pacbio_range = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "/"          },
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"          },
                             { .dict_id = { _SAM_Q3NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// {movieName}/{holeNumber}/ccs 
// example: m64136_200621_234916/18/ccs
// see: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
static SmallContainer con_pacbio_label = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "/"          }, // monotonically-increasing integer if collated
                             { .dict_id = { _SAM_Q2NAME }                            }, // usually all-the-same
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// @<MovieName>/<ZMW_number>
static SmallContainer con_pacbio_plain = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/"          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

#define PX_pacbio { "m" }

// Oxford Nanopore MinION/GridION/PromethION:
// example: af84b0c1-6945-4323-9193-d9f6f2c38f9a
static SmallContainer con_nanopore = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 8  } },
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 12 } },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
};
#define PX_nanopore { "", "-", "-", "-", "-", PX_MATE_FIXED_0_PAD }

// example: 2a228edf-218c-46b3-b1b8-3d613b8530dc_39-13665
static SmallContainer con_nanopore_rng = {
    .repeats             = 1,
    .nitems_lo           = 8,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 8  } },
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 12 } },
                             { .dict_id = { _SAM_Q5NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q6NAME },                                      }, 
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
};
#define PX_nanopore_rng { "", "-", "-", "-", "-", "_", "" }

// example: 2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template
static SmallContainer con_nanopore_ext = {
    .repeats             = 1,
    .nitems_lo           = 7,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 8  } },
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 12 } },
                             { .dict_id = { _SAM_Q5NAME }                                       },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
};
#define PX_nanopore_ext { "", "-", "-", "-", "-", "_", "" }

// example: 22:33597495-34324994_726956_727496_0:0:0_0:0:0_2963e
static SmallContainer con_bamsurgeon = {
    .repeats             = 1,
    .nitems_lo           = 9,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q5NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q6NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q7NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};


// example: SRR11215720.1_1_length=120: items 1 and 2 are usually equal, item 3 equals seq_len
static SmallContainer con_ncbi_sra_L = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."                    }, // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q1NAME }, .separator = "_"                    }, // usually sequential number
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"                    }, // usually equal the FASTQ file number (e.g 1 or 2)
                             { .dict_id = { _SAM_Q3NAME }                                      } }// length
};

#define PX_sra_len { "", "", "", "length=" }

// example: ERR811170.1
static SmallContainer con_ncbi_sra = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters (usually SRR, ERR or DRR) - usually all-the-same context. See: https://www.biostars.org/p/381527/ and https://www.ncbi.nlm.nih.gov/books/NBK569234/#search.what_do_the_different_sra_accessi
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    },  // SRA number - sometimes all-the-same but not always - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }                                      } } // usually sequential number in FASTQ
};

// example: ERR2708427.1.1
static SmallContainer con_ncbi_sra2 = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    },  // SRA number - sometimes all-the-same but not always - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }, .separator = "."                    },  // 1st number - usually sequential number in FASTQ
                             { .dict_id = { _SAM_Q3NAME }                                      } } // 2nd number - observed as FASTQ file number (i.e. 1 or 2)
};

// example: basic.1
static SmallContainer con_genozip_opt = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// example: consensus:1331
static SmallContainer con_consensus = {
    .repeats             = 1,
    .nitems_lo           = 1,
    .items               = { { .dict_id = { _SAM_Q0NAME } } }
};

#define PX_consensus { "consensus:", "" }

// example: read_1
static SmallContainer con_str_integer = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// example: adeno-reads100.fasta.000000008
static SmallContainer con_seqan = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."          },
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 9 } },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

#define PX_seqan { "", "", "", "" }

// example: "umi64163_count1" - generated by CLC Genomics Workbench
static SmallContainer con_clc_gw = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q1NAME },                           },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

#define PX_clc_gw { "umi", "count", "" }

// example: 1
static SmallContainer con_integer = {
    .repeats             = 1,
    .nitems_lo           = 2,
    .items               = { { .dict_id = { _SAM_Q0NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// example: 30cf_chr10 (possibly from wgsim simulator, not sure)
static SmallContainer con_hex_chr = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 4 } }, 
                             { .dict_id = { _SAM_Q1NAME }                                      },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }
};

#define PX_hex_chr { "", "_" }

//--------------------------------------------------------------------------------------------------------

#define QFS_MAX_EXAMPLES 5

typedef struct QnameFlavorStruct {
    QnameFlavorId id;                     // optional; required only if referenced in the code
    char name[16]; 
    char example[QFS_MAX_EXAMPLES][256];
    SeqTech tech;                         // The sequencing technology used to generate this data
    SeqTech qname1_tech;                  // QNAME2 is available only if tech==qname1_tech     
    QType only_q;                         // this QF can appear only as QNAME or only as QNAME2
    SmallContainer *con_template;         // container template - copied to con in qname_zip_initialize
    char cut_to_canonize;                 // terminate before the last character that is this, to canonoize 
    unsigned num_seps;                    // number of printable separator and prefix characters
    int integer_items[MAX_QNAME_ITEMS+1]; // array of item_i (0-based) that are expected to be integer - no leading zeros (+1 for termiantor; terminated by -1)
    int numeric_items[MAX_QNAME_ITEMS+1]; // array of item_i (0-based) that are expected to be numeric - leading zeros ok (+1 for termiantor; terminated by -1)
    int in_local[MAX_QNAME_ITEMS+1];      // int or numeric items that should be stored in local. if not set, item is segged to dict. if also in order_item, then if sorted by qname - will created a delta snip, otherwise store as local. 
    int hex_items[MAX_QNAME_ITEMS+1];     // array of item_i (0-based) that are expected to be lower case hexademical (+1 for termiantor; terminated by -1). a hex item must also be either integer or numeric.
    bool sam_qname_sorted;                // qnames are expected to be sorted, even if BAM if sorted by coordinates. This happens when SAM QNAMES have be replaced with a sorted flavor.
    int ordered_item1, ordered_item2;     // an item that may be delta'd in non-sorted files. this should be the fast-moving item which consecutive values are close
    int range_end_item1,range_end_item2;  // an item containing end of range - delta vs preceding item in container
    int seq_len_item;                     // item containing the seq_len of this read
    int barcode_item;                     // item contains barcode. if also in_local, then data is compressed with CODEC_ACGT
    int callback_item;                    // item is segged using a callback
    unsigned fixed_len;                   // fixed length (0 if it is not fixed length)
    rom px_strs[MAX_QNAME_ITEMS+1];       // prefixes of items (+1 for terminator; terminated by empty entry) 
    mSTRl (con_snip, NUM_QTYPES, 256);    // container snips (generated in qname_zip_initialize)
    #define MAX_PREFIX_LEN 30
    STRl (con_prefix, MAX_PREFIX_LEN);    // prefix of container
    SmallContainer con;                   // container
    bool is_integer[MAX_QNAME_ITEMS], is_hex[MAX_QNAME_ITEMS], is_in_local[MAX_QNAME_ITEMS], is_numeric[MAX_QNAME_ITEMS]; // indexed according to order of items in the container (NOT by order of did_i)
    bool is_mated;                        // true means qname has a /1 or /2, qfs generated with qname_genarate_qfs_with_mate()
    int barcode_item2;                    // set if we have a barcode item with a '+' separator
} QnameFlavorStruct;

static QnameFlavorStruct qf[] = { 
/*  mate    id             name             example                                       tech          qname1_tech   only_q  con_template       canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc cb  len px_strs           test_file  */
    // QNAMEs generated by sequencers
    {},  { QF_ILLUM_7gsFQ, "Illumina-gsFQ", { "ATATA-ATGCATAG|ab|A00488:61:HMLGNDSXX:4:1101:4345:1000" },   // must be before QF_ILLUM_7
                                                                                          TECH_ILLUM,   TECH_NCBI,    QNAME1, &con_illumina_7gsFQ,0,   8,  {3,5,6,7,8,-1},     {-1},           {3,5,-1},           {-1},           0,  7,8,   -1,-1, -1, -1, -1,                       /* flavor.illumina-gsFQ.fq */},
         { QF_ILLUM_7gs,   "Illumina-gs",   { "ATATA-ATGCATAG|ab|A00488:61:HMLGNDSXX:4:1101:4345:1000|1" },   
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_7gs,  0,   9,  {3,5,6,7,8,9,-1},   {-1},           {3,5,-1},           {-1},           0,  7,8,   -1,-1, -1, -1, -1,    .is_mated=true     /* flavor.illumina-gs.sam */}, // is_mated is set so that the |1 suffix is removed during canonization
    {},  { QF_ILLUM_7,     "Illumina",      { "A00488:61:HMLGNDSXX:4:1101:4345:1000" },   TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7,    0,   6,  {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, -1, -1,                       /* flavor.illumina-7.fq */ },
    {},  { QF_ILLUM_7i,    "Illumina#bc",   { "A00488:61:HMLGNDSXX:4:1101:4345:1000#CTGGGAAG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7i,   '#', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, 7,  -1,                       /* flavor.illumina#bc.sam */ },
    {},  { QF_ILLUM_7umi,  "Illumina-umi",  { "A00488:61:HMLGNDSXX:4:1101:4345:1000;umi=ACCTTCCAA" },   
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7umi, ';', 11, {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, 7,  -1, 0,  PX_illumina_7umi  /* flavor.illumina_7umi.fq */ },
    {},  { QF_ILLUM_7_2bc, "Illumina:2bc",  { "A00488:61:HMLGNDSXX:4:1101:4345:1000:TGCTGGG+ACTTTTA" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7bc2, ':', 8,  {3,4,5,6,-1},       {-1},           {3,-1},             {-1},           0,  5,6,   -1,-1, -1, 7,  -1,                       /* flavor.illumina-7-2bc.fq */ },
    {},  { QF_ILLUM_7_bc,  "Illumina:bc",   { "SDF-02:GFH-0166::1:13435:2311:1233:GTAGCCAATCA" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7bc,  ':', 7,  {3,4,5,6,-1},       {-1},           {3,-1},             {-1},           0,  5,6,   -1,-1, -1, 7,  -1,                       /* flavor.illumina-7-bc.fq, flavor.illumina.colon.sam */ },
    {},  { QF_SINGULAR,    "Singular",      { "B05:000:FC2:4:1:272670:483" },             TECH_SINGLR,  TECH_NCBI,    QANY,   &con_singular,      0,   6,  {3,4,5,6,-1},       {1,-1},         {5,6,-1},           {-1},           0,  6,-1,  -1,-1, -1, -1, -1, 0,  PX_SINGULAR       /* flavor.singualr.fq */ },
    {},  { QF_ELEMENT,     "Element",       { "PLT-16:APP-0316:UNKNOWN_FLOWCELL:1:10102:0582:0027", "PLT-03:BBS-0174:2140948523:1:10102:0293:0058" },     
                                                                                          TECH_ELEMENT, TECH_NCBI,    QANY,   &con_element,       0,   6,  {3,-1},             {5,6,-1},       {5,6,-1},           {-1},           0,  6,-1,  -1,-1, -1, -1, -1, 0,  PX_ELEMENT        /* flavor.element.fq */ },
    {},  { QF_ELEMENT_6,   "Element-6",     { "PLT-04:APP-0289::1:10102:0338:0025" },     TECH_ELEMENT, TECH_NCBI,    QANY,   &con_element_6,     0,   6,  {2,3,-1},           {4,5,-1},       {4,5,-1},           {-1},           0,  5,-1,  -1,-1, -1, -1, -1, 0,  PX_ELEMENT_6      /* flavor.element-6.fq */ },
/*  mate    id             name             example                                       tech          qname1_tech   only_q  con_template       canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc cb  len px_strs           test_file  */
    {},  { QF_BGI_varlen,  "BGI-varlen",    { "8A_V100004684L3C001R029311637", "V300022116L2C001R0012002968", "V300046476L1C001R00110001719" },          
                                                                                          TECH_BGI,     TECH_NCBI,    QANY,   &con_bgi_varlen,    0,   3,  {4,-1},             {1,2,3,-1},     {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, 0,  PX_bgi_varlen     },
    {},  { QF_BGI_r6,      "BGI-R6",        { "8A_V100004684L3C001R029011637", "V300014293BL2C001R027005967", "V300003413L4C001R016000000" },          
                                                                                          TECH_BGI,     TECH_NCBI,    QANY,   &con_bgi_R6,        0,   3,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, 0,  PX_bgi_R          },
    {},  { QF_BGI_r7,      "BGI-R7",        { "V300017009_8AL2C001R0030001805", "V300022116L2C001R0010002968", "V300014296L2C001R0013000027", "E100001117L1C001R0030000000", "E1000536L1C002R0020000005" },         
                                                                                          TECH_BGI,     TECH_NCBI,    QANY,   &con_bgi_R7,        0,   3,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, 0,  PX_bgi_R          },
    {},  { QF_BGI_rgs8FQ,  "BGI-Rgs8FQ",    { "CGGTCT-AACCT|ab|E200003777L1C001R00100888074" },         // must be before QF_BGI_r8
                                                                                          TECH_BGI,     TECH_NCBI,    QNAME1, &con_bgi_RgsFQ8,    0,   5,  {-1},               {3,4,5,6,-1},   {4,5,6,-1},         {-1},           0,  6,5,   -1,-1, -1, -1, -1, 0,  PX_bgi_RgsFQ      },
         { QF_BGI_rgs8,    "BGI-Rgs8",      { "CGGTCT-AACCT|ab|E200003777L1C001R00100888074|2" },       
                                                                                          TECH_BGI,     TECH_NCBI,    QSAM,   &con_bgi_Rgs8,      0,   6,  {-1},               {3,4,5,6,-1},   {4,5,6,-1},         {-1},           0,  6,5,   -1,-1, -1, -1, -1, 0,  PX_bgi_Rgs, .is_mated=true }, // is_mated is set so that the |1 suffix is removed during canonization
    {},  { QF_BGI_r8,      "BGI-R8",        { "V300046476L1C001R00100001719" },           TECH_BGI,     TECH_NCBI,    QANY,   &con_bgi_R8,        0,   3,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, 0,  PX_bgi_R          },
    {},  { QF_BGI_ll7,     "BGI-LL7",       { "DP8400010271TLL1C005R0511863479" },        TECH_BGI,     TECH_NCBI,    QANY,   &con_bgi_LL7,       0,   4,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, 0,  PX_bgi_LL         },
    {},  { QF_BGI_cl,      "BGI-CL",        { "CL100025298L1C002R050_244547" },           TECH_BGI,     TECH_NCBI,    QANY,   &con_bgi_CL,        0,   6,  {4,-1},             {1,2,3,-1},     {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, 0,  PX_bgi_CL         }, 
         { QF_ULTIMA_a,    "Ultima.a",      { "004733_1-X0003-0072646116" },              TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_a,      0,   3,  {1,-1},             {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 25, PX_ULTIMA_A       /* flavor.ultima-a.sam */},
         { QF_ULTIMA_a_bc, "Ultima.a-bc",   { "004733_1-X0003-0072646116_TCGTCACTCGAAAACT" },         
                                                                                          TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_a_bc,   0,   4,  {1,-1},             {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, 4,  -1, 0,  PX_ULTIMA_A_BC    /* flavor.ultima-a-bc.fq */},
         { QF_ULTIMA_b,    "Ultima.b",      { "014214_2-UFOx1-3-9871314132" },            TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_b,      0,   4,  {1,-1},             {4,-1},         {1,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_ULTIMA_B       /* flavor.ultima-b.fq  */},
         { QF_ULTIMA_b_bc, "Ultima.b_bc",   { "014214_2-UFOx1-143-9871314132_AGTAC" },    TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_b_bc,   0,   5,  {1,-1},             {4,-1},         {1,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, 5,  -1, 0,  PX_ULTIMA_B_BC    /* flavor.ultima-b-bc.fq  */},
         { QF_ULTIMA_d,    "Ultima.d",      { "014214-UFOx1-3-9871314132" },              TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_d,      0,   3,  {-1},               {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_ULTIMA_D       /* flavor.ultima-d.fq  */},
         { QF_ULTIMA_d_bc, "Ultima.d_bc",   { "014214-UFOx1-143-9871314132_AGTAC" },      TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_d_bc,   0,   4,  {-1},               {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, 4,  -1, 0,  PX_ULTIMA_D_BC    /* flavor.ultima-d-bc.fq  */},
         { QF_ULTIMA_c,    "Ultima.c",      { "V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:" },         
                                                                                          TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_c,      0,   14, {1,4,5,6,7,8,9,10,-1},{-1},         {4,7,8,10,-1},      {-1},           0,  -1,-1, -1,-1, -1, -1, 5,  0,                    /* flavor.ultima-c.fq */},
         { QF_ULTIMA_c_bc, "Ultima.c_bc",   { "V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:_GCTGCTGACA" },         
                                                                                          TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_c_bc,   0,   15, {1,4,5,6,7,8,9,10,-1},{-1},         {4,7,8,10,-1},      {-1},           0,  -1,-1, -1,-1, -1, 14, 5,  0,                    /* flavor.ultima-c-bc.fq */},
    {},  { QF_ION_TORR_3,  "IonTorrent",    { "ZEWTM:10130:07001" },                      TECH_IONTORR, TECH_NCBI,    QANY,   &con_ion_torrent_3, 0,   2,  {-1},               {-1},           {-1},               {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 17                    },
    {},  { QF_ILLUM_5i,    "Illumina-old#", { "HWI-ST550_0201:3:1101:1626:2216#ACAGTG" }, TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_5i,   '#', 5,  {1,2,3,4,-1},       {-1},           {-1},               {-1},           0,  -1,-1, -1,-1, -1, 5,  -1,                       },
    {},  { QF_ILLUM_5,     "Illumina-old",  { "SOLEXA-1GA-1_4_FC20ENL:7:258:737:870" },   TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_5,    0,   4,  {1,2,3,4,-1},       {-1},           {1,2,3,4-1},        {-1},           0,  -1,-1, -1,-1, -1, -1, -1,                       },
    {},  { QF_ILLUM_5rng,  "Illumina-oldR", { "NOVID_3053_FC625AGAAXX:6:1:1069:11483:0,84" },   
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_5rng, ':', 6,  {1,2,3,4,5,6,-1},   {-1},           {1,2,3,4,5,6-1},    {-1},           0,  -1,-1, -1,-1, 6,  -1, -1,                       },
    {},  { QF_ROCHE_454,   "Roche-454",     { "000050_1712_0767" },                       TECH_454,     TECH_NCBI,    QANY,   &con_roche_454,     0,   2,  {-1},               {0,1,2,-1},     {-1},               {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 16, PX_roche_454      },
    {},  { QF_HELICOS,     "Helicos",       { "VHE-242383071011-15-1-0-2" },              TECH_HELICOS, TECH_NCBI,    QANY,   &con_helicos,       0,   5,  {2,3,4,5,-1},       {-1},           {-1},               {-1},           0,  -1,-1, -1,-1, -1, -1, -1,                       },
    {},  { QF_PACBIO_3,    "PacBio-3",      { "0ae26d65_70722_4787" },                    TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_3,      0,   2,  {1,2,-1},           {0,-1},         {1,2,-1},           {0,-1},         0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_pacbio_3       /* test.pacbio-blasr.bam */},
    {},  { QF_PACBIO_rng,  "PacBio-Range",  { "m130802_221257_00127_c100560082550000001823094812221334_s1_p0/128361/872_4288" },
                                                                                          TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_range,  0,   4,  {1,2,3,-1},         {-1},           {1,-1},             {-1},           0,  -1,-1,  3,-1, -1, -1, -1, 0,  PX_pacbio         }, // note: tried ord1=1 for pacbio flavors, we're worse off than in_local
    {},  { QF_PACBIO_lbl,  "PacBio-Label",  { "m64136_200621_234916/18/ccs" },            TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_label,  0,   3,  {1,-1},             {-1},           {1,-1},             {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_pacbio         /* flavor.pacbio_label.aux.fq */},
    {},  { QF_PACBIO_pln,  "PacBio-Plain",  { "m64136_200621_234916/18" },                TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_plain,  0,   2,  {1,-1},             {-1},           {1,-1},             {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_pacbio         },
    {},  { QF_NANOPORE,    "Nanopore",      { "af84b0c1-6945-4323-9193-d9f6f2c38f9a" },   TECH_ONP,     TECH_NCBI,    QANY,   &con_nanopore,      0,   4,  {-1},               {0,1,2,3,4-1},  {0,1,2,3,4,-1},     {0,1,2,3,4,-1}, 0,  -1,-1, -1,-1, -1, -1, -1, 36, PX_nanopore       /* flavor.nanopore.aux.fq, flavor.nanopore.length.fq */},
    {},  { QF_NANOPORE_rng,"Nanopore-rng",  { "2a228edf-218c-46b3-b1b8-3d613b8530dc_39-13665" },
                                                                                          TECH_ONP,     TECH_NCBI,    QANY,   &con_nanopore_rng,  0,   6,  {5,6,-1},           {0,1,2,3,4,-1}, {0,1,2,3,4,5,6,-1}, {0,1,2,3,4,-1}, 0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_nanopore_rng   /* flavor.nanopore-rng.fq */}, // 14.0.31
/*  mate    id             name             example                                       tech          qname1_tech   only_q  con_template       canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc cb  len px_strs           test_file  */
    {},  { QF_NANOPORE_ext,"Nanopore-ext",  { "2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template" },
                                                                                          TECH_ONP,     TECH_NCBI,    QANY,   &con_nanopore_ext,  0,   5,  {-1},               {0,1,2,3,4,-1}, {0,1,2,3,4,-1},     {0,1,2,3,4,-1}, 0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_nanopore_ext   },
    {},  { QF_BAMSURGEON,  "BamSurgeon",    { "22:33597495-34324994_726956_727496_0:0:0_0:0:0_2963e" },   
                                                                                          TECH_UNKNOWN, NO_QNAME2,    QNAME1, &con_bamsurgeon,    0,   7,  {1,2,3,4,7,-1},     {-1},           {1,3,7,-1},         {7,-1},         0,  1,3,   2,4,   -1, -1, -1,                       },
    // NCBI QNAMEs - no mate, as mate is part of QNAME2 in this case
         { QF_SRA_L,       "NCBI_SRA_L",    { "SRR11215720.1_1_length=120" },             TECH_NCBI,    NO_QNAME2,    Q1or3,  &con_ncbi_sra_L,    0,   10, {1,2,-1},           {-1},           {2,-1},             {-1},           0,  2,-1,  -1,-1,  3, -1, -1, 0,  PX_sra_len        },
         { QF_SRA2,        "NCBI-SRA2",     { "ERR2708427.1.1" },                         TECH_NCBI,    NO_QNAME2,    Q1or3,  &con_ncbi_sra2,     0,   2,  {2,3,-1},           {-1},           {2,3,-1},           {-1},           0,  2,-1,  -1,-1, -1, -1, -1,     .is_mated=true    /* flavor.nanopore.length.fq */},
         { QF_SRA,         "NCBI-SRA",      { "SRR001666.1" },                            TECH_NCBI,    NO_QNAME2,    Q1or3,  &con_ncbi_sra,      0,   1,  {2,-1},             {-1},           {2,-1},             {-1},           0,  2,-1,  -1,-1, -1, -1, -1,                       /* flavor.sra.sam, flavor.sra.fq */},

    // QNAME2 - Illumina and Singular - no mate, as mate is part of QNAME in this case
         { QF_ILLUM_2bc,   "Illumina-2bc",  { "2:N:0:CTGAAGCT+ATAGAGGC" },                TECH_ILLUM,   TECH_ILLUM,   QNAME2, &con_illumina_2bc,  0,   4,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, 0,  PX_illumina_2bc   /* flavor.illumina_7_fq.fq */ }, 
         { QF_ILLUM_0bc,   "Illumina-0bc",  { "1:N:0:0" },                                TECH_ILLUM,   TECH_ILLUM,   QNAME2, &con_illumina_0bc,  0,   3,  {0,2,3,-1},         {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_illumina_0bc   /* special.solexa-R1.fq */  }, // v14.0.0
         { QF_ILLUM_1bc,   "Illumina-1bc",  { "2:N:0:GATATTAC" },                         TECH_ILLUM,   TECH_ILLUM,   QNAME2, &con_illumina_1bc,  0,   3,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, 0,  PX_illumina_1bc   /* flavor.illumina_1bc.fq, flavor.illumina-bc-fq.fq, flavor.illumina_0bc.fq */}, 
         { QF_SINGULR_1bc, "Singular-1bc",  { "2:N:0:NNNNNN" },                           TECH_SINGLR,  TECH_SINGLR,  QNAME2, &con_illumina_1bc,  0,   3,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, 0,  PX_illumina_1bc   /* flavor.singualr.fq */}, 
         { QF_ELEMENT_2bc, "Element-2bc",   { "1:N:0:CTACAGTG+AGTCGCTT" },                TECH_ELEMENT, TECH_ELEMENT, QNAME2, &con_illumina_2bc,  0,   4,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, 0,  PX_illumina_2bc   /* flavor.illumina_7_fq.fq */ }, 
         { QF_ELEMENT_0bc, "Element-0bc",   { "1:N:0:1" },                                TECH_ELEMENT, TECH_ELEMENT, QNAME2, &con_illumina_0bc,  0,   3,  {0,2,3,-1},         {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_illumina_0bc   /* special.solexa-R1.fq */  }, 
         { QF_ELEMENT_1bc, "Element-1bc",   { "1:N:0:CTACAGTG" },                         TECH_ELEMENT, TECH_ELEMENT, QNAME2, &con_illumina_1bc,  0,   3,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, 0,  PX_illumina_1bc   /* flavor.illumina_1bc.fq, flavor.illumina-bc-fq.fq, flavor.illumina_0bc.fq */}, 
         { QF_BGI_Rgs_2bc, "BGI-Rgs-2bc",   { "2:N:0:CTGAAGCT+ATAGAGGC" },                TECH_BGI,     TECH_BGI,     QNAME2, &con_illumina_2bc,  0,   4,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, 0,  PX_illumina_2bc   /* flavor.illumina_7_fq.fq */ }, 

    // observed as QNAME2 in NCBI (possibly with mate) and in SAM/BAM
    {},  { QF_ILLUM_X_2bc, "Illumina_X_2bc",{ "A00180:28:HC3F5DRXX:2:2110:27453:21981_1:N:0:ATTACTCGATCT+GGCTCTGA" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    Q2orSAM,&con_illumina_X_2bc,'_', 11, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1,                       /* flavor.illumina_X_2bc.fq */},
    {},  { QF_ILLUM_X_1bc, "Illumina_X_1bc",{ "ST-E00314:354:H7J2YCCXY:1:1101:7080:1450_1:N:0:NAGGCG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    Q2orSAM,&con_illumina_X_1bc,'_', 10, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1,                       /* flavor.illumina_X_1bc.sam */},
    {},  { QF_ILLUM_X_0bc, "Illumina_X_0bc",{ "A00488:61:HMLGNDSXX:4:1101:4345:1000_2:N:0" },
                                                                                          TECH_ILLUM,   TECH_NCBI,    Q2orSAM,&con_illumina_X_0bc,'_', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,-1},           {-1},           0,  5,6,   -1,-1, -1, -1, -1,                       /* flavor.illumina_X_0bc.sam */ }, // v14.0.17  
    // only possible in SAM/BAM
    {},  { QF_ILLUM_S_2bc, "Illumina_S_2bc",{ "A00180:28:HC3F5DRXX:2:2110:27453:21981 1:N:0:ATTACTCGATCT+GGCTCTGA" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_S_2bc,' ', 11, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1,                       /* flavor.illumina_S_2bc.sam */},
    {},  { QF_ILLUM_S_1bc, "Illumina_S_1bc",{ "ST-E00314:354:H7J2YCCXY:1:1101:7080:1450 1:N:0:NAGGCG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_S_1bc,' ', 10, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1,                       },
    {},  { QF_ILLUM_S_0bc, "Illumina_S_0bc",{ "A00488:61:HMLGNDSXX:4:1101:4345:1000 2:N:0" },
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_S_0bc,' ', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,-1},           {-1},           0,  5,6,   -1,-1, -1, -1, -1,                       }, // v14.0.17  
    // QNAMEs generated by various software tools
    {},  { QF_SEQAN,       "seqan",         { "adeno-reads100.fasta.000000008" },         TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_seqan,         0,   2,  {-1},               {2, -1},        {-1},               {-1},           0,  2,-1,  -1,-1, -1, -1, -1, 0,  PX_seqan         },
    {},  { QF_CLC_GW,      "CLC-GW",        { "umi64163_count1" },                        TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_clc_gw,        0,   9,  {0,1,-1},           {-1},           {0,1,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_clc_gw        },
    {},  { QF_HEX_CHR,     "hex_chr",       { "30cf_chr10" }, /* wgsim simulator? */      TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_hex_chr,       0,   1,  {-1},               {0,-1},         {0,-1},             {0,-1},         0,  -1,-1, -1,-1, -1, -1, -1, 0,  PX_hex_chr       }, // added v14
    {},  { QF_INTEGER,     "Integer",       { "123" },                                    TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_integer,       0,   0,  {0,-1},             {-1},           {0,-1},             {-1},           0,  0,-1,  -1,-1, -1, -1, -1,                       /* flavor.ncbi_srr_fq.fq */}, 
    {},  { QF_STR_INT,     "Str_Integer",   { "read_1" },   /* eg CLC */                  TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_str_integer,   0,   1,  {1,-1},             {-1},           {1,-1},             {-1},           0,  1,-1,  -1,-1, -1, -1, -1,                       },
         { QF_CONSENSUS,   "consensus",     { "consensus:23" },                           TECH_UNKNOWN, NO_QNAME2,    QSAM,   &con_consensus,     0,   10, {0,-1},             {-1},           {0,-1},             {-1},           1,  0,-1,  -1,-1, -1, -1, -1, 0,  PX_consensus      /* flavor.concensus.sam  */ },
    {},  { QF_GENOZIP_OPT, "Genozip-opt",   { "basic.1" },  /* must be last */            TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_genozip_opt,   0,   1,  {1,-1},             {-1},           {1,-1},             {-1},           1,  1,-1,  -1,-1, -1, -1, -1,                       /* flavor.Genozip-opt.fq */ },
};
#define NUM_QFs ARRAY_LEN(qf) // note: different than NUM_FLAVORS bc each flavor might have two QFs: one for mated

static QnameSegCallback qf_callbacks[NUM_FLAVORS] = { [QF_ULTIMA_c] = ultima_c_Q5NAME_cb, [QF_ULTIMA_c_bc] = ultima_c_Q5NAME_cb };

