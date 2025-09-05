// ------------------------------------------------------------------
//   qname_flavors.h
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
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

#define PX_MATE_FIXED_0_PAD (char[]){ CI0_SKIP, 0 } // add to prefix for mate item IFF preceeding item has CI0_FIXED_0_PAD
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

// Example: A00488:61:HMLGNDSXX:4:1101:4345:1000:CAGACGCGCACATACTTTTCTCACG
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

static SmallContainer con_illumina_7rbc = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = {':','r'}    },
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
static SmallContainer con_qname2_2bc = {
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
static SmallContainer con_qname2_0bc = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q3NAME }                                      } }  // barcode
};

#define PX_illumina_0bc { "", "", ":", "" } 

// Example: 2:N:0:CTGAAGCT
static SmallContainer con_qname2_1bc = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q3NAME }                                      } }  // barcode
};

#define PX_illumina_1bc { "", "", ":", "" } 

// Example: "DNBSEQT7:001:E100012314:1:002:0020000077:0020000077"
// Example: "MGI2000:001:V300053419:2:003:00100001039:00100001039"
#define CON_MGI_NEW(n)                                                                                                  \
static SmallContainer con_mgi_new##n = {                                                                                \
    .repeats             = 1,                                                                                           \
    .nitems_lo           = 10,                                                                                          \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                     },  /* Sequencer    */  \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 3  } },  /* ?            */  \
                             { .dict_id = { _SAM_Q2NAME }, .separator = ":"                     },  /* Flow cell    */  \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 1  } },  /* Lane         */  \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 3  } },  /* Column       */  \
                             { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 3  } },  /* Row          */  \
                             { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, n  } },  /* Tile         */  \
                             { .dict_id = { _SAM_Q7NAME }, .separator = { CI0_FIXED_0_PAD, 3  } },  /* Copy of Row  */  \
                             { .dict_id = { _SAM_Q8NAME }, .separator = { CI0_FIXED_0_PAD, n  } },  /* Copy of Tile */  \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } } /* Mate         */  \
};
CON_MGI_NEW(6);
CON_MGI_NEW(7);
CON_MGI_NEW(8);

#define PX_mgi_new { "", "", ":", "", ":", ":", "", ":", "", PX_MATE_FIXED_0_PAD }

static bool val_mgi_new (STRps(item)) 
{
    return str_isprefix_(STRi(item,0), _S("DNBSEQ")) || str_isprefix_(STRi(item,0), _S("MGI"));
};

// SAM + FASTQ: callback defined in qname_flavor for QF_MGI_NEW 
sSTRl(copy_q5name, 16);
sSTRl(copy_q6name, 16);
static void seg_qname_mgi_new_cb (VBlockP vb, ContextP ctx, STRp(copy))
{
    STRlast (src, ctx->did_i-2);

    if (str_issame (src, copy)) {
        if (ctx->did_i == SAM_Q7NAME)
            seg_by_ctx (vb, STRa(copy_q5name), ctx, copy_len);
        else if (ctx->did_i == SAM_Q8NAME)
            seg_by_ctx (vb, STRa(copy_q6name), ctx, copy_len);
        else
            ABORT ("unexpected ctx=%s", ctx->tag_name);
    }

    else
        seg_numeric_or_not (vb, ctx, STRa(copy), copy_len);
}

// example: SOME:2:PREFIX:L01:R001C012:0000:8199
static SmallContainer con_mgi_sap8 = {
    .repeats   = 1,
    .nitems_lo = 7,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":L"                   },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 2 } },
                   { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                   { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                   { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } 
};

#define PX_mgi_sap8 { "", "", ":R", "C", ":", ":", PX_MATE_FIXED_0_PAD }

//--------------------------------------------------------------------------------------------------------------
// MGI_E format: 27 characters: E, FlowCellSerialNumber[various], L, Lane[1], C, Column[3], R, Row[3] Tile[6-8]
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
#define CON_MGI_R(n)                    \
static SmallContainer con_mgi_R##n = {  \
    .repeats             = 1,           \
    .nitems_lo           = 6,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_MGI_R(6);
CON_MGI_R(7);
CON_MGI_R(8);

#define PX_mgi_R { "", "", "C", "R", "", PX_MATE_FIXED_0_PAD }

// example: die1_A100004684C001R029011637
#define CON_MGI_die(n)                   \
static SmallContainer con_mgi_die##n = {  \
    .repeats             = 1,            \
    .nitems_lo           = 6,            \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Die       */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "C"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_MGI_die(6);

#define PX_mgi_die { "die", "_", "", "R", "", PX_MATE_FIXED_0_PAD }

#define CON_MGI_Rgs(n) /* SAM/BAM only*/    \
static SmallContainer con_mgi_Rgs##n = {    \
    .repeats             = 1,               \
    .nitems_lo           = 8,               \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "|"                    }, /* the two parts of the barcode are correletated and hence segged together */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "|"                    },                 \
                             { .dict_id = { _SAM_Q2NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_Q7NAME }                                      } }/* number of reads merged together during the consensus analysis step */ \
}
CON_MGI_Rgs(8);

#define PX_mgi_Rgs { "", "", "", "", "C", "R", "", "|" }

// Example: CGGTCT-AACCT|ab|E200003777L1C001R00100888074
#define CON_MGI_RgsFQ(n) /* FASTQ only */ \
static SmallContainer con_mgi_RgsFQ##n = {  \
    .repeats             = 1,           \
    .nitems_lo           = 8,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "|"                    }, /* the two parts of the barcode are correletated and hence segged together */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "|"                    },                 \
                             { .dict_id = { _SAM_Q2NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_MGI_RgsFQ(8);

#define PX_mgi_RgsFQ { "", "", "", "", "C", "R", "", PX_MATE_FIXED_0_PAD }

// variant of MGI flavor where Q4NAME is a variable-length integer rather than a fixed-length zero-padded numeric
static SmallContainer con_mgi_varlen = {  \
    .repeats             = 1,           \
    .nitems_lo           = 6,           \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "L"                    }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME },                                     }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
};

#define PX_mgi_varlen { "", "", "C", "R" }

//--------------------------------------------------------------------------------------------------------------
// Same as CON_MGI_R, but the first separator is two L
//  DP8400010271TLL1C005R0511863479

#define CON_MGI_LL(n) \
static SmallContainer con_mgi_LL##n = {  \
    .repeats             = 1,            \
    .nitems_lo           = 6,            \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { 'L', 'L'           } }, /* Flow cell */ \
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, /* Lane      */ \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Column    */ \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, /* Row       */ \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, n } }, /* Tile      */ \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }/* Mate      */ \
}
CON_MGI_LL(7); // only encountered with 7 so far

#define PX_mgi_LL PX_mgi_R

//--------------------------------------------------------------------------------------------------------------------
// MGI_CL format: CL, FlowCellSerialNumber[9], L, Lane[1], C, Column[3], R, Row[3], _, variable-length number
//  CL100025298L1C002R050_244547 reported by https://en.wikipedia.org/wiki/File:BGI_seq_platform_read_name_description.png
// @CL100072652L2C001R001_12/1 for FASTQ reported by https://www.biostars.org/p/395715/
//--------------------------------------------------------------------------------------------------------------------
static SmallContainer con_mgi_CL = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "L" },                     /* Flow cell */
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },  /* Lane      */
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  /* Column    */
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  /* Row       */
                             { .dict_id = { _SAM_Q4NAME },                                     },  /* Tile      */
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } /* Mate      */
};

#define PX_mgi_CL { "CL", "", "C", "R", "_" } // "CL" is a prefix, and the 2nd L is seprator 

// ULTIMA_1 format
// Example: 004733_1-X0003-1345904491
static SmallContainer con_ultima_a = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } } 
};

#define PX_ULTIMA_A { "", "_", "", "-", PX_MATE_FIXED_0_PAD }  

// Example: 004733_1-X0003-0072646116_TCGTCACTCGAAAACT
static SmallContainer con_ultima_a_bc = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_Q4NAME }                                       },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
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
}; // no I_AM_MATE bc final separator not supported yet. TO DO: fix this

// Example: V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:_GCTGCTGACA
static SmallContainer con_ultima_c_bc = {
    .repeats   = 1,
    .nitems_lo = 16,
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
                   { .dict_id = { _SAM_QENAME }                    },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE         } }
};

// QF_ULTIMA_Z? format
// Example: 012345678_1-Z0123-0123456789
#define CON_ULTIMA_Z(n) \
static SmallContainer con_ultima_Z##n = {                                                           \
    .repeats             = 1,                                                                       \
    .nitems_lo           = 5,                                                                       \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, n  } },  \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     },  \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } },  \
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },  \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } } \
};
CON_ULTIMA_Z(9)

#define PX_ULTIMA_Z { "", "_", "", "-", PX_MATE_FIXED_0_PAD }  

// QF_ULTIMA_b? format
// Example: 014214_2-UGAv1-143-9871314132
#define CON_ULTIMA_B(n) \
static SmallContainer con_ultima_b##n = {                                                           \
    .repeats             = 1,                                                                       \
    .nitems_lo           = 6,                                                                       \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, n  } },  \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     },  \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } },  \
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-" },                      \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },  \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } } \
};
CON_ULTIMA_B(6)
CON_ULTIMA_B(9)

#define PX_ULTIMA_B { "", "_", "", "-", "", PX_MATE_FIXED_0_PAD }  

// QF_ULTIMA_b?_bc format
// Example: 014214_2-UGAv1-143-9871314132_AGTAC
#define CON_ULTIMA_B_BC(n) \
static SmallContainer con_ultima_b##n##_bc = {                                                      \
    .repeats             = 1,                                                                       \
    .nitems_lo           = 7,                                                                       \
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, n  } },  \
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"                     },  \
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5  } },  \
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-" },                      \
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },  \
                             { .dict_id = { _SAM_Q5NAME }                                       },  \
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } } \
};
CON_ULTIMA_B_BC(6)
CON_ULTIMA_B_BC(9)

#define PX_ULTIMA_B_BC { "", "_", "", "-", "", "_" }  

// QF_ULTIMA_d format
// Example: 014214-UGAv1-143-9871314132
static SmallContainer con_ultima_d = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } }, 
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
};

#define PX_ULTIMA_D { "", "-", "-", "", PX_MATE_FIXED_0_PAD }  

// QF_ULTIMA_d_bc format
// Example: 014214-UGAv1-143-9871314132_AGTAC
static SmallContainer con_ultima_d_bc = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 6  } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-" }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_Q4NAME }                                       }, 
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
};

#define PX_ULTIMA_D_BC { "", "-", "-", "", "_" }  

// example: 1
static SmallContainer con_ultima_n = {
    .repeats             = 1,
    .nitems_lo           = 2,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 10 } },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                            } }
};

#define PX_ULTIMA_N { "", PX_MATE_FIXED_0_PAD }  

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

// Example: SDF-02:GFH-0166::1:13435:2311:1233:GTAGCCAATCA
static SmallContainer con_element_bc = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_Q6NAME }, .separator = { CI0_FIXED_0_PAD, 4 } },
                   { .dict_id = { _SAM_Q7NAME },                                     },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } 
};

#define PX_ELEMENT_BC { "", "", "", "", "", "", ":", ":", "" }  

//--------------------------------------------------------------------------------------------------------------
// ION_TORRENT_3 format: 17 characters: RunID[5] : Row[5] : Column[5]
// Example: ZEWTM:10130:07001
// See: https://www.biostars.org/p/330644/
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_ion_torrent_3 = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
                   { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } }
};

#define PX_ion_torrent_3 { "", ":", ":", PX_MATE_FIXED_0_PAD } 

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

// Example: HWI-ST156_288:4:1:10000:112636:0
static SmallContainer con_illumina_6 = {
    .repeats   = 1,
    .nitems_lo = 7,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME },                           },
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
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/"              }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "/"              },
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"              },
                             { .dict_id = { _SAM_Q3NAME }                                },
                             { .dict_id = { _SAM_Q4NAME }, .separator[0] = CI0_INVISIBLE } } // used to calculate seq_len from range
}; // note: mate not yet support with CI0_INVISIBLE (TO DO)

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

static bool val_pacbio (STRps(item)) 
{
    // verify that first item is stats with something like : 64136_200621_234916 (possibly more components)
    str_split (items[0], item_lens[0], 0, '_', subitem, false);

    return n_subitems >= 3 && str_is_numeric (STRi(subitem, 1)) && str_is_numeric (STRi(subitem, 2));
};

static SmallContainer con_onso = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"                    },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q2NAME }, .separator = "-"                    },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q5NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
                   { .dict_id = { _SAM_Q6NAME }, .separator = ":"                    },
                   { .dict_id = { _SAM_Q7NAME },                                     },
                   { .dict_id = { _SAM_QmNAME }, I_AM_MATE                           } } 
};

#define PX_onso { "", "", "", "", "", "", ":", "" }

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


// example: SRR11215720.1_1_length=120: items 2 and 3 are usually equal, item 4 equals seq_len
static SmallContainer con_ncbi_sra_L = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, // three uppercase letters (usually SRR, ERR or DRR) - usually all-the-same context. See: https://www.biostars.org/p/381527/ and https://www.ncbi.nlm.nih.gov/books/NBK569234/#search.what_do_the_different_sra_accessi
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    }, // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"                    }, // usually sequential number
                             { .dict_id = { _SAM_Q3NAME }, .separator = "_"                    }, // usually equal the FASTQ file number (e.g 1 or 2)
                             { .dict_id = { _SAM_Q4NAME }                                      } }// length
};

#define PX_sra_len { "", "", "", "", "length=" }

// example: ERR811170.1
static SmallContainer con_ncbi_sra = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three uppercase letters (usually SRR, ERR or DRR) - usually all-the-same context. See: https://www.biostars.org/p/381527/ and https://www.ncbi.nlm.nih.gov/books/NBK569234/#search.what_do_the_different_sra_accessi
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    },  // SRA number - sometimes all-the-same but not always - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }                                      } } // usually sequential number in FASTQ
};

// used SRA-sra, example: "SRR001666.sra.1"
#define PX_sra_sra { "", "", "sra." }

// example: ERR811170.1_mylabel
static SmallContainer con_ncbi_sra_label = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three uppercase letters (usually SRR, ERR or DRR) - usually all-the-same context. See: https://www.biostars.org/p/381527/ and https://www.ncbi.nlm.nih.gov/books/NBK569234/#search.what_do_the_different_sra_accessi
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    },  // SRA number - sometimes all-the-same but not always - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"                    },  // usually sequential number in FASTQ
                             { .dict_id = { _SAM_Q3NAME }                                      } } // an arbitrary label
};

// example: ERR2708427.1.1
static SmallContainer con_ncbi_sra2 = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three uppercase letters - usually all-the-same
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    },  // SRA number - sometimes all-the-same but not always - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }, .separator = "."                    },  // 1st number - usually sequential number in FASTQ
                             { .dict_id = { _SAM_Q3NAME }                                      } } // 2nd number - observed as FASTQ file number (i.e. 1 or 2)
};

static bool val_sra (STRps(item)) 
{
    return item_lens[0]==3 && IS_CLETTER (items[0][0]) && IS_CLETTER (items[0][1]) && IS_CLETTER (items[0][2]) &&
           str_is_numeric (STRi(item, 1));
};

// example: basic.1
static SmallContainer con_genozip_opt = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

// example: mapped.ILLUMINA.bwa:1
static SmallContainer con_generated = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, I_AM_MATE                 } }
};

static SmallContainer con_prfx_and_int = {
    .repeats             = 1,
    .nitems_lo           = 1,
    .items               = { { .dict_id = { _SAM_Q0NAME } } }
};

#define PX_consensus { "consensus:", "" } // example: consensus:1331
#define PX_cons      { "cons"      , "" } // example: cons6532
#define PX_Sint      { "S"         , "" } // example: S2348134

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

#define PX_seqan { "", "", "", PX_MATE_FIXED_0_PAD }

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

static bool no_validate (STRps(item)) { return true; };

//--------------------------------------------------------------------------------------------------------

#define QFS_MAX_EXAMPLES 2

typedef struct QnameFlavorStruct {
    QnameFlavorId id;                     // optional; required only if referenced in the code
    char name[16]; 
    char example[QFS_MAX_EXAMPLES][80];
    SeqTech tech;                         // The sequencing technology used to generate this data
    SeqTech fq_qname1_tech;               // FASTQ only: this flavor is accepted for QNAME2 only if QNAME1.tech is this value 
    QType only_q;                         // this QF can appear only as QNAME or only as QNAME2
    SmallContainer *con_template;         // container template - copied to con in qname_zip_initialize
    bool (*validate_flavor)(STRps (item));// extra validation (beyond the qf table) for this flavor
    char cut_to_canonize;                 // terminate before the last character that is this, to canonoize 
    unsigned num_seps;                    // number of printable separator and prefix characters
    int integer_items[MAX_QNAME_ITEMS+1]; // array of item_i (0-based) that are expected to be integer - no leading zeros (+1 for termiantor; terminated by -1)
    int numeric_items[MAX_QNAME_ITEMS+1]; // array of item_i (0-based) that are expected to be numeric - leading zeros ok (+1 for termiantor; terminated by -1)
    int in_local[MAX_QNAME_ITEMS+1];      // int or numeric items that should be stored in local. if not set, item is segged to dict. if also in order_item, then if sorted by qname - will created a delta snip, otherwise store as local. 
    int hex_items[MAX_QNAME_ITEMS+1];     // array of item_i (0-based) that are expected to be lower case hexademical (+1 for termiantor; terminated by -1). a hex item must also be either integer or numeric.
    bool sam_qname_sorted;                // qnames are expected to be sorted, even if BAM if sorted by coordinates. This happens when SAM QNAMES have be replaced with a sorted flavor.
    int ordered_item1, ordered_item2;     // an item that may be delta'd in non-POS-sorted files. this should be the fast-moving item which consecutive values are close
    int range_end_item1,range_end_item2;  // an item containing end of range - delta vs preceding item in container
    int seq_len_item;                     // item containing the seq_len of this read
    int barcode_item;                     // item contains barcode. if also in_local, then data is compressed with CODEC_ACGT
    int callback_item, callback_item2;    // item is segged using a callback
    unsigned fixed_len;                   // fixed length (0 if it is not fixed length)
    rom px_strs[MAX_QNAME_ITEMS+1];       // prefixes of items (+1 for terminator; terminated by empty entry) 
    mSTRl (con_snip, NUM_QTYPES, 280);    // container snips (generated in qname_zip_initialize)
    #define MAX_PREFIX_LEN 64
    STRl (con_prefix, MAX_PREFIX_LEN);    // prefix of container
    SmallContainer con;                   // container
    bool is_integer[MAX_QNAME_ITEMS], is_hex[MAX_QNAME_ITEMS], is_in_local[MAX_QNAME_ITEMS], is_numeric[MAX_QNAME_ITEMS]; // indexed according to order of items in the container (NOT by order of did_i)
    bool is_mated;                        // true means qname has a /1 or /2 - and that mates (defined by is_first/is_last SAM flags) have opposite /1 vs /2. This field is generated with qname_genarate_qfs_with_mate()
    int barcode_item2;                    // set if we have a barcode item with a '+' separator
} QnameFlavorStruct;

static QnameFlavorStruct qf[] = { 
/*  mate    id             name             example                                       tech          fq_qname1_tech only_q con_template         validate_func  canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc RG  cb1 cb2  len px_strs           */
    // QNAMEs generated by sequencers
    {},  { QF_ILLUM_7gsFQ, "Illumina-gsFQ", { "ATATA-ATGCATAG|ab|A00488:61:HMLGNDSXX:4:1101:4345:1000" },   // must be before QF_ILLUM_7
                                                                                          TECH_ILLUM,   TECH_NCBI,    QNAME1, &con_illumina_7gsFQ, no_validate,    0,   8,  {3,5,6,7,8,-1},     {-1},           {3,5,-1},           {-1},           0,  7,8,   -1,-1, -1, -1, -1, -1,                       },
         { QF_ILLUM_7gs,   "Illumina-gs",   { "ATATA-ATGCATAG|ab|A00488:61:HMLGNDSXX:4:1101:4345:1000|1" },   
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_7gs,   no_validate,    '|', 9,  {3,5,6,7,8,9,-1},   {-1},           {3,5,-1},           {-1},           0,  7,8,   -1,-1, -1, -1, -1, -1,                       },
    {},  { QF_ILLUM_7,     "Illumina",      { "A00488:61:HMLGNDSXX:4:1101:4345:1000" },   TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7,     no_validate,    0,   6,  {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, -1, -1, -1,                       },
    {},  { QF_ILLUM_7i,    "Illumina#bc",   { "A00488:61:HMLGNDSXX:4:1101:4345:1000#CTGGGAAG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7i,    no_validate,    '#', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, 7,  -1, -1,                       },
    {},  { QF_ILLUM_7umi,  "Illumina-umi",  { "A00488:61:HMLGNDSXX:4:1101:4345:1000;umi=ACCTTCCAA" },   
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7umi,  no_validate,    ';', 11, {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, 7,  -1, -1, 0,  PX_illumina_7umi  },
    {},  { QF_ILLUM_7_2bc, "Illumina-2bc",  { "A00488:61:HMLGNDSXX:4:1101:4345:1000:GNTGTCA+GCGTTGT", }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7bc2,  no_validate,    ':', 8,  {3,4,5,6,-1},       {-1},           {3,5,6,-1},         {-1},           0,  5,6,   -1,-1, -1, 7,  -1, -1,                       },
    {},  { QF_ILLUM_7_rbc, "Illumina-rbc",  { "A00488:61:HMLGNDSXX:4:1101:4345:1000:rTGTATGTCCC" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7rbc,  no_validate,    ':', 8,  {3,4,5,6,-1},       {-1},           {3,5,6,-1},         {-1},           0,  5,6,   -1,-1, -1, 7,  -1, -1,                       },
    {},  { QF_ILLUM_7_bc,  "Illumina-bc",   { "A00488:61:HMLGNDSXX:4:1101:4345:1000:CAGACGCGCACATACTTTTCTCACG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_7bc,   no_validate,    ':', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,5,6,-1},       {-1},           0,  5,6,   -1,-1, -1, 7,  -1, -1,                       },
    {},  { QF_SINGULAR,    "Singular",      { "B05:000:FC2:4:1:272670:483" },             TECH_SINGLR,  TECH_NCBI,    QANY,   &con_singular,       no_validate,    0,   6,  {3,4,5,6,-1},       {1,-1},         {5,6,-1},           {-1},           0,  6,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_SINGULAR       },
    {},  { QF_ELEMENT,     "Element",       { "PLT-03:BBS-0174:2140948523:1:10102:0293:0058", "PLT-16:APP-0316:UNKNOWN_FLOWCELL:1:10102:0582:0027" },     
                                                                                          TECH_ELEMENT, TECH_NCBI,    QANY,   &con_element,        no_validate,    0,   6,  {3,4,-1},           {5,6,-1},       {5,6,-1},           {-1},           0,  5,6,   -1,-1, -1, -1, -1, -1, 0,  PX_ELEMENT        },
    {},  { QF_ELEMENT_bc,  "Element-bc",    { "SDF-02:GFH-0166::1:13435:2311:1233:GTAGCCAATCA", "SDF-02:GFH-0166:2140948523:1:13435:2311:1233:GTAGCCAATCA"  }, 
                                                                                          TECH_ELEMENT, TECH_NCBI,    QANY,   &con_element_bc,     no_validate,    ':', 7,  {3,4,-1},           {5,6,-1},       {5,6,-1},           {-1},           0,  5,6,   -1,-1, -1, 7,  -1, -1, 0,  PX_ELEMENT_BC     },
/*  mate    id             name             example                                       tech          fq_qname1_tech  only_q con_template        validate_func  canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc cb1 cb2 len px_strs           */
    {},  { QF_MGI_NEW6,    "MGI-NEW6",      { "DNBSEQT7:001:E100012314:1:002:002000077:002000077" }, // synthetic example, not yet observed          
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_new6,       val_mgi_new,    0,   6,  {-1},               {1,3,4,5,6,7,-1},{4,5,6,-1},        {-1},           0,  5,6,   -1,-1, -1, -1, 7,  8,  0,  PX_mgi_new        }, // 15.0.51 
    {},  { QF_MGI_NEW7,    "MGI-NEW7",      { "DNBSEQT7:001:E100012314:1:002:0020000077:0020000077" },          
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_new7,       val_mgi_new,    0,   6,  {-1},               {1,3,4,5,6,7,-1},{4,5,6,-1},        {-1},           0,  5,6,   -1,-1, -1, -1, 7,  8,  0,  PX_mgi_new        }, // 15.0.51 
    {},  { QF_MGI_NEW8,    "MGI-NEW8",      { "MGI2000:001:V300053419:2:003:00100001039:00100001039" },          
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_new8,       val_mgi_new,    0,   6,  {-1},               {1,3,4,5,6,7,-1},{4,5,6,-1},        {-1},           0,  5,6,   -1,-1, -1, -1, 7,  8,  0,  PX_mgi_new        }, // 15.0.51 
    {},  { QF_MGI_SAP8,    "MGI-SAP8",      { "SOME:2:PREFIX:L01:R001C012:0000:8199" },          
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_sap8,       no_validate,    0,   7,  {-1},               {1,2,3,4,5,-1}, {2,3,4,5,-1},       {-1},           0,  4,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_mgi_sap8       }, // 15.0.70
    {},  { QF_MGI_varlen,  "MGI-varlen",    { "8A_V100004684L3C001R029311637", "V300046476L1C001R00110001719" },          
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_varlen,     no_validate,    0,   3,  {4,-1},             {1,2,3,-1},     {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_varlen     },
    {},  { QF_MGI_r6,      "MGI-R6",        { "8A_V100004684L3C001R029011637", "V300003413L4C001R016000000" },          
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_R6,         no_validate,    0,   3,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_R          },
    {},  { QF_MGI_die6,    "MGI-die6",      { "die1_A100004684C001R029011637" },          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_die6,       no_validate,    0,   6,  {-1},               {0,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_die        }, // 15.0.67
    {},  { QF_MGI_r7,      "MGI-R7",        { "V300017009_8AL2C001R0030001805", "E100001117L1C001R0030000000" },         
                                                                                          TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_R7,         no_validate,    0,   3,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_R          },
    {},  { QF_MGI_rgs8FQ,  "MGI-Rgs8FQ",    { "CGGTCT-AACCT|ab|E200003777L1C001R00100888074" },         // must be before QF_MGI_r8
                                                                                          TECH_MGI,     TECH_NCBI,    QNAME1, &con_mgi_RgsFQ8,     no_validate,    0,   5,  {-1},               {3,4,5,6,-1},   {4,5,6,-1},         {-1},           0,  6,5,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_RgsFQ      },
         { QF_MGI_rgs8,    "MGI-Rgs8",      { "CGGTCT-AACCT|ab|E200003777L1C001R00100888074|2" },       
                                                                                          TECH_MGI,     TECH_NCBI,    QSAM,   &con_mgi_Rgs8,       no_validate,    '|', 6,  {-1},               {3,4,5,6,-1},   {4,5,6,-1},         {-1},           0,  6,5,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_Rgs,       }, 
    {},  { QF_MGI_r8,      "MGI-R8",        { "V300046476L1C001R00100001719" },           TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_R8,         no_validate,    0,   3,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_R          },
    {},  { QF_MGI_ll7,     "MGI-LL7",       { "DP8400010271TLL1C005R0511863479" },        TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_LL7,        no_validate,    0,   4,  {-1},               {1,2,3,4,-1},   {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_LL         },
    {},  { QF_MGI_cl,      "MGI-CL",        { "CL100025298L1C002R050_244547" },           TECH_MGI,     TECH_NCBI,    QANY,   &con_mgi_CL,         no_validate,    0,   6,  {4,-1},             {1,2,3,-1},     {2,3,4,-1},         {-1},           0,  4,3,   -1,-1, -1, -1, -1, -1, 0,  PX_mgi_CL         }, 
    {},  { QF_ULTIMA_a,    "Ultima-a",      { "012345_1-X0003-0072646116" },              TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_a,       no_validate,    0,   3,  {1,-1},             {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 25, PX_ULTIMA_A       },
    {},  { QF_ULTIMA_a_bc, "Ultima-a_bc",   { "012345_1-X0003-0072646116_TCGTCACTCGAAAACT" },         
                                                                                          TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_a_bc,    no_validate,    '_', 4,  {1,-1},             {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, 4,  -1, -1, 0,  PX_ULTIMA_A_BC    },
    {},  { QF_ULTIMA_b6,   "Ultima-b6",     { "012345_2-UGAv1-3-9871314132" },            TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_b6,      no_validate,    0,   4,  {1,-1},             {4,-1},         {1,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_ULTIMA_B       },
    {},  { QF_ULTIMA_b6_bc,"Ultima-b6_bc",  { "012345_2-UGAv1-143-9871314132_AGTAC" },    TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_b6_bc,   no_validate,    '_', 5,  {1,-1},             {4,-1},         {1,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, 5,  -1, -1, 0,  PX_ULTIMA_B_BC    },
    {},  { QF_ULTIMA_b9,   "Ultima-b9",     { "012345678_2-UGAv3-3-9871314132" },         TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_b9,      no_validate,    0,   4,  {1,-1},             {4,-1},         {1,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_ULTIMA_B       },
    {},  { QF_ULTIMA_Z9,   "Ultima-Z9",     { "012345678_1-Z0123-0123456789" },           TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_Z9,      no_validate,    0,   3,  {1,-1},             {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_ULTIMA_Z       },
    {},  { QF_ULTIMA_b9_bc,"Ultima-b9_bc",  { "012345678_2-UGAv3-143-9871314132_AGTAC" }, TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_b9_bc,   no_validate,    '_', 5,  {1,-1},             {4,-1},         {1,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, 5,  -1, -1, 0,  PX_ULTIMA_B_BC    },
    {},  { QF_ULTIMA_d,    "Ultima-d",      { "012345-UGAv1-3-9871314132" },              TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_d,       no_validate,    0,   3,  {-1},               {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_ULTIMA_D       },
    {},  { QF_ULTIMA_d_bc, "Ultima-d_bc",   { "012345-UGAv1-143-9871314132_AGTAC" },      TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_d_bc,    no_validate,    '_', 4,  {-1},               {3,-1},         {1,3,-1},           {-1},           0,  -1,-1, -1,-1, -1, 4,  -1, -1, 0,  PX_ULTIMA_D_BC    },
         { QF_ULTIMA_c,    "Ultima-c",      { "V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:" },         
                                                                                          TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_c,       no_validate,    ':', 14, {1,4,5,6,7,8,9,10,-1},{-1},         {4,7,8,10,-1},      {-1},           0,  -1,-1, -1,-1, -1, -1, 5,  -1, 0,                    },
    {},  { QF_ULTIMA_c_bc, "Ultima-c_bc",   { "V222:23526:::1:1:7:9831:222:1:443:N:0.99:Z0199:_GCTGCTGACA" },         
                                                                                          TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_c_bc,    no_validate,    ':', 15, {1,4,5,6,7,8,9,10,-1},{-1},         {4,7,8,10,-1},      {-1},           0,  -1,-1, -1,-1, -1, 14, 5,  -1, 0,                    },
    {},  { QF_ULTIMA_n,    "Ultima-n",      { "0871314132" },                             TECH_ULTIMA,  TECH_NCBI,    QANY,   &con_ultima_n,       no_validate,    0,   0,  {-1},               {0,-1},         {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_ULTIMA_N       },

    {},  { QF_ONSO,        "Onso"    ,      { "PSQ003:86:FB0031380-BCC:1:01001:4705:166" },   
                                                                                          TECH_ONSO,    TECH_NCBI,    QANY,   &con_onso,           no_validate,    0,   7,  {1,4,6,7,-1},       {5,-1},         {6,-1},             {-1},           0,  7,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_onso           },
    {},  { QF_ION_TORR_3,  "IonTorrent",    { "ZEWTM:00130:07001" },                      TECH_IONTORR, TECH_NCBI,    QANY,   &con_ion_torrent_3,  no_validate,    0,   2,  {-1},               {1,2,-1},       {1,2,-1},           {-1},           0,  1,2,   -1,-1, -1, -1, -1, -1, 17, PX_ion_torrent_3  },
    {},  { QF_ILLUM_5i,    "Illum-old#bc",  { "HWI-ST550_0201:3:1101:1626:2216#ACAGTG" }, TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_5i,    no_validate,    '#', 5,  {1,2,3,4,-1},       {-1},           {3,4,-1},           {-1},           0,  -1,-1, -1,-1, -1, 5,  -1, -1,                       },
    {},  { QF_ILLUM_5,     "Illum-old",     { "SOLEXA-1GA-1_4_FC20ENL:7:258:737:870" },   TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_5,     no_validate,    0,   4,  {1,2,3,4,-1},       {-1},           {1,2,3,4,-1},       {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1,                       },
    {},  { QF_ILLUM_5rng,  "Illum-oldR",    { "NOVID_3053_FC625AGAAXX:6:1:1069:11483:0,84" },   
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_5rng,  no_validate,    ':', 6,  {1,2,3,4,5,6,-1},   {-1},           {1,2,3,4,5,6,-1},   {-1},           0,  -1,-1, -1,-1, 6,  -1, -1, -1,                       },
    {},  { QF_ILLUM_6,     "Illum-old6",    { "HWI-ST156_288:4:1:10000:110537:0" },       TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_6,     no_validate,    0,   5,  {1,2,3,4,5,-1},     {-1},           {1,2,3,4,-1},       {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1,                       },
    {},  { QF_ROCHE_454,   "Roche-454",     { "000050_1712_0767" },                       TECH_LS454,   TECH_NCBI,    QANY,   &con_roche_454,      no_validate,    0,   2,  {-1},               {0,1,2,-1},     {-1},               {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 16, PX_roche_454      },
    {},  { QF_HELICOS,     "Helicos",       { "VHE-242383071011-15-1-0-2" },              TECH_HELICOS, TECH_NCBI,    QANY,   &con_helicos,        no_validate,    0,   5,  {2,3,4,5,-1},       {-1},           {-1},               {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1,                       },
    {},  { QF_PACBIO_3,    "PacBio-3",      { "0ae26d65_70722_4787" },                    TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_3,       no_validate,    0,   2,  {1,2,-1},           {0,-1},         {1,2,-1},           {0,-1},         0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_pacbio_3       },
         { QF_PACBIO_rng,  "PacBio-Range",  { "m130802_221257_00127_c100560082550000001823094812221334_s1_p0/128361/872_4288" },
                                                                                          TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_range,   val_pacbio,     0,   4,  {1,2,3,-1},         {-1},           {1,-1},             {-1},           0,  1,-1,  3,-1,  4,  -1,  4, -1, 0,  PX_pacbio         }, // note: tried ord1=1 for pacbio flavors, we're worse off than in_local
    {},  { QF_PACBIO_lbl,  "PacBio-Label",  { "m64136_200621_234916/18/ccs" },            TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_label,   val_pacbio,     0,   3,  {1,-1},             {-1},           {1,-1},             {-1},           0,  1,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_pacbio         },
    {},  { QF_PACBIO_pln,  "PacBio-Plain",  { "m64136_200621_234916/18" },                TECH_PACBIO,  TECH_NCBI,    QANY,   &con_pacbio_plain,   val_pacbio,     0, 2,  {1,-1},             {-1},           {1,-1},             {-1},           0,  1,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_pacbio         },
    {},  { QF_NANOPORE,    "Nanopore",      { "af84b0c1-6945-4323-9193-d9f6f2c38f9a" },   TECH_NANOPORE,TECH_NCBI,    QANY,   &con_nanopore,       no_validate,    0,   4,  {-1},               {0,1,2,3,4,-1}, {0,1,2,3,4,-1},     {0,1,2,3,4,-1}, 0,  -1,-1, -1,-1, -1, -1, -1, -1, 36, PX_nanopore       },
    {},  { QF_NANOPORE_rng,"Nanopore-rng",  { "2a228edf-218c-46b3-b1b8-3d613b8530dc_39-13665" },
                                                                                          TECH_NANOPORE,TECH_NCBI,    QANY,   &con_nanopore_rng,   no_validate,    '_', 6,  {5,6,-1},           {0,1,2,3,4,-1}, {0,1,2,3,4,5,6,-1}, {0,1,2,3,4,-1}, 0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_nanopore_rng   }, // 14.0.31
/*  mate    id             name             example                                       tech          fq_qname1_tech only_q con_template         validate_func  canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc cb1 cb2 len px_strs           */
    {},  { QF_NANOPORE_ext,"Nanopore-ext",  { "2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template" },
                                                                                          TECH_NANOPORE,TECH_NCBI,    QANY,   &con_nanopore_ext,   no_validate,    0,   5,  {-1},               {0,1,2,3,4,-1}, {0,1,2,3,4,-1},     {0,1,2,3,4,-1}, 0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_nanopore_ext   },
    {},  { QF_BAMSURGEON,  "BamSurgeon",    { "22:33597495-34324994_726956_727496_0:0:0_0:0:0_2963e" },   
                                                                                          TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_bamsurgeon,     no_validate,    0,   7,  {1,2,3,4,7,-1},     {-1},           {1,3,7,-1},         {7,-1},         0,  1,3,   2,4,   -1, -1, -1, -1,                       },
    // NCBI QNAMEs - no mate, as mate is part of QNAME2 in this case
         { QF_SRA_L,       "NCBI_SRA_L",    { "SRR11215720.1_1_length=120" },             TECH_NCBI,    TECH_NONE,    Q1or3,  &con_ncbi_sra_L,     val_sra,        0,   10, {2,3,4,-1},         {-1},           {2,3,4,-1},         {-1},           0,  3,-1,  -1,-1,  4, -1, -1, -1, 0,  PX_sra_len        },
         { QF_SRA2,        "NCBI-SRA2",     { "ERR2708427.1.1" },                         TECH_NCBI,    TECH_NONE,    Q1or3,  &con_ncbi_sra2,      val_sra,        0,   2,  {2,3,-1},           {-1},           {2,3,-1},           {-1},           0,  2,-1,  -1,-1, -1, -1, -1, -1,     .is_mated=true    },
         { QF_SRA,         "NCBI-SRA",      { "SRR001666.1" },                            TECH_NCBI,    TECH_NONE,    Q1or3,  &con_ncbi_sra,       val_sra,        0,   1,  {2,-1},             {-1},           {2,-1},             {-1},           0,  2,-1,  -1,-1, -1, -1, -1, -1,                       },
         { QF_SRA_label,   "NCBI-SRA-label",{ "SRR001666.1_mylabel" },                    TECH_NCBI,    TECH_NONE,    Q1or3,  &con_ncbi_sra_label, val_sra,        0,   2,  {2,-1},             {-1},           {2,-1},             {-1},           0,  2,-1,  -1,-1, -1, -1, -1, -1,                       }, // 15.0.74
         { QF_SRA_sra,     "NCBI-SRA-sra",  { "SRR001666.sra.1" },                        TECH_NCBI,    TECH_NONE,    Q1or3,  &con_ncbi_sra,       val_sra,        0,   5,  {2,-1},             {-1},           {2,-1},             {-1},           0,  2,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_sra_sra        }, // 15.0.67

    // QNAME2 - Illumina, Singular, Element... - no mate, as mate is part of QNAME in this case
         { QF_ILLUM_2bc,   "Illum-2bc",     { "2:N:0:CTGAAGCT+ATAGAGGC" },                TECH_NONE,    TECH_ANY,     QNAME2, &con_qname2_2bc,     no_validate,    0,   4,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, -1, 0,  PX_illumina_2bc   }, 
         { QF_ILLUM_0bc,   "Illum-0bc",     { "1:N:0:0" },                                TECH_NONE,    TECH_ANY,     QNAME2, &con_qname2_0bc,     no_validate,    0,   3,  {0,2,3,-1},         {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_illumina_0bc   }, // v14.0.0
         { QF_ILLUM_1bc,   "Illum-1bc",     { "2:N:0:GATATTAC", "2:N:0:NNNNNN" },         TECH_NONE,    TECH_ANY,     QNAME2, &con_qname2_1bc,     no_validate,    0,   3,  {0,2,-1},           {-1},           {0,-1},             {-1},           0,  -1,-1, -1,-1, -1, 3,  -1, -1, 0,  PX_illumina_1bc   }, 

    // observed as QNAME2 in NCBI (possibly with mate), SAM/BAM and when generated FASTQ from SAM/BAM
    {},  { QF_ILLUM_X_2bc, "Illumina_X_2bc",{ "A00180:28:HC3F5DRXX:2:2110:27453:21981_1:N:0:ATTACTCGATCT+GGCTCTGA" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_X_2bc, no_validate,    '_', 11, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1, -1,                       },
    {},  { QF_ILLUM_X_1bc, "Illumina_X_1bc",{ "ST-E00314:354:H7J2YCCXY:1:1101:7080:1450_1:N:0:NAGGCG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_X_1bc, no_validate,    '_', 10, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1, -1,                       },
    {},  { QF_ILLUM_X_0bc, "Illumina_X_0bc",{ "A00488:61:HMLGNDSXX:4:1101:4345:1000_2:N:0" },
                                                                                          TECH_ILLUM,   TECH_NCBI,    QANY,   &con_illumina_X_0bc, no_validate,    '_', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,-1},           {-1},           0,  5,6,   -1,-1, -1, -1, -1, -1,                       }, // v14.0.17  
    // only possible in SAM/BAM
    {},  { QF_ILLUM_S_2bc, "Illumina_S_2bc",{ "A00180:28:H50C3F5DRXX:2:2110:27453:21981 1:N:0:ATTACTCGATCT+GGCTCTGA" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_S_2bc, no_validate,    ' ', 11, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1, -1,                       },
    {},  { QF_ILLUM_S_1bc, "Illumina_S_1bc",{ "ST-E00314:354:H7J2YCCXY:1:1101:7080:1450 1:N:0:NAGGCG" }, 
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_S_1bc, no_validate,    ' ', 10, {1,3,4,5,6,7,9,-1}, {-1},           {1,3,7,8,9,-1},     {-1},           0,  5,6,   -1,-1, -1, 10, -1, -1,                       },
    {},  { QF_ILLUM_S_0bc, "Illumina_S_0bc",{ "A00488:61:HMLGNDSXX:4:1101:4345:1000 2:N:0" },
                                                                                          TECH_ILLUM,   TECH_NCBI,    QSAM,   &con_illumina_S_0bc, no_validate,    ' ', 7,  {1,3,4,5,6,-1},     {-1},           {1,3,-1},           {-1},           0,  5,6,   -1,-1, -1, -1, -1, -1,                       }, // v14.0.17  
/*  mate    id             name             example                                       tech          fq_qname1_tech only_q con_template         validate_func  canon #sp integer_items       numeric_items   in-local            hex_items       srt ord1,2 rng    sqln bc cb1 cb2 len px_strs           */
    // QNAMEs generated by various software tools
    {},  { QF_SEQAN,       "seqan",         { "adeno-reads100.fasta.000000008" },         TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_seqan,          no_validate,    0,   2,  {-1},               {2, -1},        {-1},               {-1},           0,  2,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_seqan          },
    {},  { QF_CLC_GW,      "CLC-GW",        { "umi64163_count1" },                        TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_clc_gw,         no_validate,    0,   9,  {0,1,-1},           {-1},           {0,1,-1},           {-1},           0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_clc_gw         },
    {},  { QF_HEX_CHR,     "hex_chr",       { "30cf_chr10" }, /* wgsim simulator? */      TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_hex_chr,        no_validate,    0,   1,  {-1},               {0,-1},         {0,-1},             {0,-1},         0,  -1,-1, -1,-1, -1, -1, -1, -1, 0,  PX_hex_chr        }, // added v14
    {},  { QF_INTEGER,     "Integer",       { "123" },                                    TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_integer,        no_validate,    0,   0,  {0,-1},             {-1},           {0,-1},             {-1},           0,  0,-1,  -1,-1, -1, -1, -1, -1,                       }, 
    {},  { QF_STR_INT,     "Str_Integer",   { "read_1" },   /* eg CLC */                  TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_str_integer,    no_validate,    0,   1,  {1,-1},             {-1},           {1,-1},             {-1},           0,  1,-1,  -1,-1, -1, -1, -1, -1,                       },
         { QF_CONSENSUS,   "consensus",     { "consensus:23" },                           TECH_CONS,    TECH_NCBI,    QANY,   &con_prfx_and_int,   no_validate,    0,   10, {0,-1},             {-1},           {0,-1},             {-1},           1,  0,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_consensus      },
         { QF_CONS,        "cons",          { "cons113" },                                TECH_CONS,    TECH_NCBI,    QANY,   &con_prfx_and_int,   no_validate,    0,   4,  {0,-1},             {-1},           {0,-1},             {-1},           1,  0,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_cons           },
         { QF_Sint,        "Sint",          { "S522414" },                                TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_prfx_and_int,   no_validate,    0,   1,  {0,-1},             {-1},           {0,-1},             {-1},           0,  0,-1,  -1,-1, -1, -1, -1, -1, 0,  PX_Sint           },
    {},  { QF_GENERATED,   "Generated",     { "mapped.ILLUMINA.bwa:1" },                  TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_generated,      no_validate,    0,   1,  {1,-1},             {-1},           {1,-1},             {-1},           1,  1,-1,  -1,-1, -1, -1, -1, -1,                       },
    {},  { QF_GENOZIP_OPT, "Genozip-opt",   { "basic.1" },  /* must be last */            TECH_UNKNOWN, TECH_NCBI,    QANY,   &con_genozip_opt,    no_validate,    0,   1,  {1,-1},             {-1},           {1,-1},             {-1},           1,  1,-1,  -1,-1, -1, -1, -1, -1,                       },
};
#define NUM_QFs ARRAY_LEN(qf) // note: different than NUM_FLAVORS bc each flavor might have two QFs: one for mated

static QnameSegCallback qf_callbacks[NUM_FLAVORS] = { 
    [QF_ULTIMA_c]    = ultima_c_Q5NAME_cb, 
    [QF_ULTIMA_c_bc] = ultima_c_Q5NAME_cb, 
    [QF_PACBIO_rng]  = seg_qname_rng2seq_len_cb,
    [QF_MGI_NEW6]    = seg_qname_mgi_new_cb,
    [QF_MGI_NEW7]    = seg_qname_mgi_new_cb,
    [QF_MGI_NEW8]    = seg_qname_mgi_new_cb,
};

