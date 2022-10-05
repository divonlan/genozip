// ------------------------------------------------------------------
//   qname_flavors.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "container.h"
#include "dict_id_gen.h"
#include "segconf.h"

//-----------------------------------------------------------------------------------------
// Illumina-7 format: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  
// Example: HWI-D00360:5:H814YADXX:1:1106:10370:52569
// Example FASTQ: A00488:61:HMLGNDSXX:4:1101:1940:1000 2:N:0:CTGAAGCT+ATAGAGGC
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
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } } 
};

static SmallContainer con_illumina_7_fq = {
    .repeats   = 1,
    .nitems_lo = 9,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q5NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q6NAME }, .separator = " "          },
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } },
                   { .dict_id = { _SAM_Q7NAME },                           } } // entire part after space is segged together 
};

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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } }/* Mate      */ \
}
CON_BGI_R(6);
CON_BGI_R(7);
CON_BGI_R(8);

#define PX_bgi_R { "", "", "C", "R", "", (char[]){CI0_SKIP} }

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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } }/* Mate      */ \
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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } } /* Mate      */
};

#define PX_bgi_CL { "CL", "", "C", "R", "_" } // "CL" is a prefix, and the 2nd L is seprator 

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
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } }
};
// static SmallContainer con_ion_torrent_3 = {
//     .repeats   = 1,
//     .nitems_lo = 4,
//     .items     = { { .dict_id = { _SAM_Q0NAME },                                     },  
//                    { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
//                    { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 5 } },
//                    { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } }
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
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
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
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// SOLEXA6_0104:3:1:1852:13550/1 1:N:0:0 
// See: test/special.solexa-R1.fq
static SmallContainer con_illumina_5_fq = {
    .repeats   = 1,
    .nitems_lo = 7,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q3NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q4NAME }, .separator = " "          },
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } },
                   { .dict_id = { _SAM_Q5NAME }                          } }  // entire part after space is segged together
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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } }
};

#define PX_roche_454 { "", "_", "_", (char[]){CI0_SKIP} }

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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

static SmallContainer con_pacbio_3 = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "_"          },
                             { .dict_id = { _SAM_Q2NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// {movieName}/{holeNumber}/ccs 
// example: m64136_200621_234916/18/ccs
// see: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
static SmallContainer con_pacbio_label = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "/"          },
                             { .dict_id = { _SAM_Q2NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// @<MovieName>/<ZMW_number>
static SmallContainer con_pacbio_plain = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "/"          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

#define PX_pacbio { "m" }

//TO DO: moving nanopore to use 0-padded hex fields - doesn't work yet
// example: af84b0c1-6945-4323-9193-d9f6f2c38f9a
static SmallContainer con_nanopore = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 8  } },
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 4  } }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 12 } },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP            } } }
};
#define PX_nanopore { "", "-", "-", "-", "-", (char[]){CI0_SKIP} }

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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP            } } }
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
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};


// example: SRR11215720.1_1_length=120: items 1 and 2 are usually equal, item 3 equals seq_len
static SmallContainer con_ncbi_sra_len = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."                    }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "_"                    }, // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }, .separator = "_"                    },
                             { .dict_id = { _SAM_Q3NAME }                                      },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP }           } }
};

#define PX_sra_len { "", "", "", "length=" }

// example: ERR811170.1
static SmallContainer con_ncbi_sra = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, // three letters (usually SRR, ERR or DRR) - usually all-the-same context. See: https://www.biostars.org/p/381527/ and https://www.ncbi.nlm.nih.gov/books/NBK569234/#search.what_do_the_different_sra_accessi
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    }, // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }                                      },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP }           } }
};

// example: ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template length=52
static SmallContainer con_ncbi_sra_fq = {
    .repeats             = 1,
    .nitems_lo           = 7,
    .items               = { { .dict_id = { _FASTQ_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_Q1NAME }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_Q2NAME }, .separator = " "                    },  // sequential number 
                             { .dict_id = { _FASTQ_QNAME2 }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_QmNAME }, .separator = { CI0_SKIP }           },
                             { .dict_id = { _FASTQ_Q3NAME }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_Q4NAME },                                     } } // length value
};

// example: ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template BOWDEN04_20151016_MN15199_FAA67113_BOWDEN04_MdC_MARC_Phase2a_4833_1_ch19_file1_strand length=52
static SmallContainer con_ncbi_sraP_fq = {
    .repeats             = 1,
    .nitems_lo           = 8,
    .items               = { { .dict_id = { _FASTQ_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_Q1NAME }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_Q2NAME }, .separator = " "                    },  // sequential number 
                             { .dict_id = { _FASTQ_QNAME2 }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_QmNAME }, .separator = { CI0_SKIP }           },
                             { .dict_id = { _FASTQ_Q3NAME }, .separator = " "                    },  // additional info 
                             { .dict_id = { _FASTQ_Q4NAME }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_Q5NAME },                                     } } // length value
};

// example: ERR2708427.1.1
static SmallContainer con_ncbi_sra2 = {
    .repeats             = 1,
    .nitems_lo           = 5,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }, .separator = "."                    },  // first number - sometimes, but not always, all-the-same - seg as integer in local
                             { .dict_id = { _SAM_Q3NAME }                                      },  // 2nd number 
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP }           } }
};

// example: ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 length=1128
static SmallContainer con_ncbi_sra2_fq = {
    .repeats             = 1,
    .nitems_lo           = 8,
    .items               = { { .dict_id = { _FASTQ_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_Q1NAME }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_Q2NAME }, .separator = "."                    },  // first number - sometimes, but not always, all-the-same - seg as integer in local
                             { .dict_id = { _FASTQ_Q3NAME }, .separator = " "                    },  // 2nd number - usually sequential
                             { .dict_id = { _FASTQ_QNAME2 }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_QmNAME }, .separator = { CI0_SKIP }           },
                             { .dict_id = { _FASTQ_Q4NAME }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_Q5NAME },                                     } } // length value
};

// example: ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 someextradata length=1128
static SmallContainer con_ncbi_sra2P_fq = {
    .repeats             = 1,
    .nitems_lo           = 9,
    .items               = { { .dict_id = { _FASTQ_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_Q1NAME }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_Q2NAME }, .separator = "."                    },  // first number - sometimes, but not always, all-the-same - seg as integer in local
                             { .dict_id = { _FASTQ_Q3NAME }, .separator = " "                    },  // 2nd number - usually sequential
                             { .dict_id = { _FASTQ_QNAME2 }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_QmNAME }, .separator = { CI0_SKIP }           },
                             { .dict_id = { _FASTQ_Q4NAME }, .separator = " "                    },  // additional info 
                             { .dict_id = { _FASTQ_Q5NAME }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_Q6NAME },                                     } } // length value
};

// example: basic.1
static SmallContainer con_genozip_opt = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// example: read_1
static SmallContainer con_str_integer = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// example: adeno-reads100.fasta.000000008
static SmallContainer con_seqan = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "."          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."          },
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 9 } },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

#define PX_seqan { "", "", "", "" }

// example: "umi64163_count1" - generated by CLC Genomics Workbench
static SmallContainer con_clc_gw = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q1NAME },                           },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

#define PX_clc_gw { "umi", "count", "" }

// example: 1
static SmallContainer con_integer = {
    .repeats             = 1,
    .nitems_lo           = 2,
    .items               = { { .dict_id = { _SAM_Q0NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// example: 30cf_chr10 (possibly from wgsim simulator, not sure)
static SmallContainer con_hex_chr = {
    .repeats             = 1,
    .nitems_lo           = 3,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q1NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

//--------------------------------------------
// FASTQ Line3 flavors

// example: ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template length=52
static SmallContainer con_ncbi_sra_l3 = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _FASTQ_T0HIRD }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_T1HIRD }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_T2HIRD }, .separator = " "                    },  // sequential number 
//                             { .dict_id = { _FASTQ_COPY_Q }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T3HIRD }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T4HIRD }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_T5HIRD },                                     } } // length value
};

// example: ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template BOWDEN04_20151016_MN15199_FAA67113_BOWDEN04_MdC_MARC_Phase2a_4833_1_ch19_file1_strand length=52
static SmallContainer con_ncbi_sraP_l3 = {
    .repeats             = 1,
    .nitems_lo           = 7,
    .items               = { { .dict_id = { _FASTQ_T0HIRD }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_T1HIRD }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_T2HIRD }, .separator = " "                    },  // sequential number 
//                             { .dict_id = { _FASTQ_COPY_Q }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T3HIRD }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T4HIRD }, .separator = " "                    },  // additional info 
                             { .dict_id = { _FASTQ_T5HIRD }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_T6HIRD },                                     } } // length value
};

// example: ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 length=1128
static SmallContainer con_ncbi_sra2_l3 = {
    .repeats             = 1,
    .nitems_lo           = 7,
    .items               = { { .dict_id = { _FASTQ_T0HIRD }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_T1HIRD }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_T2HIRD }, .separator = "."                    },  // first number - sometimes, but not always, all-the-same - seg as integer in local
                             { .dict_id = { _FASTQ_T3HIRD }, .separator = " "                    },  // 2nd number - usually sequential
//                             { .dict_id = { _FASTQ_COPY_Q }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T4HIRD }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T5HIRD }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_T6HIRD },                                     } } // length value
};

// example: ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 someextradata length=1128
static SmallContainer con_ncbi_sra2P_l3 = {
    .repeats             = 1,
    .nitems_lo           = 8,
    .items               = { { .dict_id = { _FASTQ_T0HIRD }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_T1HIRD }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_T2HIRD }, .separator = "."                    },  // first number - sometimes, but not always, all-the-same - seg as integer in local
                             { .dict_id = { _FASTQ_T3HIRD }, .separator = " "                    },  // 2nd number - usually sequential
//                             { .dict_id = { _FASTQ_COPY_Q }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T4HIRD }, .separator = " "                    },  // original qname 
                             { .dict_id = { _FASTQ_T5HIRD }, .separator = " "                    },  // additional info 
                             { .dict_id = { _FASTQ_T6HIRD }, .separator = "="                    },  // "length" (can't be a prefix, because no fixed_i)
                             { .dict_id = { _FASTQ_T7HIRD },                                     } } // length value
};


//--------------------------------------------------------------------------------------------------------

#define MAX_QNAME_ITEMS 11 // mate + 10 others (matching Q?NAME defined in sam.h, fastq.h, kraken.h) 
#define QFS_MAX_EXAMPLES 5

typedef struct QnameFlavorStruct {
    char name[16]; 
    char example[QFS_MAX_EXAMPLES][256];
    SeqTech tech;                         // The sequencing technology used to generate this data
    int fq_only;                          // this QF is only for FASTQ, not SAM/BAM/KRAKEN. FQ: 1 if /1 goes after the first string, 2 after the second space-separated string
    SmallContainer *con_template;         // container template - copied to con in qname_zip_initialize
    unsigned num_seps;                    // number of printable separator and prefix characters
    int integer_items[MAX_QNAME_ITEMS+1]; // array of item_i (0-based) that are expected to be integer - no leading zeros (+1 for termiantor; terminated by -1)
    int numeric_items[MAX_QNAME_ITEMS+1]; // array of item_i (0-based) that are expected to be numeric - leading zeros ok (+1 for termiantor; terminated by -1)
    int in_local[MAX_QNAME_ITEMS+1];      // int or numeric items that should be stored in local. if not set, item is segged to dict.
    int hex_items[MAX_QNAME_ITEMS+1];     // array of item_i (0-based) that are expected to be lower case hexademical (+1 for termiantor; terminated by -1)
    int ordered_item1, ordered_item2;     // an item that may be delta'd in non-sorted files. this should be the fast-moving item which consecutive values are close
    int range_end_item1,range_end_item2;  // an item containing end of range - delta vs preceding item in container
    int seq_len_item;                     // item containing the seq_len of this read
    unsigned fixed_len;                   // fixed length (0 if it is not fixed length)
    rom px_strs[MAX_QNAME_ITEMS+1];       // prefixes of items (+1 for terminator; terminated by empty entry) 
    int qname2;                           // embedded qname2 item (generated in qname_zip_initialize)
    STRl (con_snip,  200);                // container snips (generated in qname_zip_initialize)
    STRl (con_snip2, 200);                // same as con_snip, just with QNAME2 dict_ids
    #define MAX_PREFIX_LEN 30
    STRl (con_prefix, MAX_PREFIX_LEN);    // prefix of container
    SmallContainer con;                   // container
    bool is_int[MAX_QNAME_ITEMS], is_hex[MAX_QNAME_ITEMS], is_in_local[MAX_QNAME_ITEMS], is_numeric[MAX_QNAME_ITEMS]; // indexed according to order of items in the container (NOT by order of did_i)
} QnameFlavorStruct;

static QnameFlavorStruct qf[] = { 
/*  mate   name             example                                       tech     fq_only   con_template     #sp  integer_items   numeric_items   in-local        hex_items       ord1,2 rng    sqln len px_strs      */
    {},  { "Illumina-fastq",{ "A00488:61:HMLGNDSXX:4:1101:4345:1000 2:N:0:CTGAAGCT+ATAGAGGC" },
                                                                          TECH_ILLUM_7, 1, &con_illumina_7_fq, 7,  {1,3,4,5,6,-1}, {-1},           {1,3,-1},       {-1},           5,6,   -1,-1, -1,                   }, // mate added v14.0.0  
    {},  { "Illumina",      { "A00488:61:HMLGNDSXX:4:1101:4345:1000" },   TECH_ILLUM_7, 0, &con_illumina_7,    6,  {1,3,4,5,6,-1}, {-1},           {1,3,-1},       {-1},           5,6,   -1,-1, -1,                   },
    {},  { "BGI-R6",        { "8A_V100004684L3C001R029011637", "V300014293BL2C001R027005967", "V300003413L4C001R016000000" },          
                                                                          TECH_BGI,     0, &con_bgi_R6,        3,  {-1},           {1,2,3,4,-1},   {-1},           {-1},           1,-1,  -1,-1, -1,  0,  PX_bgi_R     },
    {},  { "BGI-R7",        { "V300017009_8AL2C001R0030001805", "V300022116L2C001R0010002968", "V300014296L2C001R0010000027", "E100001117L1C001R0030000000", "E1000536L1C002R0020000005" },         
                                                                          TECH_BGI,     0, &con_bgi_R7,        3,  {-1},           {1,2,3,4,-1},   {-1},           {-1},           1,-1,  -1,-1, -1,  0,  PX_bgi_R     },
    {},  { "BGI-R8",        { "V300046476L1C001R00100001719" },           TECH_BGI,     0, &con_bgi_R8,        3,  {-1},           {1,2,3,4,-1},   {-1},           {-1},           1,-1,  -1,-1, -1,  0,  PX_bgi_R     },
    {},  { "BGI-LL7",       { "DP8400010271TLL1C005R0511863479" },        TECH_BGI,     0, &con_bgi_LL7,       4,  {-1},           {1,2,3,4,-1},   {-1},           {-1},           1,-1,  -1,-1, -1,  0,  PX_bgi_LL    },
    {},  { "BGI-CL",        { "CL100025298L1C002R050_244547" },           TECH_BGI,     0, &con_bgi_CL,        4,  {4,-1},         {1,2,3,-1},     {-1},           {-1},           4,-1,  -1,-1, -1,  0,  PX_bgi_CL    }, 
    // {},  { "IonTorrent",    { "ZEWTM:10130:07001" },                      TECH_IONTORR, 0, &con_ion_torrent_3, 2,  {-1},           {1,2,-1},       {-1},           {-1},           -1,-1, -1,-1, -1,  17, PX_ion_torrent_3 },
    {},  { "IonTorrent",    { "ZEWTM:10130:07001" },                      TECH_IONTORR, 0, &con_ion_torrent_3, 2,  {-1},           {-1},           {-1},           {-1},           -1,-1, -1,-1, -1,  17               },
    {},  { "Illumina-o-fq", { "SOLEXA6_0104:3:1:1852:13550 1:N:0:0" },    TECH_ILLUM_5, 1, &con_illumina_5_fq, 5,  {1,2,3,4,-1},   {-1},           {-1},           {-1},           -1,-1, -1,-1, -1,                   }, // v14.0.0
    {},  { "Illumina-old#", { "HWI-ST550_0201:3:1101:1626:2216#ACAGTG" }, TECH_ILLUM_5, 0, &con_illumina_5i,   5,  {1,2,3,4,-1},   {-1},           {-1},           {-1},           -1,-1, -1,-1, -1,                   },
    {},  { "Illumina-old",  { "SOLEXA-1GA-1_4_FC20ENL:7:258:737:870" },   TECH_ILLUM_5, 0, &con_illumina_5,    4,  {1,2,3,4,-1},   {-1},           {1,2,3,4-1},    {-1},           -1,-1, -1,-1, -1,                   },
    {},  { "Roche-454",     { "000050_1712_0767" },                       TECH_454,     0, &con_roche_454,     2,  {-1},           {0,1,2,-1},     {-1},           {-1},           -1,-1, -1,-1, -1,  16, PX_roche_454 },
    {},  { "Helicos",       { "VHE-242383071011-15-1-0-2" },              TECH_HELICOS, 0, &con_helicos,       5,  {2,3,4,5,-1},   {-1},           {-1},           {-1},           -1,-1, -1,-1, -1,                   },
    {},  { "PacBio-3",      { "56cdb76f_70722_4787" },                    TECH_PACBIO,  0, &con_pacbio_3,      2,  {1,2,-1},       {-1},           {1,2,-1},       {0,-1},         -1,-1, -1,-1, -1,                   },
    {},  { "PacBio-Range",  { "m130802_221257_00127_c100560082550000001823094812221334_s1_p0/128361/872_4288" },
                                                                          TECH_PACBIO,  0, &con_pacbio_range,  4,  {1,2,3,-1},     {-1},           {-1},           {-1},           1,-1,   3,-1, -1,  0,  PX_pacbio    },
    {},  { "PacBio-Label",  { "m64136_200621_234916/18/ccs" },            TECH_PACBIO,  0, &con_pacbio_label,  3,  {1,-1},         {-1},           {-1},           {-1},           1,-1,  -1,-1, -1,  0,  PX_pacbio    },
    {},  { "PacBio-Plain",  { "m64136_200621_234916/18" },                TECH_PACBIO,  0, &con_pacbio_plain,  2,  {1,-1},         {-1},           {-1},           {-1},           1,-1,  -1,-1, -1,  0,  PX_pacbio    },
    // {},  { "Nanopore",      { "af84b0c1-6945-4323-9193-d9f6f2c38f9a" },   TECH_ONP,     0, &con_nanopore,      4,  {-1},           {-1},           {0,1,2,3,4,-1}, {-1},           -1,-1, -1,-1, -1,  36,              },
    {},  { "Nanopore",      { "af84b0c1-6945-4323-9193-d9f6f2c38f9a" },   TECH_ONP,     0, &con_nanopore,      4,  {-1},           {0,1,2,3,4-1},  {0,1,2,3,4,-1}, {0,1,2,3,4,-1}, -1,-1, -1,-1, -1,  36, PX_nanopore},
    // {},  { "Nanopore-ext",  { "2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template" },
    //                                                                       TECH_ONP,     0, &con_nanopore_ext,  5,  {-1},           {-1},           {0,1,2,3,4,-1}, {-1},           -1,-1, -1,-1, -1,                   },
    {},  { "Nanopore-ext",  { "2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template" },
                                                                          TECH_ONP,     0, &con_nanopore_ext,  5,  {-1},           {0,1,2,3,4,-1}, {0,1,2,3,4,-1}, {0,1,2,3,4,-1}, -1,-1, -1,-1, -1,  0,  PX_nanopore_ext },
    {},  { "BamSurgeon",    { "22:33597495-34324994_726956_727496_0:0:0_0:0:0_2963e" },   
                                                                          TECH_UNKNOWN, 0, &con_bamsurgeon,    7,  {1,2,3,4,-1},   {-1},           {1,3,7,-1},     {7,-1},         1,3,   2,4,   -1,                   },
    {},  { "NCBI_SRA_LEN",  { "SRR11215720.1_1_length=120" },             TECH_UNKNOWN, 0, &con_ncbi_sra_len,  10, {1,2,-1},       {-1},           {-1},           {-1},           1,-1,  2,-1,  3,   0,  PX_sra_len   },
    {},  { "NCBI-SRA2+-FQ", { "ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 someextradata length=1128" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sra2P_fq, 6,  {2,3,8,-1},     {-1},           {2,3,-1},       {-1},           3,-1,  -1,-1, 8,                    },
    {},  { "NCBI-SRA+-FQ",  { "ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template BOWDEN04_20151016_MN15199_FAA67113_BOWDEN04_MdC_MARC_Phase2a_4833_1_ch19_file1_strand length=52" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sraP_fq,  5,  {2,7,-1},       {-1},           {-1},           {-1},           2,-1,  -1,-1, 7,                    },
    {},  { "NCBI-SRA2-FQ",  { "ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 length=1128" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sra2_fq,  5,  {2,3,7,-1},     {-1},           {2,3,-1},       {-1},           3,-1,  -1,-1, 7,                    },
    {},  { "NCBI-SRA-FQ",   { "ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template length=52" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sra_fq,   4,  {2,6,-1},       {-1},           {-1},           {-1},           2,-1,  -1,-1, 6,                    },
    {},  { "NCBI-SRA2",     { "ERR2708427.1.1" },                         TECH_UNKNOWN, 0, &con_ncbi_sra2,     2,  {2,3,-1},       {-1},           {2,3,-1},       {-1},           3,-1,  -1,-1, -1,                   },
    {},  { "NCBI-SRA",      { "SRR001666.1" },                            TECH_UNKNOWN, 0, &con_ncbi_sra,      1,  {2,-1},         {-1},           {2,-1},         {-1},           2,-1,  -1,-1, -1,                   },
    {},  { "seqan",         { "adeno-reads100.fasta.000000008" },         TECH_UNKNOWN, 0, &con_seqan,         2,  {-1},           {2, -1},        {-1},           {-1},           2,-1,  -1,-1, -1,   0,  PX_seqan    },
    {},  { "CLC-GW",        { "umi64163_count1" },                        TECH_UNKNOWN, 0, &con_clc_gw,        9,  {0,1,-1},       {-1},           {0,1,-1},       {-1},           -1,-1, -1,-1, -1,   0,  PX_clc_gw   },
    {},  { "hex_chr",       { "30cf_chr10" }, /* wgsim simulator? */      TECH_UNKNOWN, 0, &con_hex_chr,       1,  {-1},           {-1},           {-1},           {0,-1},         -1,-1,  -1,-1, -1,                  }, // added v14
    {},  { "Integer",       { "123" },                                    TECH_UNKNOWN, 0, &con_integer,       0,  {0,-1},         {-1},           {0,-1},         {-1},           0,-1,  -1,-1, -1,                   }, 
    {},  { "Str_Integer",   { "read_1" },   /* eg CLC */                  TECH_UNKNOWN, 0, &con_str_integer,   1,  {1,-1},         {-1},           {1,-1},         {-1},           1,-1,  -1,-1, -1,                   },
    {},  { "Genozip-opt",   { "basic.1" },  /* must be last */            TECH_UNKNOWN, 0, &con_genozip_opt,   1,  {1,-1},         {-1},           {1,-1},         {-1},           1,-1,  -1,-1, -1,                   },

    // FASTQ Line3 QFs
         { "NCBI-SRA2+-L3", { "ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 someextradata length=1128" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sra2P_l3, 6, {2,3,7,-1},     {-1},            {2,3,-1},       {-1},           3,4,   -1,-1,  -1,                  },
         { "NCBI-SRA+-L3",  { "ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template BOWDEN04_20151016_MN15199_FAA67113_BOWDEN04_MdC_MARC_Phase2a_4833_1_ch19_file1_strand length=52" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sraP_l3,  5, {2,6,-1},       {-1},            {-1},           {-1},           2,3,   -1,-1,  -1,                  },
         { "NCBI-SRA2-L3",  { "ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 length=1128" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sra2_l3,  5, {2,3,6,-1},     {-1},            {2,3,-1},       {-1},           3,4,   -1,-1,  -1,                  },
         { "NCBI-SRA-L3",   { "ERR811170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template length=52" },
                                                                          TECH_UNKNOWN, 2, &con_ncbi_sra_l3,   4, {2,5,-1},       {-1},            {-1},           {-1},           2,3,   -1,-1,  -1,                  },
};

#define NUM_QF_L3s 4
#define NUM_QFs (ARRAY_LEN(qf) - NUM_QF_L3s)
#define FIRST_QF_L3 NUM_QFs

