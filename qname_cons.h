// ------------------------------------------------------------------
//   qname_cons.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "container.h"
#include "dict_id_gen.h"

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
                   { .dict_id = { _SAM_Q7NAME },                           }, // entire part after space is segged together
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } } 
};

//--------------------------------------------------------------------------------------------------------------
// BGI_E format: 27 characters: E, FlowCellSerialNumber[9], L, Lane[1], C, Column[3], R, Row[3] Tile[7]
// Example: E100020409L1C001R0030000234
// See: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_bgi_E = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 9 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } }, // note: CI0_FIXED_0_PAD is only useful if the field is segged as a number (which may be shorter than the field length). if it is segged as a text, the separator will have no effect.
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q4NAME }, .separator = { CI0_FIXED_0_PAD, 7 } },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } } 
};

//--------------------------------------------------------------------------------------------------------------------
// BGI_CL format: CL, FlowCellSerialNumber[9], L, Lane[1], C, Column[3], R, Row[3], _, variable-length number
//  CL100025298L1C002R050_244547 reported by https://en.wikipedia.org/wiki/File:BGI_seq_platform_read_name_description.png
// @CL100072652L2C001R001_12/1 for FASTQ reported by https://www.biostars.org/p/395715/
//--------------------------------------------------------------------------------------------------------------------
static SmallContainer con_bgi_CL = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 9 } }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = { CI0_FIXED_0_PAD, 1 } },
                             { .dict_id = { _SAM_Q2NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q3NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },
                             { .dict_id = { _SAM_Q4NAME },                                     },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP           } } } 
};

//--------------------------------------------------------------------------------------------------------------
// ION_TORRENT_3 format: 17 characters: RunID[5] : Row[5] : Column[5]
// Example: ZEWTM:10130:07001
// See: https://www.biostars.org/p/330644/
//--------------------------------------------------------------------------------------------------------------
static SmallContainer con_ion_torrent_3 = {
    .repeats   = 1,
    .nitems_lo = 4,
    .items     = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"          },  
                   { .dict_id = { _SAM_Q1NAME }, .separator = ":"          },
                   { .dict_id = { _SAM_Q2NAME },                           },
                   { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

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

// See: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats
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

// example: af84b0c1-6945-4323-9193-d9f6f2c38f9a
static SmallContainer con_nanopore = {
    .repeats             = 1,
    .nitems_lo           = 6,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "-"          },
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q4NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// example: 2a228edf-d8bc-45d4-9c96-3d613b8530dc_Basecall_2D_000_template
static SmallContainer con_nanopore_ext = {
    .repeats             = 1,
    .nitems_lo           = 7,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q1NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q2NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q3NAME }, .separator = "-"          }, 
                             { .dict_id = { _SAM_Q4NAME }, .separator = "_"          }, 
                             { .dict_id = { _SAM_Q5NAME }                            },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP } } }
};

// example: ERR1111170.1
static SmallContainer con_ncbi_sra = {
    .repeats             = 1,
    .nitems_lo           = 4,
    .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } }, // three letters - usually all-the-same
                             { .dict_id = { _SAM_Q1NAME }, .separator = "."                    }, // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q2NAME }                                      },
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP }           } }
};

// example: ERR1111170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template length=52
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

// example: ERR1111170.1 07dc4948-eb0c-45f2-9b40-a933a9bd5cf7_Basecall_2D_000_template BOWDEN04_20151016_MN15199_FAA67113_BOWDEN04_MdC_MARC_Phase2a_4833_1_ch19_file1_strand length=52
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
                             { .dict_id = { _SAM_Q2NAME }, .separator = "."                    },  // first number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _SAM_Q3NAME }                                      },  // 2nd number 
                             { .dict_id = { _SAM_QmNAME }, .separator = { CI0_SKIP }           } }
};

// example: ERR2708427.1.1 51e7525d-fa50-4b1a-ad6d-4f4ae25c1df7 length=1128
static SmallContainer con_ncbi_sra2_fq = {
    .repeats             = 1,
    .nitems_lo           = 8,
    .items               = { { .dict_id = { _FASTQ_Q0NAME }, .separator = { CI0_FIXED_0_PAD, 3 } },  // three letters - usually all-the-same
                             { .dict_id = { _FASTQ_Q1NAME }, .separator = "."                    },  // SRA number - usually all-the-same - seg as snip and not numeric
                             { .dict_id = { _FASTQ_Q2NAME }, .separator = "."                    },  // first number - usually all-the-same - seg as snip and not numeric
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
                             { .dict_id = { _FASTQ_Q2NAME }, .separator = "."                    },  // first number - usually all-the-same - seg as snip and not numeric
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
