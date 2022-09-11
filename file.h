// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "context.h"
#include "aes.h"
#include "data_types.h"
#include "mutex.h"
#include "buffer.h"

// REFERENCE file - only appears in .genozip format
#define REF_GENOZIP_   ".ref" GENOZIP_EXT

// VCF file variations
#define VCF_           ".vcf"
#define VCF_GZ_        ".vcf.gz"
#define VCF_BGZF_      ".vcf.bgz"
#define VCF_BZ2_       ".vcf.bz2"
#define VCF_XZ_        ".vcf.xz"
#define VCF_GENOZIP_   ".vcf" GENOZIP_EXT
#define DVCF_GENOZIP_  ".d.vcf" GENOZIP_EXT

// SAM file variations
#define SAM_           ".sam"
#define SAM_GZ_        ".sam.gz"
#define SAM_BGZF_      ".sam.bgz"
#define SAM_BZ2_       ".sam.bz2"
#define SAM_XZ_        ".sam.xz"
#define SAM_GENOZIP_   ".sam" GENOZIP_EXT

// BAM
#define BAM_           ".bam"    // can be bgzf-compressed or not
#define BAM_GZ_        ".bam.gz" // same as .bam - not observed in the wild, but used by test.sh
#define BAM_BGZF_      ".bam.bgz" 
#define CRAM_          ".cram"
#define BAM_GENOZIP_   ".bam.genozip" 

// BCF
#define BCF_           ".bcf"
#define BCF_GZ_        ".bcf.gz"
#define BCF_BGZF_      ".bcf.bgz"
#define BCF_GENOZIP_   ".bcf.genozip" 

// FASTQ file variations
#define FASTQ_         ".fastq"
#define FASTQ_GZ_      ".fastq.gz"
#define FASTQ_BZ2_     ".fastq.bz2"
#define FASTQ_XZ_      ".fastq.xz"
#define FASTQ_GENOZIP_ ".fastq" GENOZIP_EXT
#define FQ_            ".fq"
#define FQ_GZ_         ".fq.gz"
#define FQ_BZ2_        ".fq.bz2"
#define FQ_XZ_         ".fq.xz"
#define FQ_GENOZIP_    ".fq" GENOZIP_EXT

// FASTA file variations
#define FASTA_         ".fasta"
#define FASTA_GZ_      ".fasta.gz"
#define FASTA_BZ2_     ".fasta.bz2"
#define FASTA_XZ_      ".fasta.xz"
#define FASTA_GENOZIP_ ".fasta" GENOZIP_EXT
#define FA_            ".fa"
#define FA_GZ_         ".fa.gz"
#define FA_BZ2_        ".fa.bz2"
#define FA_XZ_         ".fa.xz"
#define FA_GENOZIP_    ".fa" GENOZIP_EXT
#define FAA_           ".faa"
#define FAA_GZ_        ".faa.gz"
#define FAA_BZ2_       ".faa.bz2"
#define FAA_XZ_        ".faa.xz"
#define FAA_GENOZIP_   ".faa" GENOZIP_EXT
#define FFN_           ".ffn"
#define FFN_GZ_        ".ffn.gz"
#define FFN_BZ2_       ".ffn.bz2"
#define FFN_XZ_        ".ffn.xz"
#define FFN_GENOZIP_   ".ffn" GENOZIP_EXT
#define FNN_           ".fnn"
#define FNN_GZ_        ".fnn.gz"
#define FNN_BZ2_       ".fnn.bz2"
#define FNN_XZ_        ".fnn.xz"
#define FNN_GENOZIP_   ".fnn" GENOZIP_EXT
#define FNA_           ".fna"
#define FNA_GZ_        ".fna.gz"
#define FNA_BZ2_       ".fna.bz2"
#define FNA_XZ_        ".fna.xz"
#define FNA_GENOZIP_   ".fna" GENOZIP_EXT
#define FRN_           ".frn"
#define FRN_GZ_        ".frn.gz"
#define FRN_BZ2_       ".frn.bz2"
#define FRN_XZ_        ".frn.xz"
#define FRN_GENOZIP_   ".frn" GENOZIP_EXT
#define FAS_           ".fas"
#define FAS_GZ_        ".fas.gz"
#define FAS_BZ2_       ".fas.bz2"
#define FAS_XZ_        ".fas.xz"
#define FAS_GENOZIP_   ".fas" GENOZIP_EXT

// GFF3 file variations (including GVF, which is a subtype of GFF3)
#define GFF3_          ".gff3"
#define GFF3_GZ_       ".gff3.gz"
#define GFF3_BZ2_      ".gff3.bz2"
#define GFF3_XZ_       ".gff3.xz"
#define GFF3_GENOZIP_  ".gff3" GENOZIP_EXT
#define GFF_           ".gff"
#define GFF_GZ_        ".gff.gz"
#define GFF_BZ2_       ".gff.bz2"
#define GFF_XZ_        ".gff.xz"
#define GFF_GENOZIP_   ".gff" GENOZIP_EXT
#define GVF_           ".gvf"
#define GVF_GZ_        ".gvf.gz"
#define GVF_BZ2_       ".gvf.bz2"
#define GVF_XZ_        ".gvf.xz"
#define GVF_GENOZIP_   ".gvf" GENOZIP_EXT

// 23andMe file variations
// note: 23andMe files come as a .txt, and therefore the user must specify --input to compress them. we have this
// made-up file extension here to avoid needing special cases throughout the code
#define ME23_          ".txt" // our made up extension - natively, 23andMe files come as a zip container containing a txt file
#define ME23_ZIP_      ".zip" 
#define ME23_GENOZIP_  ".txt" GENOZIP_EXT

// PHYLIP files
#define PHY_           ".phy"
#define PHY_GZ_        ".phy.gz"
#define PHY_BZ2_       ".phy.bz2"
#define PHY_XZ_        ".phy.xz"
#define PHY_GENOZIP_   ".phy" GENOZIP_EXT

// chain files
#define CHAIN_         ".chain"
#define CHAIN_GZ_      ".chain.gz"
#define CHAIN_BZ2_     ".chain.bz2"
#define CHAIN_XZ_      ".chain.xz"
#define CHAIN_GENOZIP_ ".chain" GENOZIP_EXT

// kraken files
#define KRAKEN_         ".kraken" // kraken files have no standard extension, we make one up to avoid irregularities
#define KRAKEN_GZ_      ".kraken.gz"
#define KRAKEN_BZ2_     ".kraken.bz2"
#define KRAKEN_XZ_      ".kraken.xz"
#define KRAKEN_GENOZIP_ ".kraken" GENOZIP_EXT

// locs files
#define LOCS_           ".locs" 
#define LOCS_GZ_        ".locs.gz"
#define LOCS_BZ2_       ".locs.bz2"
#define LOCS_XZ_        ".locs.xz"
#define LOCS_GENOZIP_   ".locs" GENOZIP_EXT

// generic files
#define GNRIC_           ""
#define GNRIC_GZ_        ".gz"
#define GNRIC_BZ2_       ".bz2"
#define GNRIC_XZ_        ".xz"
#define GNRIC_GENOZIP_   GENOZIP_EXT

typedef enum {TXT_FILE, Z_FILE} FileSupertype; 

typedef enum { UNKNOWN_FILE_TYPE, 
               REF_GENOZIP,
               VCF, VCF_GZ, VCF_BGZF, VCF_BZ2, VCF_XZ,   VCF_GENOZIP, 
               SAM, SAM_GZ, SAM_BGZF, SAM_BZ2, SAM_XZ,   SAM_GENOZIP,
               FASTQ, FASTQ_GZ, FASTQ_BZ2, FASTQ_XZ,     FASTQ_GENOZIP,
               FQ,    FQ_GZ,    FQ_BZ2,    FQ_XZ,        FQ_GENOZIP,
               FASTA, FASTA_GZ, FASTA_BZ2, FASTA_XZ,     FASTA_GENOZIP,
               FA,    FA_GZ,    FA_BZ2,    FA_XZ,        FA_GENOZIP,
               FAA,   FAA_GZ,   FAA_BZ2,   FAA_XZ,       FAA_GENOZIP,
               FFN,   FFN_GZ,   FFN_BZ2,   FFN_XZ,       FFN_GENOZIP,
               FNN,   FNN_GZ,   FNN_BZ2,   FNN_XZ,       FNN_GENOZIP,
               FNA,   FNA_GZ,   FNA_BZ2,   FNA_XZ,       FNA_GENOZIP,
               FRN,   FRN_GZ,   FRN_BZ2,   FRN_XZ,       FRN_GENOZIP,
               FAS,   FAS_GZ,   FAS_BZ2,   FAS_XZ,       FAS_GENOZIP,
               GFF3,  GFF3_GZ,  GFF3_BZ2,  GFF3_XZ,      GFF3_GENOZIP,
               GFF,   GFF_GZ,   GFF_BZ2,   GFF_XZ,       GFF_GENOZIP,
               GVF,   GVF_GZ,   GVF_BZ2,   GVF_XZ,       GVF_GENOZIP,
               ME23,  ME23_ZIP,                          ME23_GENOZIP, 
               PHY,   PHY_GZ,   PHY_BZ2,   PHY_XZ,       PHY_GENOZIP,
               CHAIN, CHAIN_GZ, CHAIN_BZ2, CHAIN_XZ,     CHAIN_GENOZIP,
               KRAKEN, KRAKEN_GZ, KRAKEN_BZ2, KRAKEN_XZ, KRAKEN_GENOZIP,
               LOCS,  LOCS_GZ,  LOCS_BZ2,  LOCS_XZ,      LOCS_GENOZIP,
               BAM,   BAM_GZ,   BAM_BGZF, CRAM,          BAM_GENOZIP,
               BCF,   BCF_GZ,   BCF_BGZF,                BCF_GENOZIP,  
               // the GNRIC row *must* be the last row, as it consists catch-all extensions (*.gz etc)
               GNRIC_GZ, GNRIC_BZ2, GNRIC_XZ,            GNRIC_GENOZIP, GNRIC, // GNRIC *must* be the very last as it is a catch-all ""
               AFTER_LAST_FILE_TYPE } FileType;

#define FILE_EXTS {"Unknown", /* order matches the FileType enum */ \
                   REF_GENOZIP_,                                                     \
                   VCF_, VCF_GZ_, VCF_BGZF_, VCF_BZ2_, VCF_XZ_, VCF_GENOZIP_,        \
                   SAM_, SAM_GZ_, SAM_BGZF_, SAM_BZ2_, SAM_XZ_, SAM_GENOZIP_,        \
                   FASTQ_, FASTQ_GZ_, FASTQ_BZ2_, FASTQ_XZ_, FASTQ_GENOZIP_,         \
                   FQ_,    FQ_GZ_,    FQ_BZ2_,    FQ_XZ_,    FQ_GENOZIP_,            \
                   FASTA_, FASTA_GZ_, FASTA_BZ2_, FASTA_XZ_, FASTA_GENOZIP_,         \
                   FA_,    FA_GZ_,    FA_BZ2_,    FA_XZ_,    FA_GENOZIP_,            \
                   FAA_,   FAA_GZ_,   FAA_BZ2_,   FAA_XZ_,   FAA_GENOZIP_,           \
                   FFN_,   FFN_GZ_,   FFN_BZ2_,   FFN_XZ_,   FFN_GENOZIP_,           \
                   FNN_,   FNN_GZ_,   FNN_BZ2_,   FNN_XZ_,   FNN_GENOZIP_,           \
                   FNA_,   FNA_GZ_,   FNA_BZ2_,   FNA_XZ_,   FNA_GENOZIP_,           \
                   FRN_,   FRN_GZ_,   FRN_BZ2_,   FRN_XZ_,   FRN_GENOZIP_,           \
                   FAS_,   FAS_GZ_,   FAS_BZ2_,   FAS_XZ_,   FAS_GENOZIP_,           \
                   GFF3_,  GFF3_GZ_,  GFF3_BZ2_,  GFF3_XZ_,  GFF3_GENOZIP_,          \
                   GFF_,   GFF_GZ_,   GFF_BZ2_,   GFF_XZ_,   GFF_GENOZIP_,           \
                   GVF_,   GVF_GZ_,   GVF_BZ2_,   GVF_XZ_,   GVF_GENOZIP_,           \
                   ME23_,  ME23_ZIP_,                        ME23_GENOZIP_,          \
                   PHY_,   PHY_GZ_,   PHY_BZ2_,   PHY_XZ_,   PHY_GENOZIP_,           \
                   CHAIN_, CHAIN_GZ_, CHAIN_BZ2_, CHAIN_XZ_, CHAIN_GENOZIP_,         \
                   KRAKEN_, KRAKEN_GZ_, KRAKEN_BZ2_, KRAKEN_XZ_, KRAKEN_GENOZIP_,    \
                   LOCS_,  LOCS_GZ_,  LOCS_BZ2_,  LOCS_XZ_,  LOCS_GENOZIP_,          \
                   BAM_,   BAM_GZ_,   BAM_BGZF_,  CRAM_,     BAM_GENOZIP_,           \
                   BCF_,   BCF_GZ_,   BCF_BGZF_,             BCF_GENOZIP_,           \
                   GNRIC_GZ_, GNRIC_BZ2_, GNRIC_XZ_,         GNRIC_GENOZIP_, GNRIC_, /* GNRIC_ is catch all */ \
                   "stdin", "stdout" }
extern rom file_exts[];

// Ordered by data_type: txt file types and their corresponding genozip file types for each data type
// first entry of each data type MUST be the default plain file
//                           { in          codec       out            }
#define TXT_IN_FT_BY_DT  { { { FASTA,      CODEC_NONE, REF_GENOZIP    }, { FASTA_GZ,  CODEC_GZ,  REF_GENOZIP    }, /* fasta for generating reference files */\
                             { FASTA_BZ2,  CODEC_BZ2,  REF_GENOZIP    }, { FASTA_XZ,  CODEC_XZ,  REF_GENOZIP    },\
                             { FNA,        CODEC_NONE, REF_GENOZIP    }, { FNA_GZ,    CODEC_GZ,  REF_GENOZIP    },\
                             { FNA_BZ2,    CODEC_BZ2,  REF_GENOZIP    }, { FNA_XZ,    CODEC_XZ,  REF_GENOZIP    },\
                             { FAS,        CODEC_NONE, REF_GENOZIP    }, { FAS_GZ,    CODEC_GZ,  REF_GENOZIP    },\
                             { FAS_BZ2,    CODEC_BZ2,  REF_GENOZIP    }, { FAS_XZ,    CODEC_XZ,  REF_GENOZIP    },\
                             { FA,         CODEC_NONE, REF_GENOZIP    }, { FA_GZ,     CODEC_GZ,  REF_GENOZIP    },\
                             { FA_BZ2,     CODEC_BZ2,  REF_GENOZIP    }, { FA_XZ,     CODEC_XZ,  REF_GENOZIP    }, { } }, \
                           { { VCF,        CODEC_NONE, VCF_GENOZIP    }, { VCF_GZ,    CODEC_GZ,  VCF_GENOZIP    }, { VCF_BGZF, CODEC_GZ,  VCF_GENOZIP  },\
                             { VCF_BZ2,    CODEC_BZ2,  VCF_GENOZIP    }, { VCF_XZ,    CODEC_XZ,  VCF_GENOZIP    }, { } },\
                           { { SAM,        CODEC_NONE, SAM_GENOZIP    }, { SAM_GZ,    CODEC_GZ,  SAM_GENOZIP    }, { SAM_BGZF, CODEC_GZ,  SAM_GENOZIP  },\
                             { SAM_BZ2,    CODEC_BZ2,  SAM_GENOZIP    }, { SAM_XZ,    CODEC_XZ,  SAM_GENOZIP    }, { }, },\
                           { { FASTQ,      CODEC_NONE, FASTQ_GENOZIP  }, { FASTQ_GZ,  CODEC_GZ,  FASTQ_GENOZIP  },\
                             { FASTQ_BZ2,  CODEC_BZ2,  FASTQ_GENOZIP  }, { FASTQ_XZ,  CODEC_XZ,  FASTQ_GENOZIP  },\
                             { FQ,         CODEC_NONE, FQ_GENOZIP     }, { FQ_GZ,     CODEC_GZ,  FQ_GENOZIP     },\
                             { FQ_BZ2,     CODEC_BZ2,  FQ_GENOZIP     }, { FQ_XZ,     CODEC_XZ,  FQ_GENOZIP     }, { } },\
                           { { FASTA,      CODEC_NONE, FASTA_GENOZIP  }, { FASTA_GZ,  CODEC_GZ,  FASTA_GENOZIP  },\
                             { FASTA_BZ2,  CODEC_BZ2,  FASTA_GENOZIP  }, { FASTA_XZ,  CODEC_XZ,  FASTA_GENOZIP  },\
                             { FAA,        CODEC_NONE, FAA_GENOZIP    }, { FAA_GZ,    CODEC_GZ,  FAA_GENOZIP    },\
                             { FAA_BZ2,    CODEC_BZ2,  FAA_GENOZIP    }, { FAA_XZ,    CODEC_XZ,  FAA_GENOZIP    },\
                             { FFN,        CODEC_NONE, FFN_GENOZIP    }, { FFN_GZ,    CODEC_GZ,  FFN_GENOZIP    },\
                             { FFN_BZ2,    CODEC_BZ2,  FFN_GENOZIP    }, { FFN_XZ,    CODEC_XZ,  FFN_GENOZIP    },\
                             { FNN,        CODEC_NONE, FNN_GENOZIP    }, { FNN_GZ,    CODEC_GZ,  FNN_GENOZIP    },\
                             { FNN_BZ2,    CODEC_BZ2,  FNN_GENOZIP    }, { FNN_XZ,    CODEC_XZ,  FNN_GENOZIP    },\
                             { FNA,        CODEC_NONE, FNA_GENOZIP    }, { FNA_GZ,    CODEC_GZ,  FNA_GENOZIP    },\
                             { FNA_BZ2,    CODEC_BZ2,  FNA_GENOZIP    }, { FNA_XZ,    CODEC_XZ,  FNA_GENOZIP    },\
                             { FRN,        CODEC_NONE, FRN_GENOZIP    }, { FRN_GZ,    CODEC_GZ,  FRN_GENOZIP    },\
                             { FRN_BZ2,    CODEC_BZ2,  FRN_GENOZIP    }, { FRN_XZ,    CODEC_XZ,  FRN_GENOZIP    },\
                             { FAS,        CODEC_NONE, FAS_GENOZIP    }, { FAS_GZ,    CODEC_GZ,  FAS_GENOZIP    },\
                             { FAS_BZ2,    CODEC_BZ2,  FAS_GENOZIP    }, { FAS_XZ,    CODEC_XZ,  FAS_GENOZIP    },\
                             { FA,         CODEC_NONE, FA_GENOZIP     }, { FA_GZ,     CODEC_GZ,  FA_GENOZIP     },\
                             { FA_BZ2,     CODEC_BZ2,  FA_GENOZIP     }, { FA_XZ,     CODEC_XZ,  FA_GENOZIP     }, { } },\
                           { { GFF3,       CODEC_NONE, GFF3_GENOZIP   }, { GFF3_GZ,   CODEC_GZ,  GFF3_GENOZIP   },\
                             { GFF3_BZ2,   CODEC_BZ2,  GFF3_GENOZIP   }, { GFF3_XZ,   CODEC_XZ,  GFF3_GENOZIP   },\
                             { GFF,        CODEC_NONE, GFF_GENOZIP    }, { GFF_GZ,    CODEC_GZ,  GFF_GENOZIP    },\
                             { GFF_BZ2,    CODEC_BZ2,  GFF_GENOZIP    }, { GFF_XZ,    CODEC_XZ,  GFF_GENOZIP    },\
                             { GVF,        CODEC_NONE, GVF_GENOZIP    }, { GVF_GZ,    CODEC_GZ,  GVF_GENOZIP    },\
                             { GVF_BZ2,    CODEC_BZ2,  GVF_GENOZIP    }, { GVF_XZ,    CODEC_XZ,  GVF_GENOZIP    }, { } },\
                           { { ME23,       CODEC_NONE, ME23_GENOZIP   }, { ME23_ZIP,  CODEC_ZIP, ME23_GENOZIP   }, { } },\
                           { { BAM,        CODEC_BGZF, BAM_GENOZIP    }, { BAM_GZ,    CODEC_GZ,  BAM_GENOZIP    }, { BAM_BGZF, CODEC_BGZF, BAM_GENOZIP  }, { CRAM, CODEC_CRAM, BAM_GENOZIP }, { } }, \
                           { { BCF,        CODEC_BCF,  BCF_GENOZIP    }, { BCF_GZ,    CODEC_BCF, BCF_GENOZIP    }, { BCF_BGZF, CODEC_BCF,  BCF_GENOZIP  }, { } }, \
                           { { GNRIC,      CODEC_NONE, GNRIC_GENOZIP  }, { GNRIC_GZ,  CODEC_GZ,  GNRIC_GENOZIP  },\
                             { GNRIC_BZ2,  CODEC_BZ2,  GNRIC_GENOZIP  }, { GNRIC_XZ,  CODEC_XZ,  GNRIC_GENOZIP  }, { } },\
                           { { PHY,        CODEC_NONE, PHY_GENOZIP    }, { PHY_GZ,    CODEC_GZ,  PHY_GENOZIP    },\
                             { PHY_BZ2,    CODEC_BZ2,  PHY_GENOZIP    }, { PHY_XZ,    CODEC_XZ,  PHY_GENOZIP    }, { } },\
                           { { CHAIN,      CODEC_NONE, CHAIN_GENOZIP  }, { CHAIN_GZ,  CODEC_GZ,  CHAIN_GENOZIP  },\
                             { CHAIN_BZ2,  CODEC_BZ2,  CHAIN_GENOZIP  }, { CHAIN_XZ,  CODEC_XZ,  CHAIN_GENOZIP  }, { } },\
                           { { KRAKEN,     CODEC_NONE, KRAKEN_GENOZIP }, { KRAKEN_GZ, CODEC_GZ,  KRAKEN_GENOZIP },\
                             { KRAKEN_BZ2, CODEC_BZ2,  KRAKEN_GENOZIP }, { KRAKEN_XZ, CODEC_XZ,  KRAKEN_GENOZIP }, { } },\
                           { { LOCS,       CODEC_NONE, LOCS_GENOZIP   }, { LOCS_GZ,   CODEC_GZ,  LOCS_GENOZIP   },\
                             { LOCS_BZ2,   CODEC_BZ2,  LOCS_GENOZIP   }, { LOCS_XZ,   CODEC_XZ,  LOCS_GENOZIP   }, { } },\
                        }

// Ordered by data_type: Supported output formats for genounzip
// plain file MUST appear first on the list - this will be the default output when redirecting 
// GZ file, if it is supported MUST be 2nd on the list - we use this type if the user outputs to eg xx.gz instead of xx.vcf.gz (see file_open_txt_write)
#define TXT_OUT_FT_BY_DT { { 0 }, /* a reference file cannot be uncompressed */  \
                           { VCF, VCF_GZ, VCF_BGZF, BCF, 0 },  \
                           { SAM, SAM_GZ, BAM, 0 },     \
                           { FASTQ, FASTQ_GZ, FQ, FQ_GZ, 0 }, \
                           { FASTA, FASTA_GZ, FA, FA_GZ, FAA, FAA_GZ, FFN, FFN_GZ, FNN, FNN_GZ, FNA, FNA_GZ, FRN, FRN_GZ, FAS, FAS_GZ, 0 },\
                           { GFF3, GFF3_GZ, GFF, GFF_GZ, GVF, GVF_GZ, 0 }, \
                           { ME23, ME23 /* no GZ */, ME23_ZIP, 0 }, \
                           { BAM }, /* There are no data_type=DT_BAM genozip files - .bam.genozip have data_type=DT_SAM */ \
                           { 0 }, /* There are no data_type=DT_BCF genozip files - .bam.genozip have data_type=DT_VCF */ \
                           { GNRIC, GNRIC_GZ, 0 }, \
                           { PHY, PHY_GZ, 0 }, \
                           { CHAIN, CHAIN_GZ, 0 }, \
                           { KRAKEN, KRAKEN_GZ, 0 }, \
                           { LOCS, LOCS_GZ, 0 }, \
                         }                        

// Ordered by data_type
#define Z_FT_BY_DT { { REF_GENOZIP, 0  },                   \
                     { VCF_GENOZIP, 0  },                   \
                     { SAM_GENOZIP, BAM_GENOZIP, 0 },       \
                     { FASTQ_GENOZIP, FQ_GENOZIP, 0 },      \
                     { FASTA_GENOZIP, FA_GENOZIP, FAA_GENOZIP, FFN_GENOZIP, FNN_GENOZIP, FNA_GENOZIP, FRN_GENOZIP, FAS_GENOZIP, 0 }, \
                     { GFF3_GENOZIP, GFF_GENOZIP, GVF_GENOZIP, 0  }, \
                     { ME23_GENOZIP, 0 },                   \
                     { 0 }, /* There are no data_type=DT_BAM genozip files - .bam.genozip have data_type=DT_SAM */ \
                     { 0 }, /* There are no data_type=DT_BCF genozip files - .bam.genozip have data_type=DT_VCF */ \
                     { GNRIC_GENOZIP, 0 },                  \
                     { PHY_GENOZIP, 0 },                    \
                     { CHAIN_GENOZIP, 0 },                  \
                     { KRAKEN_GENOZIP, 0 },                 \
                     { LOCS_GENOZIP, 0 },                   \
                   } 

#define MAX_GEN_COMP 2

typedef rom FileMode;
extern FileMode READ, WRITE, WRITEREAD;// this are pointers to static strings - so they can be compared eg "if (mode==READ)"

typedef struct File {
    void *file;
    char *name;                        // allocated by file_open(), freed by file_close()
    rom basename;                      // basename of name
    FileMode mode;
    FileSupertype supertype;            
    FileType type;
    bool is_remote;                    // true if file is downloaded from a url
    bool redirected;                   // txt_file: true if this file is redirected from stdin/stdout or a pipe
    bool is_eof;                       // we've read the entire file
    bool header_only;                  // ZIP txt_file: file has only the data-type header and no data
    DataType data_type;
    Codec source_codec;                // ZIP - txt_file: codec of txt file before redirection (eg CRAM, XZ, ZIP...)
    Codec codec;                       // ZIP - txt_file: generic codec used with this file (in PIZ we use flag.bgzf instead). If redirected - as read by txtfile (eg for cram files this is BGZF)

    // these relate to actual bytes on the disk
    int64_t disk_size;                 // 0 if not known (eg stdin or http stream). 
    int64_t disk_so_far;               // data read/write to/from "disk" (using fread/fwrite) (possibley gz/bz2 compressed)
    int64_t est_seggable_size;         // ZIP txt_file, access via txtfile_get_seggable_size(). Estimated size of txt_data in file, i.e. excluding the header. It is exact for plain files, or based on test_vb if the file has source compression

    // this relate to the textual data represented. In case of READ - only data that was picked up from the read buffer.
    int64_t txt_data_so_far_single;    // txt_file: data read (ZIP) or written (PIZ) to/from txt file so far
                                       // z_file: txt data represented in the GENOZIP data written (ZIP) or read (PIZ) to/from the genozip file so far for the current txt file
    int64_t header_size;               // txt_file ZIP: size of txt header  z_file ZIP: size of MAIN txt header
    
    //      seggable_data_so_far = txt_data_so_far_single - header_size
    int64_t seggable_data_so_far_gz_bz2; // txt_file ZIP: if source gz or bz2 compression, this is the compressed txt_data_so_far_single 
    int64_t txt_data_so_far_bind;      // z_file only: uncompressed txt data represented in the GENOZIP data written so far for all bound files
                                       // note: txt_data_so_far_single/bind accounts for txt modifications due to --optimize or --chain or compressing a Luft file
    int64_t txt_data_so_far_single_0;  // txt_file PIZ: accounting for data as it was in the original source file, as reading TxtHeader and VbHeader sections from the genozip file
                                       // z_file & ZIP only: same as txt_data_so_far_single/bind, but original sizes without modifications due to --chain/--optimize/Luft
    int64_t txt_data_so_far_bind_0;    // z_file & ZIP only: similar to txt_data_so_far_single_0, but for bound
    int64_t txt_disk_so_far_bind;      // z_file & ZIP only: compressed (with txt_file.codec - eg bgzf) txt data represented in the GENOZIP data written so far for all bound files
    int64_t num_lines;                 // z_file: number of lines in all txt files bound into this z_file
                                       // txt_file: number of lines, in source file terms, (read so far) in single txt file
    // per-component data (ZIP)
    int64_t disk_so_far_comp[MAX_NUM_COMPS];
    int64_t txt_data_so_far_bind_comp[MAX_NUM_COMPS];
    int64_t txt_data_so_far_bind_0_comp[MAX_NUM_COMPS];

    // Digest stuff - stored in z_file (ZIP & PIZ)
    Serializer digest_serializer;      // ZIP/PIZ: used for serializing VBs so they are MD5ed in order (not used for Adler32)
    DigestContext digest_ctx;          // ZIP/PIZ: z_file: digest context of txt file being compressed / reconstructed (used for MD5 and, in v9-13, for Adler32)
    Digest digest;                     // ZIP: z_file: digest of txt data read from input file  PIZ: z_file: as read from TxtHeader section (used for MD5 and, in v9-13, for Adler32)

    // Used for READING & WRITING txt files - but stored in the z_file structure for zip to support bindenation (and in the txt_file structure for piz)
    uint32_t max_lines_per_vb;         // ZIP & PIZ - in ZIP, discovered while segmenting, in PIZ - given by SectionHeaderTxtHeader
    bool piz_header_init_has_run;      // PIZ: true if we called piz_header_init to initialize (only once per outputted txt_file, even if concatenated)

    // Used for READING GENOZIP files
    uint8_t genozip_version;           // GENOZIP_FILE_FORMAT_VERSION of the genozip file being read
        
    CompIType num_components;          // PIZ z_file: set from genozip header 

    struct FlagsGenozipHeader z_flags; // z_file: genozip file flags as read from SectionHeaderGenozipHeader.h.flags
    struct FlagsTxtHeader txt_flags;   // txt_file PIZ: genozip file flags as read from SectionHeaderTxtHeader.h.flags

    Buffer vb_sections_index;          // PIZ z_file: an index into VB sections
    Buffer comp_sections_index;        // PIZ z_file: an index into Txt sections
    Section first_dict_section;        // PIZ z_file

    // Used for WRITING GENOZIP files
    uint64_t disk_at_beginning_of_this_txt_file;     // z_file: the value of disk_size when starting to read this txt file
    
    uint8_t txt_header_enc_padding[AES_BLOCKLEN-1];  // just so we can overwrite txt_header with encryption padding

    SectionHeaderTxtHeader txt_header_single;        // ZIP: store the txt header of single component in bound mode
                                                     // PIZ: z_file: header of COMP_MAIN
    uint8_t txt_header_enc_padding2[AES_BLOCKLEN-1]; // same

    bool z_closes_after_me;            // Z_FILE ZIP: this z_file will close after completing the current txt_file
    
    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    Mutex dicts_mutex;                 // this mutex protects contexts and num_contexts from concurrent adding of a new dictionary
    Did num_contexts;                  // length of populated subfield_ids and mtx_ctx;
    
    Did dict_id_to_did_i_map[65536 * 2]; // map for quick look up of did_i from dict_id : 64K for key_map, 64K for alt_map 
    Context contexts[MAX_DICTS];       // a merge of dictionaries of all VBs
    Buffer ra_buf;                     // ZIP/PIZ:  RAEntry records: ZIP: of DC_PRIMARY ; PIZ - PRIMARY or LUFT depending on flag.luft
    Buffer ra_buf_luft;                // ZIP only: RAEntry records of DC_LUFT
    
    Buffer chrom2ref_map;              // ZIP/PIZ: mapping from user file chrom to alternate chrom in reference file
    Buffer ref2chrom_map;              // ZIP: an inverted chrom2ref_map, created by ref_compress_ref

    // section list - used for READING and WRITING genozip files
    Buffer section_list_buf;           // section list to be written as the payload of the genotype header section
    Buffer section_list_vb1;           // a copy of the section_list of vb_i=1, bc it might be removed from section_list_buf due to recon plan.
    CompIType num_txts_so_far;         // ZIP z_file: number of txt files compressed into this z_file - each becomes a component
                                       // PIZ z_file: number of txt files written from this z_file - each generated from one or more components 

    // TXT file: reading
    Buffer unconsumed_txt;             // ZIP: excess data read from the txt file - moved to the next VB

    // TXT file: stuff reading and writing txt files compressed with BGZF
    Buffer unconsumed_bgzf_blocks;     // ZIP: unconsumed or partially consumed bgzf blocks - moved to the next VB
    Buffer bgzf_isizes;                // ZIP/PIZ: uncompressed size of the bgzf blocks in which this txt file is compressed (in BGEN16)
    Buffer bgzf_starts;                // ZIP: offset in txt_file of each BGZF block
    struct FlagsBgzf bgzf_flags;       // correspond to SectionHeader.flags in SEC_BGZF
    uint8_t bgzf_signature[3];         // PIZ: 3 LSB of size of source BGZF-compressed file, as passed in SectionHeaderTxtHeader.codec_info

    // TXT file: data used in --show-sex, --show-coverage and --idxstats
    Buffer coverage;
    Buffer read_count;
    Buffer unmapped_read_count;
    
    // Z_FILE: SAM/BAM SA stuff
    Buffer sag_grps;                   // Z_FILE: an SA group is a group of alignments, including the primary aligngment
    Buffer sag_gps_index;              // Z_FILE: index z_file->sag_grps by adler32(qname)
    Buffer sag_alns;                   // Z_FILE: array of {RNAME, STRAND, POS, CIGAR, NM, MAPQ} of the alignment
    Buffer sag_qnames;                 // Z_FILE
    Buffer sag_depn_index;             // Z_FILE: SAG_BY_FLAG: uniq-sorted hash(QNAME) of all depn alignments in the file
    
    union {
    Buffer sag_cigars;                 // Z_FILE: SAG_BY_SA: compressed CIGARs
    Buffer solo_data;                  // Z_FILE: SAG_BY_SOLO: solo data
    };
    Buffer sag_seq;                    // Z_FILE: bitmap of seqs in ACGT 2bit format
    Buffer sag_qual;                   // Z_FILE: compressed QUAL
    
    uint64_t saggy_near_count, mate_line_count, prim_far_count; // Z_FILE ZIP: SAM: for stats

    // Z_FILE: DVCF stuff
    Buffer rejects_report;             // Z_FILE ZIP --chain: human readable report about rejects
    Buffer apriori_tags;               // Z_FILE ZIP DVCF: used for INFO/FORMAT tag renaming. Data from command line options if --chain, or VCF header if DVCF
    
    // TXT_FILE: DVCF stuff
    uint8_t coords;                    // TXT_FILE ZIP: DVCF: Set from ##dual_coordinates and immutable thereafter
    uint64_t reject_bytes;             // ZIP DVCF: number of bytes in lines originating from ##primary_only/##luft_only, not yet assigned to a VB

    // Reconstruction plan, for reconstructing in sorted order if --sort: [0] is primary coords, [1] is luft coords
    Mutex recon_plan_mutex[2];         // TXT_FILE ZIP: VCF: protect vb_info and line_info during merging of VB data
    Buffer vb_info[2];                 // TXT_FILE ZIP: VCF: array of ZipVbInfo per VB, indexed by (vb_i-1), 0:PRIMARY, 1:LUFT
                                       // Z_FILE   ZIP: SAM: array of SamGcVbInfo for: 0:PRIM 1:DEPN
                                       // Z_FILE   PIZ: [0]: used by writer [1]: used to load SAM SA Groups - array of PlsgVbInfo - entry per PRIM vb 
    Buffer line_info[2];               // TXT_FILE ZIP: VCF: array of LineInfo per line or gapless range in txt_file
                                       //               SAM: array of uint32 - lengths of lines in PRIM/DEPN
    Buffer recon_plan;                 // TXT_FILE ZIP/PIZ: array of ReconPlanItem - order of reconstruction of ranges of lines, to achieve a sorted file. VCF: [0]=PRIM rendition [1]=LUFT rendition
                                       // Z_FILE   PIZ: plan for entire z_file, txt_file.recon_plan is assigned a portion of this plan
    Buffer recon_plan_index;           // TXT_FILE ZIP / Z_FILE PIZ: An array of BufWord, one for each VB: start and length of VB in recon_plan
    Buffer txt_header_info;            // Z_FILE   PIZ: used by writer
    uint64_t lines_written_so_far;     // TXT_FILE PIZ: number of textual lines (excluding the header) that passed all filters except downsampling, and is to be written to txt_file, or downsampled-out
    
    // Z_FILE 
    Buffer bound_txt_names;            // ZIP: Stats data: a concatenation of all bound txt_names that contributed to this genozip file
    
    // Information content stats - how many bytes and how many sections does this file have in each section type
    uint32_t num_vbs;                  // ZIP: z_file/txt_file PIZ: txt_file: number of VBs processed z_file: total VBs in file
    uint32_t num_vbs_dispatched;       // ZIP: txt_file
    uint32_t num_preproc_vbs_joined;   // PIZ: z_file
    uint32_t max_conc_writing_vbs;     // PIZ z_file: the maximal value conc_writing_vbs across all SEC_RECON_PLAN sections in this z_file
} File;

#define z_has_gencomp z_file->z_flags.has_gencomp
#define z_is_dvcf (Z_DT(DT_VCF) && z_has_gencomp)
#define z_sam_gencomp ((Z_DT(DT_SAM) || Z_DT(DT_BAM)) && z_has_gencomp) // note: is BAM file in piz are Z_DT(DT_SAM) and in zip are Z_DT(DT_BAM)

extern DataType last_z_dt; // data_type of last z_file opened

// methods
extern File *file_open (rom filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE */);
extern void file_close (FileP *file_p, bool index_txt, bool cleanup_memory /* optional */);
extern void file_write (FileP file, const void *data, unsigned len);
extern bool file_seek (File *file, int64_t offset, int whence, int soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
extern int64_t file_tell_do (File *file, bool soft_fail, FUNCLINE);
#define file_tell(file,soft_fail) file_tell_do ((file), (soft_fail), __FUNCLINE) 
extern void file_set_input_type (rom type_str);
extern void file_set_input_size (rom size_str);
extern FileType file_get_type (rom filename);
extern FileType file_get_stdin_type (void);
extern DataType file_get_data_type (FileType ft, bool is_input);
extern Codec file_get_codec_by_txt_ft (DataType dt, FileType txt_ft, bool source);
extern rom file_plain_text_ext_of_dt (DataType dt);
extern void file_get_raw_name_and_type (char *filename, char **raw_name, FileType *ft);
extern bool file_has_ext (rom filename, rom extension);
extern rom file_basename (rom filename, bool remove_exe, rom default_basename, char *basename, unsigned basename_size);
extern void file_assert_ext_decompressor (void);
extern void file_kill_external_compressors (void);
extern FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft);
extern DataType file_get_dt_by_z_ft (FileType z_ft);
extern FileType file_get_z_ft_by_dt (DataType dt);
extern rom file_plain_ext_by_dt (DataType dt);
extern void file_remove_codec_ext (char *filename, FileType ft);
extern rom ft_name (FileType ft);
extern rom file_guess_original_filename (const File *file);
extern char *file_get_fastq_pair_filename (rom fn1, rom fn2, bool test_only);

// wrapper operations for operating system files
extern void file_get_file (VBlockP vb, rom filename, BufferP buf, rom buf_name, bool add_string_terminator);
extern bool file_put_data (rom filename, const void *data, uint64_t len, mode_t mode);
extern void file_put_data_abort (void);

typedef struct { char s[64]; } PutLineFn;
extern PutLineFn file_put_line (VBlockP vb, STRp(line), rom msg);
extern bool file_exists (rom filename);
extern bool file_is_fifo (rom filename);
extern bool file_is_dir (rom filename);
extern void file_mkdir (rom dirname);
extern void file_remove (rom filename, bool fail_quietly);
extern void file_mkfifo (rom filename);
extern uint64_t file_get_size (rom filename);
extern char *file_compressible_extensions (bool plain_only);
extern char *file_make_unix_filename (char *filename);

#define FILENAME_STDIN  "(stdin)"
#define FILENAME_STDOUT "(stdout)"

#define file_printname(file) (!file             ? "(none)"       \
                           : (file)->name       ? (file)->name   \
                           : (file)->mode==READ ? FILENAME_STDIN \
                           :                      FILENAME_STDOUT)
#define txt_name file_printname(txt_file)
#define z_name   file_printname(z_file)

#define CLOSE(fd,name,quiet) do { ASSERTW (!close (fd) || (quiet),  "Warning in %s:%u: Failed to close %s: %s",  __FUNCLINE, (name), strerror(errno));} while (0)
#define FCLOSE(fp,name) do { if (fp) { ASSERTW (!fclose (fp), "Warning in %s:%u: Failed to fclose %s: %s", __FUNCLINE, (name), strerror(errno)); fp = NULL; } } while (0)
 
 // ---------------------------
// tests for compression types
// ---------------------------

static inline bool file_is_read_via_ext_decompressor(ConstFileP file) { 
    return file->supertype == TXT_FILE && (file->codec == CODEC_XZ || file->codec == CODEC_ZIP || file->codec == CODEC_BCF || file->codec == CODEC_CRAM);
}

static inline bool file_is_read_via_int_decompressor(ConstFileP file) {
    return file->supertype == TXT_FILE && (file->codec == CODEC_GZ || file->codec == CODEC_BGZF || file->codec == CODEC_BZ2);
}

static inline bool file_is_written_via_ext_compressor(ConstFileP file) {
    return file->supertype == TXT_FILE && (file->codec == CODEC_BCF || file->codec == CODEC_GZ);
}

static inline bool file_is_plain_or_ext_decompressor(ConstFileP file) {
    return file->supertype == TXT_FILE && (file->codec == CODEC_NONE || file_is_read_via_ext_decompressor(file));
}

// read the contents and a newline-separated text file and split into lines - creating char **lines, unsigned *line_lens and unsigned n_lines
#define file_split_lines(fn, name) \
    static Buffer data = EMPTY_BUFFER; \
    ASSINP0 (!data.len, "only one instance of a " name " option can be used"); \
    file_get_file (evb, fn, &data, "file_split_lines__" name, false); \
    \
    str_split_enforce (data.data, data.len, 0, '\n', line, true, (name)); \
    ASSINP (!line_lens[n_lines-1], "Expecting %s to end with a newline", (fn)); \
    \
    n_lines--; /* remove final empty line */\
    str_remove_CR (line); /* note: fopen with non-binary mode (no "b") doesn't help, because it won't remove Windows-created \r when running on Unix */ \
    str_nul_separate (line); \
    /* note: its up to the caller to free data */

// platform compatibility stuff
#ifdef _WIN32
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#ifdef __APPLE__
#define fseeko64 fseeko
#define ftello64 ftello
#endif
