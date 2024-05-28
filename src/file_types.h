// ------------------------------------------------------------------
//   file_types.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

// REFERENCE file - only appears in .genozip format
#define REF_GENOZIP_   ".ref" GENOZIP_EXT

// VCF file variations
#define VCF_           ".vcf"
#define VCF_GZ_        ".vcf.gz"
#define VCF_BGZF_      ".vcf.bgz"
#define VCF_BZ2_       ".vcf.bz2"
#define VCF_XZ_        ".vcf.xz"
#define VCF_GENOZIP_   ".vcf" GENOZIP_EXT

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
#define BAM_GENOZIP_   ".bam.genozip" 
#define CRAM_          ".cram"
#define CRAM_GENOZIP_  ".cram.genozip" 

#define DEEP_GENOZIP_  ".deep.genozip"

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
#define FASTQ_ORA_     ".fastq.ora"
#define FASTQ_GENOZIP_ ".fastq" GENOZIP_EXT
#define FQ_            ".fq"
#define FQ_GZ_         ".fq.gz"
#define FQ_BZ2_        ".fq.bz2"
#define FQ_XZ_         ".fq.xz"
#define FQ_ORA_        ".fq.ora"
#define FQ_GENOZIP_    ".fq" GENOZIP_EXT

// FASTA file variations
#define FASTA_         ".fasta"
#define FASTA_GZ_      ".fasta.gz"
#define FASTA_BZ2_     ".fasta.bz2"
#define FASTA_XZ_      ".fasta.xz"
#define FASTA_ZIP_     ".fasta.zip"
#define FASTA_GENOZIP_ ".fasta" GENOZIP_EXT
#define FA_            ".fa"
#define FA_GZ_         ".fa.gz"
#define FA_BZ2_        ".fa.bz2"
#define FA_XZ_         ".fa.xz"
#define FA_ZIP_        ".fa.zip"
#define FA_GENOZIP_    ".fa" GENOZIP_EXT
#define FAA_           ".faa"
#define FAA_GZ_        ".faa.gz"
#define FAA_BZ2_       ".faa.bz2"
#define FAA_XZ_        ".faa.xz"
#define FAA_ZIP_       ".faa.zip"
#define FAA_GENOZIP_   ".faa" GENOZIP_EXT
#define FFN_           ".ffn"
#define FFN_GZ_        ".ffn.gz"
#define FFN_BZ2_       ".ffn.bz2"
#define FFN_XZ_        ".ffn.xz"
#define FFN_ZIP_       ".ffn.zip"
#define FFN_GENOZIP_   ".ffn" GENOZIP_EXT
#define FNN_           ".fnn"
#define FNN_GZ_        ".fnn.gz"
#define FNN_BZ2_       ".fnn.bz2"
#define FNN_XZ_        ".fnn.xz"
#define FNN_ZIP_       ".fnn.zip"
#define FNN_GENOZIP_   ".fnn" GENOZIP_EXT
#define FNA_           ".fna"
#define FNA_GZ_        ".fna.gz"
#define FNA_BZ2_       ".fna.bz2"
#define FNA_XZ_        ".fna.xz"
#define FNA_ZIP_       ".fna.zip"
#define FNA_GENOZIP_   ".fna" GENOZIP_EXT
#define FRN_           ".fna"
#define FRN_GZ_        ".frn.gz"
#define FRN_BZ2_       ".frn.bz2"
#define FRN_XZ_        ".frn.xz"
#define FRN_ZIP_       ".frn.zip"
#define FRN_GENOZIP_   ".frn" GENOZIP_EXT
#define FAS_           ".fas"
#define FAS_GZ_        ".fas.gz"
#define FAS_BZ2_       ".fas.bz2"
#define FAS_XZ_        ".fas.xz"
#define FAS_ZIP_       ".fas.zip"
#define FAS_GENOZIP_   ".fas" GENOZIP_EXT
#define FSA_           ".fsa"
#define FSA_GZ_        ".fsa.gz"
#define FSA_BZ2_       ".fsa.bz2"
#define FSA_XZ_        ".fsa.xz"
#define FSA_ZIP_       ".fsa.zip"
#define FSA_GENOZIP_   ".fsa" GENOZIP_EXT // adding more FASTA extensions? also update txtheader_piz_get_filename()

// GFF2 / GFF3 file variations (including GVF, which is a subtype of GFF3 and GTF which is a subtype of GFF2)
#define GFF3_          ".gff3"
#define GFF3_GZ_       ".gff3.gz"
#define GFF3_BZ2_      ".gff3.bz2"
#define GFF3_XZ_       ".gff3.xz"
#define GFF3_GENOZIP_  ".gff3" GENOZIP_EXT
#define GFF_           ".gff" // can be GFF1 (not supported), GFF2 or GFF3
#define GFF_GZ_        ".gff.gz"
#define GFF_BZ2_       ".gff.bz2"
#define GFF_XZ_        ".gff.xz"
#define GFF_GENOZIP_   ".gff" GENOZIP_EXT
#define GVF_           ".gvf"
#define GVF_GZ_        ".gvf.gz"
#define GVF_BZ2_       ".gvf.bz2"
#define GVF_XZ_        ".gvf.xz"
#define GVF_GENOZIP_   ".gvf" GENOZIP_EXT
#define GTF_           ".gtf"
#define GTF_GZ_        ".gtf.gz"
#define GTF_BZ2_       ".gtf.bz2"
#define GTF_XZ_        ".gtf.xz"
#define GTF_GENOZIP_   ".gtf" GENOZIP_EXT

// 23andMe file variations
// note: 23andMe files come as a .txt, and therefore the user must specify --input to compress them. we have this
// made-up file extension here to avoid needing special cases throughout the code
#define ME23_          ".txt" // our made up extension - natively, 23andMe files come as a zip container containing a txt file
#define ME23_ZIP_      ".zip" 
#define ME23_GENOZIP_  ".txt" GENOZIP_EXT

// locs files
#define LOCS_           ".locs" 
#define LOCS_GZ_        ".locs.gz"
#define LOCS_BZ2_       ".locs.bz2"
#define LOCS_XZ_        ".locs.xz"
#define LOCS_GENOZIP_   ".locs" GENOZIP_EXT

// bed files
#define BED_            ".bed" 
#define BED_GZ_         ".bed.gz"
#define BED_BZ2_        ".bed.bz2"
#define BED_XZ_         ".bed.xz"
#define BED_GENOZIP_    ".bed" GENOZIP_EXT
#define TRACK_          ".track" 
#define TRACK_GZ_       ".track.gz"
#define TRACK_BZ2_      ".track.bz2"
#define TRACK_XZ_       ".track.xz"
#define TRACK_GENOZIP_  ".track" GENOZIP_EXT

// generic files
#define GNRIC_           ""
#define GNRIC_GZ_        ".gz"
#define GNRIC_BZ2_       ".bz2"
#define GNRIC_XZ_        ".xz"
#define GNRIC_ZIP_       ".zip"
#define GNRIC_GENOZIP_   GENOZIP_EXT

typedef enum {TXT_FILE, Z_FILE} FileSupertype; 

typedef enum FileType { 
    UNKNOWN_FILE_TYPE, 
    REF_GENOZIP,
    VCF, VCF_GZ, VCF_BGZF, VCF_BZ2, VCF_XZ,                 VCF_GENOZIP, 
    SAM, SAM_GZ, SAM_BGZF, SAM_BZ2, SAM_XZ,                 SAM_GENOZIP,
    FASTQ, FASTQ_GZ, FASTQ_BZ2, FASTQ_XZ,     FASTQ_ORA,    FASTQ_GENOZIP,
    FQ,    FQ_GZ,    FQ_BZ2,    FQ_XZ,        FQ_ORA,       FQ_GENOZIP,
    FASTA, FASTA_GZ, FASTA_BZ2, FASTA_XZ,     FASTA_ZIP,    FASTA_GENOZIP,
    FA,    FA_GZ,    FA_BZ2,    FA_XZ,        FA_ZIP,       FA_GENOZIP,
    FAA,   FAA_GZ,   FAA_BZ2,   FAA_XZ,       FAA_ZIP,      FAA_GENOZIP,
    FFN,   FFN_GZ,   FFN_BZ2,   FFN_XZ,       FFN_ZIP,      FFN_GENOZIP,
    FNN,   FNN_GZ,   FNN_BZ2,   FNN_XZ,       FNN_ZIP,      FNN_GENOZIP,
    FNA,   FNA_GZ,   FNA_BZ2,   FNA_XZ,       FNA_ZIP,      FNA_GENOZIP,
    FRN,   FRN_GZ,   FRN_BZ2,   FRN_XZ,       FRN_ZIP,      FRN_GENOZIP,
    FAS,   FAS_GZ,   FAS_BZ2,   FAS_XZ,       FAS_ZIP,      FAS_GENOZIP,
    FSA,   FSA_GZ,   FSA_BZ2,   FSA_XZ,       FSA_ZIP,      FSA_GENOZIP,
    GFF3,  GFF3_GZ,  GFF3_BZ2,  GFF3_XZ,                    GFF3_GENOZIP,
    GFF,   GFF_GZ,   GFF_BZ2,   GFF_XZ,                     GFF_GENOZIP,
    GVF,   GVF_GZ,   GVF_BZ2,   GVF_XZ,                     GVF_GENOZIP,
    GTF,   GTF_GZ,   GTF_BZ2,   GTF_XZ,                     GTF_GENOZIP,
    ME23,  ME23_ZIP,                                        ME23_GENOZIP, 
    LOCS,  LOCS_GZ,  LOCS_BZ2,  LOCS_XZ,                    LOCS_GENOZIP,
    BAM,   BAM_GZ,   BAM_BGZF,                              BAM_GENOZIP,
    CRAM,                                                   CRAM_GENOZIP,
    BCF,   BCF_GZ,   BCF_BGZF,                              BCF_GENOZIP,  
    BED,   BED_GZ,   BED_BZ2,   BED_XZ,                     BED_GENOZIP,
    TRACK, TRACK_GZ, TRACK_BZ2, TRACK_XZ,                   TRACK_GENOZIP,
    GNRIC_GZ, GNRIC_BZ2, GNRIC_XZ, GNRIC_ZIP, GNRIC_GENOZIP, GNRIC, // the GNRIC row *must* be the last row, as it consists catch-all extensions (*.gz etc), and the GNRIC file type *must* be last this row
    AFTER_LAST_FILE_TYPE 
} FileType;

#define FILE_EXTS {                                                                              \
    "Unknown", /* order matches the FileType enum */                                             \
                                                                         REF_GENOZIP_,           \
    VCF_,      VCF_GZ_, VCF_BGZF_, VCF_BZ2_,    VCF_XZ_,                 VCF_GENOZIP_,           \
    SAM_,      SAM_GZ_, SAM_BGZF_, SAM_BZ2_,    SAM_XZ_,                 SAM_GENOZIP_,           \
    FASTQ_,    FASTQ_GZ_,          FASTQ_BZ2_,  FASTQ_XZ_,   FASTQ_ORA_, FASTQ_GENOZIP_,         \
    FQ_,       FQ_GZ_,             FQ_BZ2_,     FQ_XZ_,      FQ_ORA_,    FQ_GENOZIP_,            \
    FASTA_,    FASTA_GZ_,          FASTA_BZ2_,  FASTA_XZ_,   FASTA_ZIP_, FASTA_GENOZIP_,         \
    FA_,       FA_GZ_,             FA_BZ2_,     FA_XZ_,      FA_ZIP_,    FA_GENOZIP_,            \
    FAA_,      FAA_GZ_,            FAA_BZ2_,    FAA_XZ_,     FAA_ZIP_,   FAA_GENOZIP_,           \
    FFN_,      FFN_GZ_,            FFN_BZ2_,    FFN_XZ_,     FFN_ZIP_,   FFN_GENOZIP_,           \
    FNN_,      FNN_GZ_,            FNN_BZ2_,    FNN_XZ_,     FNN_ZIP_,   FNN_GENOZIP_,           \
    FNA_,      FNA_GZ_,            FNA_BZ2_,    FNA_XZ_,     FNA_ZIP_,   FNA_GENOZIP_,           \
    FRN_,      FRN_GZ_,            FRN_BZ2_,    FRN_XZ_,     FRN_ZIP_,   FRN_GENOZIP_,           \
    FAS_,      FAS_GZ_,            FAS_BZ2_,    FAS_XZ_,     FAS_ZIP_,   FAS_GENOZIP_,           \
    FSA_,      FSA_GZ_,            FSA_BZ2_,    FSA_XZ_,     FSA_ZIP_,   FSA_GENOZIP_,           \
    GFF3_,     GFF3_GZ_,           GFF3_BZ2_,   GFF3_XZ_,                GFF3_GENOZIP_,          \
    GFF_,      GFF_GZ_,            GFF_BZ2_,    GFF_XZ_,                 GFF_GENOZIP_,           \
    GVF_,      GVF_GZ_,            GVF_BZ2_,    GVF_XZ_,                 GVF_GENOZIP_,           \
    GTF_,      GTF_GZ_,            GTF_BZ2_,    GTF_XZ_,                 GTF_GENOZIP_,           \
    ME23_,                                                   ME23_ZIP_,  ME23_GENOZIP_,          \
    LOCS_,     LOCS_GZ_,           LOCS_BZ2_,   LOCS_XZ_,                LOCS_GENOZIP_,          \
    BAM_,      BAM_GZ_, BAM_BGZF_,                                       BAM_GENOZIP_,           \
    CRAM_,                                                               CRAM_GENOZIP_,          \
    BCF_,      BCF_GZ_, BCF_BGZF_,                                       BCF_GENOZIP_,           \
    BED_,      BED_GZ_,            BED_BZ2_,    BED_XZ_,                 BED_GENOZIP_,           \
    TRACK_,    TRACK_GZ_,          TRACK_BZ2_,  TRACK_XZ_,               TRACK_GENOZIP_,         \
               GNRIC_GZ_,          GNRIC_BZ2_,  GNRIC_XZ_,   GNRIC_ZIP_, GNRIC_GENOZIP_, GNRIC_, /* GNRIC_ is catch all */ \
    "stdin", "stdout"                                                                            \
}
extern rom file_exts[];

// txt file types and their corresponding genozip file types for each data type
// first entry of each data type MUST be the default plain file
//                                        { in          codec       out            }
#define TXT_IN_FT_BY_DT  { /*txt_in_ft_by_dt*/ \
    [DT_REF]   = { { FASTA,      CODEC_NONE, REF_GENOZIP    }, { FASTA_GZ,  CODEC_GZ,  REF_GENOZIP    }, /* fasta for generating reference files */\
                   { FASTA_BZ2,  CODEC_BZ2,  REF_GENOZIP    }, { FASTA_XZ,  CODEC_XZ,  REF_GENOZIP    },\
                   { FNA,        CODEC_NONE, REF_GENOZIP    }, { FNA_GZ,    CODEC_GZ,  REF_GENOZIP    },\
                   { FNA_BZ2,    CODEC_BZ2,  REF_GENOZIP    }, { FNA_XZ,    CODEC_XZ,  REF_GENOZIP    },\
                   { FAS,        CODEC_NONE, REF_GENOZIP    }, { FAS_GZ,    CODEC_GZ,  REF_GENOZIP    },\
                   { FAS_BZ2,    CODEC_BZ2,  REF_GENOZIP    }, { FAS_XZ,    CODEC_XZ,  REF_GENOZIP    },\
                   { FSA,        CODEC_NONE, REF_GENOZIP    }, { FSA_GZ,    CODEC_GZ,  REF_GENOZIP    },\
                   { FSA_BZ2,    CODEC_BZ2,  REF_GENOZIP    }, { FSA_XZ,    CODEC_XZ,  REF_GENOZIP    },\
                   { FA,         CODEC_NONE, REF_GENOZIP    }, { FA_GZ,     CODEC_GZ,  REF_GENOZIP    },\
                   { FA_BZ2,     CODEC_BZ2,  REF_GENOZIP    }, { FA_XZ,     CODEC_XZ,  REF_GENOZIP    }, { } }, \
    [DT_VCF]   = { { VCF,        CODEC_NONE, VCF_GENOZIP    }, { VCF_GZ,    CODEC_GZ,  VCF_GENOZIP    }, { VCF_BGZF,  CODEC_GZ,  VCF_GENOZIP  },\
                   { VCF_BZ2,    CODEC_BZ2,  VCF_GENOZIP    }, { VCF_XZ,    CODEC_XZ,  VCF_GENOZIP    },\
                   { BCF,        CODEC_BCF,  BCF_GENOZIP    }, { BCF_GZ,    CODEC_GZ,  BCF_GENOZIP    }, { BCF_BGZF, CODEC_BGZF, BCF_GENOZIP  }, { } }, \
    [DT_SAM]   = { { SAM,        CODEC_NONE, SAM_GENOZIP    }, { SAM_GZ,    CODEC_GZ,  SAM_GENOZIP    }, { SAM_BGZF,  CODEC_GZ,  SAM_GENOZIP  },\
                   { SAM_BZ2,    CODEC_BZ2,  SAM_GENOZIP    }, { SAM_XZ,    CODEC_XZ,  SAM_GENOZIP    }, { }, },\
    [DT_BAM]   = { { BAM,        CODEC_BGZF, BAM_GENOZIP    }, { BAM_GZ,    CODEC_GZ,  BAM_GENOZIP    }, { BAM_BGZF, CODEC_BGZF, BAM_GENOZIP  }, \
                   { CRAM,       CODEC_CRAM, CRAM_GENOZIP   }, { } }, \
    [DT_FASTQ] = { { FASTQ,      CODEC_NONE, FASTQ_GENOZIP  }, { FASTQ_GZ,  CODEC_GZ,  FASTQ_GENOZIP  },\
                   { FASTQ_BZ2,  CODEC_BZ2,  FASTQ_GENOZIP  }, { FASTQ_XZ,  CODEC_XZ,  FASTQ_GENOZIP  }, { FASTQ_ORA, CODEC_ORA, FASTQ_GENOZIP  },\
                   { FQ,         CODEC_NONE, FQ_GENOZIP     }, { FQ_GZ,     CODEC_GZ,  FQ_GENOZIP     },\
                   { FQ_BZ2,     CODEC_BZ2,  FQ_GENOZIP     }, { FQ_XZ,     CODEC_XZ,  FQ_GENOZIP     }, { FQ_ORA,    CODEC_ORA, FASTQ_GENOZIP  }, { } },\
    [DT_FASTA] = { { FASTA,      CODEC_NONE, FASTA_GENOZIP  }, { FASTA_GZ,  CODEC_GZ,  FASTA_GENOZIP  },\
                   { FASTA_BZ2,  CODEC_BZ2,  FASTA_GENOZIP  }, { FASTA_XZ,  CODEC_XZ,  FASTA_GENOZIP  }, { FASTA_ZIP, CODEC_ZIP, FASTA_GENOZIP  }, \
                   { FAA,        CODEC_NONE, FAA_GENOZIP    }, { FAA_GZ,    CODEC_GZ,  FAA_GENOZIP    },\
                   { FAA_BZ2,    CODEC_BZ2,  FAA_GENOZIP    }, { FAA_XZ,    CODEC_XZ,  FAA_GENOZIP    }, { FAA_ZIP,   CODEC_ZIP, FAA_GENOZIP    }, \
                   { FFN,        CODEC_NONE, FFN_GENOZIP    }, { FFN_GZ,    CODEC_GZ,  FFN_GENOZIP    },\
                   { FFN_BZ2,    CODEC_BZ2,  FFN_GENOZIP    }, { FFN_XZ,    CODEC_XZ,  FFN_GENOZIP    }, { FFN_ZIP,   CODEC_ZIP, FFN_GENOZIP    }, \
                   { FNN,        CODEC_NONE, FNN_GENOZIP    }, { FNN_GZ,    CODEC_GZ,  FNN_GENOZIP    },\
                   { FNN_BZ2,    CODEC_BZ2,  FNN_GENOZIP    }, { FNN_XZ,    CODEC_XZ,  FNN_GENOZIP    }, { FNN_ZIP,   CODEC_ZIP, FNN_GENOZIP    }, \
                   { FNA,        CODEC_NONE, FNA_GENOZIP    }, { FNA_GZ,    CODEC_GZ,  FNA_GENOZIP    },\
                   { FNA_BZ2,    CODEC_BZ2,  FNA_GENOZIP    }, { FNA_XZ,    CODEC_XZ,  FNA_GENOZIP    }, { FNA_ZIP,   CODEC_ZIP, FNA_GENOZIP    }, \
                   { FRN,        CODEC_NONE, FRN_GENOZIP    }, { FRN_GZ,    CODEC_GZ,  FRN_GENOZIP    },\
                   { FRN_BZ2,    CODEC_BZ2,  FRN_GENOZIP    }, { FRN_XZ,    CODEC_XZ,  FRN_GENOZIP    }, { FRN_ZIP,   CODEC_ZIP, FRN_GENOZIP    }, \
                   { FAS,        CODEC_NONE, FAS_GENOZIP    }, { FAS_GZ,    CODEC_GZ,  FAS_GENOZIP    },\
                   { FAS_BZ2,    CODEC_BZ2,  FAS_GENOZIP    }, { FAS_XZ,    CODEC_XZ,  FAS_GENOZIP    }, { FAS_ZIP,   CODEC_ZIP, FAS_GENOZIP    }, \
                   { FSA,        CODEC_NONE, FSA_GENOZIP    }, { FSA_GZ,    CODEC_GZ,  FSA_GENOZIP    },\
                   { FSA_BZ2,    CODEC_BZ2,  FSA_GENOZIP    }, { FSA_XZ,    CODEC_XZ,  FSA_GENOZIP    }, { FSA_ZIP,   CODEC_ZIP, FSA_GENOZIP    }, \
                   { FA,         CODEC_NONE, FA_GENOZIP     }, { FA_GZ,     CODEC_GZ,  FA_GENOZIP     },\
                   { FA_BZ2,     CODEC_BZ2,  FA_GENOZIP     }, { FA_XZ,     CODEC_XZ,  FA_GENOZIP     }, { FA_ZIP,    CODEC_ZIP, FA_GENOZIP     }, { } },\
    [DT_GFF]   = { { GFF3,       CODEC_NONE, GFF3_GENOZIP   }, { GFF3_GZ,   CODEC_GZ,  GFF3_GENOZIP   },\
                   { GFF3_BZ2,   CODEC_BZ2,  GFF3_GENOZIP   }, { GFF3_XZ,   CODEC_XZ,  GFF3_GENOZIP   },\
                   { GFF,        CODEC_NONE, GFF_GENOZIP    }, { GFF_GZ,    CODEC_GZ,  GFF_GENOZIP    },\
                   { GFF_BZ2,    CODEC_BZ2,  GFF_GENOZIP    }, { GFF_XZ,    CODEC_XZ,  GFF_GENOZIP    },\
                   { GVF,        CODEC_NONE, GVF_GENOZIP    }, { GVF_GZ,    CODEC_GZ,  GVF_GENOZIP    },\
                   { GVF_BZ2,    CODEC_BZ2,  GVF_GENOZIP    }, { GVF_XZ,    CODEC_XZ,  GVF_GENOZIP    },\
                   { GTF,        CODEC_NONE, GTF_GENOZIP    }, { GTF_GZ,    CODEC_GZ,  GTF_GENOZIP    },\
                   { GTF_BZ2,    CODEC_BZ2,  GTF_GENOZIP    }, { GTF_XZ,    CODEC_XZ,  GTF_GENOZIP    }, { } },\
    [DT_ME23]  = { { ME23,       CODEC_NONE, ME23_GENOZIP   }, { ME23_ZIP,  CODEC_ZIP, ME23_GENOZIP   }, { } },\
    [DT_GNRIC] = { { GNRIC,      CODEC_NONE, GNRIC_GENOZIP  }, { GNRIC_GZ,  CODEC_GZ,  GNRIC_GENOZIP  },\
                   { GNRIC_BZ2,  CODEC_BZ2,  GNRIC_GENOZIP  }, { GNRIC_XZ,  CODEC_XZ,  GNRIC_GENOZIP  },\
                   { GNRIC_ZIP,  CODEC_ZIP,  GNRIC_GENOZIP  }, { } },\
    [DT_LOCS]  = { { LOCS,       CODEC_NONE, LOCS_GENOZIP   }, { LOCS_GZ,   CODEC_GZ,  LOCS_GENOZIP   },\
                   { LOCS_BZ2,   CODEC_BZ2,  LOCS_GENOZIP   }, { LOCS_XZ,   CODEC_XZ,  LOCS_GENOZIP   }, { } },\
    [DT_BED]   = { { BED,        CODEC_NONE, BED_GENOZIP    }, { BED_GZ,    CODEC_GZ,  BED_GENOZIP    },\
                   { BED_BZ2,    CODEC_BZ2,  BED_GENOZIP    }, { BED_XZ,    CODEC_XZ,  BED_GENOZIP    },\
                   { TRACK,      CODEC_NONE, TRACK_GENOZIP  }, { TRACK_GZ,  CODEC_GZ,  TRACK_GENOZIP  },\
                   { TRACK_BZ2,  CODEC_BZ2,  TRACK_GENOZIP  }, { TRACK_XZ,  CODEC_XZ,  TRACK_GENOZIP  }, { } },\
}

// Supported output formats for genounzip
// plain file MUST appear first on the list - this will be the default output when redirecting 
// GZ file, if it is supported MUST be 2nd on the list - we use this type if the user outputs to eg xx.gz instead of xx.vcf.gz (see file_open_txt_write)
#define TXT_OUT_FT_BY_DT {                                                    \
    [DT_VCF]   = { VCF, VCF_GZ, VCF_BGZF, 0 },                                \
    [DT_BCF]   = { BCF, BCF, 0 }, /* default types with or without bgzf */    \
    [DT_SAM]   = { SAM, SAM_GZ, 0 },                                          \
    [DT_BAM]   = { BAM, BAM, 0 }, /* default types with or without bgzf */    \
    [DT_CRAM]  = { CRAM, CRAM, 0 },                                           \
    [DT_FASTQ] = { FASTQ, FASTQ_GZ, FQ, FQ_GZ, 0 },                           \
    [DT_FASTA] = { FASTA, FASTA_GZ, FA, FA_GZ, FAA, FAA_GZ, FFN, FFN_GZ, FNN, FNN_GZ, FNA, FNA_GZ, FRN, FRN_GZ, FAS, FAS_GZ, FSA, FSA_GZ, 0 },\
    [DT_GFF]   = { GFF3, GFF3_GZ, GFF, GFF_GZ, GVF, GVF_GZ, GTF, GTF_GZ, 0 }, \
    [DT_ME23]  = { ME23, ME23 /* no GZ */, ME23_ZIP, 0 },                     \
    [DT_GNRIC] = { GNRIC, GNRIC_GZ, 0 },                                      \
    [DT_LOCS]  = { LOCS, LOCS_GZ, 0 },                                        \
    [DT_BED]   = { BED, BED_GZ, TRACK, TRACK_GZ, 0 },                         \
}                        

#define Z_FT_BY_DT {                                                          \
    [DT_REF]   = { REF_GENOZIP, 0  },                                         \
    [DT_VCF]   = { VCF_GENOZIP, BCF_GENOZIP, 0 },                             \
    [DT_SAM]   = { SAM_GENOZIP, BAM_GENOZIP, CRAM_GENOZIP, 0 },               \
    [DT_FASTQ] = { FASTQ_GENOZIP, FQ_GENOZIP, 0 },                            \
    [DT_FASTA] = { FASTA_GENOZIP, FA_GENOZIP, FAA_GENOZIP, FFN_GENOZIP, FNN_GENOZIP, FNA_GENOZIP, FRN_GENOZIP, FAS_GENOZIP, FSA_GENOZIP, 0 }, \
    [DT_GFF]   = { GFF3_GENOZIP, GFF_GENOZIP, GVF_GENOZIP, GTF_GENOZIP, 0  }, \
    [DT_ME23]  = { ME23_GENOZIP, 0 },                                         \
    [DT_GNRIC] = { GNRIC_GENOZIP, 0 },                                        \
    [DT_LOCS]  = { LOCS_GENOZIP, 0 },                                         \
    [DT_BED]   = { BED_GENOZIP, TRACK_GENOZIP, 0 },                           \
} 
