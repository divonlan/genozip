// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FILE_INCLUDED
#define FILE_INCLUDED

#include "genozip.h"
#include "context.h"
#include "aes.h"
#include "data_types.h"
#include "mutex.h"

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
#define CRAM_          ".cram"
#define SAM_GENOZIP_   ".sam" GENOZIP_EXT

// BAM
#define BAM_           ".bam" // can be bgzf-compressed or not
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

// GFF3 file variations (including GVF, which is a subtype of GFF3)
//#define GFF3_           ".gff3"
//#define GFF3_GZ_        ".gff3.gz"
//#define GFF3_BZ2_       ".gff3.bz2"
//#define GFF3_XZ_        ".gff3.xz"
//#define GFF3_GENOZIP_   ".gff3" GENOZIP_EXT
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

// generic files
#define GNRIC_           ""
#define GNRIC_GZ_        ".gz"
#define GNRIC_BZ2_       ".bz2"
#define GNRIC_XZ_        ".xz"
#define GNRIC_GENOZIP_   GENOZIP_EXT

typedef enum {TXT_FILE, Z_FILE} FileSupertype; 

typedef enum      { UNKNOWN_FILE_TYPE, 
                    REF_GENOZIP,
                    VCF, VCF_GZ, VCF_BGZF, VCF_BZ2, VCF_XZ, VCF_GENOZIP, 
                    SAM, SAM_GZ, SAM_BGZF, SAM_BZ2, SAM_XZ, CRAM, SAM_GENOZIP,
                    FASTQ, FASTQ_GZ, FASTQ_BZ2, FASTQ_XZ, FASTQ_GENOZIP,
                    FQ,    FQ_GZ,    FQ_BZ2,    FQ_XZ,    FQ_GENOZIP,
                    FASTA, FASTA_GZ, FASTA_BZ2, FASTA_XZ, FASTA_GENOZIP,
                    FA,    FA_GZ,    FA_BZ2,    FA_XZ,    FA_GENOZIP,
                    FAA,   FAA_GZ,   FAA_BZ2,   FAA_XZ,   FAA_GENOZIP,
                    FFN,   FFN_GZ,   FFN_BZ2,   FFN_XZ,   FFN_GENOZIP,
                    FNN,   FNN_GZ,   FNN_BZ2,   FNN_XZ,   FNN_GENOZIP,
                    FNA,   FNA_GZ,   FNA_BZ2,   FNA_XZ,   FNA_GENOZIP,
//                    GFF3,  GFF3_GZ,  GFF3_BZ2,  GFF3_XZ,  GFF3_GENOZIP,
                    GVF,   GVF_GZ,   GVF_BZ2,   GVF_XZ,   GVF_GENOZIP,
                    ME23,  ME23_ZIP,                      ME23_GENOZIP, 
                    PHY,   PHY_GZ,   PHY_BZ2,   PHY_XZ,   PHY_GENOZIP,
                    CHAIN, CHAIN_GZ, CHAIN_BZ2, CHAIN_XZ, CHAIN_GENOZIP,
                    BAM,                                  BAM_GENOZIP,
                    BCF, BCF_GZ, BCF_BGZF,                BCF_GENOZIP,  
                    // the GNRIC row *must* be the last row, as it consists catch-all extensions (*.gz etc)
                    GNRIC_GZ, GNRIC_BZ2, GNRIC_XZ,        GNRIC_GENOZIP, GNRIC, // GNRIC *must* be the very last as it is a catch-all ""
                    AFTER_LAST_FILE_TYPE } FileType;

#define FILE_EXTS {"Unknown", /* order matches the FileType enum */ \
                   REF_GENOZIP_, \
                   VCF_, VCF_GZ_, VCF_BGZF_, VCF_BZ2_, VCF_XZ_, VCF_GENOZIP_,        \
                   SAM_, SAM_GZ_, SAM_BGZF_, SAM_BZ2_, SAM_XZ_, CRAM_, SAM_GENOZIP_, \
                   FASTQ_, FASTQ_GZ_, FASTQ_BZ2_, FASTQ_XZ_, FASTQ_GENOZIP_,         \
                   FQ_,    FQ_GZ_,    FQ_BZ2_,    FQ_XZ_,    FQ_GENOZIP_,            \
                   FASTA_, FASTA_GZ_, FASTA_BZ2_, FASTA_XZ_, FASTA_GENOZIP_,         \
                   FA_,    FA_GZ_,    FA_BZ2_,    FA_XZ_,    FA_GENOZIP_,            \
                   FAA_,   FAA_GZ_,   FAA_BZ2_,   FAA_XZ_,   FAA_GENOZIP_,           \
                   FFN_,   FFN_GZ_,   FFN_BZ2_,   FFN_XZ_,   FFN_GENOZIP_,           \
                   FNN_,   FNN_GZ_,   FNN_BZ2_,   FNN_XZ_,   FNN_GENOZIP_,           \
                   FNA_,   FNA_GZ_,   FNA_BZ2_,   FNA_XZ_,   FNA_GENOZIP_,           \
                   /*GFF3_,  GFF3_GZ_,  GFF3_BZ2_,  GFF3_XZ_,  GFF3_GENOZIP_,*/      \
                   GVF_,   GVF_GZ_,   GVF_BZ2_,   GVF_XZ_,   GVF_GENOZIP_,           \
                   ME23_,  ME23_ZIP_,                        ME23_GENOZIP_,          \
                   PHY_,   PHY_GZ_,   PHY_BZ2_,   PHY_XZ_,   PHY_GENOZIP_,           \
                   CHAIN_, CHAIN_GZ_, CHAIN_BZ2_, CHAIN_XZ_, CHAIN_GENOZIP_,         \
                   BAM_,                                     BAM_GENOZIP_,           \
                   BCF_,   BCF_GZ_,   BCF_BGZF_,             BCF_GENOZIP_,           \
                   GNRIC_GZ_, GNRIC_BZ2_, GNRIC_XZ_,         GNRIC_GENOZIP_, GNRIC_, \
                   "stdin", "stdout" }
extern const char *file_exts[];

// Ordered by data_type: txt file types and their corresponding genozip file types for each data type
// first entry of each data type MUST be the default plain file
//                           { in         codec       out           }
#define TXT_IN_FT_BY_DT  { { { FASTA,     CODEC_NONE, REF_GENOZIP   }, { FASTA_GZ, CODEC_GZ,  REF_GENOZIP   },\
                             { FASTA_BZ2, CODEC_BZ2,  REF_GENOZIP   }, { FASTA_XZ, CODEC_XZ,  REF_GENOZIP   },\
                             { FA,        CODEC_NONE, REF_GENOZIP   }, { FA_GZ,    CODEC_GZ,  REF_GENOZIP   },\
                             { FA_BZ2,    CODEC_BZ2,  REF_GENOZIP   }, { FA_XZ,    CODEC_XZ,  REF_GENOZIP   }, { } }, \
                           { { VCF,       CODEC_NONE, VCF_GENOZIP   }, { VCF_GZ,   CODEC_GZ,  VCF_GENOZIP   }, { VCF_BGZF, CODEC_GZ,  VCF_GENOZIP },\
                             { VCF_BZ2,   CODEC_BZ2,  VCF_GENOZIP   }, { VCF_XZ,   CODEC_XZ,  VCF_GENOZIP   }, { } },\
                           { { SAM,       CODEC_NONE, SAM_GENOZIP   }, { SAM_GZ,   CODEC_GZ,  SAM_GENOZIP   }, { SAM_BGZF, CODEC_GZ,  SAM_GENOZIP },\
                             { SAM_BZ2,   CODEC_BZ2,  SAM_GENOZIP   }, { SAM_XZ,   CODEC_XZ,  SAM_GENOZIP   },\
                             { CRAM,      CODEC_CRAM, SAM_GENOZIP   }, { }, },\
                           { { FASTQ,     CODEC_NONE, FASTQ_GENOZIP }, { FASTQ_GZ, CODEC_GZ,  FASTQ_GENOZIP },\
                             { FASTQ_BZ2, CODEC_BZ2,  FASTQ_GENOZIP }, { FASTQ_XZ, CODEC_XZ,  FASTQ_GENOZIP },\
                             { FQ,        CODEC_NONE, FQ_GENOZIP    }, { FQ_GZ,    CODEC_GZ,  FQ_GENOZIP    },\
                             { FQ_BZ2,    CODEC_BZ2,  FQ_GENOZIP    }, { FQ_XZ,    CODEC_XZ,  FQ_GENOZIP    }, { } },\
                           { { FASTA,     CODEC_NONE, FASTA_GENOZIP }, { FASTA_GZ, CODEC_GZ,  FASTA_GENOZIP },\
                             { FASTA_BZ2, CODEC_BZ2,  FASTA_GENOZIP }, { FASTA_XZ, CODEC_XZ,  FASTA_GENOZIP },\
                             { FAA,       CODEC_NONE, FAA_GENOZIP   }, { FAA_GZ,   CODEC_GZ,  FAA_GENOZIP   },\
                             { FAA_BZ2,   CODEC_BZ2,  FAA_GENOZIP   }, { FAA_XZ,   CODEC_XZ,  FAA_GENOZIP   },\
                             { FFN,       CODEC_NONE, FFN_GENOZIP   }, { FFN_GZ,   CODEC_GZ,  FFN_GENOZIP   },\
                             { FFN_BZ2,   CODEC_BZ2,  FFN_GENOZIP   }, { FFN_XZ,   CODEC_XZ,  FFN_GENOZIP   },\
                             { FNN,       CODEC_NONE, FNN_GENOZIP   }, { FNN_GZ,   CODEC_GZ,  FNN_GENOZIP   },\
                             { FNN_BZ2,   CODEC_BZ2,  FNN_GENOZIP   }, { FNN_XZ,   CODEC_XZ,  FNN_GENOZIP   },\
                             { FNA,       CODEC_NONE, FNA_GENOZIP   }, { FNA_GZ,   CODEC_GZ,  FNA_GENOZIP   },\
                             { FNA_BZ2,   CODEC_BZ2,  FNA_GENOZIP   }, { FNA_XZ,   CODEC_XZ,  FNA_GENOZIP   },\
                             { FA,        CODEC_NONE, FA_GENOZIP    }, { FA_GZ,    CODEC_GZ,  FA_GENOZIP    },\
                             { FA_BZ2,    CODEC_BZ2,  FA_GENOZIP    }, { FA_XZ,    CODEC_XZ,  FA_GENOZIP    }, { } },\
                           {/* { GFF3,      CODEC_NONE, GFF3_GENOZIP  }, { GFF3_GZ,  CODEC_GZ,  GFF3_GENOZIP  },\
                             { GFF3_BZ2,  CODEC_BZ2,  GFF3_GENOZIP  }, { GFF3_XZ,  CODEC_XZ,  GFF3_GENOZIP  },*/ \
                             { GVF,       CODEC_NONE, GVF_GENOZIP   }, { GVF_GZ,   CODEC_GZ,  GVF_GENOZIP   },\
                             { GVF_BZ2,   CODEC_BZ2,  GVF_GENOZIP   }, { GVF_XZ,   CODEC_XZ,  GVF_GENOZIP   }, { } },\
                           { { ME23,      CODEC_NONE, ME23_GENOZIP  }, { ME23_ZIP, CODEC_ZIP, ME23_GENOZIP  }, { } },\
                           { { BAM,       CODEC_BGZF, BAM_GENOZIP   }, { } }, \
                           { { BCF,       CODEC_BCF,  BCF_GENOZIP   }, { BCF_GZ,   CODEC_BCF, BCF_GENOZIP   }, { BCF_BGZF, CODEC_BCF, BCF_GENOZIP }, { } }, \
                           { { GNRIC,     CODEC_NONE, GNRIC_GENOZIP }, { GNRIC_GZ, CODEC_GZ,  GNRIC_GENOZIP },\
                             { GNRIC_BZ2, CODEC_BZ2,  GNRIC_GENOZIP }, { GNRIC_XZ, CODEC_XZ,  GNRIC_GENOZIP }, { } },\
                           { { PHY,       CODEC_NONE, PHY_GENOZIP   }, { PHY_GZ,   CODEC_GZ,  PHY_GENOZIP   },\
                             { PHY_BZ2,   CODEC_BZ2,  PHY_GENOZIP   }, { PHY_XZ,   CODEC_XZ,  PHY_GENOZIP   }, { } },\
                           { { CHAIN,     CODEC_NONE, CHAIN_GENOZIP }, { CHAIN_GZ, CODEC_GZ,  CHAIN_GENOZIP },\
                             { CHAIN_BZ2, CODEC_BZ2,  CHAIN_GENOZIP }, { CHAIN_XZ, CODEC_XZ,  CHAIN_GENOZIP }, { } },\
                        }

// Ordered by data_type: Supported output formats for genounzip
// plain file MUST appear first on the list - this will be the default output when redirecting 
// GZ file, if it is supported MUST be 2nd on the list - we use this type if the user outputs to eg xx.gz instead of xx.vcf.gz (see file_open_txt_write)
#define TXT_OUT_FT_BY_DT { { 0 }, /* a reference file cannot be uncompressed */  \
                           { VCF, VCF_GZ, VCF_BGZF, BCF, 0 },  \
                           { SAM, SAM_GZ, BAM, 0 },     \
                           { FASTQ, FASTQ_GZ, FQ, FQ_GZ, 0 }, \
                           { FASTA, FASTA_GZ, FA, FA_GZ, FAA, FAA_GZ, FFN, FFN_GZ, FNN, FNN_GZ, FNA, FNA_GZ, 0 },\
                           { GVF, GVF_GZ, /*GFF3, GFF3_GZ,*/ 0 }, \
                           { ME23, ME23 /* no GZ */, ME23_ZIP, 0 }, \
                           { 0 }, /* There are no data_type=DT_BAM genozip files - .bam.genozip have data_type=DT_SAM */ \
                           { 0 }, /* There are no data_type=DT_BCF genozip files - .bam.genozip have data_type=DT_VCF */ \
                           { GNRIC, GNRIC_GZ, 0 }, \
                           { PHY, PHY_GZ, 0 }, \
                           { CHAIN, CHAIN_GZ, 0 }, \
                         }                        

// Ordered by data_type
#define Z_FT_BY_DT { { REF_GENOZIP, 0  },                   \
                     { VCF_GENOZIP, 0  },                   \
                     { SAM_GENOZIP, BAM_GENOZIP, 0 },       \
                     { FASTQ_GENOZIP, FQ_GENOZIP, 0 },      \
                     { FASTA_GENOZIP, FA_GENOZIP, FAA_GENOZIP, FFN_GENOZIP, FNN_GENOZIP, FNA_GENOZIP, 0 }, \
                     { GVF_GENOZIP,/* GFF3_GENOZIP,*/ 0  }, \
                     { ME23_GENOZIP, 0 },                   \
                     { 0 }, /* There are no data_type=DT_BAM genozip files - .bam.genozip have data_type=DT_SAM */ \
                     { 0 }, /* There are no data_type=DT_BCF genozip files - .bam.genozip have data_type=DT_VCF */ \
                     { GNRIC_GENOZIP, 0 },                  \
                     { PHY_GENOZIP, 0 },                    \
                     { CHAIN_GENOZIP, 0 },                  \
                   } 

typedef const char *FileMode;
extern FileMode READ, WRITE, WRITEREAD; // this are pointers to static strings - so they can be compared eg "if (mode==READ)"

typedef enum { DC_NONE, DC_PRIMARY, DC_LUFT } DualCoordinates;

// ---------------------------
// tests for compression types
// ---------------------------

#define file_is_read_via_ext_decompressor(file) (file->supertype == TXT_FILE && \
  (file->codec == CODEC_XZ || file->codec == CODEC_ZIP || file->codec == CODEC_BCF || file->codec == CODEC_CRAM))

#define file_is_read_via_int_decompressor(file) (file->supertype == TXT_FILE && \
  (file->codec == CODEC_GZ || file->codec == CODEC_BGZF || file->codec == CODEC_BZ2))

#define file_is_written_via_ext_compressor(file) (file->supertype == TXT_FILE && \
  (file->codec == CODEC_BCF || file->codec == CODEC_GZ))

#define file_is_plain_or_ext_decompressor(file) (file->supertype == TXT_FILE && \
  (file->codec == CODEC_NONE || file_is_read_via_ext_decompressor(file)))

typedef struct File {
    void *file;
    char *name;                        // allocated by file_open(), freed by file_close()
    const char *basename;              // basename of name
    FileMode mode;
    FileSupertype supertype;            
    FileType type;
    bool is_remote;                    // true if file is downloaded from a url
    bool redirected;                   // txt_file: true if this file is redirected from stdin/stdout
    bool is_eof;                       // we've read the entire file
    DataType data_type;
    Codec codec;                       // ZIP - txt_file: generic codec used with this file (in PIZ we use flag.bgzf instead)
                                       // ZIP - z_file: copy from txt_file.codec, but not for rejects file    
    // these relate to actual bytes on the disk
    int64_t disk_size;                 // 0 if not known (eg stdin or http stream). 
                                       // note: this is different from txt_data_size_single as disk_size might be compressed (gz, bz2 etc)
    int64_t disk_so_far;               // data read/write to/from "disk" (using fread/fwrite)

    // this relate to the textual data represented. In case of READ - only data that was picked up from the read buffer.
    int64_t txt_data_size_single;      // txt_file: size of the txt data. ZIP: if its a plain txt file, then its the disk_size. If not, we initially do our best to estimate the size, and update it when it becomes known.
    int64_t txt_data_so_far_single;    // txt_file: data read (ZIP) or written (PIZ) to/from txt file so far
                                       // z_file: txt data represented in the GENOZIP data written (ZIP) or read (PIZ) to/from the genozip file so far for the current VCF
    int64_t txt_data_so_far_bind;      // z_file & ZIP only: uncompressed txt data represented in the GENOZIP data written so far for all bound files
                                       // note: txt_data_so_far_single/bind accounts for txt modifications due to --optimize or --chain or compressing a Luft file
    int64_t txt_data_so_far_single_0;  // z_file & ZIP only: same as txt_data_so_far_single/bind, but original sizes without 
    int64_t txt_data_so_far_bind_0;    //      modifications due to --chain/--optimize/Luft                                   
    int64_t txt_disk_so_far_bind;      // z_file & ZIP only: compressed (with txt_file.codec - eg bgzf) txt data represented in the GENOZIP data written so far for all bound files
    int64_t num_lines;                 // z_file: number of lines in all txt files bound into this z_file
                                       // txt_file: number of lines in single txt file

    // Used for READING & WRITING txt files - but stored in the z_file structure for zip to support bindenation (and in the txt_file structure for piz)
    DigestContext digest_ctx_bound;    // md5 context of txt file. in bound mode - of the resulting bound txt file
    DigestContext digest_ctx_single;   // used only in bound mode - md5 of the single txt component
    Digest digest;                     // PIZ: as read from header: z_file: digest_bound txt_file: digest_single

    uint32_t max_lines_per_vb;         // ZIP & PIZ - in ZIP, discovered while segmenting, in PIZ - given by SectionHeaderTxtHeader

    // Used for READING GENOZIP files
    uint8_t genozip_version;           // GENOZIP_FILE_FORMAT_VERSION of the genozip file being read

    struct FlagsGenozipHeader z_flags; // genozip file flags as read from SectionHeaderGenozipHeader.h.flags
    uint32_t num_components;           // set from genozip header
    uint32_t num_copied_ref_sections;  // number of SEC_REFERENCE sections with vblock_i==0 (a result of being created with ref_copy_one_compressed_section)
     
    // Used for WRITING GENOZIP files
    uint64_t disk_at_beginning_of_this_txt_file;     // z_file: the value of disk_size when starting to read this txt file
    
    SectionHeaderTxtHeader txt_header_first;         // store the first txt header - we might need to update it at the very end;
    uint8_t txt_header_enc_padding[AES_BLOCKLEN-1];  // just so we can overwrite txt_header with encryption padding

    SectionHeaderTxtHeader txt_header_single;        // store the txt header of single component in bound mode
    uint8_t txt_header_enc_padding2[AES_BLOCKLEN-1]; // same

    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    Mutex dicts_mutex;                 // this mutex protects contexts and num_contexts from concurrent adding of a new dictionary
    DidIType num_contexts;             // length of populated subfield_ids and mtx_ctx;
    
    DidIType dict_id_to_did_i_map[65536]; // map for quick look up of did_i from dict_id 
    Context contexts[MAX_DICTS];       // a merge of dictionaries of all VBs
    Buffer ra_buf;                     // RAEntry records - in a format ready to write to disk (Big Endian etc)
    Buffer ra_min_max_by_chrom;        // max_pos of each from according to RA. An array of uint32_t indexed by chrom_word_index.
    Buffer chroms_sorted_index;        // PIZ: index into contexts[CHROM]->word_list, sorted alphabetically by snip
    
    Buffer alt_chrom_map;              // ZIP/PIZ: mapping from user file chrom to alternate chrom in reference file

    // section list - used for READING and WRITING genozip files
    Buffer section_list_buf;           // section list to be written as the payload of the genotype header section
    uint32_t num_txt_components_so_far;// ZIP/PIZ z_file

    // TXT file: stuff reading and writing txt files compressed with BGZF
    Buffer bgzf_isizes;                // of the bgzf blocks in which this txt file is compressed (in BGEN16)
    struct FlagsBgzf bgzf_flags;       // correspond to SectionHeader.flags in SEC_BGZF
    uint8_t bgzf_signature[3];         // PIZ: 3 LSB of size of source BGZF-compressed file, as passed in SectionHeaderTxtHeader.codec_info
    int32_t bzgf_passed_down_len;      // PIZ: bytes at the end of the VB too small for one bgzf block passed to the next block

    // TXT file: data used in --show-sex, --show-coverage and --idxstats
    Buffer coverage;
    Buffer read_count;
    Buffer unmapped_read_count;

    // Z_FILE: Liftover stuff
    char *rejects_file_name;            // ZIP: allocated and freed by liftover_*
    FILE *rejects_file;                 // ZIP: rejects txt file
    uint64_t rejects_disk_size;         // ZIP

    // TXT_FILE: Liftover stuff
    DualCoordinates dual_coords;       // ZIP: dual coordinate status of the TXT file. Set when reading TXT header, and immutable thereafter
    uint64_t luft_reject_bytes;        // ZIP of Luft dual coordinate file: number of bytes of that are rejected lines, not yet assigned to a VB

    // Reconstruction plan, for reconstructing in sorted order if --sort: [0] is primary coords, [1] is luft coords
    Mutex recon_plan_mutex[2];         // TXT_FILE ZIP: protect vb_info and line_info during merging of VB data
    Buffer vb_info[2];                 // TXT_FILE ZIP: array of ZipVbInfo per VB, indexed by (vb_i-1), 0:PRIMARY, 1:LUFT
                                       // Z_FILE   PIZ: array of PizVbInfo per VB, indexed by (vb_i-1), only vb_info[0] is used
    Buffer line_info[2];               // TXT_FILE ZIP: array of LineInfo per line or gapless range in txt_file
    Buffer recon_plan;                 // TXT_FILE ZIP/PIZ: array of ReconPlanItem - order of reconstruction of ranges of lines, to achieve a sorted file
                                       // Z_FILE   PIZ: plan for entire z_file, txt_file.recon_plan is assigned a portion of this plan
    Buffer comp_info;                  // Z_FILE   PIZ: array of PizCompInfo - component information
    Buffer txt_file_info;              // Z_FILE   PIZ: array of PizTxtFileInfo - txt_file information
    uint32_t lines_so_far;             // TXT_FILE PIZ: number of textual lines (excluding the header) that passed all filters except downsampling, and is to be written to txt_file, or downsampled-out

    // Z_FILE: stats data
    Buffer stats_buf, STATS_buf;       // Strings to be outputted in case of --stats or --STATS (generated during ZIP, stored in SEC_STATS)
    Buffer bound_txt_names;            // ZIP: a concatenation of all bound txt_names that contributed to this genozip file
    bool is_txt_len_frozen;            // ZIP: true if frozen_txt_len is assigned
    uint64_t frozen_txt_len[MAX_DICTS];// ZIP: copy of contexts[].txt_len, from stats_freeze_txt_len

    // Information content stats - how many bytes and how many sections does this file have in each section type
    uint32_t num_vbs;

    // Used for reading txt files
    Buffer unconsumed_txt;             // excess data read from the txt file - moved to the next VB
} File;

// methods
extern File *file_open (const char *filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE */);
extern void file_close (FileP *file_p, bool index_txt, bool cleanup_memory /* optional */);
extern void file_write (FileP file, const void *data, unsigned len);
extern bool file_seek (File *file, int64_t offset, int whence, int soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
extern uint64_t file_tell (File *file);
extern void file_set_input_type (const char *type_str);
extern void file_set_input_size (const char *size_str);
extern FileType file_get_type (const char *filename);
extern FileType file_get_stdin_type (void);
extern DataType file_get_data_type (FileType ft, bool is_input);
extern const char *file_plain_text_ext_of_dt (DataType dt);
extern void file_get_raw_name_and_type (char *filename, char **raw_name, FileType *ft);
extern bool file_has_ext (const char *filename, const char *extension);
extern const char *file_basename (const char *filename, bool remove_exe, const char *default_basename,
                                  char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */);
extern void file_assert_ext_decompressor (void);
extern void file_kill_external_compressors (void);
extern FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft);
extern DataType file_get_dt_by_z_ft (FileType z_ft);
extern FileType file_get_z_ft_by_dt (DataType dt);
extern const char *file_plain_ext_by_dt (DataType dt);
extern void file_remove_codec_ext (char *filename, FileType ft);
extern const char *ft_name (FileType ft);
extern const char *file_guess_original_filename (const File *file);
extern char *file_get_fastq_pair_filename (const char *fn1, const char *fn2, bool test_only);

// wrapper operations for operating system files
extern void file_get_file (VBlockP vb, const char *filename, Buffer *buf, const char *buf_name, bool add_string_terminator);
extern bool file_put_data (const char *filename, void *data, uint64_t len);
extern bool file_exists (const char *filename);
extern bool file_is_fifo (const char *filename);
extern bool file_is_dir (const char *filename);
extern void file_remove (const char *filename, bool fail_quietly);
extern void file_mkfifo (const char *filename);
extern uint64_t file_get_size (const char *filename);
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

#define CLOSE(fd,name,quiet) do { ASSERTW (!close (fd) || (quiet),  "Warning in %s:%u: Failed to close %s: %s",  __FUNCTION__, __LINE__, (name), strerror(errno));} while (0)
#define FCLOSE(fp,name) do { if (fp) { ASSERTW (!fclose (fp), "Warning in %s:%u: Failed to fclose %s: %s", __FUNCTION__, __LINE__, (name), strerror(errno)); fp = NULL; } } while (0)
 
// Windows compatibility stuff
#ifdef _WIN32
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#endif
