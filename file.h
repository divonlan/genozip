// ------------------------------------------------------------------
//   file.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef FILE_INCLUDED
#define FILE_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <pthread.h>
#else
#include "compatibility/visual_c_pthread.h"
#endif

#include "genozip.h"
#include "buffer.h"
#include "sections.h"
#include "move_to_front.h"
#include "aes.h"

// VCF file variations
#define VCF_           ".vcf"
#define VCF_GZ_        ".vcf.gz"
#define VCF_BGZ_       ".vcf.bgz"
#define VCF_BZ2_       ".vcf.bz2"
#define VCF_XZ_        ".vcf.xz"
#define BCF_           ".bcf"
#define BCF_GZ_        ".bcf.gz"
#define BCF_BGZ_       ".bcf.bgz"
#define VCF_GENOZIP_   ".vcf" GENOZIP_EXT

// SAM file variations
#define SAM_           ".sam"
#define SAM_GZ_        ".sam.gz"
#define SAM_BGZ_       ".sam.bgz"
#define SAM_BZ2_       ".sam.bz2"
#define SAM_XZ_        ".sam.xz"
#define BAM_           ".bam"
#define SAM_GENOZIP_   ".sam" GENOZIP_EXT

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
#define FAA_            ".faa"
#define FAA_GZ_         ".faa.gz"
#define FAA_BZ2_        ".faa.bz2"
#define FAA_XZ_         ".faa.xz"
#define FAA_GENOZIP_    ".faa" GENOZIP_EXT
#define FFN_            ".ffn"
#define FFN_GZ_         ".ffn.gz"
#define FFN_BZ2_        ".ffn.bz2"
#define FFN_XZ_         ".ffn.xz"
#define FFN_GENOZIP_    ".ffn" GENOZIP_EXT
#define FNN_            ".fnn"
#define FNN_GZ_         ".fnn.gz"
#define FNN_BZ2_        ".fnn.bz2"
#define FNN_XZ_         ".fnn.xz"
#define FNN_GENOZIP_    ".fnn" GENOZIP_EXT
#define FNA_            ".fna"
#define FNA_GZ_         ".fna.gz"
#define FNA_BZ2_        ".fna.bz2"
#define FNA_XZ_         ".fna.xz"
#define FNA_GENOZIP_    ".fna" GENOZIP_EXT

// GFF3 file variations (including GVF, which is a subtype of GFF3)
//#define GFF3_           ".gff3"
//#define GFF3_GZ_        ".gff3.gz"
//#define GFF3_BZ2_       ".gff3.bz2"
//#define GFF3_XZ_        ".gff3.xz"
//#define GFF3_GENOZIP_   ".gff3" GENOZIP_EXT
#define GVF_            ".gvf"
#define GVF_GZ_         ".gvf.gz"
#define GVF_BZ2_        ".gvf.bz2"
#define GVF_XZ_         ".gvf.xz"
#define GVF_GENOZIP_    ".gvf" GENOZIP_EXT

// 23andMe file variations
// note: 23andMe files come as a .txt, and therefore the user must specify --input to compress them. we have this
// made-up file extension here to avoid needing special cases throughout the code
#define ME23_           ".txt" // our made up extension - natively, 23andMe files come as a zip container containing a txt file
#define ME23_ZIP_       ".zip" 
#define ME23_GENOZIP_   ".txt" GENOZIP_EXT

typedef enum {TXT_FILE, Z_FILE} FileSupertype; 

typedef enum      { UNKNOWN_FILE_TYPE, 
                    VCF, VCF_GZ, VCF_BGZ, VCF_BZ2, VCF_XZ, BCF, BCF_GZ, BCF_BGZ, VCF_GENOZIP,  
                    SAM, SAM_GZ, SAM_BGZ, SAM_BZ2, SAM_XZ, BAM,                  SAM_GENOZIP,
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
                    AFTER_LAST_FILE_TYPE } FileType;

#define FILE_EXTS {"Unknown", /* order matches the FileType enum */ \
                   VCF_, VCF_GZ_, VCF_BGZ_, VCF_BZ2_, VCF_XZ_, BCF_, BCF_GZ_, BCF_BGZ_, VCF_GENOZIP_, \
                   SAM_, SAM_GZ_, SAM_BGZ_, SAM_BZ2_, SAM_XZ_, BAM_,                    SAM_GENOZIP_, \
                   FASTQ_, FASTQ_GZ_, FASTQ_BZ2_, FASTQ_XZ_, FASTQ_GENOZIP_, \
                   FQ_,    FQ_GZ_,    FQ_BZ2_,    FQ_XZ_,    FQ_GENOZIP_, \
                   FASTA_, FASTA_GZ_, FASTA_BZ2_, FASTA_XZ_, FASTA_GENOZIP_,\
                   FA_,    FA_GZ_,    FA_BZ2_,    FA_XZ_,    FA_GENOZIP_, \
                   FAA_,   FAA_GZ_,   FAA_BZ2_,   FAA_XZ_,   FAA_GENOZIP_, \
                   FFN_,   FFN_GZ_,   FFN_BZ2_,   FFN_XZ_,   FFN_GENOZIP_, \
                   FNN_,   FNN_GZ_,   FNN_BZ2_,   FNN_XZ_,   FNN_GENOZIP_, \
                   FNA_,   FNA_GZ_,   FNA_BZ2_,   FNA_XZ_,   FNA_GENOZIP_, \
                   /*GFF3_,  GFF3_GZ_,  GFF3_BZ2_,  GFF3_XZ_,  GFF3_GENOZIP_,*/ \
                   GVF_,   GVF_GZ_,   GVF_BZ2_,   GVF_XZ_,   GVF_GENOZIP_, \
                   ME23_,  ME23_ZIP_,                        ME23_GENOZIP_,\
                   "stdin", "stdout" }
extern const char *file_exts[];

// IMPORTANT: these values CANNOT BE CHANGED as they are part of the genozip file - 
// they go in SectionHeader.sec_compression_alg and also SectionHeaderTxtHeader.compression_type
#define NUM_COMPRESSION_ALGS 9
typedef enum { COMP_UNKNOWN=-1, COMP_PLN=0 /* plain - no compression */, 
               COMP_GZ=1, COMP_BZ2=2, COMP_BGZ=3, COMP_XZ=4, COMP_BCF=5, COMP_BAM=6, COMP_LZMA=7, COMP_ZIP=8 } CompressionAlg; 
#define COMPRESSED_FILE_VIEWER { "cat", "gunzip -d -c", "bzip2 -d -c", "gunzip -d -c", "xz -d -c", \
                                 "bcftools -Ov --version", "samtools view -h -OSAM", "N/A", "unzip -p" }

// txt file types and their corresponding genozip file types for each data type
// first entry of each data type MUST be the default plain file
#define TXT_IN_FT_BY_DT  { { { VCF,       COMP_PLN, VCF_GENOZIP   }, { VCF_GZ,   COMP_GZ,  VCF_GENOZIP   }, { VCF_BGZ, COMP_GZ,  VCF_GENOZIP },\
                             { VCF_BZ2,   COMP_BZ2, VCF_GENOZIP   }, { VCF_XZ,   COMP_XZ,  VCF_GENOZIP   },\
                             { BCF,       COMP_BCF, VCF_GENOZIP   }, { BCF_GZ,   COMP_BCF, VCF_GENOZIP   }, { BCF_BGZ, COMP_BCF, VCF_GENOZIP }, {0, 0, 0} },\
                           { { SAM,       COMP_PLN, SAM_GENOZIP   }, { SAM_GZ,   COMP_GZ,  SAM_GENOZIP   }, { SAM_BGZ, COMP_GZ,  SAM_GENOZIP },\
                             { SAM_BZ2,   COMP_BZ2, SAM_GENOZIP   }, { SAM_XZ,   COMP_XZ,  SAM_GENOZIP   },\
                             { BAM,       COMP_BAM, SAM_GENOZIP   },                                        { 0, 0, 0 }, },\
                           { { FASTQ,     COMP_PLN, FASTQ_GENOZIP }, { FASTQ_GZ, COMP_GZ,  FASTQ_GENOZIP },\
                             { FASTQ_BZ2, COMP_BZ2, FASTQ_GENOZIP }, { FASTQ_XZ, COMP_XZ,  FASTQ_GENOZIP },\
                             { FQ,        COMP_PLN, FQ_GENOZIP    }, { FQ_GZ,    COMP_GZ,  FQ_GENOZIP    },\
                             { FQ_BZ2,    COMP_BZ2, FQ_GENOZIP    }, { FQ_XZ,    COMP_XZ,  FQ_GENOZIP    }, { 0, 0, 0 } },\
                           { { FASTA,     COMP_PLN, FASTA_GENOZIP }, { FASTA_GZ, COMP_GZ,  FASTA_GENOZIP },\
                             { FASTA_BZ2, COMP_BZ2, FASTA_GENOZIP }, { FASTA_XZ, COMP_XZ,  FASTA_GENOZIP },\
                             { FAA,       COMP_PLN, FAA_GENOZIP   }, { FAA_GZ,   COMP_GZ,  FAA_GENOZIP   },\
                             { FAA_BZ2,   COMP_BZ2, FAA_GENOZIP   }, { FAA_XZ,   COMP_XZ,  FAA_GENOZIP   },\
                             { FFN,       COMP_PLN, FFN_GENOZIP   }, { FFN_GZ,   COMP_GZ,  FFN_GENOZIP   },\
                             { FFN_BZ2,   COMP_BZ2, FFN_GENOZIP   }, { FFN_XZ,   COMP_XZ,  FFN_GENOZIP   },\
                             { FNN,       COMP_PLN, FNN_GENOZIP   }, { FNN_GZ,   COMP_GZ,  FNN_GENOZIP   },\
                             { FNN_BZ2,   COMP_BZ2, FNN_GENOZIP   }, { FNN_XZ,   COMP_XZ,  FNN_GENOZIP   },\
                             { FNA,       COMP_PLN, FNA_GENOZIP   }, { FNA_GZ,   COMP_GZ,  FNA_GENOZIP   },\
                             { FNA_BZ2,   COMP_BZ2, FNA_GENOZIP   }, { FNA_XZ,   COMP_XZ,  FNA_GENOZIP   },\
                             { FA,        COMP_PLN, FA_GENOZIP    }, { FA_GZ,    COMP_GZ,  FA_GENOZIP    },\
                             { FA_BZ2,    COMP_BZ2, FA_GENOZIP    }, { FA_XZ,    COMP_XZ,  FA_GENOZIP    }, { 0, 0, 0 } },\
                           {/* { GFF3,      COMP_PLN, GFF3_GENOZIP  }, { GFF3_GZ,  COMP_GZ,  GFF3_GENOZIP  },\
                             { GFF3_BZ2,  COMP_BZ2, GFF3_GENOZIP  }, { GFF3_XZ,  COMP_XZ,  GFF3_GENOZIP  },*/ \
                             { GVF,       COMP_PLN, GVF_GENOZIP   }, { GVF_GZ,   COMP_GZ,  GVF_GENOZIP   },\
                             { GVF_BZ2,   COMP_BZ2, GVF_GENOZIP   }, { GVF_XZ,   COMP_XZ,  GVF_GENOZIP   }, { 0, 0, 0 } },\
                           { { ME23,      COMP_PLN, ME23_GENOZIP  }, { ME23_ZIP, COMP_ZIP, ME23_GENOZIP  }, { 0, 0, 0 } } }

// Supported output formats for genounzip
// plain file MUST appear first on the list - this will be the default output when redirecting
// GZ file, if it is supported MUST be 2nd on the list - we use this type if the user outputs to eg xx.gz instead of xx.vcf.gz
#define TXT_OUT_FT_BY_DT { { VCF,  VCF_GZ,  VCF_BGZ,  BCF, 0 },   \
                           { SAM,  SAM_GZ,  BAM, 0 },   \
                           { FASTQ, FASTQ_GZ, FQ, FQ_GZ, 0 },     \
                           { FASTA, FASTA_GZ, FA, FA_GZ, FAA, FAA_GZ, FFN, FFN_GZ, FNN, FNN_GZ, FNA, FNA_GZ, 0 },\
                           { GVF, GVF_GZ, /*GFF3, GFF3_GZ,*/ 0 }, \
                           { ME23, ME23_ZIP, 0 } }                        

// txt file types and their corresponding genozip file types for each data type
// first entry of each data type MUST be the default genozip file
#define Z_FT_BY_DT { { VCF_GENOZIP, 0  },               \
                     { SAM_GENOZIP, 0  },               \
                     { FASTQ_GENOZIP, FQ_GENOZIP, 0 },  \
                     { FASTA_GENOZIP, FA_GENOZIP, FAA_GENOZIP, FFN_GENOZIP, FNN_GENOZIP, FNA_GENOZIP, 0 }, \
                     { GVF_GENOZIP,/* GFF3_GENOZIP,*/ 0  },               \
                     { ME23_GENOZIP, 0 } } 

typedef const char *FileMode;
extern FileMode READ, WRITE; // this are pointers to static strings - so they can be compared eg "if (mode==READ)"

// ---------------------------
// tests for compression types
// ---------------------------

#define file_is_read_via_ext_decompressor(file) \
  (file->comp_alg == COMP_XZ || file->comp_alg == COMP_ZIP || file->comp_alg == COMP_BCF || file->comp_alg == COMP_BAM)

#define file_is_read_via_int_decompressor(file) \
  (file->comp_alg == COMP_GZ || file->comp_alg == COMP_BGZ || file->comp_alg == COMP_BZ2)

#define file_is_written_via_ext_compressor(file) (file->comp_alg == COMP_BCF || file->comp_alg == COMP_GZ)

#define file_is_plain_or_ext_decompressor(file) (file->comp_alg == COMP_PLN || file_is_read_via_ext_decompressor(file))

typedef struct File {
    void *file;
    char *name;                        // allocated by file_open(), freed by file_close()
    FileMode mode;
    FileSupertype supertype;            
    FileType type;
    bool is_remote;                    // true if file is downloaded from a url
    bool redirected;                   // true if this file is redirected from stdin/stdout
    DataType data_type;
    CompressionAlg comp_alg;            // txt_file: compression algorithm used with this file

    // these relate to actual bytes on the disk
    int64_t disk_size;                 // 0 if not known (eg stdin or http stream). 
                                       // note: this is different from txt_data_size_single as disk_size might be compressed (gz, bz2 etc)
    int64_t disk_so_far;               // data read/write to/from "disk" (using fread/fwrite)

    // this relate to the textual data represented. In case of READ - only data that was picked up from the read buffer.
    int64_t txt_data_size_single;      // txt_file: size of the txt data. ZIP: if its a plain txt file, then its the disk_size. If not, we initially do our best to estimate the size, and update it when it becomes known.
    int64_t txt_data_so_far_single;    // txt_file: data read (ZIP) or written (PIZ) to/from txt file so far
                                       // z_file: txt data represented in the GENOZIP data written (ZIP) or read (PIZ) to/from the genozip file so far for the current VCF
    int64_t txt_data_so_far_concat;    // z_file & ZIP only: txt data represented in the GENOZIP data written so far for all VCFs
    int64_t num_lines;                 // z_file: number of lines in all txt files concatenated into this z_file
                                       // txt_file: number of lines in single txt file

    // Used for READING & WRITING txt files - but stored in the z_file structure for zip to support concatenation (and in the txt_file structure for piz)
    Md5Context md5_ctx_concat;         // md5 context of txt file. in concat mode - of the resulting concatenated txt file
    Md5Context md5_ctx_single;         // used only in concat mode - md5 of the single txt component
    uint32_t max_lines_per_vb;         // ZIP & PIZ - in ZIP, discovered while segmenting, in PIZ - given by SectionHeaderTxtHeader

    // Used for READING GENOZIP files
    Buffer v1_next_vcf_header;         // genozip v1 only: next VCF header - used when reading in --split mode
    uint8_t genozip_version;           // GENOZIP_FILE_FORMAT_VERSION of the genozip file being read
    uint32_t num_components;           // set from genozip header

    // Used for WRITING GENOZIP files
    uint64_t disk_at_beginning_of_this_txt_file;     // z_file: the value of disk_size when starting to read this txt file
    
    SectionHeaderTxtHeader txt_header_first;         // store the first txt header - we might need to update it at the very end;
    uint8_t txt_header_enc_padding[AES_BLOCKLEN-1];  // just so we can overwrite txt_header with encryption padding

    SectionHeaderTxtHeader txt_header_single;        // store the txt header of single component in concat mode
    uint8_t txt_header_enc_padding2[AES_BLOCKLEN-1]; // same

    // dictionary information used for writing GENOZIP files - can be accessed only when holding mutex
    pthread_mutex_t dicts_mutex;
    bool dicts_mutex_initialized;      // this mutex protects contexts and num_dict_ids from concurrent adding of a new dictionary
    unsigned num_dict_ids;             // length of populated subfield_ids and mtx_ctx;
    
    uint8_t dict_id_to_did_i_map[65536]; // map for quick look up of did_i from dict_id 
    MtfContext contexts[MAX_DICTS];     // a merge of dictionaries of all VBs
    Buffer ra_buf;                     // RAEntry records - in a format ready to write to disk (Big Endian etc)
    Buffer dict_data;                  // Dictionary data accumulated from all VBs and written near the end of the file

    // section list - used for READING and WRITING genozip files
    Buffer section_list_buf;           // section list to be written as the payload of the genotype header section
    Buffer section_list_dict_buf;      // ZIP: a subset of section_list_buf - dictionaries are added here by VBs as they are being constructed
    uint32_t sl_cursor, sl_dir_cursor; // PIZ: next index into section_list for searching for sections
    uint32_t num_txt_components_so_far;

    // Information content stats - how many bytes and how many sections does this file have in each section type
    uint32_t num_vbs;

#   define READ_BUFFER_SIZE (1<<19)    // 512KB
    // Used for reading txt files
    Buffer unconsumed_txt;         // excess data read from the txt file - moved to the next VB
    
    // Used for reading genozip files
    uint32_t z_next_read, z_last_read;     // indices into read_buffer for z_file
    char read_buffer[];                  // only allocated for mode=READ files   
} File;

// globals
extern File *z_file, *txt_file; 

// methods
extern File *file_open (const char *filename, FileMode mode, FileSupertype supertype, DataType data_type /* only needed for WRITE */);
extern File *file_open_redirect (FileMode mode, FileSupertype supertype, DataType data_type);
extern bool file_open_txt (File *file);
extern void file_close (FileP *file_p, bool cleanup_memory /* optional */);
extern size_t file_write (FileP file, const void *data, unsigned len);
extern bool file_seek (File *file, int64_t offset, int whence, bool soft_fail); // SEEK_SET, SEEK_CUR or SEEK_END
extern uint64_t file_tell (File *file);
extern uint64_t file_get_size (const char *filename);
extern void file_set_input_type (const char *type_str);
extern void file_set_input_size (const char *size_str);
extern FileType file_get_stdin_type (void);
extern bool file_is_dir (const char *filename);
extern void file_remove (const char *filename, bool fail_quietly);
extern bool file_has_ext (const char *filename, const char *extension);
extern const char *file_basename (const char *filename, bool remove_exe, const char *default_basename,
                                  char *basename /* optional pre-allocated memory */, unsigned basename_size /* basename bytes */);
extern void file_get_file (VBlockP vb, const char *filename, Buffer *buf, const char *buf_name, unsigned buf_param, bool add_string_terminator);
extern void file_assert_ext_decompressor (void);
extern void file_kill_external_compressors (void);
extern FileType file_get_z_ft_by_txt_in_ft (DataType dt, FileType txt_ft);
extern const char *file_plain_ext_by_dt (DataType dt);

extern const char *ft_name (FileType ft);
extern const char *file_viewer (File *file);

#define file_printname(file) ((file)->name ? (file)->name : ((file)->mode==READ ? "(stdin)" : "(stdout)"))
#define txt_name file_printname(txt_file)
#define z_name   file_printname(z_file)

#define CLOSE(fd,name,quiet)  { ASSERTW (!close (fd) || (quiet),  "Warning in %s:%u: Failed to close %s: %s",  __FUNCTION__, __LINE__, (name), strerror(errno));}
#define FCLOSE(fp,name) { if (fp) { ASSERTW (!fclose (fp), "Warning in %s:%u: Failed to fclose %s: %s", __FUNCTION__, __LINE__, (name), strerror(errno)); fp = NULL; } }
 
// Windows compatibility stuff
#ifdef _WIN32
#define stat64  _stat64
#define fstat64 _fstat64
#else // this needs more work - there are more cases, depending if gcc is 32 or 64
#define stat64  stat
#define fstat64 fstat
#endif

#endif
