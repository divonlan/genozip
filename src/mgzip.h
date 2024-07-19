// ------------------------------------------------------------------
//   mgzip.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
 
#include "sections.h"

#define TXT_IS_PLAIN (txt_file->effective_codec == CODEC_NONE)
#define TXT_IS_BGZF  (txt_file->effective_codec == CODEC_BGZF)
#define TXT_IS_IL1M  (txt_file->effective_codec == CODEC_IL1M)
#define TXT_IS_MGZF  (txt_file->effective_codec == CODEC_MGZF)
#define TXT_IS_MGSP  (txt_file->effective_codec == CODEC_MGSP)
#define TXT_IS_EMFL  (txt_file->effective_codec == CODEC_EMFL)
#define TXT_IS_EMVL  (txt_file->effective_codec == CODEC_EMVL)
#define TXT_IS_GZ    (txt_file->effective_codec == CODEC_GZ)
#define TXT_IS_BZ2   (txt_file->effective_codec == CODEC_BZ2)

#define IS_BGZF(codec)  ((codec)==CODEC_BGZF)
#define IS_MGZF(codec)  ((codec)==CODEC_MGZF)
#define IS_MGSP(codec)  ((codec)==CODEC_MGSP)
#define IS_IL1M(codec)  ((codec)==CODEC_IL1M)
#define IS_EMFL(codec)  ((codec)==CODEC_EMFL)
#define IS_EMVL(codec)  ((codec)==CODEC_EMVL)
#define IS_GZ(codec)    ((codec)==CODEC_GZ)
#define IS_BZ2(codec)   ((codec)==CODEC_BZ2)
#define IS_NONE(codec)  ((codec)==CODEC_NONE)
#define IS_MGZIP(codec) (IS_BGZF(codec) || IS_MGZF(codec) || IS_MGSP(codec) || IS_IL1M(codec) || IS_EMFL(codec) || IS_EMVL(codec)) // multi-block gzip
#define IS_GZIP(codec)  (IS_MGZIP(codec) || IS_GZ(codec))

// note on MGSP: "gz block" in the comments below means, for MGSP, a group of gz blocks.
#define IS_IN_SYNC(codec)          (IS_MGZF(codec) || IS_MGSP(codec) || IS_EMVL(codec) || IS_EMFL(codec)) // codecs in which R1 and R2 gz blocks are guaranteed to contain whole, and precisely matching reads. Therefore, R2 gz-decompression can delegated to compute threads without further checks.
#define IS_VB_SIZE_BY_BLOCK(codec) (IS_MGZF(codec) || IS_EMVL(codec))             // codecs that are 1. variable-length 2. reads are never split between blocks 3. we use on VB per gz block
#define IS_VB_SIZE_BY_MGZIP(codec) (IS_VB_SIZE_BY_BLOCK(codec) || IS_MGSP(codec)) // like IS_VB_SIZE_BY_BLOCK, but VB can be a group of gz blocks
#define GZ_HEADER_HAS_BSIZE(codec) (IS_BGZF(codec) || IS_MGZF(codec)) // gz header contains bsize
#define IS_EXACTABLE(codec)     (IS_BGZF(codec))    // codecs for which we can we discover the library level and can reconstruct exactly  

#define TXT_IS_MGZIP            IS_MGZIP(txt_file->effective_codec)
#define TXT_IS_GZIP             IS_GZIP (txt_file->effective_codec)
#define TXT_IS_VB_SIZE_BY_BLOCK IS_VB_SIZE_BY_BLOCK(txt_file->effective_codec) 
#define TXT_IS_VB_SIZE_BY_MGZIP IS_VB_SIZE_BY_MGZIP(txt_file->effective_codec) 
#define TXT_IS_IN_SYNC          IS_IN_SYNC(txt_file->effective_codec)
#define TXT_GZ_HEADER_HAS_BSIZE GZ_HEADER_HAS_BSIZE(txt_file->effective_codec)

#define BGZF_DEFAULT_LEVEL 2 // PIZ: used if --bgzf is not specified (it is actually faster than 1 if also writing to disk)
#define BGZF_MAX_BLOCK_SIZE ((uint32_t)(64 KB)) // maximum block size of both compressed and uncompressed data of one block
#define BGZF_MAX_CHUCK_SIZE ((uint32_t)(1 MB))  // max amount we read from disk at a time  
#define BGZF_PREFIX "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00" // First 16 bytes of every BGZF block
#define BGZF_HEADER_LEN 18
#define BGZF_EOF BGZF_PREFIX "\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00" // // BGZF EOF marker is simply an empty block (note: there are multiple encoding for empty blocks, this is a specific one of them), see https://samtools.github.io/hts-specs/SAMv1.pdf section 4.1.2
#define BGZF_EOF_LEN 28

#define IL1M_HEADER "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03"
#define IL1M_ISIZE  "\x00\x00\x10\x00" // isize == 1MB in all blocks except the last

// MGI: a 32-bit version of BGZF 
#define MGZF_PREFIX_LEN 16
#define MGZF_PREFIX "\x1f\x8b\x08\x14\x00\x00\x00\x00\x00\xff\x08\x00\x49\x47\x04\x00"
#define MGZF_EOF_LEN 31
#define MGZF_EOF MGZF_PREFIX "\x1f\x00\x00\x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
#define MGZF_HEADER_LEN 29

// MGI: constant isize for all gz blocks that go into a particular VB (last block in group might slightly bigger)
#define MGSP_HEADER     IL1M_HEADER
#define MGSP_EOF_LEN 20
#define MGSP_EOF MGSP_HEADER "\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

// Element: constant isize of all gz blocks in file (last block might be smaller)
#define EMFL_HEADER_LEN 10

#define EMVL_HEADER "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"
#define EMVL_HEADER_LEN 10
#define EMVL_FIRST_BLOCK EMVL_HEADER "\x01\x00\x00\xff\xff\x00\x00\x00\x00\x00\x00\x00\x00" // EMVL files begin with this empty block and have no EOF block

// fixed length of gz header of each MGZIP codec
#define MGZIP_HEADER_LEN_BY_CODEC {     \
    [CODEC_BGZF] = BGZF_HEADER_LEN,     \
    [CODEC_MGZF] = MGZF_HEADER_LEN,     \
    [CODEC_EMFL] = EMFL_HEADER_LEN,     \
    [CODEC_IL1M] = STRLEN(IL1M_HEADER), \
    [CODEC_MGSP] = STRLEN(MGSP_HEADER), \
    [CODEC_EMVL] = STRLEN(EMVL_HEADER), \
}

// for capped-isize codecs: VB will be vb_size, and full or partial gz blocks
#define MAX_ISIZE_BY_CODEC {            \
    [CODEC_BGZF] = 64 KB,               \
    [CODEC_IL1M] = 1 MB,                \
    [CODEC_EMFL] = txt_file->max_mgzip_isize /* determined during discovery */ \
}    

typedef bool (*IsValidSize)(FileP, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb);

// "no_bsize" codecs
typedef struct { 
    IsValidSize is_valid_isize; // isize validation function
    bool valid_3_blocks_isize;  // the first 3 gz blocks are expected to have the same isize (used for discovery) 
    uint32_t max_bsize;         // an upper limit we set on compressed size (bsize) of a block based on observation (after discovery)
    bytes gz_hdr;               // fixed gz header
    uint32_t gz_hdr_len;
} NoBsizeCodecParams;

#define NO_BSIZE_CODECS_PARAMS {                                                                   \
    [CODEC_MGSP] = { mgsp_is_valid_isize, true,  4  MB, (bytes)MGSP_HEADER, STRLEN(MGSP_HEADER) }, \
    [CODEC_EMVL] = { emvl_is_valid_isize, false, 32 MB, (bytes)EMVL_HEADER, STRLEN(EMVL_HEADER) }, \
    [CODEC_IL1M] = { il1m_is_valid_isize, true,  1  MB, (bytes)IL1M_HEADER, STRLEN(IL1M_HEADER) }, \
    [CODEC_EMFL] = { emfl_is_valid_isize, true,  4  MB, NULL/*run time*/,   EMFL_HEADER_LEN     }, \
}    

typedef struct BgzfBlockPiz {
    int32_t txt_index, txt_size; // index of uncompressed block within vb->txt_data. The first block index will be negative if there is passed-down unconsumed data
} BgzfBlockPiz;

// ZIP side
typedef enum           { GZ_SUCCESS, GZ_IS_OTHER_FORMAT,  GZ_MORE_DATA, GZ_NOT_GZIP, GZ_EOF_WITHOUT_EOF_BLOCK, GZ_TRUNCATED, NUM_GZ_STATUSES } GzStatus; // file is truncated
#define GZSTATUS_NAMES {   "SUCCESS",  "IS_OTHER_FORMAT",   "MORE_DATA",  "NOT_GZIP"   "EOF_WITHOUT_EOF_BLOCK",  "TRUNCATED",                }

// data type of VBlock.gz_blocks and txt_file->unconsumed_mgzip_blocks : details of MGZIP blocks.
typedef struct GzBlockZip {
    int32_t txt_index;                    // index of uncompressed block within vb->txt_data. If there is passed-down data from previous VB/txt_header, then txt_index of the first block will be negative (see mgzip_copy_unconsumed_blocks)
    uint32_t txt_size        : 30; 
    uint32_t is_uncompressed : 1;         // true if data has been GZ-decompressed by main thread
    uint32_t is_eof          : 1;         // true if this is the last GZ-block in the file
    uint32_t compressed_index, comp_size; // index within vb->scratch
} GzBlockZip;

extern GzStatus mgzip_read_block_with_bsize (FileP file, bool discovering, Codec codec);
extern GzStatus mgzip_read_block_no_bsize (FileP file, bool discovering, Codec codec);
extern void mgzip_uncompress_vb (VBlockP vb, Codec codec);
extern void mgzip_uncompress_one_block (VBlockP vb, GzBlockZip *bb, Codec codec);
extern void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *file);
extern void mgzip_compress_mgzip_section (void);
extern void mgzip_zip_advance_index (VBlockP vb, uint32_t line_len);
extern int64_t mgzip_copy_unconsumed_blocks (VBlockP vb);
extern void mgzip_zip_init_vb (VBlockP vb);
extern void bgzf_insert_back_segconf_blocks (VBlockP vb);
extern void mgzip_return_segconf_blocks (VBlockP vb);
extern uint32_t mgzip_get_max_block_size (void);

extern void inc_disk_gz_uncomp_or_trunc_(FileP file, uint64_t inc, FUNCLINE);
#define inc_disk_gz_uncomp_or_trunc(file, inc) inc_disk_gz_uncomp_or_trunc_((file), (inc), __FUNCLINE)

// codec size validators
extern bool il1m_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb);
extern bool mgsp_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb);
extern bool emfl_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb);
extern bool emvl_is_valid_isize (FileP file, uint32_t proposed_isize, bool is_eof, bool *is_end_of_vb);

// library / level discovery
extern void bgzf_initialize_discovery (FileP file);
extern void bgzf_finalize_discovery (void);

// PIZ side
extern FlagsMgzip mgzip_piz_calculate_mgzip_flags (CompIType comp_i, Codec src_codec);
extern void bgzf_piz_set_txt_file_bgzf_info (FlagsMgzip mgzip_flags, bytes codec_info);
extern void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last);
extern void bgzf_write_finalize (void);

// misc
extern rom bgzf_library_name (MgzipLibraryType library, bool long_name);
extern rom gzstatus_name (GzStatus st);
extern void il1m_compress (void);
extern void bgzf_libdeflate_1_7_initialize (void);
extern void bgzf_sign (uint64_t disk_size, uint8_t *signature);
