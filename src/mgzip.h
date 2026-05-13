// ------------------------------------------------------------------
//   mgzip.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
 
#include "sections.h"
#include "flags.h"

// See: "Dev Docs" → "FASTQ gz formats"

#define IS_BGZF(codec)  ((codec)==CODEC_BGZF)
#define IS_MGZF(codec)  ((codec)==CODEC_MGZF)
#define IS_MGSP(codec)  ((codec)==CODEC_MGSP)
#define IS_IL1M(codec)  ((codec)==CODEC_IL1M)
#define IS_IL4M(codec)  ((codec)==CODEC_IL4M)
#define IS_EMFL(codec)  ((codec)==CODEC_EMFL)
#define IS_EMVL(codec)  ((codec)==CODEC_EMVL)
#define IS_GZBL(codec)  ((codec)==CODEC_GZBL)
#define IS_GZ(codec)    ((codec)==CODEC_GZ)
#define IS_BZ2(codec)   ((codec)==CODEC_BZ2)
#define IS_NONE(codec)  ((codec)==CODEC_NONE)
#define IS_MGZIP(codec) codec_args[codec].is_mgzip // multi-block gzip
#define IS_GZIP(codec)  (IS_MGZIP(codec) || IS_GZ(codec))

// note on MGSP: "gz block" in the comments below means, for MGSP, a group of gz blocks.
#define IS_IN_SYNC(codec)          (IS_MGZF(codec) || IS_MGSP(codec) || IS_EMVL(codec) || IS_EMFL(codec)) // codecs in which R1 and R2 gz blocks are guaranteed to contain whole, and precisely matching reads. Therefore, R2 gz-decompression can delegated to compute threads without further checks.
#define IS_VB_SIZE_BY_BLOCK(codec) (IS_MGZF(codec) || IS_EMVL(codec))             // codecs that are all of these: 1. variable-length 2. reads are never split between blocks 3. we use one VB per gz block
#define IS_VB_SIZE_BY_MGZIP(codec) (IS_VB_SIZE_BY_BLOCK(codec) || IS_MGSP(codec)) // like IS_VB_SIZE_BY_BLOCK, but VB can be a group of gz blocks
#define GZ_HEADER_HAS_BSIZE(codec) (IS_BGZF(codec) || IS_MGZF(codec)) // gz header contains bsize
#define VAR_LENGTH_NO_BSIZE(codec) (IS_GZBL(codec) || IS_EMVL(codec))
#define IS_EXACTABLE(codec)        (IS_BGZF(codec)) // codecs for which we can we discover library⁀level and can reconstruct exactly  

#define TXT_IS(c)               (txt_file->effective_codec == CODEC_##c)
#define TXT_IS_MGZIP            IS_MGZIP(txt_file->effective_codec)
#define TXT_IS_GZIP             IS_GZIP (txt_file->effective_codec)
#define TXT_IS_VB_SIZE_BY_BLOCK IS_VB_SIZE_BY_BLOCK(txt_file->effective_codec) 
#define TXT_IS_VB_SIZE_BY_MGZIP IS_VB_SIZE_BY_MGZIP(txt_file->effective_codec) 
#define TXT_IS_IN_SYNC          IS_IN_SYNC(txt_file->effective_codec)
#define TXT_GZ_HEADER_HAS_BSIZE GZ_HEADER_HAS_BSIZE(txt_file->effective_codec)
#define TXT_HAS_EMPTY_BLOCKS    (TXT_IS(BGZF) || TXT)
#define BGZF_MAX_BLOCK_SIZE ((uint32_t)(64 KB)) // maximum block size of both compressed and uncompressed data of one block
#define BGZF_MAX_CHUCK_SIZE ((uint32_t)(1 MB))  // max amount we read from disk at a time  
#define BGZF_HEADER_LEN 18
#define BGZF_PREFIX "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00" // First 16 bytes of every BGZF block (the gz header except for the last field "uint16_t bsize")
#define BGZF_PREFIX_LEN 16
#define BGZF_EOF_BSIZE_MINUS_1 "\x1b\x00"
#define BGZF_EOF_CDATA "\x03\x00" // there are many possible encodings for "empty block" in deflate. This is the one used BGZF_EOF per spec. igzip for example produces 01.00.00.FF.FF
#define BGZF_EMPTY_BLK_CRC_ISIZE "\x00\x00\x00\x00\x00\x00\x00\x00"
#define BGZF_EOF BGZF_PREFIX BGZF_EOF_BSIZE_MINUS_1 BGZF_EOF_CDATA BGZF_EMPTY_BLK_CRC_ISIZE // BGZF EOF marker is simply an empty block: SAM specification §4.1.2
#define BGZF_EOF_LEN    28

#define ILxM_PREFIX "\x1f\x8b\x08\x00\x00\x00\x00\x00" // first 8 bytes of ILxM gz header, to be followed by XFL (can be 0,2,4) and OS (3)
#define ILxM_PREFIX_LEN 8 
#define IL1M_ISIZE  "\x00\x00\x10\x00" // isize == 1MB in most blocks (except the last, and mid-blocks due to file concatenation)
#define IL4M_ISIZE  "\x00\x00\x40\x00" // isize == 4MB 

// MGI: a 32-bit version of BGZF 
#define MGZF_PREFIX "\x1f\x8b\x08\x14\x00\x00\x00\x00\x00\xff\x08\x00\x49\x47\x04\x00"
#define MGZF_PREFIX_LEN 16
#define MGZF_EOF MGZF_PREFIX "\x1f\x00\x00\x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
#define MGZF_EOF_LEN    31
#define MGZF_HEADER_LEN 29 // note: MGZF EOF header is only 21 bytes, because of zero-length comment field

// MGI: constant isize for all gz blocks that go into a particular VB (last block in group might slightly bigger)
#define MGSP_HEADER     "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03"
#define MGSP_EOF_LEN    20
#define MGSP_EOF MGSP_HEADER BGZF_EOF_CDATA BGZF_EMPTY_BLK_CRC_ISIZE

// MGI: single-block gz with non-gzip header (same as EMVL)
#define MGSB_HEADER "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"

// Element: constant isize of all gz blocks in file (last block might be smaller)
#define EMVL_HEADER "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"
#define EMVL_FIRST_BLOCK EMVL_HEADER "\x01\x00\x00\xff\xff\x00\x00\x00\x00\x00\x00\x00\x00" // EMVL files begin with this empty block and have no EOF block

typedef bool (*IsValidSize)(FileP, uint32_t proposed_isize, bool is_eof, bool discovering, bool *is_end_of_vb);

typedef struct BgzfBlockPiz {
    int32_t txt_index;  // index of uncompressed block within vb->txt_data. The first block index will be negative if there is passed-down unconsumed data
    int32_t txt_size;
    uint32_t gz_index;  // index of bgzf-compressed block within vb->comp_txt_data
    uint32_t gz_digest; // xxh3 digest of gz block (32b LSb) - inc. the header and footer of the block
} BgzfBlockPiz;

// ZIP side
typedef enum           { GZ_SUCCESS, GZ_SUCCESS_END_OF_VB, GZ_IS_OTHER_FORMAT,  GZ_MORE_DATA, GZ_NOT_GZIP, GZ_EOF_WITHOUT_EOF_BLOCK, GZ_TRUNCATED, NUM_GZ_STATUSES } GzStatus; // file is truncated
#define GZSTATUS_NAMES {   "SUCCESS",  "SUCCESS_END_OF_VB",  "IS_OTHER_FORMAT",   "MORE_DATA",  "NOT_GZIP"   "EOF_WITHOUT_EOF_BLOCK",  "TRUNCATED",                }

// data type of VBlock.gz_blocks and txt_file->unconsumed_mgzip_blocks : details of MGZIP blocks.
typedef struct GzBlockZip { // 16 bytes
    int32_t txt_index;             // index of uncompressed block within vb->txt_data. If there is passed-down data from previous VB/txt_header, then txt_index of the first block will be negative (see mgzip_copy_unconsumed_blocks)
    uint32_t txt_size        : 30; 
    uint32_t is_uncompressed : 1;  // true if data has been GZ-decompressed by main thread
    uint32_t is_eof          : 1;  // true if this is the last GZ-block in the file
    uint32_t gz_index, gz_size;    // index within vb->scratch
    uint32_t gz_digest;            // 32 LSb of xxh3 hash of the gz-data (inc. header and footer)
} GzBlockZip;

extern GzStatus mgzip_read_block_with_bsize (FileP file, bool discovering, Codec codec);
extern GzStatus mgzip_read_block_no_bsize (FileP file, bool discovering, Codec codec);
extern void mgzip_uncompress_vb (VBlockP vb, Codec codec);
extern void mgzip_uncompress_one_block (VBlockP vb, GzBlockZip *bb, Codec codec);
extern void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *file);
extern void mgzip_compress_SEC_GZ_sections (void);
extern void mgzip_zip_advance_index (VBlockP vb, uint32_t line_len);
extern int64_t mgzip_copy_unconsumed_blocks (VBlockP vb);
extern void mgzip_zip_init_vb (VBlockP vb);
extern void bgzf_insert_back_segconf_blocks (VBlockP vb);
extern void mgzip_return_segconf_blocks (VBlockP vb);
extern uint32_t mgzip_get_max_block_size (void);

extern void inc_disk_gz_uncomp_or_trunc_(FileP file, uint64_t inc, FUNCLINE);
#define inc_disk_gz_uncomp_or_trunc(file, inc) inc_disk_gz_uncomp_or_trunc_((file), (inc), __FUNCLINE)

// library / level discovery
extern void bgzf_initialize_discovery (FileP file);
extern void bgzf_finalize_discovery (void);
extern bool bgzf_read_and_uncomp_final_block (rom filename, qSTRp(uncomp));
extern void mgzip_set_is_exactable (FileP file, bool is_exactable, rom reason_why_not);
extern rom isal_error (int ret);

// PIZ side
extern FlagsMgzip mgzip_piz_calculate_mgzip_flags (CompIType comp_i, Codec src_codec);
extern void mgzip_piz_set_txt_file_info (FlagsMgzip mgzip_flags, uint32_t OLD_gz_size_3LSB);
extern void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last, bool is_txt_header);
extern void bgzf_write_finalize (void);
extern void bgzf_compress_tbi (void);

// misc
extern rom bgzf_library_name (MgzipLibraryType library, bool long_name);
extern StrText bgzf_lib_name_level (FlagsMgzip mgzip_flags);
extern rom gzstatus_name (GzStatus st);
extern void generate_il1m (void);
extern void bgzf_libdeflate_1_7_initialize (void);

extern void show_gz (rom filename);
extern void dump_gz_block (rom filename);

// private within mgzip*.c
extern rom NON_EXACT_ERROR;
extern const FlagsMgzip bgzf_recompression_levels[1+MAX_FLAG_BGZF];
