// ------------------------------------------------------------------
//   bgzf.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
 
#include "sections.h"

#define BGZF_DEFAULT_LEVEL 2 // PIZ: used if --bgzf is not specified (it is actually faster than 1 if also writing to disk)

#define BGZF_MAX_BLOCK_SIZE ((uint32_t)(64 KB)) // maximum block size of both compressed and uncompressed data of one block
#define BGZF_MAX_CHUCK_SIZE ((uint32_t)(1 MB))  // max amount we read from disk at a time  

// First 16 bytes of every BGZF block
#define BGZF_PREFIX_LEN 16
#define BGZF_PREFIX "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"

// BGZF EOF marker is simply an empty block, see https://samtools.github.io/hts-specs/SAMv1.pdf section 4.1.2
#define BGZF_EOF_LEN 28
#define BGZF_EOF BGZF_PREFIX "\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    
// maximum block size of uncompressed data, and we are going to assume that since GZIL is only 
// for FASTQ, and FASTQ is quite compressible, (GZIL_MAX_BLOCK_SIZE-GZIL_HEADER_LEN) it is an upper limit on GZIL-compressed data   
#define GZIL_MAX_BLOCK_SIZE ((uint32_t)(1 MB))  

#define GZIL_HEADER "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03"
#define GZIL_HEADER_LEN 10

#define GZIL_ISIZE "\x00\x00\x10\x00" // isize == 1MB in all blocks except the last
#define GZIL_ISIZE_LEN 4

typedef struct BgzfBlockPiz {
    int32_t txt_index, txt_size;          // index of uncompressed block within vb->txt_data. The first block index will be negative if there is passed-down unconsumed data
} BgzfBlockPiz;

// ZIP side
typedef enum           { GZ_SUCCESS, GZ_IS_GZIP_NOT_BGZF, GZ_IS_NOT_GZIL, GZ_IS_NOT_GZIP, GZ_EOF_WITHOUT_EOF_BLOCK, GZ_TRUNCATED, NUM_GZ_STATUSES } GzStatus; // file is truncated
#define GZSTATUS_NAMES {   "SUCCESS",  "IS_GZIP_NOT_BGZF",  "IS_NOT_GZIL",  "IS_NOT_GZIP",  "EOF_WITHOUT_EOF_BLOCK",  "TRUNCATED" }

// data type of VBlock.gz_blocks and txt_file->unconsumed_bgz_blocks : details of BGZF/GZIL blocks.
typedef struct GzBlockZip {
    int32_t txt_index;                    // index of uncompressed block within vb->txt_data. If there is passed-down data from previous VB/txt_header, then txt_index of the first block will be negative (see bgz_copy_unconsumed_blocks)
    uint32_t txt_size        : 30; 
    uint32_t is_decompressed : 1;         // true if data has been GZ-decompressed by main thread
    uint32_t is_eof          : 1;         // true if this is the last GZ-block in the file
    uint32_t compressed_index, comp_size; // index within vb->scratch
} GzBlockZip;

extern GzStatus bgzf_read_block (FileP file, bool discovering);
extern GzStatus gzil_read_block (FileP file, bool discovering, bool *is_eof);
extern void bgz_uncompress_vb (VBlockP vb, Codec codec);
extern void bgz_uncompress_one_block (VBlockP vb, GzBlockZip *bb, Codec codec);
extern void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *file);
extern void bgzf_compress_bgzf_section (void);
extern void bgz_zip_advance_index (VBlockP vb, uint32_t line_len);
extern int64_t bgz_copy_unconsumed_blocks (VBlockP vb);
extern void bgz_zip_init_vb (VBlockP vb);
extern void bgzf_insert_back_segconf_blocks (VBlockP vb);
extern void bgz_return_segconf_blocks (VBlockP vb);

extern void inc_disk_gz_uncomp_or_trunc_(FileP file, uint64_t inc, FUNCLINE);
#define inc_disk_gz_uncomp_or_trunc(file, inc) inc_disk_gz_uncomp_or_trunc_((file), (inc), __FUNCLINE)

// library / level discovery
extern void bgzf_initialize_discovery (FileP file);
extern void bgzf_finalize_discovery (void);

// PIZ side
extern FlagsBgzf bgzf_piz_calculate_bgzf_flags (CompIType comp_i, Codec src_codec);
extern void bgzf_piz_set_txt_file_bgzf_info (FlagsBgzf bgzf_flags, bytes codec_info);
extern void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last);
extern void bgzf_write_finalize (void);

// misc
extern rom bgzf_library_name (BgzfLibraryType library, bool long_name);
extern rom gzstatus_name (GzStatus st);
extern void gzil_compress (void);
extern void bgzf_libdeflate_1_7_initialize (void);
extern void bgzf_sign (uint64_t disk_size, uint8_t *signature);
