// ------------------------------------------------------------------
//   bgzf.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
 
#include "genozip.h"
#include "sections.h"

#define BGZF_MAX_BLOCK_SIZE 65536 // maximum block size of both compressed and uncompressed data of one block

// First 16 bytes of every BGZF block
#define BGZF_PREFIX_LEN 16
#define BGZF_PREFIX "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"

// BGZF EOF marker is simply an empty block, see https://samtools.github.io/hts-specs/SAMv1.pdf section 4.1.2
#define BGZF_EOF_LEN 28
#define BGZF_EOF BGZF_PREFIX "\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    
// data type of VBlock.bgzf_blocks
typedef struct BgzfBlockZip {
    int32_t txt_index;                    // index of uncompressed block within vb->txt_data. If there is passed-down data from previous VB/txt_header, then txt_index of the first block will be negative (see bgzf_copy_unconsumed_blocks)
    uint32_t txt_size        : 17;        // max value is BGZF_MAX_BLOCK_SIZE
    uint32_t is_decompressed : 1;         // has data been BGZF-decompressed by main thread
    uint32_t compressed_index, comp_size; // index within vb->scratch
} BgzfBlockZip;

typedef struct BgzfBlockPiz {
    int32_t txt_index, txt_size;          // index of uncompressed block within vb->txt_data. The first block index will be negative if there is passed-down unconsumed data
} BgzfBlockPiz;

extern void bgzf_libdeflate_1_7_initialize (void);
extern void bgzf_sign (uint64_t disk_size, uint8_t *signature);

//---------
// ZIP side
//---------

#define BGZF_BLOCK_SUCCESS        0
#define BGZF_BLOCK_GZIP_NOT_BGZIP -1
#define BGZF_BLOCK_IS_NOT_GZIP    -2
#define BGZF_ABRUBT_EOF           -3 // EOF without an EOF block
#define BGZF_BLOCK_TRUNCATED      -4 // likely file is truncated
extern int32_t bgzf_read_block (FileP file, uint8_t *block, uint32_t *block_size, FailType soft_fail);
extern void bgzf_uncompress_vb (VBlockP vb);
extern void bgzf_uncompress_one_block (VBlockP vb, BgzfBlockZip *bb);
extern void bgzf_reread_uncompress_vb_as_prescribed (VBlockP vb, FILE *file);
extern void bgzf_compress_bgzf_section (void);
extern void bgzf_zip_advance_index (VBlockP vb, uint32_t line_len);
extern int64_t bgzf_copy_unconsumed_blocks (VBlockP vb);
extern void bgzf_zip_init_vb (VBlockP vb);
extern void bgzf_insert_back_segconf_blocks (VBlockP vb);
extern void bgzf_return_segconf_blocks (VBlockP vb);

// library / level discovery
extern void bgzf_initialize_discovery (FileP file);
extern void bgzf_finalize_discovery (void);

#define consumed_by_prev_vb prm32[0] // bytes of the first BGZF block consumed by the prev VB or txt_header
#define current_bb_i        prm32[1] // index into vb->bgzf_blocks of first bgzf block of current line

//---------
// PIZ side
//---------

extern FlagsBgzf bgzf_piz_calculate_bgzf_flags (CompIType comp_i, Codec src_codec);
extern void bgzf_piz_set_txt_file_bgzf_info (FlagsBgzf bgzf_flags, bytes codec_info);
extern void bgzf_dispatch_compress (Dispatcher dispatcher, STRp (uncomp), CompIType comp_i, bool is_last);
extern void bgzf_write_finalize (void);

extern rom bgzf_library_name (BgzfLibraryType library, bool long_name);