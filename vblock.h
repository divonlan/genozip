// ------------------------------------------------------------------
//   vblock.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VBLOCK_INCLUDED
#define VBLOCK_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "profiler.h"
#include "aes.h"
#include "context.h"
#include "bit_array.h"
#include "data_types.h"

#define NUM_CODEC_BUFS 7       // bzlib2 compress requires 4 and decompress requires 2 ; lzma compress requires 7 and decompress 1

typedef enum { GS_READ, GS_TEST, GS_UNCOMPRESS } GrepStages;

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
#define VBLOCK_COMMON_FIELDS \
    uint32_t vblock_i;         /* number of variant block within VCF file */\
    int id;                    /* id of vb within the vb pool (-1 is the external vb) */\
    DataType data_type;        /* type of this VB */\
    \
    /* memory management  */\
    Buffer buffer_list;        /* a buffer containing an array of pointers to all buffers allocated for this VB (either by the I/O thread or its compute thread) */\
    \
    bool ready_to_dispatch;    /* line data is read, and dispatcher can dispatch this VB to a compute thread */\
    bool is_processed;         /* thread completed processing this VB - it is ready for outputting */\
    bool in_use;               /* this vb is in use */\
    \
    /* tracking lines */\
    Buffer lines;              /* An array of *DataLine* - the lines in this VB */\
    uint32_t first_line;       /* PIZ: line number in txt file (counting from 1), of this variant block */\
    uint32_t num_lines_at_1_3, num_lines_at_2_3; /* ZIP VB=1 the number of lines segmented when 1/3 + 2/3 of estimate was reached  */\
    \
    /* tracking execution */\
    uint64_t vb_position_txt_file; /* position of this VB's data in the plain text file (i.e after decompression if the txt_file is compressed) */\
    int32_t vb_data_size;      /* ZIP: actual size of txt read from file ; PIZ: expected size of decompressed txt. Might be different than original if --optimize is used. */\
    uint32_t longest_line_len; /* length of longest line of text line in this vb. calculated by seg_all_data_lines */\
    uint32_t line_i;           /* ZIP: current line in VB (0-based) being segmented PIZ: current line in txt file */\
    uint64_t line_start;       /* PIZ: position of start of line currently being reconstructed in vb->txt_data */\
    Digest digest_so_far;      /* partial calculation of MD5 up to and including this VB */ \
    uint32_t component_i;      /* PIZ: 0-based txt component within z_file that this VB belongs to */ \
    \
    bool dont_show_curr_line;  /* PIZ: line currently in reconstruction is grepped out due to --grep or --regions and should not be displayed */\
    GrepStages grep_stages;    /* PIZ: tell piz_is_skip_section what to skip in case of --grep */\
    uint8_t num_type1_subfields; \
    uint8_t num_type2_subfields; \
    RangeP range;              /* ZIP: used for compressing the reference ranges */ \
    uint32_t range_num_set_bits;  /* ZIP: I/O thread telling compute thread to how many bits are set in range.is_set */ \
    \
    /* data for dictionary compressing */ \
    char *fragment_start;        \
    uint32_t fragment_len;       \
    uint32_t fragment_num_words; \
    Context *fragment_ctx;       \
    Codec fragment_codec;        \
    \
    uint32_t refhash_layer;    /* create_ref && reading external reference: compressing/decompressing refhash */ \
    uint32_t refhash_start_in_layer;     /* create_ref && reading external reference: compressing/decompressing refhash */ \
    \
    ProfilerRec profile; \
    \
    /* bgzf - for handling bgzf-compressed files */ \
    void *gzip_compressor;     /* Handle into libdeflate compressor or decompressor, or zlib's z_stream. Pointer to codec_bufs[].data */ \
    Buffer bgzf_blocks;        /* ZIP: an array of BgzfBlockZip tracking the decompression of blocks into txt_data */\
    \
    /* random access, chrom, pos */ \
    Buffer ra_buf;             /* ZIP only: array of RAEntry - copied to z_file at the end of each vb compression, then written as a SEC_RANDOM_ACCESS section at the end of the genozip file */\
    WordIndex chrom_node_index;/* ZIP and PIZ: index and name of chrom of the current line */ \
    const char *chrom_name;    \
    unsigned chrom_name_len; \
    uint32_t seq_len;          /* PIZ - last calculated seq_len (as defined by each data_type) */\
    \
    /* regions & filters */ \
    Buffer region_ra_intersection_matrix;  /* PIZ: a byte matrix - each row represents an ra in this vb, and each column is a region specieid in the command. the cell contains 1 if this ra intersects with this region */\
    \
    /* used by --show-coverage and --show-sex */ \
    Buffer coverage; \
    \
    /* crypto stuff */\
    Buffer spiced_pw;          /* used by crypt_generate_aes_key() */\
    uint8_t aes_round_key[240];/* for 256 bit aes */\
    uint8_t aes_iv[AES_BLOCKLEN]; \
    int bi; \
    \
    /* file data */\
    Buffer z_data;             /* all headers and section data as read from disk */\
    \
    Buffer txt_data;           /* ZIP only: txt_data as read from disk - either the txt header (in evb) or the VB data lines */\
    \
    int16_t z_next_header_i;   /* next header of this VB to be encrypted or decrypted */\
    \
    Buffer z_section_headers;  /* PIZ and Pair-1 reading in ZIP-Fastq: an array of unsigned offsets of section headers within z_data */\
    \
    Buffer compressed;         /* helper buffer for writing to/from zfile: used by various functions. user must assert that its free before use, and buf_free after use. */\
    \
    /* dictionaries stuff - we use them for 1. subfields with genotype data, 2. fields 1-9 of the VCF file 3. infos within the info field */\
    DidIType num_contexts;     /* total number of dictionaries of all types */\
    Context contexts[MAX_DICTS];    \
    DidIType dict_id_to_did_i_map[65536];       /* map for quick look up of did_i from dict_id */\
    \
    /* ZIP only: reference range lookup caching */ \
    RangeP prev_range;         /* previous range returned by ref_seg_get_locked_range */ \
    uint32_t prev_range_range_i; /* range_i used to calculate previous range */ \
    WordIndex prev_range_chrom_node_index; /* chrom used to calculate previous range */ \
    \
    /* Information content stats - how many bytes does this section have more than the corresponding part of the vcf file */\
    Buffer show_headers_buf;   /* ZIP only: we collect header info, if --show-headers is requested, during compress, but show it only when the vb is written so that it appears in the same order as written to disk */\
    Buffer show_b250_buf;      /* ZIP only: for collecting b250 during generate - so we can print at onces without threads interspersing */\
    Buffer section_list_buf;   /* ZIP only: all the sections non-dictionary created in this vb. we collect them as the vb is processed, and add them to the zfile list in correct order of VBs. */\
    \
    /* Codec stuff */ \
    Buffer codec_bufs[NUM_CODEC_BUFS];   /* memory allocation for compressor so it doesn't do its own malloc/free */ \
    \
    /* used by CODEC_ACGT (For SEQ) */ \
    bool has_non_agct;         /* ZIP only */ \
    \
    /* used by CODEC_PBWT, CODEC_HAPMAT and CODEC_GTSHARK */ \
    uint32_t ht_per_line; \
    Context *ht_matrix_ctx; \
    \
    /* used by CODEC_PBWT */ \
    Context *runs_ctx, *fgrc_ctx;

typedef struct VBlock {
    VBLOCK_COMMON_FIELDS
} VBlock;

extern void vb_cleanup_memory(void);
extern VBlock *vb_get_vb (const char *task_name, uint32_t vblock_i);
extern void vb_initialize_evb(void);
extern void vb_release_vb (VBlock *vb);
extern void vb_destroy_all_vbs (void);

typedef struct {
    unsigned num_vbs; // length of array of pointers to VBlock
    unsigned num_allocated_vbs; // number of VBlocks allocated ( <= num_vbs )
    VBlock *vb[]; // variable length
} VBlockPool;
extern void vb_create_pool (unsigned num_vbs);
extern VBlockPool *vb_get_pool(void);
#endif

//-----------------------------------------
// VBlock utilities
//-----------------------------------------

// NOT thread safe, use only in execution-terminating messages
extern const char *err_vb_pos (void *vb);
