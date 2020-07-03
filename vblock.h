// ------------------------------------------------------------------
//   vblock.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VBLOCK_INCLUDED
#define VBLOCK_INCLUDED

#include "genozip.h"
#include "buffer.h"
#include "profiler.h"
#include "aes.h"
#include "context.h"
#include "bit_array.h"

#ifndef DID_I_NONE // also defined in move_to_front.h
#define DID_I_NONE   255
#endif
typedef struct SubfieldMapper {
    uint8_t num_subfields;        // (uint8_t)NIL if this mapper is not defined
    uint8_t did_i[MAX_SUBFIELDS]; // array in the order the subfields appears in FORMAT or INFO - each an index into vb->contexts[]
} SubfieldMapper;

#define MAPPER_CTX(mapper,sf) (((mapper)->did_i[(sf)] != (uint8_t)NIL) ? &vb->contexts[(mapper)->did_i[(sf)]] : NULL)

#define NUM_COMPRESS_BUFS 7   // bzlib2 compress requires 4 and decompress requires 2 ; lzma compress requires 7 and decompress 1

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
    Buffer lines;              /* An array of *DataLine* - the lines in this VB */\
    uint32_t first_line;       /* PIZ only: line number in VCF file (counting from 1), of this variant block */\
    uint32_t num_lines_at_1_3, num_lines_at_2_3; /* ZIP VB=1 the number of lines segmented when 1/3 + 2/3 of estimate was reached  */\
    \
    /* tracking execution */\
    uint64_t vb_position_txt_file; /* position of this VB's data in the plain text file (i.e after decompression if the txt_file is compressed) */\
    int32_t vb_data_size;      /* ZIP: actual size of txt read file file ; PIZ: expected size of decompressed txt. Might be different than original if --optimize is used. */\
    uint32_t vb_data_read_size;/* ZIP only: amount of data read in txtfile_read_block() (either plain VCF or gz or bz2) for this VB */\
    uint32_t longest_line_len; /* length of longest line of text line in this vb. calculated by seg_all_data_lines */\
    uint32_t line_i;           /* ZIP: current line in VB (0-based) being segmented PIZ: current line in txt file */\
    bool dont_show_curr_line;  /* PIZ only - line currently in reconstruction is grepped out due to --grep or --regions and should not be displayed */\
    GrepStages grep_stages;    /* PIZ: tell piz_is_skip_section what to skip in case of --grep */\
    uint8_t num_type1_subfields; \
    uint8_t num_type2_subfields; \
    RangeP range;              /* ZIP: used for compressing the reference ranges */ \
    uint32_t range_num_set_bits;  /* ZIP: I/O thread telling compute thread to how many bits are set in range.is_set */ \
    \
    ProfilerRec profile; \
    \
    /* random access, chrom, pos */ \
    Buffer ra_buf;             /* ZIP only: array of RAEntry - copied to z_file at the end of each vb compression, then written as a SEC_RANDOM_ACCESS section at the end of the genozip file */\
    int32_t chrom_node_index;  /* ZIP and PIZ: index and name of chrom of the current line */ \
    const char *chrom_name;    \
    unsigned chrom_name_len; \
    uint32_t seq_len;          /* PIZ - last calculated seq_len (as defined by each data_type) */\
                               \
    /* regions & filters */ \
    Buffer region_ra_intersection_matrix;  /* PIZ: a byte matrix - each row represents an ra in this vb, and each column is a region specieid in the command. the cell contains 1 if this ra intersects with this region */\
    \
    /* crypto stuff */\
    Buffer spiced_pw;  /* used by crypt_generate_aes_key() */\
    uint8_t aes_round_key[240];/* for 256 bit aes */\
    uint8_t aes_iv[AES_BLOCKLEN]; \
    int bi; \
    \
    /* file data */\
    Buffer z_data;                    /* all headers and section data as read from disk */\
    \
    Buffer txt_data;                  /* ZIP only: txt_data as read from disk - either the txt header (in evb) or the VB data lines */\
    uint32_t txt_data_next_offset;    /* we re-use txt_data memory to overlay stuff in segregate */\
    Buffer txt_data_spillover;        /* when re-using txt_data, if it is too small, we spill over to this buffer */\
    \
    int16_t z_next_header_i;          /* next header of this VB to be encrypted or decrypted */\
    \
    Buffer z_section_headers;         /* PIZ only: an array of unsigned offsets of section headers within z_data */\
    \
    Buffer compressed;                /* used by various zfile functions */\
    \
    /* dictionaries stuff - we use them for 1. subfields with genotype data, 2. fields 1-9 of the VCF file 3. infos within the info field */\
    uint32_t num_dict_ids;            /* total number of dictionaries of all types */\
    Context contexts[MAX_DICTS];    \
    uint8_t dict_id_to_did_i_map[65536];       /* map for quick look up of did_i from dict_id */\
    \
    /* ZIP only: reference range lookup caching */ \
    RangeP prev_range; /* previous range returned by ref_zip_get_locked_range */ \
    uint32_t prev_range_range_i; /* range_i used to calculate previous range */ \
    uint32_t prev_range_chrom_node_index; /* chrom used to calculate previous range */ \
    \
    /* Information content stats - how many bytes does this section have more than the corresponding part of the vcf file */\
    Buffer show_headers_buf;                   /* ZIP only: we collect header info, if --show-headers is requested, during compress, but show it only when the vb is written so that it appears in the same order as written to disk */\
    Buffer show_b250_buf;                      /* ZIP only: for collecting b250 during generate - so we can print at onces without threads interspersing */\
    Buffer section_list_buf;                   /* ZIP only: all the sections non-dictionary created in this vb. we collect them as the vb is processed, and add them to the zfile list in correct order of VBs. */\
    \
    /* Compressor stuff */ \
    Buffer compress_bufs[NUM_COMPRESS_BUFS];   /* memory allocation for compressor so it doesn't do its own malloc/free */ 

typedef struct VBlock {
    VBLOCK_COMMON_FIELDS
} VBlock;

extern void vb_cleanup_memory(void);
extern VBlock *vb_get_vb (unsigned vblock_i);
extern void vb_external_vb_initialize(void);
extern void vb_release_vb (VBlock *vb);

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
