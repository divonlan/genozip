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
#include "context_struct.h"
#include "data_types.h"
#include "sections.h"

#define NUM_CODEC_BUFS 7       // bzlib2 compress requires 4 and decompress requires 2 ; lzma compress requires 7 and decompress 1
                               // if updating, also update array in codec_alloc()

typedef enum { GS_READ, GS_TEST, GS_UNCOMPRESS } GrepStages;

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
#define VBLOCK_COMMON_FIELDS \
    uint32_t vblock_i;         /* number of variant block within VCF file */\
    int id;                    /* id of vb within the vb pool (-1 is the external vb) */\
    \
    /* compute thread stuff */ \
    ThreadId compute_thread_id;/* id of compute thread currently processing this VB */ \
    const char *compute_task;  /* task which the compute thread for this VB is performing */ \
    void (*compute_func)(VBlockP); /* compute thread entry point */\
    Mutex vb_ready_for_compute_thread; /* threads_create finished initializeing this VB */\
    \
    DataType data_type;        /* type of this VB */\
    \
    /* memory management  */\
    Buffer buffer_list;        /* a buffer containing an array of pointers to all buffers allocated for this VB (either by the main thread or its compute thread) */\
    \
    bool ready_to_dispatch;    /* line data is read, and dispatcher can dispatch this VB to a compute thread */\
    bool is_processed;         /* thread completed processing this VB - it is ready for outputting */\
    bool in_use;               /* this vb is in use */\
    \
    /* tracking lines */\
    Buffer lines;              /* ZIP: An array of *DataLine* - the lines in this VB */\
                               /* PIZ: array of (num_lines+1) x (char *) - pointer to within txt_data - start of each line. last item is AFTERENT(txt_data). */\
    uint16_t first_line;       /* PIZ: line number in source txt file (counting from 1), of this variant block */\
                               /* ZIP: used for optimize_DESC in FASTQ */ \
    uint32_t num_lines_at_1_3, num_lines_at_2_3; /* ZIP VB=1 the number of lines segmented when 1/3 + 2/3 of estimate was reached  */\
    \
    /* tracking execution */\
    uint64_t vb_position_txt_file; /* position of this VB's data in the plain text file (i.e after decompression if the txt_file is compressed) */\
    int32_t recon_size;        /* ZIP: actual size of txt read from file, modified if --optimize/--chain as expected to be reconstructed in PRIMARY coordinates */\
                               /* PIZ: expected reconstruction size in the coordinates of reconstruction */\
    int32_t recon_size_luft;   /* ZIP only: expected reconstruction size in LUFT coordinates */ \
    uint32_t recon_num_lines;  /* PIZ: expected number of non-dropped lines in the default reconstruction size in the user coordinates of reconstruction */\
                               /* ZIP: dual-coordinate files only: DC_PRIMARY + DC_BOTH lines in this VB  */ \
    uint32_t recon_num_lines_luft; /* ZIP only, dual coordinates files only: DC_LUFT + DC_BOTH lines in this VB */ \
    int32_t txt_size;          /* ZIP: original size of of text data read from the file */ \
    uint32_t longest_line_len; /* length of longest line of text line in this vb. calculated by seg_all_data_lines */\
    uint32_t line_i;           /* ZIP: current line in VB (0-based) being segmented PIZ: current line in txt file */\
    uint64_t line_start;       /* PIZ: position of start of line currently being reconstructed in vb->txt_data */\
    \
    Digest digest_so_far;      /* partial calculation of MD5 up to and including this VB */ \
    uint32_t component_i;      /* PIZ: 0-based txt component within z_file that this VB belongs to */ \
    DtTranslation translation; /* PIZ: translation to be applies to this VB */ \
    \
    const char *drop_curr_line;/* PIZ: line currently in reconstruction is to be dropped due a filter (value is filter name) */\
    uint32_t num_nondrop_lines;/* PIZ: number of lines NOT dropped as a result of drop_curr_line */\
    GrepStages grep_stages;    /* PIZ: tell piz_is_skip_section what to skip in case of --grep */\
    uint8_t num_type1_subfields; \
    uint8_t num_type2_subfields; \
    RangeP range;              /* ZIP: used for compressing the reference ranges */ \
    uint32_t range_num_set_bits;  /* ZIP: main thread telling compute thread to how many bits are set in range.is_set */ \
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
    Buffer ra_buf[2];          /* ZIP only: array of RAEntry [primary, luft] - copied to z_file at the end of each vb compression, then written as a SEC_RANDOM_ACCESS section at the end of the genozip file */\
    WordIndex chrom_node_index;/* ZIP and PIZ: index and name of chrom of the current line. Note: since v12, this is redundant with last_int (CHROM) */ \
    const char *chrom_name;    /* since v12, this redundant with last_txtx (CHROM) */ \
    unsigned chrom_name_len;   /* since v12, this redundant with last_txt_len (CHROM) */\
    uint32_t seq_len;          /* PIZ - last calculated seq_len (as defined by each data_type) */\
    \
    /* regions & filters */ \
    \
    /* used by --show-coverage and --show-sex */ \
    Buffer coverage;           /* number of bases of each contig - exluding 'S' CIGAR, excluding reads flagged as Duplicate, Seconday arnd Failed filters */ \
    Buffer read_count;         /* number of mapped reads of each contig for show-coverage/idxstats (for show-coverage - excluding reads flagged as Duplicate, Seconday arnd Failed filters) */\
    Buffer unmapped_read_count;   \
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
    /* liftover stuff */ \
    Coords vb_coords;          /* ZIP: DC_PRIMARY, DC_LUFT or DC_BOTH */ \
                               /* PIZ: DC_PRIMARY or DC_LUFT - influenced by FlagsVbHeader.coords and flag.luft */ \
    Coords line_coords;        /* Seg: coords of current line - DC_PRIMARY or DC_LUFT */ \
    Buffer liftover;           /* ZIP: map from chrom_node_index (not word index!) to entry in chain_index */ \
    uint32_t pos_aln_i;        /* ZIP: chain alignment of POS (used to compare to that of END) */\
    Buffer lo_rejects[2];      /* ZIP generating a dual-coordinates file: txt lines rejected for liftover */ \
    int32_t reject_bytes;      /* ZIP of a Luft file: number of bytes of reject data in this VB (data originating from ##primary_only/##luft_only) */ \
    bool is_rejects_vb;        /* PIZ: this is a VB of rejects variants for header ##primary_only/##luft_only */ \
    bool is_unsorted[2];       /* ZIP: line order of this VB[primary, luft] is unsorted */ \
    \
    /* Information content stats - how many bytes does this section have more than the corresponding part of the vcf file */\
    Buffer show_headers_buf;   /* ZIP only: we collect header info, if --show-headers is requested, during compress, but show it only when the vb is written so that it appears in the same order as written to disk */\
    Buffer show_b250_buf;      /* ZIP only: for collecting b250 during generate - so we can print at onces without threads interspersing */\
    Buffer section_list_buf;   /* ZIP only: all the sections non-dictionary created in this vb. we collect them as the vb is processed, and add them to the zfile list in correct order of VBs. */\
    \
    /* Codec stuff */ \
    Codec codec_using_codec_bufs; /* codec currently using codec_bufs */\
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

// VBLOCK_COMMON_LINES_ZIP needs to be at the begining of *ZipDataLine of data types that support dual-coordinates.
#define VBLOCK_COMMON_LINES_ZIP \
    WordIndex chrom_index[2];   /* Seg: enter as node_index ; Merge: convert to word_index */ \
    PosType pos[2];             /* arrays of [2] - { primary-coord, luft-coord } */ \

typedef struct ZipDataLine {
    VBLOCK_COMMON_LINES_ZIP
} ZipDataLine;

extern void vb_cleanup_memory(void);
extern VBlock *vb_get_vb (const char *task_name, uint32_t vblock_i);
extern bool vb_has_free_vb (void);
extern void vb_destroy_vb (VBlockP *vb_p);

#define EVB -1 // ID of VB used by main thread 
extern VBlockP vb_initialize_nonpool_vb(int vb_id);
extern void vb_release_vb_do (VBlock **vb_p, const char *func);
#define vb_release_vb(vb_p) vb_release_vb_do (vb_p, __FUNCTION__)
extern void vb_destroy_pool_vbs (void);

// -------------
// vb_pool stuff
// -------------
typedef struct {
    unsigned num_vbs; // length of array of pointers to VBlock
    unsigned num_allocated_vbs; // number of VBlocks allocated ( <= num_vbs )
    VBlock *vb[]; // variable length
} VBlockPool;
extern void vb_create_pool (unsigned num_vbs);
extern VBlockPool *vb_get_pool(void);
#define vb_get_from_pool(vb_i) (((vb_i) == EVB) ? evb : vb_pool->vb[vb_i])

#endif

//-----------------------------------------
// VBlock utilities
//-----------------------------------------

// NOT thread safe, use only in execution-terminating messages
extern const char *err_vb_pos (void *vb);
