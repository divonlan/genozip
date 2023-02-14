// ------------------------------------------------------------------
//   vblock.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"
#include "profiler.h"
#include "aes.h"
#include "context_struct.h"
#include "data_types.h"
#include "sections.h"

#define NUM_CODEC_BUFS 7       // bzlib2 compress requires 4 and decompress requires 2 ; lzma compress requires 7 and decompress 1
                               // if updating, also update array in codec_alloc()

// IMPORTANT: if changing fields in VBlockVCF, also update vb_release_vb
#define VBLOCK_COMMON_FIELDS \
    VBIType vblock_i;             /* number of VB within VCF file */\
    CompIType comp_i;             /* ZIP/PIZ: txt component within z_file that this VB belongs to  */ \
    bool is_eof;                  /* encountered EOF when reading this VB data from file */ \
    int32_t id;                   /* id of vb within the vb pool (-1 is the external vb) */\
    VBlockPoolP pool;             /* the VB pool to which this VB belongs */ \
    \
    /* compute thread stuff */ \
    ThreadId compute_thread_id;   /* id of compute thread currently processing this VB */ \
    rom compute_task;             /* task which the compute thread for this VB is performing */ \
    void (*compute_func)(VBlockP);/* compute thread entry point */\
    Mutex ready_for_compute;      /* threads_create finished initializeing this VB */\
    \
    DataType data_type;           /* type of this VB. In PIZ, this is the z_file data_type, NOT flag.out_dt */\
    DataType data_type_alloced;   /* type of this VB was allocated as. could be different that data_type, see vb_get_vb */\
    \
    /* memory management  */\
    Buffer buffer_list;           /* a buffer containing an array of pointers to all buffers allocated for this VB (either by the main thread or its compute thread). param=is_sorted */\
    \
    volatile DispatchStatus dispatch; /* line data is read, and dispatcher can dispatch this VB to a compute thread */\
    volatile bool is_processed;   /* thread completed processing this VB - it is ready for outputting */\
    volatile bool in_use;         /* this vb is in use */\
    \
    /* tracking lines */\
    Buffer lines;                 /* ZIP: An array of *DataLine* - the lines in this VB */\
                                  /* PIZ: array of (num_lines+1) x (char *) - pointer to within txt_data - start of each line. last item is BAFT(txt_data). */\
    BitsP is_dropped;             /* PIZ: a bitarray with a bit set is the line is marked for dropping by container_reconstruct_do */ \
    uint32_t num_lines_at_1_3, num_lines_at_2_3; /* ZIP VB=1 the number of lines segmented when 1/3 + 2/3 of estimate was reached  */\
    uint32_t debug_line_hash;     /* Seg: adler32 of line, used if Seg modifies line */\
    bool debug_line_hash_skip;    /* Seg: don't calculate debug_line_hash as line is skipped */\
    \
    /* tracking execution */\
    uint64_t vb_position_txt_file;/* ZIP/PIZ: position of this VB's data in the plain text file (i.e after decompression if the txt_file is compressed) */\
    uint64_t vb_bgzf_i;           /* ZIP: index into txt_file->bgzf_isizes of the first BGZF block of this VB */ \
    int32_t recon_size;           /* ZIP: actual size of txt if this VB is reconstructed in PRIMARY coordinates (inc. as ##primary_only in --luft) */\
                                  /* PIZ: expected reconstruction size in the coordinates of reconstruction */\
    int32_t txt_size;             /* ZIP: original size of of text data read from the file */ \
    uint32_t longest_line_len;    /* length of longest line of text line in this vb. calculated by seg_all_data_lines */\
    uint32_t sample_i;            /* ZIP/PIZ: VCF: current sample in line (0-based) */ \
    LineIType line_i;             /* ZIP/PIZ: current line in VB (0-based) being segmented/reconstructed */\
    int64_t rback_id;             /* ZIP: sequential number of current rollback point */ \
    uint32_t line_start;          /* ZIP/PIZ: position of start of line currently being segged / reconstructed in vb->txt_data */\
    uint32_t line_bgzf_uoffset;   /* ZIP: offset in uncompressed bgzf block of the start of the current line (current_bb_i) */  \
    \
    Digest digest;                /* ZIP/PIZ Adler32 (v9-13) and MD5: commulative digest up to and including this VB. Adler32 (v14+): standalone digest of this VB */ \
    Digest expected_digest;       /* PIZ: digest as transmitted in SectionHeaderVbHeader.digest */ \
    \
    DtTranslation translation;    /* PIZ: translation to be applies to this VB */ \
    union FlagsVbHeader flags;    /* ZIP: set by *_seg_finalize and consumed by zfile_compress_vb_header */ \
                                  /* PIZ: copied from SectionHeaderVbHeader.flags.vb_header */ \
    \
    rom drop_curr_line;           /* PIZ: line currently in reconstruction is to be dropped due a filter (value is filter name) */\
    uint32_t num_nondrop_lines;   /* PIZ: number of lines NOT dropped as a result of drop_curr_line */\
    uint8_t num_type1_subfields; \
    uint8_t num_type2_subfields; \
    RangeP range;                 /* ZIP: used for compressing the reference ranges. SAM PIZ: used */ \
    \
    uint32_t num_rollback_ctxs;   /* ZIP: Seg rollback contexts */ \
    ContextP rollback_ctxs[MEDIUM_CON_NITEMS]; \
    Buffer frozen_state;          /* PIZ: reconstruction state - frozen during reconstruct_peek */ \
    \
    /* data for dictionary and recon_plan compressing */ \
    char *fragment_start;        \
    uint32_t fragment_len;       \
    uint32_t fragment_num_words; \
    Context *fragment_ctx;       \
    \
    uint32_t refhash_layer;       /* create_ref && reading external reference: compressing/decompressing refhash */ \
    uint32_t refhash_start_in_layer; /* create_ref && reading external reference: compressing/decompressing refhash */ \
    \
    ProfilerRec profile; \
    \
    /* bgzf - for handling bgzf-compressed files */ \
    void *gzip_compressor;        /* Handle into libdeflate compressor or decompressor, or zlib's z_stream. Pointer to codec_bufs[].data */ \
    Buffer bgzf_blocks;           /* ZIP: an array of BgzfBlockZip tracking the decompression of bgzf blocks in scratch into txt_data.  */\
                                  /* PIZ: an array of BgzfBlockPiz */ \
    \
    /* random access, chrom, pos */ \
    Buffer ra_buf[2];             /* ZIP only: array of RAEntry [primary, luft] - copied to z_file at the end of each vb compression, then written as a SEC_RANDOM_ACCESS section at the end of the genozip file */\
    WordIndex chrom_node_index;   /* ZIP and PIZ: index and name of chrom of the current line. Note: since v12, this is redundant with last_int (CHROM) */ \
    STR(chrom_name);              /* since v12, this redundant with last_txtx/last_txt_len (CHROM) */ \
    uint32_t seq_len;             /* PIZ - last calculated seq_len (as defined by each data_type) */\
    \
    /* regions & filters */ \
    \
    /* PIZ: used by --show-coverage and --show-sex */ \
    Buffer coverage;              /* number of bases of each contig - exluding 'S' CIGAR, excluding reads flagged as Duplicate, Seconday arnd Failed filters */ \
    Buffer read_count;            /* number of mapped reads of each contig for show-coverage/idxstats (for show-coverage - excluding reads flagged as Duplicate, Seconday arnd Failed filters) */\
    Buffer unmapped_read_count;   \
    \
    /* crypto stuff */\
    Buffer spiced_pw;             /* used by crypt_generate_aes_key() */\
    int bi;                       /* used by AES */ \
    uint8_t aes_round_key[240];   /* for 256 bit aes */\
    uint8_t aes_iv[AES_BLOCKLEN]; \
    \
    /* file data */\
    Buffer z_data;                /* all headers and section data as read from disk */\
    union { \
    Buffer z_data_test;           /* ZIP: for use of codec_assign_best_codec */ \
    Buffer reread_prescription;   /* ZIP SAM/BAM DEPN: list of lines to be re-read at seg initialize */\
    }; \
    Buffer txt_data;              /* ZIP: txt_data as read from disk - either the txt header (in evb) or the VB data lines PIZ: reconstructed data */\
    Buffer z_section_headers;     /* PIZ and Pair-1 reading in ZIP-Fastq: an array of unsigned offsets of section headers within z_data */\
    Buffer scratch;               /* helper buffer: used by many functions. before usage, assert that its free, and buf_free after. */\
    int16_t z_next_header_i;      /* next header of this VB to be encrypted or decrypted */\
    \
    /* dictionaries stuff - we use them for 1. subfields with genotype data, 2. fields 1-9 of the VCF file 3. infos within the info field */\
    Did num_contexts;             /* total number of dictionaries of all types */\
    ContextArray contexts;    \
    DictIdtoDidMap d2d_map; /* map for quick look up of did_i from dict_id : 64K for key_map, 64K for alt_map */\
    \
    bool has_ctx_index; \
    ContextIndex ctx_index[MAX_DICTS]; /* PIZ VB: sorted index into contexts for binary-search lookup if d2d_map fails */\
    \
    /* reference stuff */ \
    Reference ref;                /* used by VBs created by dispatchers for uncompressing / compressing internal or external references. NOT used by VBs of the data type itself. */ \
    Buffer chrom2ref_map;         /* ZIP: mapping from user file chrom to alternate chrom in reference file (new chroms in this VB) - incides much vb->contexts[CHROM].nodes */\
    Buffer ol_chrom2ref_map;      /* ZIP: mapping from user file chrom to alternate chrom in reference file (chroms cloned) - incides much vb->contexts[CHROM].ol_nodes */\
    \
    /* reference range lookup caching */ \
    RangeP prev_range[2];         /* previous range returned by ref_seg_get_range */ \
    uint32_t prev_range_range_i;  /* range_i used to calculate previous range */ \
    WordIndex prev_range_chrom_node_index[2]; /* chrom used to calculate previous range */ \
    \
    union { \
    Buffer gencomp_lines;         /* ZIP: array of GencompLineIEntry: DVCF: lines rejected for liftover ; SAM-SA: primary/dependent lines */ \
    Buffer flusher_blocks;        /* PIZ writer vb */ \
    }; \
    \
    /* ref_iupac quick lookup */\
    ConstRangeP iupacs_last_range[2]; /* [0]=prim_ref [1]=gref */ \
    PosType64 iupacs_last_pos[2], iupacs_next_pos[2]; \
    \
    /* Information content stats - how many bytes does this section have more than the corresponding part of the vcf file */\
    Buffer show_headers_buf;      /* ZIP only: we collect header info, if --show-headers is requested, during compress, but show it only when the vb is written so that it appears in the same order as written to disk */\
    Buffer show_b250_buf;         /* ZIP only: for collecting b250 during generate - so we can print at onces without threads interspersing */\
    Buffer section_list_buf;      /* ZIP only: all the sections non-dictionary created in this vb. we collect them as the vb is processed, and add them to the zfile list in correct order of VBs. */\
    \
    /* Codec stuff */ \
    Codec codec_using_codec_bufs; /* codec currently using codec_bufs */\
    Buffer codec_bufs[NUM_CODEC_BUFS]; /* memory allocation for compressor so it doesn't do its own malloc/free */ \
    \
    /* used by CODEC_PBWT, CODEC_HAPMAT */ \
    uint32_t ht_per_line; \
    Context *ht_matrix_ctx; \
    \
    /* used by CODEC_PBWT */ \
    Context *runs_ctx, *fgrc_ctx; /* possibly diffrent did_i for different data types */\
    \
    /* copies of the values in flag, for flags that may change during the execution */\
    bool preprocessing;           /* PIZ: this VB is preprocessing, not reconstructing (SAM: loading SA Groups FASTA/FASTQ: grepping) */ \
    bool show_containers; \
    uint64_t ensure_dt_specific_fields_start_on_word_boundary; /* this ensures that the data-type-specific fields start at offset=sizeof(VBlock) - making it easier to release specific fields etc. this field is not used and remains 0. */ 

typedef struct VBlock {
    VBLOCK_COMMON_FIELDS
} VBlock;

#define in_assign_codec vb->z_data_test.prm8[0] // vb is currently in codec_assign_best_codec
#define peek_stack_level frozen_state.prm8[0]

extern void vb_cleanup_memory(void);
extern bool vb_is_valid (VBlockP vb);

extern void vb_destroy_vb_do (VBlockP *vb_p, rom func);
#define vb_destroy_vb(vb_p) vb_destroy_vb_do((vb_p), __FUNCTION__)

#define VB_ID_EVB           -1 // ID of VB used by main thread 
#define VB_ID_SEGCONF       -2 // ID of VB used by segconf_calculate
#define VB_ID_GCACHE_CREATE -3 
#define VB_ID_HCACHE_CREATE -4 
#define VB_ID_DEPN_SCAN     -5
#define VB_ID_COMPRESS_DEPN -6
#define NUM_NONPOOL_VBs      6
extern VBlockP vb_initialize_nonpool_vb(int vb_id, DataType dt, rom task);

extern void vb_release_vb_do (VBlockP *vb_p, rom task_name, rom func);
#define vb_release_vb(vb_p, task_name) vb_release_vb_do ((vb_p), (task_name), __FUNCTION__)
extern unsigned def_vb_size(DataType dt);

// -------------
// vb_pool stuff
// -------------

typedef struct VBlockPool {
    uint32_t size;              // size of memory allocated for this struct
    uint32_t num_vbs;           // length of array of pointers to VBlock
    uint32_t num_allocated_vbs; // number of VBlocks allocated ( <= num_vbs )
    uint32_t num_in_use;        // number of VBlocks currently in use ( <= num_allocated )
    VBlockP vb[];               // variable length
} VBlockPool;

extern void vb_create_pool (VBlockPoolType type);
extern VBlockPool *vb_get_pool (VBlockPoolType type, FailType soft_fail);
extern VBlockP vb_get_from_pool (VBlockPoolP pool, int32_t vb_id);
extern VBlockP vb_get_vb (VBlockPoolType type, rom task_name, VBIType vblock_i, CompIType comp_i);
extern void vb_destroy_pool_vbs (VBlockPoolType type);
extern uint32_t vb_pool_get_num_in_use (VBlockPoolType type);
extern bool vb_pool_is_full (VBlockPoolType type);
extern bool vb_pool_is_empty (VBlockPoolType type);

extern bool vb_is_processed (VBlockP vb);
extern void vb_set_is_processed (VBlockP vb);

//-----------------------------------------
// VBlock utilities
//-----------------------------------------

// NOT thread safe, use only in execution-terminating messages
extern rom err_vb_pos (void *vb);
