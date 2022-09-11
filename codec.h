// ------------------------------------------------------------------
//   codec.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

#include "lzma/7zTypes.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "data_types.h"

#define MIN_LEN_FOR_COMPRESSION 50 // less that this size, and compressed size is typically larger than uncompressed size

#define COMPRESS(func)                                                              \
    bool func (VBlockP vb,                                                          \
               ContextP ctx, /* NULL if not compressing a context */                \
               SectionHeader *header,                                               \
               rom uncompressed,         /* option 1 - compress contiguous data */  \
               uint32_t *uncompressed_len,                                          \
               LocalGetLineCB get_line_cb,  /* option 2 - call back to get lines */ \
               char *compressed, uint32_t *compressed_len/* in/out */,              \
               bool soft_fail, rom name)    

typedef COMPRESS (CodecCompress);

#define UNCOMPRESS(func)                                             \
    void func (VBlockP vb, Codec codec, uint8_t param,               \
               rom compressed, uint32_t compressed_len,              \
               BufferP uncompressed_buf, uint64_t uncompressed_len,  \
               Codec sub_codec,                                      \
               rom name)

typedef UNCOMPRESS (CodecUncompress);

typedef uint32_t CodecEstSizeFunc (Codec codec, uint64_t uncompressed_len);

#define CODEC_RECONSTRUCT(func) \
    void func (VBlockP vb, Codec codec, ContextP ctx, STRp(snip))

typedef CODEC_RECONSTRUCT (CodecReconstruct);

typedef struct {
    bool             is_simple;  // a simple codec is one that is compressed into a single section in one step
    const char       *name;
    const char       *ext;       // extensions by compression type. + if it adds to the name ; - if it replaces the extension of the uncompress name
    CodecCompress    *compress;
    CodecUncompress  *uncompress;
    CodecReconstruct *reconstruct;
    CodecEstSizeFunc *est_size;
} CodecArgs;

#define NA0 "N/A"
#define NA1 codec_compress_error
#define NA2 codec_uncompress_error
#define NA3 codec_reconstruct_error
#define NA4 codec_est_size_default
#define USE_SUBCODEC NULL

#define CODEC_ARGS { /* aligned with Codec defined in genozip.h */ \
/*  simp name    ext       compress                  uncompress                reconstruct                est_size                 */ \
    { 1, "N/A",  "+",      NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 1, "NONE", "+",      codec_none_compress,      codec_none_uncompress,    NA3,                       codec_none_est_size      }, \
    { 1, "GZ",   "+.gz",   NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 1, "BZ2",  "+.bz2",  codec_bz2_compress,       codec_bz2_uncompress,     NA3,                       NA4                      }, \
    { 1, "LZMA", "+",      codec_lzma_compress,      codec_lzma_uncompress,    NA3,                       codec_none_est_size      }, \
    { 1, "BSC",  "+",      codec_bsc_compress,       codec_bsc_uncompress,     NA3,                       codec_bsc_est_size       }, \
    { 1, "RANB", "+",      codec_RANB_compress,      codec_rans_uncompress,    NA3,                       codec_RANB_est_size      }, \
    { 1, "RANW", "+",      codec_RANW_compress,      codec_rans_uncompress,    NA3,                       codec_RANW_est_size      }, /* STRIPE */\
    { 1, "RANb", "+",      codec_RANb_compress,      codec_rans_uncompress,    NA3,                       codec_RANb_est_size      }, /* PACK */\
    { 1, "RANw", "+",      codec_RANw_compress,      codec_rans_uncompress,    NA3,                       codec_RANw_est_size      }, /* STRIPE & PACK */\
    { 0, "ACGT", "+",      codec_acgt_compress,      codec_acgt_uncompress,    NA3,                       codec_complex_est_size   }, \
    { 0, "XCGT", "+",      USE_SUBCODEC,             codec_xcgt_uncompress,    NA3,                       NA4                      }, \
    { 0, "HAPM", "+",      NA1,                      USE_SUBCODEC,             codec_hapmat_reconstruct,  NA4,                     }, /* HapMat was discontinued and replaced by PBWT. We keep it for decompressing old VCF files */ \
    { 0, "DOMQ", "+",      codec_domq_compress,      USE_SUBCODEC,             codec_domq_reconstruct,    codec_complex_est_size,  }, \
    { 0, "GTSH", "+",      NA1,                      codec_gtshark_uncompress, codec_pbwt_reconstruct,    NA4,                     }, /* gtshark discontinued in v12. keep for displaying an error */\
    { 0, "PBWT", "+",      codec_pbwt_compress,      codec_pbwt_uncompress,    codec_pbwt_reconstruct,    codec_complex_est_size   }, \
    { 1, "ARTB", "+",      codec_ARTB_compress,      codec_arith_uncompress,   NA3,                       codec_ARTB_est_size      }, \
    { 1, "ARTW", "+",      codec_ARTW_compress,      codec_arith_uncompress,   NA3,                       codec_ARTW_est_size      }, /* STRIPE */\
    { 1, "ARTb", "+",      codec_ARTb_compress,      codec_arith_uncompress,   NA3,                       codec_ARTb_est_size      }, /* PACK */\
    { 1, "ARTw", "+",      codec_ARTw_compress,      codec_arith_uncompress,   NA3,                       codec_ARTw_est_size      }, /* STRIPE & PACK */\
    { 0, "BGZF", "+.gz",   NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 0, "XZ",   "+.xz",   NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 0, "BCF",  "-.bcf",  NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 0, "BAM",  "-.bam",  NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 0, "CRAM", "-.cram", NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 0, "ZIP",  "+.zip",  NA1,                      NA2,                      NA3,                       NA4                      }, \
    { 0, "LNGR", "+",      codec_longr_compress,     USE_SUBCODEC,             codec_longr_reconstruct,   codec_longr_est_size     }, \
    { 0, "NRMQ", "+",      codec_normq_compress,     USE_SUBCODEC,             codec_normq_reconstruct,   codec_complex_est_size,  }, \
}

extern CodecArgs codec_args[NUM_CODECS];

extern CodecCompress codec_bz2_compress, codec_lzma_compress, codec_domq_compress, codec_bsc_compress, 
                     codec_none_compress, codec_acgt_compress, codec_xcgt_compress, codec_pbwt_compress, 
                     codec_RANB_compress, codec_RANW_compress, codec_RANb_compress, codec_RANw_compress, 
                     codec_ARTB_compress, codec_ARTW_compress, codec_ARTb_compress, codec_ARTw_compress,
                     codec_longr_compress, codec_normq_compress;

extern CodecUncompress codec_bz2_uncompress, codec_lzma_uncompress, codec_acgt_uncompress, codec_xcgt_uncompress,
                       codec_bsc_uncompress, codec_none_uncompress, codec_gtshark_uncompress, codec_pbwt_uncompress,
                       codec_rans_uncompress, codec_arith_uncompress;

extern CodecReconstruct codec_hapmat_reconstruct, codec_domq_reconstruct, codec_pbwt_reconstruct, 
                        codec_longr_reconstruct, codec_normq_reconstruct;

extern CodecEstSizeFunc codec_none_est_size, codec_bsc_est_size, codec_hapmat_est_size, codec_domq_est_size,
                        codec_RANB_est_size, codec_RANW_est_size, codec_RANb_est_size, codec_RANw_est_size, 
                        codec_ARTB_est_size, codec_ARTW_est_size, codec_ARTb_est_size, codec_ARTw_est_size,
                        codec_complex_est_size, codec_longr_est_size;

// non-codec-specific functions
extern void codec_initialize (void);
extern rom codec_name (Codec codec);
extern void *codec_alloc_do (VBlockP vb, uint64_t size, float grow_at_least_factor, FUNCLINE);
#define codec_alloc(vb,size,grow_at_least_factor) codec_alloc_do((vb),(size),(grow_at_least_factor), __FUNCLINE)

extern void codec_free_do (void *vb, void *addr, FUNCLINE);
#define codec_free(vb,addr) codec_free_do ((vb), (addr), __FUNCLINE)

extern void codec_free_all (VBlockP vb);
extern void codec_verify_free_all (VBlockP vb, rom op, Codec codec);
extern void codec_show_time (VBlockP vb, rom name, rom subname, Codec codec);

#define CODEC_ASSIGN_SAMPLE_SIZE 99999 // bytes (slightly better results than 50K)
extern Codec codec_assign_best_codec (VBlockP vb, ContextP ctx, BufferP non_ctx_data, SectionType st);
extern void codec_assign_best_qual_codec (VBlockP vb, Did qual_did, LocalGetLineCB callback, bool no_longr, bool maybe_revcomped);

#define TAG_NAME (ctx ? ctx->tag_name : "NoContext")

// ACGT stuff
extern const uint8_t acgt_encode[256];
extern void codec_acgt_comp_init (VBlockP vb, Did nonref_did_i);
extern void codec_acgt_reconstruct (VBlockP vb, ContextP ctx, STRp(snip));

// BSC stuff
extern void codec_bsc_initialize (void);

// HAPMAT stuff - retired, used for compressing old files
extern void codec_hapmat_piz_calculate_columns (VBlockP vb);

// DOMQ stuff
extern bool codec_domq_comp_init (VBlockP vb, Did qual_did_i, LocalGetLineCB callback);

// BZ2 stuff
extern uint64_t BZ2_consumed (void *bz_file); // a hacky addition to bzip2

// PBWT stuff
extern void codec_pbwt_seg_init (VBlockP vb, ContextP runs_ctx, ContextP fgrc_ctx);
extern void codec_pbwt_display_ht_matrix (VBlockP vb, uint32_t max_rows);

// LONGR stuff
extern void codec_longr_comp_init (VBlockP vb, Did qual_did_i);
extern void codec_longr_segconf_calculate_bins (VBlockP vb, ContextP ctx, LocalGetLineCB callback);
