// ------------------------------------------------------------------
//   codec.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "lzma/7zTypes.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "data_types.h"

#define MIN_LEN_FOR_COMPRESSION 90 // less that this size, and compressed size is typically larger than uncompressed size

#define ENANO_MIN_READ_LEN 4

typedef bool CodecCompress (VBlockP vb, 
                            SectionHeader *header,       // in / out
                            const char *uncompressed,   // option 1 - compress contiguous data
                            uint32_t *uncompressed_len, // in/out (might be modified by complex codecs)
                            LocalGetLineCB callback,    // option 2 - compress data one line at a tim
                            char *compressed, 
                            uint32_t *compressed_len,   // in/out
                            bool soft_fail);

typedef void CodecUncompress (VBlockP vb, Codec codec, uint8_t param,
                              const char *compressed, uint32_t compressed_len,
                              Buffer *uncompressed_buf, uint64_t uncompressed_len,
                              Codec sub_codec);

typedef uint32_t CodecEstSizeFunc (Codec codec, uint64_t uncompressed_len);

typedef void CodecReconstruct (VBlockP vb, Codec codec, ContextP ctx);

typedef struct {
    bool             is_simple;  // a simple codec is one that is compressed into a single section in one step
    const char       *name;
    const char       *ext;       // extensions by compression type. + if it adds to the name ; - if it replaces the extension of the uncompress name
    const char       *viewer;    // command line to view a file of this type
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
/*  simp name    ext       viewer         compress                  uncompress             reconstruct                est_size                 */ \
    { 1, "N/A",  "+",      NA0,           NA1,                      NA2,                   NA3,                       NA4                      }, \
    { 1, "NONE", "+",      "cat",         codec_none_compress,      codec_none_uncompress, NA3,                       codec_none_est_size      }, \
    { 1, "GZ",   "+.gz",   "gunzip -c",   NA1,                      NA2,                   NA3,                       NA4                      }, \
    { 1, "BZ2",  "+.bz2",  "bzip2 -d -c", codec_bz2_compress,       codec_bz2_uncompress,  NA3,                       NA4                      }, \
    { 1, "LZMA", "+",      NA0,           codec_lzma_compress,      codec_lzma_uncompress, NA3,                       NA4                      }, \
    { 1, "BSC",  "+",      NA0,           codec_bsc_compress,       codec_bsc_uncompress,  NA3,                       codec_bsc_est_size       }, \
    { 1, "RANB", "+",      NA0,           codec_RANB_compress,      codec_rans_uncompress, NA3,                       codec_RANB_est_size      }, \
    { 1, "RANW", "+",      NA0,           codec_RANW_compress,      codec_rans_uncompress, NA3,                       codec_RANW_est_size      }, \
    { 1, "RANb", "+",      NA0,           codec_RANb_compress,      codec_rans_uncompress, NA3,                       codec_RANb_est_size      }, \
    { 1, "RANw", "+",      NA0,           codec_RANw_compress,      codec_rans_uncompress, NA3,                       codec_RANw_est_size      }, \
    { 0, "ACGT", "+",      NA0,           codec_acgt_compress,      codec_acgt_uncompress, NA3,                       codec_complex_est_size   }, \
    { 0, "XCGT", "+",      NA0,           USE_SUBCODEC,             codec_xcgt_uncompress, NA3,                       NA4                      }, \
    { 0, "HAPM", "+",      NA0,           NA1,                      USE_SUBCODEC,          codec_hapmat_reconstruct,  NA4,                     }, /* HapMat was discontinued and replaced by PBWT. We keep it for decompressing old VCF files */ \
    { 0, "DOMQ", "+",      NA0,           codec_domq_compress,      USE_SUBCODEC,          codec_domq_reconstruct,    codec_complex_est_size,  }, \
    { 0, "GTSH", "+",      NA0,           NA1,                      codec_gtshark_uncompress, codec_pbwt_reconstruct, NA4,                     }, /* gtshark discontinued in v12. keep for displaying an error */\
    { 0, "PBWT", "+",      NA0,           codec_pbwt_compress,      codec_pbwt_uncompress, codec_pbwt_reconstruct,    codec_complex_est_size   }, \
    { 1, "ARTB", "+",      NA0,           codec_ARTB_compress,      codec_arith_uncompress,NA3,                       codec_ARTB_est_size      }, \
    { 1, "ARTW", "+",      NA0,           codec_ARTW_compress,      codec_arith_uncompress,NA3,                       codec_ARTW_est_size      }, \
    { 1, "ARTb", "+",      NA0,           codec_ARTb_compress,      codec_arith_uncompress,NA3,                       codec_ARTb_est_size      }, \
    { 1, "ARTw", "+",      NA0,           codec_ARTw_compress,      codec_arith_uncompress,NA3,                       codec_ARTw_est_size      }, \
    { 0, "BGZF", "+.gz",   "gunzip -c",   NA1,                      NA2,                   NA3,                       NA4                      }, \
    { 0, "XZ",   "+.xz",   "xz -d -c",    NA1,                      NA2,                   NA3,                       NA4                      }, \
    { 0, "BCF",  "-.bcf",  "bcftools view", NA1,                    NA2,                   NA3,                       NA4                      }, \
    { 0, "BAM",  "-.bam",  "samtools view -h --threads 2", NA1,     NA2,                   NA3,                       NA4                      }, \
    { 0, "CRAM", "-.cram", "samtools view -h --threads 2", NA1,     NA2,                   NA3,                       NA4                      }, \
    { 0, "ZIP",  "+.zip",  "unzip -p",    NA1,                      NA2,                   NA3,                       NA4                      }, \
    { 0, "NANO", "+",      NA0,           codec_enano_compress,     codec_enano_uncompress,NA3,                       codec_complex_est_size   }, \
}

extern CodecArgs codec_args[NUM_CODECS];

extern CodecCompress codec_bz2_compress, codec_lzma_compress, codec_domq_compress, codec_bsc_compress, 
                     codec_none_compress, codec_acgt_compress, codec_xcgt_compress, codec_pbwt_compress, 
                     codec_RANB_compress, codec_RANW_compress, codec_RANb_compress, codec_RANw_compress, 
                     codec_ARTB_compress, codec_ARTW_compress, codec_ARTb_compress, codec_ARTw_compress,
                     codec_enano_compress;

extern CodecUncompress codec_bz2_uncompress, codec_lzma_uncompress, codec_acgt_uncompress, codec_xcgt_uncompress,
                       codec_bsc_uncompress, codec_none_uncompress, codec_gtshark_uncompress, codec_pbwt_uncompress,
                       codec_rans_uncompress, codec_arith_uncompress, codec_enano_uncompress;

extern CodecReconstruct codec_hapmat_reconstruct, codec_domq_reconstruct, codec_pbwt_reconstruct;

extern CodecEstSizeFunc codec_none_est_size, codec_bsc_est_size, codec_hapmat_est_size, codec_domq_est_size,
                        codec_RANB_est_size, codec_RANW_est_size, codec_RANb_est_size, codec_RANw_est_size, 
                        codec_ARTB_est_size, codec_ARTW_est_size, codec_ARTb_est_size, codec_ARTw_est_size,
                        codec_complex_est_size;

// non-codec-specific functions
extern void codec_initialize (void);
extern const char *codec_name (Codec codec);
extern void *codec_alloc (VBlockP vb, uint64_t size, float grow_at_least_factor);
extern void codec_free (void *vb, void *addr);
extern void codec_free_all (VBlockP vb);
extern void codec_verify_free_all (VBlockP vb, const char *op, Codec codec);
extern void codec_show_time (VBlockP vb, const char *name, const char *subname, Codec codec);

#define CODEC_ASSIGN_SAMPLE_SIZE 99999 // bytes (slightly better results than 50K)
extern Codec codec_assign_best_codec (VBlockP vb, ContextP ctx, BufferP non_ctx_data, SectionType st);

// ACGT stuff
extern const uint8_t acgt_encode[256];
extern const char acgt_decode[4];
#define ACGT_DECODE(bitarr,idx) acgt_decode[bit_array_get ((bitarr), ((int64_t)(idx))*2) + (bit_array_get ((bitarr), ((int64_t)(idx))*2 + 1) << 1)]

extern void codec_acgt_comp_init (VBlockP vb);
extern void codec_acgt_reconstruct (VBlockP vb, ContextP ctx, STRp(snip));

// BSC stuff
extern void codec_bsc_initialize (void);

// HAPMAT stuff - retired, used for compressing old files
extern void codec_hapmat_piz_calculate_columns (VBlockP vb);

// DOMQ stuff
extern bool codec_domq_comp_init (VBlockP vb, DidIType qual_did_i, LocalGetLineCB callback);

// BZ2 stuff
extern uint64_t BZ2_consumed (void *bz_file); // a hacky addition to bzip2

// LZMA stuff
extern void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size);
extern void lzma_free (ISzAllocPtr alloc_stuff, void *addr);
extern const char *lzma_errstr (SRes res);
extern const char *lzma_status (ELzmaStatus status);

// PBWT stuff
extern void codec_pbwt_seg_init (VBlockP vb, ContextP runs_ctx, ContextP fgrc_ctx);
extern void codec_pbwt_display_ht_matrix (VBlockP vb, uint32_t max_rows);

// NANO stuff
extern void codec_enano_seg_init (VBlockP vb, ContextP qual_ctx);
