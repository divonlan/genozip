// ------------------------------------------------------------------
//   codec.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef CODEC_INCLUDED
#define CODEC_INCLUDED

#include "lzma/7zTypes.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "data_types.h"

typedef bool CodecCompress (VBlockP vb, 
                            Codec *codec,               // in / out
                            const char *uncompressed,   // option 1 - compress contiguous data
                            uint32_t *uncompressed_len, // in/out (might be modified by complex codecs)
                            LocalGetLineCB callback,    // option 2 - compress data one line at a tim
                            char *compressed, 
                            uint32_t *compressed_len,   // in/out
                            bool soft_fail);

typedef void CodecUncompress (VBlockP vb, Codec codec,
                              const char *compressed, uint32_t compressed_len,
                              Buffer *uncompressed_buf, uint64_t uncompressed_len,
                              Codec sub_codec);

typedef uint32_t CodecEstSizeFunc (Codec codec, uint64_t uncompressed_len);

typedef void CodecReconstruct (VBlockP vb, Codec codec, ContextP ctx);

typedef struct {
    const char       *name;
    const char       *ext; // extensions by compression type. + if it adds to the name ; - if it replaces the extension of the uncompress name
    CodecCompress    *compress;
    CodecUncompress  *uncompress;
    CodecReconstruct *reconstruct;
    CodecEstSizeFunc *est_size;
    Codec            sub_codec1, sub_codec2;  // for complex codecs that invokes another codec (for each of the two sections)
} CodecArgs;

#define NA1 codec_compress_error
#define NA2 codec_uncompress_error
#define NA3 codec_reconstruct_error
#define NA4 codec_est_size_default
#define USE_SUBCODEC NULL
#define CODEC_ARGS { /* aligned with Codec defined in genozip.h */ \
    { "NONE", "+",      codec_none_compress, codec_none_uncompress, NA3,                    codec_none_est_size }, \
    { "GZ",   "+.gz",   NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "BZ2",  "+.bz",   codec_bz2_compress,  codec_bz2_uncompress,  NA3,                    NA4                 }, \
    { "LZMA", "+",      codec_lzma_compress, codec_lzma_uncompress, NA3,                    NA4                 }, \
    { "BSC",  "+",      codec_bsc_compress,  codec_bsc_uncompress,  NA3,                    codec_bsc_est_size  }, \
    { "FFU5", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FFU6", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FFU7", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FFU8", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FFU9", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "ACGT", "+",      codec_acgt_compress, codec_acgt_uncompress, NA3,                    codec_sc1_est_size, CODEC_LZMA /* NONREF   */ }, \
    { "XCGT", "+",      codec_xcgt_compress, codec_xcgt_uncompress, NA3,                    codec_sc1_est_size, CODEC_BZ2  /* NONREF_X */ }, \
    { "HT",   "+",      codec_ht_compress,   USE_SUBCODEC,          codec_ht_reconstruct,   codec_sc1_est_size, CODEC_BZ2  /* GT_HT    */, CODEC_BZ2 /* GT_HT_INDEX */ }, \
    { "DOMQ", "+",      codec_domq_compress, USE_SUBCODEC,          codec_domq_reconstruct, codec_sc1_est_size, CODEC_BSC  /* QUAL     */, CODEC_BSC /* DOMQRUNS    */ }, \
    { "FF14", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FF15", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FF16", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FF17", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FF18", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "FF19", "+",      NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "BGZ",  "+.bgz",  NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "XZ",   "+.xz",   NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "BCF",  "-.bcf",  NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "BAM",  "-.bam",  NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "CRAM", "-.cram", NA1,                 NA2,                   NA3,                    NA4                 }, \
    { "ZIP",  "+.zip",  NA1,                 NA2,                   NA3,                    NA4                 }, \
}

extern CodecArgs codec_args[NUM_CODECS];

extern CodecCompress codec_bz2_compress, codec_lzma_compress, codec_domq_compress, codec_ht_compress, codec_bsc_compress, 
                     codec_none_compress, codec_acgt_compress, codec_xcgt_compress;

extern CodecUncompress codec_bz2_uncompress, codec_lzma_uncompress, codec_acgt_uncompress, codec_xcgt_uncompress,
                       codec_bsc_uncompress, codec_none_uncompress;

extern CodecReconstruct codec_ht_reconstruct, codec_domq_reconstruct;

extern CodecEstSizeFunc codec_none_est_size, codec_bsc_est_size, codec_ht_est_size, codec_domq_est_size;

// non-codec-specific functions
extern void codec_initialize (void);
extern const char *codec_name (Codec codec);
extern void *codec_alloc (VBlockP vb, int size, double grow_at_least_factor);
extern void codec_free (VBlockP vb, void *addr);
extern void codec_free_all (VBlockP vb);

// ACGT stuff
extern const uint8_t acgt_encode[256];
extern const char acgt_decode[4];
#define ACGT_DECODE(bitarr,idx) acgt_decode[bit_array_get ((bitarr), ((int64_t)(idx))*2) + (bit_array_get ((bitarr), ((int64_t)(idx))*2 + 1) << 1)]

extern void codec_acgt_comp_init (VBlockP vb);
extern void codec_acgt_reconstruct (VBlockP vb, ContextP ctx, const char *snip, unsigned snip_len);

// BSC stuff
extern void codec_bsc_initialize (void);

// HT stuff
extern void codec_ht_comp_init (VBlockP vb);
extern void codec_ht_piz_calculate_columns (VBlockP vb);

// DOMQ stuff
extern bool codec_domq_comp_init (VBlockP vb, LocalGetLineCB callback);

// BZ2 stuff
extern uint64_t BZ2_consumed (void *bz_file); // a hacky addition to bzip2

// LZMA stuff
extern void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size);
extern void lzma_free (ISzAllocPtr alloc_stuff, void *addr);
extern const char *lzma_errstr (SRes res);
extern const char *lzma_status (ELzmaStatus status);

#endif
