// ------------------------------------------------------------------
//   codec.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "lzma/7zTypes.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "data_types.h"

typedef bool CodecCompress (VBlockP vb, 
                            Codec codec,
                            const char *uncompressed,   // option 1 - compress contiguous data
                            uint32_t *uncompressed_len, // in/out (might be modified by complex codecs)
                            LocalGetLineCB callback,    // option 2 - compress data one line at a tim
                            char *compressed, 
                            uint32_t *compressed_len,   // in/out
                            bool soft_fail);
typedef CodecCompress (*Compressor);

extern CodecCompress codec_bz2_compress, codec_lzma_compress, codec_non_acgt_compress, codec_ht_compress, codec_bsc_compress, 
                     codec_domq_compress, codec_none_compress;

typedef void CodecUncompress (VBlockP vb, 
                               const char *compressed, uint32_t compressed_len,
                               char *uncompressed_data, uint64_t uncompressed_len,
                               Codec sub_codec);
typedef CodecUncompress (*Uncompressor);

extern CodecUncompress codec_bz2_uncompress, codec_lzma_uncompress, codec_acgt_uncompress, codec_non_acgt_uncompress,
                       codec_bsc_uncompress, codec_domq_uncompress, codec_ht_uncompress, codec_none_uncompress;

typedef uint32_t CodecEstSizeFunc (uint64_t uncompressed_len);

extern CodecEstSizeFunc codec_none_est_size, codec_bsc_est_size, codec_ht_est_size, codec_domq_est_size;

typedef struct {
    CodecCompress    *compress;
    CodecUncompress  *uncompress;
    CodecEstSizeFunc *est_size;
    Codec            sub_codec; // used if this is a novel codec that invokes another codec
} CodecArgs;

#define CODEC_ARGS { \
    /* NONE */ { codec_none_compress,     codec_none_uncompress,     codec_none_est_size    }, \
    /* GZ   */ { codec_compress_error,    codec_uncompress_error,    codec_est_size_default }, \
    /* BZ2  */ { codec_bz2_compress,      codec_bz2_uncompress,      codec_est_size_default }, \
    /* LZMA */ { codec_lzma_compress,     codec_lzma_uncompress,     codec_est_size_default }, \
    /* BSC  */ { codec_bsc_compress,      codec_bsc_uncompress,      codec_bsc_est_size     }, \
    /* FFU5 */ { codec_compress_error,    codec_uncompress_error,    codec_est_size_default }, \
    /* FFU6 */ { codec_compress_error,    codec_uncompress_error,    codec_est_size_default }, \
    /* FFU7 */ { codec_compress_error,    codec_uncompress_error,    codec_est_size_default }, \
    /* FFU8 */ { codec_compress_error,    codec_uncompress_error,    codec_est_size_default }, \
    /* FFU9 */ { codec_compress_error,    codec_uncompress_error,    codec_est_size_default }, \
    /* ACGT */ { codec_lzma_compress,     codec_acgt_uncompress,     codec_est_size_default }, \
    /* ~CGT */ { codec_non_acgt_compress, codec_non_acgt_uncompress, codec_est_size_default }, \
    /* HT   */ { codec_ht_compress,       codec_ht_uncompress,       codec_ht_est_size,     CODEC_BZ2 }, \
    /* DOMQ */ { codec_domq_compress,     codec_domq_uncompress,     codec_domq_est_size,   CODEC_BSC }, \
}

extern CodecArgs codec_args[NUM_CODECS];

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
extern void codec_acgt_pack (VBlockP vb, const char *data, uint64_t data_len, unsigned bits_consumed, bool do_lten, bool do_lten_partial_final_word);
extern void codec_acgt_pack_last_partial_word (VBlockP vb, ISeqInStream *instream);

// BSC stuff
extern void codec_bsc_initialize (void);

// HT stuff
extern void codec_ht_piz_get_one_line (VBlockP vb);
extern void codec_ht_piz_calculate_columns (VBlockP vb);

// DOMQ stuff
extern bool codec_domq_comp_init (VBlockP vb, DidIType qual_field, LocalGetLineCB callback);
extern void codec_domq_reconstruct (VBlockP vb, ContextP qual_ctx);

// BZ2 stuff
extern uint64_t BZ2_consumed (void *bz_file); // a hacky addition to bzip2

// LZMA stuff
extern void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size);
extern void lzma_free (ISzAllocPtr alloc_stuff, void *addr);
extern const char *lzma_errstr (SRes res);
extern const char *lzma_status (ELzmaStatus status);
