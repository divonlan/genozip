// ------------------------------------------------------------------
//   comp_private.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "lzma/7zTypes.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "compressor.h"

typedef bool CompressorFunc (VBlockP vb, 
                             Codec codec,
                             const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                             LocalGetLineCallback callback,                        // option 2 - compress data one line at a tim
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail);
typedef CompressorFunc (*Compressor);

extern CompressorFunc comp_bzlib_compress, comp_lzma_compress, comp_non_acgt_compress, comp_ht_compress, comp_bsc_compress, 
                      comp_none_compress;

typedef void UncompressorFunc (VBlockP vb, 
                               const char *compressed, uint32_t compressed_len,
                               char *uncompressed_data, uint64_t uncompressed_len);

extern UncompressorFunc comp_bzlib_uncompress, comp_lzma_uncompress, comp_acgt_uncompress, comp_non_acgt_uncompress,
                        comp_bsc_uncompress, comp_ht_uncompress, comp_none_uncompress;

typedef uint32_t EstSizeFunc (uint64_t uncompressed_len);

extern EstSizeFunc comp_none_est_size, comp_bsc_est_size, comp_ht_est_size;

typedef struct {
    CompressorFunc *compress;
    UncompressorFunc *uncompress;
    EstSizeFunc *est_size;
} CodecArgs;

#define CODEC_ARGS { \
    /* NONE */ { comp_none_compress,     comp_none_uncompress,     comp_none_est_size    }, \
    /* GZ   */ { comp_compress_error,    comp_uncompress_error,    comp_est_size_default }, \
    /* BZ2  */ { comp_bzlib_compress,    comp_bzlib_uncompress,    comp_est_size_default }, \
    /* LZMA */ { comp_lzma_compress,     comp_lzma_uncompress,     comp_est_size_default }, \
    /* BSC  */ { comp_bsc_compress,      comp_bsc_uncompress,      comp_bsc_est_size     }, \
    /* FFU5 */ { comp_compress_error,    comp_uncompress_error,    comp_est_size_default }, \
    /* FFU6 */ { comp_compress_error,    comp_uncompress_error,    comp_est_size_default }, \
    /* FFU7 */ { comp_compress_error,    comp_uncompress_error,    comp_est_size_default }, \
    /* FFU8 */ { comp_compress_error,    comp_uncompress_error,    comp_est_size_default }, \
    /* FFU9 */ { comp_compress_error,    comp_uncompress_error,    comp_est_size_default }, \
    /* ACGT */ { comp_lzma_compress,     comp_acgt_uncompress,     comp_est_size_default }, \
    /* ~CGT */ { comp_non_acgt_compress, comp_non_acgt_uncompress, comp_est_size_default }, \
    /* HT   */ { comp_ht_compress,       comp_ht_uncompress,       comp_ht_est_size      }, \
}

extern CodecArgs codec_args[NUM_CODECS];

typedef UncompressorFunc (*Uncompressor);

extern void *comp_alloc (VBlockP vb, int size, double grow_at_least_factor);
extern void comp_free (VBlockP vb, void *addr);
extern void comp_free_all (VBlockP vb);

extern void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size);
extern void lzma_free (ISzAllocPtr alloc_stuff, void *addr);
extern const char *lzma_errstr (SRes res);
extern const char *lzma_status (ELzmaStatus status);

extern void comp_acgt_pack (VBlockP vb, const char *data, uint64_t data_len, unsigned bits_consumed, bool do_lten, bool do_lten_partial_final_word);
extern void comp_acgt_pack_last_partial_word (VBlockP vb, ISeqInStream *instream);

extern void comp_bsc_initialize (void);
