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

extern CompressorFunc comp_compress_bzlib, comp_compress_lzma, comp_compress_non_acgt, comp_compress_ht;

typedef void UncompressorFunc (VBlockP vb, 
                               const char *compressed, uint32_t compressed_len,
                               char *uncompressed_data, uint64_t uncompressed_len);

UncompressorFunc comp_uncompress_bzlib, comp_uncompress_lzma, comp_uncompress_acgt, comp_uncompress_non_acgt;
typedef UncompressorFunc (*Uncompressor);

extern void *comp_alloc (VBlockP vb, int size, double grow_at_least_factor);
extern void comp_free (VBlockP vb, void *addr);
extern void comp_free_all (VBlockP vb);

extern void *lzma_alloc (ISzAllocPtr alloc_stuff, size_t size);
extern void lzma_free (ISzAllocPtr alloc_stuff, void *addr);
const char *lzma_errstr (SRes res);
const char *lzma_status (ELzmaStatus status);

extern void comp_acgt_pack (VBlockP vb, const char *data, uint64_t data_len, unsigned bits_consumed, bool do_lten, bool do_lten_partial_final_word);
extern void comp_acgt_pack_last_partial_word (VBlockP vb, ISeqInStream *instream);
