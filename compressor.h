// ------------------------------------------------------------------
//   compresssor.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef COMPRESSOR_INCLUDED
#define COMPRESSOR_INCLUDED

#include "genozip.h"
#include "data_types.h"

extern void comp_compress (VBlockP vb, BufferP z_data, bool is_z_file_buf,
                           SectionHeaderP header, 
                           const char *uncompressed_data, // option 1 - compress contiguous data
                           LocalGetLineCallback callback); // option 2 - compress data one line at a time

extern void comp_uncompress (VBlockP vb, CompressionAlg alg, 
                             const char *compressed_data, uint32_t compressed_data_len,
                             char *uncompressed_data, uint64_t uncompressed_len);

// a hacky addition to bzip2
extern uint64_t BZ2_consumed (void *bz_file);

typedef bool CompressorFunc (VBlockP vb, 
                             CompressionAlg alg,
                             const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                             LocalGetLineCallback callback,                        // option 2 - compress data one line at a tim
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail);
typedef CompressorFunc (*Compressor);

CompressorFunc comp_compress_bzlib, comp_compress_lzma, comp_compress_none;

extern const uint8_t acgt_encode[256];
extern const char acgt_decode[4];
#define ACGT_DECODE(bitarr,idx) acgt_decode[bit_array_get ((bitarr), ((int64_t)(idx))*2) + (bit_array_get ((bitarr), ((int64_t)(idx))*2 + 1) << 1)]
#endif