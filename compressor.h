// ------------------------------------------------------------------
//   zfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef COMPRESSOR_INCLUDED
#define COMPRESSOR_INCLUDED

#include "genozip.h"

#define NUM_COMPRESSOR_ALGS 3
typedef enum { COMPRESS_UNKNOWN=-1, // value used internally, never written to disk
               COMPRESS_NONE=0, COMPRESS_BZLIB=1, COMPRESS_LZMA=2 } CompressorAlg; // goes into SectionHeader.data_compression_alg

typedef void CompGetLineCallback (VBlockP vb, uint32_t vb_line_i, 
                                  char **line_data_1, uint32_t *line_data_len_1,
                                  char **line_data_2, uint32_t *line_data_len_2);

extern void comp_compress (VBlockP vb, BufferP z_data, bool is_z_file_buf,
                           SectionHeaderP header, 
                           const char *uncompressed_data, // option 1 - compress contiguous data
                           CompGetLineCallback callback); // option 2 - compress data one line at a time

extern void comp_uncompress (VBlockP vb, CompressorAlg alg, 
                             const char *compressed_data, uint32_t compressed_data_len,
                             BufferP uncompressed_data);

// a hacky addition to bzip2
extern uint64_t BZ2_consumed (void *bz_file);

typedef bool CompressorFunc (VBlockP vb, 
                             const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                             CompGetLineCallback callback,                        // option 2 - compress data one line at a tim
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail);
typedef CompressorFunc (*Compressor);

CompressorFunc comp_compress_bzlib, comp_compress_lzma, comp_compress_none;

#endif