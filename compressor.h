// ------------------------------------------------------------------
//   zfile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef COMPRESSOR_INCLUDED
#define COMPRESSOR_INCLUDED

// IMPORTANT: these values CANNOT BE CHANGED as they are part of the genozip file - 
// they go in SectionHeader.sec_compression_alg and also SectionHeaderTxtHeader.compression_type
#define NUM_COMPRESSION_ALGS 11 // update compressors array in comp_compress if making any change 
typedef enum { 
    COMP_UNKNOWN=-1, COMP_NONE=0 /* plain - no compression */,            
    COMP_GZ=1, COMP_BZ2=2, COMP_BGZ=3, COMP_XZ=4, COMP_BCF=5, COMP_BAM=6, COMP_LZMA=7, COMP_ZIP=8, 
    // compress a sequence of A,C,G,T nucleotides - first squeeze into 2 bits and then LZMA. It's about 25X faster and 
    // slightly better compression ratio than LZMA. Any characters that are not ACGT are stored in a complimentary 
    // COMP_NON_ACGT compression - which is \0 for ACGT locations, \1 for acgt locations and verbatim for other characters
    COMP_ACGT=9, COMP_NON_ACGT=10 
} CompressionAlg; 

#define COMPRESSED_FILE_VIEWER { "cat", "gunzip -d -c", "bzip2 -d -c", "gunzip -d -c", "xz -d -c", \
                                 "bcftools -Ov --version", "samtools view -h -OSAM --threads 2", "N/A", "unzip -p", "N/A" }

#include "genozip.h"

typedef void CompGetLineCallback (VBlockP vb, uint32_t vb_line_i, 
                                  char **line_data_1, uint32_t *line_data_len_1,
                                  char **line_data_2, uint32_t *line_data_len_2);

extern void comp_compress (VBlockP vb, BufferP z_data, bool is_z_file_buf,
                           SectionHeaderP header, 
                           const char *uncompressed_data, // option 1 - compress contiguous data
                           CompGetLineCallback callback); // option 2 - compress data one line at a time

extern void comp_uncompress (VBlockP vb, CompressionAlg alg, 
                             const char *compressed_data, uint32_t compressed_data_len,
                             char *uncompressed_data, uint64_t uncompressed_len);

// a hacky addition to bzip2
extern uint64_t BZ2_consumed (void *bz_file);

typedef bool CompressorFunc (VBlockP vb, 
                             CompressionAlg alg,
                             const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                             CompGetLineCallback callback,                        // option 2 - compress data one line at a tim
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail);
typedef CompressorFunc (*Compressor);

CompressorFunc comp_compress_bzlib, comp_compress_lzma, comp_compress_none;

extern const char actg_decode[4];
extern const uint8_t actg_encode[256];

#endif