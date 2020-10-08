// ------------------------------------------------------------------
//   compresssor.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef COMPRESSOR_INCLUDED
#define COMPRESSOR_INCLUDED

#include "genozip.h"
#include "data_types.h"

extern void comp_initialize (void);

extern void comp_compress (VBlockP vb, BufferP z_data, bool is_z_file_buf,
                           SectionHeaderP header, 
                           const char *uncompressed_data, // option 1 - compress contiguous data
                           LocalGetLineCallback callback); // option 2 - compress data one line at a time

extern void comp_uncompress (VBlockP vb, Codec codec, 
                             const char *compressed_data, uint32_t compressed_data_len,
                             char *uncompressed_data, uint64_t uncompressed_len);

extern const char *codec_name (Codec codec);

// a hacky addition to bzip2
extern uint64_t BZ2_consumed (void *bz_file);

extern const uint8_t acgt_encode[256];
extern const char acgt_decode[4];
#define ACGT_DECODE(bitarr,idx) acgt_decode[bit_array_get ((bitarr), ((int64_t)(idx))*2) + (bit_array_get ((bitarr), ((int64_t)(idx))*2 + 1) << 1)]

extern void comp_ht_piz_get_one_line (VBlockP vb);
extern void comp_ht_piz_calculate_columns (VBlockP vb);

#endif