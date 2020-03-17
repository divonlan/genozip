// ------------------------------------------------------------------
//   gtshark.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GTSHARK_INCLUDED
#define GTSHARK_INCLUDED

#include "genozip.h"
#include "sections.h"

#define GTSHARK_DEFAULT_MEMORY_PER_VB_STR     "128"
#define GTSHARK_DEFAULT_SAMPLES_PER_BLOCK_STR "16384"

extern void gtshark_compress_haplotype_data (VariantBlockP vb, ConstBufferP section_data, unsigned sb_i);

extern void gtshark_uncompress_haplotype_data (VariantBlockP vb, unsigned sb_i);

#endif

