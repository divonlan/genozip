// ------------------------------------------------------------------
//   squeeze_vcf.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SQUEEZE_INCLUDED
#define SQUEEZE_INCLUDED

#include "genozip.h"

extern unsigned squeeze_len (unsigned int len);

extern void squeeze (VBlockVCFP vb,
                     uint8_t *dst, // memory should be pre-allocated by caller
                     uint16_t *squeezed_checksum,
                     const unsigned *src, 
                     unsigned src_len);

extern void unsqueeze (VBlockVCFP vb,
                       unsigned *normal, // memory should be pre-allocated by caller
                       const uint8_t *squeezed, 
                       uint16_t squeezed_checksum,
                       unsigned normal_len);

#endif                       