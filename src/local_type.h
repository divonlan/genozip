// ------------------------------------------------------------------
//   local_type.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

// LT_* values are consistent with BAM optional 'B' types (and extend them)
typedef packed_enum { // 1 byte
    // LT values that are part of the file format - values can be added but not changed
    LT_SINGLETON  = 0,   // nul-terminated singleton snips (note: ltype=0 was called LT_TEXT until 15.0.26 and included both SINGLETONs and STRINGs)
    LT_INT8       = 1,    
    LT_UINT8      = 2,
    LT_INT16      = 3,
    LT_UINT16     = 4,
    LT_INT32      = 5,
    LT_UINT32     = 6,
    LT_INT64      = 7,   
    LT_UINT64     = 8,   
    LT_FLOAT32    = 9,   
    LT_FLOAT64    = 10,  
    LT_BLOB       = 11,  // length of data extracted is determined by vb->seq_len or provided in the LOOKUP snip (until 15.0.26 called LT_SEQUENCE)
    LT_BITMAP     = 12,  // a bitmap
    LT_CODEC      = 13,  // codec specific type with its codec specific reconstructor
    LT_UINT8_TR   = 14,  // transposed array - number of columns in original array is in param (up to 255 columns)
    LT_UINT16_TR  = 15,  // "
    LT_UINT32_TR  = 16,  // "
    LT_UINT64_TR  = 17,  // "
    LT_hex8       = 18,  // lower-case UINT8  hex
    LT_HEX8       = 19,  // upper-case UINT8  hex
    LT_hex16      = 20,  // lower-case UINT16 hex
    LT_HEX16      = 21,  // upper-case UINT16 hex
    LT_hex32      = 22,  // lower-case UINT32 hex
    LT_HEX32      = 23,  // upper-case UINT32 hex
    LT_hex64      = 24,  // lower-case UINT64 hex
    LT_HEX64      = 25,  // upper-case UINT64 hex
    LT_STRING     = 26,  // nul-terminated strings
    LT_SUPP       = 27,  // supplementary data used in reconstruction of another context. Data is not BGENed etc, and context is never reconstructed directly (15.0.40)
    LT_UINT8_PTR  = 28,  // partial transposed array - only items not-copied by VCF_COPY_SAMPLE are included
    LT_UINT16_PTR = 29,  // "
    LT_UINT32_PTR = 30,  // "
    NUM_LTYPES,          // counts LocalTypes that can appear in the Genozip file format

    // LT_DYN* - LT values that are NOT part of the file format, just used during seg 
    LT_DYN_INT,         // dynamic size local 
    LT_DYN_INT_h,       // dynamic size local - hex
    LT_DYN_INT_H,       // dynamic size local - HEX
    
    NUM_LOCAL_TYPES     // counts all LocalTypes
} LocalType;

#define IS_LT_DYN(ltype) ((ltype) == LT_DYN_INT || (ltype) == LT_DYN_INT_h || (ltype) == LT_DYN_INT_H)

typedef void BgEnBufFunc (BufferP buf, LocalType *lt); 

typedef BgEnBufFunc (*BgEnBuf);

typedef enum  {  BAM_NA=0, BAM_c, BAM_C, BAM_s, BAM_S, BAM_i, BAM_I, NUM_BAM_INT_TYPES } BamIntTypes;
#define BAM_INT_TYPES { 0,    'c',   'C',   's',   'S',   'i',   'I' }
typedef struct LocalTypeDesc { // 16 bytes
    Pointeר nameר           : 28; // 28 bit is plenty for a relative pointer within our static data
    uint32_t bam_type       : 3;
    uint32_t is_signed      : 1;

    uint32_t width          : 4;
    Pointeר file_to_nativeר : 28; // relative pointer to function

    int64_t max_int; // relevant for integer fields only. if is_signed, min_int = (-max_int-1)
} LocalTypeDesc;

extern LocalTypeDesc *lt_desc;

#define LOCALTYPE_DESC {                                                                                                          \
/*  name          bam_type  signed width file_to_native              max_int   */                                                 \
   { ר("SIN"),    0,        0,     1,    0,                          0          },                                                \
   /* 64B-alignment starts here: each entry is 16B ⇒ the 8 I/U integers fit in 2 cache lines */                                   \
   { ר("I8 "),    BAM_c,    1,     1,    ר(BGEN_deinterlace_d8_buf), INT8_MAX   },                                                \
   { ר("U8 "),    BAM_C,    0,     1,    ר(BGEN_u8_buf),             UINT8_MAX  },                                                \
   { ר("I16"),    BAM_s,    1,     2,    ר(BGEN_deinterlace_d16_buf),INT16_MAX  },                                                \
   { ר("U16"),    BAM_S,    0,     2,    ר(BGEN_u16_buf),            UINT16_MAX },                                                \
   { ר("I32"),    BAM_i,    1,     4,    ר(BGEN_deinterlace_d32_buf),INT32_MAX  },                                                \
   { ר("U32"),    BAM_I,    0,     4,    ר(BGEN_u32_buf),            UINT32_MAX },                                                \
   { ר("I64"),    0,        1,     8,    ר(BGEN_deinterlace_d64_buf),INT64_MAX  },                                                \
   { ר("U64"),    0,        0,     8,    ר(BGEN_u64_buf),            INT64_MAX  }, /* internal rep is int64_t so max is limited */\
   { ר("F32"),    0,        0,     4,    ר(BGEN_u32_buf),            0          },                                                \
   { ר("F64"),    0,        0,     8,    ר(BGEN_u64_buf),            0          },                                                \
   { ר("BLB"),    0,        0,     1,    0,                          0          },                                                \
   { ר("BMP"),    0,        0,     8,    0,                          0          },                                                \
   { ר("COD"),    0,        0,     1,    0,                          0          },                                                \
   { ר("T8 "),    0,        0,     1,    ר(BGEN_transpose_u8_buf),   UINT8_MAX  },                                                \
   { ר("T16"),    0,        0,     2,    ר(BGEN_transpose_u16_buf),  UINT16_MAX },                                                \
   { ר("T32"),    0,        0,     4,    ר(BGEN_transpose_u32_buf),  UINT32_MAX },                                                \
   { ר("N/A"),    0,        0,     8,    0,                          INT64_MAX  }, /* unused */                                   \
   { ר("h8 "),    0,        0,     1,    ר(BGEN_u8_buf),             UINT8_MAX  }, /* lower-case UINT8 hex */                     \
   { ר("H8 "),    0,        0,     1,    ר(BGEN_u8_buf),             UINT8_MAX  }, /* upper-case UINT8 hex */                     \
   { ר("h16"),    0,        0,     2,    ר(BGEN_u16_buf),            UINT16_MAX },                                                \
   { ר("H16"),    0,        0,     2,    ר(BGEN_u16_buf),            UINT16_MAX },                                                \
   { ר("h32"),    0,        0,     4,    ר(BGEN_u32_buf),            UINT32_MAX },                                                \
   { ר("H32"),    0,        0,     4,    ר(BGEN_u32_buf),            UINT32_MAX },                                                \
   { ר("h64"),    0,        0,     8,    ר(BGEN_u64_buf),            INT64_MAX  },                                                \
   { ר("H64"),    0,        0,     8,    ר(BGEN_u64_buf),            INT64_MAX  },                                                \
   { ר("STR"),    0,        0,     1,    0,                          0          },                                                \
   { ר("SUP"),    0,        0,     1,    0,                          0          },                                                \
   { ר("t8 "),    0,        0,     1,    ר(BGEN_ptranspose_u8_buf),  UINT8_MAX  },                                                \
   { ר("t16"),    0,        0,     2,    ר(BGEN_ptranspose_u16_buf), UINT16_MAX },                                                \
   { ר("t32"),    0,        0,     4,    ר(BGEN_ptranspose_u32_buf), UINT32_MAX },                                                \
   { /* NUM_LTYPES */                                                           },                                                \
   /* from here - not part of the file format, just used during seg */                                                            \
   { ר("DYN"),    0,        0,     8,    0,                          INT64_MAX  },                                                \
   { ר("DYh"),    0,        0,     8,    0,                          INT64_MAX  },                                                \
   { ר("DYH"),    0,        0,     8,    0,                          INT64_MAX  },                                                \
}

#define lt_width(ctx)       (lt_desc[(ctx)->ltype].width)
#define lt_is_signed(ltype) (lt_desc[ltype].is_signed)
#define lt_max(ltype)       (lt_desc[ltype].max_int)
#define lt_min(ltype)       (lt_is_signed(ltype) ? (-lt_max(ltype) - 1) : 0)
