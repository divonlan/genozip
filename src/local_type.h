// ------------------------------------------------------------------
//   local_type.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

// LT_* values are consistent with BAM optional 'B' types (and extend them)
typedef enum __attribute__ ((packed)) { // 1 byte
    // LT values that are part of the file format - values can be added but not changed
    LT_SINGLETON = 0,   // nul-terminated singleton snips (note: ltype=0 was called LT_TEXT until 15.0.26 and included both SINGLETONs and STRINGs)
    LT_INT8      = 1,    
    LT_UINT8     = 2,
    LT_INT16     = 3,
    LT_UINT16    = 4,
    LT_INT32     = 5,
    LT_UINT32    = 6,
    LT_INT64     = 7,   // ffu
    LT_UINT64    = 8,   // ffu
    LT_FLOAT32   = 9,   
    LT_FLOAT64   = 10,  // ffu
    LT_BLOB      = 11,  // length of data extracted is determined by vb->seq_len or provided in the LOOKUP snip (until 15.0.26 called LT_SEQUENCE)
    LT_BITMAP    = 12,  // a bitmap
    LT_CODEC     = 13,  // codec specific type with its codec specific reconstructor
    LT_UINT8_TR  = 14,  // transposed array - number of columns in original array is in param (up to 255 columns)
    LT_UINT16_TR = 15,  // "
    LT_UINT32_TR = 16,  // "
    LT_UINT64_TR = 17,  // "
    LT_hex8      = 18,  // lower-case UINT8  hex
    LT_HEX8      = 19,  // upper-case UINT8  hex
    LT_hex16     = 20,  // lower-case UINT16 hex
    LT_HEX16     = 21,  // upper-case UINT16 hex
    LT_hex32     = 22,  // lower-case UINT32 hex
    LT_HEX32     = 23,  // upper-case UINT32 hex
    LT_hex64     = 24,  // lower-case UINT64 hex
    LT_HEX64     = 25,  // upper-case UINT64 hex
    LT_STRING    = 26,  // nul-terminated strings
    NUM_LTYPES,         // counts LocalTypes that can appear in the Genozip file format

    // LT_DYN* - LT values that are NOT part of the file format, just used during seg 
    LT_DYN_INT,         // dynamic size local 
    LT_DYN_INT_h,       // dynamic size local - hex
    LT_DYN_INT_H,       // dynamic size local - HEX
    
    NUM_LOCAL_TYPES     // counts all LocalTypes
} LocalType;

#define IS_LT_DYN(ltype) ((ltype) == LT_DYN_INT || (ltype) == LT_DYN_INT_h || (ltype) == LT_DYN_INT_H)

typedef void BgEnBufFunc (BufferP buf, LocalType *lt); 

typedef BgEnBufFunc (*BgEnBuf);

typedef struct LocalTypeDesc {
    rom name;
    const char sam_type;
    unsigned width;
    bool is_signed;
    int64_t min_int, max_int; // relevant for integer fields only
    BgEnBuf file_to_native;
} LocalTypeDesc;

extern const LocalTypeDesc lt_desc[NUM_LOCAL_TYPES];
#define LOCALTYPE_DESC {                                                                              \
/*   name   sam  wid signed min_int                max_int                file_to_native           */ \
   { "SIN", 0,   1,  0,     0,                     0,                     0                        }, \
   { "I8 ", 'c', 1,  1,     -0x80LL,               0x7fLL,                BGEN_deinterlace_d8_buf  }, \
   { "U8 ", 'C', 1,  0,     0,                     0xffLL,                BGEN_u8_buf              }, \
   { "I16", 's', 2,  1,     -0x8000LL,             0x7fffLL,              BGEN_deinterlace_d16_buf }, \
   { "U16", 'S', 2,  0,     0,                     0xffffLL,              BGEN_u16_buf             }, \
   { "I32", 'i', 4,  1,     -0x80000000LL,         0x7fffffffLL,          BGEN_deinterlace_d32_buf }, \
   { "U32", 'I', 4,  0,     0,                     0xffffffffLL,          BGEN_u32_buf             }, \
   { "I64", 0,   8,  1,     -0x8000000000000000LL, 0x7fffffffffffffffLL,  BGEN_deinterlace_d64_buf }, \
   { "U64", 0,   8,  0,     0,                     0x7fffffffffffffffLL,  BGEN_u64_buf             }, /* note: our internal representation is int64_t so max is limited by that */ \
   { "F32", 'f', 4,  0,     0,                     0,                     BGEN_u32_buf             }, \
   { "F64", 0,   8,  0,     0,                     0,                     BGEN_u64_buf             }, \
   { "BLB", 0,   1,  0,     0,                     0,                     0                        }, \
   { "BMP", 0,   8,  0,     0,                     0,                     0                        }, \
   { "COD", 0,   1,  0,     0,                     0,                     0                        }, \
   { "T8 ", 0,   1,  0,     0,                     0xffLL,                BGEN_transpose_u8_buf    }, \
   { "T16", 0,   2,  0,     0,                     0xffffLL,              BGEN_transpose_u16_buf   }, \
   { "T32", 0,   4,  0,     0,                     0xffffffffLL,          BGEN_transpose_u32_buf   }, \
   { "N/A", 0,   8,  0,     0,                     0x7fffffffffffffffLL,  0                        }, /* unused - can be repurposed - used to be T64, but this was never possible in the code */ \
   { "h8",  0,   1,  0,     0,                     0xffLL,                BGEN_u8_buf              }, /* lower-case UINT8 hex */ \
   { "H8",  0,   1,  0,     0,                     0xffLL,                BGEN_u8_buf              }, /* upper-case UINT8 hex */ \
   { "h16", 0,   2,  0,     0,                     0xffffLL,              BGEN_u16_buf             }, \
   { "H16", 0,   2,  0,     0,                     0xffffLL,              BGEN_u16_buf             }, \
   { "h32", 0,   4,  0,     0,                     0xffffffffLL,          BGEN_u32_buf             }, \
   { "H32", 0,   4,  0,     0,                     0xffffffffLL,          BGEN_u32_buf             }, \
   { "h64", 0,   8,  0,     0,                     0x7fffffffffffffffLL,  BGEN_u64_buf             }, \
   { "H64", 0,   8,  0,     0,                     0x7fffffffffffffffLL,  BGEN_u64_buf             }, \
   { "STR", 0,   1,  0,     0,                     0,                     0                        }, \
   { /* NUM_LTYPES */                                                                              }, \
   /* after here - not part of the file format, just used during seg */                               \
   { "DYN", 0,   8,  0,     0x8000000000000000LL,  0x7fffffffffffffffLL,  0                        }, \
   { "DYh" ,0,   8,  0,     0x8000000000000000LL,  0x7fffffffffffffffLL,  0                        }, \
   { "DYH" ,0,   8,  0,     0x8000000000000000LL,  0x7fffffffffffffffLL,  0                        }, \
}

#define lt_width(ctx)       (lt_desc[(ctx)->ltype].width)
#define lt_min(ltype)       (lt_desc[ltype].min_int)
#define lt_max(ltype)       (lt_desc[ltype].max_int)
#define lt_is_signed(ltype) (lt_desc[ltype].is_signed)
