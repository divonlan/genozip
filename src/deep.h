// ------------------------------------------------------------------
//   deep.h
//   Copyright (C) 2023-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "qname.h"
#include "libdeflate/libdeflate.h"

#define IS_DEEPABLE(f) (!(f).secondary && !(f).supplementary)

#define num_hash_bits z_file->deep_hash.prm8[0] // number of bits of seq_hash used for hash table (i.e. hash table is of size 2^num_hash_bits)

#define deephash_issame(a,b) (!memcmp (&(a), &(b), sizeof (DeepHash)))
#define DEEPHASHf(a) (a).qname, (a).seq, (a).qual // for printf-like arguments

#define VB_start_deep_line(vb_i) (*B64(z_file->vb_start_deep_line, (vb_i)-1))

typedef struct {
    DeepHash hash;
    VBIType  vb_i;    // Initially: SAM VB of this entry. After consumption of entry by FASTQ: vb_i of FASTQ that consumed
    uint32_t line_i;  // Initially: SAM line_i within vb_i of this entry. After consumption: line_i of consuming FASTQ VB | 0x80000000 (i.e., MSb set)
    uint32_t next;    // linked list
} ZipZDeep;
#define NO_NEXT 0xffffffff

// hash of textual, forward, seq
extern uint32_t deep_seq_hash (VBlockP vb, STRp(seq), bool is_revcomp);

// hash of forward, qual. note: qual_hash=0 means "QUAL not hashed", so both crc32=0 and crc32=1 get mapped to hash=1
extern uint32_t deep_qual_hash (VBlockP vb, STRp(qual), bool is_revcomp);

// hash of qname (note: in FASTQ qname is the part of DESC up to the first whitespace)
static inline uint32_t deep_qname_hash (STRp(qname)) { return QNAME_HASH (STRa(qname), -1 ); }

//-----------------------------------------------------------------------
// PizZDeep structure:
// Part A:
// PizZDeepFlags - 0: is_seq_compressed : 0=SEQ packed 1=SEQ compressed
//                 1: is_long_seq (if long - len is 4B, if not, it is 1B)
//                 2: is_long_seq_comp  
//                 3: is_long_qual_comp
//
// Part B: exists unless segconf.deep_no_qname
// 1 Byte           : qname_len (up to 254 by BAM spec)
// qname_len bytes  : QNAME (not compressed, not nul-terminated)
// 1 or 4 bytes     : seq_len (of uncompressed SEQ)
//
// Part C:  (if is_seq_compressed=0) 
// ⌈seq_len/4⌉ bytes : 2bit packed SEQ
//
// Part C:  (if is_seq_compressed=1)
// 1 or 4 bytes     : seq_comp_len
// seq_comp_len bts : compressed SEQ
//
// Part D:
// 1 or 4 bytes     : qual_comp_len
// qual_comp_len bts: compressed QUAL
//-------------------------------------------------------------

#pragma pack(1) 
typedef struct {
    uint8_t is_seq_compressed : 1; // 1 if compressed, 0 if packed
    uint8_t is_long_seq       : 1; // is seq_len greater than 255
    uint8_t is_long_seq_comp  : 1; // is seq_comp_len greater than 255
    uint8_t is_long_qual_comp : 1; // is qual_comp_len greater than 255
    uint8_t unused            : 4;
} PizZDeepFlags;
#pragma pack() 
