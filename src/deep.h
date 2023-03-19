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

#define num_hash_bits z_file->deep_hash.prm8[0] // number of bits of seq_hash used for hash table (i.e. hash table is of size 2^num_hash_bits)

#define deephash_issame(a,b) (!memcmp (&(a), &(b), sizeof (DeepHash)))
#define DEEPHASHf(a) (a).qname, (a).seq, (a).qual // for printf-like arguments

// note: in PIZ, z_file->vb_start_deep_line is indexed by vb_idx - 0-based counter of non-DEPN VBs (unlike in ZIP which is vb_i)
#define num_deepable_sam_vbs z_file->deep_index.param

// 24 bytes
typedef struct {
    DeepHash hash;              // hashes of qname, seq, qual

    uint32_t next;              // linked list
    #define NO_NEXT 0xffffffff

    union ZipZDeepPlace {
        struct { 
            uint32_t line_i;        // conusmed=0: SAM line_i within vb_i of this entry. 
            uint32_t vb_i     : 31; // consumed=1: FASTQ line_i and vb_i of the consuming line
            uint32_t consumed : 1;
        }; // 64 bit
        uint64_t place;
    };
} ZipZDeep __attribute__((aligned (8)));

// hash of textual, forward, seq
extern uint32_t deep_seq_hash (VBlockP vb, STRp(seq), bool is_revcomp);

// hash of forward, qual. note: qual_hash=0 means "QUAL not hashed", so both crc32=0 and crc32=1 get mapped to hash=1
extern uint32_t deep_qual_hash (VBlockP vb, STRp(qual), bool is_revcomp);

// hash of qname (note: in FASTQ qname is the part of DESC up to the first whitespace)
static inline uint32_t deep_qname_hash (STRp(qname)) { return QNAME_HASH (STRa(qname), -1); }

//-----------------------------------------------------------------------
// PizZDeep structure:
// Part A:
// PizZDeepFlags - 0: is_seq_compressed : 0=SEQ packed 1=SEQ compressed
//                 1: is_long_seq (if long - seq_len is 4B, if not, it is 1B)
//                 2: is_long_seq_comp  
//                 3: is_long_qual_comp
//
// Part B: exists unless segconf.deep_no_qname
// 1 Byte           : qname_len (up to 254 by BAM spec)
// qname_len bytes  : QNAME (not compressed, not nul-terminated)
//
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

typedef enum { ZDEEP_SEQ_PACKED,        // SEQ is packed to 2-bit
               ZDEEP_SEQ_COMPRESSED,    // SEQ is compressed, because it contains a non-ACGT character 
               ZDEEP_SEQ_COMPRESSED_LONG_LEN, // SEQ is compressed, and compressed_len > 255
               ZDEEP_SEQ_PERFECT,       // SEQ is a perfect copy from the reference
               ZDEEP_SEQ_PERFECT_REV,   // SEQ is a perfect copy from the reference, just revcomp
               ZDEEP_SEQ_MIS_1,         // SEQ is a 1 mismatch short of a perfect copy from the reference
               ZDEEP_SEQ_MIS_1_REV      // SEQ is a 1 mismatch short of a perfect copy from the reference, just revcomp
} PizZDeepSeqEncoding;

// PIZ Deep: if SEQ is with a perfect CIGAR (eg 151M) a copy of the reference (possibly revcomp) with at most one mismatch,
// this tells FASTQ how to reconstruct it from the reference. Two simple cases that cover most of the reads of a typical FASTQ
typedef struct __attribute__ ((__packed__)) {
    uint32_t gpos;           // Start of this SEQ in the genome (support gpos up 1 - 0xffffffff)
    uint8_t mismatch_offset; // If SEQ has one mismatch - offset of mismatch in SEQ
    char mismatch_base;      // If SEQ has one mismatch - the different base
} PizDeepSeq;

#pragma pack(1) 
typedef struct {
    uint8_t seq_encoding      : 3; // SEQ encoding according to PizZDeepSeqEncoding
    uint8_t is_long_seq       : 1; // is seq_len greater than 255
    uint8_t is_long_qual_comp : 1; // is qual_comp_len greater than 255
    uint8_t is_qname_comp     : 1; // QNAME suffix is compressed
    uint8_t unused            : 2;
} PizZDeepFlags;
#pragma pack() 

