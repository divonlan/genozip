// ------------------------------------------------------------------
//   deep.h
//   Copyright (C) 2023-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// DATA FORMAT IN vb->deep_ent (PIZ):
// note: "1 or 5 bytes" integers: if <= 254 store in 1 byte, else 255 followed by 4 bytes
// note: "offset" means the number of non-included bases (i.e. the "gap") since the previous offset or beginning of string
//
// -- QNAME --
// 1 byte: qname_len: uncompressed length of canonical qname
// (var length) qname: huffman-compressed for files starting 15.0.65 and uncompressed for older files
//
// -- SEQ --
// 1 byte: PizDeepSeqFlags
// 1 or 5 bytes: deep_ref_consumed (if has_cigar and not --qual_only)
//               deep_seq_len (if !has_cigar or --qual_only)
// ***
// *** OPTION 1: if seq_encoding == ZDEEP_SEQ_BY_REF *** 
// 4 or 5 bytes: gpos - depending on PizDeepSeqFlags.is_long_gpos). gpos is of last reference base if SAM_FLAG.rev_comp.
// 1 or 5 bytes: num_mismatches (omitted if !PizDeepSeqFlags.has_mismatches)
// (var length) num_mismatches x (offset, base). offset is 1 or 5 bytes, base is one ASCII character. Note: mismatches are revcomped if SAM_FLAG.rev_comp
// 1 or 5 bytes: uncomp_cigar_len (omitted if !PizDeepSeqFlags.has_cigar)
// (var length) textual cigar (reversed if SAM_FLAG.rev_comp). huffman-compressed in files since 15.0.68 and uncompressed for older files.
// 1 or 5 bytes: nonref_len (omitted if !PizDeepSeqFlags.has_cigar)
// (var length) packed nonref (i.e. S,I data) (revcomped if SAM_FLAG.rev_comp) 
// ***
// *** OPTION 2: if seq_encoding == ZDEEP_SEQ_PACKED *** (if all are 'A','C','G','T')
// (var length) 2-bit packed seq (revcomped if SAM_FLAG.rev_comp) 
// ***
// *** OPTION 3: if seq_encoding == ZDEEP_SEQ_ARITH *** (if there are IUPACs)
// 1 or 5 bytes: seq_comp_len 
// (var length) Arith-compressed seq (revcomped if SAM_FLAG.rev_comp) 
//
// -- QUAL --
// *** OPTION 1: if qual_is_monochar=false
// 1 or 5 bytes: comp_len  
// (var length) huffman or arith compression (reversed if SAM_FLAG.rev_comp)
// ***
// *** OPTION 2: if qual_is_monochar=true
// 1 byte: monochar  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#pragma once

#include "qname.h"

#define deephash_issame(a,b) (!memcmp (&(a), &(b), sizeof (DeepHash)))
#define DEEPHASHf(a) (a).qname, (a).seq, (a).qual // for printf-like arguments

// note: in PIZ, z_file->vb_start_deep_line is indexed by vb_idx - 0-based counter of non-DEPN VBs (unlike in ZIP which is vb_i)
#define num_deepable_sam_vbs z_file->deep_index.prm32[0] // PIZ

// 32 bytes
typedef struct {
    DeepHash hash;                  // hashes of qname (64b), seq(32b), qual (32b)
    #define NO_NEXT 0xffffffff
    #define MAX_DEEP_ENTS 0xfffffffe// our current implementation is limited to 4G reads because the linked list is 32b
    uint32_t next;                  // linked list of entries with the with the same (hash.qname & mask)
    uint32_t seq_len;

    union ZipZDeepPlace {
        struct { 
            uint32_t line_i;        // ⸢ consumed=0: SAM line_i within vb_i of this entry. 
            uint32_t vb_i     : 30; // ⸤ consumed=1: FASTQ line_i and vb_i of the consuming line
            uint32_t consumed : 1;  // ZIP FASTQ: a FASTQ read matching this entry has been found. note: if dup=true, consumed will never be set  
            uint32_t dup      : 1;  // ZIP FASTQ: this is one of 2+ entries with the same hash values.
        }; // 64 bit
        uint64_t place;             // used to access place data as a single integer
    };
} ZipZDeep __attribute__((aligned (8))); // needs to be word-aligned so we can __atomic_compare_exchange_n of place

// hash of textual, forward, seq
extern uint32_t deep_seq_hash (VBlockP vb, STRp(seq), bool is_revcomp);

// hash of forward, qual. note: qual_hash=0 means "QUAL not hashed", so both crc32=0 and crc32=1 get mapped to hash=1
extern uint32_t deep_qual_hash (VBlockP vb, STRp(qual), bool is_revcomp);

// ZIP: hash of canonical qname (note: in FASTQ qname is the part of DESC up to the first whitespace)
extern uint64_t deep_qname_hash (VBlockP vb, QType q, STRp(qname), thool is_last, uint32_t *uncanonical_suffix_len);

typedef enum { ZDEEP_SEQ_PACKED,   // SEQ is packed to 2-bit
               ZDEEP_SEQ_BY_REF,   // SEQ is to be reconstructed from the reference
               ZDEEP_SEQ_ARITH,    // SEQ is arith rather than packed, because it contains a non-ACGT character 
} PizZDeepSeqEncoding;

// data format in vb->deep_ent: see description in sam_piz_con_item_cb

typedef struct { // 1 byte
    uint8_t seq_encoding      : 2; // SEQ encoding according to PizZDeepSeqEncoding
    uint8_t is_forward        : 1; // true if the FASTQ seq is forward relative to the reference 
    uint8_t is_long_gpos      : 1; // used 64b gpos (used if ZDEEP_SEQ_BY_REF)
    uint8_t has_cigar         : 1; // CIGAR. false if aligner was used or cigar is "perfect" (e.g. 151M) (used if ZDEEP_SEQ_BY_REF)
    uint8_t has_mismatches    : 1; // mismatches are encoded in deep_ents (used if ZDEEP_SEQ_BY_REF)
    uint8_t qual_is_monochar  : 1;
    uint8_t unused            : 1;
} PizDeepSeqFlags;

extern uint8_t num_hash_bits; 

#define DEEP_NUM_SHOW 5 // for flag.show_deep

// type of z_file->vb_start_deep_line.param
typedef enum { DEEP_PIZ_ENTS_UNINITIALIZED, DEEP_PIZ_ENTS_RECONSTRUCTING, DEEP_PIZ_ENTS_READY } DeepPizEntsStatus;