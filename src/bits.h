// ------------------------------------------------------------------
//   bits.h
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
//   Copyright claimed on additions and modifications vs public domain.

// a module for handling arrays of bits, partially based on BitArray, a public domain code located in https://github.com/noporpoise/BitArray/. 
// The unmodified license of BitArray is as follows:
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
// Statement of Purpose
// --------------------
//
// The laws of most jurisdictions throughout the world automatically confer exclusive Copyright and Related Rights (defined below) upon the creator and subsequent owner(s) (each and all, an "owner") of an original work of authorship and/or a database (each, a "Work").
//
// Certain owners wish to permanently relinquish those rights to a Work for the purpose of contributing to a commons of creative, cultural and scientific works ("Commons") that the public can reliably and without fear of later claims of infringement build upon, modify, incorporate in other works, reuse and redistribute as freely as possible in any form whatsoever and for any purposes, including without limitation commercial purposes. These owners may contribute to the Commons to promote the ideal of a free culture and the further production of creative, cultural and scientific works, or to gain reputation or greater distribution for their Work in part through the use and efforts of others.
//
// For these and/or other purposes and motivations, and without any expectation of additional consideration or compensation, the person associating CC0 with a Work (the "Affirmer"), to the extent that he or she is an owner of Copyright and Related Rights in the Work, voluntarily elects to apply CC0 to the Work and publicly distribute the Work under its terms, with knowledge of his or her Copyright and Related Rights in the Work and the meaning and intended legal effect of CC0 on those rights.
//
// -----------------------
//
// 1. *Copyright and Related Rights.* A Work made available under CC0 may be protected by copyright and related or neighboring rights ("Copyright and Related Rights"). Copyright and Related Rights include, but are not limited to, the following:
//
//     - the right to reproduce, adapt, distribute, perform, display, communicate, and translate a Work;
//     - moral rights retained by the original author(s) and/or performer(s);
//     - publicity and privacy rights pertaining to a person's image or likeness depicted in a Work;
//     - rights protecting against unfair competition in regards to a Work, subject to the limitations in paragraph 4(a), below;
//     - rights protecting the extraction, dissemination, use and reuse of data in a Work;
//     - database rights (such as those arising under Directive 96/9/EC of the European Parliament and of the Council of 11 March 1996 on the legal protection of databases, and under any national implementation thereof, including any amended or successor version of such directive); and
//     - other similar, equivalent or corresponding rights throughout the world based on applicable law or treaty, and any national implementations thereof.
//
// 2. *Waiver.* To the greatest extent permitted by, but not in contravention of, applicable law, Affirmer hereby overtly, fully, permanently, irrevocably and unconditionally waives, abandons, and surrenders all of Affirmer's Copyright and Related Rights and associated claims and causes of action, whether now known or unknown (including existing as well as future claims and causes of action), in the Work (i) in all territories worldwide, (ii) for the maximum duration provided by applicable law or treaty (including future time extensions), (iii) in any current or future medium and for any number of copies, and (iv) for any purpose whatsoever, including without limitation commercial, advertising or promotional purposes (the "Waiver"). Affirmer makes the Waiver for the benefit of each member of the public at large and to the detriment of Affirmer's heirs and successors, fully intending that such Waiver shall not be subject to revocation, rescission, cancellation, termination, or any other legal or equitable action to disrupt the quiet enjoyment of the Work by the public as contemplated by Affirmer's express Statement of Purpose.
//
// 3. *Public License Fallback.* Should any part of the Waiver for any reason be judged legally invalid or ineffective under applicable law, then the Waiver shall be preserved to the maximum extent permitted taking into account Affirmer's express Statement of Purpose. In addition, to the extent the Waiver is so judged Affirmer hereby grants to each affected person a royalty-free, non transferable, non sublicensable, non exclusive, irrevocable and unconditional license to exercise Affirmer's Copyright and Related Rights in the Work (i) in all territories worldwide, (ii) for the maximum duration provided by applicable law or treaty (including future time extensions), (iii) in any current or future medium and for any number of copies, and (iv) for any purpose whatsoever, including without limitation commercial, advertising or promotional purposes (the "License"). The License shall be deemed effective as of the date CC0 was applied by Affirmer to the Work. Should any part of the License for any reason be judged legally invalid or ineffective under applicable law, such partial invalidity or ineffectiveness shall not invalidate the remainder of the License, and in such case Affirmer hereby affirms that he or she will not (i) exercise any of his or her remaining Copyright and Related Rights in the Work or (ii) assert any associated claims and causes of action with respect to the Work, in either case contrary to Affirmer's express Statement of Purpose.
//
// 4. *Limitations and Disclaimers.*
//
//     - No trademark or patent rights held by Affirmer are waived, abandoned, surrendered, licensed or otherwise affected by this document.
//     - Affirmer offers the Work as-is and makes no representations or warranties of any kind concerning the Work, express, implied, statutory or otherwise, including without limitation warranties of title, merchantability, fitness for a particular purpose, non infringement, or the absence of latent or other defects, accuracy, or the present or absence of errors, whether or not discoverable, all to the greatest extent permissible under applicable law.
//     - Affirmer disclaims responsibility for clearing rights of other persons that may apply to the Work or any use thereof, including without limitation any person's Copyright and Related Rights in the Work. Further, Affirmer disclaims responsibility for obtaining any necessary consents, permissions or other rights required for any use of the Work.
//     - Affirmer understands and acknowledges that Creative Commons is not a party to this document and has no duty or obligation with respect to this CC0 or use of the Work.
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#pragma once

#include "genozip.h"
#include "buf_struct.h"
#include "buffer.h"

// Types 

// "The Bits invariant" that all functions need to comply with:
// 1. Top bits of final word beyond .nbits must always be 0. So:
//    A. Functions cannot pollute the top bits
//    B. Any function that changes .nwords needs to clear the unused bits in the new top word
// 2. No function should rely on these bits being 0
// 3. When writing, functions need to preserve the top unused bits being 0
// 4. Overlay: Two overlaid (small) Bits on a larger Bits cannot overlap. Each overlaid Bits must comply with these Invariants independetly.
// 5. Bits.words must be 64bit-word-aligned, including in overlaid Bits (therefore unused top bits of one overlay cannot belong to another overlay) 
typedef Buffer Bits;

typedef uint8_t word_offset_t; // Offset within a 64 bit word

// Macros 

#define BITS_IN_WORD 64

// trailing_zeros is number of least significant zeros, leading_zeros is most significant zeros.
#define trailing_zeros(x) ((x) ? (__typeof(x))__builtin_ctzll(x) : (__typeof(x))sizeof(x)*8)
#define leading_zeros(x)  ((x) ? (__typeof(x))__builtin_clzll(x) : (__typeof(x))sizeof(x)*8)

#define roundup_bits2bytes(bits)   (((bits)+7)/8)
#define roundup_bits2words64(bits) (((bits)+63)/64)
#define roundup_bits2bytes64(bits) (roundup_bits2words64(bits)*8) // number of bytes in the array of 64b words needed for bits

// Round a number up to the nearest number that is a power of two (fixed by Divon)
#define roundup2pow(x) ((__builtin_popcountll(x)==1) ? (x) : ((__typeof(x))1 << (64 - leading_zeros(x))))

// A bitmask is a value with the (nbits) lower bits set to 1, and the remaining bits set to 0: caller guarantees nbits > 0 (saves a branch)
// Caller guarantees nbits ∈ [1, bits_in_type]
#define bitmask8_(nbits)  ((uint8_t) ~(uint8_t) 0 >> (sizeof(uint8_t) *8-(nbits))) 
#define bitmask64_(nbits) ((uint64_t)~(uint64_t)0 >> (sizeof(uint64_t)*8-(nbits))) 

// Caller guarantees nbits ∈ [0, bits_in_type]  <--- nbits is allowed to be 0
#define bitmask8(nbits)   (((nbits) > 0) ? bitmask8_(nbits)  : (uint8_t) 0)
#define bitmask64(nbits)  (((nbits) > 0) ? bitmask64_(nbits) : (uint64_t)0)

// combine two words: for "1" bits in the abits, take the bit from a, and for "0" take the bit from b
#define bitmask_merge(a,b,abits) (b ^ ((a ^ b) & abits))  // equivalent to ((a & abits) | (b & ~abits))

// Bit functions on arrays - wrd:idx version
#define bitset2_get(arr,wrd,idx)       (((arr)[wrd] >> (idx)) & 0x1)
#define bitset2_get2(arr,wrd,idx)      (((arr)[wrd] >> (idx)) & 0x3) 
#define bitset2_get4(arr,wrd,idx)      (((arr)[wrd] >> (idx)) & 0xf) 
#define bitset2_set(arr,wrd,idx)       ((arr)[wrd] |= (1ULL << (idx)))
#define bitset2_del(arr,wrd,idx)       ((arr)[wrd] &= ~(1ULL << (idx)))
#define bitset2_cpy(arr,wrd,idx,bit)   ((arr)[wrd] = ((arr)[wrd] & ~(1ULL << (idx))) | ((uint64_t)(bit)  << (idx)))
#define bitset2_cpy2(arr,wrd,idx,bit2) ((arr)[wrd] = ((arr)[wrd] & ~(3ULL << (idx))) | ((uint64_t)(bit2) << (idx))) // copy 2 bits - idx must be an even number

#define bits_wrd(pos) ((pos) >> 6)
#define bits_idx(pos) ((pos) & 63)

// convert from pos version to wrd:idx version
#define bitset_op(func,arr,pos)      func((arr), bits_wrd(pos), bits_idx(pos))
#define bitset_op2(func,arr,pos,bit) func((arr), bits_wrd(pos), bits_idx(pos), (bit))

// Bit functions on arrays - pos version
#define bitset_get(arr,pos)       bitset_op (bitset2_get,  (arr), (pos))
#define bitset_get2(arr,pos)      bitset_op (bitset2_get2, (arr), (pos)) 
#define bitset_get4(arr,pos)      bitset_op (bitset2_get4, (arr), (pos)) 
#define bitset_set(arr,pos)       bitset_op (bitset2_set,  (arr), (pos))
#define bitset_del(arr,pos)       bitset_op (bitset2_del,  (arr), (pos))
#define bitset_cpy(arr,pos,bit)   bitset_op2(bitset2_cpy,  (arr), (pos), (bit))
#define bitset_cpy2(arr,pos,bit2) bitset_op2(bitset2_cpy2, (arr), (pos)/*must be even*/, (bit2)/*2 bits*/)

// clear excess bits in top used word of bitmap 
static inline void bits_clear_excess_bits_in_top_word (ConstBitsP bits) 
{
    word_offset_t top_bits = bits->nbits & 63;
    if (top_bits)
        // note: this is a read-modify-write operation, if the last word is not initialized then valgrind etc will shout upon the read
        bits->words[bits->nwords-1] &= bitmask64_(top_bits);
}

ℬ𝒾ℊℰ (extern void LTEN_bits (BitsP bits);) 
ℒ𝒾𝓉ℰ (static inline void LTEN_bits (BitsP bits) {} )

//
// Get, set, clear, assign and toggle individual bits
// Inline functions for fast access -- beware: no bounds checking. 
//
static inline bool    bits_get (ConstBitsP arr, uint64_t i) { return bitset_get (arr->words, i); }
static inline uint8_t bits_get2(ConstBitsP arr, uint64_t i) { return bitset_get2(arr->words, i); }
static inline uint8_t bits_get4(ConstBitsP arr, uint64_t i) { return bitset_get4(arr->words, i); }
static inline void bits_set    (BitsP arr, uint64_t i) { bitset_set(arr->words, i); }
static inline void bits_clear  (BitsP arr, uint64_t i) { bitset_del(arr->words, i); }
static inline void bits_assign (BitsP arr, uint64_t i, bool value) { bitset_cpy(arr->words, i, value); }
static inline void bits_assign2(BitsP arr, uint64_t i/*even number*/, uint8_t value/*0,1,2 or 3*/) { bitset_cpy2(arr->words, i, value); }

extern void bits_reverse (BitsP bits);

extern void bits_set_region (BitsP bits, uint64_t start, uint64_t len);

extern void bits_clear_region (BitsP bits, uint64_t start, uint64_t len);

extern void bits_bit_to_byte (uint8_t *restrict dst, ConstBitsP src_bits, uint64_t src_bit, uint32_t num_bits);

// create words - if the word is not aligned to the bitmap word boundaries, and hence spans 2 bitmap words, 
// we take the MSb's from the left word and the LSb's from the right word 
static inline uint64_t _bits_combined_word (uint64_t word_a, uint64_t word_b, 
                                            int shift) // must be 0-63, undefined behavior otherwise
{
#ifdef __x86_64__
    __asm__ ("shrdq %2, %1, %0" // single-cycle hardware Shift Right Double instruction
             : "+r" (word_a)
             : "r" (word_b), "cJ" ((char)shift));
    return word_a;

#else // fallback (carefully avoid branches)
    uint64_t high_part = word_a >> shift;
    uint64_t low_part  = word_b << ((64 - shift) & 63);
    
    return high_part | (low_part & -(int64_t)(shift > 0));
#endif    
}         

// gets a word - careful not to access a word beyond the end 
static inline uint64_t _get_word_(const uint64_t *words, uint64_t nwords, uint64_t start)
{
    uint64_t word_index = bits_wrd(start);
    int word_offset     = bits_idx(start);

    uint64_t word_a = words[word_index];
    uint64_t word_b = (word_index + 1 < nwords) ? words[word_index + 1] : 0ULL; // note: expected to be branchless: likely compiled to a conditional move rather than a conditional jump

    return _bits_combined_word (word_a, word_b, word_offset);    
}

// gets a word - careful not to access a word beyond the end 
static inline uint64_t _get_word (ConstBitsP bits, uint64_t start)
{
    return (_get_word_(bits->words, bits->nwords, start));
}

static inline uint64_t bits_get_wordn (ConstBitsP bits, uint64_t start, int n /* up to 64 */)
{
  ASSERT (start + n <= bits->nbits, "expecting start=%"PRIu64" + n=%d <= bits->nbits=%"PRIu64, 
          start, n, bits->nbits);
           
  return (uint64_t)(_get_word(bits, start) & bitmask64(n));
}

extern void _set_word (BitsP bits, uint64_t start, uint64_t word);

static inline void bits_set_wordn (BitsP bits, uint64_t start, uint64_t word, int n) 
{
    uint64_t w = _get_word (bits, start), m = bitmask64(n);
    _set_word (bits, start, bitmask_merge (word,w,m));
}

extern bool bits_is_fully_set (ConstBitsP bits);
extern bool bits_is_fully_clear (ConstBitsP bits);
extern uint64_t bits_num_set_bits (ConstBitsP bits);
extern uint64_t bits_num_set_bits_region (ConstBitsP bits, uint64_t start, uint64_t length); 

// Get the number of bits not set (length - hamming weight)
extern uint64_t bits_num_clear_bits (ConstBitsP bits);

// Find the index of the next bit that is set, at or after `offset`
// Returns 1 if a bit is set, otherwise 0
// Index of next set bit is stored in the integer pointed to by result
// If no next bit is set result is not changed
extern bool bits_find_next_set_bit   (ConstBitsP bits, uint64_t offset, uint64_t *result);
extern bool bits_find_next_clear_bit (ConstBitsP bits, uint64_t offset, uint64_t *result);
extern bool bits_find_prev_set_bit   (ConstBitsP bits, uint64_t offset, uint64_t *result);
extern bool bits_find_prev_clear_bit (ConstBitsP bits, uint64_t offset, uint64_t *result);
extern bool bits_find_first_set_bit  (ConstBitsP bits, uint64_t *result);
extern bool bits_find_first_clear_bit(ConstBitsP bits, uint64_t *result);
extern bool bits_find_last_set_bit   (ConstBitsP bits, uint64_t *result);
extern bool bits_find_last_clear_bit (ConstBitsP bits, uint64_t *result);

// move a range of bits to a lower index within the same array
extern void bits_sink_range (ConstBits𐤐 bits, uint64_t dstindx, uint64_t srcindx, uint64_t length, bool src_may_exceed_nbits);

// Copy bits between Bits arrays
extern void bits_copy (BitsP dst, uint64_t dstindx, ConstBitsP src, uint64_t srcindx, uint64_t length);

// revcomp one word: the order of the 32 bases in a word is swapped, and each base is complemented A⇔T C⇔G
static inline uint64_t bits_revcomp_word (uint64_t w) 
{
#if defined(__clang__) && defined(__aarch64__) // 4-5 clock cycles on ARM (bitreverse64 does not exist natively on x86)
    w = __builtin_bitreverse64(w); // reverse bit order, however this also swaps the bits within every 2bits 

    // swap back the bits within every 2bits
    w = ((w & 0xAAAAAAAAAAAAAAAAULL) >> 1) | // picks bits 1, 3, 5...
        ((w & 0x5555555555555555ULL) << 1);  // picks bits 0, 2, 4...

#else // 4 clock cycles on Intel, 7-9 on ARM
    w = ((w >> 2) & 0x3333333333333333ULL) | ((w & 0x3333333333333333ULL) << 2);  // within every nibble, swap the first 2 bits with the last 2 bits
    w = ((w >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((w & 0x0F0F0F0F0F0F0F0FULL) << 4);  // within every byte, swap the first nibble with the last nibble
    w = __builtin_bswap64(w); // reverse the 8 bytes of the word 
#endif
 
    return ~w; // complement: A(00)⇔T(11) C(01)⇔G(10)
}

extern void bits_overlay (BitsP overlaid_bits, BitsP regular_bits, uint64_t start, uint64_t nbits, rom name);

// removes flanking bits on boths sides, shrinking bits 
extern void bits_remove_flanking (BitsP bits, uint64_t lsb_flanking, uint64_t msb_flanking);

extern uint64_t bits_resize (BitsP bits, uint64_t new_num_of_bits);

extern void bits_add_bit (BitsP bits, int64_t new_bit);

extern uint32_t bits_hamming_distance (ConstBits𐤐 bits_1, ConstBits𐤐 bits_2, uint64_t index_2);

// get the number of consecutive 1s at the start of a region 
static inline uint64_t bits_get_run (ConstBitsP bits, uint64_t start, uint64_t len) 
{
    if (__builtin_expect(!len, 0)) return 0;

    // FIX 1: Overflow-safe bounds check (prevents start + len wrapping 64-bit uint)
    ASSERT(start <= bits->nbits && len <= bits->nbits - start,
           "bits_get_run out of bounds: start=%"PRIu64", len=%"PRIu64", nbits=%"PRIu64, start, len, bits->nbits);

    uint64_t word_index        = bits_wrd (start);
    uint64_t last_word_index   = bits_wrd (start + len - 1);
    word_offset_t start_offset = bits_idx (start);
    word_offset_t after_offset = bits_idx (start + len);

    // FIX 2: Use uint64_t for run calculation to avoid signed integer underflow/overflow
    uint64_t run = 0;

    for (const uint64_t *start_w = &bits->words[word_index], 
                        *last_w  = &bits->words[last_word_index], 
                        *w       = start_w;
         w <= last_w; 
         w++) {
        
        uint64_t word = ~*w; // Invert: consecutive 1s become consecutive 0s

        // handle partial last word (inject a 1 to terminate the run right after 'len')
        if (w == last_w && after_offset)
            word |= (uint64_t)1 << after_offset;

        // handle partial first word (zero out bits before 'start')
        if (w == start_w) 
            word &= ~bitmask64(start_offset); 

        uint64_t zeros = trailing_zeros(word);

        run += (w == start_w) ? (zeros - start_offset) : zeros;

        if (zeros != 64) break; // Found a 0 in original word (terminator hit)
    }

    return run;
}


// Stringification
extern StrText bits_word_to_01_string (uint64_t word, int n_bits/*0 to 64*/);
extern StrText1K bits_to_01_string (ConstBitsP bits, uint64_t start, int64_t length);
extern StrText1K bits_to_ACGT_string (ConstBitsP bits, uint64_t start_base, uint64_t n_bases);

// Debug functions
#ifdef DEBUG
extern void validate_bits (ConstBitsP bits, Caller caller);
#define DEBUG_VALIDATE_BITS(a) validate_bits((a), THIS_CODE_LINE)
#else
#define DEBUG_VALIDATE_BITS(a)
#endif
