// ------------------------------------------------------------------
//   bits.h
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
//   Copyright claimed on additions and modifications vs public domain.
//
// a module for handling arrays of 2-bit elements, partially based on Bits, a public domain code located in https://github.com/noporpoise/BitArray/. 
// The unmodified license of Bits is as follows:
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

//------------------------------------------------------------------
// Macros 
//------------------------------------------------------------------

#define WORD_SIZE 64

// trailing_zeros is number of least significant zeros
// leading_zeros is number of most significant zeros
#ifndef __GNUC__
    #define trailing_zeros(x) ({ __typeof(x) _r; _BitScanReverse64(&_r, (x)); _r; })
    #define leading_zeros(x)  ({ __typeof(x) _r; _BitScanForward64(&_r, (x)); _r; })
#else
    #define trailing_zeros(x) ((x) ? (__typeof(x))__builtin_ctzll(x) : (__typeof(x))sizeof(x)*8)
    #define leading_zeros(x)  ((x) ? (__typeof(x))__builtin_clzll(x) : (__typeof(x))sizeof(x)*8)
#endif

// Get index of top set bit. If x is 0 return nbits
#define top_set_bit(x) ((x) ? sizeof(x)*8-leading_zeros(x)-1 : sizeof(x)*8)

#define roundup_bits2bytes(bits)   (((bits)+7)/8)
#define roundup_bits2words32(bits) (((bits)+31)/32)
#define roundup_bits2words64(bits) (((bits)+63)/64)
#define roundup_bits2bytes64(bits) (roundup_bits2words64(bits)*8) // number of bytes in the array of 64b words needed for bits
#define roundup_bytes2bytes64(bytes) (((bytes)+7)  & ~(__typeof(bytes))0b111)   // rounds up bytes to nearest 8 (i.e. 64b word) 
#define roundup_bases2bytes64(bases) (((bases)+31) & ~(__typeof(bases))0b11111) // rounds up bases(=2bits) to nearest 32 (i.e. 64b word) 

// Round a number up to the nearest number that is a power of two (fixed by Divon)
#define roundup2pow(x) ((__builtin_popcountll(x)==1) ? (x) : ((__typeof(x))1 << (64 - leading_zeros(x))))

#define rot32(x,r) (((x)<<(r)) | ((x)>>(32-(r))))
#define rot64(x,r) (((x)<<(r)) | ((x)>>(64-(r))))

// A bitmask is a value with the (nbits) lower bits set to 1, and the remaining bits set to 0
#define bitmask(nbits, type) ((nbits) ? (type)~(type)0 >> (sizeof(type)*8-(nbits)) : (type)0)
#define bitmask8(nbits)  bitmask(nbits, uint8_t)
#define bitmask16(nbits) bitmask(nbits, uint16_t)
#define bitmask32(nbits) bitmask(nbits, uint32_t)
#define bitmask64(nbits) bitmask(nbits, uint64_t)

// combine two words: for "1" bits in the abits, take the bit from a, and for "0" take the bit from b
#define bitmask_merge(a,b,abits) (b ^ ((a ^ b) & abits))  // equivalent to ((a & abits) | (b & ~abits))

//
// Bit array (bitset)
//
// bitsetX_wrd(): get word for a given position
// bitsetX_idx(): get index within word for a given position

#define _TYPESHIFT(arr,word,shift) \
        ((__typeof(*(arr)))((__typeof(*(arr)))(word) << (shift)))

#define bitsetX_wrd(wrdbits,pos) ((pos) / (wrdbits))
#define bitsetX_idx(wrdbits,pos) ((pos) % (wrdbits))

#define bitset32_wrd(pos) ((pos) >> 5)
#define bitset32_idx(pos) ((pos) & 31)

#define bitset64_wrd(pos) ((pos) >> 6)
#define bitset64_idx(pos) ((pos) & 63)

//
// Bit functions on arrays
//
#define bitset2_get(arr,wrd,idx)     (((arr)[wrd] >> (idx)) & 0x1)
#define bitset2_get2(arr,wrd,idx)    (((arr)[wrd] >> (idx)) & 0x3) // divon
#define bitset2_get4(arr,wrd,idx)    (((arr)[wrd] >> (idx)) & 0xf) // divon
#define bitset2_set(arr,wrd,idx)     ((arr)[wrd] |=  _TYPESHIFT((arr),1,(idx)))
#define bitset2_del(arr,wrd,idx)     ((arr)[wrd] &=~ _TYPESHIFT((arr),1,(idx)))
#define bitset2_tgl(arr,wrd,idx)     ((arr)[wrd] ^=  _TYPESHIFT((arr),1,(idx)))
#define bitset2_or(arr,wrd,idx,bit)  ((arr)[wrd] |=  _TYPESHIFT((arr),(bit),idx))
#define bitset2_xor(arr,wrd,idx,bit) ((arr)[wrd]  = ~((arr)[wrd] ^ (~_TYPESHIFT((arr),(bit),(idx)))))
#define bitset2_and(arr,wrd,idx,bit) ((arr)[wrd] &= (_TYPESHIFT((arr),(bit),(idx)) | ~_TYPESHIFT((arr),1,(idx))))
#define bitset2_cpy(arr,wrd,idx,bit) ((arr)[wrd]  = ((arr)[wrd] &~ _TYPESHIFT((arr),1,(idx))) | _TYPESHIFT((arr),(bit),(idx)))

// divon - copy 2 bits - idx must be an even number
#define bitset2_cpy2(arr,wrd,idx,bit2) ((arr)[wrd] = ((arr)[wrd] &~ _TYPESHIFT((arr),0x3,(idx))) | _TYPESHIFT((arr),(bit2),(idx)))

//
// Auto detect size of type from pointer
//
#define bitset_wrd(arr,pos) bitsetX_wrd(sizeof(*(arr))*8,(pos))
#define bitset_idx(arr,pos) bitsetX_idx(sizeof(*(arr))*8,(pos))
#define bitset_op(func,arr,pos)      func((arr), bitset_wrd((arr),(pos)), bitset_idx((arr), (pos)))
#define bitset_op2(func,arr,pos,bit) func((arr), bitset_wrd((arr), (pos)), bitset_idx((arr), (pos)), (bit))

// Auto-detect type size: bit functions
#define bitset_get(arr,pos)       bitset_op (bitset2_get,  (arr), (pos))
#define bitset_get2(arr,pos)      bitset_op (bitset2_get2, (arr), (pos)) // divon - gets a value [0,3] of 2 bits 
#define bitset_get4(arr,pos)      bitset_op (bitset2_get4, (arr), (pos)) // divon - gets a value [0,15] of 4 bits 
#define bitset_set(arr,pos)       bitset_op (bitset2_set,  (arr), (pos))
#define bitset_del(arr,pos)       bitset_op (bitset2_del,  (arr), (pos))
#define bitset_tgl(arr,pos)       bitset_op (bitset2_tgl,  (arr), (pos))
#define bitset_or(arr,pos,bit)    bitset_op2(bitset2_or,   (arr), (pos), (bit))
#define bitset_xor(arr,pos,bit)   bitset_op2(bitset2_xor,  (arr), (pos), (bit))
#define bitset_and(arr,pos,bit)   bitset_op2(bitset2_and,  (arr), (pos), (bit))
#define bitset_cpy(arr,pos,bit)   bitset_op2(bitset2_cpy,  (arr), (pos), (bit))
#define bitset_cpy2(arr,pos,bit2) bitset_op2(bitset2_cpy2, (arr), (pos), (bit2))

// Clearing a word does not return a meaningful value
#define bitset_clear_word(arr,pos) ((arr)[bitset_wrd((arr),(pos))] = 0)

//------------------------------------------------------------------

typedef uint8_t word_offset_t; // Offset within a 64 bit word

#define BIT_INDEX_MIN 0
#define BIT_INDEX_MAX (~(uint64_t)0)

//
// Structs
//

typedef struct Bits
{
    // These fields should not be changed or added to, as they map to Buffer
    uint64_t *words;     // maps to Buffer->data
    uint64_t nbits;      // maps to Buffer->nbits
    uint64_t nwords;     // maps to Buffer->len  ; round_up (nbits / 64)
    uint64_t type  : BUF_TYPE_BITS;  // maps to Buffer->type
    uint64_t other : 64-BUF_TYPE_BITS; // used only if this is a Buffer
} Bits;

#define bits_in_top_word(nbits) ((nbits) ? bitset64_idx((nbits) - 1) + 1 : 0)

// clear excess bits in top used word of bitmap 
static inline void bits_clear_excess_bits_in_top_word (BitsP bits)
{
    if (bits->nwords) {
        word_offset_t bits_active = bits_in_top_word(bits->nbits);
        //atomic version of &= in case multiple threads do this in concurrently while working on different regions of the bitmap
        __atomic_and_fetch (&bits->words[bits->nwords-1], bitmask64(bits_active), __ATOMIC_RELAXED);
    }
}

static inline void ASSERT_excess_bits_are_0 (ConstBitsP bits) // divon
{
  ASSERT0 (!(bits->nbits & 0x3f) || !(bits->words[bits->nwords-1] & ~bitmask64 (bits->nbits & 0x3f)),
           "Expecting excess bits in top word of bit array to be 0");
}


extern Bits bits_init_do (uint64_t nbits, uint8_t *data, uint64_t data_len, bool clear, FUNCLINE);
#define bits_init(nbits, data, data_len, clear) bits_init_do ((nbits), (data), (data_len), (clear), __FUNCLINE)

extern Bits bits_alloc_do (uint64_t nbits, bool clear, FUNCLINE);
#define bits_alloc(nbits, clear) bits_alloc_do ((nbits), (clear), __FUNCLINE)

extern void bits_free (BitsP bits);

extern void LTEN_bits (BitsP bits); // divon

//
// Macros
//

//
// Get, set, clear, assign and toggle individual bits
// Inline functions for fast access -- beware: no bounds checking. 
//
static inline bool    bits_get (ConstBitsP arr, uint64_t i) { return bitset_get (arr->words, i); }
static inline uint8_t bits_get2(ConstBitsP arr, uint64_t i) { return bitset_get2(arr->words, i); }
static inline uint8_t bits_get4(ConstBitsP arr, uint64_t i) { return bitset_get4(arr->words, i); }
static inline void bits_set    (BitsP arr, uint64_t i) { bitset_set(arr->words, i); }
static inline void bits_clear  (BitsP arr, uint64_t i) { bitset_del(arr->words, i); }
static inline void bits_toggle (BitsP arr, uint64_t i) { bitset_tgl(arr->words, i); }
static inline void bits_assign (BitsP arr, uint64_t i, bool value) { bitset_cpy(arr->words, i, value); }
static inline void bits_assign2(BitsP arr, uint64_t i/*even number*/, uint8_t value/*0,1,2 or 3*/) { bitset_cpy2(arr->words, i, value); }

//
// Get, set, clear, assign and toggle individual bits
// "Safe": use assert() to check bounds
//

// Get the value of a bit (returns 0 or 1)
extern bool bits_get_bit (ConstBitsP bits, uint64_t b);
extern void bits_set_bit (BitsP bits, uint64_t b);
extern void bits_clear_bit (BitsP bits, uint64_t b);
extern void bits_toggle_bit (BitsP bits, uint64_t b);
extern void bits_assign_bit (BitsP bits, uint64_t b, bool set);

//
// Get, set, clear and toggle several bits at once
//

// Get the offsets of the set bits (for offsets start<=offset<end)
// Returns the number of bits set
// It is assumed that dst is at least of length (end-start)
extern uint64_t bits_get_bits(ConstBitsP bits, uint64_t start, uint64_t end, uint64_t *dst);

// Set multiple bits at once.
// e.g. set bits 1, 20 & 31: bits_set_bits(bits, 3, 1,20,31);
// Note: variable args are of type unsigned int
extern void bits_set_bits(BitsP bits, size_t n, ...);

// Clear multiple bits at once.
// e.g. clear bits 1, 20 & 31: bits_clear_bits(bits, 3, 1,20,31);
// Note: variable args are of type unsigned int
extern void bits_clear_bits(BitsP bits, size_t n, ...);

// Toggle multiple bits at once
// e.g. toggle bits 1, 20 & 31: bits_toggle_bits(bits, 3, 1,20,31);
// Note: variable args are of type unsigned int
extern void bits_toggle_bits(BitsP bits, size_t n, ...);

// Reverse the whole array or part of it
extern void bits_reverse (BitsP bits);
extern void bits_reverse_region (BitsP bits, uint64_t start, uint64_t len);

//
// Set, clear and toggle all bits in a region
//

// Set all the bits in a region
extern void bits_set_region (BitsP bits, uint64_t start, uint64_t len);

// Clear all the bits in a region
#define bits_clear_region(bits,start,len) bits_clear_region_do (bits, start, len, __FUNCLINE)
extern void bits_clear_region_do (BitsP bits, uint64_t start, uint64_t len, rom func, unsigned code_line);

//
// Set, clear and toggle all bits at once
//

// Set all bits in this array to 1
extern void bits_set_all(BitsP bits);

// Set all bits in this array to 0
extern void bits_clear_all(BitsP bits);

extern void bits_bit_to_byte  (uint8_t *dst, ConstBitsP src_bits, uint64_t src_bit,  uint32_t num_bits);
extern void bits_2bit_to_byte (uint8_t *dst, ConstBitsP src_bits, uint64_t base_i, uint32_t num_2bits);
extern void bits_2bit_to_ACGT (char *dst, ConstBitsP src_bits, uint64_t base_i, uint32_t num_bases);

//
// Get and set words (low level -- no bounds checking)
//
static inline uint64_t _get_word(ConstBitsP bits, uint64_t start)
{
    uint64_t word_index = bitset64_wrd(start);
    word_offset_t word_offset = bitset64_idx(start);

    uint64_t result = bits->words[word_index] >> word_offset;

    word_offset_t bits_taken = WORD_SIZE - word_offset;

    // word_offset is now the number of bits we need from the next word
    // Check the next word has at least some bits
    if(word_offset > 0 && start + bits_taken < bits->nbits)
        result |= bits->words[word_index+1] << (WORD_SIZE - word_offset);

    return result;
}

//
// Get / set a word of a given size
//

// Bits is interpreted as an integer in Little Endian - i.e. first bit in the Bits is LSb
#define BIT_ARRAY_GET_WORD(width) \
    static inline uint64_t bits_get_word ## width (ConstBitsP bits, uint64_t start) { \
        ASSERT (start + width <= bits->nbits, "expecting start(%"PRIu64") + " #width " <= bits->nbits(%"PRIu64")", start, bits->nbits); \
        return (uint ## width ## _t)_get_word(bits, start); \
    } 
BIT_ARRAY_GET_WORD(64)
BIT_ARRAY_GET_WORD(32)
BIT_ARRAY_GET_WORD(16)
BIT_ARRAY_GET_WORD(8)

extern uint64_t bits_get_wordn(ConstBitsP bits, uint64_t start, int n);

// Set 64 bits at once from a particular start position
// Doesn't extend bit array. However it is safe to TRY to set bits beyond the
// end of the array, as long as: `start` is < `nbits`
extern void bits_set_word64 (BitsP bits, uint64_t start, uint64_t word);
extern void bits_set_word32 (BitsP bits, uint64_t start, uint32_t word);
extern void bits_set_word16 (BitsP bits, uint64_t start, uint16_t word);
extern void bits_set_word8 (BitsP bits, uint64_t start, uint8_t byte);
extern void bits_set_wordn (BitsP bits, uint64_t start, uint64_t word, int n);

//
// Number of bits set
//

// Get the number of bits set (hamming weight)
extern bool bits_is_fully_set (ConstBitsP bits);
extern bool bits_is_fully_clear (ConstBitsP bits);
extern uint64_t bits_num_set_bits (ConstBitsP bits);
extern uint64_t bits_num_set_bits_region (ConstBitsP bits, uint64_t start, uint64_t length); // added by divon

// Get the number of bits not set (length - hamming weight)
extern uint64_t bits_num_clear_bits (ConstBitsP bits);

//
// Find indices of set/clear bits
//

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


//
// String and printing methods
//

// Construct a Bits from a string.
extern void bits_from_str (BitsP bits, rom bitstr);

// Construct a Bits from a substring with given on and off characters.
extern void bits_from_substr (BitsP bits, uint64_t offset, rom str, size_t len, rom on, rom off, char left_to_right);

// Takes a char array to write to.  `str` must be bits->nbits+1 in
// length. Terminates string with '\0'
extern char* bits_to_str (ConstBitsP bits, char* str);
extern char* bits_to_str_rev (ConstBitsP bits, char* str);

// Get a string representations for a given region, using given on/off characters.
extern char *bits_to_substr (ConstBitsP bits, uint64_t start, uint64_t length, char *str);

// Print this array to a file stream.  Prints '0's and '1'.  Doesn't print newline.
extern void bits_print_do (ConstBitsP bits, rom msg, FILE *file);
#define bits_print(bits) bits_print_do (bits, #bits, info_stream)

extern void bits_print_binary_word_do (uint64_t word, rom msg, FILE *file);
#define bits_print_binary_word(word) bits_print_binary_word_do (word, #word, info_stream)

extern void bits_print_substr (rom msg, ConstBitsP bits, uint64_t start, uint64_t length, FILE *file);

extern void bits_print_substr_bases (rom msg, ConstBitsP bits, uint64_t start_base, uint64_t n_bases, FILE *file);

// Print bit array as hex
extern size_t bits_print_hex (ConstBitsP bits, uint64_t start, uint64_t length, FILE* fout, char uppercase);

// Copy bits from one array to another. Destination and source can be the same bits and src/dst regions can overlap
extern void bits_copy_do (BitsP dst, uint64_t dstindx, ConstBitsP src, uint64_t srcindx, uint64_t length, rom func, unsigned code_line);
#define bits_copy(dst, dstindex, src, srcindx, length) \
    bits_copy_do ((dst), (dstindex), (src), (srcindx), (length), __FUNCLINE)

// for each 2 bits in the src array, the dst array will contain those 2 bits in the reverse
// position, as well as transform them 00->11 11->00 01->10 10->01
// works on arrays with full words
extern void bits_reverse_complement_aligned (BitsP dst, ConstBitsP src, uint64_t src_start_base, uint64_t max_num_bases);
extern void bits_reverse_complement_in_place (BitsP bits) ;

extern void bits_overlay (BitsP overlaid_bits, BitsP regular_bits, uint64_t start, uint64_t nbits);

//
// Logic operators
//

// BIT_ARRAYs can all be different or the same object
// dest array will be resized if it is too short
//
extern void bits_and (BitsP dest, ConstBitsP src1, ConstBitsP src2);
extern void bits_or  (BitsP dest, ConstBitsP src1, ConstBitsP src2);
extern void bits_xor (BitsP dest, ConstBitsP src1, ConstBitsP src2);
extern void bits_not (BitsP dest, ConstBitsP src);

extern void bits_xor_with (BitsP dst, uint64_t dst_start_bit, ConstBitsP xor_with, uint64_t xor_with_bit, uint64_t num_bits);

//
// Comparisons
//

// Note: (bits_cmp(a,b) == 0) <=> (bits_cmp_big_endian(a,b) == 0)

// comparison functions return:
//   1 iff bits1 > bits2
//   0 iff bits1 == bits2
//  -1 iff bits1 < bits2

// Compare two bit arrays by value stored, with index 0 being the Least
// Significant Bit (LSB). Arrays do not have to be the same length.
// Example: ..0101 (5) > ...0011 (3) [index 0 is LSB at right hand side]
extern int bits_cmp(ConstBitsP bits1, ConstBitsP bits2);

// Compare two bit arrays by value stored, with index 0 being the Most
// Significant Bit (MSB). Arrays do not have to be the same length.
// Example: 10.. > 01.. [index 0 is MSB at left hand side]
extern int bits_cmp_big_endian(ConstBitsP bits1, ConstBitsP bits2);

// compare bits with (bits2 << pos)
extern int bits_cmp_words(ConstBitsP bits, uint64_t pos, ConstBitsP bits2);

//
// Shift, interleave, reverse
//

// removes flanking bits on boths sides, shrinking bits (divon)
void bits_remove_flanking (BitsP bits, uint64_t lsb_flanking, uint64_t msb_flanking);

// shortens an array to a certain number of bits (divon)
void bits_truncate (BitsP bits, uint64_t new_num_of_bits);

// get number of bits in an array, excluding trailing zeros (divon)
uint64_t bits_effective_length (BitsP bits);

// create words - if the word is not aligned to the bitmap word boundaries, and hence spans 2 bitmap words, 
// we take the MSb's from the left word and the LSb's from the right word, to create shift_1 (divon)
static inline uint64_t _bits_combined_word (uint64_t word_a, uint64_t word_b, uint8_t shift)
{
    return shift ? ( (word_a >> shift)       // first  word's MSb are the LSb of our combined word
                   | (word_b << (64-shift))) // second word's LSb are the MSb of our combined word
                 : word_a;
}         

// calculate the number of bits that are different between two Bits' at arbitrary positions (divon)
// note: static inline because in the tight loop of aligner
static inline uint32_t bits_hamming_distance (ConstBitsP bits_1, uint64_t index_1, 
                                              ConstBitsP bits_2, uint64_t index_2, 
                                              uint64_t len)
{
    const uint64_t *words_1 = &bits_1->words[index_1 >> 6];
    uint8_t shift_1 = index_1 & bitmask64(6); // word 1 contributes (64-shift) most-significant bits, word 2 contribute (shift) least significant bits
    uint64_t *after_1 = bits_1->words + bits_1->nwords;

    const uint64_t *words_2 = &bits_2->words[index_2 >> 6];
    uint8_t shift_2 = index_2 & bitmask64(6); // word 1 contributes (64-shift) most-significant bits, word 2 contribute (shift) least significant bits
    uint64_t *after_2 = bits_2->words + bits_2->nwords;

    uint64_t word=0;
    uint32_t nonmatches=0; 
    uint32_t nwords = roundup_bits2words64 (len);

    for (uint32_t i=0; i < nwords; i++) {
        uint64_t word_1 = _bits_combined_word (words_1[i], (&words_1[i+1] < after_1 ? words_1[i+1] : 0), shift_1);
        uint64_t word_2 = _bits_combined_word (words_2[i], (&words_2[i+1] < after_2 ? words_2[i+1] : 0), shift_2);
        word = word_1 ^ word_2; // xor the words - resulting in 1 in each position they differ and 0 where they're equal

        // we count the number of different bits. note that identical nucleotides results in 2 equal bits while
        // different nucleotides results in 0 or 1 equal bits (a 64 bit word contains 32 nucleotides x 2 bit each)
        nonmatches += __builtin_popcountll (word); // note: expected to be _mm_popcnt_u64 with SSE4.2
    }

    // remove non-matches due to the unused part of the last word
    if (len % 64)
        nonmatches -= __builtin_popcountll (word & ~bitmask64 (len % 64));

    return nonmatches; // this is the "hamming distance" between the two Bits' - number of non-matches
}

// iterator on bases (2 bits)
#define BASE_ITER_INIT(bits, base_i, n_bases, is_fwd)       \
    uint64_t _base_i = (is_fwd) ? (base_i) : ((base_i) + (n_bases) - 1); \
    uint64_t *_word_p = &(bits)->words[_base_i >> 5];       \
    uint64_t _word = *_word_p;                              \
    int8_t _next_2bit = _base_i & 31;                       \
    bool _is_fwd __attribute__((unused)) = (is_fwd); 

#define BASE_NEXT_FWD ({                                    \
    if (_next_2bit == 0) _word = *_word_p;                  \
    uint8_t base = (_word >> (_next_2bit * 2)) & 3;         \
    if (++_next_2bit == 32) {                               \
        _next_2bit = 0;                                     \
        _word_p++;                                          \
    }                                                       \
    base; })

#define BASE_NEXT_REVCOMP ({                                \
    if (_next_2bit == 31) _word = *_word_p;                 \
    uint8_t base = 3 - ((_word >> (_next_2bit * 2)) & 3);   \
    if (--_next_2bit == -1) {                               \
        _next_2bit = 31;                                    \
        _word_p--;                                          \
    }                                                       \
    base; })

#define BASE_NEXT (_is_fwd ? BASE_NEXT_FWD : BASE_NEXT_REVCOMP)