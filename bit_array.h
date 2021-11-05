// ------------------------------------------------------------------
//   bit_array.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt
//   Copyright claimed on additions and modifications vs public domain.
//
// a module for handling arrays of 2-bit elements, partially based on BitArray, a public domain code located in https://github.com/noporpoise/BitArray/. 
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

// Round a number up to the nearest number that is a power of two (fixed by Divon)
#define roundup2pow(x) ((__builtin_popcountll(x)==1) ? (x) : ((__typeof(x))1 << (64 - leading_zeros(x))))

#define rot32(x,r) (((x)<<(r)) | ((x)>>(32-(r))))
#define rot64(x,r) (((x)<<(r)) | ((x)>>(64-(r))))

// need to check for length == 0, undefined behaviour if uint64_t >> 64 etc
#define bitmask(nbits,type) ((nbits) ? ~(type)0 >> (sizeof(type)*8-(nbits)): (type)0)
#define bitmask32(nbits) bitmask(nbits,uint32_t)
#define bitmask64(nbits) bitmask(nbits,uint64_t)

// A possibly faster way to combine two words with a mask
//#define bitmask_merge(a,b,abits) ((a & abits) | (b & ~abits))
#define bitmask_merge(a,b,abits) (b ^ ((a ^ b) & abits))

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
#define bitset2_set(arr,wrd,idx)     ((arr)[wrd] |=  _TYPESHIFT(arr,1,idx))
#define bitset2_del(arr,wrd,idx)     ((arr)[wrd] &=~ _TYPESHIFT(arr,1,idx))
#define bitset2_tgl(arr,wrd,idx)     ((arr)[wrd] ^=  _TYPESHIFT(arr,1,idx))
#define bitset2_or(arr,wrd,idx,bit)  ((arr)[wrd] |=  _TYPESHIFT(arr,bit,idx))
#define bitset2_xor(arr,wrd,idx,bit) ((arr)[wrd]  = ~((arr)[wrd] ^ (~_TYPESHIFT(arr,bit,idx))))
#define bitset2_and(arr,wrd,idx,bit) ((arr)[wrd] &= (_TYPESHIFT(arr,bit,idx) | ~_TYPESHIFT(arr,1,idx)))
#define bitset2_cpy(arr,wrd,idx,bit) ((arr)[wrd]  = ((arr)[wrd] &~ _TYPESHIFT(arr,1,idx)) | _TYPESHIFT(arr,bit,idx))

// divon - copy 2 bits - idx must be an even number
#define bitset2_cpy2(arr,wrd,idx,bit2) ((arr)[wrd] = ((arr)[wrd] &~ _TYPESHIFT(arr,0x3,idx)) | _TYPESHIFT(arr,bit2,idx))

//
// Auto detect size of type from pointer
//
#define bitset_wrd(arr,pos) bitsetX_wrd(sizeof(*(arr))*8,pos)
#define bitset_idx(arr,pos) bitsetX_idx(sizeof(*(arr))*8,pos)
#define bitset_op(func,arr,pos)      func(arr, bitset_wrd(arr,pos), bitset_idx(arr,pos))
#define bitset_op2(func,arr,pos,bit) func(arr, bitset_wrd(arr,pos), bitset_idx(arr,pos), bit)

// Auto-detect type size: bit functions
#define bitset_get(arr,pos)     bitset_op(bitset2_get, arr, pos)
#define bitset_set(arr,pos)     bitset_op(bitset2_set, arr, pos)
#define bitset_del(arr,pos)     bitset_op(bitset2_del, arr, pos)
#define bitset_tgl(arr,pos)     bitset_op(bitset2_tgl, arr, pos)
#define bitset_or(arr,pos,bit)  bitset_op2(bitset2_or, arr, pos, bit)
#define bitset_xor(arr,pos,bit) bitset_op2(bitset2_xor, arr, pos, bit)
#define bitset_and(arr,pos,bit) bitset_op2(bitset2_and, arr, pos, bit)
#define bitset_cpy(arr,pos,bit) bitset_op2(bitset2_cpy, arr, pos, bit)
#define bitset_cpy2(arr,pos,bit2) bitset_op2(bitset2_cpy2, arr, pos, bit2)

// Clearing a word does not return a meaningful value
#define bitset_clear_word(arr,pos) ((arr)[bitset_wrd(arr,pos)] = 0)


/*
 * Byteswapping
 */

/* clang uses these to check for features */
#ifndef __has_feature
#define __has_feature(x) 0
#endif

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

/* GCC versions < 4.3 do not have __builtin_bswapX() */
#if ( defined(__clang__) && !__has_builtin(__builtin_bswap64) ) ||             \
    ( !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__) &&   \
      ( (__GNUC__ < 4)  || (__GNUC__ == 4 && __GNUC_MINOR__ < 3)) )
  #define byteswap64(x) ( (((uint64_t)(x) << 56))                       | \
                          (((uint64_t)(x) << 40) & 0xff000000000000ULL) | \
                          (((uint64_t)(x) << 24) & 0xff0000000000ULL)   | \
                          (((uint64_t)(x) <<  8) & 0xff00000000ULL)     | \
                          (((uint64_t)(x) >>  8) & 0xff000000ULL)       | \
                          (((uint64_t)(x) >> 24) & 0xff0000ULL)         | \
                          (((uint64_t)(x) >> 40) & 0xff00ULL)           | \
                          (((uint64_t)(x) >> 56)) )

  #define byteswap32(x) ( (((uint32_t)(x) << 24))                       | \
                          (((uint32_t)(x) <<  8) & 0xff0000U)           | \
                          (((uint32_t)(x) >>  8) & 0xff00U)             | \
                          (((uint32_t)(x) >> 24)) )

  /* uint16_t type might be bigger than 2 bytes, so need to mask */
  #define byteswap16(x) ( (((uint16_t)(x) & 0xff) << 8) | \
                          (((uint16_t)(x) >> 8) & 0xff) )
#else
  #define byteswap64(x) __builtin_bswap64(x)
  #define byteswap32(x) __builtin_bswap64(x)
  #define byteswap16(x) __builtin_bswap64(x)
#endif

//------------------------------------------------------------------

// 64 bit words
typedef uint64_t word_t, word_addr_t, bit_index_t;
typedef uint8_t word_offset_t; // Offset within a 64 bit word

#define BIT_INDEX_MIN 0
#define BIT_INDEX_MAX (~(bit_index_t)0)

//
// Structs
//

typedef enum {BITARR_UNALLOCATED=0, BITARR_REGULAR, BITARR_OVERLAY} BitArrayType; // must be identical to BufferType

typedef struct BitArray
{
    // These fields should not be changed or added to, as they map to Buffer
    BitArrayType type;
    word_t *words;      // maps to Buffer->data
    bit_index_t nbits;  // maps to Buffer->param
    word_addr_t nwords; // maps to Buffer->len  ; round_up (nbits / 64)
} BitArray;

//
// Basics: Constructor, destructor, get length, resize
//

static inline void bit_array_clear_excess_bits_in_top_word (BitArray* bitarr) // divon
{
  if (bitarr->nbits % 64)
    bitarr->words[bitarr->nwords-1] &= bitmask64 (bitarr->nbits % 64); 
}

static inline void ASSERT_excess_bits_are_0 (BitArray* bitarr) // divon
{
  ASSERT0 (!(bitarr->nbits & 0x3f) || !(bitarr->words[bitarr->nwords-1] & ~bitmask64 (bitarr->nbits & 0x3f)),
           "Expecting excess bits in top word of bit array to be 0");
}


// Allocate using existing struct
extern BitArray bit_array_alloc_do (bit_index_t nbits, bool clear, const char *func, uint32_t code_line);
#define bit_array_alloc(nbits, clear) bit_array_alloc_do ((nbits), (clear), __FUNCTION__, __LINE__)

extern void bit_array_realloc_do (BitArray *bitarr, bit_index_t nbits, bit_index_t low_level_nbits, bool clear, const char *func, uint32_t code_line);
#define bit_array_realloc(bitarr, nbits, low_level_nbits, clear) bit_array_realloc_do ((bitarr), (nbits), (low_level_nbits), (clear), __FUNCTION__, __LINE__)

extern void bit_array_free (BitArray* bitarr);

// Get length of bit array
extern bit_index_t bit_array_length (const BitArray* bit_arr);

extern void LTEN_bit_array (BitArray* bitarr); // divon

//
// Macros
//

//
// Get, set, clear, assign and toggle individual bits
// Macros for fast access -- beware: no bounds checking
//

#define bit_array_get(arr,i)      bitset_get((arr)->words, i)
#define bit_array_set(arr,i)      do { BitArray *arr2=arr; /* avoid aliasing */ ; bitset_set((arr2)->words, i); } while(0)
#define bit_array_clear(arr,i)    bitset_del((arr)->words, i)
#define bit_array_toggle(arr,i)   bitset_tgl((arr)->words, i)
// c must be 0 or 1
#define bit_array_assign(arr,i,c) bitset_cpy((arr)->words,i,c)

// c must be 0,1,2 or 3 ; i must be even
#define bit_array_assign2(arr,i,c) bitset_cpy2((arr)->words,i,c) // divon

#define bit_array_len(arr) ((arr)->nbits)

//
// Get, set, clear, assign and toggle individual bits
// "Safe": use assert() to check bounds
//

// Get the value of a bit (returns 0 or 1)
extern char bit_array_get_bit(const BitArray* bitarr, bit_index_t b);
extern void bit_array_set_bit(BitArray* bitarr, bit_index_t b);
extern void bit_array_clear_bit(BitArray* bitarr, bit_index_t b);
extern void bit_array_toggle_bit(BitArray* bitarr, bit_index_t b);
// If char c != 0, set bit; otherwise clear bit
extern void bit_array_assign_bit(BitArray* bitarr, bit_index_t b, char c);

//
// Get, set, clear and toggle several bits at once
//

// Get the offsets of the set bits (for offsets start<=offset<end)
// Returns the number of bits set
// It is assumed that dst is at least of length (end-start)
extern bit_index_t bit_array_get_bits(const BitArray* bitarr,
                                      bit_index_t start, bit_index_t end,
                                      bit_index_t* dst);

// Set multiple bits at once.
// e.g. set bits 1, 20 & 31: bit_array_set_bits(bitarr, 3, 1,20,31);
// Note: variable args are of type unsigned int
extern void bit_array_set_bits(BitArray* bitarr, size_t n, ...);

// Clear multiple bits at once.
// e.g. clear bits 1, 20 & 31: bit_array_clear_bits(bitarr, 3, 1,20,31);
// Note: variable args are of type unsigned int
extern void bit_array_clear_bits(BitArray* bitarr, size_t n, ...);

// Toggle multiple bits at once
// e.g. toggle bits 1, 20 & 31: bit_array_toggle_bits(bitarr, 3, 1,20,31);
// Note: variable args are of type unsigned int
extern void bit_array_toggle_bits(BitArray* bitarr, size_t n, ...);

//
// Set, clear and toggle all bits in a region
//

// Set all the bits in a region
extern void bit_array_set_region(BitArray* bitarr, bit_index_t start, bit_index_t len);

// Clear all the bits in a region
#define bit_array_clear_region(bitarr,start,len) bit_array_clear_region_do (bitarr, start, len, __FUNCTION__, __LINE__)
extern void bit_array_clear_region_do (BitArray* bitarr, bit_index_t start, bit_index_t len, const char *func, unsigned code_line);

extern uint32_t bit_array_hamming_distance (const BitArray *bitarr1, bit_index_t index1, const BitArray *bitarr2, bit_index_t index2, bit_index_t len);

//
// Set, clear and toggle all bits at once
//

// Set all bits in this array to 1
extern void bit_array_set_all(BitArray* bitarr);

// Set all bits in this array to 0
extern void bit_array_clear_all(BitArray* bitarr);

//
// Get and set words (low level -- no bounds checking)
//
static inline word_t _get_word(const BitArray* bitarr, bit_index_t start)
{
    word_addr_t word_index = bitset64_wrd(start);
    word_offset_t word_offset = bitset64_idx(start);

    word_t result = bitarr->words[word_index] >> word_offset;

    word_offset_t bits_taken = WORD_SIZE - word_offset;

    // word_offset is now the number of bits we need from the next word
    // Check the next word has at least some bits
    if(word_offset > 0 && start + bits_taken < bitarr->nbits)
        result |= bitarr->words[word_index+1] << (WORD_SIZE - word_offset);

    return result;
}

//
// Get / set a word of a given size
//

// bitarray is interpreted as an integer in Little Endian - i.e. first bit in the bitarray is LSb
#define BIT_ARRAY_GET_WORD(width) \
  static inline uint64_t bit_array_get_word ## width (const BitArray* bitarr, bit_index_t start) \
  { \
    ASSERT (start + width <= bitarr->nbits, "expecting start(%"PRIu64") + " #width " <= bitarr->nbits(%"PRIu64")", start, bitarr->nbits); \
    return (uint ## width ## _t)_get_word(bitarr, start); \
  } 
BIT_ARRAY_GET_WORD(64)
BIT_ARRAY_GET_WORD(32)
BIT_ARRAY_GET_WORD(16)
BIT_ARRAY_GET_WORD(8)

extern uint64_t bit_array_get_wordn(const BitArray* bitarr, bit_index_t start, int n);

// Set 64 bits at once from a particular start position
// Doesn't extend bit array. However it is safe to TRY to set bits beyond the
// end of the array, as long as: `start` is < `bit_array_length(arr)`
extern void bit_array_set_word64(BitArray* bitarr, bit_index_t start, uint64_t word);
extern void bit_array_set_word32(BitArray* bitarr, bit_index_t start, uint32_t word);
extern void bit_array_set_word16(BitArray* bitarr, bit_index_t start, uint16_t word);
extern void bit_array_set_word8(BitArray* bitarr, bit_index_t start, uint8_t byte);
extern void bit_array_set_wordn(BitArray* bitarr, bit_index_t start, uint64_t word, int n);

//
// Number of bits set
//

// Get the number of bits set (hamming weight)
extern bit_index_t bit_array_num_bits_set(const BitArray* bitarr);
extern bit_index_t bit_array_num_bits_set_region(const BitArray* bitarr, bit_index_t start, bit_index_t length); // added by divon

// Get the number of bits not set (length - hamming weight)
extern bit_index_t bit_array_num_bits_cleared(const BitArray* bitarr);

//
// Find indices of set/clear bits
//

// Find the index of the next bit that is set, at or after `offset`
// Returns 1 if a bit is set, otherwise 0
// Index of next set bit is stored in the integer pointed to by result
// If no next bit is set result is not changed
extern char bit_array_find_next_set_bit(const BitArray* bitarr, bit_index_t offset, bit_index_t* result);

// Find the index of the next bit that is NOT set, at or after `offset`
// Returns 1 if a bit is NOT set, otherwise 0
// Index of next zero bit is stored in the integer pointed to by `result`
// If no next bit is zero, value at `result` is not changed
extern char bit_array_find_next_clear_bit(const BitArray* bitarr, bit_index_t offset,
                                 bit_index_t* result);

// Find the index of the previous bit that is set, before offset.
// Returns 1 if a bit is set, otherwise 0
// Index of previous set bit is stored in the integer pointed to by `result`
// If no previous bit is set result is not changed
extern char bit_array_find_prev_set_bit(const BitArray* bitarr, bit_index_t offset,
                                 bit_index_t* result);

// Find the index of the previous bit that is NOT set, before offset.
// Returns 1 if a bit is clear, otherwise 0
// Index of previous zero bit is stored in the integer pointed to by `result`
// If no previous bit is zero result is not changed
extern char bit_array_find_prev_clear_bit(const BitArray* bitarr, bit_index_t offset,
                                   bit_index_t* result);

// Find the index of the first bit that is set.
// Returns 1 if a bit is set, otherwise 0
// Index of first set bit is stored in the integer pointed to by `result`
// If no bit is set result is not changed
extern char bit_array_find_first_set_bit(const BitArray* bitarr, bit_index_t* result);

// Find the index of the first bit that is NOT set.
// Returns 1 if a bit is clear, otherwise 0
// Index of first zero bit is stored in the integer pointed to by `result`
// If no bit is zero result is not changed
extern char bit_array_find_first_clear_bit(const BitArray* bitarr, bit_index_t* result);

// Find the index of the last bit that is set.
// Returns 1 if a bit is set, otherwise 0
// Index of last set bit is stored in the integer pointed to by `result`
// If no bit is set result is not changed
extern char bit_array_find_last_set_bit(const BitArray* bitarr, bit_index_t* result);

// Find the index of the last bit that is NOT set.
// Returns 1 if a bit is clear, otherwise 0
// Index of last zero bit is stored in the integer pointed to by `result`
// If no bit is zero result is not changed
extern char bit_array_find_last_clear_bit(const BitArray* bitarr, bit_index_t* result);


//
// String and printing methods
//

// Construct a BitArray from a string.
extern void bit_array_from_str(BitArray* bitarr, const char* bitstr);

// Construct a BitArray from a substring with given on and off characters.
extern void bit_array_from_substr(BitArray* bitarr, bit_index_t offset,
                           const char* str, size_t len,
                           const char *on, const char *off, char left_to_right);

// Takes a char array to write to.  `str` must be bitarr->nbits+1 in
// length. Terminates string with '\0'
extern char* bit_array_to_str(const BitArray* bitarr, char* str);
extern char* bit_array_to_str_rev(const BitArray* bitarr, char* str);

// Get a string representations for a given region, using given on/off
// characters.
// Note: does not nul-terminate
extern void bit_array_to_substr (const BitArray* bitarr,
                                 bit_index_t start, bit_index_t length,
                                 char* str, char on, char off, char left_to_right);

// Print this array to a file stream.  Prints '0's and '1'.  Doesn't print
// newline.
extern void bit_array_print_do (const BitArray* bitarr, const char *msg, FILE *file);
#define bit_array_print(bitarr) bit_array_print_do (bitarr, #bitarr, info_stream)

extern void bit_array_print_binary_word_do (word_t word, const char *msg, FILE *file);
#define bit_array_print_binary_word(word) bit_array_print_binary_word_do (word, #word, info_stream)

// Print a string representations for a given region, using given on/off
// characters. Reverse prints from highest to lowest -- this is useful for
// printing binary numbers
extern void bit_array_print_substr (const char *msg, const BitArray* bitarr,
                                    bit_index_t start, bit_index_t length, FILE *file);

// Print bit array as hex
extern size_t bit_array_print_hex (const BitArray* bitarr,
                                   bit_index_t start, bit_index_t length,
                                   FILE* fout, char uppercase);

//
// Clone and copy
//

// Copy bits from one array to another
// Note: use MACRO bit_array_copy
// Destination and source can be the same bit_array and
// src/dst regions can overlap
extern void bit_array_copy(BitArray* dst, bit_index_t dstindx,
                           const BitArray* src, bit_index_t srcindx,
                           bit_index_t length);

// concatenate src at the end of dst (divon)
extern void bit_array_concat (BitArray *base, const BitArray *add, unsigned additional_concats_expected);

// for each 2 bits in the src array, the dst array will contain those 2 bits in the reverse
// position, as well as transform them 00->11 11->00 01->10 10->01
// works on arrays with full words
extern void bit_array_reverse_complement_all (BitArray *dst, const BitArray *src, bit_index_t src_start_base, bit_index_t max_num_bases);

extern void bit_array_overlay (BitArray *overlaid_bitarr, BitArray *regular_bitarr, bit_index_t start, bit_index_t nbits);

//
// Logic operators
//

// BIT_ARRAYs can all be different or the same object
// dest array will be resized if it is too short
//
extern void bit_array_and(BitArray* dest, const BitArray* src1, const BitArray* src2);
extern void bit_array_or (BitArray* dest, const BitArray* src1, const BitArray* src2);
extern void bit_array_xor(BitArray* dest, const BitArray* src1, const BitArray* src2);
extern void bit_array_not(BitArray* dest, const BitArray* src);

//extern void bit_array_or_with (BitArray *dst, bit_index_t dst_start_bit, BitArray *src , bit_index_t src_start_bit, bit_index_t len);

//
// Comparisons
//

// Note: (bit_array_cmp(a,b) == 0) <=> (bit_array_cmp_big_endian(a,b) == 0)

// comparison functions return:
//   1 iff bitarr1 > bitarr2
//   0 iff bitarr1 == bitarr2
//  -1 iff bitarr1 < bitarr2

// Compare two bit arrays by value stored, with index 0 being the Least
// Significant Bit (LSB). Arrays do not have to be the same length.
// Example: ..0101 (5) > ...0011 (3) [index 0 is LSB at right hand side]
extern int bit_array_cmp(const BitArray* bitarr1, const BitArray* bitarr2);

// Compare two bit arrays by value stored, with index 0 being the Most
// Significant Bit (MSB). Arrays do not have to be the same length.
// Example: 10.. > 01.. [index 0 is MSB at left hand side]
extern int bit_array_cmp_big_endian(const BitArray* bitarr1, const BitArray* bitarr2);

// compare bitarr with (bitarr2 << pos)
extern int bit_array_cmp_words(const BitArray *bitarr,
                        bit_index_t pos, const BitArray *bitarr2);

//
// Shift, interleave, reverse
//

// removes flanking bits on boths sides, shrinking bitarr (divon)
void bit_array_remove_flanking (BitArray* bitarr, bit_index_t lsb_flanking, bit_index_t msb_flanking);

// shortens an array to a certain number of bits (divon)
void bit_array_truncate (BitArray* bitarr, bit_index_t new_num_of_bits);

// Shift array left/right.  
void bit_array_shift_right_shrink (BitArray* bitarr, bit_index_t shift_dist);
void bit_array_shift_right(BitArray* bitarr, bit_index_t shift_dist, char fill);
void bit_array_shift_left (BitArray* bitarr, bit_index_t shift_dist, char fill);

// Cyclic shift
void bit_array_cycle_right(BitArray* bitarr, bit_index_t dist);
void bit_array_cycle_left (BitArray* bitarr, bit_index_t dist);

// Reverse the whole array or part of it
void bit_array_reverse(BitArray* bitarr);
void bit_array_reverse_region(BitArray* bitarr, bit_index_t start, bit_index_t len);

//
// Read/Write bit_array to a file
//
// File format is [8 bytes: for number of elements in array][data]
// Number of bytes of data is: (int)((nbits + 7) / 8)
//

//
// Generally useful functions
//

// Generalised 'binary to string' function
// Adds bits to the string in order of lsb to msb
// e.g. 0b11010 (26 in decimal) would come out as "01011"
char* bit_array_word2str(const void *ptr, size_t nbits, char *str);

// Same as above but in reverse
char* bit_array_word2str_rev(const void *ptr, size_t nbits, char *str);

// get number of bits in an array, excluding trailing zeros (divon)
bit_index_t bit_array_effective_length (BitArray *bitarr);

