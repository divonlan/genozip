// ------------------------------------------------------------------
//   bit_array.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//   Copyright claimed on additions and modifications vs public domain.
//
// a module for handling arrays of 2-bit elements, partially based on public domain code here: 
// https://github.com/noporpoise/BitArray/. 

#ifndef BIT_ARRAY_INCLUDED
#define BIT_ARRAY_INCLUDED

#include "genozip.h"

//------------------------------------------------------------------
// Macros 
//------------------------------------------------------------------

// trailing_zeros is number of least significant zeros
// leading_zeros is number of most significant zeros
#ifndef __GNUC__
  #define trailing_zeros(x) ({ __typeof(x) _r; _BitScanReverse64(&_r, x); _r; })
  #define leading_zeros(x) ({ __typeof(x) _r; _BitScanForward64(&_r, x); _r; })
#else
  #define trailing_zeros(x) ((x) ? (__typeof(x))__builtin_ctzll(x) : (__typeof(x))sizeof(x)*8)
  #define leading_zeros(x) ((x) ? (__typeof(x))__builtin_clzll(x) : (__typeof(x))sizeof(x)*8)
#endif

// Get index of top set bit. If x is 0 return nbits
#define top_set_bit(x) ((x) ? sizeof(x)*8-leading_zeros(x)-1 : sizeof(x)*8)

#define roundup_bits2bytes(bits)   (((bits)+7)/8)
#define roundup_bits2words32(bits) (((bits)+31)/32)
#define roundup_bits2words64(bits) (((bits)+63)/64)
#define roundup_bits2bytes64(bits) (roundup_bits2words64(bits)*8)

// Round a number up to the nearest number that is a power of two (fixed by Divon)
#define roundup2pow(x) (__builtin_popcountll(x)==1 ? (x) : (1UL << (64 - leading_zeros(x))))

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

#ifdef __cplusplus
extern "C" {
#endif

//
// Structs
//

typedef enum {BITARR_UNALLOCATED=0, BITARR_REGULAR, BITARR_OVERLAY} BitArrayType; // must be identical to BufferType

typedef struct BitArray
{
    // These fields should not be changed or added to, as they map to Buffer
    BitArrayType type;
    word_t *words;            // maps to Buffer->data
    bit_index_t num_of_bits;  // maps to Buffer->param
    word_addr_t num_of_words; // maps to Buffer->len  ; round_up (num_of_bits / 64)
} BitArray;

//
// Basics: Constructor, destructor, get length, resize
//

static inline void bit_array_clear_excess_bits_in_top_word (BitArray* bitarr) // divon
{
  if (bitarr->num_of_bits % 64)
    bitarr->words[bitarr->num_of_words-1] &= bitmask64 (bitarr->num_of_bits % 64); 
}

// Allocate using existing struct
extern BitArray* bit_array_alloc(BitArray* bitarr, bit_index_t nbits);
extern void bit_array_dealloc(BitArray* bitarr);

// Get length of bit array
extern bit_index_t bit_array_length(const BitArray* bit_arr);

extern void LTEN_bit_array (BitArray* bitarr); // divon

//
// Macros
//

//
// Get, set, clear, assign and toggle individual bits
// Macros for fast access -- beware: no bounds checking
//

#define bit_array_get(arr,i)      bitset_get((arr)->words, i)
#define bit_array_set(arr,i)      bitset_set((arr)->words, i)
#define bit_array_clear(arr,i)    bitset_del((arr)->words, i)
#define bit_array_toggle(arr,i)   bitset_tgl((arr)->words, i)
// c must be 0 or 1
#define bit_array_assign(arr,i,c) bitset_cpy((arr)->words,i,c)

// c must be 0,1,2 or 3 ; i must be even
#define bit_array_assign2(arr,i,c) bitset_cpy2((arr)->words,i,c) // divon

#define bit_array_len(arr) ((arr)->num_of_bits)

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

// Toggle all the bits in a region
extern void bit_array_toggle_region(BitArray* bitarr, bit_index_t start, bit_index_t len);

//
// Set, clear and toggle all bits at once
//

// Set all bits in this array to 1
extern void bit_array_set_all(BitArray* bitarr);

// Set all bits in this array to 0
extern void bit_array_clear_all(BitArray* bitarr);

// Set all 1 bits to 0, and all 0 bits to 1
extern void bit_array_toggle_all(BitArray* bitarr);

//
// Get / set a word of a given size
//

// First bit is in the least significant bit position
// start index must be within the range of the bit array (0 <= x < length)
extern uint64_t bit_array_get_word64(const BitArray* bitarr, bit_index_t start);
extern uint32_t bit_array_get_word32(const BitArray* bitarr, bit_index_t start);
extern uint16_t bit_array_get_word16(const BitArray* bitarr, bit_index_t start);
extern uint8_t  bit_array_get_word8(const BitArray* bitarr, bit_index_t start);
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

// Takes a char array to write to.  `str` must be bitarr->num_of_bits+1 in
// length. Terminates string with '\0'
extern char* bit_array_to_str(const BitArray* bitarr, char* str);
extern char* bit_array_to_str_rev(const BitArray* bitarr, char* str);

// Get a string representations for a given region, using given on/off
// characters.
// Note: does not nul-terminate
extern void bit_array_to_substr(const BitArray* bitarr,
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

// Copy a BitArray struct and the data it holds - returns pointer to new object
#define bit_array_dup	bit_array_clone
BitArray* bit_array_clone(const BitArray* bitarr);

// Copy bits from one array to another
// Note: use MACRO bit_array_copy
// Destination and source can be the same bit_array and
// src/dst regions can overlap
extern void bit_array_copy(BitArray* dst, bit_index_t dstindx,
                           const BitArray* src, bit_index_t srcindx,
                           bit_index_t length);

// copy all of src to dst. dst is resized to match src.
extern void bit_array_copy_all(BitArray* dst, const BitArray* src);

// for each 2 bits in the src array, the dst array will contain those 2 bits in the reverse
// position, as well as transform them 00->11 11->00 01->10 10->01
// works on arrays with full words
extern void bit_array_reverse_complement_all (BitArray *dst, const BitArray *src, bit_index_t src_start_base, bit_index_t max_num_bases);

extern void bit_array_overlay (BitArray *overlaid_bitarr, BitArray *regular_bitarr, bit_index_t start, bit_index_t num_of_bits);

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
// Number of bytes of data is: (int)((num_of_bits + 7) / 8)
//

//
// Generally useful functions
//

// Generalised 'binary to string' function
// Adds bits to the string in order of lsb to msb
// e.g. 0b11010 (26 in decimal) would come out as "01011"
char* bit_array_word2str(const void *ptr, size_t num_of_bits, char *str);

// Same as above but in reverse
char* bit_array_word2str_rev(const void *ptr, size_t num_of_bits, char *str);

#ifdef __cplusplus
}
#endif

#endif
