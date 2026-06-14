// ------------------------------------------------------------------
//   bits.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.
//   Copyright claimed on additions and modifications vs public domain.
//
// a module for handling bit arrays, partially based on: https://github.com/noporpoise/bit_array/ which says:
// "This software is in the Public Domain. That means you can do whatever you like with it. That includes being used in proprietary products 
// without attribution or restrictions. There are no warranties and there may be bugs."

#include <stdarg.h>
#include "genozip.h"
#include "endianness.h"
#include "bits.h"
#include "buffer.h"

// Used in debugging
#ifdef DEBUG

// Mostly used for debugging
typedef struct { char s[WORD_SIZE+1]; } WordStr;
static inline WordStr _print_word (uint64_t word)
{
    WordStr w;
    w.s[WORD_SIZE] = 0;

    for (word_offset_t i=0; i < WORD_SIZE; i++)
        w.s[i] = ((word >> i) & (uint64_t)0x1) == 0 ? '0' : '1';

    return w;
}

#define DEBUG_VALIDATE(a) validate_bits((a), __FUNCLINE)
static void validate_bits (ConstBitsP bits, rom file, int lineno)
{
    // Verify that its allocated
    ASSERT (bits->type != BUF_UNALLOCATED, "[%s:%i] Bits is not allocated", file, lineno);

    // Check num of words is correct
    uint64_t num_words = roundup_bits2words64 (bits->nbits);
    ASSERT (num_words == bits->nwords, "[%s:%i] num of words wrong [bits: %i, expected nwords: %i, actual nwords: %i]", 
            file, lineno, (int)bits->nbits, (int)num_words, (int)bits->nwords);

    // Check top word is masked (only if not overlayed - the unused bits of top word don't belong to this bit array and might be used eg by another bit array in genome.ref/genome.is_set)
    if (bits->type == BUF_REGULAR && bits->nwords) {
      uint64_t tw = bits->nwords - 1;
      uint64_t top_bits = bits_in_top_word (bits->nbits);

      ASSERT (top_bits == 64 ||  bits->words[tw] <= bitmask64(top_bits), "[%s:%i] Expected %i bits in top word[%i] (the rest should be 0) but word=%s\n", 
              file, lineno, (int)top_bits, (int)tw, _print_word(bits->words[tw]).s);
    }
}

#else
    #define DEBUG_VALIDATE(a)
#endif

// Reverse a word
static inline uint64_t _reverse_word (uint64_t word)
{
#if defined(__clang__) && defined(__aarch64__) // 4-5 clock cycles on ARM (bitreverse64 does not exist natively on x86)
    word = __builtin_bitreverse64 (word); // reverse bit order

#else // 5 clock cycles on Intel, 8-10 on ARM
    word = ((word >> 1) & 0x5555555555555555ULL) | ((word & 0x5555555555555555ULL) << 1); // Swap adjacent bits
    word = ((word >> 2) & 0x3333333333333333ULL) | ((word & 0x3333333333333333ULL) << 2); // Swap adjacent pairs
    word = ((word >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((word & 0x0F0F0F0F0F0F0F0FULL) << 4); // Swap nibbles 
    word = __builtin_bswap64 (word); // reverse the 8 bytes of the word 
#endif
 
    return word; // reverses bit values - effectively A(00)⇔T(11) C(01)⇔G(10)
}
    
// Wrap around
static inline uint64_t _get_word_cyclic (ConstBitsP bits, uint64_t start)
{
    uint64_t word = _get_word(bits, start);

    uint64_t bits_taken = bits->nbits - start;

    if (bits_taken < WORD_SIZE) {
        word |= (bits->words[0] << bits_taken);

        if (bits->nbits < (uint64_t)WORD_SIZE)
            // Mask word to prevent repetition of the same bits
            word = word & bitmask64(bits->nbits);
    }

  return word;
}

// Wrap around
static inline void _set_word_cyclic (BitsP bits, uint64_t start, uint64_t word)
{
    _set_word(bits, start, word);

    uint64_t bits_set = bits->nbits - start;

    if (bits_set < WORD_SIZE && start > 0) {
        word >>= bits_set;

        // Prevent overwriting the bits we've just set
        // by setting 'start' as the upper bound for the number of bits to write
        word_offset_t bits_remaining = MIN_(WORD_SIZE - bits_set, start);
        uint64_t mask = bitmask64(bits_remaining);

        bits->words[0] = bitmask_merge(word, bits->words[0], mask);
    }
}

// allocates a STANDALONE bit array
Bits bits_alloc_do (uint64_t nbits, bool clear, FUNCLINE)
{
    Bits bits = { .type   = BITS_STANDALONE, // standalone Bits that is not part of a Buffer
                  .nbits  = nbits,
                  .nwords = roundup_bits2words64(nbits),
                  .words  = buf_low_level_malloc (roundup_bits2bytes64(nbits), clear, func, code_line) };

    // zero the bits in the top word that are beyond nbits (if not already cleared)
    if (!clear) bits_clear_excess_bits_in_top_word (&bits,  false);
    
    return bits;
}

// set all the bits in a region
void bits_set_region (BitsP bits, uint64_t start, uint64_t len)
{
    if (!len) return;

    ASSERT (start + len <= bits->nbits, "Expecting: start=%"PRId64" + len=%"PRId64" <= nbits=%"PRId64"",
            start, len, bits->nbits); 

    uint64_t first_word   = bitset64_wrd(start);
    uint64_t last_word    = bitset64_wrd(start + len - 1);
    word_offset_t foffset = bitset64_idx(start);
    word_offset_t loffset = bitset64_idx(start + len - 1);

    if (__builtin_expect (first_word == last_word, false)) {
        uint64_t mask = bitmask64(len) << foffset;
        bits->words[first_word] |= mask; 
    }
    
    else {
        bits->words[first_word] |= ~bitmask64(foffset);     // set first word
        
        for (uint64_t i=first_word + 1; i < last_word; i++) // set whole words
            bits->words[i] = UINT64_MAX;

        bits->words[last_word] |= bitmask64(loffset+1);     // set last word
    }

    DEBUG_VALIDATE (bits);
}

// Clear all the bits in a region
void bits_clear_region_do (BitsP bits, uint64_t start, uint64_t len, FUNCLINE)
{
    if (!len) return; // nothing to do 

    ASSERT (start + len <= bits->nbits, "called from %s:%u: Expecting: start=%"PRId64" + len=%"PRId64" <= nbits=%"PRId64,
            func, code_line, start, len, bits->nbits); // divon fixed bug

    uint64_t first_word   = bitset64_wrd(start);
    uint64_t last_word    = bitset64_wrd(start + len - 1);
    word_offset_t foffset = bitset64_idx(start);
    word_offset_t loffset = bitset64_idx(start + len - 1);

    if (__builtin_expect (first_word == last_word, false)) {
        uint64_t mask = bitmask64(len) << foffset;
        bits->words[first_word] &= ~mask;
    }
    
    else {
        bits->words[first_word] &= bitmask64(foffset);      // clear first word

        for (uint64_t i=first_word + 1; i < last_word; i++) // clear whole words
            bits->words[i] = (uint64_t)0;

        bits->words[last_word] &= ~bitmask64(loffset+1);    // clear last word
    }

    DEBUG_VALIDATE (bits);
}

//
// Number of bits set
//

// true if all bits in bit array are set (divon)
bool bits_is_fully_set (ConstBitsP bits)
{
    if (!bits->nbits) return true; // trivially true

    // all words but last one
    for (uint64_t i=0; i < bits->nwords-1; i++)
        if (bits->words[i] != UINT64_MAX) return false;

    // last word might be partial
    word_offset_t bits_active = bits_in_top_word (bits->nbits);    
    word_offset_t num_of_bits_set = __builtin_popcountll (bits->words[bits->nwords-1] & bitmask64 (bits_active));

    return num_of_bits_set == bits_active;
}

// true if all bits in bit array are clear (divon)
bool bits_is_fully_clear (ConstBitsP bits)
{
    if (!bits->nbits) return true; // trivially true

    // all words but last one
    for (uint64_t i=0; i < bits->nwords - 1; i++)
        if (bits->words[i]) return false;

    return (bits->words[bits->nwords-1] & bitmask64(bits_in_top_word(bits->nbits))) == 0;
}

// Get the number of bits set (hamming weight)
uint64_t bits_num_set_bits (ConstBitsP bits)
{
    if (!bits->nbits) return 0;

    uint64_t num_of_bits_set = 0;

    // all words but last one
    for (uint64_t i=0; i < bits->nwords-1; i++)
        if (bits->words[i])
            num_of_bits_set += __builtin_popcountll (bits->words[i]);

    // last word might be partial
    word_offset_t bits_active = bits_in_top_word (bits->nbits);    
    num_of_bits_set += __builtin_popcountll (bits->words[bits->nwords-1] & bitmask64 (bits_active));
    
    return num_of_bits_set;
}

// added by divon
uint64_t bits_num_set_bits_region (ConstBitsP bits, uint64_t start, uint64_t length)
{
    if (length == 0) return 0;

    ASSERT (start + length <= bits->nbits, "out of range: execpting: start=%"PRIu64" + length=%"PRIu64" <= nbits=%"PRIu64, start, length, bits->nbits);

    uint64_t first_word = bitset64_wrd(start);
    uint64_t last_word  = bitset64_wrd(start+length-1);
    word_offset_t foffset = bitset64_idx(start);

    uint64_t num_of_bits_set = 0;

    if (first_word == last_word) {
        uint64_t mask = bitmask64(length) << foffset;
        num_of_bits_set += __builtin_popcountll (bits->words[first_word] & mask);
    }
    else {
        word_offset_t loffset  = bitset64_idx(start+length-1);
    
        // first word
        num_of_bits_set += __builtin_popcountll (bits->words[first_word] & ~bitmask64(foffset));

        // whole words
        for (uint64_t i = first_word + 1; i < last_word; i++)
            num_of_bits_set += __builtin_popcountll (bits->words[i]);

        // last word
        num_of_bits_set += __builtin_popcountll (bits->words[last_word] & bitmask64(loffset+1));
    }

    return num_of_bits_set;
}

// Get the number of bits not set (1 - hamming weight)
uint64_t bits_num_clear_bits (ConstBitsP bits)
{
    return bits->nbits - bits_num_set_bits(bits);
}



//
// Find indices of set/clear bits
//

// Find the index of the next bit that is set/clear, at or after `offset`
// Returns 1 if such a bit is found, otherwise 0
// Index is stored in the integer pointed to by `result`
// If no such bit is found, value at `result` is not changed
#define _next_bit_func_def(FUNC,GET) \
bool FUNC(ConstBitsP bits, uint64_t offset, uint64_t *result) \
{ \
    ASSERT (offset < bits->nbits, "expecting offset(%"PRId64") < bits->nbits(%"PRId64")", offset, bits->nbits); \
    if (bits->nbits == 0 || offset >= bits->nbits) { return false; } \
    \
    /* Find first word that is greater than zero */ \
    uint64_t i = bitset64_wrd(offset); \
    uint64_t w = GET(bits->words[i]) & ~bitmask64(bitset64_idx(offset)); \
    \
    while (1) { \
        if (w > 0) { \
            uint64_t pos = i * WORD_SIZE + trailing_zeros(w); \
            if (pos < bits->nbits) { *result = pos; return true; } \
            else { return false; } \
        } \
        i++; \
        if (i >= bits->nwords) break; \
        w = GET(bits->words[i]); \
    } \
    \
    return false; \
}

// Find the index of the previous bit that is set/clear, before `offset`.
// Returns 1 if such a bit is found, otherwise 0
// Index is stored in the integer pointed to by `result`
// If no such bit is found, value at `result` is not changed
#define _prev_bit_func_def(FUNC,GET) \
bool FUNC(ConstBitsP bits, uint64_t offset, uint64_t *result) \
{ \
    ASSERT (offset <= bits->nbits, "expecting offset=%"PRIu64" <= nbits=%"PRIu64, offset, bits->nbits); \
    if (bits->nbits == 0 || offset == 0) { return false; } \
    \
    /* Find prev word that is greater than zero */ \
    uint64_t i = bitset64_wrd(offset-1); \
    uint64_t w = GET(bits->words[i]) & bitmask64(bitset64_idx(offset-1)+1); \
    \
    if (w > 0) { *result = (i+1) * WORD_SIZE - leading_zeros(w) - 1; return true; } \
    \
    /* i is unsigned so have to use break when i == 0 */ \
    for (--i; i != UINT64_MAX; i--) { \
        w = GET(bits->words[i]); \
        if (w > 0) { \
            *result = (i+1) * WORD_SIZE - leading_zeros(w) - 1; \
            return true; \
        } \
    } \
    \
    return false; \
}

#define GET_WORD(x) (x)
#define NEG_WORD(x) (~(x))
_next_bit_func_def(bits_find_next_set_bit,  GET_WORD);
_next_bit_func_def(bits_find_next_clear_bit,NEG_WORD);
_prev_bit_func_def(bits_find_prev_set_bit,  GET_WORD);
_prev_bit_func_def(bits_find_prev_clear_bit,NEG_WORD);

// Find the index of the first bit that is set.
// Returns true if a bit is set.
// Index of first set bit is stored in the integer pointed to by result
// If no bits are set, value at `result` is not changed
bool bits_find_first_set_bit (ConstBitsP bits, uint64_t *result)
{
    return bits_find_next_set_bit (bits, 0, result);
}

// same same
bool bits_find_first_clear_bit (ConstBitsP bits, uint64_t *result)
{
    return bits_find_next_clear_bit (bits, 0, result);
}

// Find the index of the last bit that is set.
// Returns 1 if a bit is set, otherwise 0
// Index of last set bit is stored in the integer pointed to by `result`
// If no bits are set, value at `result` is not changed
bool bits_find_last_set_bit (ConstBitsP bits, uint64_t *result)
{
    return bits_find_prev_set_bit (bits, bits->nbits, result);
}

// same same
bool bits_find_last_clear_bit (ConstBitsP bits, uint64_t *result)
{
    return bits_find_prev_clear_bit (bits, bits->nbits, result);
}


//
// Strings and printing
//

// Get a string representations for a given region, using given on/off characters.
char *bits_to_substr (ConstBitsP bits, uint64_t start, uint64_t length, char *str)
{
    if (start >= bits->nbits) 
        length = 0; // nothing to output
    
    else if (start + length > bits->nbits)
        length = bits->nbits - start; // until the end of the Bits, less than the original length
      
    for (uint64_t i=0; i < length; i++) 
        str[i] = bits_get(bits, start + i) ? '1' : '0';

    str[length] = '\0';

    return str;
}

void bits_print_do (ConstBitsP bits, rom msg, FILE *file)
{
    fprintf (file, "%s (nbits=%"PRIu64"): ", msg, bits->nbits);

    for (uint64_t i=0; i < bits->nbits; i++)
        fprintf (file, "%c", bits_get(bits, i) ? '1' : '0');

    fprintf (file, "\n");
}

void bits_print_binary_word_do (uint64_t word, rom msg, FILE *file)
{
    Bits bits = { .nbits=64, .nwords=1, .words = &word };
    bits_print_do (&bits, msg, file);
}

// Print a string representations for a given region, using given on/off characters.
void bits_print_substr (rom msg, ConstBitsP bits,uint64_t start, uint64_t length, FILE *file)
{
    if (!file) file = info_stream;
    
    length = MIN_(length, bits->nbits - start);

    if (msg) fprintf (file, "%s: ", msg);
    for (uint64_t i=0; i < length; i++) 
        fputc (bits_get(bits, start + i) ? '1' : '0', file);

    fputc ('\n', file);
}

// prints a 2bit bases array - A,C,G,T bases
void bits_print_substr_bases (rom msg, ConstBitsP bits, uint64_t start_base, uint64_t n_bases, FILE *file)
{
    if (!file) file = info_stream;

    uint64_t start_bit = start_base * 2;
    uint64_t n_bits = MIN_(n_bases * 2, bits->nbits - start_bit);

    if (msg) fprintf (file, "%s: ", msg);
    for (uint64_t i=0; i < n_bits; i += 2) { 
        char b = bits_get2 (bits, start_bit + i);
        fputc (b==0 ? 'A' : b==1 ? 'C' : b==2 ? 'G' : 'T' , file);
    }

    fputc ('\n', file);
}

// destination and source may be the same bits and src/dst regions may overlap
void bits_copy_do (BitsP dst, uint64_t dstindx, ConstBitsP src, uint64_t srcindx, uint64_t length, FUNCLINE)
{
    DEBUG_VALIDATE(src);

    ASSERT (dstindx + length <= dst->nbits, "called from %s:%u dstindx(%"PRIu64") + length(%"PRIu64") > dst->nbits(%"PRIu64")", func, code_line, dstindx, length, dst->nbits);
    ASSERT (srcindx + length <= src->nbits, "called from %s:%u srcindx(%"PRIu64") + length(%"PRIu64") > src->nbits(%"PRIu64")", func, code_line, srcindx, length, src->nbits);

    // Num of full words to copy
    uint64_t num_of_full_words = length / WORD_SIZE;
    uint64_t i;

    word_offset_t bits_in_last_word = length % WORD_SIZE; 

    // if possible, work left to right (should work faster due to CPU pre-fetching etc)
    if (dst != src || srcindx > dstindx) {
        for (i=0; i < num_of_full_words; i++) {
            uint64_t word = _get_word (src, srcindx + i * WORD_SIZE);
            _set_word (dst, dstindx + i * WORD_SIZE, word);
        }

        if (bits_in_last_word > 0) {
            uint64_t src_word = _get_word (src, srcindx + i * WORD_SIZE);
            uint64_t dst_word = _get_word (dst, dstindx + i * WORD_SIZE);

            uint64_t mask = bitmask64 (bits_in_last_word);
            uint64_t word = bitmask_merge (src_word, dst_word, mask);

            _set_word (dst, dstindx+num_of_full_words * WORD_SIZE, word);
        }
    }

    else {
        // Work right to left
        for (i = 0; i < num_of_full_words; i++) {
            uint64_t word = _get_word (src, srcindx+length-(i+1)*WORD_SIZE);
            _set_word(dst, dstindx + length - (i+1) * WORD_SIZE, word);
        }

        if (bits_in_last_word > 0) {
            uint64_t src_word = _get_word (src, srcindx);
            uint64_t dst_word = _get_word (dst, dstindx);

            uint64_t mask = bitmask64 (bits_in_last_word);
            uint64_t word = bitmask_merge (src_word, dst_word, mask);
            _set_word (dst, dstindx, word);
        }
    }

    // zero top word if copy reaches it
    if (dstindx + length == dst->nbits)
        bits_clear_excess_bits_in_top_word (dst, true);

    DEBUG_VALIDATE(dst);
}

void bits_overlay (BitsP overlaid_bits, BitsP regular_bits, uint64_t start, uint64_t nbits)
{
    ASSERT (start % 64 == 0, "start=%"PRIu64" must be a multiple of 64", start);
    ASSERT (start + nbits <= regular_bits->nbits, "start(%"PRIu64") + nbits(%"PRIu64") <= regular_bits->nbits(%"PRIu64")",
            start, nbits, regular_bits->nbits);

    uint64_t word_i = start / 64;
    *overlaid_bits = (Bits){ .type   = BITS_STANDALONE ,
                             .words  = &regular_bits->words[word_i],
                             .nwords = roundup_bits2words64 (nbits),
                             .nbits  = nbits };
} 

// convert a bit array to a byte array of values 0-1
void bits_bit_to_byte (uint8_t *dst, ConstBitsP src_bits, uint64_t src_bit, uint32_t num_bits)
{
    ASSERT (src_bit + num_bits <= src_bits->nbits, "Expecting src_bit=%"PRIu64" + num_bits=%u) <= nbits=%"PRIu64,
            src_bit, num_bits, src_bits->nbits);

    uint64_t *src_word_p = &src_bits->words[src_bit >> 6];
    uint64_t src_word = *src_word_p;
    uint8_t next_bit = src_bit & 63;

    for (uint32_t i=0; i < num_bits; i++) {
        *dst++ = (src_word >> next_bit) & 1;
        
        if (++next_bit == 64) {
            next_bit = 0;
            src_word_p++;
            src_word = *src_word_p;
        }
    }
}

//
// Reverse -- coords may wrap around
//

// No bounds checking
// length cannot be zero
static void _reverse_region (BitsP bits, uint64_t start, uint64_t length)
{
    uint64_t left = start;
    uint64_t right = (start + length - WORD_SIZE) % bits->nbits;

    while(length >= 2 * WORD_SIZE) {
        // Swap entire words
        uint64_t left_word  = _get_word_cyclic (bits, left);
        uint64_t right_word = _get_word_cyclic (bits, right);

        // reverse words individually
        left_word  = _reverse_word(left_word);
        right_word = _reverse_word(right_word);

        // Swap
        _set_word_cyclic (bits, left, right_word);
        _set_word_cyclic (bits, right, left_word);

        // Update
        left = (left + WORD_SIZE) % bits->nbits;
        right = (right < WORD_SIZE ? right + bits->nbits : right) - WORD_SIZE;
        length -= 2 * WORD_SIZE;
    }

    uint64_t word, rev;

    if (length == 0) 
        return;

    else if (length > WORD_SIZE) {
        // Words overlap
        uint64_t left_word = _get_word_cyclic (bits, left);
        uint64_t right_word = _get_word_cyclic (bits, right);

        rev = _reverse_word (left_word);
        right_word = _reverse_word (right_word);

        // fill left 64 bits with right word rev
        _set_word_cyclic(bits, left, right_word);

        // Now do remaining bits (length is between 1 and 64 bits)
        left += WORD_SIZE;
        length -= WORD_SIZE;

        word = _get_word_cyclic (bits, left);
    }
    
    else {
        word = _get_word_cyclic (bits, left);
        rev = _reverse_word (word);
    }

    rev >>= WORD_SIZE - length;
    uint64_t mask = bitmask64(length);

    word = bitmask_merge(rev, word, mask);

    _set_word_cyclic (bits, left, word);
}

void bits_reverse (BitsP bits)
{
    if (bits->nbits) _reverse_region (bits, 0, bits->nbits);
    
    DEBUG_VALIDATE(bits);
}

// removes flanking bits on boths sides, shrinking bits
void bits_remove_flanking (BitsP bits, uint64_t lsb_flanking, uint64_t msb_flanking) // added by divon
{
    DEBUG_VALIDATE (bits); // catch a bug

    uint64_t cpy_length = bits->nbits - lsb_flanking;
    bits_copy (bits, 0, bits, lsb_flanking, cpy_length);

    bits->nbits -= lsb_flanking + msb_flanking;
    bits->nwords = roundup_bits2words64 (bits->nbits);   
}

// shortens an array to a certain number of bits (divon)
void bits_truncate (BitsP bits, uint64_t new_num_of_bits)
{
    ASSERT (new_num_of_bits <= bits->nbits, "expecting new_num_of_bits=%"PRIu64" to be <= bits->nbits=%"PRIu64,
            new_num_of_bits, bits->nbits);

    bits->nbits  = new_num_of_bits;
    bits->nwords = roundup_bits2words64 (new_num_of_bits);   

    bits_clear_excess_bits_in_top_word (bits, true);
}

#ifndef __LITTLE_ENDIAN__
void LTEN_bits (BitsP bits)
{
    int64_t num_words = (bits->nbits >> 6) + ((bits->nbits & 0x3f) != 0);
    
    for (uint64_t i=0; i < num_words; i++)
        bits->words[i] = LTEN64 (bits->words[i]);
}
#endif

// xor the entire dst bit array, with a subset of xor_with
void bits_xor_with (BitsP dst, ConstBitsP xor_with, uint64_t xor_with_bit)
{
    ASSERT (xor_with_bit + dst->nbits <= xor_with->nbits, 
            "Expecting xor_with_bit=%"PRId64" + dst->nbits=%"PRId64" <= xor_with->nbits=%"PRId64, 
            xor_with_bit, dst->nbits, xor_with->nbits);

    uint64_t num_of_full_words = dst->nbits / WORD_SIZE;

    // full words
    for (uint64_t word_i=0; word_i < num_of_full_words; word_i++, xor_with_bit += 64)
        dst->words[word_i] ^= _get_word (xor_with, xor_with_bit);

    // last partial word - combine xored part of the word, part of the dst word beyond the requested range
    uint8_t remaining_bits = dst->nbits & 0x3f;
    if (remaining_bits) 
        dst->words[dst->nwords-1] ^= _get_word (xor_with, xor_with_bit) & bitmask64 (remaining_bits);

    DEBUG_VALIDATE(dst);
}

// calculate the number of bits that are different between a bit array (in its entirety)
// and a forward or reverse-complemented region of another bit array
// note: static inline because in the tight loop of aligner
uint32_t bits_hamming_distance (ConstBits𐤐 bits_1, // the entire bit array (restrict-ed!)
                                ConstBits𐤐 bits_2, 
                                uint64_t index_2)  // index first bit in bits_2 
{
    const uint64_t *restrict words_1 = bits_1->words;
    const uint64_t *restrict words_2 = &bits_2->words[index_2 >> 6];
    const uint64_t *after_1 = words_1 + bits_1->nwords;
    uint64_t diff=0;
    int shift_2 = index_2 & 0b111111;
    uint32_t distance=0; // number of non-matching bits

    // if bits_1 goes beyond the end of bits_2 we can't calculate distance
    if (index_2 + bits_1->nbits > bits_2->nbits)
        return bits_1->nbits; // maximum distance = complete misfit

    uint64_t next_w2 = *words_2; 
    for (const uint64_t *w1=words_1, *w2=words_2; w1 < after_1; w1++, w2++) {
        uint64_t this_w2 = next_w2;
        next_w2 = *(w2+1);
        // note: if last bits_1 word is partial, we still calculate it in its entirety and correct later. However, in case we reach 
        // the end of bits_2 and have shift_2, we will access the word beyond the end of the bits_2->words. counting on words being a Buffer to not cause a segfault.
        uint64_t w2_value = _bits_combined_word (this_w2, next_w2, shift_2);
        diff = *w1 ^ w2_value; 
        distance += __builtin_popcountll (diff); // note: expected to be _mm_popcnt_u64 with SSE4.2
    }

    // remove distance due to the unused part of the last bits_1 word
    if (bits_1->nbits & 0b111111)
        distance -= __builtin_popcountll (diff & ~bitmask64 (bits_1->nbits & 0b111111));

    return distance; 
}