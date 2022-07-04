// ------------------------------------------------------------------
//   bits.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
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

//
// Tables of constants
//

#define assert(x) ASSERT ((x), "%s", #x)

#ifndef __GNUC__

// See http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
static uint64_t __inline windows_popcount (uint64_t w)
{
    w = w - ((w >> 1) & (uint64_t)~(uint64_t)0/3);
    w = (w & (uint64_t)~(uint64_t)0/15*3) + ((w >> 2) & (uint64_t)~(uint64_t)0/15*3);
    w = (w + (w >> 4)) & (uint64_t)~(uint64_t)0/255*15;
    w = (uint64_t)(w * ((uint64_t)~(uint64_t)0/255)) >> (sizeof(uint64_t) - 1) * 8;
}

#define POPCOUNT(x) windows_popcountl(x)
#else
#define POPCOUNT(x) (unsigned)__builtin_popcountll(x)
#endif

// word of all 1s
#define WORD_MAX  (~(uint64_t)0)

#define SET_REGION(arr,start,len)    _set_region((arr),(start),(len),FILL_REGION)
#define CLEAR_REGION(arr,start,len)  _set_region((arr),(start),(len),ZERO_REGION)
#define TOGGLE_REGION(arr,start,len) _set_region((arr),(start),(len),SWAP_REGION)

//
// Common internal functions
//

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
static inline void validate_bits (ConstBitsP arr, rom file, int lineno)
{
    // Verify that its allocated
    ASSERT (arr->type != BITARR_UNALLOCATED, "[%s:%i] Bits is not allocated", file, lineno);

    // Check num of words is correct
    uint64_t num_words = roundup_bits2words64(arr->nbits);
    ASSERT (num_words == arr->nwords, "[%s:%i] num of words wrong [bits: %i, expected nwords: %i, actual nwords: %i]", 
            file, lineno, (int)arr->nbits, (int)num_words, (int)arr->nwords);

    // Check top word is masked (only if not overlayed - the unused bits of top word don't belong to this bit array and might be used eg by another bit array in genome.ref/genome.is_set)
    if (arr->type == BITARR_REGULAR) {
      uint64_t tw = arr->nwords == 0 ? 0 : arr->nwords - 1;
      uint64_t top_bits = bits_in_top_word(arr->nbits);

      ASSERT (arr->words[tw] <= bitmask64(top_bits), "[%s:%i] Expected %i bits in top word[%i] (the rest should be 0) but word=%s\n", 
              file, lineno, (int)top_bits, (int)tw, _print_word(arr->words[tw]).s);
    }
}

#else
    #define DEBUG_VALIDATE(a)
#endif

// Reverse a word
static inline uint64_t _reverse_word (uint64_t word)
{
    static const uint64_t reverse_table[256] = {
        0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
        0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
        0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
        0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
        0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
        0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
        0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
        0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
        0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
        0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
        0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
        0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
        0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
        0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
        0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
        0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF,
    };

    uint64_t reverse = (reverse_table[(word)       & 0xff] << 56) |
                       (reverse_table[(word >>  8) & 0xff] << 48) |
                       (reverse_table[(word >> 16) & 0xff] << 40) |
                       (reverse_table[(word >> 24) & 0xff] << 32) |
                       (reverse_table[(word >> 32) & 0xff] << 24) |
                       (reverse_table[(word >> 40) & 0xff] << 16) |
                       (reverse_table[(word >> 48) & 0xff] << 8)  |
                       (reverse_table[(word >> 56) & 0xff]);

    return reverse;
}

// Set 64 bits from a particular start position
// Doesn't extend bit array
static inline void _set_word (BitsP bits, uint64_t start, uint64_t word)
{
    uint64_t word_index = bitset64_wrd(start);
    word_offset_t word_offset = bitset64_idx(start);

    if (word_offset == 0)
        bits->words[word_index] = word;
    
    else {
        bits->words[word_index] = (word << word_offset) |
                                  (bits->words[word_index] & bitmask64(word_offset));

        if (word_index+1 < bits->nwords) { // if last part of the word goes beyond nwords, we drop it

            bits->words[word_index+1] = (word >> (WORD_SIZE - word_offset)) |
                                        (bits->words[word_index+1] & (WORD_MAX << word_offset));

            // added by divon: Mask top word only if its the last word --divon
            if (word_index+2 == bits->nwords)
                bits_clear_excess_bits_in_top_word (bits);
        }
    }
}

static inline void _set_byte (BitsP bits, uint64_t start, uint8_t byte)
{
    uint64_t w = _get_word (bits, start);
    _set_word (bits, start, (w & ~(uint64_t)0xff) | byte);
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

//
// Fill a region (internal use only)
//

// FillAction is fill with 0 or 1 or toggle
typedef enum {ZERO_REGION, FILL_REGION, SWAP_REGION} FillAction;

static inline void _set_region (BitsP bits, uint64_t start, uint64_t length, FillAction action)
{
    if (length == 0) return;

    uint64_t first_word = bitset64_wrd(start);
    uint64_t last_word = bitset64_wrd(start+length-1);
    word_offset_t foffset = bitset64_idx(start);
    word_offset_t loffset = bitset64_idx(start+length-1);

    if (first_word == last_word) {
        uint64_t mask = bitmask64(length) << foffset;

        switch (action) {
            case ZERO_REGION: bits->words[first_word] &= ~mask; break;
            case FILL_REGION: bits->words[first_word] |=  mask; break;
            case SWAP_REGION: bits->words[first_word] ^=  mask; break;
        }
    }
    
    else {
        // Set first word
        switch (action) {
            case ZERO_REGION: bits->words[first_word] &=  bitmask64(foffset); break;
            case FILL_REGION: bits->words[first_word] |= ~bitmask64(foffset); break;
            case SWAP_REGION: bits->words[first_word] ^= ~bitmask64(foffset); break;
        }

        uint64_t i;

        // Set whole words
        switch (action) {
            case ZERO_REGION:
                for (i = first_word + 1; i < last_word; i++)
                    bits->words[i] = (uint64_t)0;
                break;
            case FILL_REGION:
                for (i = first_word + 1; i < last_word; i++)
                    bits->words[i] = WORD_MAX;
                break;
            case SWAP_REGION:
                for (i = first_word + 1; i < last_word; i++)
                    bits->words[i] ^= WORD_MAX;
                break;
        }

        // Set last word
        switch (action) {
            case ZERO_REGION: bits->words[last_word] &= ~bitmask64(loffset+1); break;
            case FILL_REGION: bits->words[last_word] |=  bitmask64(loffset+1); break;
            case SWAP_REGION: bits->words[last_word] ^=  bitmask64(loffset+1); break;
        }
    }
}

//
// Constructor
//

// allocates a STANDALONE bit array
Bits bits_alloc_do (uint64_t nbits, bool clear, FUNCLINE)
{
    Bits bits = { .type   = BITARR_STANDALONE, // standalone Bits that is not part of a Buffer
                  .nbits  = nbits,
                  .nwords = roundup_bits2words64(nbits),
                  .words  = buf_low_level_malloc (roundup_bits2bytes64(nbits), clear, func, code_line) };

    // zero the bits in the top word that are beyond nbits (if not already cleared)
    if (!clear) bits_clear_excess_bits_in_top_word (&bits);
    
    return bits;
}

// reallocates a Buffer-embedded bit array
void bits_realloc_do (BitsP bits, uint64_t nbits, 
                      uint64_t low_level_nbits, // if not 0, we over-allocated in anticipation of further reallocs (improves performance)
                      bool clear, FUNCLINE)
{
    ASSERT0 (bits->type == BITARR_UNALLOCATED || bits->type == BITARR_REGULAR || bits->type == BITARR_STANDALONE, 
             "bit array needs to be BITARR_UNALLOCATED or BITARR_REGULAR or BITARR_STANDALONE for realloc");

    uint64_t old_nbits = bits->nbits;

    if (!low_level_nbits) low_level_nbits = nbits;
    uint64_t nwords = roundup_bits2words64 (low_level_nbits);

    if (bits->type == BITARR_STANDALONE)
        bits->words = buf_low_level_realloc (bits->words, nwords * sizeof(uint64_t), "", func, code_line);

    else {
        BufferP buf = buf_get_buffer_from_bits (bits);
        buf_alloc_do (buf->vb, buf, nwords * sizeof(uint64_t), 1, func, code_line, NULL);
    }

    bits->nbits  = nbits;
    bits->nwords = roundup_bits2words64 (nbits); // possibly less than low_level_nbits

    if (clear && nbits >old_nbits) 
        bits_clear_region (bits, old_nbits, nbits - old_nbits);

    bits_clear_excess_bits_in_top_word (bits); // zero the bits in the top word that are beyond nbits
}

void bits_free (BitsP bits)
{
    FREE (bits->words);
    memset (bits, 0, sizeof(Bits));
}

uint64_t bits_length (ConstBitsP bit_arr)
{
    return bit_arr->nbits;
}

// If bits length < num_bits, resizes to num_bits
char bits_ensure_size (BitsP bits, uint64_t ensure_num_of_bits)
{
    ASSERT (bits->nbits >= ensure_num_of_bits, "bits of out range: bits->nbits=%"PRIu64" ensure_num_of_bits=%"PRIu64,
            bits->nbits, ensure_num_of_bits);

    return 1;
}

static inline void bits_ensure_size_critical (BitsP bits, uint64_t nbits)
{
    ASSERT (bits->nbits >= nbits, "bits of out range: bits->nbits=%"PRIu64" nbits=%"PRIu64,
            bits->nbits, nbits);
}

//
// Get, set, clear, assign and toggle individual bits
//

// Get the value of a bit (returns 0 or 1)
char bits_get_bit (ConstBitsP bits, uint64_t b)
{
    ASSERT (b < bits->nbits, "Expecting b(%"PRId64") < bits->nbits(%"PRId64")", b, bits->nbits);
    return bits_get(bits, b);
}

// set a bit (to 1) at position b
void bits_set_bit (BitsP bits, uint64_t b)
{
    assert(b < bits->nbits);
    bits_set(bits,b);
    DEBUG_VALIDATE(bits);
}

// clear a bit (to 0) at position b
void bits_clear_bit (BitsP bits, uint64_t b)
{
    assert(b < bits->nbits);
    bits_clear(bits, b);
    DEBUG_VALIDATE(bits);
}

// If bit is 0 -> 1, if bit is 1 -> 0.  AKA 'flip'
void bits_toggle_bit (BitsP bits, uint64_t b)
{
    assert(b < bits->nbits);
    bits_toggle(bits, b);
    DEBUG_VALIDATE(bits);
}

// If char c != 0, set bit; otherwise clear bit
void bits_assign_bit (BitsP bits, uint64_t b, char c)
{
    assert(b < bits->nbits);
    bits_assign(bits, b, c ? 1 : 0);
    DEBUG_VALIDATE(bits);
}

//
// Get, set, clear and toggle several bits at once
//

// Get the offsets of the set bits (for offsets start<=offset<end)
// Returns the number of bits set
// It is assumed that dst is at least of length (end-start)
uint64_t bits_get_bits (ConstBitsP bits, uint64_t start, uint64_t end, uint64_t *dst)
{
    uint64_t i, n = 0;
    assert(end <= bits->nbits);
    for (i = start; i < end; i++) 
        if (bits_get(bits, i)) dst[n++] = i;

    return n;
}

// Set multiple bits at once.
// e.g. set bits 1, 20 & 31: bits_set_bits(bits, 3, 1,20,31);
void bits_set_bits (BitsP bits, size_t n, ...)
{
    va_list argptr;
    va_start(argptr, n);

    for (size_t i = 0; i < n; i++) {
        unsigned int bit_index = va_arg(argptr, unsigned int);
        bits_set_bit(bits, bit_index);
    }

    va_end(argptr);
    DEBUG_VALIDATE(bits);
}

// Clear multiple bits at once.
// e.g. clear bits 1, 20 & 31: bits_clear_bits(bits, 3, 1,20,31);
void bits_clear_bits (BitsP bits, size_t n, ...)
{
    va_list argptr;
    va_start(argptr, n);

    for (size_t i = 0; i < n; i++) {
        unsigned int bit_index = va_arg(argptr, unsigned int);
        bits_clear_bit(bits, bit_index);
    }

    va_end(argptr);
    DEBUG_VALIDATE(bits);
}

//
// Set, clear all bits in a region
//

// Set all the bits in a region
void bits_set_region (BitsP bits, uint64_t start, uint64_t len)
{
    if (!len) return; // nothing to do 

    ASSERT (start + len - 1 <= bits->nbits, "Expecting: start(%"PRId64") + len(%"PRId64") - 1 <= bits->nbits(%"PRId64")",
            start, len, bits->nbits); // divon fixed bug

    SET_REGION (bits, start, len);
    DEBUG_VALIDATE (bits);
}


// Clear all the bits in a region
void bits_clear_region_do (BitsP bits, uint64_t start, uint64_t len, FUNCLINE)
{
    if (!len) return; // nothing to do 

    ASSERT (start + len - 1 <= bits->nbits, "called from %s:%u: Expecting: start(%"PRId64") + len(%"PRId64") - 1 <= bits->nbits(%"PRId64")",
            func, code_line, start, len, bits->nbits); // divon fixed bug

    CLEAR_REGION (bits, start, len);
    DEBUG_VALIDATE (bits);
}

//
// Set, clear all bits at once
//

// set all elements of data to one
void bits_set_all (BitsP bits)
{
    if (!bits->nwords) return; // nothing to do
    
    uint64_t num_of_bytes = bits->nwords * sizeof(uint64_t);
    memset(bits->words, 0xFF, num_of_bytes);
    bits_clear_excess_bits_in_top_word(bits);
    DEBUG_VALIDATE(bits);
}

// set all elements of data to zero
void bits_clear_all (BitsP bits)
{
    if (!bits->nwords) return; // nothing to do

    memset(bits->words, 0, bits->nwords * sizeof(uint64_t));
    DEBUG_VALIDATE(bits);
}

//
// Get a word at a time
//

uint64_t bits_get_wordn (ConstBitsP bits, uint64_t start, int n /* up to 64 */)
{
  ASSERT (start + n <= bits->nbits, "expecting start=%"PRIu64" + n=%d <= bits->nbits=%"PRIu64, 
          start, n, bits->nbits);
           
  return (uint64_t)(_get_word(bits, start) & bitmask64(n));
}

//
// Set a word at a time
//
// Doesn't extend bit array. However it is safe to TRY to set bits beyond the
// end of the array, as long as: `start` is < `bits_length(arr)`
//

void bits_set_word64 (BitsP bits, uint64_t start, uint64_t word)
{
    ASSERT (start < bits->nbits, "expecting start(%"PRIu64") < bits->nbits(%"PRIu64")", start, bits->nbits);
    _set_word(bits, start, (uint64_t)word);
}

void bits_set_word32 (BitsP bits, uint64_t start, uint32_t word)
{
    ASSERT (start < bits->nbits, "expecting start(%"PRIu64") < bits->nbits(%"PRIu64")", start, bits->nbits);
    uint64_t w = _get_word(bits, start);
    _set_word(bits, start, bitmask_merge(w, word, 0xffffffff00000000UL));
}

void bits_set_word16 (BitsP bits, uint64_t start, uint16_t word)
{
    ASSERT (start < bits->nbits, "expecting start(%"PRIu64") < bits->nbits(%"PRIu64")", start, bits->nbits);
    uint64_t w = _get_word(bits, start);
    _set_word(bits, start, bitmask_merge(w, word, 0xffffffffffff0000UL));
}

void bits_set_word8 (BitsP bits, uint64_t start, uint8_t byte)
{
    ASSERT (start < bits->nbits, "expecting start(%"PRIu64") < bits->nbits(%"PRIu64")", start, bits->nbits);
    _set_byte(bits, start, byte);
}

void bits_set_wordn(BitsP bits, uint64_t start, uint64_t word, int n)
{
    ASSERT (start < bits->nbits, "expecting start(%"PRIu64") < bits->nbits(%"PRIu64")", start, bits->nbits);
    assert(n <= 64);
    uint64_t w = _get_word(bits, start), m = bitmask64(n);
    _set_word(bits, start, bitmask_merge(word,w,m));
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
        if (bits->words[i] != 0xffffffffffffffffULL) return false;

    // last word might be partial
    word_offset_t bits_active = bits_in_top_word (bits->nbits);    
    word_offset_t num_of_bits_set = POPCOUNT (bits->words[bits->nwords-1] & bitmask64 (bits_active));

    return num_of_bits_set == bits_active;
}

// true if all bits in bit array are clear (divon)
bool bits_is_fully_clear (ConstBitsP bits)
{
    ASSERT_excess_bits_are_0 (bits);

    if (!bits->nbits) return true; // trivially true

    // all words but last one
    for (uint64_t i=0; i < bits->nwords; i++)
        if (bits->words[i]) return false;

    return true;
}

// Get the number of bits set (hamming weight)
uint64_t bits_num_set_bits (ConstBitsP bits)
{
    if (!bits->nbits) return 0;

    uint64_t num_of_bits_set = 0;

    // all words but last one
    for (uint64_t i=0; i < bits->nwords-1; i++)
        if (bits->words[i])
            num_of_bits_set += POPCOUNT (bits->words[i]);

    // last word might be partial
    word_offset_t bits_active = bits_in_top_word (bits->nbits);    
    num_of_bits_set += POPCOUNT (bits->words[bits->nwords-1] & bitmask64 (bits_active));
    
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
        num_of_bits_set += POPCOUNT (bits->words[first_word] & mask);
    }
    else {
        word_offset_t loffset  = bitset64_idx(start+length-1);
    
        // first word
        num_of_bits_set += POPCOUNT (bits->words[first_word] & ~bitmask64(foffset));

        // whole words
        for (uint64_t i = first_word + 1; i < last_word; i++)
            num_of_bits_set += POPCOUNT (bits->words[i]);

        // last word
        num_of_bits_set += POPCOUNT (bits->words[last_word] & bitmask64(loffset+1));
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
    while(1) { \
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
    assert(offset <= bits->nbits); \
    if (bits->nbits == 0 || offset == 0) { return false; } \
    \
    /* Find prev word that is greater than zero */ \
    uint64_t i = bitset64_wrd(offset-1); \
    uint64_t w = GET(bits->words[i]) & bitmask64(bitset64_idx(offset-1)+1); \
    \
    if (w > 0) { *result = (i+1) * WORD_SIZE - leading_zeros(w) - 1; return true; } \
    \
    /* i is unsigned so have to use break when i == 0 */ \
    for (--i; i != BIT_INDEX_MAX; i--) { \
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
// Returns 1 if a bit is set, otherwise 0
// Index of first set bit is stored in the integer pointed to by result
// If no bits are set, value at `result` is not changed
bool bits_find_first_set_bit (ConstBitsP bits, uint64_t *result)
{
    return bits_find_next_set_bit(bits, 0, result);
}

// same same
bool bits_find_first_clear_bit (ConstBitsP bits, uint64_t *result)
{
    return bits_find_next_clear_bit(bits, 0, result);
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

// Construct a Bits from a substring with given on and off characters.
void bits_from_substr (BitsP bits, uint64_t offset, rom str, size_t len, rom on, rom off,char left_to_right)
{
    bits_ensure_size (bits, offset + len);
    bits_clear_region (bits, offset, len);

    // Bits region is now all 0s -- just set the 1s
    size_t i;
    uint64_t j;

    for (i = 0; i < len; i++)
        if (strchr(on, str[i]) != NULL) {
            j = offset + (left_to_right ? i : len - i - 1);
            bits_set(bits, j);
        }
        else { assert(strchr(off, str[i]) != NULL); }

    DEBUG_VALIDATE(bits);
}

// From string method
void bits_from_str (BitsP bits, rom str)
{
    bits_from_substr(bits, 0, str, strlen(str), "1", "0", 1);
}

// Takes a char array to write to.  `str` must be bits->nbits+1 in length
// Terminates string with '\0'
char *bits_to_str (ConstBitsP bits, char *str)
{
    for (uint64_t i = 0; i < bits->nbits; i++)
        str[i] = bits_get(bits, i) ? '1' : '0';

    str[bits->nbits] = '\0';

    return str;
}

char *bits_to_str_rev (ConstBitsP bits, char *str)
{
    for (uint64_t i = 0; i < bits->nbits; i++)
        str[i] = bits_get(bits, bits->nbits-i-1) ? '1' : '0';

    str[bits->nbits] = '\0';

    return str;
}


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
    length = MIN_(length, bits->nbits - start);

    if (msg) fprintf (file, "%s: ", msg);
    for (uint64_t i=0; i < length; i++) 
        fputc (bits_get(bits, start + i) ? '1' : '0', file);

    fputc ('\n', file);
}

// prints a 2bit bases array - A,C,G,T bases
void bits_print_substr_bases (rom msg, ConstBitsP bits, uint64_t start, uint64_t length/*in bases*/, FILE *file)
{
    start *= 2; // now it is in bits
    length = MIN_(length * 2, bits->nbits - start);

    if (msg) fprintf (file, "%s: ", msg);
    for (uint64_t i=0; i < length; i += 2) { 
        char b = bits_get2 (bits, start + i);
        fputc (b==0 ? 'A' : b==1 ? 'C' : b==2 ? 'G' : 'T' , file);
    }

    fputc ('\n', file);
}

//
// Clone and copy
//

// destination and source may be the same bits
// and src/dst regions may overlap
static void _array_copy (BitsP dst, uint64_t dstindx, ConstBitsP src, uint64_t srcindx, uint64_t length)
{
    // Num of full words to copy
    uint64_t num_of_full_words = length / WORD_SIZE;
    uint64_t i;

    word_offset_t bits_in_last_word = length % WORD_SIZE; //bits_in_top_word(length); // divon: fixed this bug in the original library

    if (dst == src && srcindx > dstindx)
    {
        // Work left to right
        for (i=0; i < num_of_full_words; i++)
        {
            uint64_t word = _get_word(src, srcindx+i*WORD_SIZE);
            _set_word(dst, dstindx+i*WORD_SIZE, word);
        }

      if (bits_in_last_word > 0)
      {
          uint64_t src_word = _get_word(src, srcindx+i*WORD_SIZE);
          uint64_t dst_word = _get_word(dst, dstindx+i*WORD_SIZE);

          uint64_t mask = bitmask64(bits_in_last_word);
          uint64_t word = bitmask_merge(src_word, dst_word, mask);

          _set_word(dst, dstindx+num_of_full_words*WORD_SIZE, word);
      }
    }
    else
    {
        // Work right to left
        for (i = 0; i < num_of_full_words; i++)
        {
            uint64_t word = _get_word(src, srcindx+length-(i+1)*WORD_SIZE);
            _set_word(dst, dstindx+length-(i+1)*WORD_SIZE, word);
        }

        if (bits_in_last_word > 0)
        {
            uint64_t src_word = _get_word(src, srcindx);
            uint64_t dst_word = _get_word(dst, dstindx);

            uint64_t mask = bitmask64(bits_in_last_word);
            uint64_t word = bitmask_merge(src_word, dst_word, mask);
            _set_word(dst, dstindx, word);
        }
    }
}

// destination and source may be the same bits
// and src/dst regions may overlap
void bits_copy_do (BitsP dst, uint64_t dstindx, ConstBitsP src, uint64_t srcindx, uint64_t length, FUNCLINE)
{
    DEBUG_VALIDATE(src);

    ASSERT (dstindx + length <= dst->nbits, "called from %s:%u dstindx(%"PRIu64") + length(%"PRIu64") > dst->nbits(%"PRIu64")", func, code_line, dstindx, length, dst->nbits);
    ASSERT (srcindx + length <= src->nbits, "called from %s:%u srcindx(%"PRIu64") + length(%"PRIu64") > src->nbits(%"PRIu64")", func, code_line, srcindx, length, src->nbits);

    _array_copy(dst, dstindx, src, srcindx, length);

    // note: only zero top word if copy reaches it
    if (dstindx + length == dst->nbits)
        bits_clear_excess_bits_in_top_word(dst);

    DEBUG_VALIDATE(dst);
}

// concatenate src at the end of dst (divon)
void bits_concat_do (BitsP base, ConstBitsP add, unsigned additional_concats_expected, FUNCLINE) // only used for non-buffer: low-level realloc more than needed, expecting additional concats - improves perforamce
{
    uint64_t index = base->nbits;

    // IMPORTANT: if calling with a Bits that is embedded in a buffer, buf_add_to_buffer_list must be 
    // called first (eg in seg_initialize), as bits_realloc_do performs an anonymous realloc 
    bits_realloc_do (base, 
                     base->nbits + add->nbits, 
                     base->nbits + (add->nbits * (1+additional_concats_expected)), // low level allocation of anticipated additional similar-size concats - improves performance
                     false, func, code_line);
                       
//xxx fix this to support additional_concats_expected in buffer. 
/*    BufferP buf = buf_get_buffer_from_bits (base);

    uint64_t more_words =  roundup_bits2words64 (add->nbits * (1+additional_concats_expected));

    buf_alloc (buf->vb, buf, more_words, 0, uint64_t, 1, NULL);
    buf_extend_bits (buf, add->nbits);
    bits_clear_excess_bits_in_top_word (base);
*/
    bits_copy (base, index, add, 0, add->nbits);
}

void bits_overlay (BitsP overlaid_bits, BitsP regular_bits, uint64_t start, uint64_t nbits)
{
    ASSERT (start % 64 == 0, "start=%"PRIu64" must be a multiple of 64", start);
    ASSERT (start + nbits <= regular_bits->nbits, "start(%"PRIu64") + nbits(%"PRIu64") <= regular_bits->nbits(%"PRIu64")",
            start, nbits, regular_bits->nbits);

    uint64_t word_i = start / 64;
    *overlaid_bits = (Bits){ .type   = BITARR_OVERLAY,
                                   .words  = &regular_bits->words[word_i],
                                   .nwords = roundup_bits2words64 (nbits),
                                   .nbits  = nbits };
} 

// convert a bitsay to a byte array of values 0-1
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

// convert a Bits containing a series of 2-bits, to a byte array of values 0-3
void bits_base_to_byte (uint8_t *dst, ConstBitsP src_bits, uint64_t base_i, uint32_t num_bases)
{
    ASSERT (2*(base_i + num_bases) <= src_bits->nbits, "Expecting 2*(base_is=%"PRIu64" + num_bases=%u) <= nbits=%"PRIu64,
            base_i, num_bases, src_bits->nbits);

    BASE_ITER_INIT (src_bits, base_i, num_bases, true);

    for (uint32_t i=0; i < num_bases; i++) 
        *dst++ = BASE_NEXT_FWD;
}

//
// Logic operators
//

// Destination can be the same as one or both of the sources
void bits_and (BitsP dst, ConstBitsP src1, ConstBitsP src2)
{
    // Ensure dst array is big enough
    uint64_t max_bits = MAX_(src1->nbits, src2->nbits);
    bits_ensure_size_critical(dst, max_bits);

    uint64_t min_words = MIN_(src1->nwords, src2->nwords);

    for (uint64_t i = 0; i < min_words; i++)
        dst->words[i] = src1->words[i] & src2->words[i];

    // Set remaining bits to zero
    for (uint64_t i = min_words; i < dst->nwords; i++)
        dst->words[i] = (uint64_t)0;

  DEBUG_VALIDATE(dst);
}

// Destination can be the same as one or both of the sources
static void _logical_or_xor(BitsP dst, ConstBitsP src1, ConstBitsP src2,char use_xor)
{
    // Ensure dst array is big enough
    bits_ensure_size_critical(dst, MAX_(src1->nbits, src2->nbits));
    DEBUG_VALIDATE(src1);
    DEBUG_VALIDATE(src2);
    if (dst != src1 && dst != src2) DEBUG_VALIDATE(dst);

    uint64_t min_words = MIN_(src1->nwords, src2->nwords);
    uint64_t max_words = MAX_(src1->nwords, src2->nwords);

    if (use_xor) {
        for (uint64_t i = 0; i < min_words; i++)
        dst->words[i] = src1->words[i] ^ src2->words[i];
    }
    else {
        for (uint64_t i = 0; i < min_words; i++)
        dst->words[i] = src1->words[i] | src2->words[i];
    }

    // Copy remaining bits from longer src array
    if (min_words != max_words) {
        ConstBitsP longer = src1->nwords > src2->nwords ? src1 : src2;

        for (uint64_t i = min_words; i < max_words; i++)
            dst->words[i] = longer->words[i];
    }

    // Set remaining bits to zero
    size_t size = (dst->nwords - max_words) * sizeof(uint64_t);
    memset(dst->words + max_words, 0, size);

    DEBUG_VALIDATE(dst);
}

void bits_or (BitsP dst, ConstBitsP src1, ConstBitsP src2)
{
    _logical_or_xor(dst, src1, src2, 0);
}

// Destination can be the same as one or both of the sources
void bits_xor (BitsP dst, ConstBitsP src1, ConstBitsP src2)
{
    _logical_or_xor(dst, src1, src2, 1);
}

// xor the entire dst bit array, with a subset of xor_with
void bits_xor_with (BitsP dst, uint64_t dst_bit, ConstBitsP xor_with, uint64_t xor_with_bit, uint64_t num_bits)
{
    ASSERT (xor_with_bit + num_bits <= xor_with->nbits, 
            "Expecting xor_with_bit=%"PRId64" + num_bits=%"PRId64" <= xor_with->nbits=%"PRId64, xor_with_bit, num_bits, xor_with->nbits);

    ASSERT (dst_bit + num_bits <= dst->nbits, 
            "Expecting dst_bit=%"PRId64" + num_bits=%"PRId64" <= dst->nbits=%"PRId64, dst_bit, num_bits, dst->nbits);

    uint64_t num_of_full_words = num_bits / WORD_SIZE;

    // full words
    for (uint64_t i=0; i < num_of_full_words; i++, xor_with_bit += 64, dst_bit += 64) {
        uint64_t xored = _get_word (dst, dst_bit) ^ _get_word (xor_with, xor_with_bit);
        _set_word (dst, dst_bit, xored);
    }

    // last partial word - combine xored part of the word, part of the dst word beyond the requested range
    uint8_t remaining_bits = num_bits - num_of_full_words * 64;
    if (remaining_bits) {
        uint64_t dst_word = _get_word (dst, dst_bit);
        uint64_t mask     = bitmask64 (remaining_bits);
        uint64_t xored    = (dst_word ^ _get_word (xor_with, xor_with_bit)) & mask;
        uint64_t result   = xored | (dst_word & ~mask);
        _set_word (dst, dst_bit, result);
    }

    DEBUG_VALIDATE(dst);
}

// If dst is longer than src, top bits are set to 1
void bits_not (BitsP dst, ConstBitsP src)
{
    bits_ensure_size_critical(dst, src->nbits);

    for (uint64_t i = 0; i < src->nwords; i++)
        dst->words[i] = ~(src->words[i]);

    // Set remaining words to 1s
    for (uint64_t i = src->nwords; i < dst->nwords; i++)
        dst->words[i] = WORD_MAX;

    bits_clear_excess_bits_in_top_word(dst);

    DEBUG_VALIDATE(dst);
}

//
// Comparisons
//

// Compare two bit arrays by value stored, with index 0 being the Least
// Significant Bit (LSB). Arrays do not have to be the same length.
// Example: ..0101 (5) > ...0011 (3) [index 0 is LSB at right hand side]
// Sorts on length if all zeros: (0,0) < (0,0,0)
// returns:
//  >0 iff bits1 > bits2
//   0 iff bits1 == bits2
//  <0 iff bits1 < bits2
int bits_cmp (ConstBitsP bits1, ConstBitsP bits2)
{
    uint64_t word1, word2;
    uint64_t min_words = bits1->nwords;

    // i is unsigned so break when i == 0
    if (bits1->nwords > bits2->nwords) {
        min_words = bits2->nwords;
        for (uint64_t i = bits1->nwords-1; ; i--) {
            if (bits1->words[i]) return 1;
            if (i == bits2->nwords) break;
        }
    }
    else if (bits1->nwords < bits2->nwords) {
        for (uint64_t i = bits2->nwords-1; ; i--) {
            if (bits2->words[i]) return 1;
            if (i == bits1->nwords) break;
        }
    }

    if (min_words == 0) return 0;

    for (uint64_t i = min_words-1; ; i--)
    {
        word1 = bits1->words[i];
        word2 = bits2->words[i];
        if (word1 != word2) return (word1 > word2 ? 1 : -1);
        if (i == 0) break;
    }

    if (bits1->nbits == bits2->nbits) return 0;
    return bits1->nbits > bits2->nbits ? 1 : -1;
}

// Compare two bit arrays by value stored, with index 0 being the Most
// Significant Bit (MSB). Arrays do not have to be the same length.
// Example: 10.. > 01.. [index 0 is MSB at left hand side]
// Sorts on length if all zeros: (0,0) < (0,0,0)
// returns:
//  >0 iff bits1 > bits2
//   0 iff bits1 == bits2
//  <0 iff bits1 < bits2
int bits_cmp_big_endian (ConstBitsP bits1, ConstBitsP bits2)
{
    uint64_t min_words = MAX_(bits1->nwords, bits2->nwords);

    uint64_t i;
    uint64_t word1, word2;

    for (i = 0; i < min_words; i++) {
        word1 = _reverse_word(bits1->words[i]);
        word2 = _reverse_word(bits2->words[i]);
        if (word1 != word2) return (word1 > word2 ? 1 : -1);
    }

    // Check remaining words. Only one of these loops will execute
    for (; i < bits1->nwords; i++)
        if (bits1->words[i]) return 1;
    for (; i < bits2->nwords; i++)
        if (bits2->words[i]) return -1;

    if (bits1->nbits == bits2->nbits) return 0;
    return bits1->nbits > bits2->nbits ? 1 : -1;
}

// compare bits with (bits2 << pos)
// bits_cmp(bits1, bits2<<pos)
// returns:
//  >0 iff bits1 > bits2
//   0 iff bits1 == bits2
//  <0 iff bits1 < bits2
int bits_cmp_words (ConstBitsP arr1,uint64_t pos, ConstBitsP arr2)
{
    if (arr1->nbits == 0 && arr2->nbits == 0) return 0;

    uint64_t top_bit1 = 0, top_bit2 = 0;

    bool arr1_zero = !bits_find_last_set_bit(arr1, &top_bit1);
    bool arr2_zero = !bits_find_last_set_bit(arr2, &top_bit2);

    if (arr1_zero && arr2_zero) return 0;
    if (arr1_zero) return -1;
    if (arr2_zero) return 1;

    uint64_t top_bit2_offset = top_bit2 + pos;

    if (top_bit1 != top_bit2_offset) 
        return top_bit1 > top_bit2_offset ? 1 : -1;

    uint64_t i, word1, word2;
    for (i = top_bit2 / WORD_SIZE; i > 0; i--) {
        word1 = _get_word(arr1, pos + i * WORD_SIZE);
        word2 = arr2->words[i];

        if (word1 > word2) return 1;
        if (word1 < word2) return -1;
    }

    word1 = _get_word(arr1, pos);
    word2 = arr2->words[0];

    if (word1 > word2) return 1;
    if (word1 < word2) return -1;

    // return 1 if arr1[0..pos] != 0, 0 otherwise

    // Whole words
    uint64_t num_words = pos / WORD_SIZE;

    for (i = 0; i < num_words; i++)
        if (arr1->words[i] > 0) return 1;

    word_offset_t bits_remaining = pos - num_words * WORD_SIZE;

    if (arr1->words[num_words] & bitmask64(bits_remaining)) return 1;

    return 0;
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

// for each 2 bits in the src array, the dst array will contain those 2 bits in the reverse
// position, as well as transform them 00->11 11->00 01->10 10->01
// src_start_base and max_num_bases must be multiplies of 32 and src != dst: one can call this function piecemiel - eg divide it to threads
// (2) src_start_base and max_num_bases are both 0, but src and dst MAY be the same
void bits_reverse_complement_aligned (BitsP dst, ConstBitsP src, uint64_t src_start_base, uint64_t max_num_bases) 
{
    ASSERT (src != dst && src_start_base % 32 == 0 && max_num_bases % 32 == 0, "Expecting src_start_base=%"PRIu64" and max_num_bases=%"PRIu64" to be multiples of 32 and src=%p != dst=%p",
            src_start_base, max_num_bases, src, dst);

    if (!max_num_bases) max_num_bases = src->nbits / 2; // entire Bits

    static const uint64_t rev_comp_table[256] = { // 00=A 01=C 10=G 11=T
        0b11111111, 0b10111111, 0b01111111, 0b00111111, 0b11101111, 0b10101111, 0b01101111, 0b00101111, 0b11011111, 0b10011111, 0b01011111, 0b00011111, 0b11001111, 0b10001111, 0b01001111, 0b00001111,
        0b11111011, 0b10111011, 0b01111011, 0b00111011, 0b11101011, 0b10101011, 0b01101011, 0b00101011, 0b11011011, 0b10011011, 0b01011011, 0b00011011, 0b11001011, 0b10001011, 0b01001011, 0b00001011,
        0b11110111, 0b10110111, 0b01110111, 0b00110111, 0b11100111, 0b10100111, 0b01100111, 0b00100111, 0b11010111, 0b10010111, 0b01010111, 0b00010111, 0b11000111, 0b10000111, 0b01000111, 0b00000111,
        0b11110011, 0b10110011, 0b01110011, 0b00110011, 0b11100011, 0b10100011, 0b01100011, 0b00100011, 0b11010011, 0b10010011, 0b01010011, 0b00010011, 0b11000011, 0b10000011, 0b01000011, 0b00000011,
        0b11111110, 0b10111110, 0b01111110, 0b00111110, 0b11101110, 0b10101110, 0b01101110, 0b00101110, 0b11011110, 0b10011110, 0b01011110, 0b00011110, 0b11001110, 0b10001110, 0b01001110, 0b00001110,
        0b11111010, 0b10111010, 0b01111010, 0b00111010, 0b11101010, 0b10101010, 0b01101010, 0b00101010, 0b11011010, 0b10011010, 0b01011010, 0b00011010, 0b11001010, 0b10001010, 0b01001010, 0b00001010,
        0b11110110, 0b10110110, 0b01110110, 0b00110110, 0b11100110, 0b10100110, 0b01100110, 0b00100110, 0b11010110, 0b10010110, 0b01010110, 0b00010110, 0b11000110, 0b10000110, 0b01000110, 0b00000110,
        0b11110010, 0b10110010, 0b01110010, 0b00110010, 0b11100010, 0b10100010, 0b01100010, 0b00100010, 0b11010010, 0b10010010, 0b01010010, 0b00010010, 0b11000010, 0b10000010, 0b01000010, 0b00000010,
        0b11111101, 0b10111101, 0b01111101, 0b00111101, 0b11101101, 0b10101101, 0b01101101, 0b00101101, 0b11011101, 0b10011101, 0b01011101, 0b00011101, 0b11001101, 0b10001101, 0b01001101, 0b00001101,
        0b11111001, 0b10111001, 0b01111001, 0b00111001, 0b11101001, 0b10101001, 0b01101001, 0b00101001, 0b11011001, 0b10011001, 0b01011001, 0b00011001, 0b11001001, 0b10001001, 0b01001001, 0b00001001,
        0b11110101, 0b10110101, 0b01110101, 0b00110101, 0b11100101, 0b10100101, 0b01100101, 0b00100101, 0b11010101, 0b10010101, 0b01010101, 0b00010101, 0b11000101, 0b10000101, 0b01000101, 0b00000101,
        0b11110001, 0b10110001, 0b01110001, 0b00110001, 0b11100001, 0b10100001, 0b01100001, 0b00100001, 0b11010001, 0b10010001, 0b01010001, 0b00010001, 0b11000001, 0b10000001, 0b01000001, 0b00000001,
        0b11111100, 0b10111100, 0b01111100, 0b00111100, 0b11101100, 0b10101100, 0b01101100, 0b00101100, 0b11011100, 0b10011100, 0b01011100, 0b00011100, 0b11001100, 0b10001100, 0b01001100, 0b00001100, 
        0b11111000, 0b10111000, 0b01111000, 0b00111000, 0b11101000, 0b10101000, 0b01101000, 0b00101000, 0b11011000, 0b10011000, 0b01011000, 0b00011000, 0b11001000, 0b10001000, 0b01001000, 0b00001000,
        0b11110100, 0b10110100, 0b01110100, 0b00110100, 0b11100100, 0b10100100, 0b01100100, 0b00100100, 0b11010100, 0b10010100, 0b01010100, 0b00010100, 0b11000100, 0b10000100, 0b01000100, 0b00000100,
        0b11110000, 0b10110000, 0b01110000, 0b00110000, 0b11100000, 0b10100000, 0b01100000, 0b00100000, 0b11010000, 0b10010000, 0b01010000, 0b00010000, 0b11000000, 0b10000000, 0b01000000, 0b00000000
    };

    ASSERT (src->nbits == src->nwords * 64, "expecting full words, bits->nwords=%"PRIu64" and bits->num_of_bit=%"PRIu64,
            src->nwords, src->nbits);

    ASSERT0 (src->nbits == dst->nbits && src->nwords == dst->nwords, "expecting src and dst to have the same number of bits and words");

    ASSERT0 (src_start_base % 32 == 0 && max_num_bases % 32 == 0, "invalid start_base or num_bases");

    #define REV_COMP(w) ((rev_comp_table[(w)       & 0xff] << 56) | \
                         (rev_comp_table[(w >>  8) & 0xff] << 48) | \
                         (rev_comp_table[(w >> 16) & 0xff] << 40) | \
                         (rev_comp_table[(w >> 24) & 0xff] << 32) | \
                         (rev_comp_table[(w >> 32) & 0xff] << 24) | \
                         (rev_comp_table[(w >> 40) & 0xff] << 16) | \
                         (rev_comp_table[(w >> 48) & 0xff] << 8)  | \
                         (rev_comp_table[(w >> 56) & 0xff]        ))

    uint64_t after_word = MIN_(src->nwords, (src_start_base + max_num_bases) / 32); // 32 nucleotides in a word

    for (uint64_t i=src_start_base / 32; i < after_word; i++)
        dst->words[dst->nwords-1 - i] = REV_COMP (src->words[i]);
}

// reverse-complements an ACGT bit array in-place
void bits_reverse_complement_in_place (BitsP bits) 
{
    DEBUG_VALIDATE(bits);

    if (!bits->nbits) return;
    _reverse_region (bits, 0, bits->nbits);

    uint8_t *bytes = (uint8_t *)bits->words;
    uint64_t nbytes = bits->nwords * sizeof(uint64_t);

    for (uint64_t i=0; i < nbytes; i++) {
        uint8_t b = bytes[i];
        uint8_t b0 = b & 0b11;
        uint8_t b1 = (b>>2) & 0b11;
        uint8_t b2 = (b>>4) & 0b11;
        uint8_t b3 = (b>>6) & 0b11;

        #define CONVERT(x) ((x)==0b00 ? 0b11 : (x)==0b11 ? 0b00 : (x))

        // TO DO - we can do this faster with a lookup table
        bytes[i] = CONVERT(b0) | (CONVERT(b1) << 2) | (CONVERT(b2) << 4) | (CONVERT(b3) << 6);
    }

    bits_clear_excess_bits_in_top_word (bits);
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

    bits_clear_excess_bits_in_top_word (bits);
}

// get number of bits in an array, excluding trailing zeros (divon)
uint64_t bits_effective_length (BitsP bits)
{
    int64_t i;
    for (i=bits->nwords-1; i>=0; i--) 
        if (bits->words[i]) break; // break on last non-zero word

    if (i == -1) return 0;

    uint64_t last_non_zero_word = bits->words[i];

    return 64 * (1+i) - leading_zeros(last_non_zero_word);
}

// create words - if the word is not aligned to the bitmap word boundaries, and hence spans 2 bitmap words, 
// we take the MSb's from the left word and the LSb's from the right word, to create shift_1 (divon)
static inline uint64_t _bits_combined_word (uint64_t word_a, uint64_t word_b, uint8_t shift)
{
    uint64_t first_word_msb  = word_a >> shift; // first word's MSb are the LSb of our combined word
    uint64_t second_word_lsb = (word_b & bitmask64 (shift)) << (64-shift); // second word's LSb are the MSb of our combined word
    return first_word_msb | second_word_lsb; 
}

// calculate the number of bits that are different between two Bits' at arbitrary positions (divon)
uint32_t bits_hamming_distance (ConstBitsP bits_1, uint64_t index_1, 
                                ConstBitsP bits_2, uint64_t index_2, 
                                uint64_t len)
{
    const uint64_t *words_1 = &bits_1->words[index_1 >> 6];
    uint8_t shift_1 = index_1 & bitmask64(6); // word 1 contributes (64-shift) most-significant bits, word 2 contribute (shift) least significant bits

    const uint64_t *words_2 = &bits_2->words[index_2 >> 6];
    uint8_t shift_2 = index_2 & bitmask64(6); // word 1 contributes (64-shift) most-significant bits, word 2 contribute (shift) least significant bits

    uint64_t word=0;
    uint32_t nonmatches=0; 
    uint32_t nwords = roundup_bits2words64 (len);

    for (uint32_t i=0; i < nwords; i++) {

        uint64_t word_1 = shift_1 ? _bits_combined_word (words_1[i], words_1[i+1], shift_1) : words_1[i];
        uint64_t word_2 = shift_2 ? _bits_combined_word (words_2[i], words_2[i+1], shift_2) : words_2[i];
        
        word = word_1 ^ word_2; // xor the words - resulting in 1 in each position they differ and 0 where they're equal

        // we count the number of different bits. note that identical nucleotides results in 2 equal bits while
        // different nucleotides results in 0 or 1 equal bits (a 64 bit word contains 32 nucleotides x 2 bit each)
        nonmatches += __builtin_popcountll (word);
    }
    
    // remove non-matches due to the unused part of the last word
    if (len % 64)
        nonmatches -= __builtin_popcountll (word & ~bitmask64 (len % 64));

    return nonmatches; // this is the "hamming distance" between the two Bits' - number of non-matches
}

void LTEN_bits (BitsP bits)
{
#ifndef __LITTLE_ENDIAN__
    int64_t num_words = (bits->nbits >> 6) + ((bits->nbits & 0x3f) != 0);
    
    for (uint64_t i=0; i < num_words; i++)
        bits->words[i] = LTEN64 (bits->words[i]);
#endif
}
