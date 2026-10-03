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
#ifdef __x86_64__
#include <immintrin.h>
#endif
#include "genozip.h"
#include "endianness.h"
#include "bits.h"
#include "buffer.h"
#include "reference.h"

// note: word can go partially beyond nbits, in which case the excess bits are ignored 
void _set_word (BitsP bits, uint64_t start, uint64_t word)
{
    DEBUG_VALIDATE_BITS(bits);
    ASSERT (start < bits->nbits, "bits=%s: Expecting start=%"PRIu64" < nbits=%"PRIu64, unר(bits->nameר), start, bits->nbits);
    
    uint64_t *words = bits->words;
    uint64_t word_index       = bits_wrd(start);
    word_offset_t word_offset = bits_idx(start);

    if (word_offset == 0) {
        words[word_index] = word;

        // reset unused top bits to 0 if they might have been modified
        if (word_index + 1 == bits->nwords) 
            bits_clear_excess_bits_in_top_word (bits);
    }
    
    else {
        words[word_index] = (word << word_offset) |
                            (words[word_index] & bitmask64(word_offset));

        if (word_index+1 < bits->nwords) { // if last part of the word goes beyond nwords, we drop it
            words[word_index+1] = (word >> (BITS_IN_WORD - word_offset)) |
                                  (words[word_index+1] & (UINT64_MAX << word_offset));

            // reset unused top bits to 0 if they might have been modified
            if (word_index + 1 == bits->nwords || word_index + 2 == bits->nwords)
                bits_clear_excess_bits_in_top_word (bits);
        }
    }
}

// set all the bits in a region
void bits_set_region (BitsP bits, uint64_t start, uint64_t len)
{
    if (!len) return;

    ASSERT (start + len <= bits->nbits, "bits=%s: Expecting: start=%"PRId64" + len=%"PRId64" <= nbits=%"PRId64"",
            unר(bits->nameר), start, len, bits->nbits); 

    uint64_t first_word   = bits_wrd(start);
    uint64_t last_word    = bits_wrd(start + len - 1);
    word_offset_t foffset = bits_idx(start);
    word_offset_t loffset = bits_idx(start + len - 1);

    if (__builtin_expect (first_word == last_word, false)) {
        uint64_t mask = bitmask64(len) << foffset;
        bits->words[first_word] |= mask; 
    }
    
    else {
        bits->words[first_word] |= ~bitmask64(foffset);     // set first word
        memset (&bits->words[first_word + 1], 0xff, (last_word - first_word - 1) * sizeof(uint64_t));
        bits->words[last_word] |= bitmask64(loffset+1);     // set last word
    }

    DEBUG_VALIDATE_BITS (bits);
}

// Clear all the bits in a region
void bits_clear_region (BitsP bits, uint64_t start, uint64_t len)
{
    if (!len) return; // nothing to do 

    ASSERT (start + len <= bits->nbits, "bits=%s: Expecting: start=%"PRId64" + len=%"PRId64" <= nbits=%"PRId64,
            unר(bits->nameר), start, len, bits->nbits); // divon fixed bug

    uint64_t first_word   = bits_wrd(start);
    uint64_t last_word    = bits_wrd(start + len - 1);
    word_offset_t foffset = bits_idx(start);
    word_offset_t loffset = bits_idx(start + len - 1);

    if (__builtin_expect (first_word == last_word, false)) {
        uint64_t mask = bitmask64(len) << foffset;
        bits->words[first_word] &= ~mask;
    }
    
    else {
        bits->words[first_word] &= bitmask64(foffset);      // clear first word
        memset (&bits->words[first_word + 1], 0, (last_word - first_word - 1) * sizeof(uint64_t));
        bits->words[last_word] &= ~bitmask64(loffset+1);    // clear last word
    }

    DEBUG_VALIDATE_BITS (bits);
}

//
// Number of bits set
//

// true if all bits in bit array are set
bool bits_is_fully_set (ConstBitsP bits)
{
    DEBUG_VALIDATE_BITS(bits);

    if (!bits->nbits) return true; // trivially true

    uint64_t i=0;
    uint64_t full_words = bits->nbits / BITS_IN_WORD;
    word_offset_t partial_top_word_nbits = bits->nbits % BITS_IN_WORD;

#ifdef __x86_64__
    const __m256i all_set = _mm256_set1_epi64x (-1);
    
    for (; i+4 <= full_words; i+=4) {
        __m256i word = _mm256_loadu_si256 ((const __m256i *)(bits->words + i));

        if (!_mm256_testc_si256 (word, all_set))
            return false;
    }
#endif

    for (; i < full_words; i++)
        if (bits->words[i] != UINT64_MAX)
            return false;

    // test partial top word
    if (partial_top_word_nbits && 
        (bits->words[full_words] & bitmask64_(partial_top_word_nbits)) != bitmask64_(partial_top_word_nbits))
        return false;

    return true;
}

// true if all bits in bit array are clear
bool bits_is_fully_clear (ConstBitsP bits)
{
    uint64_t full_words = bits->nbits / BITS_IN_WORD;
    word_offset_t partial_top_word_nbits = bits->nbits % BITS_IN_WORD;

    if (full_words && !str_is_zero ((rom)bits->words, full_words * sizeof(uint64_t))) 
        return false;

    if (partial_top_word_nbits && (bits->words[full_words] & bitmask64_(partial_top_word_nbits))) 
        return false;

    return true;
}

// Get the number of bits set (hamming weight)
uint64_t bits_num_set_bits (ConstBitsP bits)
{
    if (__builtin_expect(!bits->nbits, 0)) return 0;

    uint64_t full_words = bits->nbits / BITS_IN_WORD;
    word_offset_t partial_top_word_nbits = bits->nbits % BITS_IN_WORD;

    uint64_t num_of_bits_set = 0;

    // full words
    for (uint64_t i=0; i < full_words; i++)
        num_of_bits_set += __builtin_popcountll (bits->words[i]);

    // last partial word
    if (partial_top_word_nbits)
        num_of_bits_set += __builtin_popcountll (bits->words[full_words] & bitmask64_(partial_top_word_nbits));
    
    return num_of_bits_set;
}

// Get the number of bits not set (1 - hamming weight)
uint64_t bits_num_clear_bits (ConstBitsP bits)
{
    return bits->nbits - bits_num_set_bits(bits);
}

uint64_t bits_num_set_bits_region (ConstBitsP bits, uint64_t start, uint64_t length)
{
    if (length == 0) return 0;

    ASSERT (start + length <= bits->nbits, "bits=%s: out of range: execpting: start=%"PRIu64" + length=%"PRIu64" <= nbits=%"PRIu64, 
            unר(bits->nameר), start, length, bits->nbits);

    uint64_t first_word = bits_wrd(start);
    uint64_t last_word  = bits_wrd(start+length-1);
    word_offset_t foffset = bits_idx(start);

    uint64_t num_of_bits_set = 0;

    if (first_word == last_word) {
        uint64_t mask = bitmask64_(length) << foffset;
        num_of_bits_set += __builtin_popcountll (bits->words[first_word] & mask);
    }
    else {
        word_offset_t loffset  = bits_idx(start+length-1);
    
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
    ASSERT (offset < bits->nbits, "bits=%s: expecting offset(%"PRId64") < bits->nbits(%"PRId64")", unר(bits->nameר), offset, bits->nbits); \
    if (bits->nbits == 0 || offset >= bits->nbits) { return false; } \
    \
    /* Find first word that is greater than zero */ \
    uint64_t i = bits_wrd(offset); \
    uint64_t w = GET(bits->words[i]) & ~bitmask64(bits_idx(offset)); \
    \
    while (1) { \
        if (w > 0) { \
            uint64_t pos = i * BITS_IN_WORD + trailing_zeros(w); \
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
    ASSERT (offset <= bits->nbits, "bits=%s: expecting offset=%"PRIu64" <= nbits=%"PRIu64, unר(bits->nameר), offset, bits->nbits); \
    if (bits->nbits == 0 || offset == 0) { return false; } \
    \
    /* Find prev word that is greater than zero */ \
    uint64_t i = bits_wrd(offset-1); \
    uint64_t w = GET(bits->words[i]) & bitmask64(bits_idx(offset-1)+1); \
    \
    if (w > 0) { *result = (i+1) * BITS_IN_WORD - leading_zeros(w) - 1; return true; } \
    \
    /* i is unsigned so have to use break when i == 0 */ \
    for (--i; i != UINT64_MAX; i--) { \
        w = GET(bits->words[i]); \
        if (w > 0) { \
            *result = (i+1) * BITS_IN_WORD - leading_zeros(w) - 1; \
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

// move a range of bits to a lower index within the same array
// note: this function writes strictly to the dstindex, length range and is guaranteed to
//       never pollute the flanking bits
void bits_sink_range (ConstBits𐤐 bits, uint64_t dstindx, uint64_t srcindx, uint64_t length,
                      bool src_may_exceed_nbits) // if true, src is allowed to go until the end of the last word, even beyond nbits
{
    if (__builtin_expect (srcindx == dstindx || !length, false)) return;

    // validate input
    uint64_t src_limit_bits = src_may_exceed_nbits ? (bits->nwords * BITS_IN_WORD) : bits->nbits;
    
    ASSERT (srcindx <= src_limit_bits && length <= src_limit_bits - srcindx, 
            "bits=%s: srcindx=%"PRIu64" + length=%"PRIu64" > src_limit_bits=%"PRIu64, unר(bits->nameר), srcindx, length, src_limit_bits);
    
    ASSERT (dstindx < srcindx, "bits=%s: dstindx(%"PRIu64") >= srcindx(%"PRIu64")", unר(bits->nameר), dstindx, srcindx);

    ASSERT (dstindx <= bits->nbits && length <= bits->nbits - dstindx, 
            "bits=%s: dstindx=%"PRIu64" + length=%"PRIu64" > nbits=%"PRIu64, unר(bits->nameר), dstindx, length, bits->nbits);
    
     // verify that bits will no longer exceed after sinking
    ASSERT (!src_may_exceed_nbits || src_limit_bits - bits->nbits <= srcindx - dstindx, 
            "Sunk bits=%s should not exceed nbits=%"PRIu64". src_limit_bits=%"PRIu64" srcindx=%"PRIu64" dstindx=%"PRIu64,
            unר(bits->nameר), bits->nbits, src_limit_bits, srcindx, dstindx);

    uint64_t *words = bits->words, 
    nwords   = bits->nwords,
    start_w  = bits_wrd(dstindx),
    end_w    = bits_wrd(dstindx + length - 1),
    dst_off  = bits_idx(dstindx), // offset of dstindx within its starting word
    shift_r  = (srcindx - dst_off) & (BITS_IN_WORD - 1), // shift relative to physical destination word 0
    shift_l  = (BITS_IN_WORD - shift_r) & 63,
    src_w    = bits_wrd(srcindx - dst_off), // Physical source word corresponding to destination word start_w
    curr_src = (src_w < nwords) ? words[src_w] : 0; // Pre-fetch initial physical source word

    for (uint64_t w=start_w; w <= end_w; w++, src_w++) {
        uint64_t next_src = (src_w + 1 < nwords) ? words[src_w + 1] : 0;

        // shift_r is constant in this loop, so branch prediction will work well here
        uint64_t aligned_src; 
        if (shift_r) aligned_src = (curr_src >> shift_r) | (next_src << shift_l);
        else         aligned_src = curr_src;

        uint64_t cur_word_bit_start = w * BITS_IN_WORD,
        write_start  = (dstindx > cur_word_bit_start) ? (dstindx - cur_word_bit_start) : 0,
        write_end    = (dstindx + length < cur_word_bit_start + BITS_IN_WORD) ? (dstindx + length - cur_word_bit_start) : BITS_IN_WORD,
        bits_to_copy = write_end - write_start,
        mask         = bitmask64(bits_to_copy) << write_start;

        words[w] = bitmask_merge(aligned_src, words[w], mask);

        curr_src = next_src;
    }

    if (!src_may_exceed_nbits) DEBUG_VALIDATE_BITS(bits);
}

// copy a range from bits to another (must be different bits)
void bits_copy (BitsP dst, uint64_t dstindx, ConstBitsP src, uint64_t srcindx, uint64_t length)
{
    DEBUG_VALIDATE_BITS(src);
    if (__builtin_expect(length == 0, 0)) return;

    ASSERT (src != dst && src->words != dst->words, "src==dst=%s unsuppported, use bits_sink_range", unר(src->nameר));
    ASSERT (dstindx + length <= dst->nbits, "dst=%s: dstindx(%"PRIu64") + length(%"PRIu64") > dst->nbits(%"PRIu64")", unר(dst->nameר), dstindx, length, dst->nbits);
    ASSERT (srcindx + length <= src->nbits, "src=%s: srcindx(%"PRIu64") + length(%"PRIu64") > src->nbits(%"PRIu64")", unר(src->nameר), srcindx, length, src->nbits);

    const uint64_t *restrict src_words = src->words;
    uint64_t *restrict dst_words = dst->words;
    uint64_t src_nwords = src->nwords;

    // STAGE 1: Align Destination to Word Boundary (Head)
    uint32_t dst_idx = bits_idx(dstindx);

    if (dst_idx != 0) {
        uint32_t head_bits = (uint32_t)MIN_((uint64_t)(64 - dst_idx), length);

        // Fetch via src_words restrict pointer + local mask
        uint64_t src_w = _get_word_(src_words, src_nwords, srcindx) & bitmask64(head_bits);
        uint64_t dst_w = dst_words[bits_wrd(dstindx)];

        uint64_t mask = bitmask64(head_bits) << dst_idx;
        dst_words[bits_wrd(dstindx)] = (dst_w & ~mask) | (src_w << dst_idx);

        srcindx += head_bits;
        dstindx += head_bits;
        length  -= head_bits;    
    }

    // STAGE 2: Direct 64-Bit Stores (Middle): dstindx is now strictly word-aligned!
    uint64_t full_words = length / BITS_IN_WORD;
    if (full_words > 0) {
        uint64_t *restrict dst_ptr = &dst_words[bits_wrd(dstindx)];

        if (srcindx % 64)
            for (uint64_t i = 0; i < full_words; i++) {
                *dst_ptr++ = _get_word_(src_words, src_nwords, srcindx);
                srcindx += BITS_IN_WORD;
            }
        else {
            memcpy (&dst_words[bits_wrd(dstindx)], 
                    &src_words[bits_wrd(srcindx)], 
                    full_words * sizeof(uint64_t));
            srcindx += full_words * BITS_IN_WORD;
        }

        dstindx += full_words * BITS_IN_WORD;
        length  %= BITS_IN_WORD;
    }

    // STAGE 3: Handle Remaining Bits (Tail)
    if (length > 0) {
        uint64_t src_w = _get_word_(src_words, src_nwords, srcindx) & bitmask64(length); // note: cannot use bits_get_wordn because it would violate restrict
        uint64_t dst_w = dst_words[bits_wrd(dstindx)];

        uint64_t mask = bitmask64(length);
        dst_words[bits_wrd(dstindx)] = (dst_w & ~mask) | src_w;

        dstindx += length;
    }

    DEBUG_VALIDATE_BITS(dst);
}

void bits_overlay (BitsP overlaid_bits, BitsP regular_bits, uint64_t start, uint64_t nbits, rom name)
{
    ASSERT (start % 64 == 0, "start=%"PRIu64" must be a multiple of 64 (name=%s)", start, name);
    ASSERT (start + nbits <= regular_bits->nbits, "start(%"PRIu64") + nbits(%"PRIu64") <= regular_bits->nbits(%"PRIu64") (name=%s)",
            start, nbits, regular_bits->nbits, name);

    uint64_t word_i = start / 64;
    *overlaid_bits = (Bits){ .type   = BITS_OVERLAY ,
                             .nbits  = nbits,
                             .words  = &regular_bits->words[word_i],
                             .nwords = roundup_bits2words64 (nbits),
                             .size   = roundup_bits2bytes64 (nbits),
                             .nameר  = ר(name) };

    // note: we can't clear top bits here, because we might overlay on a read-only shm
} 

// convert a bit array to a byte array of values 0/1
void bits_bit_to_byte (uint8_t *restrict dst, ConstBitsP src_bits, uint64_t src_bit, uint32_t num_bits)
{
    ASSERT (src_bit + num_bits <= src_bits->nbits, "src=%s: Expecting src_bit=%"PRIu64" + num_bits=%u) <= nbits=%"PRIu64,
            unר(src_bits->nameר), src_bit, num_bits, src_bits->nbits);

    while (num_bits) {
        uint64_t word = _get_word (src_bits, src_bit);
        uint32_t n = MIN_(num_bits, 64);
        const uint32_t consumed = n;

        while (n >= 8) {
#ifdef __x86_64__
            // deposit the 8(=popcount(2nd arg)) LSb of word into the bit locations indicated by the 2nd arg, which happen to be the LSb of each byte of the output word - effectively creating an uint8_t[8] array of 0/1
            *(unaligned_uint64_t *)dst = _pdep_u64 (word, 0x0101010101010101ULL); 
#else
            for (uint32_t i=0; i < 8; i++)
                dst[i] = (word >> i) & 1;
#endif
            dst += 8;
            word >>= 8;
            n -= 8;
        }

        // up to 7 tail bits
        while (n) {
            *dst++ = word & 1;
            word >>= 1;
            n--;
        }

        src_bit  += consumed;
        num_bits -= consumed;
    }
}

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

void bits_reverse (BitsP bits)
{
    if (__builtin_expect (!bits->nbits, 0)) return;

    uint64_t nwords = bits->nwords;
    uint64_t *words = bits->words;

    // Pass 1: Swap word positions and reverse bits within each word
    for (uint64_t i=0, j=nwords-1; i < j; i++, j--) {
        uint64_t tmp = words[i];
        words[i] = _reverse_word(words[j]);
        words[j] = _reverse_word(tmp);
    }

    // Reverse middle word if nwords is odd
    if (nwords & 1) 
        words[nwords / 2] = _reverse_word(words[nwords / 2]);

    // Pass 2: Shift entire range from [shift ... shift + nbits) down to [0 ... nbits)
    uint64_t shift = (BITS_IN_WORD - bits->nbits % BITS_IN_WORD) % BITS_IN_WORD;
    bits_sink_range (bits, 0, shift, bits->nbits, true); 

    // clear tops bits which were polluted when reversing the top word
    bits_clear_excess_bits_in_top_word (bits);

    DEBUG_VALIDATE_BITS(bits);
}

// resizes an array to a certain number of bits (longer or shorter) - without reallocating
uint64_t bits_resize (BitsP bits, uint64_t new_nbits)
{
    ASSERT (new_nbits <= bits->size * 8, "bits=%s: no room to extend the bitmap: nbits=%"PRIu64", num_new_bits=%"PRId64", bits->size=%"PRIu64, 
            unר(bits->nameר), bits->nbits, new_nbits, (uint64_t)bits->size);

    uint64_t next_bit = bits->nbits;

    bits->nbits  = new_nbits;
    uint64_t old_nwords = bits->nwords;
    bits->nwords = roundup_bits2words64 (new_nbits);   

    if (bits->nwords > old_nwords)
        bits->words[bits->nwords-1] = 0; // new word: zero entire word to avoid read-modify-write of an uninitialized memory
    else
        // note: these is still a possibility of a read-modify-write here (valgrind error) if we shorten bits, and the new top word was never initialized
        bits_clear_excess_bits_in_top_word (bits); // note: if shortening, abandoned data is ruined

    DEBUG_VALIDATE_BITS (bits); 
    return next_bit;
}

// removes flanking bits on boths sides, shrinking bits
void bits_remove_flanking (BitsP bits, uint64_t lsb_flanking, uint64_t msb_flanking) // added by divon
{
    ASSERT (bits->nbits >= lsb_flanking + msb_flanking, "bits=%s: Expecting nbits=&"PRIu64" >= lsb_flanking=%"PRIu64" + msb_flanking=%"PRIu64,
            unר(bits->nameר), bits->nbits, lsb_flanking, msb_flanking);

    uint64_t new_nbits = bits->nbits - lsb_flanking - msb_flanking;
    bits_sink_range (bits, 0, lsb_flanking, new_nbits, false);

    bits_resize (bits, new_nbits);
}

void bits_add_bit (BitsP bits, int64_t new_bit) 
{
    ASSERT (bits->nbits < bits->size * 8, "bits=%s: no room to extend the bitmap", unר(bits->nameר));
    bits->nbits++;     
    if (bits->nbits % 64 == 1) { // starting a new word                
        bits->nwords++;
        bits->words[bits->nwords-1] = new_bit; // LSb is as requested, other 63 unused top bits are 0
    } 
    else
        bits_assign (bits, bits->nbits-1, new_bit);  
}
  
#ifndef __LITTLE_ENDIAN__
void LTEN_bits (BitsP bits)
{
    int64_t num_words = (bits->nbits >> 6) + ((bits->nbits & 0x3f) != 0);
    
    for (uint64_t i=0; i < num_words; i++)
        bits->words[i] = LTEN64 (bits->words[i]);
}
#endif

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
    const uint64_t *after_2 = bits_2->words + bits_2->nwords;
    uint64_t diff=0;
    int shift_2 = index_2 & 63;
    uint32_t distance=0; // number of non-matching bits

    // if bits_1 goes beyond the end of bits_2 we can't calculate distance
    if (index_2 > bits_2->nbits || bits_1->nbits > (bits_2->nbits - index_2)/*remaining length in bits_2*/)
        return bits_1->nbits; // maximum distance = complete misfit

    uint64_t next_w2 = *words_2; 
    for (const uint64_t *w1=words_1, *w2=words_2; w1 < after_1; w1++, w2++) {
        uint64_t this_w2 = next_w2;
        next_w2 = (w2 + 1 < after_2) ? *(w2 + 1) : 0;
        // note: if last bits_1 word is partial, we still calculate it in its entirety and correct later. However, in case we reach 
        // the end of bits_2 and have shift_2, we will access the word beyond the end of the bits_2->words. counting on words being a Buffer to not cause a segfault.
        uint64_t w2_value = _bits_combined_word (this_w2, next_w2, shift_2);
        diff = *w1 ^ w2_value; 
        distance += __builtin_popcountll (diff); // note: expected to be _mm_popcnt_u64 with SSE4.2
    }

    // remove distance due to the unused part of the last bits_1 word
    if (bits_1->nbits & 63)
        distance -= __builtin_popcountll (diff & ~bitmask64 (bits_1->nbits & 63));

    return distance; 
}

// Stringification

StrText1K bits_to_01_string (ConstBitsP bits, uint64_t start, 
                             int64_t length) // length -1 means "to the end of the bitmap"
{
    StrText1K s;

    if (start >= bits->nbits) 
        length = 0; // nothing to output
    
    else if (length < 0 || start + length > bits->nbits)
        length = bits->nbits - start; // until the end the bits
    
    MINIMIZE (length, sizeof(s)-1);
      
    for (uint64_t i=0; i < length; i++) 
        s.s[i] = bits_get(bits, start + i) ? '1' : '0';

    s.s[length] = '\0';
    return s;
}

// prints a 2bit bases array - A,C,G,T bases
StrText1K bits_to_ACGT_string (ConstBitsP bits, uint64_t start_base, uint64_t n_bases)
{
    StrText1K s;
    MINIMIZE (n_bases, sizeof(s)-1);

    uint64_t start_bit = start_base * 2;
    uint64_t n_bits = MIN_(n_bases * 2, bits->nbits - start_bit);
    decl_acgt_decode;

    char *next = s.s;
    for (uint64_t i=0; i < n_bits; i += 2) { 
        uint8_t b = bits_get2 (bits, start_bit + i);
        *next++ = acgt_decoder[b];
    }

    *next++ = '\0';
    return s;
}

// print word as bit array (LSb is first character '0' or '1')
StrText bits_word_to_01_string (uint64_t word, int n_bits/*0 to 64*/)
{
    StrText w;
    MINIMIZE (n_bits, BITS_IN_WORD);
    w.s[n_bits] = 0;

    for (word_offset_t i=0; i < n_bits; i++)
        w.s[i] = ((word >> i) & (uint64_t)0x1) == 0 ? '0' : '1';

    return w;
}

#ifdef DEBUG
void validate_bits (ConstBitsP bits, Caller caller)
{
    // Verify that its allocated
    ASSERT (bits->type != BUF_UNALLOCATED, "[%s:%u] Bits is not allocated", CALLERf);

    // Check num of words is correct
    uint64_t num_words = roundup_bits2words64 (bits->nbits);
    ASSERT (num_words == bits->nwords, "[%s:%u] bits=%s num of words wrong, nbits=%"PRIu64", expected_nwords=%"PRIu64", actual nwords=%"PRIu64, 
            CALLERf, unר(bits->nameר), bits->nbits, num_words, bits->nwords);

    // check that words is word-aligned
    ASSERT ((uintptr_t)bits->words % sizeof(uint64_t) == 0, 
            "bits=%s bits->words=%p is not word-aligned", unר(bits->nameר), bits->words);

    // Check top word is masked (only if not overlayed - the unused bits of top word don't belong to this bit array and might be used eg by another bit array in genome.ref/genome.is_set)
    if (bits->type == BUF_REGULAR) {
        word_offset_t top_bits = bits->nbits & 63; // 0 to 63
        if (top_bits) { // not full word (and also not empty bitmap)
            uint64_t tw = bits->nwords - 1;

            ASSERT (bits->words[tw] <= bitmask64_(top_bits), "[%s:%u] bits=%s Expected top %u bits in top word_i=%"PRIu64" to be 0, but word=%s\n", 
                    CALLERf, unר(bits->nameר), 64-top_bits, tw, bits_word_to_01_string(bits->words[tw], BITS_IN_WORD).s);
        }
    }
}
#endif
