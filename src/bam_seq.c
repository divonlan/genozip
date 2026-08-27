// ------------------------------------------------------------------
//   sam_bam_seq.c - functions for handling BAM binary sequence format 
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#ifdef __x86_64__
#include <immintrin.h>
#endif

const char bam_base_codes[16] = "=ACMGRSVTWYHKDBN";

rom bam_seq_display (bytes seq, uint32_t l_seq) // caller should free memory
{
    char *str = MALLOC (l_seq + 2);

    for (uint32_t i=0; i < (l_seq+1)/2; i++) {
        str[i*2]   = bam_base_codes[seq[i] >> 4];
        str[i*2+1] = bam_base_codes[seq[i] & 0xf];
    }

    str[l_seq] = 0;
    return str;
}

// re-writes BAM format SEQ into textual SEQ
void bam_seq_to_sam (VBlockP vb, bytes𐤐 bam_seq, 
                     uint32_t seq_len,       // bases, not bytes
                     bool start_mid_byte,    // ignore first nibble of bam_seq (seq_len doesn't include the ignored nibble)
                     bool test_final_nibble, // if true, we test that the final nibble, if unused, is 0, and warn if not
                     BufferP out,            // appends to end of buffer - caller should allocate seq_len+1 (+1 for last half-byte) 
                     bool is_from_zip_cb)    // don't account for time when codec-compressing, as the codecs account for their own time
{
    START_TIMER;
        
    ASSERT (out->len32 + seq_len + 2 <= out->size, "%s: out allocation too small", LN_NAME);

    if (!seq_len) {
        BNXTc (*out) = '*';
        return;        
    }

    // we implement "start_mid_byte" by converting the redudant base too, but starting 1 character before in the buffer 
    char save = 0;
    if (start_mid_byte) {
        out->len32--;
        seq_len++;
        save = *BAFTc(*out); // this is the byte we will overwrite, and recover it later. possibly, the fence if the buffer is empty;
    }
    
    unaligned_uint16_t *restrict sam_seq = (unaligned_uint16_t *)BAFTc(*out);
    uint32_t num_bytes = (seq_len + 1) / 2;
    uint32_t i=0;

#ifdef __x86_64__  // note: AVX2 support enforced by arch_initialize   
    const __m256i bam_base_codes_m256i = _mm256_setr_epi8 (
        '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
        '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
    );

    const __m256i mask_low = _mm256_set1_epi8 (0x0F);

    // process 32 BAM bytes (64 bases) per iteration
    for (; i + 32 <= num_bytes; i += 32) {
        __m256i raw = _mm256_loadu_si256 ((const __m256i*)(bam_seq + i)); // unaligned load

        // extract high (first base) and low (second base) nibbles
        __m256i high_nibbles = _mm256_and_si256 (_mm256_srli_epi16(raw, 4), mask_low);
        __m256i low_nibbles  = _mm256_and_si256 (raw, mask_low);

        // hardware vector lookup
        __m256i ascii_first  = _mm256_shuffle_epi8 (bam_base_codes_m256i, high_nibbles);
        __m256i ascii_second = _mm256_shuffle_epi8 (bam_base_codes_m256i, low_nibbles);

        // interleave (unpacking in 128-bit lanes)
        __m256i unp_lo = _mm256_unpacklo_epi8 (ascii_first, ascii_second); // Lane0-Lo, Lane1-Lo
        __m256i unp_hi = _mm256_unpackhi_epi8 (ascii_first, ascii_second); // Lane0-Hi, Lane1-Hi

        // re-align 128-bit lanes into continuous 256-bit sequential order
        __m256i out0 = _mm256_permute2x128_si256 (unp_lo, unp_hi, 0x20); // Lane0-Lo | Lane0-Hi
        __m256i out1 = _mm256_permute2x128_si256 (unp_lo, unp_hi, 0x31); // Lane1-Lo | Lane1-Hi

        // store 64 bases (32x uint16_t words) - unaligned
        _mm256_storeu_si256 ((__m256i*)(sam_seq + i), out0);
        _mm256_storeu_si256 ((__m256i*)(sam_seq + i + 16), out1);
    }
#endif

    // remaining tail bytes (or all bytes if not using AVX2): 2 bases at a time
    for (; i < num_bytes; i++)
        sam_seq[i] =  (uint16_t)bam_base_codes[bam_seq[i] >> 4]
                   | ((uint16_t)bam_base_codes[bam_seq[i] & 0x0F] << 8);

    if (start_mid_byte) {
        *BAFTc(*out) = save;
        out->len32 += seq_len;      
        seq_len--;
    }
    else
        out->len32 += seq_len;

    ASSERTW (!test_final_nibble || !(seq_len % 2) || (*BAFTc (*out)=='='), 
             _WRN "%s: bam_seq_to_sam: expecting the unused lower 4 bits of last seq byte in an odd-length seq_len=%u to be 0, but its not. This will cause an incorrect digest",
             LN_NAME, seq_len);

    *BAFTc(*out) = 0; // nul-terminate after end of seq
    
    if (!is_from_zip_cb) COPY_TIMER(bam_seq_to_sam);
}

/*
// compare a sub-sequence to a full sequence and return true if they're the same. Sequeneces in BAM format.
bool bam_seq_has_sub_seq (bytes full_seq, uint32_t full_seq_len, 
                          bytes sub_seq,  uint32_t sub_seq_len, uint32_t start_base) // lengths are in bases, not bytes
{
    // easy case: similar byte alignment
    if (!(start_base & 1)) {
        if (memcmp (sub_seq, &full_seq[start_base/2], sub_seq_len/2)) return false; // mismatch in first even number of bases
        if (!(sub_seq_len & 1)) return true; // even number of bases

        uint8_t last_base_full = full_seq[(start_base + sub_seq_len - 1)/2] >> 4;
        uint8_t last_base_sub  = sub_seq[(sub_seq_len - 1)/2] >> 4;
        return last_base_full == last_base_sub;
    }

    // not byte-aligned
    else { 
        for (uint32_t sub_base_i=0; sub_base_i < sub_seq_len; sub_base_i++) {
    
            uint8_t base_sub = (sub_base_i % 2) ? (sub_seq[sub_base_i/2] & 15) : (sub_seq[sub_base_i/2] >> 4);
                        
            uint32_t full_base_i = start_base + sub_base_i;
            uint8_t base_full = (full_base_i % 2) ? (full_seq[full_base_i/2] & 15) : (full_seq[full_base_i/2] >> 4);

            if (base_sub != base_full) return false;
        }
        return true; // all bases are the same
    }
}
*/
static void bam_seq_copy (uint8_t *dst, bytes src, 
                          uint32_t src_start_base, uint32_t n_bases) // bases, not bytes
{
    src += src_start_base / 2;

    if (src_start_base & 1) {
        for (uint32_t i=0; i < n_bases / 2; i++) 
            *dst++ = ((src[i] & 0x0f) << 4) | (src[i+1] >> 4); 
        
        if (n_bases & 1)
            *dst = (src[n_bases / 2] & 0x0f) << 4;
    }
    else {
        memcpy (dst, src, (n_bases+1)/2);

        if (n_bases & 1)
            dst[n_bases/2] &= 0xf0; // keep the high nibble only (the first BAM base of this byte)
    }
}

// C<>G A<>T ; IUPACs: R<>Y K<>M B<>V D<>H W<>W S<>S N<>N
static void bam_seq_revcomp_in_place (uint8_t *seq, uint32_t n_bases)
{                                 // Was:  =    A    C    M    G    R    S    V    T    W    Y    H    K    D    B    N                          
                                  // Comp: =    T    G    K    C    Y    S    B    A    W    R    D    M    H    V    N
    static const uint8_t bam_comp[16] = { 0x0, 0x8, 0x4, 0xc, 0x2, 0xa, 0x6, 0xe, 0x1, 0x9, 0x5, 0xd, 0x3, 0xb, 0x7, 0xf };

    for (int32_t i=0, j=n_bases-1; i < n_bases/2; i++, j--) {
        
        uint8_t b1 = (i&1) ? (seq[i/2] & 0xf) : (seq[i/2] >> 4);
        uint8_t b1c = bam_comp[b1];

        uint8_t b2 = (j&1) ? (seq[j/2] & 0xf) : (seq[j/2] >> 4);
        uint8_t b2c = bam_comp[b2];

        seq[i/2] = (i&1) ? (b2c | (seq[i/2] & 0xf0)) : ((b2c << 4) | (seq[i/2] & 0x0f));
        seq[j/2] = (j&1) ? (b1c | (seq[j/2] & 0xf0)) : ((b1c << 4) | (seq[j/2] & 0x0f));
    }
}

// returns length in bytes 
uint32_t sam_seq_copy (char *dst, rom src, uint32_t src_start_base, uint32_t n_bases, 
                       bool revcomp, bool is_bam_format)
{
    if (!is_bam_format) {
        if (revcomp)
            str_revcomp (dst, src, n_bases);
        else
            memcpy (dst, src, n_bases);
    }
    
    else {
        bam_seq_copy ((uint8_t*)dst, (const uint8_t*)src, src_start_base, n_bases);

        if (revcomp)
            bam_seq_revcomp_in_place ((uint8_t*)dst, n_bases);
    } 

    return is_bam_format ? (n_bases + 1) / 2 : n_bases;
}

/*
// like bam_seq_has_sub_seq, but sub_seq is reverse complemented, and start_base is relative to END of full_seq
bool bam_seq_has_sub_seq_revcomp (bytes full_seq, uint32_t full_seq_len, 
                                  bytes sub_seq,  uint32_t sub_seq_len, uint32_t start_base)
{
    bam_seq_revcomp_in_place ((uint8_t*)STRa(sub_seq));

    bool same = bam_seq_has_sub_seq (STRa(full_seq), STRa(sub_seq), full_seq_len - start_base - sub_seq_len);

    // revcomp-back, if we have to
    bam_seq_revcomp_in_place ((uint8_t*)STRa(sub_seq));

    return same;
}
*/