// ------------------------------------------------------------------
//   comp_agct.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "lzma/7zTypes.h"
#include "lzma/LzmaEnc.h"
#include "lzma/LzmaDec.h"
#include "genozip.h"
#include "comp_private.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"
#include "endianness.h"

// -------------------------------------------------------------------------------------
// acgt stuff
// compress a sequence of A,C,G,T nucleotides - first squeeze into 2 bits and then LZMA.
// It's about 25X faster and slightly better compression ratio than LZMA
// -------------------------------------------------------------------------------------

// decoder of 2bit encoding of nucleotides
const char acgt_decode[4] = { 'A', 'C', 'G', 'T' };

// table to convert ASCII to ACGT encoding. A,C,G,T (lower and upper case) are encoded as 0,1,2,3 respectively, 
// and everything else (including N) is encoded as 0
const uint8_t acgt_encode[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 0
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 16
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 32
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 48
                                   0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,   // 64  A(65)->0 C(67)->1 G(71)->2
                                   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 80  T(84)->3
                                   0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,   // 96  a(97)->0 c(99)->1 g(103)->2
                                   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 112 t(116)->3
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 128
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


// table to convert ASCII to NON-ACGT encoding. The character is XORed with the entry in the table
static const uint8_t non_acgt_encode[256] = 
                                 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 0   -> XOR with 0 = stay unchanged
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 16
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 32
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 48
                                   0,'A',0,'C',0, 0, 0,'G',0, 0, 0, 0, 0, 0, 0, 0,   // 64  A(65), C(67), G(71) -> 0 (XORed with self)
                                   0, 0, 0, 0,'T',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 80  T(84)->0
                                   0, 'a'^1,0,'c'^1,0,0,0,'g'^1,0,0,0,0,0,0, 0, 0,   // 96  a(97), c(99), g(103)-> 1 (XORed with self XOR 1)
                                   0, 0, 0, 0,'t'^1,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0,   // 112 t(116)->1
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   // 128
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

// packing of an array A,G,C,T characters into a 2-bit BitArray, stored in vb->compress. Previous incomplete 64bit words in carry
// are used at the beginning of the resulting array, and any uncomplete 64bit word at the end, is stored back in carry
void comp_acgt_pack (VBlockP vb, const char *data, uint64_t data_len, unsigned bits_consumed, bool do_lten, bool do_lten_partial_final_word)
{
    START_TIMER;
    
    buf_alloc (vb, &vb->compressed, 
               (vb->compressed.len + roundup_bits2words64 (data_len * 2)) * sizeof (uint64_t), // note: len is in words
               2, buf_is_allocated (&vb->compressed) ? NULL : "compress", 0); // NULL if already allocated, to avoid overwriting param which is overlayed with BitArray->num_of_bits

    BitArray *packed = buf_get_bitarray (&vb->compressed);

    // remove consumed bits
    if (bits_consumed)
        bit_array_shift_right_shrink (packed, bits_consumed);

    // increase bit array to accomodate data
    uint64_t next_bit = packed->num_of_bits;
    packed->num_of_bits += data_len * 2;
    packed->num_of_words = roundup_bits2words64 (packed->num_of_bits);

    // pack nucleotides - each character is packed into 2 bits
    for (uint64_t i=0 ; i < data_len ; i++) {
        if (!IS_NUCLEOTIDE (data[i])) 
            vb->has_non_agct = true;

        uint8_t encoding = acgt_encode[(uint8_t)data[i]];
        bit_array_assign (packed, next_bit, encoding & 1);
        bit_array_assign (packed, next_bit + 1, (encoding & 2) >> 1);
        next_bit += 2;
    }

    // note: we store in Little Endian unlike the rest of the data that is in Big Endian, because LTEN keeps the nucleotides in their
    // original order, and improves compression ratio by about 2%
    if (do_lten)
        LTEN_bit_array (packed, do_lten_partial_final_word);

    COPY_TIMER (vb->profile.comp_acgt_pack)
}

void comp_acgt_pack_last_partial_word (VBlockP vb, ISeqInStream *instream)
{
    BitArray *packed = buf_get_bitarray (&vb->compressed);
    uint64_t bits_remaining = packed->num_of_bits - instream->bits_consumed;

    if (bits_remaining) { 
        
        ASSERT (bits_remaining >= 1 && bits_remaining <= 63, "Error in comp_acgt_pack_last_partial_word: Invalid bits_remaining%u", (uint32_t)bits_remaining);

        bit_array_shift_right_shrink (packed, instream->bits_consumed);
     
        instream->bits_consumed = 0;
        packed->num_of_bits  = packed->num_of_words = 0;
        packed->words[0]     = LTEN64 (packed->words[0]);
        instream->next_in_1  = (char *)packed->words; // last word
        instream->avail_in_1 = sizeof (uint64_t);
    }
}

static void comp_acgt_unpack (VBlockP vb, char *uncompressed_data, uint64_t uncompressed_len)
{
    BitArray *packed = buf_get_bitarray (&vb->compressed);
    packed->num_of_bits = uncompressed_len * 2;
    packed->num_of_words = roundup_bits2words64 (packed->num_of_bits);

    LTEN_bit_array (packed, true);

    bit_array_clear_excess_bits_in_top_word (packed);

    for (uint64_t i=0; i < uncompressed_len; i++) 
        uncompressed_data[i] = ACGT_DECODE(packed, i);
}

static inline void comp_non_acgt_transform (char *data, uint32_t len)
{
    for (uint32_t i=0; i < len; i++)
        data[i] ^= non_acgt_encode[(uint8_t)data[i]];
}

// returns true if successful and false if data_compressed_len is too small (but only if soft_fail is true)
// the purpose of NON-AGCT compression is to be used after ACGT, on the same data, to capture the characters that are NOT ACGT
// -- an A,C,G or T character is encoded as 0 (its encoded by CODEC_ACGT)
// -- an a,c,g or t character is encoded as 1 (its encoded by CODEC_ACGT the same as its uppercase counterpart
// -- any other character remains as is
// We modify the source data (and hence NON-AGCT compression is destructive!) and then recursively compress it with bz2.
bool comp_compress_non_acgt (VBlock *vb, Codec codec,
                             const char *uncompressed, uint32_t uncompressed_len, // option 1 - compress contiguous data
                             LocalGetLineCallback callback,                        // option 2 - compress data one line at a time
                             char *compressed, uint32_t *compressed_len /* in/out */, 
                             bool soft_fail)
{
    START_TIMER;

    // option 1 - compress contiguous data
    if (uncompressed) 
        comp_non_acgt_transform ((char*)uncompressed, uncompressed_len);
    
    // option 2 - compress data one line at a time
    else if (callback) {

        for (unsigned line_i=0; line_i < vb->lines.len; line_i++) {

            char *data_1=0, *data_2=0;
            uint32_t data_1_len=0, data_2_len=0;
            
            callback (vb, line_i, &data_1, &data_1_len, &data_2, &data_2_len);

            comp_non_acgt_transform (data_1, data_1_len);
            comp_non_acgt_transform (data_2, data_2_len);
        }
    }
    else 
        ABORT0 ("Error in comp_compress_non_acgt: neither src_data nor callback is provided");
    
    COPY_TIMER (vb->profile.comp_compress_non_acgt); // excluding bzlib
    
    // now do the compression on the non-agct data
    // note: we don't support soft-fail because the allocated amount (uncompressed_len/2) is plenty for our textual data,
    // and we can't allow re-calling of this routine as the xor will undo itself
    return comp_compress_bzlib (vb, CODEC_BZ2, uncompressed, uncompressed_len, callback, compressed, compressed_len, false);
}

static void comp_apply_non_acgt_on_top_of_acgt (char *acgt, const char *non_acgt, uint64_t len)
{
    for (uint64_t i=0; i < len; i++)
        
        // if we have a 1 - we convert the nucleotide to lower case
        if (non_acgt[i] == 1) acgt[i] += 32;

        // if its a non-0, non-1 - we copy verbatim (this is usually, but not necessarily, 'N')
        else if (non_acgt[i]) acgt[i] = non_acgt[i];
}

void comp_uncompress_non_acgt (VBlock *vb, 
                               const char *compressed, uint32_t compressed_len,
                               char *uncompressed_data, uint64_t uncompressed_len)
{
    // first - do bzip2 decoding into vb->compressed
    buf_alloc (vb, &vb->compressed, uncompressed_len, 1.5, "compressed", 0);
    comp_uncompress (vb, CODEC_BZ2, compressed, compressed_len, vb->compressed.data, uncompressed_len);

    // second, fix "compressed", already containing ACGT data, with the data from this NON_ACGT data
    comp_apply_non_acgt_on_top_of_acgt (uncompressed_data, vb->compressed.data, uncompressed_len);
}

void comp_uncompress_acgt (VBlock *vb, 
                           const char *compressed, uint32_t compressed_len,
                           char *uncompressed_data, uint64_t uncompressed_len)
{
    ISzAlloc alloc_stuff = { .Alloc = lzma_alloc, .Free = lzma_free, .vb = vb};
    ELzmaStatus status;

    SizeT compressed_len64 = (uint64_t)compressed_len - LZMA_PROPS_SIZE; // first 5 bytes in compressed stream are the encoding properties

    uint64_t bitarray_size = roundup_bits2bytes64 (uncompressed_len * 2); // 4 nucleotides per byte, rounded up to whole 64b words
    buf_alloc (vb, &vb->compressed, bitarray_size, 2, "compressed", 0);
    
    uint64_t expected_num_uncompressed_bytes = roundup_bits2bytes64 (uncompressed_len * 2); // 4 nucleotides per byte, rounded up to whole bytes
    uint64_t actual_num_uncompressed_bytes = expected_num_uncompressed_bytes;
    
    SRes ret = LzmaDecode ((uint8_t *)vb->compressed.data, &actual_num_uncompressed_bytes, 
                            (uint8_t *)compressed + LZMA_PROPS_SIZE, &compressed_len64, 
                            (uint8_t *)compressed, LZMA_PROPS_SIZE, 
                            LZMA_FINISH_END, &status, &alloc_stuff);

    ASSERT (ret == SZ_OK && status == LZMA_STATUS_FINISHED_WITH_MARK, 
            "Error: LzmaDecode failed: ret=%s status=%s", lzma_errstr (ret), lzma_status (status)); 

    ASSERT (expected_num_uncompressed_bytes == actual_num_uncompressed_bytes, "Error in comp_uncompress while decompressing ACGT: expected_num_uncompressed_bytes(%u) != actual_num_uncompressed_bytes(%u)",
            (uint32_t)expected_num_uncompressed_bytes, (uint32_t)actual_num_uncompressed_bytes);
            
    comp_acgt_unpack (vb, uncompressed_data, uncompressed_len);    
}

