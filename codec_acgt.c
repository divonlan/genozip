// ------------------------------------------------------------------
//   comp_agct.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"
#include "endianness.h"
#include "piz.h"

// -------------------------------------------------------------------------------------
// acgt stuff
// compress a sequence of A,C,G,T nucleotides - first squeeze into 2 bits and then LZMA.
// It's about 25X faster and slightly better compression ratio than LZMA
// -------------------------------------------------------------------------------------

// decoder of 2bit encoding of nucleotides
const char acgt_decode[4] = { 'A', 'C', 'G', 'T' };

// table to convert ASCII to ACGT encoding. A,C,G,T (lower and upper case) are encoded as 0,1,2,3 respectively, 
// and everything else (including N) is encoded as 0
const uint8_t acgt_encode[256] = { ['A']=0, ['C']=1, ['G']=2, ['T']=3, 
                                   ['a']=0, ['c']=1, ['g']=2, ['t']=3 }; // all others are 0

//--------------
// ZIP side
//--------------

void codec_acgt_comp_init (VBlock *vb)
{
        Context *nonref_ctx   = &vb->contexts[DTF(nonref)];
        nonref_ctx->lcodec    = CODEC_ACGT; // ACGT is better than LZMA and BSC
        nonref_ctx->ltype     = LT_SEQUENCE;

        Context *nonref_x_ctx = nonref_ctx + 1;
        nonref_x_ctx->lcodec  = CODEC_XCGT;
        nonref_x_ctx->ltype   = LT_UINT8;
}

// packing of an array A,G,C,T characters into a 2-bit BitArray, stored in vb->compressed. 
/// returns true if non-ACGT was encountered.
bool codec_acgt_pack (BitArray *packed, const char *data, uint64_t data_len)
{
    // increase bit array to accomodate data
    uint64_t next_bit    = packed->nbits;
    packed->nbits += data_len * 2;
    packed->nwords = roundup_bits2words64 (packed->nbits);
    
    bool has_non_agct    = false;

    // pack nucleotides - each character is packed into 2 bits
    for (uint64_t i=0 ; i < data_len ; i++) {
        has_non_agct = has_non_agct || !IS_NUCLEOTIDE (data[i]);
        
        uint8_t encoding = acgt_encode[(uint8_t)data[i]];
        
        bit_array_assign (packed, next_bit,   (encoding & 1)     );
        bit_array_assign (packed, next_bit+1, (encoding & 2) >> 1);
        next_bit += 2;
    }

    return has_non_agct;
}

// This function decompsoses SEQ data into two buffers:
// 1. A,G,C,T characters are packed into a 2-bit BitArray, placed in vb->compressed and then compressed with ACGT.sub_codec
// 2. NONREF_X.local is constructed to be the same length on the SEQ data, with each characer corresponding to a character in SEQ:
// -- an A,C,G or T character in SEQ is corresponds to a \0 in NONREF_X 
// -- an a,c,g or t character in SEQ is corresponds to a \1 in NONREF_X 
// -- any other character is copied from SEQ as is
// NONREF_X.local is later compressed in codec_xcgt_compress with XCGT.sub_codec
bool codec_acgt_compress (VBlock *vb, SectionHeader *header,
                          const char *uncompressed,    // option 1 - compress contiguous data
                          uint32_t *uncompressed_len,
                          LocalGetLineCB callback,     // option 2 - compress data one line at a time
                          char *compressed, uint32_t *compressed_len /* in/out */, 
                          bool soft_fail)
{
    // table to convert SEQ data to ACGT exceptions. The character is XORed with the entry in the table
    static const uint8_t acgt_exceptions[256] = { 
        ['A']='A',   ['C']='C',   ['G']='G',   ['T']='T',  // -->0 (XORed with self)
        ['a']='a'^1, ['c']='c'^1, ['g']='g'^1, ['t']='t'^1 // -->1 (XORed with self XOR 1)
    };                                                     // all others are XORed with 0 and hence remain unchanged
    
    START_TIMER;
    
    #define PACK(data,len) { if (len) vb->has_non_agct = codec_acgt_pack (packed, (data), (len)) || vb->has_non_agct; }

    Context *nonref_ctx   = &vb->contexts[DTF(nonref)];
    Context *nonref_x_ctx = nonref_ctx + 1;

    ASSERTE0 (!vb->compressed.len && !vb->compressed.param, "expecting vb->compressed to be free, but its not");

    // we will pack into vb->compressed
    buf_alloc (vb, &vb->compressed, roundup_bits2bytes64 (*uncompressed_len * 2), 1, "compress");
    BitArray *packed = buf_get_bitarray (&vb->compressed);

    // option 1 - pack contiguous data
    if (uncompressed) {
START_TIMER;        
        // overlay the NONREF.local to NONREF_X.local to avoid needing more memory, as NONREF.local is not needed after packing
        buf_set_overlayable (&nonref_ctx->local);
        buf_overlay (vb, &nonref_x_ctx->local, &nonref_ctx->local, "contexts->local");

        PACK (uncompressed, *uncompressed_len); // pack into vb->compressed

        // calculate the exception in-place in NONREF.local also overlayed to NONREF_X.local
        for (uint32_t i=0; i < *uncompressed_len; i++) \
            ((char*)uncompressed)[i] = (uint8_t)(uncompressed[i]) ^ acgt_exceptions[(uint8_t)(uncompressed[i])];
COPY_TIMER (tmp1);
    }

    // option 2 - callback to get each line
    else if (callback) {
START_TIMER;        

        buf_alloc (vb, &nonref_x_ctx->local, *uncompressed_len, CTX_GROWTH, "contexts->local");
        for (uint32_t line_i=0; line_i < vb->lines.len; line_i++) {

            char *data_1=0;
            uint32_t data_1_len=0;
            
            callback (vb, line_i, &data_1, &data_1_len, *uncompressed_len - nonref_x_ctx->local.len);

            PACK (data_1, data_1_len);

            for (uint32_t i=0; i < data_1_len; i++) 
                NEXTENT (uint8_t, nonref_x_ctx->local) = (uint8_t)(data_1[i]) ^ acgt_exceptions[(uint8_t)(data_1[i])];
        }
COPY_TIMER (tmp2);
    }
    else 
        ABORT0 ("Error in codec_acgt_compress_nonref: neither src_data nor callback is provided");

    // note: we store in Little Endian unlike the rest of the data that is in Big Endian, because LTEN keeps the nucleotides in their
    // original order, and improves compression ratio by about 2%
    LTEN_bit_array (packed);

    Codec sub_codec = codec_args[header->codec].sub_codec;
    CodecCompress *compress = codec_args[sub_codec].compress;
    uint32_t packed_uncompressed_len = packed->nwords * sizeof (word_t);

    PAUSE_TIMER; // sub-codec compresssors account for themselves
    if (flag.show_time) codec_show_time (vb, "Subcodec", vb->profile.next_subname, sub_codec);
    compress (vb, header, (char *)packed->words, &packed_uncompressed_len, NULL, compressed, compressed_len, false); // no soft fail
    RESUME_TIMER (compressor_actg);

    buf_free (&vb->compressed);

    // note: NONREF_X will be compressed after us in codec_xcgt_compress, as it is the subsequent context, and its local is now populated

    COPY_TIMER (compressor_actg); // don't include sub-codec compressor - it accounts for itself
    return true;
}

//--------------
// PIZ side
//--------------

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
void codec_xcgt_uncompress (VBlock *vb, Codec codec, uint8_t param,
                            const char *compressed, uint32_t compressed_len,
                            Buffer *uncompressed_buf, uint64_t uncompressed_len,
                            Codec sub_codec)
{
    // uncompress NONREF_X using CODEC_XCGT.sub_codec (passed to us as sub_codec)
    codec_args[sub_codec].uncompress (vb, sub_codec, param, compressed, compressed_len, uncompressed_buf, uncompressed_len, CODEC_NONE);

    const BitArray *acgt_packed = buf_get_bitarray (&vb->compressed); // data from NONREF context (2-bit per base)
    const char *acgt_x = FIRSTENT (const char, *uncompressed_buf); // data from NONREF_X context
    
    Context *nonref_ctx = &vb->contexts[DTF(nonref)];
    char *nonref = FIRSTENT (char, nonref_ctx->local); // note: local was allocated by caller ahead of comp_uncompress -> codec_acgt_uncompress of the NONREF context

    for (uint32_t i=0; i < uncompressed_len; i++) {
        if      (!acgt_x || acgt_x[i] == 0) *nonref++ = ACGT_DECODE(acgt_packed, i);      // case 0: use acgt as is - 'A', 'C', 'G' or 'T'
        else if (           acgt_x[i] == 1) *nonref++ = ACGT_DECODE(acgt_packed, i) + 32; // case 1: convert to lower case - 'a', 'c', 'g' or 't'
        else                                *nonref++ = acgt_x[i];                        // case non-0/1: use acgt_x (this is usually, but not necessarily, 'N')
    }

    buf_free (&vb->compressed)
}

// Explanation of uncompression of data compressed with the ACGT codec:
// - ACGT-compressed data is stored in two consecutive sections, NONREF which has CODEC_ACGT, and NONREF_X which has sub_codec2
// 1) NONREF contains a 2-bit representation of the bases: is is uncompressed by codec_acgt_uncompress into vb->compressed using sub_codec
// 2) NONREF_X is a character array of exceptions and is uncompressed into NONREF_X.local by codec_xcgt_uncompress
// 3) codec_xcgt_uncompress also combines vb->compressed with NONREF_X.local to recreate NONREF.local - an LT_SEQUENCE local buffer
void codec_acgt_uncompress (VBlock *vb, Codec codec, uint8_t param,
                            const char *compressed, uint32_t compressed_len,
                            Buffer *uncompressed_buf, uint64_t num_bases,
                            Codec sub_codec)
{
    ASSERTE0 (!vb->compressed.len && !vb->compressed.param, "expected vb->compressed to be free, but its not");

    uint64_t bitmap_num_bytes = roundup_bits2bytes64 (num_bases * 2); // 4 nucleotides per byte, rounded up to whole 64b words
    buf_alloc (vb, &vb->compressed, bitmap_num_bytes, 1, "compressed");    

    // uncompress bitmap using CODEC_ACGT.sub_codec (passed to us as sub_codec) into vb->compressed
    codec_args[sub_codec].uncompress (vb, sub_codec, param, compressed, compressed_len, &vb->compressed, bitmap_num_bytes, CODEC_NONE);

    // finalize bitmap structure
    BitArray *packed     = buf_get_bitarray (&vb->compressed);
    packed->nbits  = num_bases * 2;
    packed->nwords = roundup_bits2words64 (packed->nbits);

    LTEN_bit_array (packed);

    bit_array_clear_excess_bits_in_top_word (packed);
}

