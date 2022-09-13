// ------------------------------------------------------------------
//   comp_agct.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "codec.h"
#include "vblock.h"
#include "buffer.h"
#include "strings.h"
#include "endianness.h"
#include "piz.h"
#include "context.h"
#include "file.h"

// -------------------------------------------------------------------------------------
// acgt stuff
// compress a sequence of A,C,G,T nucleotides - first squeeze into 2 bits and then LZMA.
// It's about 25X faster and slightly better compression ratio than LZMA
// -------------------------------------------------------------------------------------

// table to convert ASCII to ACGT encoding. A,C,G,T (lower and upper case) are encoded as 0,1,2,3 respectively, 
// and everything else (including N) is encoded as 0
const uint8_t acgt_encode[256] = { ['A']=0, ['C']=1, ['G']=2, ['T']=3,  // all others are 0
                                   ['a']=0, ['c']=1, ['g']=2, ['t']=3, 
                                   
                                   // IUPAC codes are mapped to one of their bases: http://www.bioinformatics.org/sms/iupac.html
                                   // the base to which a IUPAC is mapped, is the lowest alphanetically of its participating bases, per 1.6.1-REF in the VCF specification: https://samtools.github.io/hts-specs/VCFv4.3.pdf
                                   // (expected by vcf_refalt_lift_snp)
                                   ['U']=3, ['R']=0, ['Y']=1, ['S']=1,
                                   ['W']=0, ['K']=2, ['M']=0, ['B']=1,
                                   ['D']=0, ['H']=0, ['V']=0, ['N']=0,

                                   ['u']=3, ['r']=0, ['y']=1, ['s']=1,
                                   ['w']=0, ['k']=2, ['m']=0, ['b']=1,
                                   ['d']=0, ['h']=0, ['v']=0, ['n']=0  }; 

// same as actg_encode, but produces the complement base. 
// Note for the IUPACs: eg B=C,G,T (the lowest of them)=> C=1 
//                   comp(B)=G,C,A (the lowest of them)=> A=0
const uint8_t acgt_encode_comp[256] = 
                                 { ['A']=3, ['C']=2, ['G']=1, ['T']=0,  // all others are 0
                                   ['a']=3, ['c']=2, ['g']=1, ['t']=0, 
                                   
                                   ['U']=0, ['R']=1, ['Y']=0, ['S']=1,
                                   ['W']=0, ['K']=0, ['M']=2, ['B']=0,
                                   ['D']=0, ['H']=0, ['V']=1, ['N']=0,

                                   ['u']=0, ['r']=1, ['y']=0, ['s']=1,
                                   ['w']=0, ['k']=0, ['m']=2, ['b']=0,
                                   ['d']=0, ['h']=0, ['v']=1, ['n']=0  }; 


//--------------
// ZIP side
//--------------

void codec_acgt_comp_init (VBlockP vb, Did nonref_did_i)
{
        Context *nonref_ctx     = CTX(nonref_did_i);
        nonref_ctx->lcodec      = CODEC_ACGT; // ACGT is better than LZMA and BSC
        nonref_ctx->ltype       = LT_SEQUENCE;

        Context *nonref_x_ctx   = nonref_ctx + 1;
        nonref_x_ctx->ltype     = LT_UINT8;
        nonref_x_ctx->local_dep = DEP_L1;     // NONREF_X.local is created with NONREF.local is compressed
        nonref_x_ctx->lcodec    = CODEC_XCGT; // prevent codec_assign_best from assigning it a different codec
}

// packing of an array A,C,G,T characters into a 2-bit Bits, stored in vb->scratch. 
static inline void codec_acgt_pack (BitsP packed, rom data, uint64_t data_len)
{
    // increase bit array to accomodate data
    uint64_t next_bit = packed->nbits;
    packed->nbits += data_len * 2;
    packed->nwords = roundup_bits2words64 (packed->nbits);
    
    // pack nucleotides - each character is packed into 2 bits
    for (uint64_t i=0 ; i < data_len ; i++, next_bit += 2)       
        bits_assign2 (packed, next_bit, acgt_encode[(uint8_t)data[i]]);
}

// This function decompsoses SEQ data into two buffers:
// 1. A,C,G,T characters are packed into a 2-bit Bits, placed in vb->scratch and then compressed with ACGT.sub_codec
// 2. NONREF_X.local is constructed to be the same length on the SEQ data, with each characer corresponding to a character in SEQ:
// -- an A,C,G or T character in SEQ is corresponds to a \0 in NONREF_X 
// -- an a,c,g or t character in SEQ is corresponds to a \1 in NONREF_X 
// -- any other character is copied from SEQ as is
// NONREF_X.local is later compressed as a normal context (codec=XCGT, subcodec=as assigned)
COMPRESS (codec_acgt_compress)
{
    // table to convert SEQ data to ACGT exceptions. The character is XORed with the entry in the table
    static const uint8_t acgt_exceptions[256] = { 
        ['A']='A',   ['C']='C',   ['G']='G',   ['T']='T',  // -->0 (XORed with self)
        ['a']='a'^1, ['c']='c'^1, ['g']='g'^1, ['t']='t'^1 // -->1 (XORed with self XOR 1)
    };                                                     // all others are XORed with 0 and hence remain unchanged
    
    START_TIMER;
    
    #define PACK(data,len) { if (len) codec_acgt_pack (packed, (data), (len)); }

    Context *nonref_ctx   = ctx;
    Context *nonref_x_ctx = nonref_ctx + 1;
    Bits *packed;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (nonref_x_ctx->local.len) {
        packed = (BitsP)&vb->scratch;
        goto compress_sub;
    }
    
    ASSERTNOTINUSE(vb->scratch);

    // we will pack into vb->scratch
    buf_alloc (vb, &vb->scratch, 0, roundup_bits2bytes64 (*uncompressed_len * 2), uint8_t, 1, "scratch");
    packed = (BitsP)&vb->scratch;

    // option 1 - pack contiguous data
    if (uncompressed) {
        // overlay the NONREF.local to NONREF_X.local to avoid needing more memory, as NONREF.local is not needed after packing
        buf_set_overlayable (&nonref_ctx->local);
        buf_overlay (vb, &nonref_x_ctx->local, &nonref_ctx->local, "contexts->local");

        PACK (uncompressed, *uncompressed_len); // pack into vb->scratch

        // calculate the exception in-place in NONREF.local also overlayed to NONREF_X.local
        for (uint32_t i=0; i < *uncompressed_len; i++) \
            ((uint8_t*)uncompressed)[i] = (uint8_t)(uncompressed[i]) ^ acgt_exceptions[(uint8_t)(uncompressed[i])];
    }

    // option 2 - callback to get each line
    else if (get_line_cb) {
        buf_alloc (vb, &nonref_x_ctx->local, 0, *uncompressed_len, uint8_t, CTX_GROWTH, "contexts->local");
        for (uint32_t line_i=0; line_i < vb->lines.len32; line_i++) {

            char *data_1=0;
            uint32_t data_1_len=0;
            
            get_line_cb (vb, ctx, line_i, &data_1, &data_1_len, *uncompressed_len - nonref_x_ctx->local.len32, NULL);

            PACK (data_1, data_1_len);

            for (uint32_t i=0; i < data_1_len; i++) 
                BNXT8 (nonref_x_ctx->local) = (uint8_t)(data_1[i]) ^ acgt_exceptions[(uint8_t)(data_1[i])];
        }
    }
    else 
        ABORT ("%s: \"%s\": neither src_data nor callback is provided", VB_NAME, name);

    // get codec for NONREF_X header->lcodec remains CODEC_XCGT, and we set subcodec to the codec discovered in assign, and set to nonref_ctx->lcode
    Codec z_lcodec = ZCTX(nonref_x_ctx->did_i)->lcodec;
    nonref_x_ctx->lcodec = z_lcodec; // possibly set by a previous VB call to codec_assign_best_codec
    nonref_x_ctx->lsubcodec_piz = codec_assign_best_codec (vb, nonref_x_ctx, NULL, SEC_LOCAL);
    if (nonref_x_ctx->lsubcodec_piz == CODEC_UNKNOWN) nonref_x_ctx->lsubcodec_piz = CODEC_NONE; // really small
    
    nonref_x_ctx->lcodec = CODEC_XCGT;

    // note: we store in Little Endian unlike the rest of the data that is in Big Endian, because LTEN keeps the nucleotides in their
    // original order, and improves compression ratio by about 2%
    LTEN_bits (packed);

    // get codec for NONREF - header->lcodec remains CODEC_ACGT, and we set subcodec to the codec discovered in assign, and set to nonref_ctx->lcodec 
    z_lcodec = ZCTX(nonref_ctx->did_i)->lcodec;
    nonref_ctx->lcodec = z_lcodec; // possibly set by a previous VB call to codec_assign_best_codec
    header->sub_codec = codec_assign_best_codec (vb, nonref_ctx, &vb->scratch, SEC_LOCAL);
    if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small
    
    compress_sub: {
        CodecCompress *compress = codec_args[header->sub_codec].compress;
        uint32_t packed_uncompressed_len = packed->nwords * sizeof (uint64_t);

        if (flag.show_time) codec_show_time (vb, "Subcodec", vb->profile.next_subname, header->sub_codec);

        PAUSE_TIMER(vb); // sub-codec compresssors account for themselves
        if (!compress (vb, ctx, header, (char *)packed->words, &packed_uncompressed_len, NULL, compressed, compressed_len, soft_fail, name)) return false;
        RESUME_TIMER (vb, compressor_actg);
    }

    buf_free (vb->scratch);
    // note: NONREF_X will be compressed after us as it is the subsequent context, and its local is now populated

    COPY_TIMER_COMPRESS (compressor_actg); // don't include sub-codec compressor - it accounts for itself
    return true;
}

//--------------
// PIZ side
//--------------

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
UNCOMPRESS (codec_xcgt_uncompress)
{
    // uncompress NONREF_X using CODEC_XCGT.sub_codec (passed to us as sub_codec)
    codec_args[sub_codec].uncompress (vb, sub_codec, param, STRa(compressed), uncompressed_buf, uncompressed_len, CODEC_NONE, name);

    const Bits *acgt_packed = (BitsP)&vb->scratch; // data from NONREF context (2-bit per base)
    rom acgt_x = B1ST (const char, *uncompressed_buf); // data from NONREF_X context
    
    Context *nonref_ctx = CTX(DTF(nonref));
    ASSERTISALLOCED (nonref_ctx->local);
    
    char *nonref = B1STc (nonref_ctx->local); // note: local was allocated by caller ahead of comp_uncompress -> codec_acgt_uncompress of the NONREF context

    for (uint32_t i=0; i < uncompressed_len; i++) {
        if      (!acgt_x || acgt_x[i] == 0) *nonref++ = base_by_idx(acgt_packed, i);      // case 0: use acgt as is - 'A', 'C', 'G' or 'T'
        else if (           acgt_x[i] == 1) *nonref++ = base_by_idx(acgt_packed, i) + 32; // case 1: convert to lower case - 'a', 'c', 'g' or 't'
        else                                *nonref++ = acgt_x[i];                        // case non-0/1: use acgt_x (this is usually, but not necessarily, 'N')
    }

    buf_free (vb->scratch);
}

// Explanation of uncompression of data compressed with the ACGT codec:
// - ACGT-compressed data is stored in two consecutive sections, NONREF which has CODEC_ACGT, and NONREF_X which has sub_codec2
// 1) NONREF contains a 2-bit representation of the bases: is is uncompressed by codec_acgt_uncompress into vb->scratch using sub_codec
// 2) NONREF_X is a character array of exceptions and is uncompressed into NONREF_X.local by codec_xcgt_uncompress
// 3) codec_xcgt_uncompress also combines vb->scratch with NONREF_X.local to recreate NONREF.local - an LT_SEQUENCE local buffer
UNCOMPRESS (codec_acgt_uncompress)
{
    ASSERTNOTINUSE (vb->scratch);

    uint64_t bitmap_num_bytes = roundup_bits2bytes64 (uncompressed_len * 2); // 4 nucleotides per byte, rounded up to whole 64b words
    buf_alloc (vb, &vb->scratch, 0, bitmap_num_bytes, char, 1, "scratch");    

    // uncompress bitmap using CODEC_ACGT.sub_codec (passed to us as sub_codec) into vb->scratch
    codec_args[sub_codec].uncompress (vb, sub_codec, param, compressed, compressed_len, &vb->scratch, bitmap_num_bytes, CODEC_NONE, name);

    // finalize bitmap structure
    Bits *packed     = (BitsP)&vb->scratch;
    packed->nbits  = uncompressed_len * 2;
    packed->nwords = roundup_bits2words64 (packed->nbits);

    LTEN_bits (packed);

    bits_clear_excess_bits_in_top_word (packed);
}

