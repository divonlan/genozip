// ------------------------------------------------------------------
//   comp_agct.c
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "codec.h"
#include "piz.h"
#include "seg.h"

// -------------------------------------------------------------------------------------
// acgt stuff
// compress a sequence of A,C,G,T nucleotides - first squeeze into 2 bits and then LZMA.
// It's about 25X faster and slightly better compression ratio than LZMA
// -------------------------------------------------------------------------------------


//--------------
// ZIP side
//--------------

void codec_acgt_seg_initialize (VBlockP vb, Did nonref_did_i,
                                bool has_x) // caller declares that sequences contain strictly only A,C,G,T (uppercase) - verified by caller
{
    ContextP nonref_ctx     = CTX(nonref_did_i);
    nonref_ctx->lcodec      = CODEC_ACGT; // ACGT is better than LZMA and BSC for "random" sequences
    nonref_ctx->ltype       = LT_BLOB;
    nonref_ctx->no_stons    = true;       // we're storing the sequencing in local, so we can't also have singletons
    nonref_ctx->flags.acgt_no_x = !has_x;

    ASSERT (ZCTX(nonref_did_i)->lcodec_hard_coded && ZCTX(nonref_did_i)->lcodec == CODEC_ACGT,
            "%s: ACGT for %s not set up in segconf finalize", VB_NAME, nonref_ctx->tag_name);

    if (has_x) {
        ContextP nonref_x_ctx   = nonref_ctx + 1;
        nonref_x_ctx->ltype     = LT_SUPP;
        nonref_x_ctx->local_dep = DEP_L1;     // NONREF_X.local is created with NONREF.local is compressed
        nonref_x_ctx->lcodec    = CODEC_XCGT; // prevent codec_assign_best from assigning it a different codec
    }
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
        bits_assign2 (packed, next_bit, acgt_encode(data[i]));
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
    alignas(64) static const uint8_t acgt_exceptions[256] = { 
        ['A']='A',   ['C']='C',   ['G']='G',   ['T']='T',  // -->0 (XORed with self)
        ['a']='a'^1, ['c']='c'^1, ['g']='g'^1, ['t']='t'^1 // -->1 (XORed with self XOR 1)
    };                                                     // all others are XORed with 0 and hence remain unchanged
    
    START_TIMER;
    
    #define PACK(data,len) { if (len) codec_acgt_pack (packed, (data), (len)); }

    ContextP nonref_ctx = ctx;
    bool has_x = !nonref_ctx->flags.acgt_no_x;

    ContextP nonref_x_ctx = has_x ? (nonref_ctx + 1) : NULL;
    BitsP packed;
    
    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (has_x && nonref_x_ctx->local.len) {
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
        if (has_x) {
            buf_set_shared (&nonref_ctx->local);
            buf_overlay (vb, &nonref_x_ctx->local, &nonref_ctx->local, C_LOCAL);
        }

        PACK (uncompressed, *uncompressed_len); // pack into vb->scratch

        // calculate the exception in-place in NONREF.local also overlayed to NONREF_X.local
        if (has_x) 
            for (uint8_t *next=(uint8_t *)uncompressed, *after=next + *uncompressed_len; next < after; next++)
                *next = *next ^ acgt_exceptions[*next];
    }

    // option 2 - callback to get each line 
    else if (get_line_cb) {
        ASSERT0 (has_x, "ACGT compression with get_line_cb is only support with has_X"); // we can easily add support if needed in the future
        
        buf_alloc (vb, &nonref_x_ctx->local, 0, *uncompressed_len, uint8_t, CTX_GROWTH, C_LOCAL);
        for_line {
            STRw0𐤐(data_1);
            get_line_cb (vb, ctx, line_i, pSTRa(data_1), *uncompressed_len - nonref_x_ctx->local.len32, NULL);

            PACK (data_1, data_1_len);

            ASSERT (nonref_x_ctx->local.len + data_1_len <= nonref_x_ctx->local.size, "nonref_x_ctx overflow: data_1_len=%u local=%.*s", 
                    data_1_len, (int)sizeof(BufDescType)-1, buf_desc (&nonref_x_ctx->local).s);
            
            for (uint8_t *restrict next=BAFT8(nonref_x_ctx->local), *after=next + data_1_len; next < after; next++, data_1++) 
                *next = (uint8_t)*data_1 ^ acgt_exceptions[(uint8_t)*data_1];

            nonref_x_ctx->local.len32 += data_1_len;
        }
    }
    else 
        ABORT ("%s: \"%s\": neither src_data nor callback is provided", VB_NAME, name);

    bits_clear_excess_bits_in_top_word (packed, false); // for good measure (V15)

    // case: no exception basess after all
    if (buf_is_zero (&nonref_x_ctx->local)) {
        has_x = false;  
        header->flags.ctx.acgt_no_x = true; 
        buf_destroy (nonref_x_ctx->local); // cannot use buf_free for an overlaid buffer 
    }

    // get codec for NONREF_X header->lcodec remains CODEC_XCGT, and we set subcodec to the codec discovered in assign, and set to nonref_ctx->lcode
    Codec z_lcodec;
    if (has_x) {
        z_lcodec = ZCTX(nonref_x_ctx->did_i)->lcodec;
        nonref_x_ctx->lcodec = z_lcodec; // possibly set by a previous VB call to codec_assign_best_codec
        PAUSE_TIMER(vb); // codec_assign_best_codec account for itself
        nonref_x_ctx->lsubcodec_piz = codec_assign_best_codec (vb, nonref_x_ctx, NULL, SEC_LOCAL);
        RESUME_TIMER (vb, compressor_acgt);
        if (nonref_x_ctx->lsubcodec_piz == CODEC_UNKNOWN) nonref_x_ctx->lsubcodec_piz = CODEC_NONE; // really small
        
        nonref_x_ctx->lcodec = CODEC_XCGT;
    }

    // note: we store in Little Endian unlike the rest of the data that is in Big Endian, because LTEN keeps the nucleotides in their
    // original order, and improves compression ratio by about 2%
    LTEN_bits (packed);

    nonref_ctx->lcodec = header->sub_codec = (vb->scratch.len32 * sizeof (uint64_t) >= MIN_LEN_FOR_COMPRESSION) ? CODEC_LZMA : CODEC_NONE;
    
    compress_sub: {
        CodecCompress *compress = codec_args[header->sub_codec].compress;
        uint32_t packed_uncompressed_len = packed->nwords * sizeof (uint64_t);

        if (flag.show_time) codec_show_time (vb, "Subcodec", vb->profile.next_subname, header->sub_codec);

        PAUSE_TIMER(vb); // sub-codec compresssors account for themselves
        if (!compress (vb, ctx, header, (char *)packed->words, &packed_uncompressed_len, NULL, compressed, compressed_len, soft_fail, name)) return false;
        RESUME_TIMER (vb, compressor_acgt);
    }

    buf_free (vb->scratch);
    // note: NONREF_X will be compressed after us as it is the subsequent context, and its local is now populated

    COPY_TIMER_COMPRESS_BY_CODEC (compressor_acgt); // don't include sub-codec compressor - it accounts for itself
    return true;
}

//--------------
// PIZ side
//--------------

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
UNCOMPRESS (codec_xcgt_uncompress)
{
    ContextP nonref_ctx = ctx - 1;
    ASSERTISALLOCED (nonref_ctx->local);

    // uncompress NONREF_X using CODEC_XCGT.sub_codec (passed to us as sub_codec)
    codec_args[sub_codec].uncompress (vb, ctx, sub_codec, param, STRa(compressed), uncompressed_buf, uncompressed_len, CODEC_NONE, name);

    START_TIMER;
    ConstBitsP packed = (BitsP)&nonref_ctx->packed;    // data from NONREF context (2-bit per base)
    rom acgt_x = B1ST (const char, *uncompressed_buf); // data from NONREF_X context (possibly NULL)
        
    char *nonref = B1STc (nonref_ctx->local); // note: local was allocated by caller ahead of comp_uncompress -> codec_acgt_uncompress of the NONREF context

    decl_acgt_decode;
    // note: in case of !acgt_x branch predictor would be right always, so very fast
    for (uint32_t i=0; i < uncompressed_len; i++) {
        if      (!acgt_x || acgt_x[i] == 0) *nonref++ = base_by_idx(packed, i);      // case 0: use acgt as is - 'A', 'C', 'G' or 'T'
        else if (           acgt_x[i] == 1) *nonref++ = base_by_idx(packed, i) + 32; // case 1: convert to lower case - 'a', 'c', 'g' or 't'
        else                                *nonref++ = acgt_x[i];                   // case non-0/1: use acgt_x (this is usually, but not necessarily, 'N')
    }

    buf_free (vb->scratch);
    COPY_TIMER (compressor_xcgt);
}

// Explanation of uncompression of data compressed with the ACGT codec:
// - ACGT-compressed data is stored in two consecutive sections, NONREF which has CODEC_ACGT, and NONREF_X which has sub_codec2
// 1) NONREF contains a 2-bit representation of the bases: is is uncompressed by codec_acgt_uncompress into vb->scratch using sub_codec
// 2) NONREF_X is a character array of exceptions and is uncompressed into NONREF_X.local by codec_xcgt_uncompress
// 3) codec_xcgt_uncompress also combines vb->scratch with NONREF_X.local to recreate NONREF.local - an LT_BLOB local buffer
UNCOMPRESS (codec_acgt_uncompress)
{
    BufferP packed_buf = ctx->flags.acgt_no_x ? &vb->scratch : &ctx->packed;
    ASSERTNOTINUSE (*packed_buf);

    uint64_t bitmap_num_bytes = roundup_bits2bytes64 (uncompressed_len * 2); // 4 nucleotides per byte, rounded up to whole 64b words
    buf_alloc (vb, packed_buf, 0, bitmap_num_bytes, char, 1, "packed");    

    // uncompress bitmap using CODEC_ACGT.sub_codec (passed to us as sub_codec) into vb->scratch
    codec_args[sub_codec].uncompress (vb, ctx, sub_codec, param, compressed, compressed_len, packed_buf, bitmap_num_bytes, CODEC_NONE, name);

    // finalize bitmap structure
    START_TIMER;
    BitsP packed   = (BitsP)packed_buf;
    packed->nbits  = uncompressed_len * 2;
    packed->nwords = roundup_bits2words64 (packed->nbits);

    LTEN_bits (packed);

    bits_clear_excess_bits_in_top_word (packed, false);

    // decode here if no X. If there's X we decode in codec_xcgt_uncompress (acgt_no_x added in 15.0.13)
    if (ctx->flags.acgt_no_x) {
        char *nonref = B1STc (ctx->local); // note: local was allocated by caller ahead of comp_uncompress -> codec_acgt_uncompress of the NONREF context

        decl_acgt_decode;
        for (uint32_t i=0; i < uncompressed_len; i++) 
            *nonref++ = base_by_idx(packed, i);

        buf_free (*packed_buf);
    }

    COPY_TIMER (compressor_acgt);
}

