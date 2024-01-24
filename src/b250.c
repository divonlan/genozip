// ------------------------------------------------------------------
//   b250.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "b250.h"
#include "context.h"
#include "codec.h"

// single-length encoding (up to 15.0.37)
// Format on data in Context.b250: Each entry is either a single-byte special-code value 0xFA-0xFF, OR a 1, 2 or 4 big-endian integer.
// The number of bytes is determined by Context.b250_size transmitted via SectionHeaderCtx.b250_size, its selection is done separately for each VB.
// In FASTQ, for b250-pairing, b250 of the pair is stored in Context.pair and its size in Context.pair_b250_size
#define BASE250_EMPTY_SF   0xFA // 250 empty string
#define BASE250_MISSING_SF 0xFB // 251 container item missing, remove preceding separator
#define BASE250_ONE_UP     0xFC // 252 value is one higher than previous value. 
#define BASE250_MOST_FREQ0 0xFD // 253 this translates to 0,1,2 representing the most frequent values (according to vb_i=1 sorting).
#define BASE250_MOST_FREQ1 0xFE // 254
#define BASE250_MOST_FREQ2 0xFF // 255

#define B250_MAX_WI_1BYTE  (0xFA         - 1) // maximal value in Context.b250 of non-special-code entries, for b250_size to be set to 1 byte
#define B250_MAX_WI_2BYTES ((0xFA << 8)  - 1) // same, for 2 bytes
#define B250_MAX_WI_3BYTES ((0xFA << 16) - 1) // same, for 3 bytes

// var-length encording (since 15.0.38)
#define VARL_MIN_1B 0
#define VARL_MAX_1B 126 // 127 reserved for WORD_INDEX_ONE_UP
#define VARL_MIN_2B (VARL_MAX_1B + 1)
#define VARL_MAX_2B (VARL_MIN_2B + (1<<14) - 1 - 2) // 3ffe, 3fff reserved for WORD_INDEX_EMPTY/MISSING 
#define VARL_MIN_3B (VARL_MAX_2B + 1)
#define VARL_MAX_3B (VARL_MIN_3B + (1<<21) - 1)
#define VARL_MIN_4B (VARL_MAX_3B + 1)  
#define VARL_MAX_4B ((1<<29) - 1)  // // 4B can represent any word_index from 0 - node_indices before conversion use this regardless of value - so its range spans 0 to 2^29-1

#define VARL_EMPTY   0xBFFE // 2B representation of WORD_INDEX_EMPTY
#define VARL_MISSING 0xBFFF // 2B representation of WORD_INDEX_MISSING
#define VARL_ONE_UP  127

#define VARL_BYTES(msb) ((((msb)>>7)==0)?1 : (((msb)>>6)==0b10)?2 : (((msb)>>5)==0b110)?3 : 4)

// ------
// ZIP
// ------

// caller guarantees that n∈[1-4]
static inline void memcpy4 (uint8_t *dst, uint8_t *src, uint_fast8_t n)
{
    if      (n == 1) { dst[0]=src[0]; }
    else if (n == 2) { dst[0]=src[0]; dst[1]=src[1]; }
    else if (n == 3) { dst[0]=src[0]; dst[1]=src[1]; dst[2]=src[2]; }
    else             { dst[0]=src[0]; dst[1]=src[1]; dst[2]=src[2]; dst[3]=src[3]; }
}

// SEG (before conversion): MSB containing type is last byte
WordIndex b250_seg_get_wi (const uint8_t *msb_p) // pointer to msb - last byte of the word
{
    uint8_t msb = *msb_p;

    if ((msb >> 7) == 0)
        return msb;

    else if ((msb >> 6) == 0b10) {
        uint_fast16_t word = GET_UINT16 (msb_p-1);
        return word == VARL_EMPTY   ? WORD_INDEX_EMPTY
             : word == VARL_MISSING ? WORD_INDEX_MISSING
             :                        (WordIndex)(word & 0x3fff) + VARL_MIN_2B;  
    }

    else if ((msb >> 5) == 0b110)
        return (GET_UINT24 (msb_p-2) & 0x1fffff) + VARL_MIN_3B;

    else 
        return (GET_UINT32 (msb_p-3) & 0x1fffffff);
}

// returns number of bytes written
static inline uint32_t b250_set_wi (uint8_t *dst, // begining of writing if piz_format=false, end of writing if true 
                                    WordIndex wi, bool piz_format)
{
    uint32_t enc; // encoding
    uint_fast8_t enc_len;

    switch (wi) {
        case VARL_MIN_1B ... VARL_MAX_1B : enc = wi;                                   enc_len = 1; break; // most significant bit = 0
        case VARL_MIN_2B ... VARL_MAX_2B : enc = (0b10UL  << 14) | (wi - VARL_MIN_2B); enc_len = 2; break; // most significant 2 bits = 10
        case VARL_MIN_3B ... VARL_MAX_3B : enc = (0b110UL << 21) | (wi - VARL_MIN_3B); enc_len = 3; break; // most significant 3 bits = 110
        case VARL_MIN_4B ... VARL_MAX_4B : enc = (0b111UL << 29) | wi;                 enc_len = 4; break; // most significant 3 bits = 111. note that value is relative to 0, not VARL_MIN_4B
        case WORD_INDEX_ONE_UP           : enc = VARL_ONE_UP;                          enc_len = 1; break; 
        case WORD_INDEX_EMPTY            : enc = VARL_EMPTY;                           enc_len = 2; break; 
        case WORD_INDEX_MISSING          : enc = VARL_MISSING;                         enc_len = 2; break; 

        default                          : ABORT ("wi=%d out of range [-4..-2,0..%u]", wi, VARL_MAX_4B);
    }

    // convert to piz format: big endian, so that MSB, carrying the type, is the first byte.
    // example for enc_len=2: [n=0x00004411 bytes=(11,44,00,00)] ; after BGEN=[n=0x11440000, bytes=(00,00,44,11)] ; after shift=[n=0x00001144, bytes=(44,11,00,00)] -> encode 44,11.
    if (piz_format) {
        enc = BGEN32 (enc) >> ((4 - enc_len) * 8); 
        dst -= enc_len - 1;
    }

    memcpy4 (dst, (uint8_t *)&enc, enc_len);

    return enc_len;
}

void b250_seg_append (VBlockP vb, ContextP ctx, WordIndex node_index)
{
    #define AT_LEAST(did_i) ((uint64_t)(10.0 + (((did_i) < MAX_NUM_PREDEFINED) ? segconf.b250_per_line[did_i] * (float)(vb->lines.len32) : 0)))

    // case: context is currently all_the_same - (ie count>=1, but only one actual entry)...
    if (ctx->flags.all_the_same) {
        
        // case: this node_index causes the context to no longer be all-the-same
        if (b250_seg_get_last (ctx) != node_index) { 
            int word_len = ctx->b250.len32;

            buf_alloc (vb, &ctx->b250, word_len * (ctx->b250.count+1), AT_LEAST(ctx->did_i), char, CTX_GROWTH, CTX_TAG_B250); // add 1 more, meaning a total of len+1
        
            // add (count-1) copies of first_node to the one already existing
            for (unsigned i=0; i < ctx->b250.count - 1; i++) {
                memcpy4 (BAFT8(ctx->b250), B1ST8(ctx->b250), word_len);
                ctx->b250.len32 += word_len;
            }

            ctx->flags.all_the_same = false; // no longer all_the_same 

            goto append_do;
        }
    }

    // b250 is not all_the_same: either it is empty, or it contains multiple different node_index values
    else append_do: {
        // case: this is the first entry - its all_the_same in a trivial way (still, we mark it so b250_zip_generate_section can potentially eliminate it)
        if (!ctx->b250.count) ctx->flags.all_the_same = true;

        buf_alloc (vb, &ctx->b250, 4, ctx->flags.all_the_same ? 0 : AT_LEAST(ctx->did_i), char, CTX_GROWTH, CTX_TAG_B250);

        // Logic: during seg, we encode variable length integers with the last byte (the MSB in little endian)
        // indicating the type. This is we can identify the type in b250_seg_remove_last. Then, in b250_zip_generate_one,
        // which switch the order of the bytes, to allow forward traversal during piz.

        // new nodes always get 4B b250 - it might be shortened later, after conversion
        if (node_index >= (WordIndex)ctx->ol_nodes.len32) {
            PUT_UINT32 (BAFT8(ctx->b250), (0b111UL << 29) | node_index); // MSb 111 = 32bit. Note that this is last byte of the Little Endian representation
            ctx->b250.len32 += 4;
        }

        else
            ctx->b250.len32 += b250_set_wi (BAFT8(ctx->b250), node_index, false);
    }

    ctx->b250.count++; // counts number of words in this b250    
}

void b250_seg_remove_last (VBlockP vb, ContextP ctx, WordIndex node_index/*optional*/)
{
    ctx_decrement_count (vb, ctx, node_index != WORD_INDEX_NONE ? node_index : b250_seg_get_last (ctx));
    
    if (!ctx->flags.all_the_same || ctx->b250.count == 1) {
        int remove_bytes = VARL_BYTES (*BLST8(ctx->b250));
        ASSERT (ctx->b250.len32 >= remove_bytes, "Cannot remove %u bytes from %s.b250 because its len=%u", remove_bytes, ctx->tag_name, ctx->b250.len32);

        ctx->b250.len32 -= remove_bytes;
    }

    // update all_the_same - true, if len is down to 1, and false if it is down to 0
    ctx->b250.count--;
    if (ctx->b250.count <= 1) ctx->flags.all_the_same = (bool)ctx->b250.count;
}

static inline uint_fast8_t get_converted_wi (VBlockP vb, ContextP ctx, const uint8_t *msb_p, WordIndex needs_conversion_threadshold, WordIndex *converted_wi)
{
    uint_fast8_t orig_wi_len = VARL_BYTES (*msb_p);
    WordIndex orig_wi = b250_seg_get_wi (msb_p);

    if (orig_wi >= needs_conversion_threadshold) 
        *converted_wi = node_index_to_word_index (vb, ctx, orig_wi); 

    else 
        *converted_wi = orig_wi;

    return orig_wi_len;
}

// we convert the b250 data to PIZ VARL format, by making the following modifications:
// 1. change words to big endian, so that the MSB, used to determine the word length, is the first byte
// 2. every node_index new to this VB (which is always 4B) is converted to a word_index of the appropriate length
// 3. in case a word_index is +1 the previous word index, it is replaced with VARL_ONE_UP
bool b250_zip_generate (VBlockP vb, ContextP ctx)
{
    START_TIMER;
    bool ret = true;

    ASSERT (ctx->dict_id.num, "tag_name=%s did_i=%u: ctx->dict_id=0 despite ctx->b250 containing data", ctx->tag_name, (unsigned)(ctx - vb->contexts));
    ASSERT (ctx->nodes_converted, "expecting nodes of %s to be converted", ctx->tag_name); // still index/len

    bool show = flag.show_b250 || dict_id_typeless (ctx->dict_id).num == flag.dict_id_show_one_b250.num;
    if (show) bufprintf (vb, &vb->show_b250_buf, "%s %s: ", VB_NAME, ctx->tag_name);

    // case: all-the-same b250 survived dropping (in ctx_drop_all_the_same) - we just shorten it to one entry
    if (ctx->flags.all_the_same && ctx->b250.count > 1) { 
        if (flag.debug_generate) 
            iprintf ("%s: %s[%u].b250 is \"all_the_same\" - shortened b250 from count=%"PRIu64" to 1\n", 
                     VB_NAME, ctx->tag_name, ctx->did_i, ctx->b250.count);
        
        ctx->b250.count = 1;
    }
    else
        if (flag.debug_generate) 
            iprintf ("%s: %s[%u].b250 count=%"PRIu64" all_the_same=%s no_stons=%s\n", 
                     VB_NAME, ctx->tag_name, ctx->did_i, ctx->b250.count, TF(ctx->flags.all_the_same), TF(ctx->no_stons));

    // determine size of word_index elements
    ctx->b250_size = B250_VARL; // since 15.0.38

    // we modify in-place. Note that the converted b250 will be of length <= the original, this is
    // because words might be shortened (due to ONE_UP and ni->wi) but never lengthened
    const uint8_t *first = B1ST8(ctx->b250);
    const uint8_t *src   = BLST8(ctx->b250); // point to the *last* byte of the next original wi to be read
    uint8_t *dst         = BLST8(ctx->b250); // point to the *last* byte of the next converted wi to be written

    WordIndex converted_wi, prev_converted_wi = WORD_INDEX_NONE;
    uint_fast8_t orig_wi_len, prev_orig_wi_len=0;
    WordIndex needs_conversion_threadshold = ctx->ol_nodes.len32;
    
    // scan backwards as type is in MSB which is the last byte in each yet-to-be-converted b250
    while (src >= first) {
        if (prev_converted_wi != WORD_INDEX_NONE) {
            converted_wi     = prev_converted_wi;
            orig_wi_len      = prev_orig_wi_len;
        }
        else  // first iteration
            orig_wi_len = get_converted_wi (vb, ctx, src, needs_conversion_threadshold, &converted_wi);

        // get previous b250 (that will be converted in the next iteration, since we are scanning backwards)
        if (src - orig_wi_len >= first) 
            prev_orig_wi_len = get_converted_wi (vb, ctx, src - orig_wi_len, needs_conversion_threadshold, &prev_converted_wi);

        if (prev_converted_wi >= 0 && converted_wi >= 0 && converted_wi == prev_converted_wi + 1) 
            converted_wi = WORD_INDEX_ONE_UP;

        src -= orig_wi_len;
        dst -= b250_set_wi (dst, converted_wi, true);
    } 

    ASSERT (src + 1 == first, "%s: src in backward scan exceeded start of %s.b250 array", VB_NAME, ctx->tag_name);
    ASSERT (dst + 1 >= first, "%s: dst in backward scan exceeded start of %s.b250 array", VB_NAME, ctx->tag_name);

    // now, the b250 buffer is possibly shortened from its start, adjust buffer fields
    uint32_t shortened_by = (dst + 1) - first;

    // shift start of buffer in memory. will be reset in buf_free. similar to buffer partial overlay logic.
    ctx->b250.data  += shortened_by;
    ctx->b250.len32 -= shortened_by;
    ctx->b250.size  -= shortened_by;

    // in case we are using "pair identical", drop this section if it is an R2 section identical to its R1 counterpart
    if (is_fastq_pair_2 (vb) && fastq_zip_use_pair_identical (ctx->dict_id) && buf_issame (&ctx->b250, &ctx->b250R1, 1)) {
        ctx->b250.len = 0; 
        
        if (flag.debug_generate) iprintf ("%s: %s[%u].b250 dropped because it is an R2 section which is identical to its R1 counterpart\n", VB_NAME, ctx->tag_name, ctx->did_i);

        ret = false;
    }
    
    // xxx - print b250 if show

    if (show && ret) {
        bufprintf (vb, &vb->show_b250_buf, "%s", "\n");
        iprintf ("%.*s", STRfb(vb->show_b250_buf));
    }

    if (show) buf_free (vb->show_b250_buf);
    COPY_TIMER (b250_zip_generate); // codec_assign measures its own time

    if (ret) codec_assign_best_codec (vb, ctx, NULL, SEC_B250);

    return ret;
}

// ------
// PIZ
// ------

WordIndex b250_piz_decode (bytes *b, bool advance, B250Size b250_size, rom ctx_name)
{
    ASSERT (*b, "*b is NULL in ctx=%s", ctx_name);

    #define RETURN(res,n) ({ if (advance) { *b += (n); } return (WordIndex)(res); })

    // case: files starting 15.0.38
    if (b250_size == B250_VARL) { 
        uint8_t msb = (*b)[0];
        
        if ((msb >> 7) == 0) // 1 byte (7 bit)
            RETURN (((msb == VARL_ONE_UP) ? WORD_INDEX_ONE_UP : msb), 1);

        else if ((msb >> 6) == 0b10) { // 2 bytes (14 bits)
            uint16_t word = BGEN16 (GET_UINT16 (*b));
            RETURN ((word == VARL_EMPTY   ? WORD_INDEX_EMPTY
                   : word == VARL_MISSING ? WORD_INDEX_MISSING
                   :                        (word & 0x3fff) + VARL_MIN_2B), 2);  
        }

        else if ((msb >> 5) == 0b110) { // 3 bytes (21 bits)
            WordIndex wi = (BGEN24 (GET_UINT24 (*b)) & 0x1fffff) + VARL_MIN_3B; // if embedding this expression directly into the RETURN macro the & operation seems to be ignored. I can't figure out why.
            RETURN (wi, 3);
}
        else { // 4 bytes (29 bits)
            WordIndex wi = BGEN32 (GET_UINT32 (*b)) & 0x1fffffff;
            RETURN (wi, 4);
        }
    }

    // case: files up to 15.0.37
    else switch ((*b)[0]) {
        case BASE250_MOST_FREQ0 : RETURN (0, 1);
        case BASE250_MOST_FREQ1 : RETURN (1, 1);
        case BASE250_MOST_FREQ2 : RETURN (2, 1);
        case BASE250_ONE_UP     : RETURN (WORD_INDEX_ONE_UP,  1);
        case BASE250_EMPTY_SF   : RETURN (WORD_INDEX_EMPTY,   1);
        case BASE250_MISSING_SF : RETURN (WORD_INDEX_MISSING, 1);
        default /* 0 ... 249 */ : {
            WordIndex value;
            switch (b250_size) {
                case B250_BYTES_1: 
                    value = (*b)[0]; 
                    RETURN (value, 1);
                case B250_BYTES_2:
                    value = ((uint32_t)(*b)[0] << 8) | (uint32_t)(*b)[1]; // careful not to use BGEN as string might not be aligned to word boundary
                    RETURN (value, 2);
                case B250_BYTES_3:
                    value = ((uint32_t)(*b)[0] << 16) | ((uint32_t)(*b)[1] << 8) | (uint32_t)(*b)[2]; 
                    RETURN (value, 3);
                case B250_BYTES_4:
                    value = ((uint32_t)(*b)[0] << 24) | ((uint32_t)(*b)[1] << 16) | ((uint32_t)(*b)[2] << 8) | (uint32_t)(*b)[3]; 
                    RETURN (value, 4);
                default:
                    ABORT ("Invalid b250_size=%u", b250_size);
            }
        }
        #undef RETURN
    }
}
