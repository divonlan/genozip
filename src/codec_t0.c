// ------------------------------------------------------------------
//   codec_t0.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// codec for Ultima's T0:Z

#include "reconstruct.h"
#include "compressor.h"

//--------------
// ZIP side
//--------------

bool codec_t0_data_is_a_fit_for_t0 (VBlockP vb)
{
    #define NUM_HPS_IN_SAMPLE 100  // test at least this number of HPs (or as many as are available)
    #define SUCCESS_CRITERION 0.85 // at least 85% of tested homopolymers are condensable

    ContextP t0_ctx = CTX(OPTION_t0_Z);
    uint32_t count[2] = {}; // [0]=count non-condensable homopolymers [1]=condensable

    for (LineIType line_i=0; line_i < vb->lines.len32 && (count[false] + count[true] < NUM_HPS_IN_SAMPLE); line_i++) {   
        STRw(t0);
        sam_zip_t0 (vb, t0_ctx, line_i, pSTRa (t0), CALLBACK_NO_SIZE_LIMIT, NULL);
        
        if (!t0_len) continue;

        STRw(seq);
        sam_zip_seq (vb, NULL, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, NULL);

        for (unsigned i=0; i < seq_len; i++) {
            unsigned hp_len = homopolymer_len (STRa(seq), i);

            if (hp_len > 1) {
                count[str_is_monochar (&t0[i], hp_len)]++;
                i += hp_len - 1; // skip to end of homopolymer
            }
        }
    }

    bool success = (count[true] || count[false]) && 
                   ((double)count[true] / (double)(count[true] + count[false])) > SUCCESS_CRITERION; 
    
    if (flag.show_qual) 
        printf ("T0: breakdown of first %u homopolymers in %10s: [condensable]=%u [non-condensable]=%u success=%s\n",
                 count[true] + count[false], VB_NAME, count[true], count[false], TF(success));

    return success;
}

void codec_t0_comp_init (VBlockP vb)
{
    ContextP t0_ctx = CTX(OPTION_t0_Z);
    t0_ctx->ltype   = LT_CODEC;
    t0_ctx->lcodec  = CODEC_T0;
}

COMPRESS (codec_t0_compress)
{
    START_TIMER;

    // case: this is our second entry, after soft-failing. Just continue from where we stopped
    if (!soft_fail) goto do_compress;

    // first pass - condese t0 of homopolymers by removing redundant scores
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        STRw(t0);
        get_line_cb (vb, ctx, line_i, pSTRa (t0), CALLBACK_NO_SIZE_LIMIT, NULL);

        if (!t0_len) continue;

        char *next_condensed = t0; // copy in-place;

        STRw(seq);
        sam_zip_seq (vb, NULL, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, NULL);

        for (uint32_t i=0; i < seq_len; i++) {
            unsigned hp_len = homopolymer_len (STRa(seq), i);  
            
            if (hp_len > 1) {
                // test condensability (rarely, homopolymer t0 may be non-condensable)
                if (str_is_monochar (&t0[i], hp_len)) 
                    *next_condensed++ = t0[i];

                else {
                    *next_condensed++ = t0[i] | 0x80; // "not-condensed" flag

                    for (unsigned h=1; h < hp_len; h++)
                        *next_condensed++ = t0[i + h];
                }

                i += hp_len - 1; // skip to end of homopolymer
            }
            
            else
                *next_condensed++ = t0[i];
        }

        uint32_t condensed_len = next_condensed - t0;

        if (condensed_len != t0_len) {
            // shrink t0_len because sub_codec uses the callback to get the condensed t0
            sam_ultima_update_t0_len (vb, line_i, condensed_len); 
            ctx->local.len32 -= t0_len - condensed_len;
        }
    }

    // second pass - compress t0 data with redundant homopolymer t0 scores removed
    ctx->lcodec = CODEC_UNKNOWN;
    header->sub_codec = codec_assign_best_codec (vb, ctx, NULL, SEC_LOCAL); 
    if (header->sub_codec == CODEC_UNKNOWN) header->sub_codec = CODEC_NONE; // really small

do_compress: ({});
    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = ctx->local.len32;

    // make sure we have enough memory
    uint32_t min_required_compressed_len = codec_args[header->sub_codec].est_size (header->sub_codec, ctx->local.len);
    if (*compressed_len < min_required_compressed_len) {
        if (soft_fail) return false; // call me again with more memory
        ABORT ("%s: Compressing %s with %s need %u bytes, but allocated only %u", 
               VB_NAME, ctx->tag_name, codec_name(header->sub_codec), min_required_compressed_len, *compressed_len);
    }

    COPY_TIMER_COMPRESS (compressor_t0); // don't account for sub-codec compressor, it accounts for itself

    return compress (vb, ctx, header, 0, uncompressed_len, get_line_cb, STRa(compressed), false, name);
}

//--------------
// PIZ side
//--------------

CODEC_RECONSTRUCT (codec_t0_reconstruct)
{
    START_TIMER;
 
    rom seq = sam_piz_get_textual_seq (vb);

    ASSPIZ (len == vb->seq_len, "expecting len=%u == vb->seq_len=%u", len, vb->seq_len);

    // reconsutrct by expanding homopolymers 
    char *next_recon = BAFTtxt;
    char *next_condensed = Bc(ctx->local, ctx->next_local);

    for (uint32_t i=0; i < len; i++) {
        unsigned hp_len = homopolymer_len (seq, len, i);  

        if (hp_len > 1) {
            // case: non-condensable
            if (*next_condensed & 0x80) {
                *next_recon++ = *next_condensed++ & 0x7f;
                for (unsigned h=1; h < hp_len; h++)
                    *next_recon++ = *next_condensed++;
            }

            else { 
                char c = *next_condensed++;
                for (unsigned h=0; h < hp_len; h++)
                    *next_recon++ = c;
            }

            i += hp_len - 1; // skip to end of homopolymer
        }
        
        else
            *next_recon++ = *next_condensed++;
    }

    if (reconstruct) Ltxt += len;
    
    ctx->next_local = BNUM (ctx->local, next_condensed);

    COPY_TIMER(codec_t0_reconstruct);
}
