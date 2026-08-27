// ------------------------------------------------------------------
//   codec_longr.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// This codec is for quality scores of Nanopore data, based on https://pubmed.ncbi.nlm.nih.gov/32470109/. This file contains 
// only original Genozip code, and the algorithm itself, derived from ENano source code, is located in codec_longr_alg.c

#include "codec.h"
#include "reconstruct.h"
#include "compressor.h"
#include "context.h"
#include "piz.h"
#include "stats.h"

#include "codec_longr_alg.c" // separate source file for this, as it derived from external code with a different license

//--------------
// PIZ side
//--------------

static bool codec_longr_recon_one_read (LongrState *state, STRp(seq), bool is_rev,
                                        bytes sorted_qual, uint32_t *next_of_chan, char *recon)
{
    codec_longr_alg_init_read (state, STRa(seq), is_rev);
    uint8_t prev_q = 0;

    #define RECON_ONE_QUAL                              \
        codec_longr_update_state (state, b, q, prev_q); \
        if (q == 255) return false; /* 255+'!' == ' ' == missing qual */ \
        prev_q = q;                                     \
        recon[i] = q  + '!';

    if (!is_rev) // separate loops to save one "if" in the tight loop
        for (uint32_t i=0; i < seq_len; i++) {        
            uint8_t b = acgt_encode(codec_longr_next_base (STRa(seq), i));
            uint8_t q = sorted_qual[next_of_chan[state->chan.channel.n]++];
            RECON_ONE_QUAL;
        }
    else
        for (int32_t i=seq_len-1; i >= 0; i--) {        
            uint8_t b = acgt_encode_comp (codec_longr_next_base_rev (STRa(seq), i));
            uint8_t q = sorted_qual[next_of_chan[state->chan.channel.n]++];
            RECON_ONE_QUAL;
        }

    return true; 
}

// order of decompression: lens_ctx is decompressed, then baseq_ctx is decompressed with its codec, and then this function is called
// as a subcodec for baseq_ctx. 
// This function converts lens_ctx to be "next_of_chan" for each channel, and initializes LongrState
static void codec_longr_reconstruct_init (VBlockP vb, ContextP lens_ctx, ContextP values_ctx)
{
    // we adjust the buffer here, since it didn't get adjusted in piz_uncompress_all_ctxs because its ltype is LT_CODEC
    lens_ctx->local.len /= sizeof (uint32_t); 
    BGEN_u32_buf (&lens_ctx->local, NULL);

    ARRAY (uint32_t, next_of_chan, lens_ctx->local);

    // transform len array to next array
    uint32_t next=0;
    for (uint32_t chan=0; chan < next_of_chan_len; chan++) {
        uint32_t len = next_of_chan[chan];
        next_of_chan[chan] = next;
        next += len;
    }
    
    // retrieve the global value-to-bin mapper from SEC_COUNTS and store it in values_ctx.value_to_bin
    buf_alloc (vb, &values_ctx->value_to_bin, 0, 256, uint8_t, 0, "value_to_bin");
    values_ctx->value_to_bin.len = 256;

    // until 13.0.11, stored in SEC_COUNTS of lens_ctx, and after in values_ctx
    BufferP value_to_bin = ZCTX(values_ctx->did_i)->counts.len ? &ZCTX(values_ctx->did_i)->counts 
                                                               : &ZCTX(lens_ctx->did_i)->counts;

    ASSERTISALLOCED (*value_to_bin);
    
    ARRAY (uint64_t, value_to_bin_src, *value_to_bin); 
    ARRAY (uint8_t,  value_to_bin_dst, values_ctx->value_to_bin);

    for (int i=0; i < 256; i++) 
        value_to_bin_dst[i] = value_to_bin_src[i]; // uint64 -> uint8

    // initialize longr state - stored in lens_ctx.longr_state
    buf_alloc_zero (vb, &lens_ctx->longr_state, 1, 0, LongrState, 0, C_"longr_state"); 
    codec_longr_alg_init (B1ST (LongrState, lens_ctx->longr_state));

    lens_ctx->is_initialized = true;
}

// When reconstructing a QUAL field on a specific line, piz calls the LT_CODEC reconstructor for CODEC_LONGR, 
// codec_longr_reconstruct, which combines data from the local buffers of lens, values and SQBITMAP to reconstruct the original QUAL field.
CODEC_RECONSTRUCT (codec_longr_reconstruct)
{
    START_TIMER;

    ContextP lens_ctx   = ctx;
    ContextP values_ctx = ctx + 1;
    
    if (!lens_ctx->is_initialized) 
        codec_longr_reconstruct_init (vb, lens_ctx, values_ctx);

    bool is_rev = !VB_DT(FASTQ) /* SAM or BAM */ && sam_is_last_flags_rev_comp(vb);

    ARRAY (uint8_t, sorted_qual, values_ctx->local);
    ARRAY (uint32_t, next_of_chan, lens_ctx->local);
    LongrState *state = B1ST (LongrState, lens_ctx->longr_state);
    state->value_to_bin = B1ST8 (values_ctx->value_to_bin);

    rom seq = VB_DT(SAM) ? sam_piz_get_textual_seq(vb) : last_txtx (vb, CTX(FASTQ_SQBITMAP)); 

    // case: Deep, and len is only the trimmed suffix as the rest if copied from SAM (see fastq_special_deep_copy_QUAL)
    if (flag.deep && len < vb->seq_len)
        seq += (vb->seq_len - len); // advance seq to the trimmed part too
    else
        ASSPIZ (len == vb->seq_len, "expecting len=%u == vb->seq_len=%u", len, vb->seq_len);
    
    if (codec_longr_recon_one_read (state, seq, len, is_rev, sorted_qual, next_of_chan, BAFTtxt)) {
        if (reconstruct) Ltxt += len;
    } else // missing qual
        sam_reconstruct_missing_quality (vb, reconstruct);

    COPY_TIMER(codec_longr_reconstruct);
}

