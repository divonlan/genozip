// ------------------------------------------------------------------
//   sam_dragen.c
//   Copyright (C) 2023-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"

void sam_dragen_seg_initialize (VBlockSAMP vb)
{
    seg_mux_init (vb, OPTION_sd_f, SAM_SPECIAL_sd, false, dragen_sd);
}

static int sd_channel_i (int seq_len, int as)
{
    return (seq_len == as);
}

void sam_dragen_seg_sd_f (VBlockSAMP vb, ZipDataLineSAMP dl, STRp(sd), ValueType numeric, unsigned add_bytes)
{
    if (dl->AS) { 
        int channel_i = sd_channel_i (dl->SEQ.len, dl->AS);
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_sd_f, (MultiplexerP)&vb->mux_dragen_sd, channel_i);

        sam_seg_float_as_snip (vb, channel_ctx, STRa(sd), numeric, add_bytes);

        seg_by_did (VB, STRa(vb->mux_dragen_sd.snip), OPTION_sd_f, 0); // de-multiplexer
    }

    // no AS in line (not expected in Dragen) or AS=0 (not expected either)
    else
        sam_seg_float_as_snip (vb, CTX(OPTION_sd_f), STRa(sd), numeric, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_sd)
{
    ContextP as_ctx;
    ASSPIZ0 (ctx_has_value_in_line (vb, _OPTION_AS_i, &as_ctx), "AS:i was not reconstructed for this line");

    int channel_i = sd_channel_i (vb->seq_len, as_ctx->last_value.i);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}
