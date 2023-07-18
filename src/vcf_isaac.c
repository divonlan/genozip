// ------------------------------------------------------------------
//   sam_isaac.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "vcf_private.h"
#include "reconstruct.h"

void vcf_isaac_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, FORMAT_GQ, FORMAT_GQX, DID_EOL); 

    // GQX stuff
    seg_mux_init (VB, CTX(FORMAT_GQX), 3, VCF_SPECIAL_DEMUX_GQX, false, (MultiplexerP)&vb->mux_GQX, "012"); 
    seg_by_did (VB, STRa(vb->mux_GQX.snip), FORMAT_GQX, 0); // all the same
}

void vcf_seg_FORMAT_GQX (VBlockVCFP vb, ContextP ctx, STRp(gqx))
{
    Allele *ht_data = B(Allele, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line + vb->ploidy * vb->sample_i);

    int channel_i = (!ctx_encountered (VB, FORMAT_GT) || *ht_data == '.') ? 0
                  : !ctx_encountered (VB, FORMAT_GQ)                      ? 1
                  :                                                         2;

    ContextP chan_ctx = seg_mux_get_channel_ctx (VB, FORMAT_GQX, (MultiplexerP)&vb->mux_GQX, channel_i);
    
    if (channel_i == 2) 
        seg_delta_vs_other (VB, chan_ctx, CTX(FORMAT_GQ), STRa(gqx));

    else
        seg_by_ctx (VB, STRa(gqx), chan_ctx, gqx_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_DEMUX_GQX)
{
    rom gt = last_txt (VB, FORMAT_GT);

    int channel_i = (!ctx_encountered (VB, FORMAT_GT) || *gt == '.') ? 0
                  : !ctx_encountered (VB, FORMAT_GQ)                 ? 1
                  :                                                    2;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}
