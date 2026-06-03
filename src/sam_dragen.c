// ------------------------------------------------------------------
//   sam_dragen.c
//   Copyright (C) 2023-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"

void sam_dragen_seg_initialize (VBlockSAMP vb)
{
    seg_mux_init (vb, OPTION_sd_f, SAM_SPECIAL_sd, false, dragen_sd);

    ctx_consolidate_stats (VB, OPTION_ga_Z, OPTION_ga_CONTIG, OPTION_ga_POS, OPTION_ga_STRAND, OPTION_ga_CIGAR, OPTION_ga_MAPQ, OPTION_ga_NM, DID_EOL);

    ctx_set_no_stons (VB, OPTION_ga_STRAND, DID_EOL);
}

static int sd_channel_i (int seq_len, int as)
{
    return (seq_len == as);
}

void sam_dragen_seg_sd_f (VBlockSAMP vb, ZipDataLineSAM𐤐 dl, STRp(sd), ValueType numeric, unsigned add_bytes)
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

// ga:Z is: "contig, pos, strand, CIGAR, mapQ, NM;". contig is not necessarily (and usually not) in the SAM header.
// Example: ga:Z:Edico_decoy_GRCh38_chrUn_JTFH01000594v1_decoy,1210,-,17S76M1D57M,0,12;
// Example: ga:Z:Edico_decoy_GRCh38_hub_3267197_GCA_009914755.4_hub_3267197_assembly_CP068263.2:0-15676637,5217078,+,4S147M,0,13; 
// See: https://help.dragen.illumina.com/product-guides/dragen-v4.5/dragen-host-software
void sam_dragen_seg_ga_Z (VBlockSAMP vb, STRp(ga), unsigned add_bytes)
{
    decl_ctx (OPTION_ga_Z);

    static const Container(6) container_ga = { 
        .nitems_lo = 6,
        .repeats   = 1,          
        .items     = { { .dict_id = { _OPTION_ga_CONTIG }, .separator = {','} }, // usually, not in the RNAME dict (extra contigs Illumina provides beyond the user's FASTA reference)
                       { .dict_id = { _OPTION_ga_POS    }, .separator = {','} },  
                       { .dict_id = { _OPTION_ga_STRAND }, .separator = {','} },  
                       { .dict_id = { _OPTION_ga_CIGAR  }, .separator = {','} },  
                       { .dict_id = { _OPTION_ga_MAPQ   }, .separator = {','} },  
                       { .dict_id = { _OPTION_ga_NM     }, .separator = {';'} } } 
    };

    static SegCallback callbacks[6] = { [SA_POS]   = sam_seg_0A_pos_cb, 
                                        [SA_CIGAR] = sam_seg_0A_cigar_cb, 
                                        [SA_MAPQ]  = sam_seg_0A_mapq_cb   };
     
    seg_struct (VB, ctx,(ContainerP)&container_ga, STRa(ga), callbacks, add_bytes, true);

    seg_set_last_txt (VB, ctx, STRa(ga));
}
