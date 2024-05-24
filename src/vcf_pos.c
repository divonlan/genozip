// ------------------------------------------------------------------
//   vcf_pos.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "reconstruct.h"
#include "random_access.h"

void vcf_seg_pos (VBlockVCFP vb, ZipDataLineVCF *dl, STRp(pos_str))
{
    PosType64 pos;

    if (segconf.vcf_is_gvcf) {
        // note: using a multiplexor for distinguising END and POS, while keeping them as an alias
        // has the advantage the that delta=1 snip common in GVCF will be the same self-delta snip
        // regardless if previous line has an END or not (i.e. delta might be against the END or the POS, but these will result in the same snip)
        ContextP subctx = seg_mux_get_channel_ctx (VB, VCF_POS, (MultiplexerP)&vb->mux_POS, 0); // goes into channel_i=0: "this is POS"

        pos = dl->pos = seg_pos_field (VB, subctx->did_i, VCF_POS, 0, '.', STRa(pos_str), 0, pos_str_len+1);
        ctx_set_last_value (VB, CTX(VCF_POS), pos);
    
        seg_by_did (VB, STRa(vb->mux_POS.snip), VCF_POS, 0); // de-multiplexer
    }

    else 
        pos = dl->pos = seg_pos_field (VB, VCF_POS, VCF_POS, 0, '.', STRa(pos_str), 0, pos_str_len+1);
    
    if (pos == 0 && !(*pos_str == '.' && pos_str_len == 1)) // POS == 0 - invalid value return from seg_pos_field
        WARN_ONCE ("FYI: invalid POS=%"PRId64" value in chrom=%.*s vb_i=%u vb_line_i=%d: line will be compressed, but not indexed", 
                    pos, vb->chrom_name_len, vb->chrom_name, vb->vblock_i, vb->line_i);
            
    if (pos) random_access_update_pos (VB, VCF_POS);

    set_last_txt_(VCF_POS, pos_str, pos_str_len); // consumed by vcf_seg_FORMAT_PS, vcf_seg_ILLUMINA_POS
}

// --------
// INFO/END
// --------

void vcf_seg_INFO_END (VBlockVCFP vb, ContextP end_ctx, STRp(end_str)) // end_ctx is INFO_END, despite being an alias
{
    // END is an alias of POS
    if (segconf.vcf_is_gvcf) {
        ContextP subctx = seg_mux_get_channel_ctx (VB, VCF_POS, (MultiplexerP)&vb->mux_POS, 1); // goes into channel_i=1: "this is END"

        PosType64 end = seg_pos_field (VB, subctx->did_i, VCF_POS, SPF_BAD_SNIPS_TOO | SPF_ZERO_IS_BAD, 0, STRa(end_str), 0, end_str_len);
        ctx_set_last_value (VB, CTX(VCF_POS), end); // END is an alias of POS

        seg_by_did (VB, STRa(vb->mux_POS.snip), VCF_POS, 0); // de-multiplexer
    }
    
    else  
        seg_pos_field (VB, VCF_POS, VCF_POS, SPF_BAD_SNIPS_TOO | SPF_ZERO_IS_BAD, 0, STRa(end_str), 0, end_str_len);

    ctx_set_encountered (VB, end_ctx); // encountered END (but segged into POS and set last_value of POS as its an alias)
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_END)
{
    int channel_i = (vb->con_stack_len == 2); // 1 if INFO_END (stack=TOPLEVEL->INFO), and 0 if VCF_POS (stack=TOPLEVEL)

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// used for FORMAT/PS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_COPYPOS)
{
    if (!reconstruct) return NO_NEW_VALUE;
    
    bool has_end = CTX(INFO_END)->last_end_line_i == vb->line_i; // true if INFO/END was encountered

    ContextP pos_ctx = CTX (VCF_POS);
    int64_t pos;

    if (has_end) {
        int64_t end   = pos_ctx->last_value.i;
        int64_t delta = pos_ctx->last_delta;
        pos = end - delta;
    }
    else
        pos = pos_ctx->last_value.i;

    if (snip_len)
        pos += atoi (snip); // add optional delta (since 13.0.5)

    RECONSTRUCT_INT (pos); 
    return NO_NEW_VALUE;
}


