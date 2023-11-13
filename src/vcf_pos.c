// ------------------------------------------------------------------
//   vcf_pos.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
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
    if (vb->line_coords == DC_PRIMARY) {
        PosType64 pos;

        if (segconf.vcf_is_gvcf) {
            // note: using a multiplexor for distinguising END and POS, while keeping them as an alias
            // has the advantage the that delta=1 snip common in GVCF will be the same self-delta snip
            // regardless if previous line has an END or not (i.e. delta might be against the END or the POS, but these will result in the same snip)
            ContextP subctx = seg_mux_get_channel_ctx (VB, VCF_POS, (MultiplexerP)&vb->mux_POS, 0); // goes into channel_i=0: "this is POS"

            pos = dl->pos[0] = seg_pos_field (VB, subctx->did_i, VCF_POS, 0, '.', STRa(pos_str), 0, pos_str_len+1);
            ctx_set_last_value (VB, CTX(VCF_POS), pos);
        
            seg_by_did (VB, STRa(vb->mux_POS.snip), VCF_POS, 0); // de-multiplexer
        }

        else 
            pos = dl->pos[0] = seg_pos_field (VB, VCF_POS, VCF_POS, 0, '.', STRa(pos_str), 0, pos_str_len+1);
        
        if (pos == 0 && !(*pos_str == '.' && pos_str_len == 1)) // POS == 0 - invalid value return from seg_pos_field
            WARN_ONCE ("FYI: invalid POS=%"PRId64" value in chrom=%.*s vb_i=%u vb_line_i=%d: line will be compressed, but not indexed", 
                       pos, vb->chrom_name_len, vb->chrom_name, vb->vblock_i, vb->line_i);
                
        if (pos) random_access_update_pos (VB, 0, VCF_POS);
    }

    else { // LUFT
        dl->pos[1] = seg_pos_field (VB, VCF_oPOS, VCF_oPOS, 0, '.', STRa(pos_str), 0, pos_str_len);
        if (dl->pos[1]) random_access_update_pos (VB, 1, VCF_oPOS);
        CTX(vb->vb_coords==DC_LUFT ? VCF_oPOS : VCF_POS)->txt_len++; // account for the tab - in oPOS in the ##luft_only VB and in POS (on behalf on the primary POS) if this is a Dual-coord line (we will rollback accounting later if its not)
    }

    set_last_txt_(VCF_POS, pos_str, pos_str_len); // consumed by vcf_seg_FORMAT_PS, vcf_seg_ILLUMINA_POS
}

// --------
// INFO/END
// --------

static void vcf_seg_INFO_END_liftover (VBlockVCFP vb, ContextP end_ctx, STRp(end_str))
{
    bool is_xstrand = (vb->last_index (VCF_oXSTRAND) > 0); // set in vcf_lo_seg_generate_INFO_DVCF
    PosType64 aln_last_pos = chain_get_aln_prim_last_pos (vb->pos_aln_i); 
    PosType64 end = vb->last_int (INFO_END); 

    // case: we don't yet handle END translation in case of a reverse strand
    if (is_xstrand)            
        REJECT_SUBFIELD (LO_INFO, end_ctx, ".\tVariant with INFO/END and chain file alignment with a negative strand%s", "");

    // case: END goes beyond the end of the chain file alignment and its a <DEL>
    else if (vb->is_del_sv && end > aln_last_pos) {

        // case: END goes beyond end of alignment
        PosType64 gap_after = chain_get_aln_gap_after (vb->pos_aln_i);
        
        // case: END falls in the gap after - <DEL> is still valid but translated END needs to be closer to POS to avoid gap - 
        // we don't yet do this
        if (end <= aln_last_pos + gap_after)
            REJECT_SUBFIELD (LO_INFO, end_ctx, ".\t<DEL> variant: INFO/END=%.*s is in the gap after the end of the chain file alignment", end_str_len, end_str);

        // case: END falls on beyond the gap (next alignment or beyond) - this variant cannot be lifted
        else
            REJECT_SUBFIELD (LO_INFO, end_ctx, ".\t<DEL> variant: INFO/END=%.*s is beyond the end of the chain file alignment and also beyond the gap after the alignment", end_str_len, end_str);
    }

    // case: END goes beyond the end of the chain file alignment and its NOT a <DEL>
    else if (!vb->is_del_sv && end > aln_last_pos) 
        REJECT_SUBFIELD (LO_INFO, end_ctx, ".\tPOS and INFO/END=%.*s are not on the same chain file alignment", end_str_len, end_str);

    // case: invalid value. since we use SPF_UNLIMITED_DELTA, any integer value should succeed
    else if (!CTX(INFO_END)->last_delta)
        REJECT_SUBFIELD (LO_INFO, end_ctx, ".\tINFO/END=%.*s has an invalid value", end_str_len, end_str);        

} 

void vcf_seg_INFO_END (VBlockVCFP vb, ContextP end_ctx, STRp(end_str)) // end_ctx is INFO_END, despite being an alias
{
    bool is_liftover = chain_is_loaded && LO_IS_OK (last_ostatus);

    // END is an alias of POS
    if (segconf.vcf_is_gvcf) {
        ContextP subctx = seg_mux_get_channel_ctx (VB, VCF_POS, (MultiplexerP)&vb->mux_POS, 1); // goes into channel_i=1: "this is END"

        PosType64 end = seg_pos_field (VB, subctx->did_i, VCF_POS, SPF_BAD_SNIPS_TOO | SPF_ZERO_IS_BAD | (is_liftover ? SPF_UNLIMITED_DELTA : 0), 0, STRa(end_str), 0, end_str_len);
        ctx_set_last_value (VB, CTX(VCF_POS), end); // END is an alias of POS

        seg_by_did (VB, STRa(vb->mux_POS.snip), VCF_POS, 0); // de-multiplexer
    }
    
    else  
        seg_pos_field (VB, VCF_POS, VCF_POS, SPF_BAD_SNIPS_TOO | SPF_ZERO_IS_BAD | SPF_UNLIMITED_DELTA, 0, STRa(end_str), 0, end_str_len);

    // add end_delta to dl for sorting. it is used only in case chrom and pos are identical
    DATA_LINE (vb->line_i)->end_delta = vb->last_delta (INFO_END);

    // case --chain: if we have lifted-over POS (as primary POS field or in INFO/LIFTBACK), 
    // check that lifting-over of END is delta-encoded and is lifted over to the same, non-xstrand, Chain alignment, and reject if not
    if (is_liftover) 
        vcf_seg_INFO_END_liftover (vb, end_ctx, STRa(end_str));
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_END)
{
    int channel_i = (vb->con_stack_len == 2); // 1 if INFO_END (stack=TOPLEVEL->INFO), and 0 if VCF_POS (stack=TOPLEVEL)

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// END data resides in POS (its an alias), but has a different translator as its a different container item. 
// For END, we didn't add an oPOS entry, because we can't consume it when showing Primary. Instead, we do delta arithmetic.
// returns true if successful (return value used only if validate_only)
TRANSLATOR_FUNC (vcf_piz_luft_END)
{
    // ZIP liftover validation: postpone to vcf_seg_INFO_END
    if (validate_only) return true; 

    PosType64 translated_end;
    ContextP pos_ctx  = CTX (VCF_POS);
    ContextP opos_ctx = CTX (VCF_oPOS);
    
    // ZIP liftback (POS is always before END, because we seg INFO/LIFTOVER first)
    if (IS_ZIP && VB_VCF->line_coords == DC_LUFT) { // liftback
        PosType64 oend;
        if (!str_get_int_range64 (STRa(recon), 0, MAX_POS, &oend))
            return false;

        translated_end = pos_ctx->last_value.i + (oend - opos_ctx->last_value.i);
    }

    // PIZ liftover: we have already reconstructed oPOS and POS (as an item in VCF_TOPLUFT)
    else {
        translated_end = opos_ctx->last_value.i + pos_ctx->last_delta; ; // delta for generated this END value (END - POS)

        CTX(INFO_END)->last_end_line_i = vb->line_i; // so vcf_piz_special_COPYPOS knows that END was reconstructed
    }

    // re-reconstruct END
    Ltxt -= recon_len; 
    RECONSTRUCT_INT (translated_end);
    
    return true;    
}

// Called to reconstruct the POS subfield of INFO/LIFTBACK, handling the possibility of INFO/END
// Also used for FORMAT/PS
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

    RECONSTRUCT_INT (pos); // vcf_piz_luft_END makes sure it always contains the value of POS, not END
    return NO_NEW_VALUE;
}


