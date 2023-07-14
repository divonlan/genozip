// ------------------------------------------------------------------
//   vcf_giab.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "context.h"
#include "strings.h"
#include "dict_id.h"
#include "reconstruct.h"
#include "lookback.h"

void vcf_giab_seg_initialize (VBlockVCFP vb)
{
    seg_mux_init (VB, CTX(FORMAT_IGT), 2, VCF_SPECIAL_MUX_BY_SAMPLE_I,  false, (MultiplexerP)&vb->mux_IGT, NULL);
    seg_mux_init (VB, CTX(FORMAT_IPS), 2, VCF_SPECIAL_MUX_BY_IGT_PHASE, false, (MultiplexerP)&vb->mux_IPS, NULL);
}

//--------------------------------------------------------------------------------
// FORMAT/IPS: <ID=IPS,Number=1,Type=String,Description="Phase set for IGT">
//--------------------------------------------------------------------------------

// TO DO: store current_phase (only in vb->sample_i==0): in Seg, can store in last_txt, 
// in Piz, to avoid needing to deal with txt_data insertions, store in a "special_store" context-specific union of ol_nodes.
// Seg: if same as current_phase, seg SPECIAL_COPY_STORED. If not, check if same as pos,ref,alt and seg SPECIAL_IPS,
// with a parameter of whether ref and alt appear rev-comped.
// (no need for lookback as IPS is expected to always be either equal to pos/ref/alt, or to previous phase)
void vcf_seg_FORMAT_IPS (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP ctx, STRp(ips))
{
    STRlast (igt, FORMAT_IGT);
    bool is_phased = (ctx_encountered (VB, FORMAT_IGT) && igt_len == 3 && igt[1] == '|');

    ContextP channel_ctx = 
        seg_mux_get_channel_ctx (VB, FORMAT_IPS, (MultiplexerP)&vb->mux_IPS, is_phased);

    seg_by_ctx (VB, STRa(ips), channel_ctx, ips_len); // note: for channel_i=0 - expeceted to be '.'

    seg_by_ctx (VB, STRa(vb->mux_IPS.snip), ctx, 0);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_IGT_PHASE)
{    
    STRlast (igt, FORMAT_IGT);
    bool is_phased = (ctx_encountered (VB, FORMAT_IGT) && igt_len == 3 && igt[1] == '|');

    HasNewValue ret = reconstruct_demultiplex (vb, ctx, STRa(snip), is_phased, new_value, reconstruct);

    return ret;
}    

//--------------------------------------------------------------------------------
// FORMAT/IGT: <ID=IGT,Number=1,Type=String,Description="Original input genotype">
//--------------------------------------------------------------------------------

void vcf_seg_FORMAT_IGT (VBlockVCFP vb, ContextP ctx, STRp(igt))
{
    seg_set_last_txt (VB, ctx, STRa(igt)); // consumed by vcf_seg_FORMAT_IPS

    if (!ctx_encountered (VB, FORMAT_GT) || vcf_num_samples != 3 || igt_len > 3 ||
        igt_len != CTX(FORMAT_GT)->gt_prev_ploidy * 2 - 1) {
        seg_by_ctx (VB, STRa(igt), ctx, igt_len);
        return;
    }

    ContextP channel_ctx = 
        seg_mux_get_channel_ctx (VB, FORMAT_IGT, (MultiplexerP)&vb->mux_IGT, (vb->sample_i > 0));

    if (vb->sample_i == 0) {
        Allele *gt = B(Allele, CTX(FORMAT_GT_HT)->local, vb->line_i * vb->ht_per_line + vb->ploidy * vb->sample_i);
        
        // case monoploid: predicting GT=IGT
        if (str_is_1char (igt, *gt))
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_IGT }, 2, channel_ctx, igt_len);     

        // case: unphased diploid IGT: predicting GT=IGT, possibly flipping the order so that the smaller allele appears first
        else if (igt[1] == '/' && (   (gt[0] <= gt[1] && igt[0] == gt[0] && igt[2] == gt[1])
                                   || (gt[0] >  gt[1] && igt[0] == gt[1] && igt[2] == gt[0])))
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_IGT, '/' }, 3, channel_ctx, igt_len);     

        // case phased diploid IGT: predicting IGT alleles = GT alleles
        else if (igt[1] == '|' && igt[0] == gt[0] && igt[2] == gt[1])
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_IGT, '|' }, 3, channel_ctx, igt_len);     

        // case phased diploid IGT: predicting IGT alleles = GT alleles, but in reverse order
        else if (igt[1] == '|' && igt[0] == gt[1] && igt[2] == gt[0])
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_IGT, ':' }, 3, channel_ctx, igt_len);     

        else
            goto fallback;
    }

    else fallback:
        seg_by_ctx (VB, STRa(igt), channel_ctx, igt_len); // note: for channel_i=1 - expeceted to be '.'

    seg_by_ctx (VB, STRa(vb->mux_IGT.snip), ctx, 0);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_SAMPLE_I)
{    
    int channel_i = (vb->sample_i > 0);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}    

SPECIAL_RECONSTRUCTOR (vcf_piz_special_IGT)
{
    STRlast(gt, FORMAT_GT);

    if (gt_len == 1)
        RECONSTRUCT1 (*gt);

    else if ((snip[0] == '/' && gt[0] <= gt[2]) || snip[0] == '|') {
        RECONSTRUCT1 (gt[0]);
        RECONSTRUCT1 (snip[0]);
        RECONSTRUCT1 (gt[2]);
    }

    else { // flip allele order
        RECONSTRUCT1 (gt[2]);
        RECONSTRUCT1 (snip[0]=='/' ? '/' : '|');
        RECONSTRUCT1 (gt[0]);
    }

    return NO_NEW_VALUE;
}

//--------------------------------------------------------------------------------------------------------
// FORMAT/ADALL: <ID=ADALL,Number=R,Type=Integer,Description="Net allele depths across all datasets">
//--------------------------------------------------------------------------------------------------------

// Sepcial treatment for item 0
void vcf_seg_ADALL_items (VBlockVCFP vb, ContextP ctx, STRps(item), ContextP *item_ctxs, const int64_t *values)
{
    for (unsigned i=0; i < n_items; i++) 
        if (i==0 || i==1) {
            if (!vb->mux_ADALL[i].num_channels)
                seg_mux_init (VB, item_ctxs[i], 4, VCF_SPECIAL_MUX_BY_DOSAGE, false, (MultiplexerP)&vb->mux_ADALL[i], "0123");
            
            vcf_seg_FORMAT_mux_by_dosage (vb, item_ctxs[i], STRi(item, i), &vb->mux_ADALL[i]);
        }
        else 
            seg_by_ctx (VB, STRi(item, i), item_ctxs[i], item_lens[i]);
}
