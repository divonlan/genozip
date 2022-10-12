// ------------------------------------------------------------------
//   vcf_illum_gtyping.c
//   Copyright (C) 2022-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "seg.h"
#include "piz.h"
#include "context.h"
#include "reconstruct.h"
#include "file.h"
#include "bits.h"

sSTRl(copy_CHROM_snip, 30);
sSTRl(copy_POS_snip, 30);

void vcf_illum_gtyping_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, INFO_ALLELE_A, INFO_ILLUMINA_POS, DID_EOL);
    seg_id_field_init (CTX(INFO_refSNP));

    seg_mux_init (VB, CTX(FORMAT_BAF), 4, VCF_SPECIAL_MUX_BY_ADJ_DOSAGE, true, (MultiplexerP)&vb->mux_BAF, NULL);
    seg_mux_init (VB, CTX(FORMAT_X),   4, VCF_SPECIAL_MUX_BY_ADJ_DOSAGE, true, (MultiplexerP)&vb->mux_X,   NULL);
    seg_mux_init (VB, CTX(FORMAT_Y),   4, VCF_SPECIAL_MUX_BY_ADJ_DOSAGE, true, (MultiplexerP)&vb->mux_Y,   NULL);

    seg_prepare_snip_other (SNIP_COPY, _VCF_CHROM, false, 0, copy_CHROM_snip);
    seg_prepare_snip_other (SNIP_COPY, _VCF_POS, false,   0, copy_POS_snip);

    ctx_set_ltype (VB, LT_SEQUENCE, INFO_PROBE_A, INFO_PROBE_B, DID_EOL);

    // create ILLUMINA_STRAND nodes as we will refer to them by index (nodes will be deleted in vcf_zip_after_vbs if not used)
    ctx_create_node (VB, INFO_ILLUMINA_STRAND, cSTR("TOP"));   // word_index=0
    ctx_create_node (VB, INFO_ILLUMINA_STRAND, cSTR("BOT"));   // word_index=1
    ctx_create_node (VB, INFO_ILLUMINA_STRAND, cSTR("PLUS"));  // word_index=2
    ctx_create_node (VB, INFO_ILLUMINA_STRAND, cSTR("MINUS")); // word_index=3
    ctx_set_store (VB, STORE_INDEX, INFO_ILLUMINA_STRAND, DID_EOL);
}

// <ID=ILLUMINA_STRAND,Number=1,Type=String,Description="Probe strand">
void vcf_seg_ILLUMINA_STRAND (VBlockVCFP vb, ContextP ctx, STRp(strand))
{
    WordIndex wi = seg_by_ctx (VB, STRa(strand), ctx, strand_len); // nodes pre-created in vcf_illum_gtyping_initialize
    ctx_set_last_value (VB, ctx, (int64_t)wi);
}

// <ID=ILLUMINA_CHR,Number=1,Type=String,Description="Chromosome in Illumina manifest">
void vcf_seg_ILLUMINA_CHR (VBlockVCFP vb, ContextP ctx, STRp(chr))
{
    if (str_issame (chr, vb->chrom_name)) {
        seg_by_ctx (VB, STRa(copy_CHROM_snip), ctx, chr_len);
        ctx_set_encountered (VB, ctx); // encoutered means the CHR matches CHROM
    }
    else
        seg_by_ctx (VB, STRa(chr), ctx, chr_len);
}

// <ID=ILLUMINA_POS,Number=1,Type=Integer,Description="Position in Illumina manifest">
void vcf_seg_ILLUMINA_POS (VBlockVCFP vb, ContextP ctx, STRp(pos))
{
    if (str_issame_(STRa(pos), STRtxtw(CTX(VCF_POS)->last_txt)))
        seg_by_ctx (VB, STRa(copy_POS_snip), ctx, pos_len);

    else
        seg_by_ctx (VB, STRa(pos), ctx, pos_len);

    PosType illumina_pos;
    if (str_get_int (STRa(pos), &illumina_pos))
        ctx_set_last_value (VB, ctx, illumina_pos);
}

// <ID=PROBE_A,Number=1,Type=String,Description="Probe base pair sequence">
// The method here captures ~85% of lines. TO DO: improve  
void vcf_seg_PROBE_A (VBlockVCFP vb, ContextP ctx, STRp(probe))
{
    START_TIMER;

    RefLock lock = REFLOCK_NONE;

    // if we have a reference, we use it 
    if (IS_REF_LOADED_ZIP && !z_is_dvcf       &&           // to do: test this for DVCF
        ctx_has_value (VB, INFO_ILLUMINA_POS) &&           // we go by ILLUMINA_POS, not POS
        ctx_encountered_in_line (VB, INFO_ILLUMINA_CHR) && // verified that CHR is the same as CHROM, so vb->chrom_node_index is correct
        ctx_has_value (VB, INFO_ILLUMINA_STRAND)) {

        WordIndex strand = CTX(INFO_ILLUMINA_STRAND)->last_value.i;
        PosType pos = vb->last_int(INFO_ILLUMINA_POS);
        
        Range *range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos - probe_len, probe_len*2 + 1, 
                                          WORD_INDEX_NONE, NULL, (IS_REF_EXT_STORE ? &lock : NULL));
        if (!range) goto fallback;
        
        // test_fwd:
        bool is_rev = false;
        PosType probe_pos = pos - probe_len + (strand==2 || strand==3);
        for (uint32_t i=0; i < probe_len; i++)
            if (probe[i] != REFp (probe_pos + i)) 
                goto test_rev;
        goto done;

        test_rev:
        is_rev = true;
        probe_pos = pos + 1 - (strand==2 || strand==3);
        for (uint32_t i=0; i < probe_len; i++)
            if (COMPLEM[(uint8_t)(probe[probe_len - i - 1])] != REFp (probe_pos + i))    
                goto fallback;

        done:
        if (IS_REF_EXT_STORE)
            bits_set_region (&range->is_set, (probe_pos - range->first_pos), probe_len);

        SNIPi3 (SNIP_SPECIAL, VCF_SPECIAL_PROBE_A, '0'+is_rev, probe_len);
        seg_by_ctx (VB, STRa(snip), ctx, probe_len);
    }

    else fallback: 
        seg_add_to_local_fixed (VB, ctx, STRa(probe), LOOKUP_WITH_LENGTH, probe_len);

    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE

    seg_set_last_txt (VB, ctx, STRa(probe));
    ctx_set_encountered (VB, ctx);

    COPY_TIMER (vcf_seg_PROBE_A);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PROBE_A)
{    
    if (!reconstruct) goto done;

    ConstRangeP range = ref_piz_get_range (vb, gref, false);

    WordIndex strand = CTX(INFO_ILLUMINA_STRAND)->last_value.i;
    PosType pos = vb->last_int(INFO_ILLUMINA_POS);
    
    uint32_t probe_len;
    str_get_uint32 (snip+1, snip_len-1, &probe_len);
    char *probe = BAFTtxt;
    
    // forward
    if (*snip == '0') {
        PosType probe_pos = pos - probe_len + (strand==2 || strand==3);

        for (uint32_t i=0; i < probe_len; i++)
            probe[i] = REFp (probe_pos+i);
    }

    // reverse complement
    else {
        PosType probe_pos = pos + 1 - (strand==2 || strand==3);

        for (uint32_t i=0; i < probe_len; i++)
            probe[probe_len - i - 1] = COMPLEM[(uint8_t)REFp (probe_pos + i)];
    }

    vb->txt_data.len32 += probe_len;

    done:
    return NO_NEW_VALUE;
}

// <ID=PROBE_B,Number=1,Type=String,Description="Probe base pair sequence; not missing for strand-ambiguous SNPs">
void vcf_seg_PROBE_B (VBlockVCFP vb, ContextP ctx, STRp(seq))
{
    // predicted: . or same as PROBE_A except for the final base 
    if (seq_len > 1 && ctx_encountered_in_line (VB, INFO_PROBE_A) &&  
        str_issame_(STRa(seq)-1, STRtxtw(CTX(INFO_PROBE_A)->last_txt)-1))
        
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_PROBE_B, seq[seq_len-1] }, 3, ctx, seq_len);

    // another real sequence - store verbatim
    else if (seq_len > 1)
        seg_add_to_local_fixed (VB, ctx, STRa(seq), LOOKUP_WITH_LENGTH, seq_len);

    // likely a '.'
    else
        seg_by_ctx (VB, STRa(seq), ctx, seq_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PROBE_B)
{    
    if (reconstruct) {
        RECONSTRUCT (Btxt (CTX(INFO_PROBE_A)->last_txt.index), CTX(INFO_PROBE_A)->last_txt.len - 1);
        RECONSTRUCT1 (snip[0]);
    }

    return NO_NEW_VALUE;
}

// <ID=ALLELE_A,Number=1,Type=String,Description="A allele">
// Expected to be either REF followed by a '*' or ALT
void vcf_seg_ALLELE_A (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    // short cut for common case of SNP - 'R'
    if ((value_len == 2 && vb->main_ref_len ==1 && value[0] == vb->main_ref[0] && value[1] == '*') || // short cut for common case of a SNP
        (value_len > 2 && value_len == vb->main_ref_len + 1 && str_issame_(STRa(vb->main_ref), STRa(value)-1) && value[value_len-1] == '*')) {
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_ALLELE_A, 'R' }, 3, ctx, value_len);
        ctx_set_last_value (VB, ctx, (int64_t)'R'); 
    }

    // short cut for common case of SNP - 'A'
    else if ((value_len == 1 && vb->main_alt_len == 1 && value[0] == vb->main_alt[0]) || // short cut for common case of a SNP
             (value_len > 1 && str_issame_(STRa(vb->main_alt), STRa(value)))) {
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_ALLELE_A, 'A' }, 3, ctx, value_len);
        ctx_set_last_value (VB, ctx, (int64_t)'A'); 
    }

    else // for example, ALT with multi-alleles
        seg_by_ctx (VB, STRa(value), ctx, value_len);
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_ALLELE_A)
{    
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (reconstruct)
        switch (*snip) {
            case 'R' : RECONSTRUCT (vb->main_ref, vb->main_ref_len); RECONSTRUCT1 ('*'); break;
            case 'A' : RECONSTRUCT (vb->main_alt, vb->main_alt_len);                     break;
            default  : ASSPIZ (false, "Invalid snip=%.*s", STRf(snip));
        }

    new_value->i = snip[0];
    return HAS_NEW_VALUE;
}    

void vcf_seg_ALLELE_B (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    if (!ctx_encountered_in_line (VB, INFO_ALLELE_A)) // not encounted, or SPECIAL not used
        goto fallback;

    // prediction ALLELE_B will the other allele out of REF or ALT, than ALLELE_A
    if ((CTX(INFO_ALLELE_A)->last_value.i == 'R' && value_len == 1 && value[0] == vb->main_alt[0]) ||
        (CTX(INFO_ALLELE_A)->last_value.i == 'A' && value_len == 2 && value[0] == vb->main_ref[0] && value[1] == '*'))

        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_ALLELE_B }, 2, ctx, value_len);

    else fallback:
        seg_by_ctx (VB, STRa(value), ctx, value_len);
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_ALLELE_B)
{    
    VBlockVCFP vb = (VBlockVCFP)vb_;

    if (reconstruct)
        switch (CTX(INFO_ALLELE_A)->last_value.i) {
            case 'A' : RECONSTRUCT1 (vb->main_ref[0]); RECONSTRUCT1 ('*'); break;
            case 'R' : RECONSTRUCT1 (vb->main_alt[0]);                     break;
            default  : ASSPIZ (false, "Invalid INFO_ALLELE_A.last_value=%"PRId64, CTX(INFO_ALLELE_A)->last_value.i);
        }

    return NO_NEW_VALUE;
}    

static int vcf_seg_adjust_channel_i (VBlockP vb, int channel_i)
{
    char allele_a = ctx_has_value (vb, INFO_ALLELE_A) ? CTX(INFO_ALLELE_A)->last_value.i : 0;

    if (allele_a == 'A' && (channel_i==0 || channel_i==2))
        return 2 - channel_i; // if homozygot and the ALLELE_A/B are switched vs REF/ALT, we switch the channel

    else if (!allele_a && channel_i==2) // without ALLELE_A, both homozygots are assigned channel=0
        return 0;

    return channel_i;
}

void vcf_seg_mux_by_adjusted_dosage (VBlockVCFP vb, ContextP ctx, STRp(baf), const DosageMultiplexer *mux) 
{
    if (z_is_dvcf || !ctx_encountered_in_line (VB, FORMAT_GT))  // note: if a line fails this test, likely all other lines do too
        seg_by_ctx (VB, STRa(baf), ctx, baf_len); 

    else {
        int channel_i = vcf_seg_adjust_channel_i (VB, vcf_seg_get_mux_channel_i (vb, true));
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, (MultiplexerP)mux, channel_i);

        seg_by_ctx (VB, STRa(baf), channel_ctx, baf_len);
        seg_by_ctx (VB, STRa(mux->snip), ctx, 0);
    }
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_ADJ_DOSAGE)
{    
    int channel_i = vcf_seg_adjust_channel_i (vb, vcf_piz_get_mux_channel_i (VB));

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}    

