// ------------------------------------------------------------------
//   vcf_platypus.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "vcf_private.h"
#include "reconstruct.h"

sSTRl (tcr_snip, 32);

void vcf_platypus_zip_initialize (void)
{
    seg_prepare_minus_snip (VCF, _INFO_TC, _INFO_TCF, tcr_snip);
}

void vcf_platypus_seg_initialize (VBlockVCFP vb)
{
    ctx_set_store (VB, STORE_INT, INFO_MGOF, FORMAT_GOF, INFO_WS, INFO_WE,
                   INFO_TC, INFO_TCF, INFO_TCR, INFO_PP, INFO_HP, INFO_TR,
                   DID_EOL);

    ctx_set_ltype (VB, LT_DYN_INT, INFO_HP, DID_EOL);
}

// ---------------------------------------------------------
// INFO/SC
// Genomic sequence 10 bases either side of variant position
// ---------------------------------------------------------

void vcf_seg_playpus_INFO_SC (VBlockVCFP vb, ContextP ctx, STRp(seq))
{
    decl_acgt_decode;
    RefLock lock = REFLOCK_NONE;

    if (seq_len != 21 || !flag.reference || segconf.running ||
        !str_is_ACGT (STRa(seq), NULL))  // reference doesn't support N or IUPACs
        goto fallback;

    PosType64 pos = DATA_LINE(vb->line_i)->pos[0] - 10; // 10 before to 10 after

    RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, seq_len, WORD_INDEX_NONE, 
                                      (IS_REF_EXT_STORE ? &lock : NULL));
    
    if (!range || pos < range->first_pos || pos + 20 > range->last_pos)
        goto fallback;

    // verify that SC matches the reference
    for (int i=0; i < seq_len; i++)
        if (seq[i] != REFp (pos + i)) 
            goto fallback;
    
    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_PLATYPUS_SC }, 2, ctx, seq_len);

    if (IS_REF_EXT_STORE)
        bits_set_region (&range->is_set, pos - range->first_pos, seq_len);

    ref_unlock (gref, &lock); 
    return;

fallback: 
    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
    seg_by_ctx (VB, STRa(seq), ctx, seq_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PLATYPUS_SC)
{    
    if (!reconstruct) return NO_NEW_VALUE;

    ConstRangeP range = ref_piz_get_range (vb, gref, HARD_FAIL);

    PosType64 pos = CTX(VCF_POS)->last_value.i - 10;
    
    char *next = BAFTtxt;
    decl_acgt_decode;
    
    for (uint32_t i=0; i < 21; i++)
        *next++ = REFp (pos++);

    Ltxt += 21;

    return NO_NEW_VALUE;
}

// ---------------------------------------------------------
// INFO/HP
// Homopolymer run length around variant locus
// ---------------------------------------------------------

#define HP_MAX_SPAN 20 // empirically observed maximum value

static int64_t HP_prediction (ConstRangeP range, PosType64 pos)
{
    decl_acgt_decode;
    int64_t hp_prediction[2] = {};

    char bases[2] = { REFp(pos - 1), REFp(pos + 1) };

    for (int side=0; side < 2; side++) {
        for (int i=1; REFp(pos - i) == bases[side] && i <= HP_MAX_SPAN; i++) hp_prediction[side]++;
        for (int i=1; REFp(pos + i) == bases[side] && i <= HP_MAX_SPAN; i++) hp_prediction[side]++;
    }

    return MAX_(hp_prediction[0], hp_prediction[1]);
} 

void vcf_seg_playpus_INFO_HP (VBlockVCFP vb, ContextP ctx, STRp(hp_str))
{
    RefLock lock = REFLOCK_NONE;
    int64_t hp;

    if (!str_get_int (STRa(hp_str), &hp)) {
        seg_by_ctx (VB, STRa(hp_str), ctx, hp_str_len);
        return;
    }

    if (!flag.reference || segconf.running) goto fallback;

    PosType64 pos = DATA_LINE(vb->line_i)->pos[0]; 

    RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos - HP_MAX_SPAN, 2*HP_MAX_SPAN+1, WORD_INDEX_NONE, 
                                      (IS_REF_EXT_STORE ? &lock : NULL));
    
    if (!range || pos - HP_MAX_SPAN < range->first_pos || pos + HP_MAX_SPAN > range->last_pos)
        goto fallback;

    if (hp == HP_prediction (range, pos)) {
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_PLATYPUS_HP }, 2, ctx, hp_str_len);

        if (IS_REF_EXT_STORE)
            bits_set_region (&range->is_set, pos - HP_MAX_SPAN, HP_MAX_SPAN * 2 + 1);
    }

    else fallback:
        seg_integer (VB, ctx, hp, true, hp_str_len);

    ref_unlock (gref, &lock); 
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_PLATYPUS_HP)
{    
    ConstRangeP range = ref_piz_get_range (vb, gref, HARD_FAIL);

    new_value->i = HP_prediction (range, CTX(VCF_POS)->last_value.i);

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ---------------------------------------------------------
// FORMAT/GOF - Goodness of fit value
// INFO/MGOF  - Worst goodness-of-fit value reported across all samples
// ---------------------------------------------------------

void vcf_seg_platypus_FORMAT_GOF (VBlockVCFP vb, ContextP ctx, STRp(gof_str))
{
    int64_t gof;

    if (ctx_has_value_in_line_(VB, CTX(INFO_MGOF)) && str_get_int (STRa(gof_str), &gof) &&
        gof == CTX(INFO_MGOF)->last_value.i)

        seg_delta_vs_other_do (VB, ctx, CTX(INFO_MGOF), 0, 0, gof, -1, gof_str_len); // always delta 0. TO DO: more effecient to use SNIP_COPY

    else 
        seg_integer_or_not (VB, ctx, STRa(gof_str), gof_str_len);
}

// ---------------------------------------------------------
// INFO/WS - Starting position of calling window
// INFO/WE - End position of calling window
// ---------------------------------------------------------

void vcf_seg_playpus_INFO_WS_WE (VBlockVCFP vb, ContextP ctx, STRp(value))
{
    ContextP other_ctx = CTX(ctx->did_i == INFO_WS ? INFO_WE : INFO_WS);

    if (ctx_has_value_in_line_(VB, other_ctx))
        seg_delta_vs_other_do (VB, ctx, other_ctx, STRa(value), 0, 128, value_len);
    else
        seg_delta_vs_other_do (VB, ctx, CTX(VCF_POS), STRa(value), 0, 128, value_len);
}

// ---------------------------------------------------------
// INFO/TC  - Total coverage at this locus
// INFO/TCF - Total forward strand coverage at this locus
// INFO/TCR - Total reverse strand coverage at this locus
// ---------------------------------------------------------

void vcf_seg_playpus_INFO_TCR (VBlockVCFP vb, ContextP ctx, STRp(tcr_str))
{
    int64_t tcr;

    if (ctx_has_value_in_line_(VB, CTX(INFO_TC)) && ctx_has_value_in_line_(VB, CTX(INFO_TCF))
        && str_get_int (STRa(tcr_str), &tcr) && 
        (tcr == CTX(INFO_TC)->last_value.i - CTX(INFO_TCF)->last_value.i))

        seg_by_ctx (VB, STRa(tcr_snip), ctx, tcr_str_len);


    else 
        seg_integer_or_not (VB, ctx, STRa(tcr_str), tcr_str_len);
}
