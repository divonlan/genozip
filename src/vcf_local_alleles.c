// ------------------------------------------------------------------
//   vcf_local_alleles.c
//   Copyright (C) 2024-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

static void vcf_LAA_predict (VBlockVCFP vb, qSTRp (pred))
{
    STRlast(gt, FORMAT_GT);

    if (gt_len > *pred_len) {
        *pred_len = 0;
        return;
    }    
    
    char *out = pred;

    char sep = memchr (gt, '/', gt_len) ? '/' : '|';
    str_split (gt, gt_len, vb->ploidy, sep, ht, false);

    // 0/1 → 0 ; 0/0 → ∅ ; ./. → ∅
    int first_ht = 0;
    while (first_ht < n_hts && ht_lens[first_ht]==1 && (hts[first_ht][0] == '0' || hts[first_ht][0] == '.'))
        first_ht++;

    // ∅ → .
    if (first_ht == n_hts) { 
        pred[0] = '.'; 
        *pred_len = 1;
    }

    // remove redundant alleles: 23/23 → 23 ; 2/500/2 → 2,500
    else {
        // copy alleles one by one to out, skipping 2nd+ occurance of the same allele
        for (int i=first_ht; i < n_hts; i++) {
            for (int j=first_ht; j < i; j++)
                if (str_issame_(STRi(ht,i), STRi(ht,j))) 
                    goto skip_i; // case: the allele on ht=i was already encountered in an earlier ht=j

            memcpy (out, hts[i], ht_lens[i]);
            out += ht_lens[i];
            *out++ = ',';
            skip_i: {}
        }
        out--; // remove final ,

        *pred_len = out - pred;
    }
}

void vcf_seg_FORMAT_LAA (VBlockVCFP vb, ContextP ctx, STRp(laa))
{
    STRlic(prediction, 100);
    vcf_LAA_predict (vb, qSTRa (prediction));

    if (str_issame (laa, prediction))
        seg_special0 (VB, VCF_SPECIAL_LAA, ctx, laa_len);
    else
        seg_by_ctx (VB, STRa(laa), ctx, laa_len);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_LAA)
{
    if (reconstruct) {
        STRlic(prediction, 100);
        vcf_LAA_predict (VB_VCF, qSTRa (prediction));

        RECONSTRUCT_str (prediction);
    }

    return NO_NEW_VALUE;
}

