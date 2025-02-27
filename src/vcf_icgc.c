// ------------------------------------------------------------------
//   vcf_vagrent.c
//   Copyright (C) 2022-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"

// ##INFO=<ID=mutation,Number=1,Type=String,Description="Somatic mutation definition">
// SNP: REF/ALT=C T         mutation=C>T
// INS: REF/ALT=T TTTTTTGT  mutation=->TTTTTGT
// DEL: REF/ALT=TA T        mutation=A>-
// other: REF/ALT=TA AC     mutation=TA>AC
void vcf_seg_INFO_mutation (VBlockVCFP vb, ContextP ctx, STRp(mut))
{
    thool keep_delins_anchor = unknown;
  
    if (IS_PERIOD(vb->ALT) && 
        mut_len == vb->REF_len * 2 + 1 && mut[vb->REF_len] == '>' &&
        !memcmp (mut, vb->REF, vb->REF_len) &&
        !memcmp (mut + vb->REF_len + 1, vb->REF, vb->REF_len))
        keep_delins_anchor = no; // ALT='.'
  
    else 
    if (mut_len == 3 && mut[1] == '>' && VT0(SNP) && 
        vb->REF[0] == mut[0] && vb->ALT[0] == mut[2])
        keep_delins_anchor = no; // SNP

    else 
    if (VT0(INS) && 
        mut_len == vb->ALT_len + 1 && mut[0] == '-' && mut[1] == '>' &&
        !memcmp (mut+2, vb->ALT + 1, vb->ALT_len - 1))
        keep_delins_anchor = no; // Insertion

    else 
    if (VT0(DEL) &&
        mut_len == vb->REF_len + 1 && mut[mut_len-1] == '-' && mut[mut_len-2] == '>' &&
        !memcmp (mut, vb->REF + 1, vb->REF_len - 1))
        keep_delins_anchor = no; // Deletion
        
    else 
    if (vb->REF_len + 1 + vb->ALT_len == mut_len && mut[vb->REF_len] == '>' &&
        !memcmp (mut, vb->REF, vb->REF_len) &&
        !memcmp (mut + vb->REF_len + 1, vb->ALT, vb->ALT_len))
        keep_delins_anchor = (vb->REF[0] == vb->ALT[0]); // Other

    if (keep_delins_anchor != unknown)
        seg_special1 (VB, VCF_SPECIAL_mutation, '0' + keep_delins_anchor, ctx, mut_len);
    else {
        seg_by_ctx (VB, STRa(mut), ctx, mut_len);
}
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_mutation)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    bool keep_delins_anchor = snip[0] - '0';

    // case: SNP
    // case: ALT=.
    // if (!keep_delins_anchor && refalt[refalt_len-1] == '.' && refalt[refalt_len-2] == '\t') {
    if (!keep_delins_anchor && VT0(NO_ALT)) {
        RECONSTRUCT_str (vb->REF);
        RECONSTRUCT1 ('>');
        RECONSTRUCT_str (vb->REF);
    }

    else if (!keep_delins_anchor && VT0(SNP)) {
        RECONSTRUCT1 (*vb->REF);
        RECONSTRUCT1 ('>');
        RECONSTRUCT1 (*ALTi(0)->alt);
    }

    // case: Insertion (left-anchored)
    else if (!keep_delins_anchor && VT0(INS)) {
        RECONSTRUCT ("->", 2);
        RECONSTRUCT (ALTi(0)->alt + 1, ALTi(0)->alt_len - 1);
    }

    // case: Deletion (left-anchored)
    else if (!keep_delins_anchor && VT0(DEL)) {
        RECONSTRUCT (vb->REF+1, vb->REF_len-1);
        RECONSTRUCT (">-", 2);
    }

    // case: Other
    else {
        RECONSTRUCT_str (vb->REF);
        RECONSTRUCT1 ('>');
        RECONSTRUCT_str (ALTi(0)->alt);
    }

    return NO_NEW_VALUE;
}
