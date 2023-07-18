// ------------------------------------------------------------------
//   vcf_vagrent.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "piz.h"
#include "reconstruct.h"

// ##INFO=<ID=mutation,Number=1,Type=String,Description="Somatic mutation definition">
// SNP: REF/ALT=C T         mutation=C>T
// INS: REF/ALT=T TTTTTTGT  mutation=->TTTTTGT
// DEL: REF/ALT=TA T        mutation=A>-
// other: REF/ALT=TA AC     mutation=TA>AC
void vcf_seg_INFO_mutation (VBlockVCFP vb, ContextP ctx, STRp(mut))
{
    thool keep_delins_anchor = unknown;
  
    if (vb->main_alt_len == 1 && vb->main_alt[0] == '.' && 
        mut_len == vb->main_ref_len * 2 + 1 && mut[vb->main_ref_len] == '>' &&
        !memcmp (mut, vb->main_ref, vb->main_ref_len) &&
        !memcmp (mut + vb->main_ref_len + 1, vb->main_ref, vb->main_ref_len))
        keep_delins_anchor = no; // ALT='.'
  
    else 
    if (mut_len == 3 && mut[1] == '>' && vb->main_ref_len == 1 && vb->main_alt_len == 1 && 
        vb->main_ref[0] == mut[0] && vb->main_alt[0] == mut[2])
        keep_delins_anchor = no; // SNP

    else 
    if (vb->main_ref_len == 1 && vb->main_alt_len > 1 && 
        vb->main_ref[0] == vb->main_alt[0] && // left anchored
        mut_len == vb->main_alt_len + 1 && mut[0] == '-' && mut[1] == '>' &&
        !memcmp (mut+2, vb->main_alt + 1, vb->main_alt_len - 1))
        keep_delins_anchor = no; // Insertion

    else 
    if (vb->main_ref_len > 1 && vb->main_alt_len == 1 &&
        vb->main_ref[0] == vb->main_alt[0] && // left anchored
        mut_len == vb->main_ref_len + 1 && mut[mut_len-1] == '-' && mut[mut_len-2] == '>' &&
        !memcmp (mut, vb->main_ref + 1, vb->main_ref_len - 1))
        keep_delins_anchor = no; // Deletion
        
    else 
    if (vb->main_ref_len + 1 + vb->main_alt_len == mut_len && mut[vb->main_ref_len] == '>' &&
        !memcmp (mut, vb->main_ref, vb->main_ref_len) &&
        !memcmp (mut + vb->main_ref_len + 1, vb->main_alt, vb->main_alt_len))
        keep_delins_anchor = (vb->main_ref[0] == vb->main_alt[0]); // Other

    if (keep_delins_anchor != unknown)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_mutation, '0' + keep_delins_anchor }, 3, ctx, mut_len);
    else {
        seg_by_ctx (VB, STRa(mut), ctx, mut_len);
}
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_mutation)
{
    STRlast (refalt, VCF_REFALT);
    bool keep_delins_anchor = snip[0] - '0';

    // case: SNP
    // case: ALT=.
    if (!keep_delins_anchor && refalt[refalt_len-1] == '.' && refalt[refalt_len-2] == '\t') {
        RECONSTRUCT (refalt, refalt_len-2);
        RECONSTRUCT1 ('>');
        RECONSTRUCT (refalt, refalt_len-2);
    }

    else if (!keep_delins_anchor && refalt_len == 3) {
        RECONSTRUCT1 (refalt[0]);
        RECONSTRUCT1 ('>');
        RECONSTRUCT1 (refalt[2]);
    }

    // case: Insertion (left-anchored)
    else if (!keep_delins_anchor && refalt[1] == '\t' && refalt[0] == refalt[2]) {
        RECONSTRUCT ("->", 2);
        RECONSTRUCT (refalt + 3, refalt_len - 3);
    }

    // case: Deletion (left-anchored)
    else if (!keep_delins_anchor && refalt[refalt_len-2] == '\t' && refalt[0] == refalt[refalt_len-1]) {
        RECONSTRUCT (refalt + 1, refalt_len - 3);
        RECONSTRUCT (">-", 2);
    }

    // case: Other
    else {
        char *start = BAFTtxt;
        RECONSTRUCT (refalt, refalt_len);
        str_replace_letter (start, refalt_len, '\t', '>');
    }

    return NO_NEW_VALUE;
}
