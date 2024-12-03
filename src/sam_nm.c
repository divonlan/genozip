// ------------------------------------------------------------------
//   sam_nm.c
//   Copyright (C) 2021-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"

// ----------------------------
// NM:i "Number of differences"
// ----------------------------

static inline void sam_seg_NM_get_prediction (VBlockSAMP vb, ZipDataLineSAMP dl, bool with_indels,
                                              int32_t *predicted_by_SEQ, int32_t *predicted_by_MD)
{
    int32_t indels_count = with_indels ? (vb->deletions + vb->insertions) : 0;

    *predicted_by_SEQ = (vb->mismatch_bases_by_SEQ != -1          && !vb->cigar_missing && !vb->seq_missing && !dl->FLAG.unmapped) ? 
        vb->mismatch_bases_by_SEQ + indels_count : -1;
    
    *predicted_by_MD  = (has_MD && vb->mismatch_bases_by_MD != -1 && !vb->cigar_missing && !vb->seq_missing && !dl->FLAG.unmapped) ? 
        vb->mismatch_bases_by_MD  + indels_count : -1;
}

// Two variations:
// 1) Integer NM per SAM specification https://samtools.github.io/hts-specs/SAMtags.pdf: "Number of differences (mismatches plus inserted and deleted bases) 
// between the sequence and reference, counting only (case-insensitive) A, C, G and T bases in sequence and reference as potential matches, with everything
// else being a mismatch. Note this means that ambiguity codes in both sequence and reference that match each other, such as ‘N’ in both, or compatible 
// codes such as ‘A’ and ‘R’, are still counted as mismatches. The special sequence base ‘=’ will always be considered to be a match, even if the reference 
// is ambiguous at that point. Alignment reference skips, padding, soft and hard clipping (‘N’, ‘P’, ‘S’ and ‘H’ CIGAR operations) do not count as mismatches,
// but insertions and deletions count as one mismatch per base."
// Note: we observed cases (eg PacBio data with bwa-sw) that NM is slightly different than expected, potentially
// seggable with a delta. However, the added entropy to b250 outweighs the benefit, and we're better off without delta.
// 2) Binary NM: 0 if sequence fully matches the reference when aligning according to CIGAR, 1 is not.
void sam_seg_NM_i (VBlockSAMP vb, ZipDataLineSAMP dl, SamNMType nm, unsigned add_bytes)
{
    START_TIMER;

    decl_ctx (OPTION_NM_i);

    // check for evidence that NM is integer (not binary) - usually, this happens in segconf, but it can also happen later in the execution
    if (nm > 1 && !segconf.NM_is_integer) 
        segconf.NM_is_integer = true; 

    // make a copy, so it remains consistent throughout the if statements below, even if segconf.NM_is_integer is set by another thread in their midst
    bool NM_is_integer = segconf.NM_is_integer; 

    if (segconf_running) {
        if (has_MD && vb->idx_MD_Z > vb->idx_NM_i) segconf.NM_after_MD = false; // we found evidence that sometimes NM is before MD
        goto no_special;
    }

    // possible already segged - from sam_seg_SA_Z
    if (ctx_has_value_in_line_(vb, ctx)) return;

    ctx_set_last_value (VB, ctx, (int64_t)nm);

    int32_t predicted_by_SEQ, predicted_by_MD;
    sam_seg_NM_get_prediction (vb, dl, true, &predicted_by_SEQ, &predicted_by_MD);

    if (nm < 0) goto no_special; // invalid nm value

    // method 0: if we use 'X' in CIGAR, NM:i can be derived from the sum of 'X' bases
    else if (segconf.CIGAR_has_eqx && nm == vb->mismatch_bases_by_CIGAR + vb->deletions + vb->insertions) // observed in pbmm2
        seg_special1 (VB, SAM_SPECIAL_NM, 'x', ctx, add_bytes); // 15.0.66 

    // method 1: if we have MD:Z, we use prediction of number of mismatches derived by analyzing it. This is almost always correct, 
    // but the downside is that reconstruction takes longer due to the need to peek MD:Z. Therefore, we limit it to certain cases.
    else if (NM_is_integer && nm == predicted_by_MD && 
               (segconf.NM_after_MD            || // case 1: MD is reconstructed before NM so peek is fast
                IS_REF_INTERNAL                || // case 2: prediction against SEQ performs poorly
                predicted_by_SEQ != nm         || // case 3: rare cases in which prediction by SEQ is wrong with an external reference.
                flag.best))                       // case 4: the user request the best method
        seg_special1 (VB, SAM_SPECIAL_NM, 'm', ctx, add_bytes);  // 'm' type since v14

    // method 2: copy from SA Group. DEPN or PRIM line. Note: in DEPN, nm already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    // note: if SA_NM_by_CIGAR_X=true then NM(=mismatches+insertions+deletions) != SA_NM(=mismatches) 
    else if (sam_seg_has_sag_by_SA (vb) && !segconf.SA_NM_by_CIGAR_X) 
        sam_seg_against_sa_group (vb, ctx, add_bytes); 

    // method 3: use prediction of the number of mismatches derived by comparing SEQ to a reference.
    // this is usually, but surprisingly not always, correct for an external reference, and often correct for an internal one.
    else if (NM_is_integer && nm == predicted_by_SEQ)
        seg_special1 (VB, SAM_SPECIAL_NM, 'i', ctx, add_bytes); 

    // case nm is a binary 0/1 rather than an integer. We use prediction against SEQ. TO DO: Support prediction against MD:Z too.
    else if (!NM_is_integer && predicted_by_SEQ != -1 && (nm > 0) == (predicted_by_SEQ > 0)) 
        seg_special1 (VB, SAM_SPECIAL_NM, 'b', ctx, add_bytes); 

    else no_special: 
        seg_integer (VB, ctx, nm, true, add_bytes);

    // in PRIM with SA, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
    if (IS_PRIM(vb) && !segconf.SA_NM_by_CIGAR_X && sam_seg_has_sag_by_SA (vb)) {
        seg_integer_as_snip (VB, OPTION_SA_NM, nm, 0);  // note: for PRIM lines without SA:Z and nm:i, we seg "0" into OPTION_SA_NM in sam_seg_sag_stuff

        // count NM field contribution to OPTION_SA_NM, so sam_stats_reallocate can allocate the z_data between NM and SA:Z
        CTX(OPTION_SA_NM)->counts.count += add_bytes; 
    }

    COPY_TIMER (sam_seg_NM_i);
}

// -------------------------------------------------
// XM:i "Number of mismatches in the alignment" 
// appears in bwa, dragen, bowtie2, novoalign...
// -------------------------------------------------

void sam_seg_XM_i (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t xm, int16_t idx, unsigned add_bytes)
{
    decl_ctx (OPTION_XM_i);
    int32_t predicted_by_SEQ, predicted_by_MD;
    sam_seg_NM_get_prediction (vb, dl, false, &predicted_by_SEQ, &predicted_by_MD);

    // method 1: if we have MD:Z, we use prediction of number of mismatches derived by analyzing it. This is almost always correct, 
    // but the downside is that reconstruction takes longer due to the need to peek MD:Z. Therefore, we limit it to certain cases.
    if (xm == predicted_by_MD && 
            ((vb->idx_MD_Z >= 0 && idx > vb->idx_MD_Z) || // case 1: MD is reconstructed before XM so peek is fast
            IS_REF_INTERNAL                            || // case 2: prediction against SEQ performs poorly
            predicted_by_SEQ != xm                     || // case 3: rare cases in which prediction by SEQ is wrong with an external reference.
            flag.best))                                   // case 4: the user request the best method
        seg_special1 (VB, SAM_SPECIAL_NM, 'M', ctx, add_bytes);  // 'm' type since v14

    // method 2: use prediction of the number of mismatches derived by comparing SEQ to a reference.
    // this is usually, but surprisingly not always, correct for an external reference, and often correct for an internal one.
    else if (xm == predicted_by_SEQ)
        seg_special1 (VB, SAM_SPECIAL_NM, 'I', ctx, add_bytes); 

    else 
        seg_integer (VB, ctx, xm, true, add_bytes);
}

// --------------------------------------------------------------------------------------------------------------
// nM:i: (STAR) the number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate.
// --------------------------------------------------------------------------------------------------------------

void sam_seg_STAR_nM (VBlockSAMP vb, ZipDataLineSAMP dl, SamNMType nm, unsigned add_bytes)
{
    ContextP ctx = CTX (OPTION_nM_i);

    // case: in paired files, its expected to be the same value as the mate
    if (segconf.is_paired && !IS_DEPN(vb)) 
        sam_seg_buddied_i_fields (vb, dl, OPTION_nM_i, nm, &dl->nM, (MultiplexerP)&vb->mux_nM, STRa(copy_mate_nM_snip), add_bytes);

    // else: in non-paired files, we use the same method as NM:i
    else {
        if (segconf_running) {
            if (has_MD && !ctx_encountered_in_line (VB, OPTION_MD_Z)) segconf.nM_after_MD = false; // we found evidence that sometimes nM is before MD
            goto fallback;
        }

        int32_t predicted_by_SEQ, predicted_by_MD;
        sam_seg_NM_get_prediction (vb, dl, true, &predicted_by_SEQ, &predicted_by_MD);

        if (nm < 0) goto fallback; // invalid nm value

        // method 1: if we have MD:Z, we use prediction of number of mismatches derived by analyzing it. This is almost always correct, 
        // but the downside is that reconstruction takes longer due to the need to peek MD:Z. Therefore, we limit it to certain cases.
        else if (nm == predicted_by_MD && 
                    (segconf.nM_after_MD            || // case 1: MD is reconstructed before nM so peek is fast
                     IS_REF_INTERNAL                || // case 2: prediction against SEQ performs poorly
                     predicted_by_SEQ != nm         || // case 3: rare cases in which prediction by SEQ is wrong with an external reference.
                     flag.best))                       // case 4: the user request the best method
            seg_special1 (VB, SAM_SPECIAL_NM, 'm', ctx, add_bytes);  // 'm' type since v14

        // method 3: use prediction of the number of mismatches derived by comparing SEQ to a reference.
        // this is usually, but surprisingly not always, correct for an external reference, and often correct for an internal one.
        else if (nm == predicted_by_SEQ)
            seg_special1 (VB, SAM_SPECIAL_NM, 'i', ctx, add_bytes); 

        else fallback: 
            seg_integer (VB, ctx, nm, true, add_bytes);
    }
}

static inline int64_t sam_piz_NM_get_mismatches_by_MD (VBlockSAMP vb)
{
    STR(MD);
    reconstruct_peek (VB, CTX(OPTION_MD_Z), pSTRa(MD));

    // count mismatches
    bool deletion=false;
    int32_t mismatches = 0;
    for (int i=0; i < MD_len; i++) 
        switch (MD[i]) {
            case '^'         : deletion = true;  break;
            case '0' ... '9' : deletion = false; break;
            default          : if (!deletion) mismatches++;
        }

    return mismatches;
}

// used for NM:i, nM:i and XM:i
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_NM)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    // note: 'M' and 'I' introduced 15.0.25
    int32_t indels_count = (*snip!='M' && *snip!='I') ? (vb->deletions + vb->insertions) : 0;

    switch (*snip) {
        case 'M' : case 'm' : new_value->i = sam_piz_NM_get_mismatches_by_MD(vb) + indels_count; break;
        case 'I' : case 'i' : new_value->i = vb->mismatch_bases_by_SEQ + indels_count;           break;
        case 'b' :            new_value->i = (vb->mismatch_bases_by_SEQ + indels_count) > 0;     break;
        case 'x' :            new_value->i = vb->mismatch_bases_by_CIGAR + indels_count;         break; // 15.0.66
        default  :            ABORT_PIZ ("unrecognized opcode '%c'. %s", *snip, genozip_update_msg());
    }

    if (reconstruct) // will be false if BAM, reconstruction is done by translator based on new_value set here
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

