// ------------------------------------------------------------------
//   sam_fields.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "sam_private.h"
#include "seg.h"
#include "codec.h"
#include "piz.h"
#include "reconstruct.h"
#include "context.h"
#include "container.h"
#include "chrom.h"
#include "optimize.h"
#include "strings.h"
#include "segconf.h"

static const StoreType aux_field_store_flag[256] = {
    ['c']=STORE_INT, ['C']=STORE_INT, 
    ['s']=STORE_INT, ['S']=STORE_INT,
    ['i']=STORE_INT, ['I']=STORE_INT,
    ['f']=STORE_FLOAT
};

static const LocalType aux_field_to_ltype[256] = {
    ['c']=LT_INT8,   ['C']=LT_UINT8, 
    ['s']=LT_INT16,  ['S']=LT_UINT16,
    ['i']=LT_INT32,  ['I']=LT_UINT32,
    ['f']=LT_FLOAT32
};

const char aux_sep_by_type[2][256] = { { // compressing from SAM
    ['c']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['C']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
    ['s']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['S']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['i']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['I']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['f']=CI0_NATIVE_NEXT | CI0_TRANS_NOR,                                        // compressing SAM - a float is stored as text, and when piz with translate to BAM - is not reconstructed, instead - translated
    ['Z']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, ['H']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, // reconstruct text and then \t separator if SAM and \0 if BAM 
    ['A']=CI0_NATIVE_NEXT,                                                        // reconstruct character and then \t separator if SAM and no separator for BAM
    ['B']=CI0_NATIVE_NEXT                                                         // reconstruct array and then \t separator if SAM and no separator for BAM
}, 
{ // compressing from BAM
    ['c']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['C']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // reconstruct number and \t separator is SAM, and don't reconstruct anything if BAM (reconstruction will be done by translator)
    ['s']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['S']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['i']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, ['I']=CI0_NATIVE_NEXT | CI0_TRANS_NOR, // -"-
    ['f']=CI0_NATIVE_NEXT,                                                        // compressing SAM - a float is stored as a SPECIAL, and the special reconstructor handles the SAM and BAM reconstructing
    ['Z']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, ['H']=CI0_NATIVE_NEXT | CI0_TRANS_NUL, // reconstruct text and then \t separator if SAM and \0 if BAM 
    ['A']=CI0_NATIVE_NEXT,                                                        // reconstruct character and then \t separator if SAM and no separator for BAM
    ['B']=CI0_NATIVE_NEXT                                                         // reconstruct array and then \t separator if SAM and no separator for BAM
} };

#define DICT_ID_ARRAY(dict_id) (DictId){ .id = { (dict_id).id[0], (dict_id).id[1], '_','A','R','R','A','Y' } } // DTYPE_2

// gets integer field value - may appear before or after this field. 
// return false if not found or not integer or out of [min_value,max_value]
bool sam_seg_peek_int_field (VBlockSAMP vb, Did did_i, int16_t idx, int32_t min_value, int32_t max_value, 
                             bool set_last_value, // set last value if not already set (only if returning true)
                             int32_t *value)      // optional out
{
    decl_ctx (did_i);
    int32_t my_value;

    if (idx == -1)
        return false; // this tag doesn't appear in this alignment

    else if (ctx_has_value_in_line_(VB, ctx)) {
        if (ctx->last_value.i < min_value || ctx->last_value.i > max_value) 
            return false; // note: we test 64 bit variable before truncating to 32 bit

        else {
            if (value) *value = ctx->last_value.i;
            return true;
        }
    }

    else if (idx != -1 && sam_seg_get_aux_int (vb, idx, value ? value : &my_value, IS_BAM_ZIP, min_value, max_value, SOFT_FAIL)) {
        if (set_last_value)
            ctx_set_last_value (VB, ctx, (int64_t)(value ? *value : my_value));
        return true;
    }

    else
        return false;
}

//--------------
// FLAG
//--------------

SamFlagStr sam_dis_FLAG (SamFlags f)
{
    SamFlagStr result;
    char *next = result.s;

    if (f.multi_segs) next += sprintf (next, "MultiSeg,");
    if (f.is_aligned)     next += sprintf (next, "Aligned,");
    if (f.unmapped)       next += sprintf (next, "Unmapped,");
    if (f.next_unmapped)  next += sprintf (next, "NextUnmapped,");
    if (f.rev_comp)       next += sprintf (next, "Revcomp,");
    if (f.next_rev_comp)  next += sprintf (next, "NextRevcomp,");
    if (f.is_first)       next += sprintf (next, "First,");
    if (f.is_last)        next += sprintf (next, "Last,");
    if (f.secondary)      next += sprintf (next, "Secondary,");
    if (f.filtered)       next += sprintf (next, "Filtered,");
    if (f.duplicate)      next += sprintf (next, "Duplicate,");
    if (f.supplementary)  next += sprintf (next, "Supplementary,");
    
    if (next != result.s) *(next-1) = 0; // removal comma separator from final item

    return result;
}

void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAM *dl, unsigned add_bytes)
{    
    // note: for simplicity, we only mux in MAIN.
    // note: we can only mux FLAG by buddy, not mate, bc sam_piz_special_DEMUX_BY_MATE needs FLAG to determine mate
    // note: in collated files, we're better off without this.
    bool do_mux = (sam_is_main_vb && segconf.is_paired && !segconf.is_collated);
    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, SAM_FLAG, (MultiplexerP)&vb->mux_FLAG, sam_has_mate || sam_has_saggy) 
                                  : CTX(SAM_FLAG);

    // case: PRIM line: 
    if (sam_is_prim_vb) {
        if (!IS_SAG_SA) goto normal; // NH-based SA Groups: FLAG will be used in preprocessing to set revcomp, and then also reconstructed normally

        // case: SAG_SA: we store FLAG.rev_comp in OPTION_SA_STRAND instead of FLAG
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~SAM_FLAG_REV_COMP, add_bytes);

        ContextP sa_strand_ctx = CTX(OPTION_SA_STRAND);
        seg_by_ctx (VB, dl->FLAG.rev_comp ? "-" : "+", 1, sa_strand_ctx, 0); 

        // count FLAG field contribution to OPTION_SA_CIGAR, so sam_stats_reallocate can allocate the z_data between CIGAR and SA:Z
        sa_strand_ctx->counts.count++; // contributed z_data due to a single-byte + or -
    }

    // case: DEPN line with SA Group: we know some flags from SA Groups, so we store FLAG without them (reduces FLAG entropy)
    else if (sam_is_depn_vb && sam_seg_has_sag_by_nonSA(vb)) // First, Last, Multi are the same for the group (see sam_sa_seg_depn_find_sagroup_noSA)
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~(SAM_FLAG_MULTI_SEG | SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST), add_bytes);

    else if (sam_is_depn_vb && sam_seg_has_sag_by_SA(vb)) // in SA, we can get the revcomp from the SA alignment (see sam_sa_seg_depn_find_sagroup_SAtag)
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~(SAM_FLAG_MULTI_SEG | SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST | SAM_FLAG_REV_COMP), add_bytes);
    
    // case: depn line in main VB
    // tested this - if snip is not identical to mate's case, the the compression worsens, if snip
    // is identical, we will need to figure out if this is depn from another source - such of a H in CIGAR
    // not doing - too much effort for too small benefit

    // case: we can retrieve the FLAG from this line's mate
    #define SAME_AS_MATE (SAM_FLAG_MULTI_SEG | SAM_FLAG_IS_ALIGNED | SAM_FLAG_SECONDARY | SAM_FLAG_FILTERED | \
                          SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY)
    else if (sam_has_mate && !sam_is_prim_vb && // TO DO: improve sam_load_groups_add_flags to allow FLAG mating in prim VB (bug 620)
        ({ ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); 
           (dl->FLAG.value & SAME_AS_MATE) == (mate_dl->FLAG.value & SAME_AS_MATE) &&
            dl->FLAG.unmapped              == mate_dl->FLAG.next_unmapped          &&
            dl->FLAG.next_unmapped         == mate_dl->FLAG.unmapped               &&
            dl->FLAG.rev_comp              == mate_dl->FLAG.next_rev_comp          &&
            dl->FLAG.next_rev_comp         == mate_dl->FLAG.rev_comp               &&
            dl->FLAG.is_first              == mate_dl->FLAG.is_last                &&
            dl->FLAG.is_last               == mate_dl->FLAG.is_first; }))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_MATE_FLAG }, 2, channel_ctx, add_bytes); // added 12.0.41

    // case: normal snip
    else normal:
        seg_integer_as_text_do (VB, channel_ctx, dl->FLAG.value, add_bytes);

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_FLAG.snip), SAM_FLAG, 0); // de-multiplexor

    // first pairing test (another test is in sam_seg_finalize_segconf)
    if (segconf.running) {
        if (dl->FLAG.is_last && !dl->FLAG.is_first)
            segconf.is_paired = true;
        
        if (sam_is_depn (dl->FLAG))
            segconf.sam_has_depn = true;
    }
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_MATE_FLAG)
{
    if (!VER(14)) sam_piz_set_buddy_v13 (vb); 

    SamFlags flag = { .value = history64(SAM_FLAG, VB_SAM->mate_line_i) }; 

    // flip fields to opposite mate
    SWAPbit (flag.unmapped, flag.next_unmapped);
    SWAPbit (flag.rev_comp, flag.next_rev_comp);
    SWAPbit (flag.is_first, flag.is_last);

    new_value->i = flag.value;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE; 
}

// -------------------------------------------------------------------------------------------------------------------------------------------
// U2:Z "Phred probability of the 2nd call being wrong conditional on the best being wrong" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// -------------------------------------------------------------------------------------------------------------------------------------------

// callback function for compress to get data of one line
COMPRESSOR_CALLBACK (sam_zip_U2)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);

    *line_data_len = MIN_(maximum_size, dl->U2.len);

    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->U2.index);

    if (flag.optimize_QUAL)
        optimize_phred_quality_string (*line_data, *line_data_len);

    if (is_rev) *is_rev = dl->FLAG.rev_comp;
}

// ---------
// BD and BI
// ---------

static void sam_seg_BD_BI_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(field), DictId dict_id, unsigned add_bytes)
{
    bool is_bi = (dict_id.num == _OPTION_BI_Z);
    Context *this_ctx  = is_bi ? CTX(OPTION_BI_Z) : CTX (OPTION_BD_Z);

    if (field_len != dl->SEQ.len) {
        seg_by_ctx (VB, STRa(field), this_ctx, field_len);
        return;
    }
    
    dl->BD_BI[is_bi] = TXTWORD (field);

    CTX(OPTION_BD_BI)->txt_len += add_bytes; 
    CTX(OPTION_BD_BI)->local.len32 += field_len;

    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, SAM_SPECIAL_BDBI }), 2, this_ctx, 0);
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
COMPRESSOR_CALLBACK (sam_zip_BD_BI)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    rom bd = dl->BD_BI[0].index ? Btxt (dl->BD_BI[0].index) : NULL;
    rom bi = dl->BD_BI[1].index ? Btxt (dl->BD_BI[1].index) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    ASSERT (bd && bi, "%s/%u: A line has one of the BD:Z/BI:Z pair - Genozip can only compress lines that have either both BD:Z and BI:Z or neither", VB_NAME, vb_line_i); 
    
    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->SEQ.len * 2);
    if (is_rev) *is_rev = dl->FLAG.rev_comp;

    if (!line_data) return; // only length was requested

    buf_alloc_exact (vb, VB_SAM->interlaced, dl->SEQ.len * 2, uint8_t, "interlaced");

    uint8_t *next = B1ST8 (VB_SAM->interlaced);
    for (uint32_t i=0; i < dl->SEQ.len; i++) {
        *next++ = bd[i];
        *next++ = bi[i] - bd[i];
    }

    *line_data = B1STc (VB_SAM->interlaced);
}   

// BD and BI - reconstruct from BD_BI context which contains interlaced BD and BI data. 
SPECIAL_RECONSTRUCTOR (sam_piz_special_BD_BI)
{
    if (!vb->seq_len || !reconstruct) goto done;

    Context *bdbi_ctx = CTX(OPTION_BD_BI);

    // note: bd and bi use their own next_local to retrieve data from bdbi_ctx. the actual index
    // in bdbi_ctx.local is calculated given the interlacing
    ASSPIZ (ctx->next_local + vb->seq_len * 2 <= bdbi_ctx->local.len, "Error reading: unexpected end of %s data. Expecting ctx->next_local=%u + vb->seq_len=%u * 2 <= bdbi_ctx->local.len=%"PRIu64, 
            dis_dict_id (bdbi_ctx->dict_id).s, ctx->next_local, vb->seq_len, bdbi_ctx->local.len);

    char *dst        = BAFTtxt;
    rom src          = Bc (bdbi_ctx->local, ctx->next_local * 2);
    uint32_t seq_len = vb->seq_len; // automatic var for effeciency

    if (ctx->dict_id.num == _OPTION_BD_Z)
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src;
    else
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src + *(src+1);
    
    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;

done:
    return NO_NEW_VALUE;
}

// ----------------------------
// NM:i "Number of differences"
// ----------------------------

static inline void sam_seg_NM_get_prediction (VBlockSAMP vb, ZipDataLineSAM *dl, 
                                              int32_t *predicted_by_SEQ, int32_t *predicted_by_MD)
{
    *predicted_by_SEQ = (vb->mismatch_bases_by_SEQ != -1          && !vb->cigar_missing && !vb->seq_missing && !dl->FLAG.unmapped) ? 
        vb->mismatch_bases_by_SEQ + vb->deletions + vb->insertions : -1;
    
    *predicted_by_MD  = (has_MD && vb->mismatch_bases_by_MD != -1 && !vb->cigar_missing && !vb->seq_missing && !dl->FLAG.unmapped) ? 
        vb->mismatch_bases_by_MD  + vb->deletions + vb->insertions : -1;
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
void sam_seg_NM_i (VBlockSAMP vb, ZipDataLineSAM *dl, SamNMType nm, unsigned add_bytes)
{
    START_TIMER;

    ContextP ctx = CTX (OPTION_NM_i);

    // check for evidence that NM is integer (not binary) - usually, this happens in segconf, but it can also happen later in the execution
    if (nm > 1 && !segconf.NM_is_integer) 
        segconf.NM_is_integer = true; 

    // make a copy, so it remains consistent throughout the if statements below, even if segconf.NM_is_integer is set by another thread in their midst
    bool NM_is_integer = segconf.NM_is_integer; 

    if (segconf.running) {
        if (has_MD && vb->idx_MD_Z > vb->idx_NM_i) segconf.NM_after_MD = false; // we found evidence that sometimes NM is before MD
        goto no_special;
    }

    // possible already segged - from sam_seg_SA_Z
    if (ctx_has_value_in_line_(vb, ctx)) return;

    ctx_set_last_value (VB, ctx, (int64_t)nm);

    int32_t predicted_by_SEQ, predicted_by_MD;
    sam_seg_NM_get_prediction (vb, dl, &predicted_by_SEQ, &predicted_by_MD);

    if (nm < 0) goto no_special; // invalid nm value

    // method 1: if we have MD:Z, we use prediction of number of mismatches derived by analyzing it. This is almost always correct, 
    // but the downside is that reconstruction takes longer due to the need to peek MD:Z. Therefore, we limit it to certain cases.
    else if (NM_is_integer && nm == predicted_by_MD && 
               (segconf.NM_after_MD            || // case 1: MD is reconstructed before NM so peek is fast
                IS_REF_INTERNAL || // case 2: prediction against SEQ performs poorly
                predicted_by_SEQ != nm         || // case 3: rare cases in which prediction by SEQ is wrong with an external reference.
                flag.best))                       // case 4: the user request the best method
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'm'}, 3, OPTION_NM_i, add_bytes);  // 'm' type since v14

    // method 2: copy from SA Group. DEPN or PRIM line. Note: in DEPN, nm already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    else if (sam_seg_has_sag_by_SA (vb)) 
        sam_seg_against_sa_group (vb, ctx, add_bytes); 

    // method 3: use prediction of the number of mismatches derived by comparing SEQ to a reference.
    // this is usually, but surprisingly not always, correct for an external reference, and often correct for an internal one.
    else if (NM_is_integer && nm == predicted_by_SEQ)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'i'}, 3, OPTION_NM_i, add_bytes); 

    // case nm is a binary 0/1 rather than an integer. We use prediction against SEQ. TO DO: Support prediction against MD:Z too.
    else if (!NM_is_integer && predicted_by_SEQ != -1 && (nm > 0) == (predicted_by_SEQ > 0)) 
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'b'}, 3, OPTION_NM_i, add_bytes); 

    else no_special: 
        seg_integer (VB, ctx, nm, true, add_bytes);

    // in PRIM with SA, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
    if (sam_is_prim_vb && sam_seg_has_sag_by_SA (vb)) {
        seg_integer_as_text (VB, OPTION_SA_NM, nm, 0);  // note: for PRIM lines without SA:Z and nm:i, we seg "0" into OPTION_SA_NM in sam_seg_sag_stuff

        // count NM field contribution to OPTION_SA_NM, so sam_stats_reallocate can allocate the z_data between NM and SA:Z
        CTX(OPTION_SA_NM)->counts.count += add_bytes; 
    }

    COPY_TIMER (sam_seg_NM_i);
}

// --------------------------------------------------------------------------------------------------------------
// nM:i: (STAR) the number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate.
// --------------------------------------------------------------------------------------------------------------

void sam_seg_nM_i (VBlockSAMP vb, ZipDataLineSAM *dl, SamNMType nm, unsigned add_bytes)
{
    ContextP ctx = CTX (OPTION_nM_i);

    // case: in paired files, its expected to be the same value as the mate
    if (segconf.is_paired && !sam_is_depn_vb) 
        sam_seg_buddied_i_fields (vb, dl, OPTION_nM_i, nm, &dl->nM, (MultiplexerP)&vb->mux_nM, STRa(copy_mate_nM_snip), add_bytes);

    // else: in non-paired files, we use the same method as NM:i
    else {
        if (segconf.running) {
            if (has_MD && !ctx_encountered_in_line (VB, OPTION_MD_Z)) segconf.nM_after_MD = false; // we found evidence that sometimes nM is before MD
            goto fallback;
        }

        int32_t predicted_by_SEQ, predicted_by_MD;
        sam_seg_NM_get_prediction (vb, dl, &predicted_by_SEQ, &predicted_by_MD);

        if (nm < 0) goto fallback; // invalid nm value

        // method 1: if we have MD:Z, we use prediction of number of mismatches derived by analyzing it. This is almost always correct, 
        // but the downside is that reconstruction takes longer due to the need to peek MD:Z. Therefore, we limit it to certain cases.
        else if (nm == predicted_by_MD && 
                    (segconf.nM_after_MD            || // case 1: MD is reconstructed before nM so peek is fast
                     IS_REF_INTERNAL || // case 2: prediction against SEQ performs poorly
                     predicted_by_SEQ != nm         || // case 3: rare cases in which prediction by SEQ is wrong with an external reference.
                     flag.best))                       // case 4: the user request the best method
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'm'}, 3, OPTION_nM_i, add_bytes);  // 'm' type since v14

        // method 3: use prediction of the number of mismatches derived by comparing SEQ to a reference.
        // this is usually, but surprisingly not always, correct for an external reference, and often correct for an internal one.
        else if (nm == predicted_by_SEQ)
            seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'i'}, 3, OPTION_nM_i, add_bytes); 

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

SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_NM)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    switch (*snip) {
        case 'm' : new_value->i = sam_piz_NM_get_mismatches_by_MD(vb) + vb->deletions + vb->insertions; break;
        case 'i' : new_value->i = vb->mismatch_bases_by_SEQ + vb->deletions + vb->insertions;           break;
        case 'b' : new_value->i = (vb->mismatch_bases_by_SEQ + vb->deletions + vb->insertions) > 0;     break;
        default  : ASSPIZ (false, "unrecognized opcode '%c'", *snip);
    }

    if (reconstruct) // will be false if BAM, reconstruction is done by translator based on new_value set here
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// SM:i: Template-independent mapping quality
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_SM_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t SM, unsigned add_bytes)
{
    decl_ctx (OPTION_SM_i);

    if (SM >= 0 && SM <= 255 && 
        SM != 254 &&           // note: 254 is a valid, but highly improbable value - we use 254 for "copy from MAPQ" so a actual 254 is segged as an exception
        !(SM && !dl->MAPQ)) {  // we're expecting XM=0 if MAPQ=0
        
        // store value in local (254 means "copy from MAPQ"), except if MAPQ=0 - we know SM=0 so no need to store
        if (dl->MAPQ) {
            uint8_t SM8 = (SM == dl->MAPQ) ? 254 : SM;
            seg_integer_fixed (VB, ctx, &SM8, false, 0);
        }

        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SM }, 2, ctx, add_bytes); // this usually is an all-the-same
    }

    else 
        seg_integer_as_text_do (VB, ctx, SM, add_bytes); 

    ctx_set_last_value (VB, ctx, SM);
}    

SPECIAL_RECONSTRUCTOR (sam_piz_special_SM)
{
    uint8_t MAPQ = CTX(SAM_MAPQ)->last_value.i;

    if (!MAPQ)
        new_value->i = 0;

    else {
        uint8_t value = NEXTLOCAL(uint8_t, ctx);
        new_value->i = (value==254) ? MAPQ : value;
    }

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// AM:i The smallest template-independent mapping quality in the template
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_AM_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t AM, unsigned add_bytes)
{
    decl_ctx (OPTION_AM_i);

    // note: currently we only support for this algorithm AM appearing after SM. Easily fixable if ever needed.
    // AM is often one of 3 options: 0, =SM =MAPQ-SM. If SM=0 then AM is expected to be 0.
    if (has(SM_i) && 
        AM >= 0 && AM <= 255 &&   // valid value
        AM != 253 && AM != 254) { // note: 253,254 are valid, but highly improbable values

        int32_t SM;
        if (!sam_seg_peek_int_field (vb, OPTION_SM_i, vb->idx_SM_i, 0, 255, false, &SM)) goto fallback;

        if (!SM && AM) goto fallback; // usually, AM=0 if SM=0
        
        if (SM) { // no need to store if SM=0, as we know AM=0
            uint8_t AM8 = (AM == SM)                  ? 253 // copy SM (note: AM is not 0, since SM is not 0)
                        : (AM && AM == dl->MAPQ - SM) ? 254 // note subtelty: substraction in uint8_t arithmetics. we are careful to do the same in sam_piz_special_AM.
                        :                               AM;

            seg_integer_fixed (VB, ctx, &AM8, false, 0);
        }
        
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_AM }, 2, ctx, add_bytes);
    }

    else fallback:
        seg_integer_as_text_do (VB, ctx, AM, add_bytes);    
}    

SPECIAL_RECONSTRUCTOR (sam_piz_special_AM)
{
    uint8_t MAPQ = CTX(SAM_MAPQ)->last_value.i;

    uint8_t SM = reconstruct_peek (vb, CTX(OPTION_SM_i), 0, 0).i;

    if (!SM)
        new_value->i = 0;

    else {
        uint8_t value = NEXTLOCAL(uint8_t, ctx);
        new_value->i = (value==253) ? SM 
                     : (value==254) ? MAPQ - SM
                     :                value;
    }

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// UQ:i Phred likelihood of the segment, conditional on the mapping being correct
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_UQ_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t UQ, unsigned add_bytes)
{
    int32_t other; // either AS or NM

    // in Novoalign, usually UQ==AS
    if (MP(NOVOALIGN) && sam_seg_peek_int_field (vb, OPTION_AS_i, vb->idx_AS_i, 0, 1000, true, &other) && other == UQ)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_UQ, '1' }, 3, OPTION_UQ_i, add_bytes);

    // In GATK produced data, in many cases (~half) UQ==2*NM
    else if (!MP(NOVOALIGN) && sam_seg_peek_int_field (vb, OPTION_NM_i, vb->idx_NM_i, 0, 0x3fffffff, 
             false, // bc seg_NM_isam_seg_NM_i inteprets would interpret a set value as "segged already"
             &other) && other*2 == UQ)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_UQ, '2' }, 3, OPTION_UQ_i, add_bytes);

    else 
        seg_integer (VB, CTX(OPTION_UQ_i), UQ, true, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_UQ)
{
    ASSPIZ0 (*snip=='1' || *snip=='2', "Unrecognized snip - upgrade to the most recent version of Genozip");

    int64_t other = reconstruct_peek (vb, CTX(*snip=='1' ? OPTION_AS_i : OPTION_NM_i), 0, 0).i; 

    new_value->i = other * (*snip - '0');

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// NH:i Number of reported alignments that contain the query in the current record
// ----------------------------------------------------------------------------------------------------------

static inline void sam_seg_NH_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t nh, unsigned add_bytes)
{
    decl_ctx (OPTION_NH_i);

    dl->NH = nh; 
    ctx_set_encountered (VB, ctx);  // needed by sam_seg_HI_i

    // build SAG structure in VB, to be later ingested into z_file->sag_*
    if (sam_is_depn_vb && vb->sag && vb->sag->num_alns == nh) {
        sam_seg_against_sa_group (vb, ctx, add_bytes);
        add_bytes = 0; // if this is prim, we seg again below, but we already added the bytes
    }

    else
        sam_seg_buddied_i_fields (vb, dl, OPTION_NH_i, nh, &dl->NH, (MultiplexerP)&vb->mux_NH, 
                                  STRa (copy_buddy_NH_snip), add_bytes);
}

// ----------------------------------------------------------------------------------------------------------
// HI:i Query hit index (a number [1,NH])
// ----------------------------------------------------------------------------------------------------------

SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_BY_BUDDY_MAP)
{
    int channel_i = (last_flags.unmapped?2 : (sam_has_mate || sam_has_saggy)?1 : 0);

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

static inline void sam_seg_HI_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t hi, unsigned add_bytes)
{
    decl_ctx (OPTION_HI_i);

    if (segconf.is_collated) {
        int64_t last = ctx->last_value.i;
        int prediction = (ctx_encountered_in_line (VB, OPTION_NH_i) && dl->NH==0) ? 0         // NH=0 -> HI=0
                       : (!vb->line_i || last == (dl-1)->NH)                      ? 1         // first in new SA Group 
                       :                                                            last + 1; // next HI in current SA Group

        if (prediction == hi)
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_HI }, 2, ctx, add_bytes);

        else  
            seg_integer (VB, ctx, hi, true, add_bytes);
    }
    else  // not collated
        seg_integer (VB, ctx, hi, false, add_bytes);

    ctx_set_last_value (VB, ctx, hi);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_HI)
{
    ContextP nh_ctx = CTX(OPTION_NH_i);
    int64_t last_nh = nh_ctx->last_value.i;
    int64_t last_hi = ctx->last_value.i;

    new_value->i = (ctx_encountered_in_line (VB, OPTION_NH_i) && last_nh == 0)            ? 0 // NH=0 -> HI=0
                 : (!vb->line_i || last_hi == *B(int64_t, nh_ctx->history, vb->line_i-1)) ? 1 // first in new SA Group
                 :                                                                          last_hi + 1; // next HI in current SA Group

    if (reconstruct)
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------
// CP:i Leftmost coordinate of the next hit
// ----------------------------------------

static inline void sam_seg_CP_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t cp, unsigned add_bytes)
{
    if (sam_is_depn_vb && vb->sag && IS_SAG_CC && vb->cc_aln->pos == cp) 
        sam_seg_against_sa_group (vb, CTX(OPTION_CP_i), add_bytes);

    // if there's no gencomp - CP tends to compress well by delta, as subsequent lines also tend to have subsequent secondaries
    else if (sam_is_main_vb && segconf.is_sorted)
        seg_pos_field (VB, OPTION_CP_i, OPTION_CP_i, 0, 0, 0, 0, cp, add_bytes); 

    else {
        uint32_t cp32 = cp;
        seg_integer_fixed (VB, CTX(OPTION_CP_i), &cp32, true, add_bytes);
    }
}

// -------------------------
// OA:Z "Original alignment"
// -------------------------

static bool sam_seg_OA_rname_cb (VBlockP vb, ContextP ctx, STRp(oa_rname), uint32_t repeat)
{
    if (str_issame(oa_rname, vb->chrom_name))
        seg_by_ctx (vb, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_OA_RNAME }, 2, ctx, oa_rname_len); 

    else
        chrom_seg_ex (vb, OPTION_OA_RNAME, STRa(oa_rname), 0, NULL, oa_rname_len, true, NULL);
    
    return true; // segged successfully
}

// copy from RNAME (note: we always recon textual, unlike RNAME that is reconstructed binary in BAM)
SPECIAL_RECONSTRUCTOR (sam_piz_special_OA_RNAME)
{
    if (!reconstruct) goto done;

    ContextP rname_ctx = CTX(SAM_RNAME);
    ASSPIZ0 (ctx_has_value_in_line_(vb, rname_ctx), "No value for RNAME in line");

    ctx_get_snip_by_word_index (rname_ctx, rname_ctx->last_value.i, snip);
    
    RECONSTRUCT_snip;

done:
    return NO_NEW_VALUE;
}

static bool sam_seg_OA_pos_cb (VBlockP vb, ContextP ctx, STRp(oa_pos), uint32_t repeat)
{
    seg_pos_field (vb, OPTION_OA_POS, SAM_POS, 0, 0, STRa(oa_pos), 0, oa_pos_len);
    return true; // segged successfully
}

// OA is: (rname, pos, strand, CIGAR, mapQ, NM ;)+ . NM is optional (but its , is not)
// Example OA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
// See: https://samtools.github.io/hts-specs/SAMtags.pdf
static void sam_seg_OA_Z (VBlockSAMP vb, STRp(field), unsigned add_bytes)
{
    static const MediumContainer container_OA = { .nitems_lo = 6,          
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_OA_RNAME  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_POS    }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_STRAND }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_CIGAR  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_MAPQ   }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_NM     },                    } } };

    SegCallback callbacks[6] = { [SA_RNAME]=sam_seg_OA_rname_cb, [SA_POS]=sam_seg_OA_pos_cb, [SA_CIGAR]=sam_seg_0A_cigar_cb, [SA_MAPQ]=sam_seg_0A_mapq_cb };
     
    seg_array_of_struct (VB, CTX(OPTION_OA_Z), container_OA, STRa(field), callbacks, add_bytes);
}

// ----------------------------------------------------------------------------------------------
// AS:i "Alignment score generated by aligner" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// ----------------------------------------------------------------------------------------------

// AS has a value set (at least as set by BWA and IonTorrent TMAP) of at most vb->ref_consumed, and often equal to it. we modify
// it to be new_value=(value-ref_consumed) 
static inline void sam_seg_AS_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t as, unsigned add_bytes)
{
    START_TIMER;

    if (segconf.running) 
        if (ABS((int32_t)2*vb->ref_consumed - as) < 10) segconf.AS_is_2ref_consumed++;
    
    ctx_set_last_value (VB, CTX (OPTION_AS_i), as);

    // depn VB - seg against prim line AS (sag_has_AS determines if its beneficial to do so)
    if (sam_is_depn_vb && segconf.sag_has_AS && vb->sag) 
        sam_seg_against_sa_group_int (vb, CTX(OPTION_AS_i), (int64_t)vb->sag->as - as, add_bytes);

    // in bowtie2-like data, we might be able to copy from mate
    else if (segconf.is_bowtie2) {
        ASSERT (as >= MIN_AS_i && as <= MAX_AS_i, "%s: AS=%"PRId64" is out of range [%d,%d]", LN_NAME, as, MIN_AS_i, MAX_AS_i);    
        
        ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_AS_i, (MultiplexerP)&vb->mux_AS, sam_has_mate);

        if (sam_has_mate && segconf.has[OPTION_YS_i] && mate_dl->YS == as) 
            seg_by_ctx (VB, STRa(copy_mate_YS_snip), channel_ctx, add_bytes);
        
        // case: AS:i tends to be close to 2 X ref_consumed
        else if (segconf.AS_is_2ref_consumed && !segconf.running && ABS((int32_t)2*vb->ref_consumed - as) < 10) {
            SNIPi3 (SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED, 'x', (int32_t)2*vb->ref_consumed - as);
            seg_by_ctx (VB, STRa(snip), channel_ctx, add_bytes); 
        }

        else
            // TODO: AS prediction, see bug 520
            //seg_integer_as_text_do (VB, channel_ctx, as, add_bytes);    
            seg_integer (VB, channel_ctx, as, true, add_bytes);

        seg_by_did (VB, STRa(vb->mux_AS.snip), OPTION_AS_i, 0); // de-multiplexor
    }

    // case: in STAR paired files, AS expected to be the same value as its mate's
    else if (MP(STAR) && segconf.is_paired && !sam_is_depn_vb) 
        sam_seg_buddied_i_fields (vb, dl, OPTION_AS_i, as, &dl->AS, (MultiplexerP)&vb->mux_AS, STRa(copy_mate_AS_snip), add_bytes);

    // not bowtie2: store a special snip with delta from ref_consumed
    else if (ABS((int32_t)vb->ref_consumed-as) < 10) {
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED, (int32_t)vb->ref_consumed-as);
        seg_by_did (VB, STRa(snip), OPTION_AS_i, add_bytes); 
    }

    // if we have minimap2-produced ms:i, AS is close to it
    else if (segconf.sam_ms_type == ms_MINIMAP2 && sam_seg_peek_int_field (vb, OPTION_ms_i, vb->idx_ms_i, -0x8000000, 0x7fffffff, true/*needed for delta*/, NULL)) 
        seg_delta_vs_other (VB, CTX(OPTION_AS_i), CTX(OPTION_ms_i), NULL, add_bytes);

    else
        seg_integer (VB, CTX (OPTION_AS_i), as, true, add_bytes);

    dl->AS = as;

    COPY_TIMER(sam_seg_AS_i);
}

// reconstruct seq_len or (seq_len-snip)
// Note: This is used by AS:i fields in files compressed up to 12.0.37
SPECIAL_RECONSTRUCTOR (sam_piz_special_AS_old)
{
    new_value->i = (int32_t)vb->seq_len - atoi (snip); // seq_len if snip=""
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// reconstruct ref_consumed or (ref_consumed-snip)
SPECIAL_RECONSTRUCTOR (sam_piz_special_REF_CONSUMED)
{
    if (snip_len && *snip == 'x')
        new_value->i = 2 * (int32_t)VB_SAM->ref_consumed - atoi (snip+1); // snip="" is 0
    else
        new_value->i =     (int32_t)VB_SAM->ref_consumed - atoi (snip); 

    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// Prediction based on SEQ.len and soft_clip. Use for:
// qs:i (pacvio) : per-read: 0-based start of query in the ZMW read (absent in CCS)
// qe:i (pacbio) : per-read: 0-based end of query in the ZMW read (absent in CCS)
// XS:i (blasr)  : read alignment start position without counting previous soft clips (1 based)
// XE:i (blasr)  : read alignment end position without counting previous soft clips (1 based)
// XL:i (blasr)  : aligned read length
// XQ:i (blasr)  : query read length
// QS:i (ngmlr)  : query start
// QE:i (ngmlr)  : query end
// XR:i (ngmlr)  : query length minus soft clips
// ----------------------------------------------------------------------------------------------

static inline int32_t sam_SEQ_END_prediction (rom op, int32_t seq_len, int32_t soft0, int32_t soft1) // all signed
{
    return (op[0]=='+' ? seq_len : 0)
         + (op[1]=='+' ? 1       : 0)
         + (op[2]=='+' ? soft0   : 0)    
         - (op[2]=='-' ? soft0   : 0)    
         + (op[3]=='+' ? soft1   : 0)    
         - (op[3]=='-' ? soft1   : 0);
} 

static inline void sam_seg_SEQ_END (VBlockSAMP vb, ZipDataLineSAM *dl, ContextP ctx, int64_t value, rom op, unsigned add_bytes)
{
    int32_t prediction = sam_SEQ_END_prediction (op, dl->SEQ.len, vb->soft_clip[0], vb->soft_clip[1]);

    if (value == prediction)
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SEQ_LEN, op[0], op[1], op[2], op[3] }, 6, ctx, add_bytes);

    else 
        seg_integer (VB, ctx, value, true, add_bytes);
}

// reconstruct a modified dl->SEQ.len based on a 4-parameter snip
SPECIAL_RECONSTRUCTOR (sam_piz_special_SEQ_LEN)
{
    new_value->i = sam_SEQ_END_prediction (snip, CTX(SAM_SQBITMAP)->last_txt.len, VB_SAM->soft_clip[0], VB_SAM->soft_clip[1]);

    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// MAPQ and MQ:i (mate MAPQ)
// ----------------------------------------------------------------------------------------------

// We seg against a previous mate line's MQ if one exists, but not if this is a single-MAPQ-value file
void sam_seg_MAPQ (VBlockSAMP vb, ZipDataLineSAM *dl, unsigned add_bytes)
{
    if (segconf.running && dl->MAPQ) {
        if (!segconf.MAPQ_value) 
            segconf.MAPQ_value = dl->MAPQ;
        else if (segconf.MAPQ_value != dl->MAPQ) 
            segconf.MAPQ_has_single_value = false;
    }

    ctx_set_last_value (VB, CTX(SAM_MAPQ), (int64_t)dl->MAPQ);

    bool do_mux = sam_is_main_vb && segconf.is_paired; // for simplicity. To do: also for prim/depn components
    int channel_i = sam_has_mate?1 : sam_has_prim?2 : 0;
    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, SAM_MAPQ, (MultiplexerP)&vb->mux_MAPQ, channel_i) 
                                  : CTX(SAM_MAPQ);

    ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

    if (sam_seg_has_sag_by_SA (vb)) {
        sam_seg_against_sa_group (vb, channel_ctx, add_bytes);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (sam_is_prim_vb) {
            ContextP sa_mapq_ctx = CTX(OPTION_SA_MAPQ);
            seg_integer_fixed (VB, sa_mapq_ctx, &dl->MAPQ, false, 0);

            // count MAPQ field contribution to OPTION_SA_MAPQ, so sam_stats_reallocate can allocate the z_data between MAPQ and SA:Z
            sa_mapq_ctx->counts.count += add_bytes; 
        }
    }

    // case: seg against mate
    else if (!segconf.MAPQ_has_single_value && segconf.has[OPTION_MQ_i] &&
             sam_has_mate && dl->MAPQ == mate_dl->MQ)
        seg_by_ctx (VB, STRa(copy_mate_MQ_snip), channel_ctx, add_bytes); // copy MQ from earlier-line mate 

    // // case: seg against predicted alignment in prim line SA:Z 
    // else if (!segconf.MAPQ_has_single_value && zip_has_prim && 
    //          dl->MAPQ && dl->MAPQ != (dl-1)->MAPQ &&
    //          sam_seg_is_item_predicted_by_prim_SA (vb, SA_MAPQ, dl->MAPQ)) 
    //     seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_PRIM, '0'+SA_MAPQ }, 3, channel_ctx, add_bytes); 

    else {
        channel_ctx->ltype = LT_UINT8;
        seg_integer_fixed (VB, channel_ctx, &dl->MAPQ, true, add_bytes);
    }

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_MAPQ.snip), SAM_MAPQ, 0); // de-multiplexor
}

// MQ:i Mapping quality of the mate/next segment
// Seg against mate if we have one, or else against MAPQ as it is often very similar
static inline void sam_seg_MQ_i (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t mq, unsigned add_bytes)
{
    ASSERT (mq >=0 && mq <= 255, "%s: Invalid MQ:i=%"PRId64": expecting an integer [0,255]", LN_NAME, mq);
    dl->MQ = mq; 
    
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_MQ_i, (MultiplexerP)&vb->mux_MQ, sam_has_mate);

    if (sam_has_mate && 
        dl->MQ == DATA_LINE (vb->mate_line_i)->MAPQ) 
        seg_by_ctx (VB, STRa(copy_mate_MAPQ_snip), channel_ctx, add_bytes); // copy MAPQ from earlier-line mate 

    else 
        seg_delta_vs_other_do (VB, channel_ctx, CTX(SAM_MAPQ), 0, 0, mq, -1, add_bytes);

    seg_by_did (VB, STRa(vb->mux_MQ.snip), OPTION_MQ_i, 0); // de-multiplexor
}

// RG:Z, PG:Z, PU:Z, LB:Z, RX:Z and others: we predict that the value will be the same as the buddy
void sam_seg_buddied_Z_fields (VBlockSAMP vb, ZipDataLineSAM *dl, MatedZFields f, STRp(value), 
                               SegBuddiedCallback seg_cb, // optional
                               unsigned add_bytes)
{
    ContextP ctx = CTX(buddied_Z_dids[f]);

    ZipDataLineSAM *buddy_dl = sam_has_mate  ? DATA_LINE (vb->mate_line_i)  // mate before saggy, consistent with sam_piz_special_COPY_BUDDY
                             : sam_has_saggy ? DATA_LINE (vb->saggy_line_i) // if collated, it will likely be the same as prev line, so better not use special 
                             : NULL;

    bool do_mux = segconf.is_paired && !segconf.is_collated; // muxing just on prim if not paired is counter-beneficial

    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, buddied_Z_dids[f], (MultiplexerP)&vb->mux_mated_z_fields[f], !!buddy_dl)
                                  : ctx;

    if (do_mux && buddy_dl && str_issame_(STRa(value), STRtxtw (buddy_dl->mated_z_fields[f]))) 
        seg_by_ctx (VB, STRi(copy_buddy_Z_snip,f), channel_ctx, add_bytes); 

    else if (seg_cb && !segconf.running) 
        seg_cb (vb, channel_ctx, STRa(value), add_bytes);

    else 
        seg_by_ctx (VB, STRa(value), channel_ctx, add_bytes); 

    if (do_mux)
        seg_by_did (VB, STRa(vb->mux_mated_z_fields[f].snip), buddied_Z_dids[f], 0); // de-multiplexor

    dl->mated_z_fields[f] = TXTWORD(value);
}

void sam_seg_buddied_i_fields (VBlockSAMP vb, ZipDataLineSAM *dl, Did did_i, 
                               int64_t my_value, 
                               int32_t *dl_value, // pointer to my dl value (value will be set)
                               MultiplexerP mux,
                               STRp(copy_snip),   // buddy_type embedded in snip needs to be consistent with mux->special_code
                               unsigned add_bytes)
{
    decl_ctx (did_i);

    ASSERT (ctx->flags.store_per_line || ctx->flags.spl_custom || segconf.running,  
            "%s: expecting ctx=%s to have store_per_line=true or spl_custom=true", LN_NAME, ctx->tag_name);

    // BAM spec permits values up to 0xffffffff, and SAM is unlimited, however for code covenience we limit
    // values segged with this method to int32_t. If this is ever an issue, it can be solved.
    ASSERT (my_value >= -0x80000000LL && my_value <= 0x7fffffffLL, "%s: Value of %s is %"PRId64", outside the supported range by Genozip of [%d,%d]",
            LN_NAME, ctx->tag_name, my_value, -0x80000000, 0x7fffffff);

    #define by_mate      (mux->special_code == SAM_SPECIAL_DEMUX_BY_MATE)
    #define by_buddy     (mux->special_code == SAM_SPECIAL_DEMUX_BY_BUDDY)
    #define by_buddy_map (mux->special_code == SAM_SPECIAL_DEMUX_BY_BUDDY_MAP)
    #define by_mate_prim (mux->special_code == SAM_SPECIAL_DEMUX_BY_MATE_PRIM)

    if (!segconf.running && segconf.is_paired && !segconf.is_collated) {
        int channel_i = by_mate      ? sam_has_mate
                      : by_buddy     ? (sam_has_mate || sam_has_saggy)
                      : by_buddy_map ? (dl->FLAG.unmapped?2 : (sam_has_mate || sam_has_saggy)?1 : 0)
                      : by_mate_prim ? (sam_has_mate?1 : sam_has_prim?2 : 0)
                      : -1; // invalid

        ASSERT (channel_i >= 0, "mux of %s has unsupported special code %u", CTX(did_i)->tag_name, mux->special_code);

        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, did_i, mux, channel_i);

        // calculating mate value's pointer. Note: this only works if ZipDataLineSAM is NOT packed
        int32_t *mate_value  = &dl_value[(vb->mate_line_i  - vb->line_i) * sizeof (ZipDataLineSAM) / sizeof (int32_t)];
        int32_t *saggy_value = &dl_value[(vb->saggy_line_i - vb->line_i) * sizeof (ZipDataLineSAM) / sizeof (int32_t)];

        // note: we don't check if our mate has this field. the math still works, because if it doesn't have 
        // this field, the dl value in Seg will be 0, as will be the History value in recon.

        if (by_buddy_map && dl->FLAG.unmapped)
            seg_integer_as_text_do (VB, channel_ctx, my_value, add_bytes); // seg as snip, as this will likely become all-the-same

        else if (sam_has_mate && my_value == *mate_value)    
            seg_by_ctx (VB, STRa(copy_snip), channel_ctx, add_bytes); // note: prior to v14, we stored ms qual_score in QUAL history, not in ms:i history 
        
        else if (((by_buddy && sam_has_saggy) || (by_mate_prim && sam_has_prim)) && my_value == *saggy_value)
            seg_by_ctx (VB, STRa(copy_snip), channel_ctx, add_bytes); 

        else 
            seg_integer (VB, channel_ctx, my_value, true, add_bytes);    

        seg_by_ctx (VB, MUX_SNIP(mux), MUX_SNIP_LEN(mux), ctx, 0); // de-multiplexor
    }
    else
        seg_integer (VB, ctx, my_value, true, add_bytes);        
}

// E2 - SEQ data. Currently broken. To do: fix (bug 403)
/*static void sam_seg_E2_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->SEQ.len, field, "E2 tag without a SEQ"); 
    ASSINP (field_len == dl->SEQ.len, 
            "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
            txt_name, dl->SEQ.len, field_len, field_len, field);

    PosType32 this_pos = vb->last_int(SAM_POS);

    sam_seg_SEQ (vb, OPTION_E2_Z, (char *)STRa(field), this_pos, vb->last_cigar, vb->ref_consumed, vb->ref_and_seq_consumed, 0, field_len, // remove const bc SEQ data is actually going to be modified
                        vb->last_cigar, add_bytes); 
}*/

// U2 - QUAL data
static void sam_seg_U2_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->SEQ.len, field, "U2 tag without a SEQ"); 
    ASSINP (field_len == dl->SEQ.len, 
            "%s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. U2=%.*s",
            LN_NAME, dl->SEQ.len, field_len, STRf(field));

    dl->U2 = TXTWORD (field);
    CTX(OPTION_U2_Z)->txt_len   += add_bytes;
    CTX(OPTION_U2_Z)->local.len += field_len;
}

static inline unsigned sam_seg_aux_add_bytes (char type, unsigned value_len, bool is_bam)
{
    if (!is_bam || type=='Z' || type=='H')
        return value_len + 1; // +1 for \0 (BAM Z/H) or \t (SAM)

    else
        return aux_width[(uint8_t)type]; // BAM

    // note: this will return 0 for 'B'
}

static inline SmallContainer *sam_seg_array_field_get_con (VBlockSAMP vb, Context *con_ctx, uint8_t type, bool has_callback,
                                                           ContextP *elem_ctx) // out
{
    // case: cached with correct type
    if (con_ctx->con_cache.param == type) {
        SmallContainer *con = B1ST (SmallContainer, con_ctx->con_cache);
        *elem_ctx = ctx_get_ctx (vb, con->items[1].dict_id);
        return con; // note: we return a pointer into con_cache- sam_seg_array_field may modify the number of repeats. that's ok.
    }

    buf_alloc (vb, &con_ctx->con_cache, 0, 1, SmallContainer, 0, "con_cache");

    SmallContainer *con = B1ST (SmallContainer, con_ctx->con_cache);
    *con = (SmallContainer){ .nitems_lo = 2, 
                             .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repsep and obsoleting this flag 
                             .repsep    = {0,0}, 
                             .items     = { { .translator = SAM2BAM_ARRAY_SELF  },  // item[0] is translator-only item - to translate the Container itself in case of reconstructing BAM 
                                            { .separator  = {0, ','}            } } // item[1] is actual array item
                           };            
    
    // prepare context where array elements will go in
    con->items[1].dict_id      = DICT_ID_ARRAY (con_ctx->dict_id);
    con->items[1].translator   = aux_field_translator (type); // instructions on how to transform array items if reconstructing as BAM (array[0] is the subtype of the array)
    con->items[1].separator[0] = aux_sep_by_type[IS_BAM_ZIP][type];
    
    *elem_ctx = ctx_get_ctx (vb, con->items[1].dict_id);
    (*elem_ctx)->st_did_i = con_ctx->did_i;
    
    StoreType store_type = aux_field_store_flag[type];
    (*elem_ctx)->flags.store = store_type;

    ASSERT (store_type, "%s: Invalid type \"%c\" in array of %s", LN_NAME, type, con_ctx->tag_name);

    if (store_type == STORE_INT) {
        (*elem_ctx)->ltype = aux_field_to_ltype[type];
        (*elem_ctx)->local_is_lten = true; // we store in local in LTEN (as in BAM) and *not* in machine endianity
    }

    con_ctx->con_cache.param = type;

    return con;
}

// an array - all elements go into a single item context, multiple repeats. items are segged as dynamic integers or floats, or a callback is called to seg them.
typedef void (*ArrayItemCallback) (VBlockSAMP vb, ContextP ctx, void *cb_param, void *array, uint32_t array_len);
static void sam_seg_array_field (VBlockSAMP vb, DictId dict_id, uint8_t type, 
                                 rom array, int/*signed*/ array_len, // SAM: comma separated array ; BAM : arrays original width and machine endianity
                                 ArrayItemCallback callback, void *cb_param) // optional - call back for each item to seg the item
{   
    // prepare array container - a single item, with number of repeats of array element. array type is stored as a prefix
    Context *con_ctx = ctx_get_ctx (vb, dict_id), *elem_ctx;
    SmallContainer *con = sam_seg_array_field_get_con (vb, con_ctx, type, !!callback, &elem_ctx);

    int width = aux_width[type];
    bool is_bam = IS_BAM_ZIP;

    int array_bytes = is_bam ? (width * array_len) : array_len;
    elem_ctx->txt_len += array_bytes;

    con->repeats = is_bam ? array_len : (1 + str_count_char (STRa(array), ','));
    ASSERT (con->repeats < CONTAINER_MAX_REPEATS, "%s: array has too many elements, more than %u", LN_NAME, CONTAINER_MAX_REPEATS);

    bool is_int = (elem_ctx->flags.store == STORE_INT);
    if (is_int) {
        elem_ctx->local.len *= width; // len will be calculated in bytes in this function
        buf_alloc (vb, &elem_ctx->local, con->repeats * width, con->repeats * width * 50, char, CTX_GROWTH, "contexts->local"); // careful not * line.len - we can get OOM
    }

    ASSERT (is_int || (type=='f' && !callback), "%s: Type not supported for SAM/BAM arrays '%c'(%u) in ctx=%s",
            LN_NAME, type, type, con_ctx->tag_name);

    char *local_start = BAFTc (elem_ctx->local);

    if (is_bam) {        

        if (is_int) 
            buf_add (&elem_ctx->local, array, array_bytes); // LTEN, not machine endianity (local_is_lten is set)

        else // FLOAT
            // TO DO: we can use SNIP_LOOKUP like seg_float_or_not and store the data in local (bug 500)
            for (uint32_t i=0; i < array_len; i++, array += width) {
                SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_FLOAT, GET_UINT32(array));
                seg_by_ctx (VB, STRa(snip), elem_ctx, 0); // TODO: seg in local and update SPECIAL to get from local if no snip
            }
    }

    else { // SAM
        // note: we're not using str_split on array, because the number of elements can be very large (eg one per base in PacBio ip:B) - possibly stack overflow
        for (uint32_t i=0; i < con->repeats; i++) { // str_len will be -1 after last number

            rom snip = array;
            for (; array_len && *array != ','; array++, array_len--) {};

            unsigned snip_len = (unsigned)(array - snip);

            if (is_int) { 
                int64_t value;
                ASSERT (str_get_int (STRa(snip), &value), "%s: Invalid array: \"%.*s\", expecting an integer in an array element of %s", 
                        LN_NAME, STRf(snip), con_ctx->tag_name);
                value = LTEN64 (value); // consistent with BAM and avoid a condition before the memcpy below 
                memcpy (&local_start[i*width], &value, width);
            }
            else 
                seg_float_or_not (VB, elem_ctx, STRa(snip), 0);
            
            array_len--; // skip comma
            array++;
        }

        if (is_int) elem_ctx->local.len += con->repeats * width;
    }

    if (callback)
        callback (vb, elem_ctx, cb_param, local_start, con->repeats);

    if (is_int)
       elem_ctx->local.len /= width; // return len back to counting in units of ltype
 
    // add bytes here in case of BAM - all to main field
    unsigned container_add_bytes = is_bam ? (4/*count*/ + 1/*type*/) : (2/*type - eg "i,"*/ + 1/*\t or \n*/);
    container_seg (vb, con_ctx, (ContainerP)con, ((char[]){ CON_PX_SEP, type, ',', CON_PX_SEP }), 4, container_add_bytes);
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the AUX dictionary entry
DictId sam_seg_aux_field (VBlockSAMP vb, ZipDataLineSAM *dl, bool is_bam, 
                          rom tag, char bam_type, char array_subtype, 
                          STRp(value), ValueType numeric) // two options 
{
    char sam_type = sam_seg_bam_type_to_sam_type (bam_type);
    char dict_name[6] = { tag[0], tag[1], ':', sam_type, ':', array_subtype }; // last 2 are ignored if not array
    DictId dict_id = dict_id_make (dict_name, (sam_type=='B' ? 6 : 4), DTYPE_SAM_AUX); // match dict_id as declared in #pragma GENDICT

    unsigned add_bytes = sam_seg_aux_add_bytes (bam_type, value_len, is_bam);

    #define COND0(condition, seg) if (condition) { seg; break; } else  
    #define COND(condition,  seg) if (condition) { seg; break; } else goto fallback; 

    if (segconf.running)
        segconf.has[ctx_get_ctx (VB, dict_id)->did_i]++;

    switch (dict_id.num) {

        // ---------------------
        // Standard fields
        // ---------------------
        case _OPTION_SA_Z: sam_seg_SA_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_OA_Z: sam_seg_OA_Z (vb, STRa(value), add_bytes); break;

        case _OPTION_OC_Z: sam_seg_other_CIGAR (vb, CTX(OPTION_OC_Z), STRa(value), true, add_bytes); break;

        case _OPTION_OP_i: seg_pos_field (VB, OPTION_OP_i, SAM_POS, 0, 0, 0, 0, numeric.i, add_bytes); break;

        case _OPTION_OQ_Z: sam_seg_other_qual (vb, &dl->OQ, OPTION_OQ_Z, STRa(value), true, add_bytes); break;
        case _OPTION_TQ_Z: sam_seg_other_qual (vb, &dl->TQ, OPTION_TQ_Z, STRa(value), false, add_bytes); break;

        //case _OPTION_CC_Z: see unused code, bug 609

        case _OPTION_MC_Z: sam_cigar_seg_MC_Z (vb, dl, STRa(value), add_bytes); break;

        // note: we cannot seg MD using SPECIAL if we're segging the line against SA Group, because
        // the MD alg requires the SQBITMAP.
        case _OPTION_MD_Z: sam_seg_MD_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_NM_i: COND (!MP(BLASR), sam_seg_NM_i (vb, dl, (SamNMType)numeric.i, add_bytes)); // blasr uses NM:i in a non-standard way

        case _OPTION_nM_i: COND (MP(STAR), sam_seg_nM_i (vb, dl, (SamNMType)numeric.i, add_bytes)); break;

        case _OPTION_BD_Z:
        case _OPTION_BI_Z: sam_seg_BD_BI_Z (vb, dl, STRa(value), dict_id, add_bytes); break;
        
        case _OPTION_RG_Z: COND (segconf.sam_multi_RG, sam_seg_buddied_Z_fields (vb, dl, MATED_RG, STRa(value), 0, add_bytes));

        case _OPTION_PG_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_PG, STRa(value), 0, add_bytes); break;
        case _OPTION_PU_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_PU, STRa(value), 0, add_bytes); break;
        case _OPTION_LB_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_LB, STRa(value), 0, add_bytes); break;
        case _OPTION_OX_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_OX, STRa(value), 0, add_bytes); break;
        case _OPTION_MI_Z: sam_seg_buddied_Z_fields (vb, dl, MATED_MI, STRa(value), 0, add_bytes); break;

        case _OPTION_CY_Z: dl->dont_compress_CY = sam_seg_barcode_qual (vb, dl, OPTION_CY_Z, SOLO_CY, segconf.n_CR_CB_CY_seps, STRa(value), qSTRa(segconf.CY_con_snip), &segconf.CY_con, add_bytes); break;
        case _OPTION_QT_Z: dl->dont_compress_QT = sam_seg_barcode_qual (vb, dl, OPTION_QT_Z, SOLO_QT, segconf.n_BC_QT_seps,    STRa(value), qSTRa(segconf.QT_con_snip), &segconf.QT_con, add_bytes); break;

        case _OPTION_CR_Z: sam_seg_CR_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_CB_Z: sam_seg_CB_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_UR_Z: // alias of RX (cellranger)
        case _OPTION_RX_Z: sam_seg_RX_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_UB_Z: // alias of BX (cellranger)
        case _OPTION_BX_Z: sam_seg_BX_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_UY_Z: // alias of QX (cellranger)
        case _OPTION_QX_Z: sam_seg_QX_Z (vb, dl, STRa(value), add_bytes); break;
        case _OPTION_BC_Z: sam_seg_BC_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_2Y_Z: COND (segconf.has_cellranger, sam_seg_other_qual (vb, &dl->_2Y, OPTION_2Y_Z, STRa(value), false, add_bytes));
        
        case _OPTION_GR_Z: COND (segconf.has_cellranger, sam_seg_GR_Z (vb, dl, STRa(value), add_bytes));
        case _OPTION_GY_Z: COND (segconf.has_cellranger, sam_seg_GY_Z (vb, dl, STRa(value), add_bytes));
        case _OPTION_2R_Z: COND (segconf.has_cellranger, sam_seg_other_seq (vb, dl, OPTION_2R_Z, STRa(value), add_bytes));
        case _OPTION_TR_Z: COND (segconf.has_cellranger || MP(LONGRANGER), sam_seg_other_seq (vb, dl, OPTION_TR_Z, STRa(value), add_bytes));

        case _OPTION_TX_Z: COND (segconf.has_cellranger, sam_seg_TX_AN_Z (vb, dl, OPTION_TX_Z, STRa(value), add_bytes));
        case _OPTION_AN_Z: COND (segconf.has_cellranger, sam_seg_TX_AN_Z (vb, dl, OPTION_AN_Z, STRa(value), add_bytes));

        case _OPTION_GN_Z: sam_seg_gene_name_id (vb, dl, OPTION_GN_Z, STRa(value), add_bytes); break;
        case _OPTION_GX_Z: sam_seg_gene_name_id (vb, dl, OPTION_GX_Z, STRa(value), add_bytes); break;
        case _OPTION_gn_Z: sam_seg_gene_name_id (vb, dl, OPTION_gn_Z, STRa(value), add_bytes); break;
        case _OPTION_gx_Z: sam_seg_gene_name_id (vb, dl, OPTION_gx_Z, STRa(value), add_bytes); break;

        
        //case _OPTION_E2: sam_seg_E2_field (vb, dl, STRa(value), add_bytes); // BROKEN. To do: fix.

        case _OPTION_U2_Z: sam_seg_U2_Z (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_AS_i: sam_seg_AS_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_MQ_i: sam_seg_MQ_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_SM_i: sam_seg_SM_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_AM_i: sam_seg_AM_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_HI_i: sam_seg_HI_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_NH_i: sam_seg_NH_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_CP_i: sam_seg_CP_i (vb, dl, numeric.i, add_bytes); break;
        case _OPTION_UQ_i: sam_seg_UQ_i (vb, dl, numeric.i, add_bytes); break;

        // ---------------------
        // Non-standard fields
        // ---------------------
        case _OPTION_XA_Z: COND (segconf.sam_has_BWA_XA_Z != no, sam_seg_BWA_XA_Z (vb, STRa(value), add_bytes)); // yes or unknown
        
        case _OPTION_XS_i: COND0 (segconf.sam_has_BWA_XS_i, sam_seg_BWA_XS_i (vb, dl, OPTION_XS_i, numeric.i, add_bytes))
                           COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XS_i), numeric.i, "0+00", add_bytes)); 

        case _OPTION_XC_i: COND (segconf.sam_has_BWA_XC_i, sam_seg_BWA_XC_i (vb, dl, numeric.i, add_bytes)); 

        case _OPTION_XT_A: COND (segconf.sam_has_BWA_XT_A, sam_seg_BWA_XT_A (vb, value[0], add_bytes));

        case _OPTION_XM_i: COND0 (segconf.sam_has_BWA_XM_i, sam_seg_BWA_XM_i (vb, numeric, add_bytes))
                           COND (MP(TMAP), sam_seg_TMAP_XM_i (vb, numeric, add_bytes));

        case _OPTION_X1_i: COND (segconf.sam_has_BWA_X01_i, sam_seg_BWA_X1_i (vb, numeric.i, add_bytes)); break;

        case _OPTION_ZS_i: COND (MP(HISAT2), sam_seg_BWA_XS_i (vb, dl, OPTION_ZS_i, numeric.i, add_bytes)); 

        case _OPTION_YS_i: COND (sam_has_bowtie2_YS_i(), sam_seg_bowtie2_YS_i (vb, dl, numeric, add_bytes));

        // case _OPTION_XR_Z: COND (MP(BISMARK) || MP(DRAGEN), sam_seg_bismark_XR_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XO_Z: COND (MP(BSSEEKER2), sam_seg_bsseeker2_XO_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XG_Z: COND0 (MP(BISMARK) || MP(DRAGEN), sam_seg_bismark_XG_Z (vb, dl, STRa(value), add_bytes))
                           COND (MP(BSSEEKER2), sam_seg_bsseeker2_XG_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XM_Z: COND0 (MP(BISMARK) || MP(DRAGEN), sam_seg_bismark_XM_Z (vb, dl, OPTION_XM_Z, SAM_SPECIAL_BISMARK_XM, STRa(value), add_bytes))
                           COND (MP(BSSEEKER2), sam_seg_bsseeker2_XM_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_XB_A: COND (MP(GEM3), sam_seg_gem3_XB_A (vb, dl, STRa(value), add_bytes));

        case _OPTION_XB_Z: COND (MP(BSBOLT), sam_seg_bsbolt_XB (vb, dl, STRa(value), add_bytes));

        case _OPTION_YS_Z: COND (MP(BSBOLT), sam_seg_bsbolt_YS_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_mc_i: COND (segconf.is_biobambam2_sort, sam_seg_mc_i (vb, numeric.i, add_bytes));

        case _OPTION_ms_i: COND (segconf.sam_ms_type == ms_BIOBAMBAM && !flag.optimize_QUAL, sam_seg_ms_i (vb, dl, numeric.i, add_bytes)); // ms:i produced by biobambam or samtools

        case _OPTION_s1_i: COND (segconf.is_minimap2, sam_seg_s1_i (vb, dl, numeric.i, add_bytes));

        case _OPTION_cm_i: COND (segconf.is_minimap2, sam_seg_cm_i (vb, dl, numeric.i, add_bytes));

        // tx:i: - we seg this as a primary field SAM_TAX_ID
        case _OPTION_tx_i: seg_by_did (VB, taxid_redirection_snip, taxid_redirection_snip_len, OPTION_tx_i, add_bytes); break;

        case _OPTION_Z5_i: seg_pos_field (VB, OPTION_Z5_i, SAM_PNEXT, 0, 0, 0, 0, numeric.i, add_bytes); break;

        case _OPTION_Zs_Z: COND (MP(HISAT2), sam_seg_HISAT2_Zs_Z (vb, STRa(value), add_bytes));

        case _OPTION_ZM_B_s: COND (MP(TMAP) && flag.optimize_ZM, sam_seg_array_field (vb, _OPTION_ZM_B_s, array_subtype, STRa(value), sam_optimize_TMAP_ZM, 0));

        case _OPTION_YH_Z: COND (MP(NOVOALIGN), seg_add_to_local_text (VB, CTX(OPTION_YH_Z), STRa(value), LOOKUP_NONE, add_bytes)); break;
        case _OPTION_YQ_Z: COND (MP(NOVOALIGN), seg_add_to_local_text (VB, CTX(OPTION_YQ_Z), STRa(value), LOOKUP_NONE, add_bytes)); break;

        case _OPTION_qs_i: COND (TECH(PACBIO), sam_seg_SEQ_END (vb, dl, CTX(OPTION_qs_i), numeric.i, "0+00", add_bytes));
        case _OPTION_qe_i: COND (TECH(PACBIO), sam_seg_SEQ_END (vb, dl, CTX(OPTION_qe_i), numeric.i, "++00", add_bytes));

        case _OPTION_QS_i: COND (MP(NGMLR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_QS_i), numeric.i, "00+0", add_bytes));
        case _OPTION_QE_i: COND (MP(NGMLR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_QE_i), numeric.i, "+00-", add_bytes));
        case _OPTION_XR_i: COND (MP(NGMLR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XR_i), numeric.i, "+0--", add_bytes));
        case _OPTION_XE_i: COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XE_i), numeric.i, "++00", add_bytes));
        case _OPTION_XQ_i: COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XQ_i), numeric.i, "+000", add_bytes));
        case _OPTION_XL_i: COND (MP(BLASR), sam_seg_SEQ_END (vb, dl, CTX(OPTION_XL_i), numeric.i, "+0--", add_bytes));

        case _OPTION_FI_i: COND (MP(BLASR), sam_seg_blasr_FI_i (vb, dl, numeric.i, add_bytes));

        case _OPTION_fx_Z: COND (segconf.has_cellranger, sam_seg_fx_Z (vb, dl, STRa(value), add_bytes));

        case _OPTION_dq_Z: COND0 (segconf.use_pacbio_iqsqdq, sam_seg_pacbio_xq (vb, dl, OPTION_dq_Z, &dl->dq, STRa(value), add_bytes))
                           COND (TECH(PACBIO), seg_add_to_local_text (VB, CTX(OPTION_dq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes)); 
        
        case _OPTION_iq_Z: COND0 (segconf.use_pacbio_iqsqdq, sam_seg_pacbio_xq (vb, dl, OPTION_iq_Z, &dl->iq, STRa(value), add_bytes))
                           COND (TECH(PACBIO), seg_add_to_local_text (VB, CTX(OPTION_iq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes)); 
        
        case _OPTION_sq_Z: COND0 (segconf.use_pacbio_iqsqdq, sam_seg_pacbio_xq (vb, dl, OPTION_sq_Z, &dl->sq, STRa(value), add_bytes))
                           COND (TECH(PACBIO), seg_add_to_local_text (VB, CTX(OPTION_sq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes)); 
        
        case _OPTION_dt_Z: COND (TECH(PACBIO), seg_add_to_local_text (VB, CTX(OPTION_dt_Z), STRa(value), LOOKUP_SIMPLE, add_bytes));
        case _OPTION_mq_Z: COND (TECH(PACBIO), seg_add_to_local_text (VB, CTX(OPTION_mq_Z), STRa(value), LOOKUP_SIMPLE, add_bytes));
        case _OPTION_st_Z: COND (TECH(PACBIO), seg_add_to_local_text (VB, CTX(OPTION_st_Z), STRa(value), LOOKUP_SIMPLE, add_bytes));
        
        default: fallback:
            
            // all types of integer
            if (sam_type == 'i') {
                ContextP ctx = ctx_get_ctx (vb, dict_id);

                if (!ctx->is_initialized) {
                    ctx->ltype = LT_DYN_INT;
                    ctx->flags.store = STORE_INT; // needs this be reconstructable as BAM
                    ctx->is_initialized = true;
                }

                seg_integer (VB, ctx, numeric.i, true, add_bytes);
            }

            else if (sam_type == 'f') {
                ContextP ctx = ctx_get_ctx (vb, dict_id);

                if (!ctx->is_initialized) {
                    ctx->flags.store = STORE_FLOAT; // needs this be reconstructable as BAM
                    ctx->is_initialized = true;
                }

                if (is_bam) {
                    SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_FLOAT, numeric.i);
                    seg_by_ctx (VB, STRa(snip), ctx, add_bytes); 
                }
                else
                    seg_float_or_not (VB, ctx, STRa(value), add_bytes);
            }

            // Numeric array
            else if (sam_type == 'B') 
                sam_seg_array_field (vb, dict_id, array_subtype, STRa(value), NULL, NULL);

            // Z,H,A - normal snips in their own dictionary
            else        
                seg_by_dict_id (VB, STRa(value), dict_id, add_bytes); 
    }
    
    return dict_id;
}
