// ------------------------------------------------------------------
//   sam_fields.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

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
#include "profiler.h"
#include "lookback.h"
#include "endianness.h"

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

//--------------
// FLAG
//--------------

SamFlagStr sam_dis_FLAG (SamFlags f)
{
    SamFlagStr result;
    char *next = result.s;

    if (f.bits.multi_segments) next += sprintf (next, "MultiSeg,");
    if (f.bits.is_aligned)     next += sprintf (next, "Aligned,");
    if (f.bits.unmapped)       next += sprintf (next, "Unmapped,");
    if (f.bits.next_unmapped)  next += sprintf (next, "NextUnmapped,");
    if (f.bits.rev_comp)       next += sprintf (next, "Revcomp,");
    if (f.bits.next_rev_comp)  next += sprintf (next, "NextRevcomp,");
    if (f.bits.is_first)       next += sprintf (next, "First,");
    if (f.bits.is_last)        next += sprintf (next, "Last,");
    if (f.bits.secondary)      next += sprintf (next, "Secondary,");
    if (f.bits.filtered)       next += sprintf (next, "Filtered,");
    if (f.bits.duplicate)      next += sprintf (next, "Duplicate,");
    if (f.bits.supplementary)  next += sprintf (next, "Supplementary,");
    
    if (next != result.s) *(next-1) = 0; // removal comma separator from final item

    return result;
}

void sam_seg_FLAG (VBlockSAMP vb, ZipDataLineSAM *dl, unsigned add_bytes)
{    
    // note: for simplicity, we only mux in MAIN.
    // note: we can only mux FLAG by buddy, not mate, bc sam_piz_special_DEMUX_BY_MATE needs FLAG to determine mate
    // note: in collated files, we're better off without this.
    bool do_mux = (sam_is_main_vb && segconf.sam_is_paired && !segconf.sam_is_collated);
    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_FLAG, zip_has_mate || zip_has_prim) 
                                  : CTX(SAM_FLAG);

    // case: PRIM line: we store FLAG.rev_comp in OPTION_SA_STRAND instead of FLAG
    if (sam_is_prim_vb) {
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~SAM_FLAG_REV_COMP, add_bytes);

        ContextP sa_strand_ctx = CTX(OPTION_SA_STRAND);
        seg_by_ctx (VB, dl->FLAG.bits.rev_comp ? "-" : "+", 1, sa_strand_ctx, 0); 

        // count FLAG field contribution to OPTION_SA_CIGAR, so sam_stats_reallocate can allocate the z_data between CIGAR and SA:Z
        sa_strand_ctx->counts.count++; // contributed z_data due to a single-byte + or -
    }

    // case: DEPN line with SA Group: we know some flags from SA Groups, so we store FLAG without them (reduces FLAG entropy)
    else if (sam_is_depn_vb && sam_seg_has_SA_Group(vb))
        #define SA_GROUP_FLAGS (SAM_FLAG_MULTI_SEGMENTS | SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST | SAM_FLAG_REV_COMP)
        sam_seg_against_sa_group_int (vb, CTX(SAM_FLAG), dl->FLAG.value & ~SA_GROUP_FLAGS, add_bytes);
    
    // case: depn line in main VB
    // tested this - if snip is not identical to mate's case, the the compression worsens, if snip
    // is identical, we will need to figure out if this is depn from another source - such of a H in CIGAR
    // not doing - too much effort for too small benefit

    // case: we can retrieve the FLAG from this line's mate
    #define SAME_AS_MATE (SAM_FLAG_MULTI_SEGMENTS | SAM_FLAG_IS_ALIGNED | SAM_FLAG_SECONDARY | SAM_FLAG_FILTERED | \
                        SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY)
    else if (zip_has_mate &&
        ({ ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); 
           (dl->FLAG.value & SAME_AS_MATE) == (mate_dl->FLAG.value & SAME_AS_MATE) &&
           dl->FLAG.bits.unmapped          == mate_dl->FLAG.bits.next_unmapped     &&
           dl->FLAG.bits.next_unmapped     == mate_dl->FLAG.bits.unmapped          &&
           dl->FLAG.bits.rev_comp          == mate_dl->FLAG.bits.next_rev_comp     &&
           dl->FLAG.bits.next_rev_comp     == mate_dl->FLAG.bits.rev_comp          &&
           dl->FLAG.bits.is_first          == mate_dl->FLAG.bits.is_last           &&
           dl->FLAG.bits.is_last           == mate_dl->FLAG.bits.is_first; }))
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_BUDDY_FLAG }, 2, channel_ctx, add_bytes); // added 12.0.41

    // case: normal snip
    else 
        seg_integer_as_text_do (VB, channel_ctx, dl->FLAG.value, add_bytes);

    if (do_mux)
        seg_by_did_i (VB, STRa(vb->mux_FLAG.snip), SAM_FLAG, 0); // de-multiplexor

    if (segconf.running && dl->FLAG.bits.is_last && !dl->FLAG.bits.is_first)
        segconf.sam_is_paired = true;
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_BUDDY_FLAG)
{
    if (vb->buddy_line_i == NO_LINE) 
        reconstruct_set_buddy (vb); 

    SamFlags flag = { .value = *B(int64_t, CTX(SAM_FLAG)->history, vb->buddy_line_i) }; 

    // flip fields to opposite mate
    SWAPbit (flag.bits.unmapped, flag.bits.next_unmapped);
    SWAPbit (flag.bits.rev_comp, flag.bits.next_rev_comp);
    SWAPbit (flag.bits.is_first, flag.bits.is_last);

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

    *line_data = Bc (vb->txt_data, dl->U2.index);

    if (flag.optimize_QUAL)
        optimize_phred_quality_string (*line_data, *line_data_len);

    if (is_rev) *is_rev = dl->FLAG.bits.rev_comp;
}

// ---------
// BD and BI
// ---------

static void sam_seg_BD_BI_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(field), DictId dict_id, unsigned add_bytes)
{
    bool is_bi = (dict_id.num == _OPTION_BI_Z);
    Context *this_ctx  = is_bi ? CTX(OPTION_BI_Z) : CTX (OPTION_BD_Z);

    if (field_len != dl->SEQ.len) {
        seg_by_ctx (VB, STRa(field), this_ctx, field_len);
        return;
    }
    
    dl->BD_BI[is_bi] = TXTWORD (field);

    CTX(OPTION_BD_BI)->txt_len += add_bytes; 

    if (!dl->BD_BI[!is_bi].index) // the first of BD and BI increments local.len, so it is incremented even if just one of BD/BI appears
        CTX(OPTION_BD_BI)->local.len += field_len * 2;

    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, SAM_SPECIAL_BDBI }), 2, this_ctx, 0);
}

// callback function for compress to get BD_BI data of one line: this is an
// interlaced line containing a character from BD followed by a character from BI - since these two fields are correlated
COMPRESSOR_CALLBACK (sam_zip_BD_BI)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    rom bd = dl->BD_BI[0].index ? Bc (vb->txt_data, dl->BD_BI[0].index) : NULL;
    rom bi = dl->BD_BI[1].index ? Bc (vb->txt_data, dl->BD_BI[1].index) : NULL;
    
    if (!bd && !bi) return; // no BD or BI on this line

    ASSERT (bd && bi, "%s: A line has one of the BD:Z/BI:Z pair - Genozip can only compress lines that have either both BD:Z and BI:Z or neither", LN_NAME); 
    
    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->SEQ.len * 2);
    if (is_rev) *is_rev = dl->FLAG.bits.rev_comp;

    if (!line_data) return; // only length was requested

    buf_alloc (vb, &VB_SAM->bd_bi_line, 0, dl->SEQ.len * 2, uint8_t, 2, "bd_bi_line");

    // calculate character-wise delta
    for (unsigned i=0; i < dl->SEQ.len; i++) {
        *B8 (VB_SAM->bd_bi_line, i*2    ) = bd[i];
        *B8 (VB_SAM->bd_bi_line, i*2 + 1) = bi[i] - (bd ? bd[i] : 0);
    }

    *line_data = B1STc (VB_SAM->bd_bi_line);
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
void sam_seg_NM_field (VBlockSAMP vb, ZipDataLineSAM *dl, SamNMType NM, unsigned add_bytes)
{
    ContextP ctx = CTX (OPTION_NM_i);

    if (segconf.running) {
        if (NM > 1) segconf.NM_is_integer = true; // we found evidence of integer NM
        if (has_MD && vb->idx_MD_Z > vb->idx_NM_i) segconf.NM_after_MD = false; // we found evidence that sometimes NM is before MD
        goto no_special;
    }

    // possible already segged - from sam_SA_Z_seg
    if (ctx_has_value_in_line_(vb, ctx)) return;

    ctx_set_last_value (VB, ctx, (int64_t)NM);

    int32_t predicted_by_SEQ = (vb->mismatch_bases_by_SEQ != -1          && !vb->cigar_missing && !vb->seq_missing && !dl->FLAG.bits.unmapped) ? 
        vb->mismatch_bases_by_SEQ + vb->deletions + vb->insertions : -1;
    
    int32_t predicted_by_MD  = (has_MD && vb->mismatch_bases_by_MD != -1 && !vb->cigar_missing && !vb->seq_missing && !dl->FLAG.bits.unmapped) ? 
        vb->mismatch_bases_by_MD  + vb->deletions + vb->insertions : -1;

    if (NM < 0) goto no_special; // invalid NM value

    // method 1: if we have MD:Z, we use prediction of number of mismatches derived by analyzing it. This is almost always correct, 
    // but the downside is that reconstruction takes longer due to the need to peek MD:Z. Therefore, we limit it to certain cases.
    else if (segconf.NM_is_integer && NM == predicted_by_MD && 
               (segconf.NM_after_MD            || // case 1: MD is reconstructed before NM so peek is fast
                flag.reference == REF_INTERNAL || // case 2: prediction against SEQ performs poorly
                predicted_by_SEQ != NM         || // case 3: rare cases in which prediction by SEQ is wrong with an external reference.
                flag.best))                       // case 4: the user request the best method
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'm'}, 3, OPTION_NM_i, add_bytes);  // 'm' type since v14

    // method 2: copy from SA Group. DEPN or PRIM line. Note: in DEPN, nm already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    else if (sam_seg_has_SA_Group(vb)) 
        sam_seg_against_sa_group (vb, ctx, add_bytes); 

    // method 3: use prediction of the number of mismatches derived by comparing SEQ to a reference.
    // this is usually, but surprisingly not always, correct for an external reference, and often correct for an internal one.
    else if (segconf.NM_is_integer && NM == predicted_by_SEQ)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'i'}, 3, OPTION_NM_i, add_bytes); 

    // case NM is a binary 0/1 rather than an integer. We use prediction against SEQ. TO DO: Support prediction against MD:Z too.
    else if (!segconf.NM_is_integer && predicted_by_SEQ != -1 && (NM > 0) == (predicted_by_SEQ > 0)) 
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_NM, 'b'}, 3, OPTION_NM_i, add_bytes); 

    else no_special: 
        seg_integer (VB, ctx, NM, true, add_bytes);

    // in PRIM with SA, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
    if (sam_is_prim_vb && sam_seg_has_SA_Group(vb)) {
        seg_integer_as_text (VB, OPTION_SA_NM, NM, 0);  // note: for PRIM lines without SA:Z and NM:i, we seg "0" into OPTION_SA_NM in sam_seg_sa_group_stuff

        // count NM field contribution to OPTION_SA_NM, so sam_stats_reallocate can allocate the z_data between NM and SA:Z
        CTX(OPTION_SA_NM)->counts.count += add_bytes; 
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
// XC:i (bwa) Undocumented: usually seq_len minus the final soft-clip (right if forward and left if rev-comp) 
// ----------------------------------------------------------------------------------------------------------
void sam_seg_XC_field (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t XC, unsigned add_bytes)
{
    ContextP ctx = CTX(OPTION_XC_i);
    int64_t prediction = dl->SEQ.len - vb->soft_clip[!dl->FLAG.bits.rev_comp || dl->FLAG.bits.unmapped];

    if (XC == prediction) 
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_XC }, 2, ctx, add_bytes);

    // case: prediction failed: cases: 1. not bwa 2. rev_comp != next_rev_comp
    else if (ABS(XC - prediction) < 1000) {  // if its not to big - seg a delta to avoid creating a local section for rare cases
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_XC, XC - prediction);
        seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
    }

    else {
        ctx->dynamic_size_local = true;
        ctx->flags.store = STORE_INT; 
        seg_integer (VB, ctx, XC, true, add_bytes);
    }
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_XC)
{
    int64_t delta=0;
    if (snip_len)  // we have a delta
        ASSPIZ (str_get_int (STRa(snip), &delta), "Invalid delta=\"%.*s\"", STRf(snip));

    new_value->i = vb->seq_len - VB_SAM->soft_clip[!last_flags.bits.rev_comp || last_flags.bits.unmapped] + delta;
    
    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// XT:A (bwa) Unique/Repeat/N/Mate-sw
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_XT_field (VBlockSAMP vb, char XT, unsigned add_bytes)
{
    bool X0_is_1 = CTX(OPTION_X0_i)->last_value.i == 1; // undefined unless has_X0
    
    // predict based on the existance and value of X0
    char prediction = !has_X0 ? 'M'  // Mate Smith-Waterman used for mapping. Note: could also be 'N', but we set our prediction to 'M' as it is very much more common
                    : X0_is_1 ? 'U'  // Unique mapping
                    :           'R'; // Repeat (i.e. not unique)

    if (prediction == XT)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_XT }, 2, OPTION_XT_A, add_bytes);

    else
        seg_by_did_i (VB, (char[]){ XT }, 1, OPTION_XT_A, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_XT)
{
    if (reconstruct) {
        char prediction = !container_does_contain (ctx, (DictId)_OPTION_X0_i)   ? 'M'
                        : reconstruct_peek (vb, CTX(OPTION_X0_i), 0, 0).i == 1 ? 'U'
                        :                                                        'R';
        RECONSTRUCT1 (prediction);
    }

    return NO_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------------------
// SM: Template-independent mapping quality
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_SM_field (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t SM, unsigned add_bytes)
{
    ContextP ctx = CTX(OPTION_SM_i);

    if (SM >= 0 && SM <= 255 && 
        SM != 254 &&           // note: 254 is a valid, but highly improbable value - we use 254 for "copy from MAPQ" so a actual 254 is segged as an exception
        !(SM && !dl->MAPQ)) {  // we're expecting XM=0 if MAPQ=0
        
        // store value in local (254 means "copy from MAPQ"), except if MAPQ=0 - we know SM=0 so no need to store
        if (dl->MAPQ) {
            uint8_t SM8 = (SM == dl->MAPQ) ? 254 : SM;
            seg_add_to_local_nonresizeable (VB, ctx, &SM8, false, 0);
        }

        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SM }, 2, ctx, add_bytes); // this usually is an all-the-same
    }

    else 
        seg_integer_as_text_do (VB, ctx, SM, add_bytes); 

    dl->SM = SM;
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
// AM: The smallest template-independent mapping quality in the template
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_AM_field (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t AM, unsigned add_bytes)
{
    ContextP ctx = CTX(OPTION_AM_i);

    // note: currently we only support for this algorithm AM appearing after SM. Easily fixable if ever needed.
    // AM is often one of 3 options: 0, =SM =MAPQ-SM. If SM=0 then AM is expected to be 0.
    if (AM >= 0 && AM <= 255 &&   // valid value
        AM != 253 && AM != 254 && // note: 253,254 are valid, but highly improbable values
        !(!dl->SM && AM)) {       // usually, AM=0 if SM=0
        
        if (dl->SM) { // no need to store if SM=0, as we know AM=0
            uint8_t AM8 = (AM == dl->SM)                  ? 253 // copy SM (note: AM is not 0, since SM is not 0)
                        : (AM && AM == dl->MAPQ - dl->SM) ? 254 // note subtelty: substraction in uint8_t arithmetics. we are careful to do the same in sam_piz_special_AM.
                        :                                   AM;

            seg_add_to_local_nonresizeable (VB, ctx, &AM8, false, 0);
        }
        
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_AM }, 2, ctx, add_bytes);
    }

    else  
        seg_integer_as_text_do (VB, ctx, AM, add_bytes);    
}    

SPECIAL_RECONSTRUCTOR (sam_piz_special_AM)
{
    uint8_t MAPQ = CTX(SAM_MAPQ)->last_value.i;
    uint8_t SM   = CTX(OPTION_SM_i)->last_value.i;

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
// X1:A (bwa) Number of suboptimal hits found by BWA
// ----------------------------------------------------------------------------------------------------------
static void sam_seg_X1_field (VBlockSAMP vb, int64_t X1, unsigned add_bytes)
{
    if (!has_X0) goto fallback; // need X0 to calculate prediction

    // in some files, where XA:Z includes sub-optimal alignments, this is usually true: 
    // (X0 + X1 = 1 + XA.repeats), and (X0 >= 1)
    int64_t prediction = has_XA ? CTX(OPTION_XA_Z)->last_value.i + 1 - CTX(OPTION_X0_i)->last_value.i : 0;

    if (prediction == X1)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_X1 }, 2, OPTION_X1_i, add_bytes);

    // note: rarely, X1 can be very large (tens of thousands)
    else 
        fallback:
        seg_integer (VB, CTX(OPTION_X1_i), X1, true, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_X1)
{
    new_value->i = 0; // default if no XA

    if (container_does_contain (ctx, (DictId)_OPTION_XA_Z)) {
        int64_t XA_repeats = reconstruct_peek_repeats (vb, CTX(OPTION_XA_Z), ';');
        int64_t X0 = reconstruct_peek (vb, CTX(OPTION_X0_i), 0, 0).i;
        
        new_value->i = XA_repeats + 1 - X0;
    }

    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// -------------------------------------------------------
// XA:Z "Alternative alignments" (a BWA-backtrack feature)
// -------------------------------------------------------

// Lookup buffer
#define lookback_buf zip_lookback_buf // we store the previous rname, pos, strand in their ctx->zip_lookback_buf buffer
#define lookback_value last_value.i

// Seg lookback callback for POS item of XA - seg the POS relative to lookback pos, and replace the already segged RNAME with a lookback if there is one.
static void sam_seg_XA_pos (VBlockP vb, STRp(pos_str), uint32_t rep)
{
    #define MAX_POS_DISTANCE 10000 // the smaller it is, the more we search for a better XA - the delta-pos will be better, but lookback worse. 10K works well.
    START_TIMER;

    ContextP rname_ctx = CTX (OPTION_XA_RNAME);
    ContextP pos_ctx   = CTX (OPTION_XA_POS);
    ContextP lb_ctx    = CTX (OPTION_XA_LOOKBACK);

    // look back for a node with this index and a similar POS - we use word_index to store the original rname_node_index, pos
    WordIndex rname_index = LASTb250(rname_ctx);
    int64_t lookback = 0;
    SamPosType pos = -MAX_POS_DISTANCE; // initial to "invalid pos" - value chosen so we can storeit in poses, in lieu of an invalid non-integer value, without a future pos being considered close to it

    if (!segconf.running && str_get_int_range32 (STRa(pos_str), 0, MAX_POS_SAM, &pos)) {

        int64_t iterator = -1;
        while ((lookback = lookback_get_next (vb, lb_ctx, rname_ctx, rname_index, &iterator))) {

            SamPosType lookback_pos = lookback_get_value (vb, lb_ctx, pos_ctx, lookback).i;

            // case: we found a lookback - same rname and close enough pos
            if (ABS (pos-lookback_pos) < MAX_POS_DISTANCE) {
            //    if (ABS (pos-lookback_pos) <= ABS(CTX(SAM_TLEN)->last_value.i)) { <-- better POS deltas but bigger index - its a wash

                // replace rname with lookback
                rname_ctx->b250.len--;
                ctx_decrement_count (vb, rname_ctx, rname_index);
                
                seg_by_ctx (vb, STRa(XA_lookback_snip), rname_ctx, 0); // add_bytes=0 bc we already added them when we segged rname the first time

                // seg pos as a delta
                SNIP(48);
                memcpy (snip, XA_lookback_snip, XA_lookback_snip_len);
                snip_len = XA_lookback_snip_len + str_int (pos - lookback_pos, &snip[XA_lookback_snip_len]);
                seg_by_ctx (vb, STRa(snip), pos_ctx, pos_str_len);
                
                break;
            }
        }
    }

    if (!lookback) 
        seg_integer_or_not (vb, pos_ctx, STRa(pos_str), pos_str_len);

    seg_add_to_local_resizable (vb, CTX(OPTION_XA_LOOKBACK), lookback, 0);

    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_RNAME, false, (ValueType){.i = rname_index }, true);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_POS,   false, (ValueType){.i = pos }, false);
    
    CTX(OPTION_XA_Z)->lookback_value = lookback; // for use when segging STRAND

    COPY_TIMER (sam_seg_XA_pos);
}

// Seg lookback callback for STRAND item of XA
static void sam_seg_XA_strand (VBlockP vb, WordIndex strand_index)
{
    ContextP strand_ctx = CTX (OPTION_XA_STRAND);
    int64_t lookback = CTX(OPTION_XA_Z)->lookback_value; // calculated in sam_seg_XA_pos

    if (lookback && lookback_get_index (vb, CTX(OPTION_XA_LOOKBACK), strand_ctx, lookback) == strand_index) 
        seg_by_ctx (vb, STRa(XA_lookback_snip), strand_ctx, 1);
    else
        seg_known_node_index (vb, strand_ctx, strand_index, 1); 

    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_STRAND, false, (ValueType){ .i = strand_index }, true);
}

// split the pos strand-pos string, eg "-10000" to strand "-" and pos "10000"
static bool seg_XA_strand_pos_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    if (field_len < 2 || (field[0] != '+' && field[0] != '-'))  
        return false; // invalid XA format - expecting pos to begin with the strand

    if (segconf.sam_is_sorted && !segconf.running) {
        sam_seg_XA_pos (vb, &field[1], field_len-1, rep);
        sam_seg_XA_strand (vb, *field == '+'); // index is 0 if '-' and 1 if '+' (set in sam_seg_0X_initialize)
    }
    // case: for a collated (or otherwise unsorted) file, we just seg normally (also: in segconf.running)
    else {
        seg_by_did_i (VB, field, 1, OPTION_XA_STRAND, 1);
        seg_integer_or_not (vb, CTX(OPTION_XA_POS), &field[1], field_len-1, field_len-1);
    }
    seg_by_ctx (VB, xa_strand_pos_snip, xa_strand_pos_snip_len, ctx, 0); // pre-created constant container

    return true; // segged successfully
}

static bool seg_XA_lookup_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    return true; // "segged successfully" - do nothing, we will seg it in seg_XA_strand_pos_cb
}

// XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
// Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
// See: http://bio-bwa.sourceforge.net/bwa.shtml
// Note that a different XA was observed in IonXpress data: "XA:Z:map4-1". This will be just segged as fallback
static void sam_seg_BWA_XA_field (VBlockSAMP vb, STRp(xa))
{
    static const MediumContainer container_XA = {
        .repeats      = 0, 
        .nitems_lo    = 5, 
        .filter_items = true, 
        .repsep       = {';'}, // including last item
        .items        = { { .dict_id = { _OPTION_XA_LOOKBACK   }, .separator = { CI0_INVISIBLE } }, 
                          { .dict_id = { _OPTION_XA_RNAME      }, .separator = {','} }, 
                          { .dict_id = { _OPTION_XA_STRAND_POS }, .separator = {','} },
                          { .dict_id = { _OPTION_XA_CIGAR      }, .separator = {','} }, // we don't mix the prirmary CIGAR field as the primary has a SNIP_SPECIAL
                          { .dict_id = { _OPTION_XA_NM         },                    } }  };

    static const SegCallback callbacks[5] = { [0]=seg_XA_lookup_cb, [1]=chrom_seg_cb, [2]=seg_XA_strand_pos_cb, [3]=sam_seg_0A_cigar_cb };
    int32_t repeats = seg_array_of_struct (VB, CTX(OPTION_XA_Z), container_XA, STRa(xa), callbacks);
    CTX(OPTION_XA_Z)->txt_len++; // +1 for '\t' in SAM or '\0' in BAM 

    // case: we failed to seg as a container - flush lookbacks (rare condition, and complicated to rollback given the round-robin and unlimited repeats)
    if (repeats == -1) {
        lookback_flush (VB, CTX(OPTION_XA_RNAME));
        lookback_flush (VB, CTX(OPTION_XA_POS));
        lookback_flush (VB, CTX(OPTION_XA_STRAND));
    }
}

static int sam_seg_which_XA (STRp(xa))
{
    unsigned semicolons;
    if ((semicolons = str_count_char (STRa(xa), ';')) >= 1 && str_count_char (STRa(xa), ',') == semicolons * 3) 
        return XA_BWA;

    if (memchr (xa, '-', xa_len))
        return XA_IONTORRENT;

    else
        return XA_UNKNOWN;
}

// this is called for XA that are a container, but not for invalid XA that are segged as a simple snip
void sam_piz_XA_field_insert_lookback (VBlockP vb)
{
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_RNAME,  true, (ValueType){ .i = 0 }, true);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_STRAND, true, (ValueType){ .i = 0 }, true);
    lookback_insert (vb, OPTION_XA_LOOKBACK, OPTION_XA_POS,    true, (ValueType){ .i = 0 }, false);
}

// -------------------------
// OA:Z "Original alignment"
// -------------------------

// OA is: (rname, pos, strand, CIGAR, mapQ, NM ;)+ . NM is optional (but its , is not)
// Example OA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
// See: https://samtools.github.io/hts-specs/SAMtags.pdf
static void sam_seg_OA_field (VBlockSAMP vb, STRp(field))
{
    static const MediumContainer container_OA = { .nitems_lo = 6,          
                                                  .repsep    = { ';' }, // including on last repeat    
                                                  .items     = { { .dict_id = { _OPTION_OA_RNAME  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_POS    }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_STRAND }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_CIGAR  }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_MAPQ   }, .separator = {','} },  
                                                                 { .dict_id = { _OPTION_OA_NM     },                    } } };

    SegCallback callbacks[6] = { [SA_RNAME]=chrom_seg_cb, [SA_POS]=seg_pos_field_cb, [SA_CIGAR]=sam_seg_0A_cigar_cb, [SA_MAPQ]=sam_seg_0A_mapq_cb };
     
    seg_array_of_struct (VB, CTX(OPTION_OA_Z), container_OA, STRa(field), callbacks);

    CTX(OPTION_OA_Z)->txt_len++; // 1 for \t in SAM and \0 in BAM 
}

// -------------------------------------------------------------------------------------------------------------------
// XM:i Case 1: BWA: "Number of mismatches in the alignment" 
//      Case 2: IonTorrent TMAP: "The target length, that is, the number of reference bases spanned by the alignment." 
// -------------------------------------------------------------------------------------------------------------------

static void sam_seg_XM_field (VBlockSAMP vb, ValueType XM, unsigned add_bytes)
{
    ContextP NM_ctx;
    
    // check for BWA case - XM is similar to NM - in our test file > 99% identical to NM.
    if (ctx_has_value_in_line (vb, _OPTION_NM_i, &NM_ctx) && XM.i == NM_ctx->last_value.i) 
        seg_by_did_i (VB, STRa(copy_NM_snip), OPTION_XM_i, add_bytes); // copy from NM
    
    // check IonTorrent TMAP case - XM is supposed to be ref_consumed
    else if (XM.i == vb->ref_consumed)                              
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED }, 2, OPTION_XM_i, add_bytes);
        
    else
        seg_integer (VB, CTX(OPTION_XM_i), XM.i, true, add_bytes);
}

// ----------------------------------------------------------------------------------------------
// AS:i "Alignment score generated by aligner" (https://samtools.github.io/hts-specs/SAMtags.pdf)
// ----------------------------------------------------------------------------------------------

// AS has a value set (at least as set by BWA and IonTorrent TMAP) of at most vb->ref_consumed, and often equal to it. we modify
// it to be new_value=(value-ref_consumed) 
static inline void sam_seg_AS_field (VBlockSAMP vb, ZipDataLineSAM *dl, ValueType AS, unsigned add_bytes)
{
    ctx_set_last_value (VB, CTX (OPTION_AS_i), AS);

    // in bowtie2, we might be able to copy from mate
    if (segconf.sam_bowtie2) {
        ASSERT (AS.i >= MIN_AS_i && AS.i <= MAX_AS_i, "%s: AS=%"PRId64" is out of range [%d,%d]", LN_NAME, AS.i, MIN_AS_i, MAX_AS_i);    
        
        ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_AS, zip_has_mate);

        if (zip_has_mate && mate_dl->YS == AS.i) 
            seg_by_ctx (VB, STRa(copy_buddy_YS_snip), channel_ctx, add_bytes);
        else
            // TODO: AS prediction, see bug 520
            seg_integer_as_text_do (VB, channel_ctx, AS.i, add_bytes);    

        seg_by_did_i (VB, STRa(vb->mux_AS.snip), OPTION_AS_i, 0); // de-multiplexor

        dl->AS = AS.i;
    }

    // not bowtie2: store a special snip with delta from ref_consumed
    else {
        SNIPi2 (SNIP_SPECIAL, SAM_SPECIAL_REF_CONSUMED, (int32_t)vb->ref_consumed-AS.i);
        seg_by_ctx (VB, STRa(snip), CTX (OPTION_AS_i), add_bytes); 
    }
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
    new_value->i = (int32_t)VB_SAM->ref_consumed - atoi (snip); // ref_consumed if snip=""
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return HAS_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// XS:i Suboptimal alignment score 
// XS:i BS-Seeker2: 1 when read is recognized as not fully converted by bisulfite treatment, or else 0
// ----------------------------------------------------------------------------------------------
static inline int sam_XS_get_mux_channel (uint8_t MAPQ)
{
    switch (MAPQ) {
        case 0         : return 0;
        case 1  ... 10 : return 1;
        case 11 ... 59 : return 2;
        default        : return 3;
    }    
}

static inline void sam_seg_XS_field (VBlockSAMP vb, ZipDataLineSAM *dl, ValueType XS, unsigned add_bytes)
{
    // "Suboptimal alignment score" - multiplex by MAPQ and (sometimes) delta vs AS
    if (!segconf.has_bsseeker2 && ctx_has_value_in_line_(VB, CTX(OPTION_AS_i)) && XS.i >= -10000 && XS.i < 10000) {

        int channel_i = sam_XS_get_mux_channel (dl->MAPQ);
        ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_XS, channel_i);

        if (channel_i == 3) {
            channel_ctx->dynamic_size_local = true;
            seg_integer (VB, channel_ctx, XS.i, false, add_bytes);
        }

        else if (XS.i) {
            ctx_set_last_value (VB, channel_ctx, XS); // needed for seg_delta_vs_other_do
            seg_delta_vs_other_do (VB, channel_ctx, CTX(OPTION_AS_i), NULL, 0, -1, add_bytes);
        }
        else
            seg_by_ctx (VB, "0", 1, channel_ctx, add_bytes);

        seg_by_did_i (VB, STRa(vb->mux_XS.snip), OPTION_XS_i, 0);
    }
    
    // BSSeeker2 or delta doesn't make sense - store as snip
    else
        seg_integer_as_text_do (VB, CTX(OPTION_XS_i), XS.i, add_bytes); // seg as text to not prevent singletons
}

// v14: De-multiplex XS by MAPQ
SPECIAL_RECONSTRUCTOR (sam_piz_special_XS)
{
    int channel_i = sam_XS_get_mux_channel (CTX(SAM_MAPQ)->last_value.i);
    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

// ----------------------------------------------------------------------------------------------
// YS:i mate alignment score (bowtie2 only)
// ----------------------------------------------------------------------------------------------

static inline void sam_seg_YS_field (VBlockSAMP vb, ZipDataLineSAM *dl, ValueType YS, unsigned add_bytes)
{
    ASSERT (YS.i >= MIN_AS_i && YS.i <= MAX_AS_i, "%s: YS=%"PRId64" is out of range [%d,%d]", LN_NAME, YS.i, MIN_AS_i, MAX_AS_i);    

    ctx_set_last_value (VB, CTX (OPTION_YS_i), YS);
    dl->YS = YS.i;

    ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_YS, zip_has_mate);

    if (zip_has_mate && mate_dl->AS == YS.i) 
        seg_by_ctx (VB, STRa(copy_buddy_AS_snip), channel_ctx, add_bytes);

    else 
        seg_integer_as_text_do (VB, channel_ctx, YS.i, add_bytes);    

    seg_by_did_i (VB, STRa(vb->mux_YS.snip), OPTION_YS_i, 0); // de-multiplexor
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

    bool do_mux = sam_is_main_vb && segconf.sam_is_paired; // for simplicity. To do: also for prim/depn components
    int channel_i = zip_has_mate?1 : zip_has_prim?2 : 0;
    ContextP channel_ctx = do_mux ? seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_MAPQ, channel_i) 
                                  : CTX(SAM_MAPQ);

    ZipDataLineSAM *mate_dl = DATA_LINE (vb->mate_line_i); // an invalid pointer if mate_line_i is -1

    if (sam_seg_has_SA_Group(vb)) {
        sam_seg_against_sa_group (vb, channel_ctx, add_bytes);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (sam_is_prim_vb) {
            ContextP sa_mapq_ctx = CTX(OPTION_SA_MAPQ);
            seg_add_to_local_nonresizeable (VB, sa_mapq_ctx, &dl->MAPQ, false, 0);

            // count MAPQ field contribution to OPTION_SA_MAPQ, so sam_stats_reallocate can allocate the z_data between MAPQ and SA:Z
            sa_mapq_ctx->counts.count += add_bytes; 
        }
    }

    // case: seg against mate
    else if (!segconf.MAPQ_has_single_value && segconf.has[OPTION_MQ_i] &&
             zip_has_mate && dl->MAPQ == mate_dl->MQ)
        seg_by_ctx (VB, STRa(copy_buddy_MQ_snip), channel_ctx, add_bytes); // copy MQ from earlier-line mate 

    // // case: seg against predicted alignment in prim line SA:Z 
    // else if (!segconf.MAPQ_has_single_value && zip_has_prim && 
    //          dl->MAPQ && dl->MAPQ != (dl-1)->MAPQ &&
    //          sam_seg_is_item_predicted_by_prim_SA (vb, SA_MAPQ, dl->MAPQ)) 
    //     seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_COPY_PRIM, '0'+SA_MAPQ }, 3, channel_ctx, add_bytes); 

    else {
        channel_ctx->ltype = LT_UINT8;
        seg_add_to_local_nonresizeable (VB, channel_ctx, &dl->MAPQ, true, add_bytes);
    }

//xx if (zip_has_prim) 
// printf ("mapq=%d pred=%d\n", dl->MAPQ, sam_seg_is_item_predicted_by_prim_SA (vb, SA_MAPQ, dl->MAPQ));
    if (do_mux)
        seg_by_did_i (VB, STRa(vb->mux_MAPQ.snip), SAM_MAPQ, 0); // de-multiplexor
}

// MQ:i Mapping quality of the mate/next segment
// Seg against mate if we have one, or else against MAPQ as it is often very similar
static inline void sam_seg_MQ_field (VBlockSAMP vb, ZipDataLineSAM *dl, ValueType MQ, unsigned add_bytes)
{
    segconf_set_has (OPTION_MQ_i);

    ASSERT (MQ.i >=0 && MQ.i <= 255, "%s: Invalid MQ:i=%"PRId64": expecting an integer [0,255]", LN_NAME, MQ.i);
    dl->MQ = MQ.i; 
    
    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_MQ, zip_has_mate);

    if (zip_has_mate && 
        dl->MQ == DATA_LINE (vb->mate_line_i)->MAPQ) 
        seg_by_ctx (VB, STRa(copy_buddy_MAPQ_snip), channel_ctx, add_bytes); // copy MAPQ from earlier-line mate 

    else {
        int64_t delta = MQ.i - (int64_t)dl->MAPQ;
        SNIP(32);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, _SAM_MAPQ, true, delta, snip);
        seg_by_ctx (VB, STRa(snip), channel_ctx, add_bytes);
    }

    seg_by_did_i (VB, STRa(vb->mux_MQ.snip), OPTION_MQ_i, 0); // de-multiplexor
}

// mc:i: (output of bamsormadup and other biobambam tools - mc in small letters) 
// appears to be a pos value usually close to PNEXT, but it is -1 is POS=PNEXT.
// from bamsort manual: "adddupmarksupport=<0|1>: add information required for streaming duplicate marking in the aux fields MS and MC.
// Input is assumed to be collated by query name. This option is ignored unless fixmates=1. By default it is disabled."
// https://github.com/gt1/biobambam2/blob/master/src/programs/bamsort.cpp says: "biobambam used MC as a mate coordinate tag which now has a clash
// with the official SAM format spec.  New biobambam version uses mc."
// ms="MateBaseScore" - sum all characters in QUAL of the mate, where the value of each character is its ASCII minus 33 (i.e. the Phred score)
// mc="MateCoordinate"
static inline void sam_seg_mc_field (VBlockSAMP vb, ValueType mc, unsigned add_bytes)
{
    // if snip is "-1", store as simple snip
    if (mc.i == -1)
        seg_by_did_i (VB, "-1", 2, OPTION_mc_i, add_bytes);
    
    // delta vs PNEXT
    else
        seg_pos_field (VB, OPTION_mc_i, SAM_PNEXT, SPF_BAD_SNIPS_TOO, 0, 0, 0, mc.i, add_bytes);
}

// multi-RG: possibly a special snip for copying RG from our buddy (mate or prim)
static inline void sam_seg_RG_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(rg), unsigned add_bytes)
{
    ZipDataLineSAM *buddy_dl = zip_has_mate ? DATA_LINE (vb->mate_line_i) 
                             : zip_has_prim ? DATA_LINE (vb->prim_line_i) // if not sorted, it will likely be the same as prev line, so better not use special 
                             : NULL;

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, (MultiplexerP)&vb->mux_RG, !!buddy_dl);

    if (segconf.sam_multi_RG && buddy_dl && str_issame_(STRa(rg), STRtxtw (buddy_dl->RG)))
        seg_by_ctx (VB, STRa(copy_buddy_RG_snip), channel_ctx, add_bytes); 
    else {
        channel_ctx->no_stons = true; // RGs are defined in the SAM header - there can't be too many of them
        seg_by_ctx (VB, STRa(rg), channel_ctx, add_bytes);    
    }

    seg_by_did_i (VB, STRa(vb->mux_RG.snip), OPTION_RG_Z, 0); // de-multiplexor

    dl->RG = TXTWORD(rg);
}

// optimization for Ion Torrent flow signal (ZM) - negative values become zero, positives are rounded to the nearest 10
static void sam_optimize_ZM (VBlockSAMP vb, ContextP ctx, void *cb_param, void *array_, uint32_t array_len)
{
    int16_t *array = (int16_t *)array_;

    for (uint32_t i=0; i < array_len; i++)
        if (array[i] >= 0) 
#ifdef __BIG_ENDIAN__
            array[i] = LTEN16 (((LTEN16(array[i]) + 5) / 10) * 10);
#else            
            array[i] = ((array[i] + 5) / 10) * 10;
#endif
        else array[i] = 0;
}

// E2 - SEQ data. Currently broken. To do: fix.
/*static void sam_seg_E2_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->SEQ.len, field, "E2 tag without a SEQ"); 
    ASSINP (field_len == dl->SEQ.len, 
            "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
            txt_name, dl->SEQ.len, field_len, field_len, field);

    SamPosType this_pos = vb->last_int(SAM_POS);

    sam_seg_SEQ (vb, OPTION_E2_Z, (char *)STRa(field), this_pos, vb->last_cigar, vb->ref_consumed, vb->ref_and_seq_consumed, 0, field_len, // remove const bc SEQ data is actually going to be modified
                        vb->last_cigar, add_bytes); 
}*/

// U2 - QUAL data (note: U2 doesn't have a context - it shares with QUAL)
static void sam_seg_U2_field (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(field), unsigned add_bytes)
{
    ASSSEG0 (dl->SEQ.len, field, "U2 tag without a SEQ"); 
    ASSINP (field_len == dl->SEQ.len, 
            "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. U2=%.*s",
            txt_name, dl->SEQ.len, field_len, field_len, field);

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
                             .drop_final_item_sep_of_final_repeat = true, // TODO - get rid of this flag and move to making the seperators to be repeat seperators as they should have been, using drop_final_repeat_sep and obsoleting this flag 
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

    #define SEG_COND(condition, seg) if (condition) { seg; break; } else goto fallback; 

    switch (dict_id.num) {

        case _OPTION_SA_Z: sam_SA_Z_seg (vb, dl, STRa(value)); break;

        case _OPTION_OA_Z: sam_seg_OA_field (vb, STRa(value)); break;

        case _OPTION_XA_Z: 
            if (!segconf.XA_type) segconf.XA_type = sam_seg_which_XA (STRa(value)); // doesn't matter if several VBs set this in a parallel, even if conflicting results. The last one wins.

            SEG_COND (segconf.XA_type == XA_BWA, sam_seg_BWA_XA_field (vb, STRa(value))); 
        
        case _OPTION_MC_Z: sam_cigar_seg_MC (vb, dl, STRa(value), add_bytes); break;

        // note: we cannot seg MD using SPECIAL if we're segging the line against SA Group, because
        // the MD alg requires the SQBITMAP.
        case _OPTION_MD_Z: sam_MD_Z_seg (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_NM_i: sam_seg_NM_field (vb, dl, (SamNMType)numeric.i, add_bytes); break;

        case _OPTION_BD_Z:
        case _OPTION_BI_Z: sam_seg_BD_BI_field (vb, dl, STRa(value), dict_id, add_bytes); break;
        
        case _OPTION_AS_i: sam_seg_AS_field (vb, dl, numeric, add_bytes); break;

        case _OPTION_XS_i: sam_seg_XS_field (vb, dl, numeric, add_bytes); break;

        case _OPTION_YS_i: SEG_COND (segconf.sam_bowtie2, sam_seg_YS_field (vb, dl, numeric, add_bytes));

        case _OPTION_XM_i: sam_seg_XM_field (vb, numeric, add_bytes); break;

        case _OPTION_MQ_i: sam_seg_MQ_field (vb, dl, numeric, add_bytes); break;

        case _OPTION_XO_Z: SEG_COND (segconf.has_bsseeker2, sam_seg_XO_Z_field (vb, dl, STRa(value), add_bytes));

        case _OPTION_XG_Z: SEG_COND (segconf.has_bsseeker2, sam_seg_XG_Z_field (vb, dl, STRa(value), add_bytes));

        case _OPTION_XM_Z: SEG_COND (segconf.has_bsseeker2, sam_seg_XM_Z_field (vb, dl, STRa(value), add_bytes));

        case _OPTION_XC_i: sam_seg_XC_field (vb, dl, numeric.i, add_bytes); break;

        case _OPTION_XT_A: sam_seg_XT_field (vb, value[0], add_bytes); break;

        case _OPTION_X1_i: sam_seg_X1_field (vb, numeric.i, add_bytes); break;

        case _OPTION_SM_i: sam_seg_SM_field (vb, dl, numeric.i, add_bytes); break;

        case _OPTION_AM_i: sam_seg_AM_field (vb, dl, numeric.i, add_bytes); break;

        case _OPTION_mc_i: sam_seg_mc_field (vb, numeric, add_bytes); break;

        case _OPTION_ms_i: segconf_set_has (OPTION_ms_i);
                           SEG_COND (segconf.sam_ms_type == ms_BIOBAMBAM, sam_seg_ms_field (vb, numeric, add_bytes));

        case _OPTION_RG_Z: SEG_COND (segconf.sam_multi_RG, sam_seg_RG_field (vb, dl, STRa(value), add_bytes));

        // tx:i: - we seg this as a primary field SAM_TAX_ID
        case _OPTION_tx_i: seg_by_did_i (VB, taxid_redirection_snip, taxid_redirection_snip_len, OPTION_tx_i, add_bytes); break;

        //case _OPTION_E2: sam_seg_E2_field (vb, dl, STRa(value), add_bytes); // BROKEN. To do: fix.

        case _OPTION_U2_Z: sam_seg_U2_field (vb, dl, STRa(value), add_bytes); break;

        case _OPTION_Z5_i: seg_pos_field (VB, OPTION_Z5_i, SAM_PNEXT, 0, 0, 0, 0, numeric.i, add_bytes); break;

        case _OPTION_CP_i: SEG_COND (segconf.sam_is_sorted, seg_pos_field (VB, OPTION_CP_i, OPTION_CP_i, 0, 0, 0, 0, numeric.i, add_bytes)); 

        case _OPTION_ZM_B_s : sam_seg_array_field (vb, (DictId)_OPTION_ZM_B_s, array_subtype, STRa(value), 
                                                   (segconf.tech == TECH_IONTORR && flag.optimize_ZM) ? sam_optimize_ZM : NULL, 0);
                              break;

        default: fallback:
            
            // all types of integer
            if (sam_type == 'i') {
                ContextP ctx = ctx_get_ctx (vb, dict_id);
                ctx->dynamic_size_local = true;
                ctx->flags.store = STORE_INT; // needs this be reconstructable as BAM
                seg_integer (VB, ctx, numeric.i, true, add_bytes);
            }

            else if (sam_type == 'f') {
                ContextP ctx = ctx_get_ctx (vb, dict_id);
                ctx->flags.store = STORE_FLOAT; // needs this be reconstructable as BAM

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
