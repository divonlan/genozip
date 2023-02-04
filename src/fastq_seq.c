// ------------------------------------------------------------------
//   fastq_seq.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "seg.h"
#include "piz.h"
#include "deep.h"
#include "aligner.h"
#include "refhash.h"
#include "coverage.h"

static void fastq_get_pair_1_gpos_strand (VBlockFASTQP vb, PosType *pair_gpos, bool *pair_is_forward)
{
    ContextP bitmap_ctx = CTX(FASTQ_SQBITMAP);

    // case: we are pair-1 OR we are pair-2, but this line in pair-1 is not aligned
    if (!bitmap_ctx->pair_b250 || 
        ({ STR(snip); ctx_get_next_snip (VB, bitmap_ctx, true, pSTRa(snip)); *snip != SNIP_LOOKUP; })) { // note: SNIP_LOOKUP in case of aligned, SNIP_SPECIAL in case of unaligned

        *pair_gpos = NO_GPOS;
        *pair_is_forward = 0;
    }

    // case: we are pair-2, and the corresponding line in pair-1 is aligned: get its gpos and is_forward
    else {
        ContextP gpos_ctx   = CTX(FASTQ_GPOS);
        ContextP strand_ctx = CTX(FASTQ_STRAND);

        ASSERT (gpos_ctx->pair.next < gpos_ctx->pair.len32, "%s: not enough data GPOS.pair (len=%u)", LN_NAME, gpos_ctx->pair.len32); 

        ASSERT (gpos_ctx->pair.next < strand_ctx->pair.nbits, "%s: cannot get pair_1 STRAND bit because pair_1 strand bitarray has only %u bits",
                LN_NAME, (unsigned)strand_ctx->pair.nbits);

        *pair_gpos = (PosType)*B32 (gpos_ctx->pair, gpos_ctx->pair.next); 
        *pair_is_forward = bits_get ((BitsP)&strand_ctx->pair, gpos_ctx->pair.next); 
        gpos_ctx->pair.next++;
    }
}

static inline bool seq_len_by_qname (VBlockFASTQP vb, uint32_t seq_len)
{
    // TO DO: add seq_len_by_qname also to the case of aligned sequence ^ 
    return segconf.qname_seq_len_dict_id.num &&  // QNAME flavor has "length=""
           seq_len == ECTX(segconf.qname_seq_len_dict_id)->last_value.i; // length is equal seq_len
}

void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool deep)
{
    // note: monobase reads don't work well with deep - too much contention. In illumina, we observe many monobase N reads,
    // and in for R2, also many monobase G reads. These seem to be an artifact of the technology rather than true sequence values.
    if (dl->monobase) {
        SNIPi3 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, seq[0], seq_len_by_qname (vb, seq_len) ? 0 : seq_len);
        seg_by_ctx (VB, STRa(snip), CTX(FASTQ_SQBITMAP), 0); 
        CTX(FASTQ_NONREF)->txt_len += seq_len; // account for the txt data in NONREF        
// printf ("xxx ZIP: monobase: %s %c\n", LN_NAME, seq[0]);
        return;
    }
            
    if (deep) { 
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_copy_deep }, 2, CTX(FASTQ_SQBITMAP), seq_len); 
        return;
    }

    // get pair-1 gpos and is_forward, but only if we are pair-2 and the corresponding pair-1 line is aligned
    PosType pair_gpos; 
    bool pair_is_forward;
    fastq_get_pair_1_gpos_strand (vb, &pair_gpos, &pair_is_forward); // does nothing if we are pair-1
               
    // case: aligned - lookup from SQBITMAP
    MappingType aln_res;
    if (flag.aligner_available && ((aln_res = aligner_seg_seq (VB, CTX(FASTQ_SQBITMAP), STRa(seq), true, (vb->comp_i == FQ_COMP_R2), pair_gpos, pair_is_forward)))) 
        seg_lookup_with_length (VB, CTX(FASTQ_SQBITMAP), (aln_res==MAPPING_PERFECT ? -(int32_t)seq_len : (int32_t)seq_len), seq_len);

    // case: not aligned - just add data to NONREF
    else {
        if (flag.show_aligner && !segconf.running) iprintf ("%s: unaligned\n", LN_NAME);

        Context *nonref_ctx = CTX(FASTQ_NONREF);
        buf_alloc (VB, &nonref_ctx->local, seq_len + 3, vb->txt_data.len / 64, char, CTX_GROWTH, "contexts->local"); 
        buf_add (&nonref_ctx->local, seq, seq_len);

        // case: we don't need to consume pair-1 gpos (bc we are pair-1, or pair-1 was not aligned): look up from NONREF
        SNIPi3 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, ' '/*copy from NONREF*/, seq_len_by_qname (vb, seq_len) ? 0 : seq_len);
        seg_by_ctx (VB, STRa(snip), CTX(FASTQ_SQBITMAP), 0); // note: FASTQ_SQBITMAP is always segged whether read is aligned or not
        CTX(FASTQ_NONREF)->txt_len += seq_len; // account for the txt data in NONREF
    }
}

COMPRESSOR_CALLBACK (fastq_zip_seq) 
{
    STRtxtset (*line_data, DATA_LINE (vb_line_i)->seq);
    if (is_rev) *is_rev = 0;
}

//------------------------
// Pair-2 GPOS (ZIP & PIZ)
//------------------------

void fastq_seg_pair2_gpos (VBlockP vb, PosType pair1_pos, PosType pair2_gpos)
{
    #define MAX_GPOS_DELTA 1000 // paired reads are usually with a delta less than 300 - so this is more than enough

    ContextP ctx = CTX(FASTQ_GPOS);

    // case: we are pair-2 ; pair-1 is aligned ; delta is small enough : store a delta
    PosType gpos_delta = pair2_gpos - pair1_pos; 
    if (pair1_pos != NO_GPOS && pair2_gpos != NO_GPOS &&
        gpos_delta <= MAX_GPOS_DELTA && gpos_delta >= -MAX_GPOS_DELTA) {
    
        SNIPi2(SNIP_SPECIAL, FASTQ_SPECIAL_PAIR2_GPOS, gpos_delta);
        seg_by_ctx (VB, STRa(snip), ctx, 0);
    }

    // case: otherwise, store verbatim
    else {
        BNXT32 (ctx->local) = (uint32_t)pair2_gpos;
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_PAIR2_GPOS }, 2, ctx, 0); // lookup from local and advance pair.next to consume gpos
    }
}

// note: in files up to v13, this is triggereed by v13_SNIP_FASTQ_PAIR2_GPOS
SPECIAL_RECONSTRUCTOR (fastq_special_PAIR2_GPOS)
{
    // case: no delta
    if (!snip_len) {
        if (CTX(FASTQ_SQBITMAP)->pair1_is_aligned == PAIR1_ALIGNED)
            ctx->pair.next++; // we didn't use this pair value 

        new_value->i = NEXTLOCAL(uint32_t, ctx);
    }

    // case: pair-1 is aligned, and GPOS is a delta vs pair-1
    else {
        int64_t pair_value = (int64_t)(VER(14) ? *B32 (ctx->pair, ctx->pair.next++) // starting v14, only aligned lines have GPOS
                                               : *B32 (ctx->pair, vb->line_i));     // up to v13, all lines segged GPOS (possibly NO_GPOS value, but not in this case, since we have a delta)
        int64_t delta = (int64_t)strtoull (snip, NULL, 10 /* base 10 */); 
        new_value->i = pair_value + delta; // just sets value, doesn't reconstruct

        ASSPIZ (pair_value != NO_GPOS, "pair_value=NO_GPOS - not expected as we have delta=%d", (int)delta);
    }
 
    return HAS_NEW_VALUE;
}

// can only be called before fastq_special_PAIR2_GPOS, because it inquires GPOS.pair.next
bool fastq_piz_get_pair2_is_forward (VBlockP vb)
{
    ContextP ctx = CTX(FASTQ_STRAND);

    // case: paired read (including all reads in up to v13) - diff vs pair1
    if (CTX(FASTQ_SQBITMAP)->pair1_is_aligned == PAIR1_ALIGNED) { // always true for files up to v13 and all lines had a is_forward value

        bool is_forward_pair_1 = VER(14) ? bits_get ((BitsP)&ctx->pair, CTX(FASTQ_GPOS)->pair.next) // since v14, gpos_ctx->pair.next is an iterator for both gpos and strand, and is incremented in fastq_special_PAIR2_GPOS
                                         : bits_get ((BitsP)&ctx->pair, vb->line_i);                // up to v13, all lines had strand, which was 0 if unmapped

        return NEXTLOCALBIT (ctx) ? is_forward_pair_1 : !is_forward_pair_1;
    }

    // case: unpaired read - just take bit
    else
        return NEXTLOCALBIT (ctx);
}

// called when reconstructing R2, to check if R1 is aligned
static bool fastq_piz_R1_test_aligned (VBlockP vb)
{
    ContextP ctx = CTX(FASTQ_SQBITMAP);

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (ctx->pair1_is_aligned == PAIR1_ALIGNED_UNKNOWN) { // case: not already set by fastq_special_mate_lookup
        STR(snip);
        ctx_get_next_snip (vb, ctx, true, pSTRa(snip));
        ctx->pair1_is_aligned = (snip_len && *snip == SNIP_LOOKUP) ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;
    }
    
    return (ctx->pair1_is_aligned == PAIR1_ALIGNED);
}

static void fastq_update_coverage_aligned (VBlockFASTQP vb)
{
    PosType gpos;
    Context *gpos_ctx = CTX(FASTQ_GPOS);

    if (vb->comp_i == FQ_COMP_R1) 
        gpos = NEXTLOCAL (uint32_t, gpos_ctx);

    else { // pair-2
        reconstruct_from_ctx (VB, FASTQ_GPOS, 0, false); // calls fastq_special_PAIR2_GPOS
        gpos = gpos_ctx->last_value.i;
    }
    
    ASSPIZ0 (gpos != NO_GPOS, "expecting a GPOS, because sequence is aligned");

    WordIndex ref_index = ref_contig_get_by_gpos (gref, gpos, NULL);
    ASSPIZ0 (ref_index != WORD_INDEX_NONE, "expecting ref_index, because sequence is aligned");

    if (flag.show_coverage || flag.show_sex)
        *B64 (vb->coverage, ref_index) += vb->seq_len;

    if (flag.show_coverage || flag.idxstats)
        (*B64 (vb->read_count, ref_index))++;
}

// PIZ: aligned SEQ reconstruction - called by reconstructing FASTQ_SQBITMAP which is a LOOKUP (either directly, or via fastq_special_mate_lookup)
void fastq_recon_aligned_SEQ (VBlockP vb_, ContextP bitmap_ctx, STRp(seq_len_str), ReconType reconstruct)
{
    VBlockFASTQP vb = (VBlockFASTQP )vb_;

    fastq_piz_R1_test_aligned (VB); // set pair1_is_aligned
    
    // v14: perfect alignment is expressed by a negative seq_len
    bool perfect_alignment = (seq_len_str[0] == '-');
    if (perfect_alignment) { seq_len_str++; seq_len_str_len--; }

    ASSERT (str_get_int_range32 (STRa(seq_len_str), 0, 0x7fffffff, (int32_t*)&vb->seq_len), "seq_len_str=\"%.*s\" out range [0,0x7fffffff]", STRf(seq_len_str));

    // just update coverage
    if (flag.collect_coverage) 
        fastq_update_coverage_aligned (vb);

    // --qual-only: only set vb->seq_len without reconstructing
    else if (flag.qual_only) {}

    // normal reconstruction
    else 
        aligner_reconstruct_seq (vb_, bitmap_ctx, vb->seq_len, vb->comp_i == FQ_COMP_R2, perfect_alignment, reconstruct);
}

// PIZ: SEQ reconstruction - in case of unaligned sequence 
SPECIAL_RECONSTRUCTOR (fastq_special_unaligned_SEQ)
{
    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (vb->comp_i == FQ_COMP_R2)
        if (fastq_piz_R1_test_aligned (vb) || !VER(14)) // up to v13, even non-aligned reads had a GPOS entry
            CTX(FASTQ_GPOS)->pair.next++; // gpos_ctx->pair.next is an iterator for both gpos and strand

    // just update coverage (unaligned)
    if (flag.collect_coverage) {
        if (flag.show_coverage || flag.show_sex)
            *(BAFT64 (vb->coverage) - NUM_COVER_TYPES + CVR_UNMAPPED) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*(BAFT64 (vb->read_count) - NUM_COVER_TYPES + CVR_UNMAPPED))++;
    }

    else {
        if (flag.show_aligner) iprintf ("%s: unaligned\n", LN_NAME);
        
        char monobase = 0;
        if (VER(15)) { 
            if (snip[0] != ' ') monobase = snip[0];
            STRinc (snip);
        }

        if (snip_len==1 && *snip=='0') 
            vb->seq_len = ECTX(segconf.qname_seq_len_dict_id)->last_value.i;
        else
            vb->seq_len = atoi(snip);

        // case: take seq_len from DESC item with length=
        if (!monobase) 
            reconstruct_from_local_sequence (vb, CTX(FASTQ_NONREF), 0, 0, reconstruct);

        else {
// printf ("xxx PIZ: monobase: %s %c\n", LN_NAME, monobase);

            char *start = BAFTtxt;
            uint32_t len = vb->seq_len;
            
            for (uint32_t i=0; i < len; i++)
                start[i] = monobase;
            
            vb->txt_data.len32 += len;
        }
    }

    return NO_NEW_VALUE;
}

