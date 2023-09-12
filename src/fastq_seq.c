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
    
#define SEQ_LEN_BY_QNAME 0x7fffffff

static void fastq_get_pair_1_gpos_strand (VBlockFASTQP vb, PosType64 *pair_gpos, bool *pair_is_forward)
{
    declare_seq_contexts;

    STR(snip);
    ctx_get_next_snip (VB, bitmap_ctx, true, pSTRa(snip));
    if (*snip != SNIP_LOOKUP) return; // pair-1 is unaligned

    ASSERT (gpos_ctx->localR1.next < gpos_ctx->localR1.len32, "%s: not enough data GPOS.localR1 (len=%u)", LN_NAME, gpos_ctx->localR1.len32); 

    ASSERT (gpos_ctx->localR1.next < strand_ctx->localR1.nbits, "%s: cannot get pair_1 STRAND bit because pair_1 strand bits has only %u bits",
            LN_NAME, (unsigned)strand_ctx->localR1.nbits);

    // the corresponding line in pair-1 is aligned: get its gpos and is_forward
    *pair_gpos = (PosType64)*B32(gpos_ctx->localR1, gpos_ctx->localR1.next); 
    *pair_is_forward = bits_get ((BitsP)&strand_ctx->localR1, gpos_ctx->localR1.next); 
    gpos_ctx->localR1.next++; // iterator for both GPOS and STRAND
}

static inline bool seq_len_by_qname (VBlockFASTQP vb, uint32_t seq_len)
{
    return segconf.seq_len_dict_id.num &&  // QNAME or FASTQ_AUX have a length (eg "length=")
           seq_len == ECTX(segconf.seq_len_dict_id)->last_value.i; // length is equal seq_len
}

void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ *dl, STRp(seq), bool deep)
{
    START_TIMER;

    declare_seq_contexts;

    // get pair-1 gpos and is_forward, but only if we are pair-2 and the corresponding pair-1 line is aligned
    PosType64 pair_gpos = NO_GPOS; 
    bool pair_is_forward = false;
    if (vb->pair_vb_i) // R2
        fastq_get_pair_1_gpos_strand (vb, &pair_gpos, &pair_is_forward); // advance iterators even if we don't need the pair data

    // note: monochar reads don't work well with deep - too much contention. In illumina, we observe many monochar N reads,
    // and in for R2, also many monochar G reads. These seem to be an artifact of the technology rather than true sequence values.
    if (dl->monochar) {
        SNIPi3 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, seq[0], seq_len_by_qname (vb, seq_len) ? 0 : seq_len);
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, 0); 
        nonref_ctx->txt_len += seq_len; // account for the txt data in NONREF        

        goto done;
    }
            
    if (deep) { 
        fastq_deep_seg_SEQ (vb, dl, STRa(seq), bitmap_ctx, nonref_ctx);
        goto done;
    }
               
    bool aligner_ok = flag.aligner_available && !segconf.is_long_reads && !segconf.running;

    // case: aligner - lookup from SQBITMAP
    MappingType aln_res;
    if (aligner_ok && ((aln_res = aligner_seg_seq (VB, STRa(seq), true, (vb->pair_vb_i > 0), pair_gpos, pair_is_forward)))) {
    
        int32_t pseudo_seq_len = seq_len_by_qname (vb, seq_len) ? SEQ_LEN_BY_QNAME : seq_len;    

        seg_lookup_with_length (VB, bitmap_ctx, (aln_res==MAPPING_PERFECT ? -(int32_t)pseudo_seq_len : (int32_t)pseudo_seq_len), seq_len);
    }

    // case: not aligned - just add data to NONREF
    else {
        if (flag.show_aligner && !segconf.running) iprintf ("%s: unaligned\n", LN_NAME);

        seg_add_to_local_fixed (VB, nonref_ctx, STRa(seq), LOOKUP_NONE, seq_len);

        // case: we don't need to consume pair-1 gpos (bc we are pair-1, or pair-1 was not aligned): look up from NONREF
        SNIPi3 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, ' '/*copy from NONREF*/, seq_len_by_qname (vb, seq_len) ? 0 : seq_len);
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, 0); // note: FASTQ_SQBITMAP is always segged whether read is aligned or not
    }

done:
    if (seq_len > vb->longest_seq_len) vb->longest_seq_len = seq_len;

    COPY_TIMER (fastq_seg_SEQ);
}

// used by QUAL codecs: LONGR and HOMP
COMPRESSOR_CALLBACK (fastq_zip_seq) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);
    bool trimmed = flag.deep && (dl->seq.len > dl->sam_seq_len);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = trimmed ? MIN_(maximum_size, dl->seq.len - dl->sam_seq_len) // compress only trimmed bases, other bases will be copied from Deep
                    :           MIN_(maximum_size, dl->seq.len);
    
    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->seq.index) + (trimmed ? dl->sam_seq_len : 0);

    if (is_rev) *is_rev = 0;
}

//------------------------
// Pair-2 GPOS (ZIP & PIZ)
//------------------------

void fastq_seg_pair2_gpos (VBlockP vb, PosType64 pair1_pos, PosType64 pair2_gpos)
{
    #define MAX_GPOS_DELTA 1000 // paired reads are usually with a delta less than 300 - so this is more than enough
    declare_seq_contexts;

    // case: we are pair-2 ; pair-1 is aligned ; delta is small enough : store a delta
    PosType64 gpos_delta = pair2_gpos - pair1_pos; 
    if (pair1_pos != NO_GPOS && pair2_gpos != NO_GPOS && ABS(gpos_delta) <= MAX_GPOS_DELTA) {
        int16_t gpos_delta16 = gpos_delta;
        seg_integer_fixed (VB, gpos_d_ctx, &gpos_delta16, false, 0);
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_PAIR2_GPOS, 'D' }, 3, gpos_ctx, 0); // lookup from local and advance localR1.next to consume gpos
    }

    // case: otherwise, store verbatim
    else {
        BNXT32 (gpos_ctx->local) = (uint32_t)pair2_gpos;
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_PAIR2_GPOS }, 2, gpos_ctx, 0); // lookup from local and advance localR1.next to consume gpos
    }
}

// note: in files up to v13, this is triggereed by v13_SNIP_FASTQ_PAIR2_GPOS
SPECIAL_RECONSTRUCTOR (fastq_special_PAIR2_GPOS)
{
    declare_seq_contexts;

    // case: no delta
    if (!snip_len) {
        if (bitmap_ctx->pair1_is_aligned == PAIR1_ALIGNED)
            gpos_ctx->localR1.next++; // we didn't use this pair value 

        new_value->i = NEXTLOCAL(uint32_t, gpos_ctx);
    }

    // case: pair-1 is aligned, and GPOS is a delta vs pair-1. 
    else {   
        ASSPIZ (gpos_ctx->localR1.next < gpos_ctx->localR1.len, "gpos_ctx->localR1.next=%"PRId64" overflow: gpos_ctx->localR1.len=%u",
                gpos_ctx->localR1.next, gpos_ctx->localR1.len32);

        int64_t pair_value = (int64_t)(VER(14) ? *B32 (gpos_ctx->localR1, gpos_ctx->localR1.next++) // starting v14, only aligned lines have GPOS
                                               : *B32 (gpos_ctx->localR1, vb->line_i));                // up to v13, all lines segged GPOS (possibly NO_GPOS value, but not in this case, since we have a delta)
        
        int16_t delta = VER(15) ? reconstruct_from_local_int (VB, gpos_d_ctx, 0, RECON_OFF) // starting v15, delta is stored in FASTQ_GPOS_DELTA.local
                                : (int64_t)strtoull (snip, NULL, 10 /* base 10 */);         // up to v14, delta was encoded in the snip

        new_value->i = pair_value + delta; // just sets value, doesn't reconstruct

        ASSPIZ (pair_value != NO_GPOS, "pair_value=NO_GPOS - not expected as we have delta=%d", (int)delta);
    }
 
    return HAS_NEW_VALUE;
}

// can only be called before fastq_special_PAIR2_GPOS, because it inquires GPOS.localR1.next
bool fastq_piz_get_pair2_is_forward (VBlockP vb)
{
    declare_seq_contexts;

    // defect 2023-02-11: until 14.0.30, we allowed dropping SQBITMAP.b250 sections if all the same (since 14.0.31, we 
    // set no_drop_b250). Due to the defect, if b250 is dropped, we always segged is_forward verbatim and not as a diff to R1.
    bool defect_2023_02_11 = !bitmap_ctx->b250R1.len32 && EXACT_VER(14); // can only happen up to 14.0.30

    // case: paired read (including all reads in up to v13) - diff vs pair1
    if (!defect_2023_02_11 && bitmap_ctx->pair1_is_aligned == PAIR1_ALIGNED) { // always true for files up to v13 and all lines had a is_forward value
        ASSPIZ (!VER(14) || gpos_ctx->localR1.next < strand_ctx->localR1.nbits, "gpos_ctx->localR1.next=%"PRId64" overflow: strand_ctx->localR1.nbits=%"PRIu64,
                gpos_ctx->localR1.next, strand_ctx->localR1.nbits);

        bool is_forward_pair_1 = VER(14) ? bits_get ((BitsP)&strand_ctx->localR1, gpos_ctx->localR1.next) // since v14, gpos_ctx->localR1.next is an iterator for both gpos and strand, and is incremented in fastq_special_PAIR2_GPOS
                                         : bits_get ((BitsP)&strand_ctx->localR1, vb->line_i);            // up to v13, all lines had strand, which was 0 if unmapped

        return NEXTLOCALBIT (strand_ctx) ? is_forward_pair_1 : !is_forward_pair_1;
    }

    // case: unpaired read - just take bit
    else
        return NEXTLOCALBIT (strand_ctx);
}

// called when reconstructing R2, to check if R1 is aligned
bool fastq_piz_R1_test_aligned (VBlockFASTQP vb)
{
    declare_seq_contexts;

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (bitmap_ctx->pair1_is_aligned == PAIR1_ALIGNED_UNKNOWN) { // case: not already set by fastq_special_mate_lookup
        STR(snip);
        ctx_get_next_snip (VB, bitmap_ctx, true, pSTRa(snip));
        bitmap_ctx->pair1_is_aligned = (snip_len && *snip == SNIP_LOOKUP) ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;
    }
    
    return (bitmap_ctx->pair1_is_aligned == PAIR1_ALIGNED);
}

static void fastq_update_coverage_aligned (VBlockFASTQP vb)
{
    declare_seq_contexts;
    PosType64 gpos;

    if (vb->comp_i == FQ_COMP_R1) 
        gpos = NEXTLOCAL (uint32_t, gpos_ctx);

    else { // pair-2
        reconstruct_from_ctx (VB, FASTQ_GPOS, 0, false); // calls fastq_special_PAIR2_GPOS
        gpos = gpos_ctx->last_value.i;
    }
    
    ASSPIZ0 (gpos != NO_GPOS, "expecting a GPOS, because sequence is aligned");

    WordIndex ref_index = ref_contig_get_by_gpos (gref, gpos, 0, NULL);
    ASSPIZ0 (ref_index != WORD_INDEX_NONE, "expecting ref_index, because sequence is aligned");

    if (flag.show_coverage || flag.show_sex)
        *B64 (vb->coverage, ref_index) += vb->seq_len;

    if (flag.show_coverage || flag.idxstats)
        (*B64 (vb->read_count, ref_index))++;
}

// PIZ: reconstruct_seq callback: aligned SEQ reconstruction - called by reconstructing FASTQ_SQBITMAP which is a LOOKUP (either directly, or via fastq_special_mate_lookup)
void fastq_recon_aligned_SEQ (VBlockP vb_, STRp(seq_len_str), ReconType reconstruct)
{
    VBlockFASTQP vb = (VBlockFASTQP )vb_;
    declare_seq_contexts;

    if (vb->pair_vb_i) // R2 
        fastq_piz_R1_test_aligned (vb); // set pair1_is_aligned
    
    // v14: perfect alignment is expressed by a negative seq_len
    bool perfect_alignment = (seq_len_str[0] == '-');
    if (perfect_alignment) { seq_len_str++; seq_len_str_len--; }

    ASSERT (str_get_int_range32 (STRa(seq_len_str), 0, 0x7fffffff, (int32_t*)&vb->seq_len), "seq_len_str=\"%.*s\" out range [0,0x7fffffff]", STRf(seq_len_str));

    if (vb->seq_len == SEQ_LEN_BY_QNAME) // introduced v15
        vb->seq_len = reconstruct_peek_by_dict_id (VB, segconf.seq_len_dict_id, 0, 0).i; // peek, since length can come from either line1 or line3

    // just update coverage
    if (flag.collect_coverage) 
        fastq_update_coverage_aligned (vb);

    // --qual-only: only set vb->seq_len without reconstructing
    else if (flag.qual_only) {}

    // normal reconstruction
    else 
        aligner_reconstruct_seq (vb_, vb->seq_len, vb->pair_vb_i > 0, perfect_alignment, reconstruct, NULL, NULL, NULL);
}

// PIZ: SEQ reconstruction - in case of unaligned sequence 
SPECIAL_RECONSTRUCTOR (fastq_special_unaligned_SEQ)
{
    declare_seq_contexts;

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (VB_FASTQ->pair_vb_i) // R2
        if (fastq_piz_R1_test_aligned (VB_FASTQ) || !VER(14)) // up to v13, even non-aligned reads had a GPOS entry
            gpos_ctx->localR1.next++; // gpos_ctx->localR1.next is an iterator for both gpos and strand

    // just update coverage (unaligned)
    if (flag.collect_coverage) {
        if (flag.show_coverage || flag.show_sex)
            *(BAFT64 (vb->coverage) - NUM_COVER_TYPES + CVR_UNMAPPED) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*(BAFT64 (vb->read_count) - NUM_COVER_TYPES + CVR_UNMAPPED))++;
    }

    else {
        if (flag.show_aligner) iprintf ("%s: unaligned\n", LN_NAME);
        
        char monochar = 0;
        if (VER(15)) { 
            if (snip[0] != ' ') monochar = snip[0];
            STRinc (snip, 1);
        }
        
        if (snip_len==1 && *snip=='0') 
            vb->seq_len = reconstruct_peek_by_dict_id (vb, segconf.seq_len_dict_id, 0, 0).i; // peek, since length can come from either line1 or line3
        else
            vb->seq_len = atoi(snip);

        // --qual-only: only set vb->seq_len without reconstructing
        if (flag.qual_only) {}

        // case: take seq_len from DESC item with length=
        else if (!monochar) 
            reconstruct_from_local_sequence (vb, nonref_ctx, vb->seq_len, reconstruct);

        else {
            memset (BAFTtxt, monochar, vb->seq_len);
            Ltxt += vb->seq_len;
        }
    }

    return NO_NEW_VALUE;
}

