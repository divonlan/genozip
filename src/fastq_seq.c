// ------------------------------------------------------------------
//   fastq_seq.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fastq_private.h"
#include "deep.h"
#include "aligner.h"
#include "refhash.h"
#include "coverage.h"
    
#define SEQ_LEN_BY_QNAME 0x7fffffff
#define NONBIO_EXCESS_ALIGNED '@'
#define NONBIO_CONTAINERIZED  '^'

static void fastq_get_pair_1_gpos_strand (VBlockFASTQP vb, PosType64 *gpos_R1, bool *is_forward_R1)
{
    declare_seq_contexts;

    STR(snip);
    ctx_get_next_snip (VB, bitmap_ctx, true, pSTRa(snip));
    if ((snip_len >= 1 && snip[0] == SNIP_LOOKUP) || // pair-1 has gpos from our aligner 
        (snip_len >= 2 && snip[0] == SNIP_SPECIAL && snip[1] == FASTQ_SPECIAL_SEQ_by_bamass)) { // pair-1 has gpos from bamass

        ASSERT (gpos_ctx->localR1.next < gpos_ctx->localR1.len32, "%s: not enough data GPOS.localR1 (len=%u)", LN_NAME, gpos_ctx->localR1.len32); 

        ASSERT (gpos_ctx->localR1.next < strand_ctx->localR1.nbits, "%s: cannot get pair_1 STRAND bit because pair_1 strand bits has only %u bits",
                LN_NAME, (unsigned)strand_ctx->localR1.nbits);

        // the corresponding line in pair-1 is aligned: get its gpos and is_forward
        *is_forward_R1 = bits_get ((BitsP)&strand_ctx->localR1, gpos_ctx->localR1.next); 
        *gpos_R1 = reconstruct_from_pair_int (vb, gpos_ctx); // also increments gpos_ctx->localR1.next: iterator for both GPOS and STRAND
    }
}

static inline bool seq_len_by_qname (VBlockFASTQP vb, uint32_t seq_len)
{
    return segconf.seq_len_dict_id.num &&  // QNAME or FASTQ_AUX have a length (eg "length=")
           seq_len == ECTX(segconf.seq_len_dict_id)->last_value.i; // length is equal seq_len
}

void fastq_seg_SEQ (VBlockFASTQP vb, ZipDataLineFASTQ𐤐  dl, STRp(seq), bool deeped)
{
    START_TIMER;

    declare_seq_contexts;
    MappingType mapping_type;

    // get pair-1 gpos and is_forward, but only if we are pair-2 and the corresponding pair-1 line is aligned
    PosType64 gpos_R1 = NO_GPOS; 
    bool is_forward_R1 = false;

    if (seq_len == 28 && vb->line_i==0 && segconf_running && !IS_R2)
        segconf.nonbio_type = NONBIO_10xGen; // optimistically - we will decide finally based on aligner results

    bool aligner_ok = flag.aligner_available && seq_len && !segconf.is_long_reads && 
                      (!segconf_running || IS_NONBIO(10xGen));
    bool is_excess_aligned;
    bool am_i_R2 = false;

    // case: R2 in paired file
    if (vb->R1_vb_i) {
        fastq_get_pair_1_gpos_strand (vb, &gpos_R1, &is_forward_R1); // advance iterators even if we don't need the pair data
        am_i_R2 = true;
    }

    // case: R2 in interleaved file
    else if (segconf.is_interleaved && vb->line_i % 2) {
        if (ctx_encountered_in_prev_line (VB, FASTQ_GPOS)) { // note: encountered only if segged with aligner
            gpos_R1 = gpos_ctx->last_value.i;
            is_forward_R1 = strand_ctx->last_value.i;
        }
        am_i_R2 = true;
    }

    // ---------------------------------
    //            SEG methods 
    // ---------------------------------

    if (deeped && flag.deep) 
        fastq_deep_seg_SEQ (vb, dl, STRa(seq), bitmap_ctx, nonref_ctx);

    // case: Parse Biosciences or SPLiT-seq non-biological read (added 15.0.83)
    else if (IS_NONBIO(Parse)) {
        if (!fastq_parse_seg_SEQ (vb, dl, STRa(seq), gpos_R1, is_forward_R1, &is_excess_aligned)) goto verbatim;

        // note: length is not needed, because it can deduced from the nonbio container
        seg_special1 (VB, FASTQ_SPECIAL_unaligned_SEQ, is_excess_aligned ? NONBIO_EXCESS_ALIGNED : NONBIO_CONTAINERIZED, bitmap_ctx, 0);        
        vb->num_nonbio++;
    }

    else if (IS_NONBIO(10xGen) && vb->comp_i == FQ_COMP_R1) {
        if (!fastq_10xGen_seg_SEQ (vb, STRa(seq))) goto verbatim;

        seg_special1 (VB, FASTQ_SPECIAL_unaligned_SEQ, NONBIO_CONTAINERIZED, bitmap_ctx, 0);        
        vb->num_nonbio++;
    }

    else if (deeped && flag.bam_assist) {
        fastq_bamass_seg_CIGAR (vb);
        
        mapping_type = fastq_bamass_seg_SEQ (vb, dl, STRa(seq), am_i_R2, gpos_R1, is_forward_R1);
        
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_SEQ_by_bamass, '0'+MT(PERFECT) }, 3, FASTQ_SQBITMAP, seq_len); 
        vb->num_bamass++;
    }

    else if (!seq_len) { // empty reads supported since 15.0.81
        seg_special1 (VB, FASTQ_SPECIAL_unaligned_SEQ, '*', bitmap_ctx, 0);        

        vb->num_empty_read++;
    }

    else if (dl->monochar) { // note: seq_len=0 is considered monochar
        SNIPi3 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, seq[0], seq_len_by_qname (vb, seq_len) ? 0 : seq_len); 
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, 0); 
        nonref_ctx->txt_len += seq_len; // account for the txt data in NONREF        
        
        vb->num_monochar++;
    }
                           
    // case: aligner - lookup from SQBITMAP
    else if (aligner_ok && ((mapping_type = aligner_seg_seq (VB, STRa(seq), am_i_R2, gpos_R1, is_forward_R1)))) {
        int32_t pseudo_seq_len = seq_len_by_qname (vb, seq_len) ? SEQ_LEN_BY_QNAME : seq_len;    

        STRl(snip, 32) = 0;
        snip[snip_len++] = SNIP_LOOKUP;
        if (MT(PERFECT) || MT(PERFECT_SPLICED)) snip[snip_len++] = ALIGNED_PERFECT;
        if (MT(SPLICED) || MT(PERFECT_SPLICED)) snip[snip_len++] = ALIGNED_SPLICED; // 15.0.83
        snip_len += str_int (pseudo_seq_len, &snip[snip_len]);
        
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, seq_len);
    }

    // case: not aligned - just add data to NONREF verbatim
    else verbatim: {
        if (flag.show_aligner && !segconf_running) iprintf ("%s: verbatim\n", LN_NAME);

        seg_add_to_local_fixed (VB, nonref_ctx, STRa(seq), LOOKUP_NONE, seq_len);

        // case: we don't need to consume pair-1 gpos (bc we are pair-1, or pair-1 was not aligned): look up from NONREF
        SNIPi3 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, ' '/*copy from NONREF*/, seq_len_by_qname (vb, seq_len) ? 0 : seq_len);
        seg_by_ctx (VB, STRa(snip), bitmap_ctx, 0); // note: FASTQ_SQBITMAP is always segged whether read is aligned or not
    
        vb->num_verbatim++; // for stats

        // test: if pratically all non-bamassed reads are unalignable, then turn off the aligner (test once every 16K reads)
        if (flag.bam_assist && flag.aligner_available && vb->line_i % 16384 == 16383)
            fastq_bamass_consider_stopping_aligner (vb);
    }

    MAXIMIZE (vb->longest_seq_len, seq_len);

    COPY_TIMER (fastq_seg_SEQ);
}

uint32_t fastq_zip_get_seq_len (VBlockP vb, uint32_t line_i) 
{ 
    return DATA_LINE (line_i)->sam_seq_len;
}

// used by SEQ-dependent QUAL codecs: LONGR, PACB, SMUX and HOMP
COMPRESSOR_CALLBACK (fastq_zip_seq) 
{
    ZipDataLineFASTQ𐤐  dl = DATA_LINE (vb_line_i);
    bool trimmed = flag.deep && (dl->seq.len > dl->sam_seq_len);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    // note: see bug 1208
    *line_data_len  = trimmed ? MIN_(maximum_size, dl->seq.len - dl->sam_seq_len) // compress only trimmed bases, other bases will be copied from Deep
                    :           MIN_(maximum_size, dl->seq.len);
    
    if (__builtin_expect (!line_data, false)) return; // only lengths were requested

    *line_data = Btxt (dl->seq.index) + (trimmed ? dl->sam_seq_len : 0);

    if (is_rev) *is_rev = 0;
}

//-----------------------------------------------------
// Pair-2 GPOS, and Single-file interleaved (ZIP & PIZ)
//-----------------------------------------------------

void fastq_seg_gpos_R2 (VBlockP vb, PosType64 gpos_R1, PosType64 gpos_R2, bool is_forward_R2)
{
    declare_seq_contexts;
    ContextP my_gpos_ctx = segconf.is_interleaved ? gpos_r2_ctx : gpos_ctx;
    
    #define MAX_GPOS_DELTA (1 MB - 1)

    // case: we are pair-2 ; pair-1 is aligned ; delta is small enough : store a delta
    PosType64 gpos_Δ = is_forward_R2 ? (gpos_R1 - gpos_R2) : (gpos_R2 - gpos_R1); 

    if (gpos_R1 != NO_GPOS && gpos_R2 != NO_GPOS && ABS(gpos_Δ) <= MAX_GPOS_DELTA) {
        seg_integer (VB, gpos_Δ_ctx, gpos_Δ, false, 0);
        seg_special1 (VB, FASTQ_SPECIAL_PAIR2_GPOS, segconf.is_interleaved ? 'I' : 'D', my_gpos_ctx, 0); // lookup from local and advance localR1.next to consume gpos

        ctx_set_last_value (vb, gpos_Δ_ctx, gpos_Δ); // consumed in aligner_seg_seq

        VB_FASTQ->num_r2_gpos_Δ++;
    }

    // case: otherwise, store verbatim
    else {
        seg_integer (VB, my_gpos_ctx, gpos_R2, false, 0);
        seg_special0 (VB, FASTQ_SPECIAL_PAIR2_GPOS, my_gpos_ctx, 0); 
    }
}

// note: in files up to v13, this is triggered by v13_SNIP_FASTQ_PAIR2_GPOS
// since 15.0.58, also used for R2 in single-file interleaved
SPECIAL_RECONSTRUCTOR (fastq_special_PAIR2_GPOS)
{
    declare_seq_contexts;

    // case: no delta
    if (!snip_len) {
        if (bitmap_ctx->r1_is_aligned == PAIR1_ALIGNED)
            ctx->localR1.next++; // we didn't use this pair value 

        new_value->i = reconstruct_from_local_int (vb, ctx, 0, RECON_OFF);
    }

    // case: pair-1 is aligned, and GPOS is a delta vs pair-1. 
    else {   
        PosType64 gpos_r1;

        // case: delta vs pair-1 data (snip is "D" since v15 or a textual integer up to v14)
        if (*snip != 'I') { 
            ASSPIZ (ctx->localR1.next < ctx->localR1.len, "gpos_ctx->localR1.next=%"PRId64" overflow: gpos_ctx->localR1.len=%u",
                    ctx->localR1.next, ctx->localR1.len32);

            gpos_r1 = (PosType64)(VER(14) ? reconstruct_from_pair_int (VB_FASTQ, ctx) // starting v14, only aligned lines have GPOS. starting 15.0.69 gpos is dyn_int.
                                          : *B32 (ctx->localR1, vb->line_i));   // up to v13, all lines segged GPOS (possibly NO_GPOS value, but not in this case, since we have a delta)
        }

        // case: delta vs previous line in interleaved file (snip is "I") since 15.0.58
        else 
            gpos_r1 = gpos_ctx->last_value.i;

        int64_t delta = VER(15) ? reconstruct_from_local_int (vb, gpos_Δ_ctx, 0, RECON_OFF) // starting v15, delta is stored in FASTQ_GPOS_DELTA.local
                                : (int64_t)strtoull (snip, NULL, 10 /* base 10 */);         // up to v14, delta was encoded in the snip

        if (strand_ctx->last_value.i/*=is_forward_R2*/ && VER2(15,83))
            delta = -delta;
            
        new_value->i = gpos_r1 + delta; // just sets value, doesn't reconstruct

        ASSPIZ (gpos_r1 != NO_GPOS, "gpos_r1=NO_GPOS - not expected as we have delta=%d", (int)delta);
    }
 
    return HAS_NEW_VALUE;
}

// can only be called before fastq_special_PAIR2_GPOS, because it inquires GPOS.localR1.next
bool fastq_piz_get_r2_is_forward (VBlockP vb)
{
    declare_seq_contexts;

    // defect 2023-02-11: until 14.0.30, we allowed dropping SQBITMAP.b250 sections if all the same (since 14.0.31, we 
    // set no_drop_b250). Due to the defect, if b250 is dropped, we always segged is_forward verbatim and not as a diff to R1.
    bool defect_2023_02_11 = !bitmap_ctx->b250R1.len32 && EXACT_VER(14); // can only happen up to 14.0.30

    // case: paired read (including all reads in up to v13) - diff vs pair1
    if (!defect_2023_02_11 && bitmap_ctx->r1_is_aligned == PAIR1_ALIGNED) { // always true for files up to v13 and all lines had a is_forward value
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

bool fastq_piz_get_interleaved_r2_is_forward (VBlockP vb)
{
    declare_seq_contexts;

    // case: paired read - diff vs r1
    if (bitmap_ctx->r1_is_aligned == PAIR1_ALIGNED) {
        bool is_forward_pair_1 = strand_ctx->last_value.i;
        return NEXTLOCALBIT (strand_r2_ctx) ? is_forward_pair_1 : !is_forward_pair_1;
    }

    // case: unpaired read - just take bit
    else
        return NEXTLOCALBIT (strand_r2_ctx);
}

// called when reconstructing R2, to check if R1 is aligned
bool fastq_piz_R1_test_aligned (VBlockFASTQP vb)
{
    declare_seq_contexts;

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (bitmap_ctx->r1_is_aligned == PAIR1_ALIGNED_UNKNOWN) { // case: not already set by fastq_special_mate_lookup
        STR(snip);
        ctx_get_next_snip (VB, bitmap_ctx, true, pSTRa(snip)); // paired file - get snip from R1

        // note: bitmap_ctx is segged with LOOKUP if aligned and SNIP_SPECIAL if not aligned 
        // bitmap_ctx->r1_is_aligned = (snip_len && *snip == SNIP_LOOKUP) ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;
        bitmap_ctx->r1_is_aligned = ((snip_len >= 1 && snip[0] == SNIP_LOOKUP) || // pair-1 has gpos from our aligner 
                                     (snip_len >= 2 && snip[0] == SNIP_SPECIAL && snip[1] == FASTQ_SPECIAL_SEQ_by_bamass)) 
                                  ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;
    }
    
    return (bitmap_ctx->r1_is_aligned == PAIR1_ALIGNED);
}

void fastq_update_coverage_aligned (VBlockFASTQP vb)
{
    declare_seq_contexts;
    PosType64 gpos;

    if (vb->comp_i == FQ_COMP_R1) 
        gpos = NEXTLOCAL (uint32_t, gpos_ctx);

    else { // pair-2
        reconstruct_from_ctx (VB, FASTQ_GPOS, 0, false); // calls fastq_special_PAIR2_GPOS
        gpos = gpos_ctx->last_value.i;
    }

    // TO DO: interleaved files, see aligner_recon_get_gpos_and_fwd
    
    ASSPIZ0 (gpos != NO_GPOS, "expecting a GPOS, because sequence is aligned");

    WordIndex ref_index = ref_contig_get_by_gpos (gpos, 0, NULL, true); // if gpos is in a gap between to contigs, it means that bulk of seq is on the next contig while its beginning is in the gap
    ASSPIZ0 (ref_index != WORD_INDEX_NONE, "expecting ref_index, because sequence is aligned");

    if (flag.show_coverage)
        *B64 (vb->coverage, ref_index) += vb->seq_len;

    if (flag.show_coverage || flag.idxstats)
        (*B64 (vb->read_count, ref_index))++;
}

// PIZ: reconstruct_seq callback: aligned SEQ reconstruction - called by reconstructing FASTQ_SQBITMAP which is a LOOKUP (either directly, or via fastq_special_mate_lookup)
void fastq_recon_aligned_SEQ (VBlockP vb_, STRp(snip), ReconType reconstruct)
{
    VBlockFASTQP vb = (VBlockFASTQP )vb_;
    declare_seq_contexts;

    if (vb->R1_vb_i) // R2 
        fastq_piz_R1_test_aligned (vb); // set r1_is_aligned
    
    // v14: perfect alignment is expressed by a negative seq_len
    bool perfect_alignment = (snip[0] == ALIGNED_PERFECT);
    if (perfect_alignment) STRinc(snip, 1);

    bool spliced_alignment = (snip[0] == ALIGNED_SPLICED);
    if (spliced_alignment) STRinc(snip, 1);

    ASSERT (str_get_int_range32 (STRa(snip), 0, 0x7fffffff, (int32_t*)&vb->seq_len), "Invalid seq_len=\"%.*s\"", STRf(snip));

    if (vb->seq_len == SEQ_LEN_BY_QNAME) // introduced v15
        vb->seq_len = reconstruct_peek_by_dict_id (VB, segconf.seq_len_dict_id, 0, 0).i; // peek, since length can come from either line1 or line3

    // just update coverage
    if (flag.collect_coverage) 
        fastq_update_coverage_aligned (vb);

    // --qual-only: only set vb->seq_len without reconstructing
    else if (flag.qual_only) {}

    else if (spliced_alignment) {
        int64_t junction = reconstruct_from_local_int (VB, junction_ctx, 0, false); // negative means first segment uses ref2
        ctx_set_last_value (VB, junction_ctx, junction); // consumed by aligner_recon_get_gpos_and_fwd

        // reconstruct first and then second segment
        aligner_reconstruct_seq (VB, junction/*seq_len*/,    vb->R1_vb_i > 0, false, perfect_alignment, reconstruct, 0, NULL, NULL, NULL);
        PosType64 G1 = gpos_ctx->last_value.i;

        aligner_reconstruct_seq (VB, vb->seq_len - junction, vb->R1_vb_i > 0, true,  perfect_alignment, reconstruct, 0, NULL, NULL, NULL);
        PosType64 G2 = gpos_ctx->last_value.i;
        bool fwd     = strand_ctx->last_value.i;

        // set GPOS.last_value as it would have been for non-spliced. consumed if interleaved - 2nd read's gpos is delta vs first
        gpos_ctx->last_value.i = fwd ? G1 : (G2 - gpos_gap_ctx->last_value.i);
    }

    // normal reconstruction
    else 
        aligner_reconstruct_seq (VB, vb->seq_len, vb->R1_vb_i > 0, false, perfect_alignment, reconstruct, 0, NULL, NULL, NULL);
}

// Used in R2: used for pair-assisted b250 reconstruction. copy parallel b250 snip from R1. 
// Since v14, used for only for SQBITMAP. For files up to v14, also used for all QNAME subfield contexts
SPECIAL_RECONSTRUCTOR (fastq_special_mate_lookup)
{
    ASSPIZ (ctx->b250R1.len32, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2.", ctx->tag_name);
            
    ctx_get_next_snip (vb, ctx, true, pSTRa(snip));

    if (ctx->did_i == FASTQ_SQBITMAP && VER(14))
        ctx->r1_is_aligned = (snip_len && *snip == SNIP_LOOKUP) ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;

    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE /* we can't cache pair items */, STRa(snip), reconstruct, __FUNCLINE); // might include delta etc - works because in --pair, ALL the snips in a context are FASTQ_SPECIAL_mate_lookup

    return NO_NEW_VALUE; // last_value already set (if needed) in reconstruct_one_snip
}

// PIZ: SEQ reconstruction - in case of unaligned sequence 
SPECIAL_RECONSTRUCTOR (fastq_special_unaligned_SEQ)
{
    declare_seq_contexts;

    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (VB_FASTQ->R1_vb_i && snip[0] != NONBIO_EXCESS_ALIGNED) // R2
        if (fastq_piz_R1_test_aligned (VB_FASTQ) || !VER(14)) // up to v13, even non-aligned reads had a GPOS entry
            gpos_ctx->localR1.next++; // gpos_ctx->localR1.next is an iterator for both gpos and strand

    // just update coverage (unaligned)
    if (flag.collect_coverage) {
        if (flag.show_coverage)
            *(BAFT64 (vb->coverage) - NUM_COVER_TYPES + CVR_UNMAPPED) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*(BAFT64 (vb->read_count) - NUM_COVER_TYPES + CVR_UNMAPPED))++;
    }

    else {
        // case: non-biological (containerized) sequence
        if (VER(15) && (snip[0] == NONBIO_EXCESS_ALIGNED || snip[0] == NONBIO_CONTAINERIZED)) {
            uint32_t len_before = Ltxt;
            reconstruct_from_ctx (vb, FASTQ_NONBIO, 0, true); // always reconstruct, so we can calculate seq_len
            vb->seq_len = Ltxt - len_before;

            if (!reconstruct) Ltxt = len_before;
            goto done;
        }

        if (flag.show_aligner) iprintf ("%s: unaligned\n", LN_NAME);
        
        if (VER(15) && snip[0] == '*') { // empty reads supported since 15.0.81
            vb->seq_len = 0;
            goto done;
        }

        char monochar = 0;

        if (VER(15)) { 
            if (snip[0] != ' ') monochar = snip[0];
            STRinc (snip, 1);
        }
        
        if (IS_CHAR0(snip)) // '0' means "no seq_len_dict_id", it does not mean seq_len=0 
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

    done: return NO_NEW_VALUE;
}

