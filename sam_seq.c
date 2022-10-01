// ------------------------------------------------------------------
//   sam_seq.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "random_access.h"
#include "aligner.h"
#include "regions.h"
#include "codec.h"
#include "refhash.h"

//---------------
// Shared ZIP/PIZ
//---------------

#define REF(idx) ref_base_by_idx (range, (idx))

#define piz_is_set codec_bufs[0] // bytemap

// called when SEQ in a prim or supp/sec line is segged against ref against prim: 
// Sets vb->mismatch_bases_by_SEQ. ZIP: also sets vb->md_verified and. PIZ: also returns sqbitmap of this line in line_sqbitmap 
static bool sam_analyze_copied_SEQ (VBlockSAMP vb, STRp(seq), const SamPosType pos, bool is_revcomp,
                                    uint32_t ref_consumed, uint32_t ref_and_seq_consumed,
                                    BufferP line_sqbitmap)
{
    RefLock lock = REFLOCK_NONE;

    BitsP bitmap = (BitsP)line_sqbitmap;
    uint32_t bit_i=0;

    if (vb->cigar_missing) goto fail;
    
    ASSERT (!line_sqbitmap->len32, "%s: line_sqbitmap is in use", LN_NAME);
    buf_alloc_bits (VB, line_sqbitmap, ref_and_seq_consumed, line_sqbitmap->name ? line_sqbitmap->name : "line_sqbitmap");        
    bits_set_region (bitmap, 0, ref_and_seq_consumed); // we initialize all the bits to "set", and clear as needed.

    ConstRangeP range = IS_ZIP ? ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, seq, 
                                                    IS_REF_EXTERNAL ? NULL : &lock)
                               : ref_piz_get_range (VB, gref, true);
    if (!range) goto fail; // can happen in ZIP/REF_INTERNAL due to contig cache contention ; in ZIP/PIZ with an external reference - contig is not in the reference

    uint32_t pos_index = pos - range->first_pos;
    uint32_t next_ref  = pos_index;

    uint32_t i=0, n=0;
    BamCigarOpType op = BC_NONE;

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
        
    while (i < seq_len || next_ref < pos_index + ref_consumed) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_consumed, "i or next_ref are out of range");

        if (vb->binary_cigar.next < vb->binary_cigar.len) {
            BamCigarOp *next_op = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next++);
            op = next_op->op;
            n  = next_op->n;
        }

        switch (op) {

            case BC_M: case BC_E: case BC_X: // alignment match or sequence match or mismatch

                ASSERT (n > 0 && n <= (seq_len - i), "%s: CIGAR implies seq_len longer than actual seq_len=%u (line=%s n=%u). CIGAR=\"%.*s\"", 
                        LN_NAME, seq_len, line_name(VB).s, n, STRfb(vb->textual_cigar));

                uint32_t start_i = i;
                for (; n && next_ref < pos_index + ref_consumed; n--, next_ref++, i++) {

                    // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                    uint32_t idx = (next_ref >= range_len) ? next_ref - range_len : next_ref; 

                    // case: we don't have the value for the REF site, therefore we can't seg MD against REF (its not worth it to add the REF site just for MD)
                    if (IS_ZIP && (flag.reference & REF_STORED) && !ref_is_nucleotide_set (range, idx)) 
                        goto fail;

                    // case our seq is identical to the reference at this site
                    else if (IS_NUCLEOTIDE (seq[i]) && seq[i] == REF(idx)) 
                        bit_i++; // bit remains set. 
 
                    // case: ref is set to a different value - update the bitmap
                    else {
                        bits_clear (bitmap, bit_i++);
                        vb->mismatch_bases_by_SEQ++;
                    } 
                }
                ref_and_seq_consumed -= (i - start_i); // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap

                break; // end if 'M', '=', 'X'

            // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
            case BC_I: case BC_S: 
                ASSERT (n > 0 && n <= (seq_len - i), "%s: CIGAR \"%s\" (n_ops=%u) implies seq_len longer than actual seq_len=%u", LN_NAME, vb->last_cigar ? vb->last_cigar : "", vb->binary_cigar.len32, seq_len);
                i += n;
                n = 0;
                break;

            // for Deletion or Skipping - we move the next_ref ahead
            case BC_D: case BC_N: {
                unsigned ref_consumed_skip = (IS_REF_INTERNAL ? MIN_(n, range_len - next_ref)
                                                                             : n);
                next_ref += ref_consumed_skip;
                n        -= ref_consumed_skip;
                break;
            }

            // Hard clippping (H) or padding (P) - nothing much to do
            case BC_H: case BC_P: 
                n = 0;
                break;

            case BC_NONE:
                ABORT ("%s: End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                        "(cigar=%s pos=%d seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u op.n=%u range=[%.*s %"PRId64"-%"PRId64"])",
                        LN_NAME, pos_index + ref_consumed - next_ref, seq_len-i, vb->last_cigar ? vb->last_cigar : "", pos, seq_len,
                        ref_consumed, next_ref, pos_index, n, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            default:
                ABORT ("%s: Invalid CIGAR op=%u", LN_NAME, op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_consumed && n) break;
    }

    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE

    // an error here can indicate that CIGAR is inconsistent with sequence
    ASSERT (i == seq_len, "%s: expecting i(%u) == seq_len(%u) pos=%d range=[%.*s %"PRId64"-%"PRId64"] (possibly reason: inconsistency between seq_len and CIGAR=\"%s\")", 
            LN_NAME, i, seq_len, pos, STRf(range->chrom_name), range->first_pos, range->last_pos, (vb->last_cigar ? vb->last_cigar : ""));

    // ZIP: final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (IS_ZIP) {
        sam_MD_Z_verify_due_to_seq (vb, STRa(seq), pos, bitmap, 0); // sets vb->md_verified
        buf_free (*line_sqbitmap);
    }
    
    return true;

fail:
    // case ZIP: if we cannot verify against the reference, the MD:Z is not verified, and we don't use the SPECIAL for segging
    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE
    buf_free (*line_sqbitmap);
    vb->md_verified = false;
    vb->mismatch_bases_by_SEQ = -1; // therefore NM cannot seg against mismatch_bases_by_SEQ
    
    return false;
}

//---------
// ZIP
//---------

void sam_seg_SEQ_initialize (VBlockSAMP vb)
{
    Context *bitmap_ctx = CTX(SAM_SQBITMAP);
    Context *strand_ctx = CTX(SAM_STRAND);
    Context *gpos_ctx   = CTX(SAM_GPOS);
    Context *nonref_ctx = CTX(SAM_NONREF);
    Context *seqmis_ctx = CTX(SAM_SEQMIS_A); // 4 contexts

    bitmap_ctx->ltype = LT_BITMAP;
    strand_ctx->ltype = LT_BITMAP;
    gpos_ctx->ltype   = LT_UINT32;

    // MAIN: we may seg depn lines against in-VB prim lines
    if (sam_is_main_vb)
        bitmap_ctx->flags.store_per_line = true; // 14.0.0

    // initial allocations, these might grow during segging if needed
    int factor = segconf.sam_is_unmapped ? 1 : 32; // if mapped, we allocate assuming 1 out of 32 lines is unmapped
    buf_alloc (vb, &bitmap_ctx->local, 1, vb->txt_data.len / (4 * factor), uint8_t, 0, "contexts->local"); 
    buf_alloc (vb, &strand_ctx->local, 1, roundup_bits2bytes64 (vb->lines.len / factor), uint8_t, 0, "contexts->local"); 
    buf_alloc (vb, &gpos_ctx->local, 1, vb->lines.len / factor, uint32_t, CTX_GROWTH, "contexts->local"); 

    if (!segconf.sam_is_unmapped || flag.aligner_available) {
        buf_alloc (vb, &nonref_ctx->local, 0, vb->txt_data.len / 64, char, 0, "contexts->local");

        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, 1, vb->txt_data.len / 128, char, 0, "contexts->local"); 
    }
    else // we store seq vertabim - no reference and no CIGAR
        buf_alloc (vb, &nonref_ctx->local, 0, vb->txt_data.len / 3, char, 0, "contexts->local");
}

// align nonref_ctx->local to a 4-character boundary. this is because CODEC_ACGT squeezes every 4 characters into a byte,
// before compressing it with LZMA. In sorted SAM, we want subsequent identical sequences to have the same byte alignment
// so that LZMA can catch their identicality. (note: as of v13, we assign a codec rather than hard-coding LZMA, but it is usually LZMA anyway)
static inline void sam_seg_SEQ_pad_nonref (VBlockSAMP vb)
{
    ContextP ctx = CTX(SAM_NONREF);
    uint32_t add_chars = (4 - (ctx->local.len32 & 3)) & 3;
    if (add_chars) buf_add (&ctx->local, "AAA", add_chars); // add 1 to 3 As
}

// SEG: called from sam_seg_SEQ_vs_ref to handle an M segment of a bisulfite-treated SEQ
static void sam_seg_bisulfite_M (VBlockSAMP vb, 
                                 RangeP range, int32_t range_len, int32_t idx,
                                 STRp(Mseg), uint32_t M_i, // M segment of a SEQ string
                                 BitsP bitmap, uint32_t first_bit,
                                 Context *seqmis_ctx) // 4 contexts
{
    START_TIMER;

    char bisulfite = vb->bisulfite_strand;
    bool MD_NM_by_unconverted = segconf.MD_NM_by_unconverted;
    idx = RR_IDX (idx);

    // iterating on an M segmenet of a SEQ
    for (uint32_t i=0; i < Mseg_len; i++) {
    
        char base = Mseg[i];
        char unconverted_ref_base = REF(idx);
        
        // C->T or G->A conversion
        char converted_ref_base = (bisulfite == unconverted_ref_base) ? ((bisulfite == 'C') ? 'T' : 'A') 
                                                                      : unconverted_ref_base;

        // if needed update unconverted_bitmap - bitmap vs unconverted reference
        if (MD_NM_by_unconverted) {
            vb->unconverted_bitmap.nbits++;
        
            // mismatch vs unconverted
            if (base != unconverted_ref_base) { 
                bits_clear ((BitsP)&vb->unconverted_bitmap, vb->unconverted_bitmap.nbits - 1);
                vb->mismatch_bases_by_SEQ++; // used by NM:i and MD:Z
            }
        }
        
        // case: seq base mismatches converted reference base
        if (base != converted_ref_base) {
            uint8_t ref_base_2bit = acgt_encode[(uint8_t)unconverted_ref_base];

            BNXTc (seqmis_ctx[ref_base_2bit].local) = base;
            bits_clear (bitmap, first_bit + i);

            if (!MD_NM_by_unconverted)
                vb->mismatch_bases_by_SEQ++; // used by NM:i and MD:Z - in case we are calculate MD and NM vs converted reference
        }

        // set methylation call. note: bisulfite is not supported in REF_INTERNAL, so we know range is the entire contig
        if (segconf.sam_predict_meth_call &&
            (bisulfite == unconverted_ref_base && (base == converted_ref_base || base == unconverted_ref_base)))

            sam_bismark_zip_update_meth_call (vb, range, range_len, idx, (bisulfite == base), M_i, Mseg_len, i);

        idx = RR_IDX(idx+1); // round robin around contig (useful in case of circular contigs)
    }

    COPY_TIMER(sam_seg_bisulfite_M);
}

// Creates a bitmap from seq data - exactly one bit per base that consumes reference (i.e. not bases with CIGAR I/S)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF (unless we can use an aligner)
// - Edge case: no CIGAR (it is "*") - we just store the sequence in SAM_NONREF
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
static MappingType sam_seg_SEQ_vs_ref (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(seq), const SamPosType pos, bool is_revcomp,
                                       uint32_t ref_consumed, uint32_t ref_and_seq_consumed, bool no_lock)
{
    START_TIMER;

    Context *bitmap_ctx = CTX(SAM_SQBITMAP);
    Context *nonref_ctx = CTX(SAM_NONREF);
    Context *seqmis_ctx = CTX(SAM_SEQMIS_A); // 4 contexts

    BitsP bitmap = (BitsP)&bitmap_ctx->local;
    uint32_t bitmap_start = bitmap_ctx->next_local;        

    ASSERTW (seq_len < 100000 || segconf.running || segconf.is_long_reads, 
                "%s: Warning: sam_seg_SEQ: seq_len=%u is suspiciously high and might indicate a bug", LN_NAME, seq_len);

    // we don't need to lock if the entire ref_consumed of this read was already is_set by previous reads (speed optimization)
    no_lock = IS_REF_EXTERNAL || ((vb->chrom_node_index == vb->consec_is_set_chrom) && (pos >= vb->consec_is_set_pos) && (pos + ref_consumed <= vb->consec_is_set_pos + vb->consec_is_set_len));

    RefLock lock = REFLOCK_NONE;
    Range *range = vb->cigar_missing ? NULL : ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, seq, (no_lock ? NULL : &lock));

    // Cases where we don't consider the refernce and just copy the seq as-is
    // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
    // 2. (loaded:) case contig doesn't exist in the reference, or POS is out of range of contig (observed in the wild with chrM)    
    if (!range) {
        COPY_TIMER (sam_seg_SEQ_vs_ref);
        return MAPPING_NO_MAPPING; // seg as a verbatim copy or use aligner
    }

    int64_t missing_bits = (int64_t)ref_and_seq_consumed - ((int64_t)bitmap_ctx->local.nbits - (int64_t)bitmap_ctx->next_local);
    if (missing_bits > 0) {
        buf_alloc_do (VB, &bitmap_ctx->local, roundup_bits2bytes64 (bitmap_ctx->local.nbits + seq_len), CTX_GROWTH, __FUNCLINE, NULL); 
        buf_extend_bits (&bitmap_ctx->local, missing_bits);
    }
    
    bits_set_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // we initialize all the bits to "set", and clear as needed.
    
    if (vb->bisulfite_strand) {
        if (segconf.sam_predict_meth_call) {
            buf_alloc_exact (vb, vb->meth_call, seq_len, char, "meth_call");
            memset (vb->meth_call.data, '.', vb->meth_call.len32);
        }

        if (segconf.MD_NM_by_unconverted) {
            buf_alloc_bits (vb, &vb->unconverted_bitmap, ref_and_seq_consumed, "unconverted_bitmap");
            bits_set_all ((BitsP)&vb->unconverted_bitmap); // we initialize all the bits to "set", and clear as needed.
            vb->unconverted_bitmap.nbits = 0; // we will grow it back to ref_and_seq_consumed as we go along
        }
    }

    for (int i=0; i < 4; i++)
        buf_alloc (vb, &seqmis_ctx[i].local, ref_and_seq_consumed, 0, char, CTX_GROWTH, "contexts->local"); 

    buf_alloc (vb, &nonref_ctx->local, seq_len + 3, 0, uint8_t, CTX_GROWTH, "contexts->local"); 

    bitmap_ctx->local_num_words++;

    uint32_t pos_index = pos - range->first_pos;
    uint32_t next_ref  = pos_index;

    uint32_t i=0, n=0;
    BamCigarOpType op = BC_NONE;

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    
    if (IS_REF_EXT_STORE && !no_lock) 
        bits_set_region (&range->is_set, pos_index, ref_consumed); // we will need this ref to reconstruct

    bool has_D_N = false;
    while (i < seq_len || next_ref < pos_index + ref_consumed) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_consumed, "i or next_ref are out of range");

        if (vb->binary_cigar.next < vb->binary_cigar.len) {
            BamCigarOp *next_op = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next++);
            op = next_op->op;
            n  = next_op->n;
        }

        switch (op) {

            case BC_M: case BC_E: case BC_X: // alignment match or sequence match or mismatch

                ASSERT (n > 0 && n <= (seq_len - i), 
                        "%s: CIGAR implies seq_len longer than actual seq_len=%u (n=%u ). CIGAR=\"%s\"", 
                        LN_NAME, seq_len, n, vb->last_cigar ? vb->last_cigar : "");

                uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
                uint32_t save_n = n;

                if (vb->bisulfite_strand) {
                    sam_seg_bisulfite_M (vb, range, range_len, next_ref, &seq[i], n, i, bitmap, bit_i, seqmis_ctx);
                    next_ref += n;
                    i += n;
                    n = 0;
                }
                
                else while (n && next_ref < pos_index + ref_consumed) {

                    bool normal_base = IS_NUCLEOTIDE (seq[i]);

                    // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                    uint32_t idx = (next_ref >= range_len) ? next_ref - range_len : next_ref; 

                    // case: we have not yet set a value for this site - we set it now. note: in ZIP, is_set means that the site
                    // will be needed for pizzing. With REF_INTERNAL, this is equivalent to saying we have set the ref value for the site
                    if (IS_REF_INTERNAL && !no_lock && !ref_is_nucleotide_set (range, idx)) { 

                        // note: in case this is a non-normal base (eg N), set the reference to an arbitrarily to 'A' as we 
                        // we will store this non-normal base in seqmis_ctx multiplexed by the reference base (i.e. in seqmis_ctx['A']).
                        ref_set_nucleotide (range, idx, normal_base ? seq[i] : 'A');
                        
                        bits_set (&range->is_set, idx); // we will need this ref to reconstruct

                        if (normal_base) 
                            bit_i++; 
                        else
                            goto mismatch;
                    }

                    // case our seq is identical to the reference at this site
                    else if (normal_base && seq[i] == REF(idx)) 
                        bit_i++; // bit remains set. 

                     // case: ref is set to a different value - we store our value in seqmis_ctx
                    else mismatch: {
                        uint8_t ref_base_2bit = bits_get2 (&range->ref, idx * 2);

                        BNXTc (seqmis_ctx[ref_base_2bit].local) = seq[i];
                        bits_clear (bitmap, bit_i++);
                        vb->mismatch_bases_by_SEQ++; // used by NM:i and MD:Z
                    } 
  
                    n--;
                    next_ref++;
                    i++;
                }
                
                ref_and_seq_consumed   -= save_n - n; // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap
                bitmap_ctx->next_local += save_n - n;

                break; // end if 'M', '=', 'X'

            // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
            case BC_I: case BC_S: 
                ASSSEG (n > 0 && n <= (seq_len - i), seq,
                        "CIGAR %s implies seq_len longer than actual seq_len=%u", vb->last_cigar ? vb->last_cigar : "", seq_len);

                buf_add (&nonref_ctx->local, &seq[i], n);
                i += n;
                n = 0;
                break;

            // for Deletion or Skipping - we move the next_ref ahead
            case BC_D: case BC_N: {
                unsigned ref_consumed_skip = (IS_REF_INTERNAL ? MIN_(n, range_len - next_ref)
                                                                             : n);
                next_ref += ref_consumed_skip;
                n        -= ref_consumed_skip;
                has_D_N  =  true;
                break;
            }

            // Hard clippping (H) or padding (P) - nothing much to do
            case BC_H: case BC_P: 
                n = 0;
                break;

            case BC_NONE:
                ASSSEG (false, vb->last_cigar, "End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                        "(cigar=%s pos=%d seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u op.n=%u range=[%.*s %"PRId64"-%"PRId64"])",
                        pos_index + ref_consumed - next_ref, seq_len-i, vb->last_cigar, pos, seq_len,
                        ref_consumed, next_ref, pos_index, n, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            default:
                ASSSEG (false, vb->last_cigar, "Invalid CIGAR op=%u", op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_consumed && n) break;
    }

    ref_unlock (gref, &lock); // does nothing if REFLOCK_NONE

    // an error here can indicate that CIGAR is inconsistent with sequence
    ASSERT (i == seq_len, "%s: expecting i(%u) == seq_len(%u) pos=%d range=[%.*s %"PRId64"-%"PRId64"] (possibly reason: inconsistency between seq_len and CIGAR=\"%s\")", 
            LN_NAME, i, seq_len, pos, STRf(range->chrom_name), range->first_pos, range->last_pos, (vb->last_cigar ? vb->last_cigar : ""));

    // if we set the entire consecutive reference range covered by this read - extend consec_is_set 
    if (segconf.is_sorted && !no_lock && (IS_REF_EXT_STORE || !has_D_N/*REF_INTERNAL*/)) {
        // case: current region is not consecutive - start a new region
        if (vb->consec_is_set_chrom != vb->chrom_node_index || vb->consec_is_set_pos + vb->consec_is_set_len < pos) {
            vb->consec_is_set_chrom = vb->chrom_node_index;
            vb->consec_is_set_pos   = pos;
            vb->consec_is_set_len   = ref_consumed;
        }

        // case: extend current region
        else if (vb->consec_is_set_pos + vb->consec_is_set_len < pos + ref_consumed)
            vb->consec_is_set_len = (pos + ref_consumed) - vb->consec_is_set_pos;
    }
    
    // verify - does MD:Z correctly reflect matches and mismatches of M/X/=
    bool use_un = segconf.MD_NM_by_unconverted && vb->bisulfite_strand;
    sam_MD_Z_verify_due_to_seq (vb, STRa(seq), pos, 
                                use_un ? (BitsP)&vb->unconverted_bitmap : bitmap, 
                                use_un ? 0                              : bitmap_start);

    // note: change of logic in v14: up to v13, we use to align every recursion level, it seems that this was a bug
    sam_seg_SEQ_pad_nonref (vb);

    COPY_TIMER (sam_seg_SEQ_vs_ref);

    if (vb->mismatch_bases_by_SEQ != 0 || vb->bisulfite_strand) // we don't support MAPPING_PERFECT with bisulfite (bug 649)
        return MAPPING_ALIGNED;
    
    else {
        bitmap_ctx->next_local -= vb->ref_and_seq_consumed; // note: we truncate bitmap, if needed, in sam_seg_finalize
        return MAPPING_PERFECT;
    } 
}

// verify that the depn line is identical (modulo rev_comp and hard_clips) to the prim line
static bool sam_seg_verify_saggy_line_SEQ (VBlockSAMP vb, ZipDataLineSAM *my_dl, ZipDataLineSAM *prim_dl, rom my_seq/*textual*/)
{
    START_TIMER;

    bool xstrand = (my_dl->FLAG.rev_comp != prim_dl->FLAG.rev_comp); // primary and dependent are on opposite strands
    uint32_t seq_len = my_dl->SEQ.len;

    rom prim_seq;
     // TO DO: compare BAM native-to-native without converting to textual
    if (IS_BAM_ZIP) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc (vb, &vb->scratch, 0, prim_dl->SEQ.len, char, 0, "scratch");
        bam_seq_to_sam (vb, (bytes)STRtxtw (prim_dl->SEQ), false, false, &vb->scratch);
        prim_seq = Bc(vb->scratch, vb->hard_clip[xstrand]);
    }
    else
        prim_seq = Btxt (prim_dl->SEQ.index) + vb->hard_clip[xstrand];

    bool success = true;
    if (!xstrand)
        success = !memcmp (my_seq, prim_seq, seq_len);

    else  {
        rom my_next = my_seq;
        rom prim_next = prim_seq + seq_len - 1;
        for (uint32_t i=0; i < seq_len; i++)
            if (*my_next++ != COMPLEM[(uint8_t)*prim_next--]) {
                success = false;
                break;
            }
    }

    if (IS_BAM_ZIP)
        buf_free (vb->scratch);

    COPY_TIMER(sam_seg_verify_saggy_line_SEQ);
    return success;
}

void sam_seg_SEQ (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(textual_seq), unsigned add_bytes)
{
    START_TIMER;

    bool force_sequence = false; // true if reconstructor would normally diff, but in this case it should reconstruct from SQBITMAP.local
    bool force_verbatim = false; // true if reconstructor would normally be vs_ref, but it verbatim
    bool aligner_used = false;
    bool perfect = false;
    bool vs_prim = false;
                            
    ZipDataLineSAM *saggy_dl;
    bool unmapped = dl->FLAG.unmapped || vb->cigar_missing || !dl->POS || str_issame_(STRa(vb->chrom_name), "*", 1);

    // case segconf: we created the contexts for segconf_set_vb_size accounting. that's enough - actually segging will mark is_set and break sam_seg_MD_Z_analyze.
    if (segconf.running) {
        segconf.longest_seq_len = MAX_(segconf.longest_seq_len, textual_seq_len);

        if (!vb->line_i) segconf_mark_as_used (VB, 5, SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND);

        if ((vb->bisulfite_strand=='G') != dl->FLAG.rev_comp)
            segconf.bs_strand_not_by_rev_comp = true; // note: in files that bisulfite_strand is determined by rev_comp, it is true for all lines

        return; 
    }

    // case: unmapped line and we have refhash: align to reference
    if (unmapped && flag.aligner_available) {
        use_aligner:
        switch (aligner_seg_seq (VB, CTX(SAM_SQBITMAP), STRa(textual_seq), true, false, NO_GPOS, false)) {
            case MAPPING_NO_MAPPING : force_verbatim = true; goto add_seq_verbatim; 
            case MAPPING_PERFECT    : perfect        = true; // fallthrough
            case MAPPING_ALIGNED    : aligner_used   = true; break;
            default                 : ABORT0 ("bad value");
        }

        sam_seg_SEQ_pad_nonref (vb);
        vb->md_verified = false;    
        vb->mismatch_bases_by_SEQ = -1; // we can't seg NM:i with special
    }

    // case: unmapped line, no refhash: just store the sequence in nonref without an indication in the bitmap
    else if (unmapped && !flag.aligner_available) { 
        add_seq_verbatim:        
        buf_add_more (VB, &CTX(SAM_NONREF)->local, STRa(textual_seq), "contexts->local");
        sam_seg_SEQ_pad_nonref (vb);
        vb->md_verified = false;    
        vb->mismatch_bases_by_SEQ = -1; 
        force_verbatim = true;    
    }

    // case: no SEQ. we already add '-' to CIGAR - no data added here
    else if (textual_seq[0] == '*') {
        vb->md_verified = false;    
        vb->mismatch_bases_by_SEQ = -1; 
        vb->seq_missing = dl->no_seq = true;
    }

    // in case of DEPN line with SEQ confirmed by sam_seg_depn_is_subseq_of_prim to be a subset of the sag sequence, 
    else if (sam_is_depn_vb && vb->sag) // depn 
        vs_prim = true;

    // case: DEPN line (or STAR line - see sam_seg_buddy) vs same-VB prim
    else if (({ saggy_dl = DATA_LINE (vb->saggy_line_i) ; sam_has_saggy && !saggy_dl->hard_clip[0] && !saggy_dl->hard_clip[1]; })) {
        
        // note: these conditions must generate "force_sequence" because PIZ cannot check for them
        if (saggy_dl->no_seq || 
            (textual_seq_len + vb->hard_clip[0] + vb->hard_clip[1] != saggy_dl->SEQ.len) ||
            !sam_seg_verify_saggy_line_SEQ (vb, dl, saggy_dl, textual_seq)) {

            force_sequence = true;    
            goto vs_sequence;
        }
        
        vs_prim = true;
    }
    
    else vs_sequence: 
        switch (sam_seg_SEQ_vs_ref (vb, dl, STRa(textual_seq), dl->POS, dl->FLAG.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed, false)) {
            case MAPPING_NO_MAPPING : if (flag.aligner_available) goto use_aligner;
                                      else { force_verbatim = true; goto add_seq_verbatim; };
            case MAPPING_PERFECT    : perfect = true; break;
            default                 : break;
        }

    // set vb->md_verified and vb->mismatch_bases_by_SEQ needed by MD:Z and NM:i
    bool force_no_analyze_depn_SEQ = false;
    if (vs_prim && segconf.has_MD_or_NM && !vb->cigar_missing && 
        !sam_analyze_copied_SEQ (vb, STRa(textual_seq), dl->POS, dl->FLAG.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed, &vb->codec_bufs[0]))
            force_no_analyze_depn_SEQ = true;

    // case: PRIM - SQBITMAP context is reconstructed when loading SA Groups (preprocessing) - invoking SPECIAL_SEQ and reconstructing the sequence
    //       and SQBITMAP is reconstructed again during recon - but this SPECIAL_SEQ copies from SA Groups
    // case: MAIN/DEPN - SPECIAL_SEQ is invoked by SQBITMAP - it either reconstructs sequence or diffs vs prim.
    //       (prim SEQ may be taken from SA Groups or from earlier in the current VB (the latter only in MAIN VBs)
    #define SEQ_SNIP_FORCE_SEQ           0
    #define SEQ_SNIP_ALIGNER_USER        1
    #define SEQ_SNIP_PERFECT             2
    #define SEQ_SNIP_NO_ANALYZE_DEPN_SEQ 3
    #define SEQ_SNIP_BISULFATE           4
    #define SEQ_FORCE_VERBATIM           5

    seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SEQ, '0'+force_sequence, '0'+aligner_used, '0'+perfect, '0'+force_no_analyze_depn_SEQ, 
                                (!vb->bisulfite_strand) ? '0'
                              : (vb->bisulfite_strand && !segconf.bs_strand_not_by_rev_comp && (vb->bisulfite_strand=='G') == dl->FLAG.rev_comp) ? '*' // bisulfite_strand is predicted by rev_comp
                              :                           vb->bisulfite_strand,
                              '0'+force_verbatim },
                8, SAM_SQBITMAP, add_bytes); 

    COPY_TIMER (sam_seg_SEQ);
}

// converts native SAM/BAM format to 2bit ACGT - if soft_fail, returns false if any base is not A,C,G or T 
bool sam_sa_native_to_acgt (VBlockSAMP vb, Bits *packed, uint64_t next_bit, STRp(seq), bool bam_format, bool revcomp, bool soft_fail)
{
    if (bam_format) {
        if (!revcomp)
            for (uint32_t i=0; i < seq_len; i++, next_bit += 2) {
                uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
                switch (b) {
                    case 0b0001 : bits_assign2 (packed, next_bit, 0); break;
                    case 0b0010 : bits_assign2 (packed, next_bit, 1); break;
                    case 0b0100 : bits_assign2 (packed, next_bit, 2); break;
                    case 0b1000 : bits_assign2 (packed, next_bit, 3); break;
                    default     : if (soft_fail) return false;
                                  ABORT ("%s: Unexpected base: '%c' i=%u seq_len=%u", LN_NAME, bam_base_codes[b], i, seq_len);
                }
            }    
        else
            for (int32_t i=seq_len-1; i >= 0; i--, next_bit += 2) {
                uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
                switch (b) {
                    case 0b0001 : bits_assign2 (packed, next_bit, 3); break;
                    case 0b0010 : bits_assign2 (packed, next_bit, 2); break;
                    case 0b0100 : bits_assign2 (packed, next_bit, 1); break;
                    case 0b1000 : bits_assign2 (packed, next_bit, 0); break;
                    default     : if (soft_fail) return false;
                                  ABORT ("%s: Unexpected base: '%c' i=%u seq_len=%u", LN_NAME, bam_base_codes[b], i, seq_len);
                }
            }    
    }

    else { // SAM
        if (!revcomp)
            for (uint32_t i=0; i < seq_len; i++, next_bit += 2) 
                switch (seq[i]) {
                    case 'A' : bits_assign2 (packed, next_bit, 0); break;
                    case 'C' : bits_assign2 (packed, next_bit, 1); break;
                    case 'G' : bits_assign2 (packed, next_bit, 2); break;
                    case 'T' : bits_assign2 (packed, next_bit, 3); break;
                    default  : if (soft_fail) return false;
                               ABORT ("%s: Unexpected base: '%c'(ASCII %u) i=%u seq_len=%u", LN_NAME, seq[i], (uint8_t)seq[i], i, seq_len);
                }
        else
            for (int32_t i=seq_len-1; i >= 0; i--, next_bit += 2) 
                switch (seq[i]) {
                    case 'A' : bits_assign2 (packed, next_bit, 3); break;
                    case 'C' : bits_assign2 (packed, next_bit, 2); break;
                    case 'G' : bits_assign2 (packed, next_bit, 1); break;
                    case 'T' : bits_assign2 (packed, next_bit, 0); break;
                    default  : if (soft_fail) return false;
                               ABORT ("%s: Unexpected base: '%c'(ASCII %u) i=%u seq_len=%u", LN_NAME, seq[i], (uint8_t)seq[i], i, seq_len);
                }
    }

    bits_clear_excess_bits_in_top_word (packed);

    return true;
}

// setting ref bases by analyze_* functions
rom ERR_ANALYZE_RANGE_NOT_AVAILABLE = "Range not available";
rom ERR_ANALYZE_DEPN_NOT_IN_REF     = "Depn alignment base not is_set in reference";
rom ERR_ANALYZE_INCORRECT_REF_BASE  = "incorrect reference base";
rom sam_seg_analyze_set_one_ref_base (VBlockSAMP vb, bool is_depn, SamPosType pos, char base, 
                                      uint32_t ref_consumed, // remaining ref_consumed starting at pos                                      
                                      RangeP *range_p, RefLock *lock)
{
    // case: pos is beyond the existing range
    if ((*range_p) && (*range_p)->last_pos < pos) {
        ref_unlock (gref, lock);
        *range_p = NULL;
    }

    // get range (and lock it if needed)
    if (! *range_p) {
        *range_p = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, NULL, 
                                      (IS_REF_EXTERNAL || is_depn) ? NULL : lock);
        if (! *range_p) return ERR_ANALYZE_RANGE_NOT_AVAILABLE; // cannot access this range in the reference
    }

    uint32_t pos_index = pos - (*range_p)->first_pos; // index within range

    // case: depn line, but we haven't set the reference data for it. since we don't need the reference data for reconstructing
    // SEQ (it is copied from prim), we don't store it just for MD:Z - so we won't use SPECIAL for this MD:Z
    if (is_depn && (flag.reference & REF_STORED) && !ref_is_nucleotide_set (*range_p, pos_index))
        return ERR_ANALYZE_DEPN_NOT_IN_REF;

    bool internal_pos_is_populated = IS_REF_INTERNAL && ref_is_nucleotide_set (*range_p, pos_index);

    // case: reference already contains a base - but unfortunately it is not "base" - so MD is not reconstractable from reference
    if (((flag.reference & REF_ZIP_LOADED) || internal_pos_is_populated) && (base != ref_base_by_pos (*range_p, pos))) 
        return ERR_ANALYZE_INCORRECT_REF_BASE; // encountered in the wild when the reference base is a IUPAC

    // case: reference is not set yet - set it now (extra careful never to set anything in depn, as range is not locked)
    if (!is_depn && IS_REF_INTERNAL && !internal_pos_is_populated)
        ref_set_nucleotide (*range_p, pos_index, base);

    // set is_set - we will need this base in the reference to reconstruct MD
    if (!is_depn && flag.reference & REF_STORED)
        bits_set (&(*range_p)->is_set, pos_index); // we will need this ref to reconstruct

    return NULL; // success
}

// pack seq (into 2bit ACGT format) of each PRIM line separately 
void sam_zip_prim_ingest_vb_pack_seq (VBlockSAMP vb, Sag *vb_grps, uint32_t vb_grps_len,
                                      BufferP underlying_buf, BufferP packed_seq_buf, bool is_bam_format)
{
    uint32_t total_seq_len=0;
    for (uint32_t grp_i=0; grp_i < vb_grps_len; grp_i++) 
        total_seq_len += vb_grps[grp_i].seq_len;

    // allocate memory but don't extend the bitmap yet
    underlying_buf->len = roundup_bits2bytes64 (total_seq_len * 2);
    buf_alloc (vb, underlying_buf, underlying_buf->len, 0, char, 0, "z_data");
    buf_set_overlayable (underlying_buf);
    buf_overlay (vb, packed_seq_buf, underlying_buf, "packed_seq_buf");

    for (uint32_t vb_grp_i=0; vb_grp_i < vb_grps_len; vb_grp_i++) {
        
        Sag *vb_grp = &vb_grps[vb_grp_i];

        uint64_t next_bit = packed_seq_buf->nbits;
        buf_extend_bits (packed_seq_buf, vb_grp->seq_len * 2); // extend now 
        Bits *sag_seq = (BitsP)packed_seq_buf;

        sam_sa_native_to_acgt (vb, sag_seq, next_bit, Bc (vb->txt_data, vb_grp->seq), vb_grp->seq_len, is_bam_format, false, false);
        vb_grp->seq = next_bit / 2; // update from an index into txt_data to an index (bases not bits) into sag_seq
    }
}

// used by codec_longr_compress
COMPRESSOR_CALLBACK (sam_zip_seq) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    *line_data_len = dl->SEQ.len;

    if (VB_DT(SAM)) 
        *line_data = Btxt (dl->SEQ.index);
    else { // BAM
        VB_SAM->textual_seq.len = 0; // reset
        buf_alloc (vb, &VB_SAM->textual_seq, 0, dl->SEQ.len+1 /* +1 for last half-byte */, char, 1.5, "textual_seq");

        bam_seq_to_sam (VB_SAM, (uint8_t *)Btxt(dl->SEQ.index), dl->SEQ.len, false, false, &VB_SAM->textual_seq);
        *line_data = B1STc (VB_SAM->textual_seq);
    }

    if (is_rev) *is_rev = dl->FLAG.rev_comp;
}

//---------
// PIZ
//---------

// get a bytemap of ref_consumed values. returns NULL if no range exists. get 2 bases before and after
// in case needed for methylation calling.
static inline uint8_t *sam_reconstruct_SEQ_get_ref_bytemap (VBlockSAMP vb, ContextP bitmap_ctx, bool v14, SamPosType pos, bool predict_meth_call)
{
    // note: in an edge case, when all is_set bits are zero, the range might not even be written to the file
    bool uses_ref_data = vb->ref_consumed && 
                            (v14 || // v14: we always have is_set for all the ref data (mismatches use it for SEQMIS)
                             bits_num_set_bits_region ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local, vb->ref_and_seq_consumed) > 0);

    vb->range = uses_ref_data ? (RangeP)ref_piz_get_range (VB, gref, true) : NULL; 
    if (!vb->range) return NULL; 

    SamPosType range_len = vb->range ? (vb->range->last_pos - vb->range->first_pos + 1) : 0;
    uint32_t num_bases = vb->ref_consumed + (predict_meth_call ? 4 : 0); // 2 bases before and after in case needed for methylation

    ASSERTNOTINUSE (vb->scratch);
    ARRAY_alloc (uint8_t, ref, num_bases, false, vb->scratch, vb, "scratch");

#ifdef DEBUG
    ASSERTNOTINUSE (vb->piz_is_set);
    ARRAY_alloc (uint8_t, is_set, num_bases, false, vb->piz_is_set, vb, "codec_bufs[0]");
#endif

    int32_t idx = RR_IDX ((int32_t)(pos - vb->range->first_pos) - (predict_meth_call ? 2 : 0)); // two bases before pos in case needed for methylation 

    uint32_t num_ref_bases = MIN_(num_bases, range_len - idx);
    bits_base_to_byte (ref, &vb->range->ref, idx, num_ref_bases); // entries with is_set=0 will be garbage
#ifdef DEBUG
    if (!IS_REF_EXTERNAL)
        bits_bit_to_byte (is_set, &vb->range->is_set, idx, num_ref_bases); 
#endif
 
    // if ref_consumed goes beyond end of range, take the rest from the beginning of range (i.e. circling around)
    if (num_ref_bases < num_bases) {
        bits_base_to_byte (&ref[num_ref_bases], &vb->range->ref, 0, num_bases - num_ref_bases); 
#ifdef DEBUG
        if (!IS_REF_EXTERNAL)
            bits_bit_to_byte (&is_set[num_ref_bases], &vb->range->is_set, 0, num_bases - num_ref_bases); 
#endif
    }
    
    return ref + (predict_meth_call ? 2 : 0); // start of ref consumed
}

// PIZ: SEQ reconstruction 
void sam_reconstruct_SEQ_vs_ref (VBlockP vb_, ContextP bitmap_ctx, STRp(snip), bool reconstruct)
{
    START_TIMER;

    #ifdef DEBUG
        #define verify_is_set(base_i) ASSPIZ (IS_REF_EXTERNAL || *Bc(vb->piz_is_set, (base_i) + (predict_meth_call ? 2 : 0)) == 1, "Expecting POS=%u + base_i=%u to have is_set=1", (SamPosType)vb->contexts[SAM_POS].last_value.i, (base_i))
    #else
        #define verify_is_set(base_i) 
    #endif

    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (!bitmap_ctx->is_loaded) return; // if case we need to skip the SEQ field (for the entire VB)

    bool v14              = VER(14);
    Context *nonref_ctx   = CTX(SAM_NONREF);
    Context *seqmis_ctx   = CTX(SAM_SEQMIS_A);
    rom nonref            = Bc (nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    rom nonref_start      = nonref;
    const SamPosType pos  = vb->last_int(SAM_POS);
    unsigned seq_consumed=0, ref_consumed=0;

    // note: prim alignments (loaded in preprocessing) are always mapped - enforced by sam_seg_is_gc_line
    bool unmapped = v14 ? (!vb->preprocessing && (last_flags.unmapped || vb->cigar_missing || !pos || IS_ASTERISK (vb->chrom_name)))
                        : (!pos || IS_ASTERISK (vb->chrom_name)); // the criterion seg used in up to v13
                                                  
    bool aligner_used = v14 ? (snip && snip[SEQ_SNIP_ALIGNER_USER] == '1') // starting v14, seg tells us explicitly
                            : (unmapped && z_file->z_flags.aligner); 

    bool is_perfect = v14 && snip && snip[SEQ_SNIP_PERFECT] == '1';

    bool force_verbatim = v14 && snip && snip[SEQ_FORCE_VERBATIM] == '1';

    char bisulfite = vb->bisulfite_strand; // automatic var for efficiency in tight loop
    bool predict_meth_call = bisulfite && segconf.sam_predict_meth_call;

    bool MD_NM_by_unconverted = v14 && bisulfite && segconf.MD_NM_by_unconverted; // we calculate MD:Z and NM:i vs unconverted reference

    if (predict_meth_call && vb->ref_and_seq_consumed) {
        buf_alloc_exact (vb, vb->meth_call, vb->seq_len, char, "meth_call");
        memset (vb->meth_call.data, '.', vb->meth_call.len32);
    }
    
    // case: unmapped, segged against reference using our aligner
    if (aligner_used) {
        aligner_reconstruct_seq (VB, bitmap_ctx, vb->seq_len, false, is_perfect, reconstruct);
        nonref_ctx->next_local = ROUNDUP4 (nonref_ctx->next_local);
        return; // aligner accounts for time separately, so we don't account for the time here
    }

    // case: unmapped, and copied verbatim
    if (unmapped || force_verbatim) unmapped: {
        if (reconstruct) RECONSTRUCT (nonref, vb->seq_len); 
        nonref_ctx->next_local += ROUNDUP4 (vb->seq_len);
        goto done;
    }

    // case: missing sequence - sequence is '*' - just reconstruct the '*' (but only if it is the primary SEQ field,
    // which we can tell by not being encountered earlier on the line - bc in v8 (at least) E2 was also stored in the same SAM_SQBITMAP)
    if (vb->seq_missing && !ctx_encountered_in_line (VB, bitmap_ctx->did_i)) {
        if (reconstruct) RECONSTRUCT1 ('*');
        goto done;
    }

    ASSERT0 (is_perfect || bitmap_ctx->local.len32 || nonref_ctx->local.len32, "No SEQ data, perhaps sections were skipped?");

    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence

    uint8_t *ref = sam_reconstruct_SEQ_get_ref_bytemap (vb, bitmap_ctx, v14, pos, predict_meth_call);
    if (v14 && !ref) goto unmapped; // starting v14, a missing range is another case of unmapped

    const BamCigarOp *cigar = vb->binary_cigar.len32 ? B1ST(BamCigarOp, vb->binary_cigar) 
                                                     : &(BamCigarOp){ .n = vb->seq_len, .op = BC_S }; // up to v13, alignments with CIGAR=* but with a sequence where segged this way
                                                     
    BamCigarOp op = {};
    bool consumes_query=false, consumes_reference=false;
    uint32_t save_n;
    char *recon = BAFTtxt;

    while (seq_consumed < vb->seq_len || ref_consumed < vb->ref_consumed) {
        
        if (!op.n) {
            op = *cigar++;
            save_n = op.n;

            switch (op.op) {
                case BC_M : case BC_X : case BC_E : consumes_query=true;  consumes_reference=true;  break;
                case BC_I : case BC_S :             consumes_query=true;  consumes_reference=false; break;
                case BC_D : case BC_N :             consumes_query=false; consumes_reference=true;  break;
                case BC_H : case BC_P :             consumes_query=false; consumes_reference=false; break;
                default   : ASSPIZ (false, "Invalid op=%d", op.op);
            }
        }

        // shortcut if perfect (= no mismatches)
        if (is_perfect && consumes_query && consumes_reference) {
            if (reconstruct) 
                for (int i=0; i < op.n; i++) {
                    verify_is_set (ref_consumed + i);
                    *recon++ = acgt_decode (ref[ref_consumed + i]);
                }

            seq_consumed       += op.n;
            ref_consumed       += op.n;
            op.n = 0;
        }
        
        else {
            if (consumes_query) {
            
                if (consumes_reference) {
                            
                    if (NEXTLOCALBIT (bitmap_ctx)) { // copy from reference (or from converted reference if bisulfite)
                        verify_is_set (ref_consumed);
                        if (reconstruct) {
                            char ref_base = acgt_decode(ref[ref_consumed]);
                            
                            // case: bisulfite data, we segged against the converted reference (C->T or G->A) - converted bases represent unmethylated bases (methylated bases remain unconverted)
                            if (ref_base == bisulfite) {
                                ref_base = (bisulfite == 'C') ? 'T' : 'A'; // convert reference base C->T or G->A
                                
                                if (MD_NM_by_unconverted) {
                                    vb->mismatch_bases_by_SEQ++; // this is actually a mismatch vs the unconverted reference
                                    bits_clear ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local-1); // adjust sqbitmap for use for reconstructing MD:Z
                                }

                                if (predict_meth_call)
                                    sam_bismark_piz_update_meth_call (vb, ref, ref_consumed, seq_consumed, save_n-op.n, save_n, cigar, bisulfite, false); // converted therefore unmethylated
                            }

                            *recon++ = ref_base;
                        } 
                    }
                    else {
                        vb->mismatch_bases_by_SEQ++;

                        if (v14) {
                            verify_is_set (ref_consumed); 
                            // note: a "not enough data" error here is often an indication that pos is wrong
                            ContextP ctx = &seqmis_ctx[ref[ref_consumed]];
                            char base = *Bc(ctx->local, ctx->next_local++);
                            *recon++ = base;

                            if (bisulfite && base == acgt_decode(ref[ref_consumed])) { // uncoverted base = metyhlated

                                if (segconf.sam_predict_meth_call)
                                    sam_bismark_piz_update_meth_call (vb, ref, ref_consumed, seq_consumed, save_n-op.n, save_n, cigar, bisulfite, true); 

                                if (MD_NM_by_unconverted) {
                                    vb->mismatch_bases_by_SEQ--; // this is actually a match vs the unconverted referenc, undo the mismatch counter
                                    bits_set ((BitsP)&bitmap_ctx->local, bitmap_ctx->next_local-1); // adjust sqbitmap for use for reconstructing MD:Z                        
                                }
                            }
                        }
                        else
                            goto nonref;
                    }
                }
                                
                else nonref: {
                    if (reconstruct) *recon++ = (*nonref);
                    nonref++;
                }

                seq_consumed++;
            }

            if (consumes_reference) 
                ref_consumed++;

            op.n--;
        }
    }

    vb->txt_data.len32 += vb->seq_len;

    ASSPIZ (seq_consumed == vb->seq_len,      "expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSPIZ (ref_consumed == vb->ref_consumed, "expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    nonref_ctx->next_local += ROUNDUP4 (nonref - nonref_start);

    buf_free (vb->scratch); // allocated by sam_reconstruct_SEQ_get_ref_bytemap
#ifdef DEBUG
    buf_free (vb->piz_is_set); 
#endif

done:
    COPY_TIMER (sam_reconstruct_SEQ_vs_ref);
}


// reconstruct from a 2bit array - start_base and seq_len are in bases (not in bits)
static void reconstruct_SEQ_acgt (VBlockSAMP vb, BitsP seq_2bits, uint64_t start_base, uint32_t seq_len, bool revcomp)
{
    // Reconstruct as SAM. If BAM, we later translate to BAM. To do: reconstruct directly as BAM if needed, bug 530
    char *next = BAFTc(vb->txt_data);
    
    if (!revcomp)
        for (uint32_t i=0; i < seq_len; i++) {
            uint8_t b = bits_get2 (seq_2bits, (start_base + i)*2);
            *next++ = b==0 ? 'A' : b==1 ? 'C' : b==2 ? 'G' : 'T';
        }
    
    else
        for (int32_t i=seq_len-1; i >= 0; i--) {
            uint8_t b = bits_get2 (seq_2bits, (start_base + i)*2);
            *next++ = b==3 ? 'A' : b==2 ? 'C' : b==1 ? 'G' : 'T';
        }

    vb->txt_data.len32 += seq_len;
}

static void reconstruct_SEQ_copy_sag_prim (VBlockSAMP vb, ContextP ctx, BitsP prim, uint64_t prim_start_base, uint32_t prim_seq_len, bool xstrand, bool reconstruct, bool force_no_analyze_depn_SEQ)
{
    START_TIMER;
    
    rom seq = BAFTtxt;
    uint32_t depn_seq_len = prim_seq_len - vb->hard_clip[0] - vb->hard_clip[1]; // depn sequence is a sub-sequence of the prim sequence
    uint64_t depn_start_base = prim_start_base + (xstrand ? (prim_seq_len - depn_seq_len - vb->hard_clip[0]) 
                                                          : vb->hard_clip[0]);
    if (reconstruct && ctx->is_loaded) {
        reconstruct_SEQ_acgt (vb, prim, depn_start_base, depn_seq_len, xstrand);

        // in needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i
        // NOTE: if SEQ is not loaded or !reconstruct - it is expected that NM and MD also don't reconstruct
        if (segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing && !force_no_analyze_depn_SEQ) 
            sam_analyze_copied_SEQ (vb, seq, depn_seq_len, CTX(SAM_POS)->last_value.i, last_flags.rev_comp,
                                    vb->ref_consumed, vb->ref_and_seq_consumed, &ctx->line_sqbitmap);        
    }

    COPY_TIMER (reconstruct_SEQ_copy_sag_prim);
}

static void reconstruct_SEQ_copy_saggy (VBlockSAMP vb, ContextP ctx, bool reconstruct, bool force_no_analyze_depn_SEQ)
{
    if (!reconstruct || !ctx->is_loaded) return;

    SamFlags saggy_flags = (SamFlags){ .value = *B(int64_t, CTX(SAM_FLAG)->history, vb->saggy_line_i) };
    bool xstrand = saggy_flags.rev_comp != last_flags.rev_comp;

    HistoryWord word = *B(HistoryWord, ctx->history, vb->saggy_line_i); // SEQ is always stored as LookupTxtData or LookupPerLine
    uint32_t seq_len = word.len - vb->hard_clip[0] - vb->hard_clip[1]; // our sequence is a sub-sequence of the saggh sequence

    char *seq = BAFTtxt;

    rom saggy_seq = ((word.lookup == LookupTxtData) ? Btxt(word.index) : Bc(ctx->per_line, word.index));

    if (!IS_RECON_BAM) {
        if (xstrand)
            str_revcomp_in_out (seq, saggy_seq + vb->hard_clip[1], seq_len);
        else
            memcpy (seq, saggy_seq + vb->hard_clip[0], seq_len);
    
        vb->txt_data.len32 += seq_len;
    }

    else {
        uint32_t clip = vb->hard_clip[xstrand];
        bam_seq_to_sam (vb, (uint8_t*)saggy_seq + clip/2, seq_len, clip%2, false, &vb->txt_data);

        if (xstrand)
            str_revcomp_in_out (seq, seq, seq_len);
    }

    // in needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i
    // NOTE: if SEQ is not loaded or !reconstruct - it is expected that NM and MD also don't reconstruct
    if (segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing && !force_no_analyze_depn_SEQ) 
        sam_analyze_copied_SEQ (vb, STRa(seq), CTX(SAM_POS)->last_value.i, last_flags.rev_comp,
                                vb->ref_consumed, vb->ref_and_seq_consumed, &ctx->line_sqbitmap);        
}

// PIZ of a SEQ in a MAIN and DEPN components - this is an all-the-same context - SEQ all lines in a DEPN component
// come here. We multiplex between SAGroup-based reconstruction vs normal SEQ reconstruction
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_SEQ)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    vb->bisulfite_strand = (!VER(14) || snip[SEQ_SNIP_BISULFATE]=='0') ? 0
                         : snip[SEQ_SNIP_BISULFATE] == '*'             ? "CG"[last_flags.rev_comp]
                         :                                               snip[SEQ_SNIP_BISULFATE];
    // case: reconstructor would normally diff, but in this case it should reconstruct from SQBITMAP.local    
    if (snip[SEQ_SNIP_FORCE_SEQ] == '1') 
        goto vs_seq__aligner__verbatim;

    // case: PRIM component - copy from SA Group 
    else if (sam_is_prim_vb && !vb->preprocessing) {
        rom seq = BAFTtxt;
        reconstruct_SEQ_acgt (vb, (BitsP)&z_file->sag_seq, vb->sag->seq, vb->sag->seq_len, false);

        // if needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i
        if (segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing) 
            sam_analyze_copied_SEQ (vb, seq, vb->sag->seq_len, CTX(SAM_POS)->last_value.i, last_flags.rev_comp,
                                    vb->ref_consumed, vb->ref_and_seq_consumed, &ctx->line_sqbitmap);        
    }

    else if (snip[SEQ_SNIP_ALIGNER_USER] == '1' || snip[SEQ_FORCE_VERBATIM] == '1')
        goto vs_seq__aligner__verbatim;

    // case: DEPN component - line with sag - copy from sag loaded from prim VB
    else if (sam_is_depn_vb && SAM_PIZ_HAS_SAG) {
        const Sag *grp = vb->sag;
        ASSPIZ (grp->seq_len > vb->hard_clip[0] + vb->hard_clip[1], "grp_i=%u aln_i=%"PRIu64" grp->seq_len=%u <= hard_clip[LEFT]=%u + hard_clip[RIGHT]=%u", 
                ZGRP_I(grp), ZALN_I(vb->sa_aln), grp->seq_len, vb->hard_clip[0], vb->hard_clip[1]);

        reconstruct_SEQ_copy_sag_prim (vb, ctx, (BitsP)&z_file->sag_seq, STRa(grp->seq), 
                                       (last_flags.rev_comp != grp->revcomp), reconstruct, snip[SEQ_SNIP_NO_ANALYZE_DEPN_SEQ]-'0');
    }

    // case: copy from saggy line in this VB
    else if (sam_has_saggy && !B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->saggy_line_i)->hard_clip[0] 
                           && !B(CigarAnalItem, CTX(SAM_CIGAR)->cigar_anal_history, vb->saggy_line_i)->hard_clip[1])
        reconstruct_SEQ_copy_saggy (vb, ctx, reconstruct, snip[SEQ_SNIP_NO_ANALYZE_DEPN_SEQ]-'0');

    // case: reconstruct sequence directly
    else 
        vs_seq__aligner__verbatim:
        sam_reconstruct_SEQ_vs_ref (VB, ctx, STRa(snip), reconstruct);

    return NO_NEW_VALUE;
}

// SAM-to-BAM translator: translate SAM ASCII sequence characters to BAM's 4-bit characters:
TRANSLATOR_FUNC (sam_piz_sam2bam_SEQ)
{
    // the characters "=ACMGRSVTWYHKDBN" are mapped to BAM 0->15, in this matrix we add 0x80 as a validity bit. All other characters are 0x00 - invalid
    static const uint8_t sam2bam_seq_map[256] = { ['=']=0x80, ['A']=0x81, ['C']=0x82, ['M']=0x83, ['G']=0x84, ['R']=0x85, ['S']=0x86, ['V']=0x87, 
                                                  ['T']=0x88, ['W']=0x89, ['Y']=0x8a, ['H']=0x8b, ['K']=0x8c, ['D']=0x8d, ['B']=0x8e, ['N']=0x8f };

    BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)Bc (vb->txt_data, vb->line_start);
    uint32_t l_seq = LTEN32 (alignment->l_seq);

    // backward compatability note: prior to v14 the condition was:
    // if (CTX(SAM_QUAL)->lcodec == CODEC_LONGR) ...
    // Since v14, it is determined by a flag. Since this flag is 0 in V<=13, earlier files will always
    // store textual_seq, including if CODEC_LONGR

    // case downstream contexts need access to the textual_seq: copy it.
    // Examples: LONGR codec for QUAL; sam_piz_special_BSSEEKER2_XM
    if (!CTX(SAM_SQBITMAP)->flags.no_textual_seq) {
        buf_alloc (vb, &VB_SAM->textual_seq, 0, recon_len, char, 0, "textual_seq");
        memcpy (B1STc (VB_SAM->textual_seq), recon, recon_len);
        VB_SAM->textual_seq.len = recon_len;
    }

    if (vb->drop_curr_line) return 0; // sequence was not reconstructed - nothing to translate
    
    // if l_seq=0, just remove the '*'
    if (!l_seq) {
        vb->txt_data.len32--;
        return 0;
    }

    // if l_seq is odd, 0 the next byte that will be half of our last result byte
    if (l_seq % 2) *BAFTtxt = 0; 

    uint8_t *seq_before=(uint8_t *)recon, *seq_after=(uint8_t *)recon; 
    for (uint32_t i=0; i < (l_seq+1)/2; i++, seq_after++, seq_before += 2) {
        uint8_t base[2] = { sam2bam_seq_map[(uint8_t)seq_before[0]], sam2bam_seq_map[(uint8_t)seq_before[1]] };
        
        // check for invalid characters - issue warning (only once per execution), and make then into an 'N'
        for (unsigned b=0; b < 2; b++)
            if (!base[b] && !(b==1 && (i+1)*2 > l_seq)) {
                char printable[MIN_(1000,l_seq)+1]; // +1 for \0
                WARN_ONCE ("Warning: when converting SAM sequence data to BAM (QNAME=%.*s): invalid character encountered, it will be converted as 'N': '%c' (ASCII %u) (this warning will appear only once). SEQ(first 1000 bases)=\"%s\" seq_i=%u", 
                           vb->last_txt_len(SAM_QNAME), last_txt (vb, SAM_QNAME), base[b], base[b], str_to_printable (recon, MIN_(1000,l_seq), printable), i*2+b);
                base[b] = 0x0f;
            }

        *seq_after = (base[0] << 4) | (base[1] & 0x0f);
    }

    if (l_seq&1) seq_after[-1] &= 0xf0; // if sequence has an odd number of bases, make the final half-byte 0
    
    vb->txt_data.len32 = vb->txt_data.len32 - l_seq + (l_seq+1)/2;

    return 0;
}

// SAM-to-FASTQ translator: reverse-complement the sequence if needed, and drop if "*"
TRANSLATOR_FUNC (sam_piz_sam2fastq_SEQ)
{  
    uint16_t sam_flag = (uint16_t)vb->last_int(SAM_FLAG);
    
    // case: SEQ is "*" - don't show this fastq record
    if (recon_len==1 && *recon == '*') 
        vb->drop_curr_line = "no_seq";

    // case: this sequence is reverse complemented - reverse-complement it
    else if (sam_flag & SAM_FLAG_REV_COMP) 
        str_revcomp (STRa(recon));

    return 0;
}
