// ------------------------------------------------------------------
//   sam_seq.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

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

// called when SEQ in a prim or supp/sec line is segged against ref against prim: 
// Sets vb->mismatch_bases_by_SEQ. ZIP: also sets vb->md_verified and. PIZ: also returns sqbitmap of this line in line_sqbitmap 
static bool sam_analyze_copied_SEQ (VBlockSAMP vb, STRp(seq), const SamPosType pos, bool is_revcomp,
                                    uint32_t ref_consumed, uint32_t ref_and_seq_consumed,
                                    unsigned recursion_level, uint32_t level_0_seq_len, uint32_t bit_i/*0 for recursion level 0*/,
                                    BufferP line_sqbitmap)
{
    RefLock lock = REFLOCK_NONE;
    uint32_t save_ref_and_seq_consumed = ref_and_seq_consumed;

    ASSERT (recursion_level < 500, "%s: excess recursion recursion_level=%u seq_len=%u", // a large number of recursion calls can happen if CIGAR=9M10910N86M3274690N30M1S as observed with the STAR aligner https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            LN_NAME, recursion_level, seq_len);

    BitArrayP bitmap = (BitArrayP)line_sqbitmap;

    if (vb->cigar_missing) goto fail;

    if (!recursion_level) {
        ASSERT (!line_sqbitmap->len32, "%s: line_sqbitmap is in use", LN_NAME);
        buf_alloc_bitarr (VB, line_sqbitmap, ref_and_seq_consumed, line_sqbitmap->name ? line_sqbitmap->name : "line_sqbitmap");        
        bit_array_clear_excess_bits_in_top_word (bitmap);
        bit_array_set_region (bitmap, 0, ref_and_seq_consumed); // we initialize all the bits to "set", and clear as needed.
    }

    ConstRangeP range = (command == ZIP) ? ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, seq, 
                                                                     flag.reference == REF_EXTERNAL ? NULL : &lock)
                                         : ref_piz_get_range (VB, gref, false);
    if (!range) goto fail; // can only happen in ZIP

    uint32_t pos_index = pos - range->first_pos;
    uint32_t next_ref  = pos_index;

    uint32_t i=0, n=0;
    BamCigarOpType op = BC_NONE;

    uint32_t ref_len_this_level = (flag.reference == REF_INTERNAL ? MIN_(ref_consumed, range->last_pos - pos + 1)
                                                                  : ref_consumed); // possibly going around the end of the chromosome in case of a circular chromosome                                   

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
        
    while (i < seq_len || next_ref < pos_index + ref_len_this_level) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_len_this_level, "i or next_ref are out of range");

        if (vb->binary_cigar.next < vb->binary_cigar.len) {
            BamCigarOp *next_op = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next++);
            op = next_op->op;
            n  = next_op->n;
        }

        switch (op) {

            case BC_M: case BC_E: case BC_X: // alignment match or sequence match or mismatch

                ASSERT (n > 0 && n <= (seq_len - i), 
                        "%s: CIGAR implies seq_len longer than actual seq_len=%u (line=%s n=%u recursion_level=%u level_0_seq_len=%u). CIGAR=\"%s\"", 
                        LN_NAME, seq_len, line_name(VB).s, n, recursion_level, level_0_seq_len, vb->last_cigar);

                uint32_t start_i = i;
                for (; n && next_ref < pos_index + ref_len_this_level; n--, next_ref++, i++) {

                    // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                    uint32_t actual_next_ref = (next_ref >= range_len) ? next_ref - range_len : next_ref; 

                    // case: we don't have the value for the REF site, therefore we can't seg MD against REF (its not worth it to add the REF site just for MD)
                    if (command == ZIP && (flag.reference & REF_STORED) && !ref_is_nucleotide_set (range, actual_next_ref)) 
                        goto fail;

                    // case our seq is identical to the reference at this site
                    else if (IS_NUCLEOTIDE (seq[i]) && seq[i] == ref_base_by_idx (range, actual_next_ref)) 
                        bit_i++; // bit remains set. 
 
                    // case: ref is set to a different value - update the bitmap
                    else {
                        bit_array_clear (bitmap, bit_i++);
                        vb->mismatch_bases_by_SEQ++;
                    } 
                }
                ref_and_seq_consumed -= (i - start_i); // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap

                break; // end if 'M', '=', 'X'

            // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
            case BC_I: case BC_S: 
                ASSERT (n > 0 && n <= (seq_len - i), "%s: CIGAR %s implies seq_len longer than actual seq_len=%u", LN_NAME, vb->last_cigar, seq_len);
                i += n;
                n = 0;
                break;

            // for Deletion or Skipping - we move the next_ref ahead
            case BC_D: case BC_N: {
                unsigned ref_consumed_skip = (flag.reference == REF_INTERNAL ? MIN_(n, range_len - next_ref)
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
                        "(cigar=%s pos=%d recursion_level=%u level_0_seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u ref_len_this_level=%u op.n=%u range=[%.*s %"PRId64"-%"PRId64"])",
                        LN_NAME, pos_index + ref_len_this_level - next_ref, seq_len-i, vb->last_cigar, pos, recursion_level, level_0_seq_len,
                        ref_consumed, next_ref, pos_index, ref_len_this_level, n, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            default:
                ABORT ("%s: Invalid CIGAR op=%u", LN_NAME, op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_len_this_level && n) break;
    }

    ref_unlock (gref, lock); // does nothing if REFLOCK_NONE

    uint32_t this_seq_last_pos = pos + (next_ref - pos_index) - 1;

    // in REF_INTERNAL, the sequence can flow over to the next range as each range is 1M bases. this cannot happen
    // in REF_EXTERNAL as each range is the entire contig
    ASSERT (flag.reference == REF_INTERNAL || i == seq_len, "%s: expecting i(%u) == seq_len(%u) pos=%d range=[%.*s %"PRId64"-%"PRId64"] (cigar=%s recursion_level=%u level_0_seq_len=%u)", 
            LN_NAME, i, seq_len, pos, STRf(range->chrom_name), range->first_pos, range->last_pos, vb->last_cigar, recursion_level, level_0_seq_len);

    // case: we have reached the end of the current reference range, but we still have sequence left - 
    // call recursively with remaining sequence and next reference range 
    if (i < seq_len) {

        ASSERT (this_seq_last_pos <= MAX_POS_SAM, "%s: POS=%d and the consumed reference implied by CIGAR=\"%s\", exceeding MAX_POS=%"PRId64
                " (next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                LN_NAME, pos, vb->last_cigar, MAX_POS_SAM, next_ref, pos_index, ref_len_this_level, n, 
                STRf(range->chrom_name), range->first_pos, range->last_pos);

        // if current op is not exhausted, next recursion will be starting from current op
        uint32_t save_op_n=0;
        if (n) {
            vb->binary_cigar.next--; 
            save_op_n = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next)->n;
            B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next)->n = n;
        }

        bool success = sam_analyze_copied_SEQ (vb, seq + i, seq_len - i, range->last_pos + 1, 
                                            is_revcomp, ref_consumed - ref_len_this_level, ref_and_seq_consumed,
                                            recursion_level + 1, level_0_seq_len, bit_i, line_sqbitmap); 

        if (n) B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next)->n = save_op_n; // restore;
    
        if (!success) goto fail;
    }

    // ZIP: final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (!recursion_level && command == ZIP) {
        sam_MD_Z_verify_due_to_seq (vb, STRa(seq), pos, bitmap, 0); // sets vb->md_verified
        buf_free (*line_sqbitmap);
    }
    
    return true;

fail:
    // case PIZ: this function is called if required by reconstruction, it's a bug if it doesn't succeed.
    ASSPIZ (command != PIZ, "failed to reconstruct MD because could not get bitmap for chrom=\"%.*s\" pos=%u ref_and_seq_consumed=%u", STRf(vb->chrom_name), pos, save_ref_and_seq_consumed);

    // case ZIP: if we cannot verify against the reference, the MD:Z is not verified, and we don't use the SPECIAL for segging
    ref_unlock (gref, lock); // does nothing if REFLOCK_NONE
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
    Context *seqsa_ctx  = CTX(SAM_SEQSA); 

    bitmap_ctx->ltype = LT_BITMAP;
    strand_ctx->ltype = LT_BITMAP;
    gpos_ctx->ltype   = LT_UINT32;

    // depn sequences - may occur in MAIN and DEPN components
    if (sam_is_main_vb || sam_is_depn_vb) {
        seqsa_ctx->ltype = LT_BITMAP; // bitwise-xor of prim vs depn sequence

        // initialize so we can use bit_array_realloc
        BufferP xor_buf = &seqsa_ctx->local;
        xor_buf->vb     = VB; 
        xor_buf->name   = "contexts->local";
        buf_add_to_buffer_list_(vb, xor_buf, "local"); 
    }

    // MAIN: we may seg depn lines against in-VB prim lines
    if (sam_is_main_vb)
        bitmap_ctx->flags.store_per_line = true; // 14.0.0

    // initial allocations, these might grow during segging if needed
    int factor = segconf.sam_is_unmapped ? 1 : 32; // if mapped, we allocate assuming 1 out of 32 lines is unmapped
    buf_alloc (vb, &bitmap_ctx->local, 1, vb->txt_data.len / (4 * factor), uint8_t, 0, "contexts->local"); 
    buf_alloc (vb, &strand_ctx->local, 1, roundup_bits2bytes64 (vb->lines.len / factor), uint8_t, 0, "contexts->local"); 
    buf_alloc (vb, &gpos_ctx->local, 1, vb->lines.len / factor, uint32_t, CTX_GROWTH, "contexts->local"); 

    if (!segconf.sam_is_unmapped || flag.aligner_available) {
        buf_alloc (vb, &nonref_ctx->local, 0, vb->txt_data.len / 64, char, 0, "local");

        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, 1, vb->txt_data.len / 128, char, 0, "contexts->local"); 
    }
    else // we store seq vertabim - no reference and no CIGAR
        buf_alloc (vb, &nonref_ctx->local, 0, vb->txt_data.len / 3, char, 0, "local");
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

// Creates a bitmap from seq data - exactly one bit per base that consumes reference (i.e. not bases with CIGAR I/S)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF (unless we can use an aligner)
// - Edge case: no CIGAR (it is "*") - we just store the sequence in SAM_NONREF
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
static MappingType sam_seg_SEQ_vs_ref (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(seq), const SamPosType pos, bool is_revcomp,
                                       uint32_t ref_consumed, uint32_t ref_and_seq_consumed,
                                       unsigned recursion_level, uint32_t level_0_seq_len)
{
    Context *bitmap_ctx = CTX(SAM_SQBITMAP);
    Context *nonref_ctx = CTX(SAM_NONREF);
    Context *seqmis_ctx = CTX(SAM_SEQMIS_A); // 4 contexts
    uint32_t save_ref_and_seq_consumed = ref_and_seq_consumed;

    ASSERT (recursion_level < 500, "excess recursion recursion_level=%u seq_len=%u", // a large number of recursion calls can happen if CIGAR=9M10910N86M3274690N30M1S as observed with the STAR aligner https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            recursion_level, seq_len);

    BitArrayP bitmap = (BitArrayP)&bitmap_ctx->local;
    uint32_t bitmap_start = bitmap_ctx->next_local;        

    if (!recursion_level) {

        ASSERTW (seq_len < 100000 || segconf.running || segconf.is_long_reads, 
                 "%s: Warning: sam_seg_SEQ: seq_len=%u is suspiciously high and might indicate a bug", LN_NAME, seq_len);

        int64_t missing_bits = (int64_t)ref_and_seq_consumed - ((int64_t)bitmap_ctx->local.nbits - (int64_t)bitmap_ctx->next_local);
        if (missing_bits > 0) {
            buf_alloc_do (VB, &bitmap_ctx->local, roundup_bits2bytes64 (bitmap_ctx->local.nbits + seq_len), CTX_GROWTH, __FUNCLINE, NULL); 
            buf_extend_bits (&bitmap_ctx->local, missing_bits);
        }
        
        bit_array_set_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // we initialize all the bits to "set", and clear as needed.
        
        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, ref_and_seq_consumed, 0, char, CTX_GROWTH, "contexts->local"); 

        buf_alloc (vb, &nonref_ctx->local, seq_len + 3, 0, uint8_t, CTX_GROWTH, "contexts->local"); 
    
        bitmap_ctx->local_num_words++;
    }

    RefLock lock = REFLOCK_NONE;
    Range *range = vb->cigar_missing ? NULL 
                                     : ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, seq, 
                                                                 flag.reference == REF_EXTERNAL ? NULL : &lock);

    // Cases where we don't consider the refernce and just copy the seq as-is
    // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
    // 2. (loaded:) case contig doesn't exist in the reference, or POS is out of range of contig (observed in the wild with chrM)    
    if (!range) {
        vb->md_verified = false; // we can't seg MD:Z with special
        vb->mismatch_bases_by_SEQ = -1; // we can't seg NM:i with special

        // case: we're in the initial call - just return and seg as a verbatim copy
        if (!recursion_level) return MAPPING_NO_MAPPING;

        // case: recursive call - we already segged partially against the reference, so we will continue that    
        else {
            // entire remaining seq is added to NONREF
            buf_add_more (VB, &nonref_ctx->local, STRa(seq), "contexts->local");
            bit_array_clear_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // note: vb->ref_and_seq_consumed==0 if cigar="*"
            bitmap_ctx->next_local += ref_and_seq_consumed;

            goto done;
        }    
    }

    uint32_t pos_index     = pos - range->first_pos;
    uint32_t next_ref      = pos_index;

    uint32_t i=0, n=0;
    BamCigarOpType op = BC_NONE;

    uint32_t ref_len_this_level = (flag.reference == REF_INTERNAL ? MIN_(ref_consumed, range->last_pos - pos + 1)
                                                                  : ref_consumed); // possibly going around the end of the chromosome in case of a circular chromosome                                   

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    
    if (flag.reference == REF_EXT_STORE) {
        uint32_t overflow = (pos_index + ref_len_this_level > range_len) ? (pos_index + ref_len_this_level - range_len) : 0;
        bit_array_set_region (&range->is_set, pos_index, ref_len_this_level - overflow); // we will need this ref to reconstruct
        if (overflow) // can only happen with external reference, expected only with circular chromosomes
            bit_array_set_region (&range->is_set, 0, overflow); // round robin to beginning
    }

    while (i < seq_len || next_ref < pos_index + ref_len_this_level) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_len_this_level, "i or next_ref are out of range");

        if (vb->binary_cigar.next < vb->binary_cigar.len) {
            BamCigarOp *next_op = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next++);
            op = next_op->op;
            n  = next_op->n;
        }

        switch (op) {

            case BC_M: case BC_E: case BC_X: // alignment match or sequence match or mismatch

                ASSERT (n > 0 && n <= (seq_len - i), 
                        "%s: CIGAR implies seq_len longer than actual seq_len=%u (n=%u recursion_level=%u level_0_seq_len=%u). CIGAR=\"%s\"", 
                        LN_NAME, seq_len, n, recursion_level, level_0_seq_len, vb->last_cigar);

                uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
                uint32_t start_i = i;
                while (n && next_ref < pos_index + ref_len_this_level) {

                    // when we have an X we don't enter it into our internal ref, and we wait for a read with a = or M for that site,
                    // as we assume that on average, more reads will have the reference base, leading to better compression
                
                    bool normal_base = IS_NUCLEOTIDE (seq[i]);

                    // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                    uint32_t actual_next_ref = (next_ref >= range_len) ? next_ref - range_len : next_ref; 

                    // case: we have not yet set a value for this site - we set it now. note: in ZIP, is_set means that the site
                    // will be needed for pizzing. With REF_INTERNAL, this is equivalent to saying we have set the ref value for the site
                    if (flag.reference == REF_INTERNAL && !ref_is_nucleotide_set (range, actual_next_ref)) { 

                        // note: in case this is a non-normal base (eg N), set the reference to an arbitrarily to 'A' as we 
                        // we will store this non-normal base in seqmis_ctx multiplexed by the reference base (i.e. in seqmis_ctx['A']).
                        ref_set_nucleotide (range, actual_next_ref, normal_base ? seq[i] : 'A');
                        
                        bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct

                        if (normal_base) 
                            bit_i++; 
                        else
                            goto mismatch;
                    }

                    // case our seq is identical to the reference at this site
                    else if (normal_base && seq[i] == ref_base_by_idx (range, actual_next_ref)) 
                        bit_i++; // bit remains set. 
 
                    // case: ref is set to a different value - we store our value in nonref_ctx
                    else mismatch: {
                        uint8_t ref_base_2bit = bit_array_get2 (&range->ref, actual_next_ref * 2);
                        BNXTc (seqmis_ctx[ref_base_2bit].local) = seq[i];
                        bit_array_clear (bitmap, bit_i++);
                        vb->mismatch_bases_by_SEQ++;
                    } 

                    n--;
                    next_ref++;
                    i++;
                }
                ref_and_seq_consumed -= (i - start_i); // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap
                bitmap_ctx->next_local = bit_i;

                break; // end if 'M', '=', 'X'

            // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
            case BC_I: case BC_S: 
                ASSSEG (n > 0 && n <= (seq_len - i), seq,
                        "CIGAR %s implies seq_len longer than actual seq_len=%u", vb->last_cigar, seq_len);

                buf_add (&nonref_ctx->local, &seq[i], n);
                i += n;
                n = 0;
                break;

            // for Deletion or Skipping - we move the next_ref ahead
            case BC_D: case BC_N: {
                unsigned ref_consumed_skip = (flag.reference == REF_INTERNAL ? MIN_(n, range_len - next_ref)
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
                ASSSEG (false, vb->last_cigar, "End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                        "(cigar=%s pos=%d recursion_level=%u level_0_seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u ref_len_this_level=%u op.n=%u range=[%.*s %"PRId64"-%"PRId64"])",
                        pos_index + ref_len_this_level - next_ref, seq_len-i, vb->last_cigar, pos, recursion_level, level_0_seq_len,
                        ref_consumed, next_ref, pos_index, ref_len_this_level, n, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            default:
                ASSSEG (false, vb->last_cigar, "Invalid CIGAR op=%u", op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_len_this_level && n) break;
    }

    ref_unlock (gref, lock); // does nothing if REFLOCK_NONE

    uint32_t this_seq_last_pos = pos + (next_ref - pos_index) - 1;

    // in REF_INTERNAL, the sequence can flow over to the next range as each range is 1M bases. this cannot happen
    // in REF_EXTERNAL as each range is the entire contig
    ASSERT (flag.reference == REF_INTERNAL || i == seq_len, "expecting i(%u) == seq_len(%u) pos=%d range=[%.*s %"PRId64"-%"PRId64"] (cigar=%s recursion_level=%u level_0_seq_len=%u)", 
            i, seq_len, pos, STRf(range->chrom_name), range->first_pos, range->last_pos, vb->last_cigar, recursion_level, level_0_seq_len);

    // case: we have reached the end of the current reference range, but we still have sequence left - 
    // call recursively with remaining sequence and next reference range 
    if (i < seq_len) {

        ASSSEG (this_seq_last_pos <= MAX_POS_SAM, vb->last_cigar, "POS=%d and the consumed reference implied by CIGAR=\"%s\", exceeding MAX_POS=%"PRId64
                " (next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                pos, vb->last_cigar, MAX_POS_SAM, next_ref, pos_index, ref_len_this_level, n, 
                STRf(range->chrom_name), range->first_pos, range->last_pos);

        // if current op is not exhausted, next recursion will be starting from current op
        uint32_t save_op_n=0;
        if (n) {
            vb->binary_cigar.next--; 
            save_op_n = B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next)->n;
            B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next)->n = n;
        }

        sam_seg_SEQ_vs_ref (vb, dl, seq + i, seq_len - i, range->last_pos + 1, 
                            is_revcomp, ref_consumed - ref_len_this_level, ref_and_seq_consumed,
                            recursion_level + 1, level_0_seq_len); // ignore return value

        if (n) B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next)->n = save_op_n; // restore;
    }
    
    // final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (!recursion_level)
        sam_MD_Z_verify_due_to_seq (vb, STRa(seq), pos, bitmap, bitmap_start);
        //xxx BitArrayP M_is_ref = (BitArrayP)&vb->md_M_is_ref;

        // bool bitmap_matches_MD = vb->md_verified && !bit_array_hamming_distance (M_is_ref, 0, bitmap, bitmap_start, M_is_ref->nbits);

        // if (flag.show_wrong_md && vb->md_verified && !bitmap_matches_MD) {
        //     iprintf ("%s RNAME=%.*s POS=%d CIGAR=%s MD=%.*s SEQ=%.*s\n", 
        //             LN_NAME, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(VB, OPTION_MD_Z), STRf(seq));
        //     bit_array_print_substr ("SEQ match to ref", bitmap, bitmap_start, M_is_ref->nbits, info_stream);
        //     bit_array_print_substr ("MD implied match", M_is_ref, 0, M_is_ref->nbits, info_stream); 
        // }

        // vb->md_verified = bitmap_matches_MD;

done:
    if (recursion_level) return 0; // recursive calling ignores return value

    // note: change of logic in v14 vs v13: in v13, we use to align every recursion level, it seems that this was a bug
    sam_seg_SEQ_pad_nonref (vb);

    if (vb->mismatch_bases_by_SEQ != 0) 
        return MAPPING_ALIGNED;
    
    else {
        bitmap_ctx->local.nbits -= save_ref_and_seq_consumed; // we don't use the bitmap if there are no mismatches
        bitmap_ctx->next_local  -= save_ref_and_seq_consumed;
        return MAPPING_PERFECT;
    } 
}

// add the diff between our seq and prim to seq to SAM_SEQSA.local. 
static bool sam_seg_depn_SEQ (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(seq), bool is_revcomp, // textual (SAM) format
                              BufferP prim_acgt, uint64_t prim_seq_index, bool prim_revcomp,
                              bool soft_fail, bool *force_no_analyze_depn_SEQ/*in/out*/)
{
    // extend SEQSA
    BitArray *xor_bits = buf_get_bitarray (&CTX(SAM_SEQSA)->local);
    uint64_t start_xor = xor_bits->nbits;
    bit_array_realloc (xor_bits, start_xor + seq_len * 2, 0, false);

    // SEQSA = DEPN (converted to ACGT format, possibly with revcomp)
    bool xstrand = (is_revcomp != prim_revcomp); // primary and dependent are on opposite strands
    if (!sam_sa_native_to_acgt (vb, xor_bits, start_xor, STRa(seq), false, xstrand, soft_fail)) // If PRIM in SA_Grp: we already verified that it is ACGT in sam_sa_seg_depn_find_sagroup
        return false; // note: the function ABORTs if not soft_fail

    // SEQSA ^= PRIM 
    BitArray *sa_seq = buf_get_bitarray (prim_acgt);
        
    uint64_t prim_start_base = prim_seq_index + vb->hard_clip[xstrand];

    bit_array_xor_with (xor_bits, start_xor, sa_seq, 2*prim_start_base, 2*seq_len);

    // set vb->md_verified and vb->mismatch_bases_by_SEQ needed by MD:Z and NM:i
    if (segconf.has_MD_or_NM && !vb->cigar_missing && 
        !sam_analyze_copied_SEQ (vb, STRa(seq), dl->POS, dl->FLAG.bits.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed, 0, seq_len, 0, &vb->codec_bufs[0]))
        *force_no_analyze_depn_SEQ = true;

    return true;
}

void sam_seg_SEQ (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(seq), unsigned add_bytes)
{
    START_TIMER;

    bool force_sequence = false; // true if reconstructor would normally diff, but in this case it should reconstruct from SQBITMAP.local
    bool aligner_used = false;
    bool perfect = false;

    ZipDataLineSAM *prim_dl;
    bool unmapped = dl->FLAG.bits.unmapped || vb->cigar_missing || !dl->POS || str_issame_(STRa(vb->chrom_name), "*", 1);
    bool force_no_analyze_depn_SEQ = false;

    // case segconf: we created the contexts for segconf_set_vb_size accounting. that's enough - actually segging will mark is_set and break sam_seg_MD_Z_analyze.
    if (segconf.running) {
        segconf.longest_seq_len = MAX_(segconf.longest_seq_len, seq_len);

        if (has_MD || has_NM) segconf.has_MD_or_NM = true;

        if (!vb->line_i) segconf_mark_as_used (VB, 5, SAM_SQBITMAP, SAM_NONREF, SAM_NONREF_X, SAM_GPOS, SAM_STRAND);

        return; 
    }

    // case: unmapped line and we have refhash: align to reference
    if (unmapped && flag.aligner_available) {
    
        switch (aligner_seg_seq (VB, CTX(SAM_SQBITMAP), STRa(seq), true, false, NO_GPOS, false)) {
            case MAPPING_NO_MAPPING : goto add_seq_verbatim;
            case MAPPING_PERFECT    : perfect      = true; // fallthrough
            case MAPPING_ALIGNED    : aligner_used = true; break;
            default                 : ABORT0 ("bad value");
        }

        sam_seg_SEQ_pad_nonref (vb);
        vb->md_verified = false;    
    }

    // case: unmapped line, no refhash: just store the sequence in nonref without an indication in the bitmap
    else if (unmapped && !flag.aligner_available) {
        add_seq_verbatim:        
        buf_add_more (VB, &CTX(SAM_NONREF)->local, STRa(seq), "contexts->local");
        sam_seg_SEQ_pad_nonref (vb);
        vb->md_verified = false;    
    }

    // case: no SEQ. we already add '-' to CIGAR - no data added here
    else if (seq[0] == '*') {
        vb->md_verified = false;    
        vb->seq_missing = true;
    }

    // in case of DEPN line with SEQ confirmed by sam_sa_seg_depn_find_sagroup to be a subset of the SA Group, 
    // no need to Seg - sam_piz_special_SEQ will copy from SA Group
    else if (sam_is_depn_vb && vb->sa_grp) 
        sam_seg_depn_SEQ (vb, dl, STRa(seq), dl->FLAG.bits.rev_comp, &z_file->sa_seq, vb->sa_grp->seq, vb->sa_grp->revcomp, false, &force_no_analyze_depn_SEQ);

    // case: DEPN line vs same-VB prim
    else if (({ prim_dl = DATA_LINE (vb->prim_line_i) ; zip_has_prim; })) {
        
        if (prim_dl->QUAL.len != prim_dl->SEQ.len || !dl->SEQ.len) {
            force_sequence = true;
            goto vs_ref;
        }

        ASSERTNOTINUSE (vb->scratch);
        BitArrayP prim = buf_alloc_bitarr (vb, &vb->scratch, prim_dl->SEQ.len * 2, "scratch");
        
        // convert prim to acgt and compare to depn. fails if either is not all A,C,G,T
        bool success = sam_sa_native_to_acgt (vb, prim, 0, STRtxtw(prim_dl->SEQ), IS_BAM_ZIP, false, true) &&
                       sam_seg_depn_SEQ (vb, dl, STRa(seq), dl->FLAG.bits.rev_comp, &vb->scratch, 0, prim_dl->FLAG.bits.rev_comp, true, &force_no_analyze_depn_SEQ);

        buf_free (vb->scratch);

        if (!success) {
            force_sequence = true;
            goto vs_ref;
        }
    }

    else vs_ref:         
        switch (sam_seg_SEQ_vs_ref (vb, dl, STRa(seq), dl->POS, dl->FLAG.bits.rev_comp, vb->ref_consumed, vb->ref_and_seq_consumed, 0, seq_len)) {
            case MAPPING_NO_MAPPING : goto add_seq_verbatim;
            case MAPPING_PERFECT    : perfect = true; // fallthrough
            default                 : break;
        }

    // case: PRIM - SQBITMAP context is reconstructed when loading SA Groups (preprocessing) - invoking SPECIAL_SEQ and reconstructing the sequence
    //       and SQBITMAP is reconstructed again during recon - but this SPECIAL_SEQ copies from SA Groups
    // case: MAIN/DEPN - SPECIAL_SEQ is invoked by SQBITMAP - it either reconstructs sequence or diffs vs prim.
    //       (prim SEQ may be taken from SA Groups or from earlier in the current VB (the latter only in MAIN VBs)
    seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SEQ, '0'+force_sequence, '0'+aligner_used, '0'+perfect, '0'+force_no_analyze_depn_SEQ }, 6, 
                  SAM_SQBITMAP, add_bytes); 

    COPY_TIMER (sam_seg_SEQ);
}

// converts native SAM/BAM format to 2bit ACGT - if soft_fail, returns false if any base is not A,C,G or T 
bool sam_sa_native_to_acgt (VBlockSAMP vb, BitArray *packed, uint64_t next_bit, STRp(seq), bool bam_format, bool revcomp, bool soft_fail)
{
    if (bam_format) {
        if (!revcomp)
            for (uint32_t i=0; i < seq_len; i++, next_bit += 2) {
                uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
                switch (b) {
                    case 0b0001 : bit_array_assign2 (packed, next_bit, 0); break;
                    case 0b0010 : bit_array_assign2 (packed, next_bit, 1); break;
                    case 0b0100 : bit_array_assign2 (packed, next_bit, 2); break;
                    case 0b1000 : bit_array_assign2 (packed, next_bit, 3); break;
                    default     : if (soft_fail) return false;
                                  ABORT ("%s: Unexpected base: '%c' i=%u seq_len=%u", LN_NAME, bam_base_codes[b], i, seq_len);
                }
            }    
        else
            for (int32_t i=seq_len-1; i >= 0; i--, next_bit += 2) {
                uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
                switch (b) {
                    case 0b0001 : bit_array_assign2 (packed, next_bit, 3); break;
                    case 0b0010 : bit_array_assign2 (packed, next_bit, 2); break;
                    case 0b0100 : bit_array_assign2 (packed, next_bit, 1); break;
                    case 0b1000 : bit_array_assign2 (packed, next_bit, 0); break;
                    default     : if (soft_fail) return false;
                                  ABORT ("%s: Unexpected base: '%c' i=%u seq_len=%u", LN_NAME, bam_base_codes[b], i, seq_len);
                }
            }    
    }

    else { // SAM
        if (!revcomp)
            for (uint32_t i=0; i < seq_len; i++, next_bit += 2) 
                switch (seq[i]) {
                    case 'A' : bit_array_assign2 (packed, next_bit, 0); break;
                    case 'C' : bit_array_assign2 (packed, next_bit, 1); break;
                    case 'G' : bit_array_assign2 (packed, next_bit, 2); break;
                    case 'T' : bit_array_assign2 (packed, next_bit, 3); break;
                    default  : if (soft_fail) return false;
                               ABORT ("%s: Unexpected base: '%c'(ASCII %u) i=%u seq_len=%u", LN_NAME, seq[i], (uint8_t)seq[i], i, seq_len);
                }
        else
            for (int32_t i=seq_len-1; i >= 0; i--, next_bit += 2) 
                switch (seq[i]) {
                    case 'A' : bit_array_assign2 (packed, next_bit, 3); break;
                    case 'C' : bit_array_assign2 (packed, next_bit, 2); break;
                    case 'G' : bit_array_assign2 (packed, next_bit, 1); break;
                    case 'T' : bit_array_assign2 (packed, next_bit, 0); break;
                    default  : if (soft_fail) return false;
                               ABORT ("%s: Unexpected base: '%c'(ASCII %u) i=%u seq_len=%u", LN_NAME, seq[i], (uint8_t)seq[i], i, seq_len);
                }
    }

    bit_array_clear_excess_bits_in_top_word (packed);

    return true;
}

// pack seq (into 2bit ACGT format) of each PRIM line separately 
void sam_zip_prim_ingest_vb_pack_seq (VBlockSAMP vb, SAGroupType *vb_grps, uint32_t vb_grps_len,
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
        
        SAGroupType *vb_grp = &vb_grps[vb_grp_i];

        uint64_t next_bit = packed_seq_buf->nbits;
        buf_extend_bits (packed_seq_buf, vb_grp->seq_len * 2); // extend now 
        BitArray *sa_seq = buf_get_bitarray (packed_seq_buf);

        sam_sa_native_to_acgt (vb, sa_seq, next_bit, Bc (vb->txt_data, vb_grp->seq), vb_grp->seq_len, is_bam_format, false, false);
        vb_grp->seq = next_bit / 2; // update from an index into txt_data to an index (bases not bits) into sa_seq
    }
}

// used by codec_longr_compress
COMPRESSOR_CALLBACK (sam_zip_seq) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    *line_data_len = dl->SEQ.len;

    if (VB_DT(DT_SAM)) 
        *line_data = Btxt (dl->SEQ.index);
    else { // BAM
        VB_SAM->textual_seq.len = 0; // reset
        buf_alloc (vb, &VB_SAM->textual_seq, 0, dl->SEQ.len+1 /* +1 for last half-byte */, char, 1.5, "textual_seq");

        bam_seq_to_sam (VB_SAM, (uint8_t *)Btxt(dl->SEQ.index), dl->SEQ.len, false, false, &VB_SAM->textual_seq);
        *line_data = B1STc (VB_SAM->textual_seq);
    }

    if (is_rev) *is_rev = dl->FLAG.bits.rev_comp;
}

//---------
// PIZ
//---------

static inline bool v13_sam_recon_SEQ_consume_query 
    (VBlockSAMP vb, ContextP bitmap_ctx, ConstRangeP range, SamPosType range_len,
     unsigned seq_consumed, unsigned ref_consumed, SamPosType pos, bool consumes_reference)
{
    if (consumes_reference && NEXTLOCALBIT (bitmap_ctx)) /* copy from reference */ {

        if (!vb->drop_curr_line) { // note: if this line is excluded with --regions, then the reference section covering it might not be loaded
            uint32_t idx = ((pos - range->first_pos) + ref_consumed) % range_len; // circle around (this can only happen when compressed with an external reference)

            if (!ref_is_nucleotide_set (range, idx) &&
                    (!flag.regions || regions_is_site_included (VB))) { // if this line is not included, then possibly its reference range is not loaded. we complete consumption (resulting in bad reconstruction) and drop the line in container_reconstruct_do

                ref_print_is_set (range, pos + ref_consumed, stderr);
                ASSPIZ (false, "Error in sam_reconstruct_SEQ: reference is not set: chrom=%u \"%.*s\" pos=%u range=[%"PRId64"-%"PRId64"]"
                        " (cigar=%s seq_start_pos=%u ref_consumed=%u seq_consumed=%u)",
                        range->chrom, range->chrom_name_len, range->chrom_name, pos + ref_consumed, 
                        range->first_pos, range->last_pos, vb->last_cigar, pos, ref_consumed, seq_consumed);
            }

            char ref = ref_base_by_idx (range, idx);
            RECONSTRUCT1 (ref); 
        }
        return true;
    }
    else
        return false;
}

// PIZ: SEQ reconstruction 
void sam_reconstruct_SEQ (VBlockP vb_, Context *bitmap_ctx, STRp(snip), bool reconstruct)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (!bitmap_ctx->is_loaded) return; // if case we need to skip the SEQ field (for the entire VB)

    Context *nonref_ctx   = CTX(SAM_NONREF);
    Context *seqmis_ctx   = CTX(SAM_SEQMIS_A);
    rom nonref            = Bc (nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    rom nonref_start      = nonref;
    const SamPosType pos  = vb->last_int(SAM_POS);
    ConstRangeP range     = NULL;
    unsigned seq_consumed=0, ref_consumed=0;

    bool unmapped = VER(14) ? (last_flags.bits.unmapped || vb->cigar_missing || !pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*'))
                            : (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')); // the criterion seg used in up to v13
                                                  
    bool aligner_used = VER(14) ? (snip && snip[1] == '1') // starting v14, seg tells us explicitly
                                : (unmapped && z_file->z_flags.aligner); 

    bool is_perfect = VER(14) && snip && snip[2] == '1';

    // case: unmapped, segged against reference using our aligner
    if (aligner_used) {
        aligner_reconstruct_seq (VB, bitmap_ctx, vb->seq_len, false, is_perfect, reconstruct);
        nonref_ctx->next_local = ROUNDUP4 (nonref_ctx->next_local);
        return;
    }

    // case: unmapped, and copied verbatim
    if (unmapped) unmapped: {
        if (reconstruct) RECONSTRUCT (nonref, vb->seq_len); 
        nonref_ctx->next_local += ROUNDUP4 (vb->seq_len);
        return;
    }

    // case: missing sequence - sequence is '*' - just reconstruct the '*' (but only if it is the primary SEQ field,
    // which we can tell by not being encountered earlier on the line - bc in v8 (at least) E2 was also stored in the same SAM_SQBITMAP)
    if (vb->seq_missing && !ctx_encountered_in_line (VB, bitmap_ctx->did_i)) {
        if (reconstruct) RECONSTRUCT1 ('*');
        return;
    }

    ASSERT0 (bitmap_ctx->local.len || nonref_ctx->local.len, "No SEQ data, perhaps sections were skipped?");

    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence

    if (vb->ref_consumed) {
        if (!(range = ref_piz_get_range (VB, gref, true))) {
            if (VER(14)) goto unmapped; // starting v14, a missing range is another case of unmapped
            ref_display_all_ranges (gref);
            ASSPIZ0 (false, "range is NULL");
        }
    }

    vb->range = (RangeP)range; // store in vb->range too, for sam_piz_special_MD

    SamPosType range_len = range ? (range->last_pos - range->first_pos + 1) : 0;

    const BamCigarOp *cigar = vb->binary_cigar.len32 ? B1ST(BamCigarOp, vb->binary_cigar) 
                                                     : &(BamCigarOp){ .n = vb->seq_len, .op = BC_S }; // up to v13, alignments with CIGAR=* but with a sequence where segged this way
                                                     
    BamCigarOp op = {};
    bool consumes_query=false, consumes_reference=false;

    while (seq_consumed < vb->seq_len || ref_consumed < vb->ref_consumed) {
        
        if (!op.n) {
            op = *cigar++;
            switch (op.op) {
                case BC_M : case BC_X : case BC_E : consumes_query=true;  consumes_reference=true;  break;
                case BC_I : case BC_S :             consumes_query=true;  consumes_reference=false; break;
                case BC_D : case BC_N :             consumes_query=false; consumes_reference=true;  break;
                case BC_H : case BC_P :             consumes_query=false; consumes_reference=false; break;
                default   : ASSPIZ (false, "Invalid op=%d", op.op);
            }
        }

        if (consumes_query) {
            
            if (!VER(14) && v13_sam_recon_SEQ_consume_query (vb, bitmap_ctx, range, range_len, seq_consumed, ref_consumed, pos, consumes_reference))
                {} // done

            else if (VER(14) && consumes_reference) {
                                
                uint32_t idx = (pos - range->first_pos) + ref_consumed; 
                if (idx >= range_len) idx -= range_len; // circle around (this can only happen when compressed with an external reference)

                if (!ref_is_nucleotide_set (range, idx) &&
                    (!flag.regions || regions_is_site_included (VB))) { 

                    ref_print_is_set (range, pos + ref_consumed, stderr);
                    ASSPIZ (false, "Unexpectedly, reference at this locus has is_set=false: chrom=%u \"%.*s\" pos=%d range=[%"PRId64"-%"PRId64"]"
                            " (cigar=%s seq_start_pos=%d ref_consumed=%u seq_consumed=%u)",
                            range->chrom, STRf(range->chrom_name), pos + ref_consumed, 
                            range->first_pos, range->last_pos, B1STc (vb->textual_cigar), pos, ref_consumed, seq_consumed);
                }

                if (is_perfect || NEXTLOCALBIT (bitmap_ctx)) /* copy from reference */ {
                    char ref = ref_base_by_idx (range, idx);
                    if (reconstruct) RECONSTRUCT1 (ref); 
                }
                else {
                    vb->mismatch_bases_by_SEQ++;

                    if (VER(14)) {
                        uint8_t ref_base_2bit = bit_array_get2 (&range->ref, idx * 2);
                        RECONSTRUCT_NEXT (&seqmis_ctx[ref_base_2bit], 1);
                    }
                    else goto nonref;
                }
            }
            else nonref: {
                if (reconstruct) RECONSTRUCT1 (*nonref);
                nonref++;
            }

            seq_consumed++;
        }

        if (consumes_reference) 
            ref_consumed++;

        op.n--;
    }

    ASSPIZ (seq_consumed == vb->seq_len,      "expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSPIZ (ref_consumed == vb->ref_consumed, "expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    nonref_ctx->next_local += ROUNDUP4 (nonref - nonref_start);
}


// reconstruct from a 2bit array - start_base and seq_len are in bases (not in bits)
static void reconstruct_SEQ_acgt (VBlockSAMP vb, BitArrayP seq_2bits, uint64_t start_base, uint32_t seq_len, bool revcomp)
{
    // Reconstruct as SAM. If BAM, we later translate to BAM. To do: reconstruct directly as BAM if needed, bug 530
    char *next = BAFTc(vb->txt_data);
    
    if (!revcomp)
        for (uint32_t i=0; i < seq_len; i++) {
            uint8_t b = bit_array_get2 (seq_2bits, (start_base + i)*2);
            *next++ = b==0 ? 'A' : b==1 ? 'C' : b==2 ? 'G' : 'T';
        }
    
    else
        for (int32_t i=seq_len-1; i >= 0; i--) {
            uint8_t b = bit_array_get2 (seq_2bits, (start_base + i)*2);
            *next++ = b==3 ? 'A' : b==2 ? 'C' : b==1 ? 'G' : 'T';
        }

    vb->txt_data.len32 += seq_len;
}

// reconstruct as diff vs primary
// note: SEQ can only be reconstructed once, because the reconstruction is destructive to SAM_SEQSA.local (the xor bitarray)
static void reconstruct_SEQ_depn (VBlockSAMP vb, ContextP ctx, BitArrayP prim, uint64_t prim_start_base, uint32_t prim_seq_len, bool xstrand, bool reconstruct, bool force_no_analyze_depn_SEQ)
{
    ContextP xor_bits_ctx = CTX(SAM_SEQSA); // initially, contains XOR diff

    uint32_t depn_seq_len = prim_seq_len - vb->hard_clip[0] - vb->hard_clip[1]; // depn sequence is a sub-sequence of the prim sequence

    if (reconstruct && ctx->is_loaded) {
        BitArray *xor_bits = buf_get_bitarray (&xor_bits_ctx->local); // XOR of DEPN and PRIM 

        // xor_bits ^= prim  --> regenerates depn (still revcomp'ed if xstrand)
        bit_array_xor_with (xor_bits, 2 * xor_bits_ctx->next_local, prim,   
                            2 * (prim_start_base + (xstrand ? (prim_seq_len - depn_seq_len - vb->hard_clip[0]) 
                                                            : vb->hard_clip[0])),
                            2 * depn_seq_len);
        
        // Reconstruct as SAM. If BAM, we later translate to BAM. To do: reconstruct directly as BAM if needed, bug 530
        rom seq = BAFTtxt;
        reconstruct_SEQ_acgt (vb, xor_bits, xor_bits_ctx->next_local, depn_seq_len, xstrand);

        // in needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i
        // NOTE: if SEQ is not loaded or !reconstruct - it is expected that NM and MD also don't reconstruct
        if (segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing && !force_no_analyze_depn_SEQ) 
            sam_analyze_copied_SEQ (vb, seq, depn_seq_len, CTX(SAM_POS)->last_value.i, last_flags.bits.rev_comp,
                                    vb->ref_consumed, vb->ref_and_seq_consumed, 0, depn_seq_len, 0, &ctx->line_sqbitmap);        
    }

    xor_bits_ctx->next_local += depn_seq_len; // in bases
}

// PIZ of a SEQ in a MAIN and DEPN components - this is an all-the-same context - SEQ all lines in a DEPN component
// come here. We multiplex between SAGroup-based reconstruction vs normal SEQ reconstruction
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_SEQ)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    // case: reconstructor would normally diff, but in this case it should reconstruct from SQBITMAP.local    
    if (snip[0] == '1') 
        goto force_sequence;

    // case: PRIM component - copy from SA Group 
    else if (sam_is_prim_vb && !vb->preprocessing) {
        sam_piz_set_sa_grp (vb);

        rom seq = BAFTtxt;
        reconstruct_SEQ_acgt (vb, (BitArrayP)&z_file->sa_seq, vb->sa_grp->seq, vb->sa_grp->seq_len, false);

        // if needed, get vb->mismatch_bases_by_SEQ and line_sqbitmap for use in reconstructing MD:Z and NM:i
        if (segconf.has_MD_or_NM && !vb->cigar_missing && !vb->seq_missing) 
            sam_analyze_copied_SEQ (vb, seq, vb->sa_grp->seq_len, CTX(SAM_POS)->last_value.i, last_flags.bits.rev_comp,
                                    vb->ref_consumed, vb->ref_and_seq_consumed, 0, vb->sa_grp->seq_len, 0, &ctx->line_sqbitmap);        
    }

    // case: DEPN component - line with SA Group - diff vs prim sequence
    else if (sam_is_depn_vb && SAM_PIZ_HAS_SA_GROUP) {
        const SAGroupType *grp = vb->sa_grp;
        ASSPIZ (grp->seq_len > vb->hard_clip[0] + vb->hard_clip[1], "grp_i=%u aln_i=%"PRIu64" grp->seq_len=%u <= hard_clip[LEFT]=%u + hard_clip[RIGHT]=%u", 
                ZGRP_I(grp), ZALN_I(vb->sa_aln), grp->seq_len, vb->hard_clip[0], vb->hard_clip[1]);

        reconstruct_SEQ_depn (vb, ctx, buf_get_bitarray (&z_file->sa_seq), STRa(grp->seq), 
                              (vb->sa_aln->revcomp != grp->revcomp), reconstruct, snip[3]-'0');
    }

    // case: MAIN component - depn line - has buddy against prim line in this VB - diff against prim line
    else if (sam_is_main_vb && (last_flags.bits.secondary || last_flags.bits.supplementary) &&
             piz_has_buddy) {

        SamFlags prim_flags = (SamFlags){ .value = *B(int64_t, CTX(SAM_FLAG)->history, vb->buddy_line_i) };

        HistoryWord word = *B(HistoryWord, ctx->history, vb->buddy_line_i); // SEQ is always stored as LookupTxtData or LookupPerLine
        rom prim_seq = (word.lookup == LookupTxtData) ? Btxt(word.index) : Bc(ctx->per_line, word.index);
        uint32_t prim_seq_len = word.len;

        // get prim in actg format
        ASSERTNOTINUSE (vb->scratch);
        BitArrayP prim = buf_alloc_bitarr (vb, &vb->scratch, prim_seq_len * 2, "scratch");
        sam_sa_native_to_acgt (vb, prim, 0, STRa(prim_seq), flag.out_dt==DT_BAM, false, false);
        reconstruct_SEQ_depn (vb, ctx, prim, 0, prim_seq_len, prim_flags.bits.rev_comp != last_flags.bits.rev_comp, reconstruct, snip[3]-'0');
    
        buf_free (vb->scratch);
    }
    // case: reconstruct sequence directly
    else 
        force_sequence:
        sam_reconstruct_SEQ (VB, ctx, STRa(snip), reconstruct);

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
        vb->txt_data.len--;
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
    
    vb->txt_data.len = vb->txt_data.len - l_seq + (l_seq+1)/2;

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




