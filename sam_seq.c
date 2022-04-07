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

//---------
// ZIP
//---------

static void sam_seg_depn_SEQ (VBlockSAMP vb, STRp(seq), bool is_revcomp); // forward

void sam_seg_SEQ_initialize (VBlockSAMP vb)
{
    // if PRIM - TOPLEVEL reconstructs SAM_SEQSA which copies SEQ from the in-memory SA Group. 
    // SAM_SQBITMAP is still segged, but it is consumed to for loading the in-memory SA Groups.
    if (sam_is_prim_vb)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PRIM_SEQ }, 2, SAM_SEQSA, 0);

    // SEQ: in DEPN component - we seg an all-the-same SPECIAL that either reconstructs SEQ from 
    // SQBITMAP.local or from sam_piz_special_pull_from_SAGROUP. 
    else if (sam_is_depn_vb)
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_DEPN_SEQ }, 2, SAM_SQBITMAP, 0); // seg and not just create_node because we have a local section too - so we need a b250 section not just dict
}

// Creates a bitmap from seq data - exactly one bit per base that consumes reference (i.e. not bases with CIGAR I/S)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF (unless we can use an aligner)
// - Edge case: no CIGAR (it is "*") - we just store the sequence in SAM_NONREF
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
void sam_seg_SEQ (VBlockSAMP vb, DidIType bitmap_did, STRp(seq), const PosType pos, rom cigar, bool is_revcomp,
                  uint32_t ref_consumed, uint32_t ref_and_seq_consumed,
                  unsigned recursion_level, uint32_t level_0_seq_len, rom level_0_cigar, unsigned add_bytes)
{
    START_TIMER;

    Context *bitmap_ctx = CTX(bitmap_did);
    Context *nonref_ctx = bitmap_ctx + 1;

    ASSERT (recursion_level < 500, "excess recursion recursion_level=%u seq_len=%u level_0_cigar=%s", // a large number of recursion calls can happen if CIGAR=9M10910N86M3274690N30M1S as observed with the STAR aligner https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            recursion_level, seq_len, level_0_cigar);

    bitmap_ctx->txt_len += add_bytes; // byte counts for --show-sections

    // in case of DEPN line with SEQ confirmed by sam_sa_seg_depn_find_sagroup to be a subset of the SA Group, 
    // no need to Seg - sam_piz_special_DEPN_SEQ (called by all-the-same SAM_SPECIAL_DEPN_SEQ) will copy from SA Group
    // note: in PRIM, we seg, and it will be conusmed when loading SA Group. When reconstructing Toplevel goes to
    // all-the-same SAM_SEQSA for copying from SA Group.
    if (sam_is_depn_vb && vb->sa_grp) {
        sam_seg_depn_SEQ(vb, STRa(seq), is_revcomp);
        goto done;
    }

    // case: PRIM - TOPLEVEL reconstructs SAM_SEQSA which copies SEQ from the in-memory SA Group. 
    // SAM_SQBITMAP is still segged too, but it is consumed to for loading the in-memory SA Groups.
    if (sam_is_prim_vb) 
        seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PRIM_SEQ }, 2, SAM_SEQSA, 0);

    // for unaligned lines, if we have refhash loaded, use the aligner instead of CIGAR-based segmenting
    if (!pos && segconf.sam_use_aligner) {
        aligner_seg_seq (VB, bitmap_ctx, STRa(seq));
        goto align_nonref_local;
    }

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);
    uint32_t bitmap_start = bitmap_ctx->next_local;        

    if (!recursion_level) {

        ASSERTW (seq_len < 100000 || segconf.running || segconf.is_long_reads, 
                 "Warning: sam_seg_SEQ: seq_len=%u is suspiciously high and might indicate a bug. vb=%u line_i=%d", seq_len, vb->vblock_i, vb->line_i);

        buf_alloc (vb, &bitmap_ctx->local, roundup_bits2words64 (ref_and_seq_consumed), 
                   (segconf.is_long_reads ? vb->lines.len/2 : vb->lines.len) * (MIN_(ref_and_seq_consumed, segconf.sam_seq_len)) / 64,  // this formula was chosen after extensive trial & error - easy to get wrong, esp for long reads
                   uint64_t /* len is in 64b words in a bitmap */, CTX_GROWTH, "contexts->local"); 

        buf_extend_bits (&bitmap_ctx->local, ref_and_seq_consumed); 
        bit_array_set_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // we initiaze all the bits to "set", and clear as needed.

        buf_alloc (vb, &nonref_ctx->local, seq_len + 3, 
                   vb->lines.len * (MIN_(seq_len, segconf.sam_seq_len + 1)) / (segconf.is_long_reads ? 50 : 20), // this one too - it was chosen after extensive trial & error - easy to get wrong, esp for long reads
                   uint8_t, CTX_GROWTH, "contexts->local"); 
    
        if (segconf.running) {
            if (seq_len > segconf.longest_seq_len) segconf.longest_seq_len = seq_len;
            return; // case segconf: we created the contexts for segconf_set_vb_size accounting. that's enough - actually segging will mark is_set and break sam_seg_MD_Z_analyze.
        }

        bitmap_ctx->local_num_words++;
    }

    // we can't compare to the reference if it is unaligned: we store the seqeuence in nonref without an indication in the bitmap
    if (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')) {
        buf_add (&nonref_ctx->local, seq, seq_len);
        goto align_nonref_local; 
    }

    if (seq[0] == '*') {
        vb->md_verified = false;    
        vb->seq_missing = true;
        goto done; // we already handled a missing seq (SEQ="*") by adding a '-' to CIGAR - no data added here
    }

    bool no_cigar = cigar[0] == '*' && cigar[1] == 0; // there's no CIGAR. 

    RefLock lock = REFLOCK_NONE;
    Range *range = no_cigar ? NULL : ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, seq, 
                                                               flag.reference == REF_EXTERNAL ? NULL : &lock);

    // Cases where we don't consider the refernce and just copy the seq as-is
    if (!range) { // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
                  // 2. (loaded:) case contig doesn't exist in the reference, or POS is out of range of contig (observed in the wild with chrM)
                  // 3. no cigar - the sequence is not aligned to the reference even if we have RNAME and POS (and its length can exceed the reference contig)
        buf_add (&nonref_ctx->local, seq, seq_len);
        
        bit_array_clear_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // note: vb->ref_and_seq_consumed==0 if cigar="*"
        bitmap_ctx->next_local += ref_and_seq_consumed;
        vb->mismatch_bases     += ref_and_seq_consumed;
        vb->md_verified = false; // we can't use MD special alg if there is no range

        random_access_update_last_pos (VB, 0, pos + ref_consumed - 1);

        // note: in case of a missing range (which can be the first range in this seq, or a subsequent range), we zero the entire remaining bitmap.
        // this is because, absent a reference, we don't know how much ref is consumed by this missing range.
        goto align_nonref_local; 
    }    

    uint32_t pos_index     = pos - range->first_pos;
    uint32_t next_ref      = pos_index;
    rom next_cigar = cigar;

    uint32_t i=0;
    int subcigar_len=0;
    char cigar_op;

    uint32_t ref_len_this_level = (flag.reference == REF_INTERNAL ? MIN_(ref_consumed, range->last_pos - pos + 1)
                                                                  : ref_consumed); // possibly going around the end of the chromosome in case of a circular chromosome                                   

    uint32_t range_len = (range->last_pos - range->first_pos + 1);
    
    while (i < seq_len || next_ref < pos_index + ref_len_this_level) {

        ASSERT0 (i <= seq_len && next_ref <= pos_index + ref_len_this_level, "i or next_ref are out of range");

        subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
        cigar_op = *(next_cigar++);

        if (cigar_op == 'M' || cigar_op == '=' || cigar_op == 'X') { // alignment match or sequence match or mismatch

            ASSERT (subcigar_len > 0 && subcigar_len <= (seq_len - i), 
                    "CIGAR implies seq_len longer than actual seq_len=%u (recursion_level=%u level0: cigar=%s seq_len=%u). CIGAR=\"%s\"", 
                    seq_len, recursion_level, level_0_cigar, level_0_seq_len, cigar);

            uint32_t bit_i = bitmap_ctx->next_local; // copy to automatic variable for performance
            uint32_t start_i = i;
            while (subcigar_len && next_ref < pos_index + ref_len_this_level) {

                // when we have an X we don't enter it into our internal ref, and we wait for a read with a = or M for that site,
                // as we assume that on average, more reads will have the reference base, leading to better compression
            
                bool normal_base = IS_NUCLEOTIDE (seq[i]);

                // circle around to beginning of chrom if out of range (can only happen with external reference, expected only with circular chromosomes) 
                uint32_t actual_next_ref = next_ref % range_len; 

                // case: we have not yet set a value for this site - we set it now. note: in ZIP, is_set means that the site
                // will be needed for pizzing. With REF_INTERNAL, this is equivalent to saying we have set the ref value for the site
                if (flag.reference == REF_INTERNAL && range && normal_base 
                    && !ref_is_nucleotide_set (range, actual_next_ref)) { 

                    ref_set_nucleotide (range, actual_next_ref, seq[i]);
                    bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                    bit_i++; // bit remains set. cannot increment inside the macro
                }

                // case our seq is identical to the reference at this site
                else if (range && normal_base && seq[i] == ref_base_by_idx (range, actual_next_ref)) {
                    bit_i++; // bit remains set. 

                    if (flag.reference == REF_EXT_STORE)
                        bit_array_set (&range->is_set, actual_next_ref); // we will need this ref to reconstruct
                }
                
                // case: ref is set to a different value - we store our value in nonref_ctx
                else {
                    BNXTc (nonref_ctx->local) = seq[i];
                    bit_array_clear (bitmap, bit_i); bit_i++;
                    vb->mismatch_bases++;
                } 

                subcigar_len--;
                next_ref++;
                i++;
            }
            ref_and_seq_consumed -= (i - start_i); // update in case a range in a subsequent recursion level is missing and we need to clear the bitmap
            bitmap_ctx->next_local = bit_i;
        } // end if 'M', '=', 'X'

        // for Insertion or Soft clipping - this SEQ segment doesn't align with the reference - we leave it as is 
        else if (cigar_op == 'I' || cigar_op == 'S') {

            ASSSEG (subcigar_len > 0 && subcigar_len <= (seq_len - i), seq,
                    "CIGAR %s implies seq_len longer than actual seq_len=%u", cigar, seq_len);

            buf_add (&nonref_ctx->local, &seq[i], subcigar_len);
            i += subcigar_len;
            subcigar_len = 0;
        }

        // for Deletion or Skipping - we move the next_ref ahead
        else if (cigar_op == 'D' || cigar_op == 'N') {
            unsigned ref_consumed_skip = (flag.reference == REF_INTERNAL ? MIN_(subcigar_len, range_len - next_ref)
                                                                         : subcigar_len);
            next_ref     += ref_consumed_skip;
            subcigar_len -= ref_consumed_skip;
        }

        // Hard clippping (H) or padding (P) - nothing much to do
        else if (cigar_op == 'H' || cigar_op == 'P') 
            subcigar_len = 0;

        else {
            ASSSEG (cigar_op, vb->last_cigar, "End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                    "(cigar=%s pos=%"PRId64" recursion_level=%u level_0_cigar=%s level_0_seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                    pos_index + ref_len_this_level - next_ref, seq_len-i,   cigar, pos, recursion_level, level_0_cigar, level_0_seq_len,
                    ref_consumed, next_ref, pos_index, ref_len_this_level, subcigar_len, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            ASSSEG (false, vb->last_cigar, "Invalid CIGAR op: '%c' (ASCII %u)", cigar_op, cigar_op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_len_this_level && subcigar_len) break;
    }

    ref_unlock (gref, lock); // does nothing if REFLOCK_NONE

    uint32_t this_seq_last_pos = pos + (next_ref - pos_index) - 1;

    // in REF_INTERNAL, the sequence can flow over to the next range as each range is 1M bases. this cannot happen
    // in REF_EXTERNAL as each range is the entire contig
    ASSERT (flag.reference == REF_INTERNAL || i == seq_len, "expecting i(%u) == seq_len(%u) pos=%"PRId64" range=[%.*s %"PRId64"-%"PRId64"] (cigar=%s recursion_level=%u level0: cigar=%s seq_len=%u)", 
            i, seq_len, pos, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos, cigar, recursion_level, level_0_cigar, level_0_seq_len);

    // case: we have reached the end of the current reference range, but we still have sequence left - 
    // call recursively with remaining sequence and next reference range 
    if (i < seq_len) {

        ASSSEG (this_seq_last_pos <= MAX_POS, cigar, "POS=%"PRId64" and the consumed reference implied by CIGAR=\"%s\", exceeding MAX_POS=%"PRId64
                " (next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                pos, cigar, MAX_POS, next_ref, pos_index, ref_len_this_level, subcigar_len, 
                range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);

        char updated_cigar[strlen (next_cigar) + 20];
        if (subcigar_len) sprintf (updated_cigar, "%u%c%s", subcigar_len, cigar_op, next_cigar);

        sam_seg_SEQ (vb, bitmap_did, seq + i, seq_len - i, range->last_pos + 1, subcigar_len ? updated_cigar : next_cigar, 
                     is_revcomp, ref_consumed - ref_len_this_level, ref_and_seq_consumed,
                     recursion_level + 1, level_0_seq_len, level_0_cigar, 0);
    }
    else { // update RA of the VB with last pos of this line as implied by the CIGAR string
        if (this_seq_last_pos <= range->last_pos) // always the case in INTERNAL and non-circular EXTERNAL 
            random_access_update_last_pos (VB, 0, this_seq_last_pos);

        else  // we circled back to the beginning for the chromosome - i.e. this VB RA is the entire chromosome
            random_access_update_to_entire_chrom (VB, 0, range->first_pos, range->last_pos);
    }

    // final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (!recursion_level) {
        BitArray *M_is_ref = buf_get_bitarray (&vb->md_M_is_ref);

        bool bitmap_matches_MD = vb->md_verified && !bit_array_hamming_distance (M_is_ref, 0, bitmap, bitmap_start, M_is_ref->nbits);

        if (flag.show_wrong_md && vb->md_verified && !bitmap_matches_MD) {
            iprintf ("vb=%u line=%d RNAME=%.*s POS=%"PRId64" CIGAR=%s MD=%.*s SEQ=%.*s\n", 
                    vb->vblock_i, vb->line_i, STRf(vb->chrom_name), pos, vb->last_cigar, vb->last_txt_len(OPTION_MD_Z), last_txt(vb, OPTION_MD_Z), STRf(seq));
            bit_array_print_substr ("SEQ match to ref", bitmap, bitmap_start, M_is_ref->nbits, info_stream);
            bit_array_print_substr ("MD implied match", M_is_ref, 0, M_is_ref->nbits, info_stream); 
        }

        vb->md_verified = bitmap_matches_MD;
    }

align_nonref_local: {
    // we align nonref_ctx->local to a 4-character boundary. this is because CODEC_ACGT squeezes every 4 characters into a byte,
    // before compressing it with LZMA. In sorted SAM, we want subsequent identical sequences to have the same byte alignment
    // so that LZMA can catch their identicality. (note: as of v13, we assign a codec rather than hard-coding LZMA, but it is usually LZMA anyway)
    uint64_t add_chars = (4 - (nonref_ctx->local.len & 3)) & 3;
    if (add_chars) buf_add (&nonref_ctx->local, "AAA", add_chars); // add 1 to 3 As
}
done:
    COPY_TIMER (sam_seg_SEQ);
}


// converts native SAM/BAM format to 2bit ACGT - if soft_fail, returns false if any base is not A,C,G or T 
bool sam_sa_native_to_actg (VBlockSAMP vb, BitArray *packed, uint64_t next_bit, STRp(seq), bool bam_format)
{
    if (bam_format) 
        for (uint32_t i=0; i < seq_len; i++, next_bit += 2) {
            uint8_t b = (!(i&1)) ? (((uint8_t*)seq)[i>>1] >> 4) : (((uint8_t*)seq)[i>>1] & 0xf);
            switch (b) {
                case 0b0001 : bit_array_assign2 (packed, next_bit, 0); break;
                case 0b0010 : bit_array_assign2 (packed, next_bit, 1); break;
                case 0b0100 : bit_array_assign2 (packed, next_bit, 2); break;
                case 0b1000 : bit_array_assign2 (packed, next_bit, 3); break;
                default     : ABORT ("Unexpected base: '%c' vb=%u i=%u seq_len=%u", bam_base_codes[b], vb->vblock_i, i, seq_len);
            }
        }    

    else  // SAM
        for (uint32_t i=0; i < seq_len; i++, next_bit += 2) 
            switch (seq[i]) {
                case 'A' : bit_array_assign2 (packed, next_bit, 0); break;
                case 'C' : bit_array_assign2 (packed, next_bit, 1); break;
                case 'G' : bit_array_assign2 (packed, next_bit, 2); break;
                case 'T' : bit_array_assign2 (packed, next_bit, 3); break;
                default  : ABORT ("Unexpected base: '%c'(ASCII %u) vb=%u i=%u seq_len=%u", seq[i], (uint8_t)seq[i], vb->vblock_i, i, seq_len);
            }

    bit_array_clear_excess_bits_in_top_word (packed);

    return true;
}

// pack seq (into 2bit ACGT format) of each PRIM line separately 
void sam_zip_prim_ingest_vb_pack_seq (VBlockSAMP vb, SAGroupType *vb_grps, uint32_t vb_grps_len,
                                      Buffer *underlying_buf, Buffer *packed_seq_buf, bool is_bam_format)
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

        sam_sa_native_to_actg (vb, sa_seq, next_bit, Bc (vb->txt_data, vb_grp->seq), vb_grp->seq_len, is_bam_format);
        vb_grp->seq = next_bit / 2; // update from an index into txt_data to an index (bases not bits) into sa_seq
    }
}

// add the diff between our seq and prim to seq to SAM_SEQSA.local. SAM_SPECIAL_DEPN_SEQ is implicitly
// segged to SAM_SQBITMAP - this is an all-the-same context initialized in sam_seg_initialize. 
static void sam_seg_depn_SEQ (VBlockSAMP vb, STRp(seq), bool is_revcomp) // textual (SAM) format
{
    bool xstrand = (is_revcomp != vb->sa_grp->revcomp); // primary and dependent are on opposite strands

    // convert to ACGT format
    ASSERTNOTINUSE (vb->scratch);
    BitArray *depn = buf_alloc_bitarr (vb, &vb->scratch, seq_len * 2, "scratch");

    sam_sa_native_to_actg (vb, depn, 0, STRa(seq), false); // we already verified that it is ACGT in sam_sa_seg_depn_find_sagroup
    
    if (xstrand)
        bit_array_reverse_complement_in_place (depn);

    // copy subset of prim seq (TODO: this is wasteful, develop bit_array_xor to xor bitarr ranges directly)
    BitArray *sa_seq = buf_get_bitarray (&z_file->sa_seq);
    
    ASSERTNOTINUSE (vb->codec_bufs[0]);
    BitArray *prim = buf_alloc_bitarr (vb, &vb->codec_bufs[0], seq_len * 2, "codec_bufs[0]");
    
    uint64_t prim_start_base = xstrand ? vb->sa_grp->seq + vb->hard_clip[1]
                                       : vb->sa_grp->seq + vb->hard_clip[0];

    bit_array_copy (prim, 0, sa_seq, 2 * prim_start_base, 2 * seq_len);

    bit_array_xor (depn, depn, prim);

    Buffer *seqdepn_buf = &CTX(SAM_SEQSA)->local;
    seqdepn_buf->vb     = VB; // needed by bit_array_concat
    seqdepn_buf->name   = "contexts->local";
    BitArray *seqdepn   = buf_get_bitarray (seqdepn_buf);
    
    bit_array_concat (seqdepn, depn, 0);

    // in DEPN component - we seg an all-the-same SPECIAL that either reconstructs SEQ from 
    // SQBITMAP.local or from sam_piz_special_pull_from_SAGROUP. 
    seg_by_did_i (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_DEPN_SEQ }, 2, SAM_SQBITMAP, 0); // seg and not just create_node because we have a local section too - so we need a b250 section not just dict

    buf_free (vb->scratch);
    buf_free (vb->codec_bufs[0]);
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

// PIZ: SEQ reconstruction 
void sam_reconstruct_SEQ (VBlockP vb_, Context *bitmap_ctx, rom unused, unsigned unused2)
{
#define ROUNDUP_TO_NEAREST_4(x) ((uint32_t)(x) + 3) & ~((uint32_t)0x3)

    VBlockSAMP vb = (VBlockSAMP)vb_;

    if (!bitmap_ctx->is_loaded) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx   = CTX(SAM_NONREF);
    rom nonref            = Bc (nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    rom nonref_start      = nonref;
    unsigned subcigar_len = 0;
    char cigar_op         = 0;
    const PosType pos     = vb->last_int(SAM_POS);
    ConstRangeP range     = NULL;
    unsigned seq_consumed=0, ref_consumed=0;

    // case: unaligned sequence - pos is 0 
    if (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')) {
        // case: compressed with a reference, using our aligner
        if (z_file->z_flags.aligner) {
            aligner_reconstruct_seq (VB, bitmap_ctx, vb->seq_len, false);
            nonref_ctx->next_local = ROUNDUP_TO_NEAREST_4 (nonref_ctx->next_local);
        }
        // case: no reference was used - in this case, the sequence is not encoded in the bitmap at all. we just copy it from NONREF
        else {
            RECONSTRUCT (nonref, vb->seq_len); 
            nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (vb->seq_len);
        }
        return;
    }

    // case: missing sequence - sequence is '*' - just reconstruct the '*' (but only if it is the primary SEQ field,
    // which we can tell by not being encountered earlier on the line - bc in v8 (at least) E2 was also stored in the same SAM_SQBITMAP)
    if (vb->seq_missing && !ctx_encountered_in_line (VB, bitmap_ctx->did_i)) {
        RECONSTRUCT1 ('*');
        return;
    }

    rom next_cigar = B1STc(vb->textual_cigar); 
    
    ASSERTNOTEMPTY (bitmap_ctx->local); // make sure its not skipped

    // in an edge case, when all is_set bits are zero, the range might not even be written to the file
    bool uses_ref_data = vb->ref_consumed && 
                         (bit_array_num_bits_set_region (buf_get_bitarray (&bitmap_ctx->local), bitmap_ctx->next_local, vb->ref_and_seq_consumed) > 0);
    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence

    range = uses_ref_data ? ref_piz_get_range (VB, gref, true) : NULL; 
    vb->range = (RangeP)range; // store is vb->range too, for sam_piz_special_MD

    PosType range_len = range ? (range->last_pos - range->first_pos + 1) : 0;

    while (seq_consumed < vb->seq_len || ref_consumed < vb->ref_consumed) {
        
        if (!subcigar_len) {
            subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
            cigar_op = cigar_lookup_sam[(uint8_t)*(next_cigar++)];
            ASSERT (cigar_op, "Invalid CIGAR op while reconstructing line %d: '%c' (ASCII %u)", vb->line_i, *(next_cigar-1), *(next_cigar-1));
            cigar_op &= 0x0f; // remove validity bit
        }

        if (cigar_op & CIGAR_CONSUMES_QUERY) {

            if ((cigar_op & CIGAR_CONSUMES_REFERENCE) && NEXTLOCALBIT (bitmap_ctx)) /* copy from reference */ {

                if (!vb->drop_curr_line) { // note: if this line is excluded with --regions, then the reference section covering it might not be loaded
                    if (!range) {
                        ref_display_all_ranges (gref);
                        ASSPIZ0 (false, "range is NULL");
                    }
                    
                    uint32_t idx = ((pos - range->first_pos) + ref_consumed) % range_len; // circle around (this can only happen when compressed with an external reference)

                    if (!ref_is_nucleotide_set (range, idx) &&
                         (!flag.regions || regions_is_site_included (VB))) { // if this line is not included, then possibly its reference range is not loaded. we complete consumption (resulting in bad reconstruction) and drop the line in container_reconstruct_do

                        ref_print_is_set (range, pos + ref_consumed, stderr);
                        ASSPIZ (false, "reference is not set: chrom=%u \"%.*s\" pos=%"PRId64" range=[%"PRId64"-%"PRId64"]"
                                " (cigar=%s seq_start_pos=%"PRId64" ref_consumed=%u seq_consumed=%u)",
                                range->chrom, range->chrom_name_len, range->chrom_name, pos + ref_consumed, 
                                range->first_pos, range->last_pos, B1STc (vb->textual_cigar), pos, ref_consumed, seq_consumed);
                    }

                    char ref = ref_base_by_idx (range, idx);
                    RECONSTRUCT1 (ref); 
                }
            }
            else {
                vb->mismatch_bases += (cigar_op & CIGAR_CONSUMES_REFERENCE) > 0;
                RECONSTRUCT1 (*nonref++);
            }

            seq_consumed++;
        }

        if (cigar_op & CIGAR_CONSUMES_REFERENCE) 
            ref_consumed++;

        subcigar_len--;
    }

    ASSPIZ (seq_consumed == vb->seq_len,      "expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSPIZ (ref_consumed == vb->ref_consumed, "expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (nonref - nonref_start);
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
    if (l_seq % 2) *BAFTc (vb->txt_data) = 0; 

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
        str_revcomp (recon, recon_len);

    return 0;
}

// reconstruct from a 2bit array - start_base and seq_len are in bases (not in bits)
static void reconstruct_2bit_sequence (VBlockP vb, BitArrayP seq_2bits, uint64_t start_base, uint32_t seq_len)
{
    // Reconstruct as SAM. If BAM, we later translate to BAM. To do: reconstruct directly as BAM if needed, bug 530
    char *next = BAFTc(vb->txt_data);
    
    for (uint32_t i=0; i < seq_len; i++) {
        uint8_t b = bit_array_get2 (seq_2bits, (start_base + i)*2);
        *next++ = b==0 ? 'A' : b==1 ? 'C' : b==2 ? 'G' : 'T';
    }

    vb->txt_data.len += seq_len;
}

// PIZ of a SEQ in a DEPN component - this is an all-the-same contexts - SEQ all lines in a DEPN component
// come here. We multiplex between SAGroup-based reconstruction vs normal SEQ reconstruction
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_DEPN_SEQ)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    // case: normal reconstruction, not against PRIM
    if (!SAM_PIZ_HAS_SA_GROUP) {
        sam_reconstruct_SEQ (VB, ctx, 0, 0);
        return false;
    }

    const SAGroupType *grp = vb->sa_grp;

    ASSPIZ (grp->seq_len > vb->hard_clip[0] + vb->hard_clip[1], "grp_i=%u aln_i=%"PRIu64" grp->seq_len=%u <= hard_clip[LEFT]=%u + hard_clip[RIGHT]=%u", 
            ZGRP_I(grp), ZALN_I(vb->sa_aln), grp->seq_len, vb->hard_clip[0], vb->hard_clip[1]);

    uint32_t depn_seq_len = grp->seq_len - vb->hard_clip[0] - vb->hard_clip[1]; // depn sequence is a sub-sequence of the prim sequence

    ContextP seqdepn_ctx = CTX(SAM_SEQSA); // initially, contains XOR diff
    const BitArray *seqdepn = buf_get_bitarray (&seqdepn_ctx->local);

    if (!reconstruct || !ctx->is_loaded) 
        goto done;

    BitArray *all_prim = buf_get_bitarray (&z_file->sa_seq); // SEQ of entire PRIM component

    // get prim (with flanking regions removed according to depn's hard clipping)
    ASSERTNOTINUSE (vb->codec_bufs[0]);
    BitArray *prim = buf_alloc_bitarr (vb, &vb->codec_bufs[0], depn_seq_len * 2, "codec_bufs[0]");
    bool xstrand = vb->sa_aln->revcomp != grp->revcomp; // primary and dependent are on opposite strands

    // remove prim hard-clipped flanking regions
    if (!xstrand)
        bit_array_copy (prim, 0, all_prim,  2 * (grp->seq + vb->hard_clip[0]), 2 * depn_seq_len);
    else
        bit_array_copy (prim, 0, all_prim, 2 * (grp->seq + grp->seq_len - depn_seq_len - vb->hard_clip[0]), 2 * depn_seq_len);

    bit_array_truncate (prim, 2 * depn_seq_len);

    // build sequence from PRIM sequence XOR diff
    ASSERTNOTINUSE (vb->scratch);
    ASSPIZ (2 * (seqdepn_ctx->next_local + depn_seq_len) <= seqdepn->nbits, 
            "SEQ bitmap out of range in grp_i=%u grp_seq=(%"PRIu64", %u) hard-clip=(%u,%u): seqdepn.next_local=%u + depn_seq_len=%u > seqdepn(xor).bases=%"PRIu64, 
            ZGRP_I(grp), grp->seq, grp->seq_len, vb->hard_clip[0], vb->hard_clip[1], seqdepn_ctx->next_local, depn_seq_len, seqdepn->nbits/2);

    BitArray *depn = buf_alloc_bitarr (vb, &vb->scratch, depn_seq_len * 2, "scratch");
    bit_array_copy (depn, 0, seqdepn, 2 * seqdepn_ctx->next_local, 2 * depn_seq_len); // this is prim XOR depn 
    bit_array_xor (depn, depn, prim); // this is now depn in actg format

    // reverse complement if needed
    if (xstrand) 
        bit_array_reverse_complement_in_place (depn);

    // Reconstruct as SAM. If BAM, we later translate to BAM. To do: reconstruct directly as BAM if needed, bug 530
    reconstruct_2bit_sequence (VB, depn, 0, depn_seq_len);

    buf_free (vb->scratch);
    buf_free (vb->codec_bufs[0]);

done:
    seqdepn_ctx->next_local += depn_seq_len; // in bases
    return false; // no new value
}

// PIZ of a SEQ in a PRIM component (an all-the-same context) - copy SEQ from the in-memory SA Group
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_PRIM_SEQ)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    sam_piz_set_sa_grp (vb);

    BitArray *all_prim = buf_get_bitarray (&z_file->sa_seq); // SEQ of entire PRIM component    
    reconstruct_2bit_sequence (VB, all_prim, vb->sa_grp->seq, vb->sa_grp->seq_len);
    
    return false; // no new value
}

