// ------------------------------------------------------------------
//   sam_seq.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
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

// Creates a bitmap from seq data - exactly one bit per base that is mapped to the reference (e.g. not for INSERT bases)
// - Normal SEQ: tracking CIGAR, we compare the sequence to the reference, indicating in a SAM_SQBITMAP whether this
//   base in the same as the reference or not. In case of REF_INTERNAL, if the base is not already in the reference, we add it.
//   bases that differ from the reference are stored in SAM_NONREF
// - Edge case: no POS (i.e. unaligned read) - we just store the sequence in SAM_NONREF
// - Edge case: no CIGAR (it is "*") - we just treat it as an M and compare to the reference
// - Edge case: no SEQ (it is "*") - we '*' in SAM_NONREF and indicate "different from reference" in the bitmap. We store a
//   single entry, regardless of the number of entries indicated by CIGAR
void sam_seg_SEQ (VBlockSAM *vb, DidIType bitmap_did, STRp(seq), const PosType pos, const char *cigar, 
                  uint32_t ref_consumed, uint32_t ref_and_seq_consumed,
                  unsigned recursion_level, uint32_t level_0_seq_len, const char *level_0_cigar, unsigned add_bytes)
{
    START_TIMER;

    Context *bitmap_ctx = CTX(bitmap_did);
    Context *nonref_ctx = bitmap_ctx + 1;

    ASSERT (recursion_level < 500, "excess recursion recursion_level=%u seq_len=%u level_0_cigar=%s", // a large number of recursion calls can happen if CIGAR=9M10910N86M3274690N30M1S as observed with the STAR aligner https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            recursion_level, seq_len, level_0_cigar);

    bitmap_ctx->txt_len += add_bytes; // byte counts for --show-sections

    // for unaligned lines, if we have refhash loaded, use the aligner instead of CIGAR-based segmenting
    if (!pos && segconf.sam_use_aligner) {
        aligner_seg_seq (VB, bitmap_ctx, STRa(seq));
        goto align_nonref_local;
    }

    BitArray *bitmap = buf_get_bitarray (&bitmap_ctx->local);
    uint32_t bitmap_start = bitmap_ctx->next_local;        

    if (!recursion_level) {

        ASSERT (seq_len < 100000 || segconf_is_long_reads(), 
                 "Warning: sam_seg_SEQ: seq_len=%u is suspiciously high and might indicate a bug. vb=%u line_i=%"PRIu64, seq_len, vb->vblock_i, vb->line_i);

        buf_alloc (vb, &bitmap_ctx->local, roundup_bits2words64 (ref_and_seq_consumed), 
                   (segconf_is_long_reads() ? vb->lines.len/2 : vb->lines.len) * (MIN_(ref_and_seq_consumed, segconf.sam_seq_len)) / 64,  // this formula was chosen after extensive trial & error - easy to get wrong, esp for long reads
                   uint64_t /* len is in 64b words in a bitmap */, CTX_GROWTH, "contexts->local"); 

        buf_extend_bits (&bitmap_ctx->local, ref_and_seq_consumed); 
        bit_array_set_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // we initiaze all the bits to "set", and clear as needed.

        buf_alloc (vb, &nonref_ctx->local, seq_len + 3, 
                   vb->lines.len * (MIN_(seq_len, segconf.sam_seq_len)) / (segconf_is_long_reads() ? 50 : 20), // this one too - it was chosen after extensive trial & error - easy to get wrong, esp for long reads
                   uint8_t, CTX_GROWTH, "contexts->local"); 
    
        if (segconf.running) return; // case segconf: we created the contexts for segconf_set_vb_size accounting. that's enough - actually segging will mark is_set and break sam_md_analyze.
    }

    // we can't compare to the reference if it is unaligned: we store the seqeuence in nonref without an indication in the bitmap
    if (!pos || (vb->chrom_name_len==1 && vb->chrom_name[0]=='*')) {
        buf_add (&nonref_ctx->local, seq, seq_len);
        goto align_nonref_local; 
    }

    if (seq[0] == '*') {
        vb->md_verified = false;    
        goto done; // we already handled a missing seq (SEQ="*") by adding a '-' to CIGAR - no data added here
    }

    bool no_cigar = cigar[0] == '*' && cigar[1] == 0; // there's no CIGAR. 

    RefLock lock;
    Range *range = no_cigar ? NULL : ref_seg_get_locked_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), pos, ref_consumed, WORD_INDEX_NONE, seq, &lock);

    // Cases where we don't consider the refernce and just copy the seq as-is
    if (!range) { // 1. (denovo:) this contig defined in @SQ went beyond the maximum genome size of 4B and is thus ignored
                  // 2. (loaded:) case contig doesn't exist in the reference, or POS is out of range of contig (observed in the wild with chrM)
                  // 3. no cigar - the sequence is not aligned to the reference even if we have RNAME and POS (and its length can exceed the reference contig)
        buf_add (&nonref_ctx->local, seq, seq_len);
        
        bit_array_clear_region (bitmap, bitmap_ctx->next_local, ref_and_seq_consumed); // note: vb->ref_and_seq_consumed==0 if cigar="*"
        bitmap_ctx->next_local += ref_and_seq_consumed;
        vb->mismatch_bases     += ref_and_seq_consumed;

        random_access_update_last_pos (VB, DC_PRIMARY, pos + ref_consumed - 1);

        // note: in case of a missing range (which can be the first range in this seq, or a subsequent range), we zero the entire remaining bitmap.
        // this is because, absent a reference, we don't know how much ref is consumed by this missing range.
        goto align_nonref_local; 
    }    

    uint32_t pos_index     = pos - range->first_pos;
    uint32_t next_ref      = pos_index;
    const char *next_cigar = cigar;

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
                    "CIGAR %s implies seq_len longer than actual seq_len=%u (recursion_level=%u level0: cigar=%s seq_len=%u)", 
                    cigar, seq_len, recursion_level, level_0_cigar, level_0_seq_len);

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
                    NEXTENT (char, nonref_ctx->local) = seq[i];
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
        else if (cigar_op == 'H' || cigar_op == 'P') {
            subcigar_len = 0;
        }

        else {
            ASSSEG (cigar_op, vb->last_cigar, "Error in sam_seg_SEQ: End of CIGAR reached but we still have %u reference and %u sequence bases to consume"
                    "(cigar=%s pos=%"PRId64" recursion_level=%u level_0_cigar=%s level_0_seq_len=%u) (ref_consumed=%d next_ref=%u pos_index=%u ref_len_this_level=%u subcigar_len=%u range=[%.*s %"PRId64"-%"PRId64"])",
                    pos_index + ref_len_this_level - next_ref, seq_len-i,   cigar, pos, recursion_level, level_0_cigar, level_0_seq_len,
                    ref_consumed, next_ref, pos_index, ref_len_this_level, subcigar_len, range->chrom_name_len, range->chrom_name, range->first_pos, range->last_pos);        

            ASSSEG (false, vb->last_cigar, "Invalid CIGAR op: '%c' (ASCII %u)", cigar_op, cigar_op);        
        }

        // case: we're at the end of the reference AND we want more of it
        if (next_ref == pos_index + ref_len_this_level && subcigar_len) break;
    }

    if (range) ref_unlock (gref, lock);       

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
                           ref_consumed - ref_len_this_level, ref_and_seq_consumed,
                           recursion_level + 1, level_0_seq_len, level_0_cigar, 0);
    }
    else { // update RA of the VB with last pos of this line as implied by the CIGAR string
        if (this_seq_last_pos <= range->last_pos) // always the case in INTERNAL and non-circular EXTERNAL 
            random_access_update_last_pos (VB, DC_PRIMARY, this_seq_last_pos);

        else  // we circled back to the beginning for the chromosome - i.e. this VB RA is the entire chromosome
            random_access_update_to_entire_chrom (VB, DC_PRIMARY, range->first_pos, range->last_pos);
    }

    // final verification step - does MD:Z correctly reflect matches and mismatches of M/X/=
    if (!recursion_level) {
        BitArray *M_is_ref = buf_get_bitarray (&vb->md_M_is_ref);

        bool bitmap_matches_MD = vb->md_verified && !bit_array_hamming_distance (M_is_ref, 0, bitmap, bitmap_start, M_is_ref->nbits);

        if (flag.show_wrong_md && vb->md_verified && !bitmap_matches_MD) {
            iprintf ("vb=%u line=%"PRIu64" RNAME=%.*s POS=%"PRId64" CIGAR=%s MD=%.*s SEQ=%.*s\n", 
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

// used by codec_enano_compress
COMPRESSOR_CALLBACK (sam_zip_seq) 
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    *line_data_len = dl->seq_len;
    *line_data     = ENT (char, vb->txt_data, dl->seq_data_start);

    if (is_rev) *is_rev = dl->FLAG.bits.rev_comp;
}

//---------
// PIZ
//---------

// PIZ: SEQ reconstruction 
void sam_reconstruct_SEQ (VBlock *vb_, Context *bitmap_ctx, const char *unused, unsigned unused2)
{
#define ROUNDUP_TO_NEAREST_4(x) ((uint32_t)(x) + 3) & ~((uint32_t)0x3)

    VBlockSAMP vb = (VBlockSAMP)vb_;
    ASSERT0 (bitmap_ctx && bitmap_ctx->did_i == SAM_SQBITMAP, "context is not SAM_SQBITMAP");

    if (piz_is_skip_section (vb, SEC_LOCAL, bitmap_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx      = CTX(SAM_NONREF);
    const char *nonref       = ENT (const char, nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    const char *nonref_start = nonref;
    unsigned subcigar_len    = 0;
    char cigar_op            = 0;
    const PosType pos        = vb->last_int(SAM_POS);
    const Range *range       = NULL;
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

    // case: missing sequence - sequence is '*' (which zip marked with a '-' in the cigar) - just reconstruct the '*'
    if (*vb->last_cigar == '-') {
        RECONSTRUCT1 ('*');
        vb->last_cigar++; // skip the '-' so it doesn't affect a subsequent E2 on the same line
        return;
    }

    const char *next_cigar = vb->last_cigar; // don't change vb->last_cigar as we may still need it, eg if we have an E2 optional field
    
    // in an edge case, when all is_set bits are zero, the range might not even be written to the file
    bool uses_ref_data = vb->ref_consumed && 
                         (bit_array_num_bits_set_region (buf_get_bitarray (&bitmap_ctx->local), bitmap_ctx->next_local, vb->ref_and_seq_consumed) > 0);
    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence

    range = uses_ref_data ? ref_piz_get_range (vb_, gref, pos, vb->ref_consumed) : NULL; 
    vb->range = (RangeP)range; // store is vb->range too, for sam_piz_special_MD

    PosType range_len = range ? (range->last_pos - range->first_pos + 1) : 0;

    while (seq_consumed < vb->seq_len || ref_consumed < vb->ref_consumed) {
        
        if (!subcigar_len) {
            subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
            cigar_op = cigar_lookup_sam[(uint8_t)*(next_cigar++)];
            ASSERT (cigar_op, "Invalid CIGAR op while reconstructing line %"PRIu64": '%c' (ASCII %u)", vb->line_i, *(next_cigar-1), *(next_cigar-1));
            cigar_op &= 0x0f; // remove validity bit
        }

        if (cigar_op & CIGAR_CONSUMES_QUERY) {

            if ((cigar_op & CIGAR_CONSUMES_REFERENCE) && NEXTLOCALBIT (bitmap_ctx)) /* copy from reference */ {

                if (!vb->drop_curr_line) { // note: if this line is excluded with --regions, then the reference section covering it might not be loaded
                    uint32_t idx = ((pos - range->first_pos) + ref_consumed) % range_len; // circle around (this can only happen when compressed with an external reference)

                    if (!ref_is_nucleotide_set (range, idx) &&
                         (!flag.regions || regions_is_site_included (VB))) { // if this line is not included, then possibly its reference range is not loaded. we complete consumption (resulting in bad reconstruction) and drop the line in container_reconstruct_do

                        ref_print_is_set (range, pos + ref_consumed, stderr);
                        ABORT ("Error in sam_reconstruct_SEQ: while reconstructing line %"PRIu64" (vb_i=%u: last_txt_line=%"PRIu64" num_lines=%"PRIu64"): reference is not set: chrom=%u \"%.*s\" pos=%"PRId64" range=[%"PRId64"-%"PRId64"]"
                               " (cigar=%s seq_start_pos=%"PRId64" ref_consumed=%u seq_consumed=%u)",
                               vb->line_i, vb->vblock_i, (vb->first_line + vb->lines.len - 1), vb->lines.len, 
                               range->chrom, range->chrom_name_len, range->chrom_name, pos + ref_consumed, 
                               range->first_pos, range->last_pos, vb->last_cigar, pos, ref_consumed, seq_consumed);
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

    ASSERT (seq_consumed == vb->seq_len,      "expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSERT (ref_consumed == vb->ref_consumed, "expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (nonref - nonref_start);
}

// SAM-to-BAM translator: translate SAM ASCII sequence characters to BAM's 4-bit characters:
TRANSLATOR_FUNC (sam_piz_sam2bam_SEQ)
{
    // the characters "=ACMGRSVTWYHKDBN" are mapped to BAM 0->15, in this matrix we add 0x80 as a validity bit. All other characters are 0x00 - invalid
    static const uint8_t sam2bam_seq_map[256] = { ['=']=0x80, ['A']=0x81, ['C']=0x82, ['M']=0x83, ['G']=0x84, ['R']=0x85, ['S']=0x86, ['V']=0x87, 
                                                  ['T']=0x88, ['W']=0x89, ['Y']=0x8a, ['H']=0x8b, ['K']=0x8c, ['D']=0x8d, ['B']=0x8e, ['N']=0x8f };
    
    if (vb->drop_curr_line) return 0; // sequence was not reconstructed - nothing to translate
    
    BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
    uint32_t l_seq = LTEN32 (alignment->l_seq);

    // if l_seq=0, just remove the '*'
    if (!l_seq) {
        vb->txt_data.len--;
        return 0;
    }

    // if l_seq is odd, 0 the next byte that will be half of our last result byte
    if (l_seq % 2) *AFTERENT (char, vb->txt_data) = 0; 

    uint8_t *seq_before=(uint8_t *)recon, *seq_after=(uint8_t *)recon; 
    for (uint32_t i=0; i < (l_seq+1)/2; i++, seq_after++, seq_before += 2) {
        uint8_t base[2] = { sam2bam_seq_map[(uint8_t)seq_before[0]], sam2bam_seq_map[(uint8_t)seq_before[1]] };
        
        // check for invalid characters - issue warning (only once per execution), and make then into an 'N'
        for (unsigned b=0; b < 2; b++)
            if (!base[b] && !(b==1 && (i+1)*2 > l_seq)) {
                char printable[l_seq*2+1];
                WARN_ONCE ("Warning: when converting SAM sequence data to BAM (QNAME=%.*s): invalid character encodered, it will be converted as 'N': '%c' (ASCII %u) (this warning will appear only once). SEQ=%s seq_i=%u", 
                           vb->last_txt_len(SAM_QNAME), last_txt (vb, SAM_QNAME), base[b], base[b], str_to_printable (recon, l_seq, printable), i*2+b);
                base[b] = 0x0f;
            }

        *seq_after = (base[0] << 4) | (base[1] & 0x0f);
    }

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

