// ------------------------------------------------------------------
//   sam_shared.c
//   Copyright (C) 2020-2024 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "sam_private.h"
#include "strings.h"
#include "file.h"
#include "endianness.h"

const uint8_t aux_width[256] = { ['i']=4, ['I']=4, ['s']=2, ['S']=2, ['c']=1, ['C']=1, ['A']=1, ['f']=4 };

// table of valid cigar_op as defined in https://samtools.github.io/hts-specs/SAMv1.pdf
const uint8_t cigar_lookup_sam[256] = { // note: bit 4 (0x10) is set for all valid values
    ['0'...'9']=0x11,                   // digits
    ['H']=0x10, ['P']=0x10,             // consume neither
    ['I']=0x12, ['S']=0x12, ['*']=0x12, // consume query (seq) only. Note: '*' is when CIGAR is "151*" - alignment with no CIGAR but a SEQ
    ['D']=0x14, ['N']=0x14,             // consume reference only
    ['M']=0x16, ['=']=0x16, ['X']=0x16  // consume both query and reference
};

const uint8_t cigar_lookup_bam[16] = {  // note: bit 4 (0x10) is set for all valid values
    [5/*H*/]=0x10, [6/*P*/]=0x10,       // consume neither
    [1/*I*/]=0x12, [4/*S*/]=0x12,       // consume query (seq) only.
    [2/*D*/]=0x14, [3/*N*/]=0x14,       // consume reference only
    [0/*M*/]=0x16, [7/*=*/]=0x16, [8/*X*/]=0x16  // consume both query and reference
};

unsigned sam_vb_size (DataType dt) { return sizeof (VBlockSAM); }
unsigned sam_vb_zip_dl_size (void) { return sizeof (ZipDataLineSAM); }

// initialization of the line
void sam_reset_line (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    ASSERT (VB_DT(SAM) || VB_DT(BAM), "VB has wrong data type: %s", dt_name (vb->data_type));
    
    vb->textual_cigar.len = vb->binary_cigar.len = vb->binary_cigar.next = 0;
    vb->textual_seq.len = 0;
    vb->qual_missing = vb->seq_missing = vb->seq_is_monochar = vb->cigar_missing = vb->line_not_deepable = false;
    vb->XG.len = 0;
    vb->seq_len = 0;
    vb->ref_consumed = vb->ref_and_seq_consumed = 0;
    vb->soft_clip[0] = vb->soft_clip[1] = 0;
    vb->hard_clip[0] = vb->hard_clip[1] = 0;
    vb->deletions = vb->insertions = 0;
    vb->mismatch_bases_by_SEQ = vb->mismatch_bases_by_MD = 0;
    vb->saggy_line_i = vb->mate_line_i = NO_LINE;
    vb->saggy_is_prim = false;
    vb->meth_call.len32 = 0;
    vb->bisulfite_strand = 0;
    
    if (IS_PIZ) {
        vb->textual_seq_str = NULL;
        vb->aux_con = NULL;
        vb->chrom_node_index = WORD_INDEX_NONE;
        vb->chrom_name = "";
        vb->chrom_name_len = 0;
        vb->range = NULL;
        vb->deep_seq = (PizDeepSeq){}; 
        vb->seq_is_monochar = false;
        CTX(SAM_SQBITMAP)->line_sqbitmap.nbits = CTX(SAM_SQBITMAP)->line_sqbitmap.nwords = 0;
        CTX(OPTION_tp_B_ARR)->tp = (struct ctx_tp){};

        // make sure we have enough room for this line if translating. 
        // note: having this allocation here allows us to keep vb->translation.factor relatively small to avoid over-allocation
        if (!vb->translation.is_src_dt && buf_user_count (&vb->txt_data) == 1) // not overlaid which happens loading sag (sam_load_groups_add_solo_data())
            buf_alloc (vb, &vb->txt_data, vb->longest_line_len * 5, 0, char, 1.15, "txt_data");
    }

    else { // ZIP
        memset (&vb->first_idx, 0xff, (char*)&vb->after_idx - (char*)&vb->first_idx); // set all idx's to -1
        vb->md_verified = false;
        vb->auxs = NULL;
        vb->aux_lens = NULL;
        vb->n_auxs = 0;
        vb->md_M_is_ref.nbits = vb->md_M_is_ref.nwords = 0;
        vb->unconverted_bitmap.nbits = vb->unconverted_bitmap.nwords = 0;
    }
}

// calculate bin given an alignment covering [first_pos_0,last_pos_0) (0-based positions, half-closed, half-open)
// code adapted from https://samtools.github.io/hts-specs/SAMv1.pdf section 5.3
uint16_t bam_reg2bin (int32_t first_pos, int32_t last_pos)
{
    int32_t first_pos_0 = first_pos - 1; // -1 to make it 0-based
    int32_t last_pos_0  = last_pos  - 1; // -1 to make it 0-based

    // Note: I found actual files where the bin was calculated by "last_pos_0=last_pos-2" but other actuals files are
    // according to the formula above (or maybe I am missing something?)
    // THIS SHOULD NEVER BE CORRECTED - as SAM_SPECIAL_BIN snips out there in the wild rely on this formula (even if incorrect)
    // to reconstruct (in sam_cigar_special_CIGAR)

    if (first_pos_0>>14 == last_pos_0>>14) return ((1<<15)-1)/7 + (first_pos_0>>14);
    if (first_pos_0>>17 == last_pos_0>>17) return ((1<<12)-1)/7 + (first_pos_0>>17);
    if (first_pos_0>>20 == last_pos_0>>20) return ((1<<9 )-1)/7 + (first_pos_0>>20);
    if (first_pos_0>>23 == last_pos_0>>23) return ((1<<6 )-1)/7 + (first_pos_0>>23);
    if (first_pos_0>>26 == last_pos_0>>26) return ((1<<3 )-1)/7 + (first_pos_0>>26);
    return 0;
}

DisFlagsStr sam_dis_flags (SamFlags f)
{
    DisFlagsStr s;
    sprintf (s.s, "multi_segs=%u is_aligned=%u unmapped=%u:%u revcomp=%u:%u mate#=%u:%u secondary=%u failfilt=%u dup=%u supp=%u",
             f.multi_segs, f.is_aligned, f.unmapped, f.next_unmapped, f.rev_comp, f.next_rev_comp, 
             f.is_first, f.is_last, f.secondary, f.filtered, f.duplicate, f.supplementary);
    return s;
}

rom buddy_type_name (BuddyType bt)
{
    return bt==BUDDY_EITHER?"EITHER" : bt==BUDDY_MATE?"MATE" : bt==BUDDY_SAGGY?"SAGGY" : "INVALID_BUDDY_TYPE";    
}
