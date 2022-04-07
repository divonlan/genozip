// ------------------------------------------------------------------
//   sam_shared.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

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

void sam_vb_release_vb (VBlockSAMP vb)
{
    vb->plsg_i = 0;
    vb->last_cigar = NULL;
    vb->ref_consumed = vb->ref_and_seq_consumed = vb->soft_clip[0] = vb->soft_clip[1] = 0;
    vb->mismatch_bases = vb->hard_clip[0] = vb->hard_clip[1] = 0;
    vb->a_bases = vb->x_bases = vb->y_bases = 0;
    vb->a_index = vb->x_index = vb->y_index = 0;
    vb->md_verified = 0;
    vb->qual_codec_no_longr = false;
    vb->seq_missing = vb->cigar_missing = vb->check_for_gc = 0;
    vb->qual_missing = 0;
    vb->sa_grp = 0;
    vb->sa_aln = 0;
    vb->comp_qual_len = vb->comp_cigars_len = 0;
    vb->NM = 0;
    vb->NM_len = 0;
    vb->XG_inc_S = 0;
    vb->first_grp_i = 0;
    vb->sa_grp_line_i = 0;

    buf_free (vb->bd_bi_line);
    buf_free (vb->XG);
    buf_free (vb->textual_cigar);
    buf_free (vb->binary_cigar);
    buf_free (vb->textual_seq);
    buf_free (vb->md_M_is_ref);
    buf_free (vb->sa_groups);
    buf_free (vb->sa_alns);
    buf_free (vb->sa_prim_cigars);
    buf_free (vb->qname_hash);
    buf_free (vb->buddy_textual_cigars);
}

void sam_vb_destroy_vb (VBlockSAMP vb)
{
    buf_destroy (vb->bd_bi_line);
    buf_destroy (vb->XG);
    buf_destroy (vb->textual_cigar);
    buf_destroy (vb->binary_cigar);
    buf_destroy (vb->textual_seq);
    buf_destroy (vb->md_M_is_ref);
    buf_destroy (vb->sa_groups);
    buf_destroy (vb->sa_alns);
    buf_destroy (vb->sa_prim_cigars);
    buf_destroy (vb->qname_hash);
    buf_destroy (vb->buddy_textual_cigars);
}

// initialization of the line
void sam_reset_line (VBlockP vb_)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    vb->textual_cigar.len = VB_SAM->binary_cigar.len = 0;
    vb->seq_missing = vb->cigar_missing = false;
    vb->qual_missing = QUAL_NOT_MISSING;
    vb->XG.len = 0;
    vb->buddy_line_i = NO_BUDDY; 
    vb->ref_consumed         = 0;
    vb->ref_and_seq_consumed = 0;
    vb->mismatch_bases       = 0;
    vb->soft_clip[0] = vb->soft_clip[1] = 0;
    vb->hard_clip[0] = vb->hard_clip[1] = 0;

    if (command == PIZ) {
        vb->chrom_node_index = WORD_INDEX_NONE;
        vb->chrom_name = "";
        vb->chrom_name_len = 0;

        if (!vb->preprocessing) 
            buf_alloc_zero (vb, &CTX(SAM_CIGAR)->piz_ctx_specific_buf, 0, vb->lines.len, uint32_t, 0, "piz_ctx_specific_buf"); // initialize to exactly one per line.
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

DisFlagsStr sam_dis_flags (SamFlags flags)
{
    struct SamFlagsBits f = flags.bits;
    DisFlagsStr s;
    sprintf (s.s, "multi_segments=%u is_aligned=%u unmapped=%u:%u revcomp=%u:%u mate#=%u:%u secondary=%u failfilt=%u dup=%u supp=%u",
             f.multi_segments, f.is_aligned, f.unmapped, f.next_unmapped, f.rev_comp, f.next_rev_comp, 
             f.is_first, f.is_last, f.secondary, f.filtered, f.duplicate, f.supplementary);
    return s;
}
