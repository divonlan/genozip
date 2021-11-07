// ------------------------------------------------------------------
//   sam_shared.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
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

void sam_vb_release_vb (VBlockSAM *vb)
{
    vb->last_cigar = NULL;
    vb->ref_consumed = vb->ref_and_seq_consumed = vb->soft_clip = vb->mismatch_bases = 0;
    vb->a_bases = vb->x_bases = vb->y_bases = 0;
    vb->a_index = vb->x_index = vb->y_index = 0;
    vb->md_verified = 0;
    
    buf_free (&vb->bd_bi_line);
    buf_free (&vb->textual_cigar);
    buf_free (&vb->binary_cigar);
    buf_free (&vb->textual_seq);
    buf_free (&vb->md_M_is_ref);
    buf_free (&vb->qname_hash);
    buf_free (&vb->buddy_textual_cigars);
}

void sam_vb_destroy_vb (VBlockSAM *vb)
{
    buf_destroy (&vb->bd_bi_line);
    buf_destroy (&vb->textual_cigar);
    buf_destroy (&vb->textual_seq);
    buf_destroy (&vb->md_M_is_ref);
    buf_destroy (&vb->qname_hash);
    buf_destroy (&vb->buddy_textual_cigars);
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
