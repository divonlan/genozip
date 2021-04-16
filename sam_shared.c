// ------------------------------------------------------------------
//   sam_shared.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "sam_private.h"
#include "strings.h"
#include "file.h"

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

unsigned sam_vb_size (void) { return sizeof (VBlockSAM); }
unsigned sam_vb_zip_dl_size (void) { return sizeof (ZipDataLineSAM); }

void sam_vb_release_vb (VBlockSAM *vb)
{
    vb->last_cigar = NULL;
    vb->ref_consumed = vb->ref_and_seq_consumed = 0;
    vb->a_bases = vb->x_bases = vb->y_bases = 0;
    vb->a_index = vb->x_index = vb->y_index = 0;
    vb->soft_clip = 0;
    buf_free (&vb->bd_bi_line);
    buf_free (&vb->textual_cigar);
    buf_free (&vb->textual_seq);
    buf_free (&vb->textual_opt);
}

void sam_vb_destroy_vb (VBlockSAM *vb)
{
    buf_destroy (&vb->bd_bi_line);
    buf_destroy (&vb->textual_cigar);
    buf_destroy (&vb->textual_seq);
    buf_destroy (&vb->textual_opt);
}

// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
void sam_analyze_cigar (VBlockSAMP vb, const char *cigar, unsigned cigar_len, 
                        unsigned *seq_consumed, unsigned *ref_consumed, unsigned *seq_and_ref, unsigned *soft_clip) // optional outs
{
    if (seq_consumed) *seq_consumed = 0;
    if (ref_consumed) *ref_consumed = 0;
    if (seq_and_ref)  *seq_and_ref  = 0;
    if (soft_clip)    *soft_clip    = 0;

    ASSERT (cigar[0] != '*' || cigar_len == 1, "Invalid CIGAR: %.*s", cigar_len, cigar); // a CIGAR start with '*' must have 1 character

    // ZIP case: if the CIGAR is "*", later sam_seg_cigar_field uses the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*", so this condition is always false
    if (cigar[0] == '*') return;

    // PIZ case: CIGAR string starts with '-' (indicating missing SEQ) - just skip the '-' for now
    if (*cigar == '-') {
        cigar++;
        cigar_len--;
    }

    // if we're reconstructing a BAM, we will create the BAM cigar data in textual_cigar. 
    bool bam_piz = (command == PIZ && flag.out_dt == DT_BAM);
    if (bam_piz) buf_alloc_old (vb, &vb->textual_cigar, cigar_len/2 /* max possible n_cigar_op */ * sizeof(uint32_t), 2, "textual_cigar");

    unsigned n=0;
    for (unsigned i=0; i < cigar_len; i++) {

        char c = cigar[i];
        char lookup = cigar_lookup_sam[(uint8_t)c];

        ASSINP (lookup, "Invalid CIGAR in %s: invalid operation %c. CIGAR=%.*s", txt_name, cigar[i], cigar_len, cigar);
        lookup &= 0x0f; // remove validity bit

        if (lookup == CIGAR_DIGIT) 
            n = n*10 + (c - '0');
        
        else {
            ASSINP (n, "Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            
            if ((lookup & CIGAR_CONSUMES_QUERY)     && seq_consumed) *seq_consumed += n;
            if ((lookup & CIGAR_CONSUMES_REFERENCE) && ref_consumed) *ref_consumed += n;
            if ((lookup & CIGAR_CONSUMES_QUERY) && (lookup & CIGAR_CONSUMES_REFERENCE) && seq_and_ref) *seq_and_ref += n;
            if (soft_clip && c == 'S') *soft_clip += n;

            // note: piz: in case of eg "151*" - *seq_consumed will be updated to the length, but textual_cigar will be empty
            if (bam_piz && c != '*') { 
                // convert character CIGAR op to BAM cigar field op  "MIDNSHP=X" -> 012345678 
                static const uint8_t cigar_char_to_op[256] = { ['M']=0, ['I']=1, ['D']=2, ['N']=3, ['S']=4, 
                                                               ['H']=5, ['P']=6, ['=']=7, ['X']=8           }; 
                uint32_t bam_cigar = (n << 4) | (uint32_t)cigar_char_to_op[(uint8_t)c];
                NEXTENT (uint32_t, vb->textual_cigar) = LTEN32 (bam_cigar);
            }
            n = 0;
        }
    }                          

    ASSINP (!n, "Invalid CIGAR in %s: expecting it to end with an operation character. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    ASSINP (!seq_consumed || *seq_consumed, "Invalid CIGAR in %s: CIGAR implies 0-length SEQ. CIGAR=%.*s", txt_name, cigar_len, cigar);
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
    // to reconstruct (in sam_piz_special_CIGAR)

    if (first_pos_0>>14 == last_pos_0>>14) return ((1<<15)-1)/7 + (first_pos_0>>14);
    if (first_pos_0>>17 == last_pos_0>>17) return ((1<<12)-1)/7 + (first_pos_0>>17);
    if (first_pos_0>>20 == last_pos_0>>20) return ((1<<9 )-1)/7 + (first_pos_0>>20);
    if (first_pos_0>>23 == last_pos_0>>23) return ((1<<6 )-1)/7 + (first_pos_0>>23);
    if (first_pos_0>>26 == last_pos_0>>26) return ((1<<3 )-1)/7 + (first_pos_0>>26);
    return 0;
}
