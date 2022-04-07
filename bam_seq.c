// ------------------------------------------------------------------
//   sam_bam_seq.c - functions for handling BAM binary sequence format 
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"

// the characters "=ACMGRSVTWYHKDBN" are mapped to BAM 0->15, in this matrix we add 0x80 as a validity bit. All other characters are 0x00 - invalid
static const uint8_t sam2bam_seq_map[256] = { ['=']=0x80, ['A']=0x81, ['C']=0x82, ['M']=0x83, ['G']=0x84, ['R']=0x85, ['S']=0x86, ['V']=0x87, 
                                              ['T']=0x88, ['W']=0x89, ['Y']=0x8a, ['H']=0x8b, ['K']=0x8c, ['D']=0x8d, ['B']=0x8e, ['N']=0x8f };

const char bam_base_codes[16] = "=ACMGRSVTWYHKDBN";

rom bam_seq_display (bytes seq, uint32_t l_seq) // caller should free memory
{
    char *str = MALLOC (l_seq + 2);

    for (uint32_t i=0; i < (l_seq+1)/2; i++) {
        str[i*2]   = bam_base_codes[seq[i] >> 4];
        str[i*2+1] = bam_base_codes[seq[i] & 0xf];
    }

    str[l_seq] = 0;
    return str;
}

// called from sam_zip_prim_ingest_vb, somewhat similar to sam_piz_sam2bam_SEQ
void sam_seq_to_bam (STRp (seq_sam), Buffer *seq_bam_buf)
{
    uint8_t *seq_bam = BAFT8 (*seq_bam_buf);
    uint32_t seq_bam_len = (seq_sam_len+1)/2;

    for (uint32_t i=0; i < seq_bam_len; i++, seq_bam++, seq_sam += 2) {
        uint8_t base[2] = { sam2bam_seq_map[(uint8_t)seq_sam[0]], sam2bam_seq_map[(uint8_t)seq_sam[1]] };
        
        // check for invalid characters 
        for (unsigned b=0; b < 2; b++)
            if (!base[b] && !(b==1 && (i+1)*2 > seq_sam_len)) {
                char printable[MIN_(1000,seq_sam_len)+1]; // +1 for \0
                ASSINP (false, "Invalid base: invalid character encountered in sequence: '%c' (ASCII %u). position %u SEQ(first 1000 bases)=\"%s\"", 
                        base[b], base[b], i*2+b, str_to_printable (seq_sam, MIN_(1000,seq_sam_len), printable));
                base[b] = 0x0f;
            }

        *seq_bam = (base[0] << 4) | (base[1] & 0x0f);
    }

    // if number of bases is odd, zero the last, unused, base
    if (seq_bam_len & 1) 
        seq_bam[seq_bam_len-1] &= 0xf0;

    seq_bam_buf->len += seq_bam_len;
}

// re-writes BAM format SEQ into textual SEQ
void bam_seq_to_sam (VBlockSAMP vb, bytes bam_seq, uint32_t seq_len /*bases, not bytes*/, 
                     bool start_mid_byte,    // ignore first nibble of bam_seq (seq_len doesn't include the ignored nibble)
                     bool test_final_nibble, // if true, we test that the final nibble, if unused, is 0, and warn if not
                     Buffer *out)            // appends to end of buffer - caller should allocate
{
    if (!seq_len) {
        BNXTc (*out) = '*';
        return;        
    }

    // we implement "start_mid_byte" by converting the redudant base too, but starting 1 character before in the buffer 
    char save = 0;
    if (start_mid_byte) {
        out->len--;
        seq_len++;
        save = *BAFTc (*out); // this is the byte we will overwrite, and recover it later. possibly, the fence if the buffer is empty;
    }
    
    char *next = BAFTc (*out);

    for (uint32_t i=0; i < (seq_len+1) / 2; i++) {
        *next++ = bam_base_codes[*(uint8_t*)bam_seq >> 4];
        *next++ = bam_base_codes[*(uint8_t*)bam_seq & 0xf];
        bam_seq++;
    }

    if (start_mid_byte) {
        *BAFTc(*out) = save;
        out->len += seq_len;      
        seq_len--;
    }
    else
        out->len += seq_len;

    ASSERTW (!test_final_nibble || !(seq_len % 2) || (*BAFTc (*out)=='='), 
             "Warning in bam_seq_to_sam vb=%u: expecting the unused lower 4 bits of last seq byte in an odd-length seq_len=%u to be 0, but its not. This will cause an incorrect MD5",
             vb->vblock_i, seq_len);
}

void bam_seq_copy (VBlockSAMP vb, bytes bam_seq, uint32_t seq_len /*bases, not bytes*/, 
                   bool start_mid_byte, // ignore first nibble of bam_seq (seq_len doesn't include the ignored nibble)
                   Buffer *out_buf)     // appends to end of buffer (byte-aligned) - caller should allocate
{
    uint8_t *out = BAFT(uint8_t, *out_buf);

    if (start_mid_byte) {
        for (uint32_t i=0; i < seq_len / 2; i++) 
            *out++ = ((bam_seq[i] & 0x0f) << 4) | (bam_seq[i+1] >> 4); 
        
        if (seq_len & 1)
            *out++ = (bam_seq[seq_len / 2] & 0x0f) << 4;
    }
    else {
        memcpy (out, bam_seq, (seq_len+1)/2);

        if (vb->seq_len & 1)
            out[seq_len/2] &= 0xf0; // keep the high nibble only (the first BAM base of this byte)
    }

    out_buf->len += (seq_len+1)/2;
}

/*
// compare a sub-sequence to a full sequence and return true if they're the same. Sequeneces in BAM format.
bool bam_seq_has_sub_seq (bytes full_seq, uint32_t full_seq_len, 
                          bytes sub_seq,  uint32_t sub_seq_len, uint32_t start_base) // lengths are in bases, not bytes
{
    // easy case: similar byte alignment
    if (!(start_base & 1)) {
        if (memcmp (sub_seq, &full_seq[start_base/2], sub_seq_len/2)) return false; // mismatch in first even number of bases
        if (!(sub_seq_len & 1)) return true; // even number of bases

        uint8_t last_base_full = full_seq[(start_base + sub_seq_len - 1)/2] >> 4;
        uint8_t last_base_sub  = sub_seq[(sub_seq_len - 1)/2] >> 4;
        return last_base_full == last_base_sub;
    }

    // not byte-aligned
    else { 
        for (uint32_t sub_base_i=0; sub_base_i < sub_seq_len; sub_base_i++) {
    
            uint8_t base_sub = (sub_base_i % 2) ? (sub_seq[sub_base_i/2] & 15) : (sub_seq[sub_base_i/2] >> 4);
                        
            uint32_t full_base_i = start_base + sub_base_i;
            uint8_t base_full = (full_base_i % 2) ? (full_seq[full_base_i/2] & 15) : (full_seq[full_base_i/2] >> 4);

            if (base_sub != base_full) return false;
        }
        return true; // all bases are the same
    }
}

// C<>G A<>T ; IUPACs: R<>Y K<>M B<>V D<>H W<>W S<>S N<>N
static void bam_seq_revcomp_in_place (uint8_t *seq, uint32_t seq_len)
{                                 // Was:  =    A    C    M    G    R    S    V    T    W    Y    H    K    D    B    N                          
                                  // Comp: =    T    G    K    C    Y    S    B    A    W    R    D    M    H    V    N
    static const uint8_t bam_comp[16] = { 0x0, 0x8, 0x4, 0xc, 0x2, 0xa, 0x6, 0xe, 0x1, 0x9, 0x5, 0xd, 0x3, 0xb, 0x7, 0xf };

    for (int32_t i=0, j=seq_len-1; i < seq_len/2; i++, j--) {
        uint8_t b1 = (i&1) ? (seq[i/2] & 0xf) : (seq[i/2] >> 4);
        uint8_t b1c = bam_comp[b1];

        uint8_t b2 = (j&1) ? (seq[j/2] & 0xf) : (seq[j/2] >> 4);
        uint8_t b2c = bam_comp[b2];

        seq[i/2] = (i&1) ? (b2c | (seq[i/2] & 0xf0)) : ((b2c << 4) | (seq[i/2] & 0xf));
        seq[j/2] = (j&1) ? (b1c | (seq[j/2] & 0xf0)) : ((b1c << 4) | (seq[j/2] & 0xf));
    }
}

// like bam_seq_has_sub_seq, but sub_seq is reverse complemented, and start_base is relative to END of full_seq
bool bam_seq_has_sub_seq_revcomp (bytes full_seq, uint32_t full_seq_len, 
                                  bytes sub_seq,  uint32_t sub_seq_len, uint32_t start_base)
{
    bam_seq_revcomp_in_place ((uint8_t*)STRa(sub_seq));

    bool same = bam_seq_has_sub_seq (STRa(full_seq), STRa(sub_seq), full_seq_len - start_base - sub_seq_len);

    // revcomp-back, if we have to
    bam_seq_revcomp_in_place ((uint8_t*)STRa(sub_seq));

    return same;
}
*/