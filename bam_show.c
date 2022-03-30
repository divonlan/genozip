// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "genozip.h"
#include "buffer.h"
#include "vblock.h"
#include "sam_private.h"
#include "strings.h"

static void bam_show_one_aux (STRp(aux))
{
    ASSERT (aux_len > 3, "\nExpecting aux_len=%u > 3", aux_len);

    iprintf ("%c%c:%c:", aux[0], aux[1], aux[2]);

    rom next_field = &aux[3];
    switch (aux[2]) {
        case 'c' : iprintf ("%d ", NEXT_UINT8);   break;
        case 'C' : iprintf ("%u ", NEXT_UINT8);   break;
        case 's' : iprintf ("%d ", NEXT_UINT16);  break;
        case 'S' : iprintf ("%u ", NEXT_UINT16);  break;
        case 'i' : iprintf ("%d ", NEXT_UINT32);  break;
        case 'I' : iprintf ("%u ", NEXT_UINT32);  break;
        case 'f' : iprintf ("%f ", NEXT_FLOAT32); break;
        case 'Z' : iprintf ("%s ", next_field); next_field += strlen (next_field)+1 ; break;
        
        case 'H' : { 
            int len = strlen(next_field);
            char hex[len];
            iprintf ("%s ", str_to_hex ((bytes)next_field, len, hex)); next_field += len+1 ; 
            break;
        }
        
        case 'B' : {
            uint8_t type   = NEXT_UINT8;
            uint32_t count = NEXT_UINT32;
            
            for (uint32_t i=0; i < count; i++) 
                switch (type) {
                    #define SEP (i==count-1 ? ' ' : ',')
                    case 'c' : iprintf ("%d%c", NEXT_UINT8,   SEP); break;
                    case 'C' : iprintf ("%u%c", NEXT_UINT8,   SEP); break;
                    case 's' : iprintf ("%d%c", NEXT_UINT16,  SEP); break;
                    case 'S' : iprintf ("%u%c", NEXT_UINT16,  SEP); break;
                    case 'i' : iprintf ("%d%c", NEXT_UINT32,  SEP); break;
                    case 'I' : iprintf ("%u%c", NEXT_UINT32,  SEP); break;
                    case 'f' : iprintf ("%f%c", NEXT_FLOAT32, SEP); break;
                }
        }
    }
}

rom bam_show_line (VBlockSAMP vb, rom alignment, uint32_t remaining_txt_len)   
{

    rom next_field = alignment;
    uint32_t block_size = NEXT_UINT32;
    rom after = alignment + block_size + sizeof (uint32_t);

    if (segconf.running) return after;

    iprintf ("block_size=%d ", block_size);

    // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
    ASSERT (block_size + 4 >= sizeof (BAMAlignmentFixed) && block_size + 4 <= remaining_txt_len, 
            "vb=%u line_i=%"PRIu64" (block_size+4)=%u is out of range - too small, or goes beyond end of txt data: remaining_txt_len=%u",
            vb->vblock_i, vb->line_i, block_size+4, remaining_txt_len);


    iprintf ("rname=%d ", NEXT_UINT32);
    iprintf ("pos=%u ", 1 + NEXT_UINT32);
    uint8_t l_read_name  = NEXT_UINT8;               
    iprintf ("l_qname=%u ", l_read_name);    
    iprintf ("mapq=%u ", NEXT_UINT8);
    iprintf ("bin=%u ", NEXT_UINT16);
    uint16_t n_cigar_op  = NEXT_UINT16;
    iprintf ("n_cigar_op=%u ", n_cigar_op);    
    iprintf ("flag=%u ", NEXT_UINT16);
    uint32_t l_seq       = NEXT_UINT32;
    iprintf ("l_seq=%u ", l_seq);    
    iprintf ("rnext=%d ", NEXT_UINT32);
    iprintf ("pnext=%u ", 1+NEXT_UINT32);
    iprintf ("tlen=%u ", NEXT_UINT32);
    iprintf ("qname=%.*s cigar=", l_read_name-1, next_field); next_field += l_read_name; // note: l_read_name includes \0

    for (int i=0; i < n_cigar_op; i++) {
        BamCigarOp op = (BamCigarOp){ .value = NEXT_UINT32 };
        iprintf ("%u%c", op.n, cigar_op_to_char[op.op]);
    }
    
    rom textual_seq = bam_seq_display ((bytes)next_field, l_seq);
    iprintf (" seq=%s ", textual_seq); next_field += (l_seq+1)/2; FREE(textual_seq);
    
    rom textual_qual = bam_qual_display ((bytes)next_field, l_seq);
    iprintf ("qual=%s ", textual_qual); next_field += l_seq; FREE (textual_qual);

    // split auxillary fields
    rom auxs[MAX_FIELDS]; 
    uint32_t aux_lens[MAX_FIELDS];
    uint32_t n_auxs = bam_split_aux (vb, next_field, after, auxs, aux_lens);
    for (int i=0; i < n_auxs; i++)
        bam_show_one_aux (auxs[i], aux_lens[i]);

    iprint0 ("\n");

    return after;
}
