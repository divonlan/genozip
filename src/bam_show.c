// ------------------------------------------------------------------
//   sam_bam.c
//   Copyright (C) 2020-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"

static void bam_show_one_aux (STRp(aux))
{
    ASSERT (aux_len > 3, "\nExpecting aux_len=%u > 3", aux_len);

    iprintf ("%c%c:%c:", aux[0], aux[1], aux[2]);

    rom next_field = &aux[3];
    switch (aux[2]) {
        case 'c' : iprintf ("%d ", (int8_t) NEXT_UINT8);   break;
        case 'C' : iprintf ("%u ",          NEXT_UINT8);   break;
        case 'A' : iprintf ("%c ",          NEXT_UINT8);   break;
        case 's' : iprintf ("%d ", (int16_t)NEXT_UINT16);  break;
        case 'S' : iprintf ("%u ",          NEXT_UINT16);  break;
        case 'i' : iprintf ("%d ", (int32_t)NEXT_UINT32);  break;
        case 'I' : iprintf ("%u ",          NEXT_UINT32);  break;
        case 'f' : iprintf ("%f ",          NEXT_FLOAT32); break;
        case 'H' :
        case 'Z' : iprintf ("%s ", next_field); next_field += strlen (next_field)+1 ; break;
        
        case 'B' : {
            uint8_t type   = NEXT_UINT8;
            uint32_t count = NEXT_UINT32;
            
            for (uint32_t i=0; i < count; i++) 
                switch (type) {
                    #define SEP (i==count-1 ? ' ' : ',')
                    case 'c' : iprintf ("%d%c", (int8_t) NEXT_UINT8,   SEP); break;
                    case 'C' : iprintf ("%u%c",          NEXT_UINT8,   SEP); break;
                    case 's' : iprintf ("%d%c", (int16_t)NEXT_UINT16,  SEP); break;
                    case 'S' : iprintf ("%u%c",          NEXT_UINT16,  SEP); break;
                    case 'i' : iprintf ("%d%c", (int32_t)NEXT_UINT32,  SEP); break;
                    case 'I' : iprintf ("%u%c",          NEXT_UINT32,  SEP); break;
                    case 'f' : iprintf ("%f%c",          NEXT_FLOAT32, SEP); break;
                }
        }
    }
}

rom bam_show_line (VBlockSAMP vb, rom alignment, uint32_t remaining_txt_len)   
{
    rom next_field = alignment;
    uint32_t block_size = NEXT_UINT32;
    rom after = alignment + block_size + sizeof (uint32_t);

    if (segconf_running) return after;

    iprintf ("block_size=%d ", block_size);

    // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
    ASSERT (block_size + 4 >= sizeof (BAMAlignmentFixed) && block_size + 4 <= remaining_txt_len, 
            "%s: (block_size+4)=%u is out of range - too small, or goes beyond end of txt data: remaining_txt_len=%u",
            LN_NAME, block_size+4, remaining_txt_len);


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

    iprintf ("qname=\"%s\" cigar=\"", str_to_printable_(next_field, MIN_(l_read_name-1,512)).s); // restrict to 512 in case l_read_name is corrupted
    next_field += l_read_name; // note: l_read_name includes \0

    for (int i=0; i < n_cigar_op; i++) {
        BamCigarOp op = (BamCigarOp){ .value = NEXT_UINT32 };
        iprintf ("%u%c", op.n, cigar_op_to_char[op.op]);
    }
    
    rom textual_seq = bam_seq_display ((bytes)next_field, l_seq);
    iprintf ("\" seq=\"%s\" ", textual_seq); next_field += (l_seq+1)/2; FREE(textual_seq);
    
    rom textual_qual = bam_qual_display ((bytes)next_field, l_seq);
    iprintf ("qual=\"%s\" ", textual_qual); next_field += l_seq; FREE (textual_qual);

    // split auxillary fields
    STR_ARRAY (aux, MAX_FIELDS) = bam_split_aux (vb, alignment, next_field, after, auxs, aux_lens);
    for (int i=0; i < n_auxs; i++)
        bam_show_one_aux (auxs[i], aux_lens[i]);

    iprint0 ("\n");

    return after;
}

// similar as bam_show_line, but for ASSSEG
rom bam_assseg_line (VBlockP vb)
{
    static Buffer show_buf = { .name = "bam_show_buf" };

    DO_ONCE { // if multiple concurrent threads - one prints and the other threads stall
        if (vb->line_start >= Ltxt) return "Invalid line_start";

        rom alignment = Btxt (vb->line_start);
        uint32_t remaining_txt_len = Ltxt - vb->line_start;

        rom next_field = alignment;
        uint32_t block_size = NEXT_UINT32;
        rom after = alignment + block_size + sizeof (uint32_t);

        bufprintf (vb, &show_buf, "block_size=%d ", block_size);

        // a non-sensical block_size might indicate an false-positive identification of a BAM alignment in bam_unconsumed
        if (block_size + 4 < sizeof (BAMAlignmentFixed) || block_size + 4 > remaining_txt_len)
            return "Cannot show Alignment - it has an invalid block_size"; 

        bufprintf (vb, &show_buf, "rname=%d ", NEXT_UINT32);
        bufprintf (vb, &show_buf, "pos=%u ", 1 + NEXT_UINT32);
        uint8_t l_read_name  = NEXT_UINT8;               
        bufprintf (vb, &show_buf, "l_qname=%u ", l_read_name);    
        bufprintf (vb, &show_buf, "mapq=%u ", NEXT_UINT8);
        bufprintf (vb, &show_buf, "bin=%u ", NEXT_UINT16);
        uint16_t n_cigar_op  = NEXT_UINT16;
        bufprintf (vb, &show_buf, "n_cigar_op=%u ", n_cigar_op);    
        bufprintf (vb, &show_buf, "flag=%u ", NEXT_UINT16);
        uint32_t l_seq       = NEXT_UINT32;
        bufprintf (vb, &show_buf, "l_seq=%u ", l_seq);    
        bufprintf (vb, &show_buf, "rnext=%d ", NEXT_UINT32);
        bufprintf (vb, &show_buf, "pnext=%u ", 1+NEXT_UINT32);
        bufprintf (vb, &show_buf, "tlen=%u ", NEXT_UINT32);

        bufprintf (vb, &show_buf, "qname=\"%s\" cigar=\"", str_to_printable_(next_field, MIN_(l_read_name-1,512)).s); // restrict to 512 in case l_read_name is corrupted
        next_field += l_read_name; // note: l_read_name includes \0

        for (int i=0; i < n_cigar_op; i++) {
            BamCigarOp op = (BamCigarOp){ .value = NEXT_UINT32 };
            bufprintf (vb, &show_buf, "%u%c", op.n, cigar_op_to_char[op.op]);
        }

        // use bufprint0, not bufprintf, for seq and qual as the latter is limited in length
        bufprint0 (vb, &show_buf, "\" seq=");
        rom textual_seq = bam_seq_display ((bytes)next_field, l_seq);
        bufprint0 (vb, &show_buf, textual_seq); 
        next_field += (l_seq+1)/2; 
        FREE (textual_seq);
        
        bufprint0 (vb, &show_buf, "\" qual=\""); 
        rom textual_qual = bam_qual_display ((bytes)next_field, l_seq);
        bufprint0 (vb, &show_buf, textual_qual);
        bufprint0 (vb, &show_buf, "\"");
        next_field += l_seq; 
        FREE (textual_qual);

        // split auxillary fields
        STR_ARRAY (aux, MAX_FIELDS) = bam_split_aux (VB_SAM, alignment, next_field, after, auxs, aux_lens);
        for (int i=0; i < n_auxs; i++)
            bam_show_one_aux (auxs[i], aux_lens[i]);

        return B1STc(show_buf);
    }
    else {
        while (1); // stall 
        return "";
    }
}
