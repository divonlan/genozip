// ------------------------------------------------------------------
//   piz_fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "txtfile.h"
#include "header.h"
#include "vblock.h"
#include "base250.h"
#include "move_to_front.h"
#include "file.h"
#include "endianness.h"
#include "piz.h"
#include "sections.h"
#include "dict_id.h"

static const char *eol[2] = { "\n", "\r\n"};
static const unsigned eol_len[2] = { 1, 2 };

static void piz_fastq_reconstruct_vb (VBlockFAST *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        uint32_t snip_len;
        const char *snip;

        uint32_t txt_line_i = vb->first_line + vb_line_i;

        // metadata looks like this - "X151YXX" - the 4 X,Y characters specify whether each row has a \r (Y=has)
        // and the number is the seq_len=qual_len
        RECONSTRUCT_FROM_DICT (FAST_LINEMETA);
        const char *md = snip;

        // description line
        LOAD_SNIP (FAST_DESC);
        piz_reconstruct_compound_field ((VBlockP)vb, &vb->desc_mapper, eol[md[0]-'X'], eol_len[md[0]-'X'], 
                                        snip, snip_len, txt_line_i);

        // sequence line
        uint32_t seq_len = atoi (&md[4]); // numeric string terminated by dictionary's \t separator
        piz_reconstruct_seq_qual ((VBlockP)vb, seq_len, &vb->seq_data, &vb->next_seq, SEC_SEQ_DATA, txt_line_i);
        buf_add (&vb->txt_data, eol[md[1]-'X'], eol_len[md[1]-'X']); // end of line

        // + line
        buf_add (&vb->txt_data, md[2]-'X' ? "+\n" : "+\r\n", eol_len[md[2]-'X'] + 1);

        // quality line
        piz_reconstruct_seq_qual ((VBlockP)vb, seq_len, &vb->qual_data, &vb->next_qual, SEC_QUAL_DATA, txt_line_i);
        buf_add (&vb->txt_data, eol[md[3]-'X'], eol_len[md[3]-'X']); // end of line
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

static void piz_fasta_reconstruct_vb (VBlockFAST *vb)
{
    START_TIMER;

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        uint32_t snip_len;
        const char *snip;

        uint32_t txt_line_i = vb->first_line + vb_line_i;

        // metadata looks like this - "X>" (desc line), "X;" (comment line) "X123" (sequence line)
        // X, Y characters specify whether each row has a \r (Y=has)
        RECONSTRUCT_FROM_DICT (FAST_LINEMETA);
        const char *md = snip;
        bool has_13 = md[0] - 'X';

        switch (md[1]) {
            case '>': // description line 
                LOAD_SNIP (FAST_DESC);
                piz_reconstruct_compound_field ((VBlockP)vb, &vb->desc_mapper, eol[has_13], eol_len[has_13], 
                                                snip, snip_len, txt_line_i);
                break;

            case ';': // comment line
                RECONSTRUCT_FROM_BUF (vb->comment_data, vb->next_comment, "COMMENT", '\n', eol[has_13], eol_len[has_13]);
                break;

            default: { // sequence line
                uint32_t seq_len = atoi (&md[1]); // numeric string terminated by dictionary's \t separator
                piz_reconstruct_seq_qual ((VBlockP)vb, seq_len, &vb->seq_data, &vb->next_seq, SEC_SEQ_DATA, txt_line_i);
                buf_add (&vb->txt_data, eol[has_13], eol_len[has_13]); // end of line
            }
        }
    }

    COPY_TIMER(vb->profile.piz_reconstruct_vb);
}

static void piz_fast_uncompress_all_sections (VBlockFAST *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->num_lines        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    unsigned section_i=1;

    // uncompress the fields     
    piz_uncompress_fields ((VBlockP)vb, section_index, &section_i);
    
    // DESC subfields
    piz_uncompress_compound_field ((VBlockP)vb, SEC_FAST_DESC_B250, SEC_FAST_DESC_SF_B250, &vb->desc_mapper, &section_i);

    // SEQ    
    SectionHeader *seq_header = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
    zfile_uncompress_section ((VBlockP)vb, seq_header, &vb->seq_data, "seq_data", SEC_SEQ_DATA);    

    // QUAL (FASTQ only)
    if (z_file->data_type == DT_FASTQ) {
        SectionHeader *qual_header = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
        if (!flag_strip) zfile_uncompress_section ((VBlockP)vb, qual_header, &vb->qual_data, "qual_data", SEC_QUAL_DATA);    
    }

    // COMMENT (FASTA only)
    else {
        SectionHeader *comment_header = (SectionHeader *)(vb->z_data.data + section_index[section_i++]);
        if (!flag_strip) zfile_uncompress_section ((VBlockP)vb, comment_header, &vb->comment_data, "comment_data", SEC_FASTA_COMMENT_DATA);    
    }
}

void piz_fast_uncompress_one_vb (VBlock *vb_)
{
    START_TIMER;

    VBlockFAST *vb = (VBlockFAST *)vb_;

    piz_fast_uncompress_all_sections ((VBlockFASTP)vb);

    if (z_file->data_type == DT_FASTQ)
        piz_fastq_reconstruct_vb ((VBlockFASTP)vb);
    else
        piz_fasta_reconstruct_vb ((VBlockFASTP)vb);

    vb->is_processed = true; // tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway
    
    COPY_TIMER (vb->profile.compute);
}
