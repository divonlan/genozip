// ------------------------------------------------------------------
//   piz_fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "txtfile.h"
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

// called by I/O thread in piz_fast_read_one_vb, in case of --grep, to decompress and reconstruct the desc line, to 
// see if this vb is included. 
bool piz_fast_test_grep (VBlockFAST *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->lines.len        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    // we only need room for one line for now 
    buf_alloc (vb, &vb->txt_data, vb->longest_line_len, 1.1, "txt_data", vb->vblock_i);

    // uncompress & map desc field (filtered by zfile_is_skip_section)
    vb->grep_stages = GS_TEST; // tell zfile_is_skip_section to skip decompressing sections not needed for determining the grep
    piz_uncompress_all_ctxs ((VBlockP)vb);
    vb->grep_stages = GS_UNCOMPRESS; // during uncompress in the compute thread, uncompress only what was not already uncompressed here

    piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper);

    // reconstruct each description line and check for string matching with flag_grep
    bool found = false, match = false;

    MtfContext *desc_ctx = &vb->mtf_ctx[FAST_DESC];
    desc_ctx->iterator.next_b250 = FIRSTENT (uint8_t, desc_ctx->b250); 

    uint32_t txt_line_i = vb->data_type == DT_FASTQ ? 4 * vb->first_line : vb->first_line;

    while (desc_ctx->iterator.next_b250 < AFTERENT (uint8_t, desc_ctx->b250)) {
        DECLARE_SNIP;
        LOAD_SNIP (FAST_DESC);
        piz_reconstruct_compound_field ((VBlockP)vb, &vb->desc_mapper, 0, 0, snip, snip_len, txt_line_i);

        *AFTERENT (char, vb->txt_data) = 0; // terminate the desc string

        match = !!strstr (vb->txt_data.data, flag_grep);

        vb->txt_data.len = 0; // reset

        if (match) { // 
            found = true; // we've found a match to the grepped string
            if (vb->data_type == DT_FASTQ) break; // for FASTA, we need to go until the last line, for FASTQ, we can break here
        }

        if (vb->data_type == DT_FASTQ) txt_line_i += 4; // note: for FASTA we have no idea what txt line we're on, because we're only tracking DESC lines
    }

    // last FASTA - carry over whether its grepped to the next VB - in case next VB starts not from the description line
    // similarly, note whether the previous VB ended with a grepped sequence. If previous VB didn't have any description
    // i.e the entire VB was a sequence that started in an earlier VB - the grep status of the easier VB is carried forward
    if (vb->data_type == DT_FASTA) {
        static bool fasta_prev_vb_last_line_was_grepped = false;

        // if no match was found for this VB, but we have one carried over from previous VB then include this VB anyway
        if (!found) found = fasta_prev_vb_last_line_was_grepped;

        vb->fasta_prev_vb_last_line_was_grepped = fasta_prev_vb_last_line_was_grepped;
        if (desc_ctx->b250.len) fasta_prev_vb_last_line_was_grepped = match; // update to the last description, IF this VB contained any description
    }

    // reset iterators - piz_fastq_reconstruct_vb will use them again 
    mtf_init_iterator (desc_ctx);
    for (unsigned sf_i=0; sf_i < vb->desc_mapper.num_subfields; sf_i++)
        mtf_init_iterator (&vb->mtf_ctx[vb->desc_mapper.did_i[sf_i]]);

    return found; // no match found
}

void piz_fastq_reconstruct_vb (VBlockFAST *vb)
{
    if (!flag_grep) piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper); // it not already done during grep

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);
    
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_line_i = 4 * (vb->first_line + vb_line_i); // each vb line is a fastq record which is 4 txt lines
        
        uint32_t txt_data_start_line = vb->txt_data.len;

        // metadata looks like this - "X151YXX" - the 4 X,Y characters specify whether each row has a \r (Y=has)
        // and the number is the seq_len=qual_len
        DECLARE_SNIP;
        LOAD_SNIP (FAST_LINEMETA);
        const char *md = snip;

        // description line
        LOAD_SNIP (FAST_DESC);
        piz_reconstruct_compound_field ((VBlockP)vb, &vb->desc_mapper, eol[md[0]-'X'], eol_len[md[0]-'X'], 
                                        snip, snip_len, txt_line_i);
        *AFTERENT (char, vb->txt_data) = 0; // null-terminate for sec, in case of grep

        bool grepped_out = false;
        // case: we're grepping, and this line doesn't match
        if (flag_grep && !strstr (&vb->txt_data.data[txt_data_start_line], flag_grep)) { 
            vb->txt_data.len = txt_data_start_line; // rollback
            grepped_out = true;
        }

        if (flag_header_one) continue; // this is invoked by --header-only (re-written to flag_header_one in piz_read_global_area)

        // sequence line
        uint32_t seq_len = atoi (&md[4]); // numeric string terminated by dictionary's SNIP_SEP
        piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[FAST_SEQ], seq_len, txt_line_i, grepped_out);
        if (!grepped_out) RECONSTRUCT (eol[md[1]-'X'], eol_len[md[1]-'X']); // end of line

        // + line
        if (!grepped_out) RECONSTRUCT (md[2]-'X' ? "+\r\n" : "+\n", eol_len[md[2]-'X'] + 1);

        // quality line
        piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[FASTQ_QUAL], seq_len, txt_line_i, grepped_out);
        
        if (!grepped_out) RECONSTRUCT (eol[md[3]-'X'], eol_len[md[3]-'X']); // end of line
    }
}

void piz_fasta_reconstruct_vb (VBlockFAST *vb)
{
    if (!flag_grep) piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper); // it not already done during grep

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size, 1.1, "txt_data", vb->vblock_i);

    bool grepped_out = flag_grep && !vb->fasta_prev_vb_last_line_was_grepped; // if we're continuing the previous VB's sequence, we obey its grep status

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_line_i = vb->first_line + vb_line_i;

        // metadata looks like this - "X>" (desc line), "X;" (comment line) "X123" (sequence line)
        // X, Y characters specify whether each row has a \r (Y=has)
        DECLARE_SNIP;
        LOAD_SNIP (FAST_LINEMETA);
        const char *md = snip;
        bool has_13 = md[0] - 'X';

        uint32_t txt_data_start_line = vb->txt_data.len;

        switch (md[1]) {
            case '>': // description line 
                LOAD_SNIP (FAST_DESC);
                piz_reconstruct_compound_field ((VBlockP)vb, &vb->desc_mapper, eol[has_13], eol_len[has_13], 
                                                snip, snip_len, txt_line_i);
                *AFTERENT (char, vb->txt_data) = 0; // terminate the desc string - for strstr below

                vb->last_line = FASTA_DESC;

                // case: we're grepping, and this line doesn't match
                if (flag_grep && !strstr (&vb->txt_data.data[txt_data_start_line], flag_grep)) { 
                    vb->txt_data.len = txt_data_start_line; // rollback
                    grepped_out = true;
                }
                else grepped_out = false;
                break;

            case ';': // comment line
                if (!flag_header_one && !flag_grep)
                    piz_reconstruct_from_ctx ((VBlockP)vb, FASTA_COMMENT, eol[has_13], eol_len[has_13], txt_line_i); 
                    //RECONSTRUCT_FROM_BUF (vb->comment_data, vb->next_comment, "COMMENT", eol[has_13], eol_len[has_13]);

                //if (flag_header_one && !snip_len)
                //    vb->txt_data.len -= eol_len[has_13]; // don't show empty lines in --header-only mode

                vb->last_line = FASTA_COMMENT;
                break;

            default:  // sequence line
                if (!flag_header_one) { // this is invoked by --header-only (re-written to flag_header_one in piz_read_global_area)

                    if (flag_fasta_sequential && vb->last_line == FASTA_SEQ && vb->txt_data.len >= 2)
                        vb->txt_data.len -= 1 + (vb->txt_data.data[vb->txt_data.len-2]=='\r');

                    uint32_t seq_len = atoi (&md[1]); // numeric string terminated by dictionary's SNIP_SEP separator
                    piz_reconstruct_seq_qual ((VBlockP)vb, &vb->mtf_ctx[FASTA_SEQ], seq_len, txt_line_i, grepped_out);
                    if (!grepped_out) RECONSTRUCT (eol[has_13], eol_len[has_13]); // end of line
                    vb->last_line = FASTA_SEQ;
                }
        }
    }
}

bool piz_fast_read_one_vb (VBlock *vb, SectionListEntry *sl)
{ 
    // if we're grepping we we uncompress and reconstruct the DESC from the I/O thread, and terminate here if this VB is to be skipped
    if (flag_grep && !piz_fast_test_grep ((VBlockFASTP)vb)) return false; 

    return true;
}