// ------------------------------------------------------------------
//   seg_fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "header.h"
#include "file.h"
#include "strings.h"

#define DATA_LINE(vb,i) (&((ZipDataLineFAST *)((vb)->data_lines))[(i)])

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *seg_fastq_data_line (VBlock *vb_,   
                                 const char *field_start_line,     // index in vb->txt_data where this line starts
                                 uint32_t vb_line_i) // line within this vb (starting from 0)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;
    ZipDataLineFAST *dl = DATA_LINE (vb, vb_line_i);

    const char *next_field, *field_start=field_start_line;
    unsigned field_len=0;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n
    char metadata[30];   // the miscellanious non-description, non-seq parts of the line

    int32_t len = AFTERENT (char, vb->txt_data) - field_start_line;

    // the leading @ - just verify it (it will be included in D0ESC subfield)
    ASSERT (*field_start == '@', "%s: Invalid FASTQ file format in %s: expecting description line to start with @ but it starts with %c. vb_i=%u vb_line_i=%u",
            global_cmd, txt_name, *field_start, vb->vblock_i, vb_line_i);

    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG (<-- this is Illumina format)
    // See here for details of Illumina subfields: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    next_field = seg_get_next_item (field_start, &len, true, false, false, vb_line_i, &field_len, &separator, &has_13, "DESC");

    // we segment it using / | : and " " as separators. 
    seg_compound_field ((VBlockP)vb, &vb->mtf_ctx[FAST_DESC], field_start, field_len, &vb->desc_mapper,
                        dict_id_fast_desc_sf (dict_id_make ("D0ESC", 5)), ' ', 
                        SEC_FAST_DESC_B250, SEC_FAST_DESC_SF_B250, vb_line_i);
    metadata[0] = 'X' + has_13;
    vb->txt_section_bytes[SEC_FAST_DESC_B250] += has_13;

    // SEQ - just get the whole line
    const char *seq_start = next_field;
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, true, false, false, vb_line_i, &dl->seq_len, &separator, &has_13, "SEQ");
    vb->txt_section_bytes[SEC_SEQ_DATA] += dl->seq_len + 1 + has_13;
    
    metadata[1] = 'X' + has_13;

    // next line is expected to be a "+"
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, true, false, false, vb_line_i, &field_len, &separator, &has_13, "+");
    ASSERT (*field_start=='+' && field_len==1, "%s: Invalid FASTQ file format in %s: expecting middle line to be a \"+\" (with no spaces) but it is \"%.*s\". vb_i=%u vb_line_i=%u",
            global_cmd, txt_name, field_len, field_start, vb->vblock_i, vb_line_i);

    metadata[2] = 'X' + has_13;    
    vb->txt_section_bytes[SEC_FAST_LINEMETA_B250] += 2 + has_13; // + \n and sometimes \r

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = next_field - vb->txt_data.data;
    unsigned qual_len;
    field_start = next_field;
    next_field = seg_get_next_item (next_field, &len, true, false, false, vb_line_i, &qual_len, &separator, &has_13, "QUAL");
    vb->txt_section_bytes[SEC_QUAL_DATA] += dl->seq_len + 1 + has_13;
    metadata[3] = 'X' + has_13;    

    ASSERT (qual_len == dl->seq_len, "%s: Invalid FASTQ file format in %s: sequence_len=%u and quality_len=%u. Expecting them to be the same. vb_i=%u vb_line_i=%u\nSEQ=%.*s\nQUAL==%.*s",
            global_cmd, txt_name, dl->seq_len, qual_len, vb->vblock_i, vb_line_i, dl->seq_len, seq_start, qual_len, field_start);

    // METADATA - eg XYXX151 - the 4 X,Y characters specify whether each row has a \r (Y=has), and the number is the seq_len=qual_len
    
    // add the seq_len to the end of the metadata
    unsigned seq_len_str_len;
    str_uint (dl->seq_len, &metadata[4], &seq_len_str_len);

    seg_one_field (vb, metadata, 4 + seq_len_str_len, vb_line_i, FAST_LINEMETA);

    return next_field;
}

// Fasta format(s): https://en.wikipedia.org/wiki/FASTA_format
// concept: we segment each line separately, and for each line, we store an element in TEMPLATE about it. The
// Metadata elements are:
// > - description line - this (1) any line starting with > or (2) the first line starting with ; at the start 
//     of a file or after a sequence
//     the descrition line data is further segmented and stored in the DESC dictionary and D0SEC subfields
// ; - a comment line - any other line that starts with a ; or an empty line
//     the comment data (which can be empty for an empty line) is stored in a data buffer (not dictionary)
//     note: if a comment line is the first line in a VB - it will be segmented as a description. No harm done.
// 123 - a sequence line - any line that's not a description of sequence line - store its length
// these ^ are preceded by a 'Y' if the line has a Windows-style \r\n line ending or 'X' if not
const char *seg_fasta_data_line (VBlock *vb_,   
                                 const char *line_start, // index in vb->txt_data where this line starts
                                 uint32_t vb_line_i) // line within this vb (starting from 0)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;	

    // get entire line
    unsigned line_len;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n
    int32_t remaing_vb_txt_len = AFTERENT (char, vb->txt_data) - line_start;
    const char *next_field = seg_get_next_item (line_start, &remaing_vb_txt_len, true, false, false, vb_line_i, &line_len, &separator, &has_13, "DESC");
    
    // case: description line - we segment it to its components
    if (*line_start == '>' || (*line_start == ';' && vb->last_line == FASTA_LINE_SEQ)) {
        // we segment using / | : and " " as separators. 
        seg_compound_field ((VBlockP)vb, &vb->mtf_ctx[FAST_DESC], line_start, line_len, &vb->desc_mapper,
                            dict_id_fast_desc_sf (dict_id_make ("D0ESC", 5)), ' ', 
                            SEC_FAST_DESC_B250, SEC_FAST_DESC_SF_B250, vb_line_i);
        
        static const char *desc_md[2] = { "X>", "Y>" };
        seg_one_field (vb, desc_md[has_13], 2, vb_line_i, FAST_LINEMETA);
    }

    // case: comment line - stored in the comment buffer
    else if (*line_start == ';' || !line_len) {
        seg_add_to_data_buf (vb_, &vb->comment_data, SEC_FASTA_COMMENT_DATA, line_start, line_len, '\n', line_len + 1 + has_13); 

        static const char *comment_md[2] = { "X;", "Y;" };
        seg_one_field (vb, comment_md[has_13], 2, vb_line_i, FAST_LINEMETA);
    }

    // case: sequence line
    else {
        DATA_LINE (vb, vb_line_i)->seq_data_start = line_start - vb->txt_data.data;
        DATA_LINE (vb, vb_line_i)->seq_len        = line_len;
        vb->txt_section_bytes[SEC_SEQ_DATA] += line_len + 1 + has_13;

        char seq_md[30];
        unsigned seq_len_len;
        seq_md[0] = 'X' + has_13;
        str_uint (line_len, &seq_md[1], &seq_len_len);
        seg_one_field (vb, seq_md, seq_len_len + 1, vb_line_i, FAST_LINEMETA);
    }

    return next_field;
}
