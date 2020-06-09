// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fast_private.h"
#include "seg.h"
#include "move_to_front.h"
#include "file.h"
#include "strings.h"
#include "piz.h"
#include "optimize.h"
#include "dict_id.h"

void fastq_seg_initialize (VBlockFAST *vb)
{
    // thread safety: this will be initialized by vb_i=1, while it holds a mutex in zip_compress_one_vb
    static bool structured_initialized = false;
    if (!structured_initialized) {
        seg_initialize_compound_structured ((VBlockP)vb, "D?ESC", &structured_DESC); 
        structured_initialized = true;
    }

    vb->contexts[FASTQ_SEQ].flags  = CTX_FL_LOCAL_LZMA;
    vb->contexts[FASTQ_SEQ].ltype  = CTX_LT_SEQUENCE;
    vb->contexts[FASTQ_QUAL].ltype = CTX_LT_SEQUENCE;
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *fastq_seg_txt_line (VBlockFAST *vb, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    ZipDataLineFAST *dl = DATA_LINE (vb->line_i);

    const char *next_field, *field_start=field_start_line;
    unsigned field_len=0;
    char separator;
/*
    MiniStructured st = { .repeats=1, .num_items=1, .flags=0, .repsep={0,0},
                          .items = { { .dict_id = dict_id_FASTQ_DESC, .did_i = DID_I_NONE, .seperator={0, '\n'} },
                                     { .dict_id = dict_id_FASTQ_SEQ,  .did_i = DID_I_NONE, .seperator={0, '\n'} },
                                     { .dict_id = dict_id_FASTQ_PLUS, .did_i = DID_I_NONE, .seperator={0, '\n'} },
                                     { .dict_id = dict_id_FASTQ_QUAL, .did_i = DID_I_NONE, .seperator={0, 'n'} } } };
*/
    int32_t len = (int32_t)(AFTERENT (char, vb->txt_data) - field_start_line);

    // the leading @ - just verify it (it will be included in D0ESC subfield)
    ASSSEG (*field_start == '@', field_start, "%s: Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            global_cmd, *field_start);

    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG (<-- this is Illumina format)
    // See here for details of Illumina subfields: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    next_field = seg_get_next_line (vb, field_start, &len, &field_len, has_13, "DESC");
 
    // we segment it using / | : and " " as separators. 
    seg_compound_field ((VBlockP)vb, &vb->contexts[FASTQ_DESC], field_start, field_len, &vb->desc_mapper, structured_DESC, true, 0);
    SEG_EOL (FASTQ_E1L, true);

    // SEQ - just get the whole line
    const char *seq_start = next_field;
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (vb, next_field, &len, true, false, false, &dl->seq_len, &separator, has_13, "SEQ");
    vb->contexts[FASTQ_SEQ].local.len += dl->seq_len;
    
    // Add LOOKUP snip with seq_len
    char snip[10];
    snip[0] = SNIP_LOOKUP;
    unsigned seq_len_str_len = str_int (dl->seq_len, &snip[1]);
    seg_by_did_i (vb, snip, 1 + seq_len_str_len, FASTQ_SEQ, dl->seq_len); 

    SEG_EOL (FASTQ_E1L, true);

    // PLUS - next line is expected to be a "+"
    GET_LAST_ITEM ("+");
    ASSSEG (*field_start=='+' && field_len==1, field_start, "%s: Invalid FASTQ file format: expecting middle line to be a \"+\" (with no spaces) but it is \"%.*s\"",
            global_cmd, field_len, field_start);

    seg_by_did_i (vb, "+", 1, FASTQ_PLUS, 1);
    SEG_EOL (FASTQ_E1L, true);

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = next_field - vb->txt_data.data;
    GET_LAST_ITEM ("QUAL");
    vb->contexts[FASTQ_QUAL].local.len += dl->seq_len;
    vb->contexts[FASTQ_QUAL].txt_len   += dl->seq_len;

    // End Of Line    
    SEG_EOL (FASTQ_E1L, true);

    ASSSEG (field_len == dl->seq_len, field_start, "%s: Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ=%.*s\nQUAL==%.*s",
            global_cmd, dl->seq_len, field_len, dl->seq_len, seq_start, field_len, field_start);
 
    return next_field;
}

// callback function for compress to get data of one line (called by comp_compress_bzlib)
void fastq_zip_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
                                          char **line_qual_data, uint32_t *line_qual_len, // out
                                          char **unused_data,   uint32_t *unused_len) 
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);
     
    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);
    *line_qual_len  = dl->seq_len;
    *unused_data    = NULL;
    *unused_len     = 0;

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag_optimize_QUAL) optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

// returns true if section is to be skipped reading / uncompressing
bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT sections

    // note that piz_read_global_area rewrites --header-only as flag_header_one
    if (flag_header_one && 
        (dict_id.num == dict_id_fields[FASTQ_SEQ] || dict_id.num == dict_id_fields[FASTQ_QUAL] || dict_id.num == dict_id_fields[FASTQ_PLUS]))
        return true;
        
    // when grepping by I/O thread - skipping all sections but DESC
    if (vb && flag_grep && (vb->grep_stages == GS_TEST) && 
        dict_id.num != dict_id_fields[FASTQ_DESC] && !dict_id_is_fast_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if (vb && flag_grep && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == dict_id_fields[FASTQ_DESC] || dict_id_is_fast_desc_sf (dict_id)))
        return true;

    return false;
}

void fastq_piz_reconstruct_vb (VBlockFAST *vb)
{
    if (!flag_grep) piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper); // it not already done during grep

    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        vb->line_i = 4 * (vb->first_line + vb_line_i); // each vb line is a fastq record which is 4 txt lines
        vb->dont_show_curr_line = false; // might become true due --regions or --grep
        
        uint32_t txt_data_start_line = vb->txt_data.len;

        piz_reconstruct_from_ctx (vb, FASTQ_DESC, 0);
        piz_reconstruct_from_ctx (vb, FASTQ_E1L,  0);

        // minor bug here: since all the EOL fields are aliases, in case header_one, we dont consume the
        // missing line EOLs and they will show wrong ones to the remaining lines. minor bc in virtually all files,
        // the EOL is identical in all lines
        if (!flag_header_one) { // note that piz_read_global_area rewrites --header-only as flag_header_one
            piz_reconstruct_from_ctx (vb, FASTQ_SEQ,  0);
            piz_reconstruct_from_ctx (vb, FASTQ_E2L,  0);
            piz_reconstruct_from_ctx (vb, FASTQ_PLUS, 0);
            piz_reconstruct_from_ctx (vb, FASTQ_E3L,  0);
            piz_reconstruct_from_ctx (vb, FASTQ_QUAL, 0);
            piz_reconstruct_from_ctx (vb, FASTQ_E4L,  0);
        }

        // case: we're grepping, and this line doesn't match
        *AFTERENT (char, vb->txt_data) = 0; // for strstr
        if (vb->dont_show_curr_line || (flag_grep && !strstr (ENT (char, vb->txt_data, txt_data_start_line), flag_grep)))
            vb->txt_data.len = txt_data_start_line; // rollback
    }
}
