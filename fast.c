// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "file.h"
#include "strings.h"
#include "endianness.h"
#include "piz.h"
#include "dict_id.h"
#include "sections.h"
#include "zip.h"
#include "optimize.h"

typedef struct {
    uint32_t seq_data_start, qual_data_start; // start within vb->txt_data
    uint32_t seq_len;                         // length within vb->txt_data (in case of FASTQ, this length is also applies to quality, per FASTQ spec)
} ZipDataLineFAST;

// IMPORTANT: if changing fields in VBlockFAST, also update vb_fastq_release_vb and vb_fastq_destroy_vb
typedef struct VBlockFAST { // for FASTA and FASTQ
    VBLOCK_COMMON_FIELDS
    SubfieldMapper desc_mapper; // FASTA and FASTQ - ZIP & PIZ
    bool fasta_prev_vb_last_line_was_grepped; // FASTA-only -PIZ: whether previous VB's last line was grepped successfully. 

    // note: last_line is initialized to FASTA_LINE_SEQ (=0) so that a ; line as the first line of the VB is interpreted as a description, not a comment
    enum { FASTA_LINE_SEQ, FASTA_LINE_DESC, FASTA_LINE_COMMENT } last_line; // FASTA ZIP
} VBlockFAST;

#define DATA_LINE(i) ENT (ZipDataLineFAST, vb->lines, i)

static Structured structured_DESC;

unsigned fast_vb_size (void) { return sizeof (VBlockFAST); }
unsigned fast_vb_zip_dl_size (void) { return sizeof (ZipDataLineFAST); }

void fast_vb_release_vb (VBlock *vb_)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;

    vb->last_line = 0;
    vb->fasta_prev_vb_last_line_was_grepped = 0;
    
    memset (&vb->desc_mapper, 0, sizeof (vb->desc_mapper));
}

void fast_seg_initialize (VBlockFAST *vb)
{
    // thread safety: this will be initialized by vb_i=1, while it holds a mutex in zip_compress_one_vb
    static bool structured_initialized = false;
    if (!structured_initialized) {
        seg_initialize_compound_structured ((VBlockP)vb, "D?ESC", &structured_DESC); 
        structured_initialized = true;
    }

    if (vb->data_type == DT_FASTQ) {
        vb->mtf_ctx[FASTQ_SEQ].flags  = CTX_FL_LOCAL_LZMA;
        vb->mtf_ctx[FASTQ_SEQ].ltype  = CTX_LT_SEQUENCE;
        vb->mtf_ctx[FASTQ_QUAL].ltype = CTX_LT_SEQUENCE;
    }
    else {
        MtfContext *ctx = mtf_get_ctx (vb, (DictIdType)dict_id_FASTA_SEQ);
        ctx->flags  = CTX_FL_LOCAL_LZMA;
        ctx->ltype  = CTX_LT_SEQUENCE;
    }
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *fastq_seg_txt_line (VBlockFAST *vb, const char *field_start_line, bool *has_special_eol)     // index in vb->txt_data where this line starts
{
    ZipDataLineFAST *dl = DATA_LINE (vb->line_i);

    const char *next_field, *field_start=field_start_line;
    unsigned field_len=0;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n

    int32_t len = (int32_t)(AFTERENT (char, vb->txt_data) - field_start_line);

    // the leading @ - just verify it (it will be included in D0ESC subfield)
    ASSSEG (*field_start == '@', field_start, "%s: Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            global_cmd, *field_start);

    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG (<-- this is Illumina format)
    // See here for details of Illumina subfields: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    next_field = seg_get_next_line (vb, field_start, &len, &field_len, &has_13, "DESC");
 
    // we segment it using / | : and " " as separators. 
    seg_compound_field ((VBlockP)vb, &vb->mtf_ctx[FASTQ_DESC], field_start, field_len, &vb->desc_mapper, structured_DESC, true, 0);
    SEG_EOL (FASTQ_E1L);

    // SEQ - just get the whole line
    const char *seq_start = next_field;
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (vb, next_field, &len, true, false, false, &dl->seq_len, &separator, &has_13, "SEQ");
    vb->mtf_ctx[FASTQ_SEQ].local.len += dl->seq_len;
    
    // Add LOOKUP snip with seq_len
    char snip[10];
    snip[0] = SNIP_LOOKUP;
    unsigned seq_len_str_len = str_int (dl->seq_len, &snip[1]);
    seg_by_did_i (vb, snip, 1 + seq_len_str_len, FASTQ_SEQ, dl->seq_len); // we account for the '+' \n and maybe \r

    SEG_EOL (FASTQ_E1L);

    // PLUS - next line is expected to be a "+"
    field_start = next_field;
    next_field = seg_get_next_item (vb, field_start, &len, true, false, false, &field_len, &separator, &has_13, "+");
    ASSSEG (*field_start=='+' && field_len==1, field_start, "%s: Invalid FASTQ file format: expecting middle line to be a \"+\" (with no spaces) but it is \"%.*s\"",
            global_cmd, field_len, field_start);

    seg_by_did_i (vb, "+", 1, FASTQ_PLUS, 1);
    SEG_EOL (FASTQ_E1L);

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = next_field - vb->txt_data.data;
    unsigned qual_len;
    field_start = next_field;
    next_field = seg_get_next_item (vb, next_field, &len, true, false, false, &qual_len, &separator, &has_13, "QUAL");
    vb->mtf_ctx[FASTQ_QUAL].local.len += dl->seq_len;
    vb->mtf_ctx[FASTQ_QUAL].txt_len   += dl->seq_len;
    SEG_EOL (FASTQ_E1L);

    ASSSEG (qual_len == dl->seq_len, field_start, "%s: Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ=%.*s\nQUAL==%.*s",
            global_cmd, dl->seq_len, qual_len, dl->seq_len, seq_start, qual_len, field_start);
 
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
const char *fasta_seg_txt_line (VBlockFAST *vb, const char *line_start, bool *has_special_eol) // index in vb->txt_data where this line starts
{
    // get entire line
    unsigned line_len;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n
    int32_t remaining_vb_txt_len = AFTERENT (char, vb->txt_data) - line_start;
    const char *next_field = seg_get_next_line (vb, line_start, &remaining_vb_txt_len, &line_len, &has_13, "FASTA line");
    char redirect_snip[100];
    unsigned redirect_snip_len;

    // case: description line - we segment it to its components
    // note: we store the DESC structured in its own ctx rather than just directly in LINEMETA, to make it easier to grep
    if (*line_start == '>' || (*line_start == ';' && vb->last_line == FASTA_LINE_SEQ)) {
        // we segment using / | : and " " as separators. 
        seg_compound_field ((VBlockP)vb, mtf_get_ctx (vb, (DictIdType)dict_id_FASTA_DESC), 
                            line_start, line_len, &vb->desc_mapper, structured_DESC, true, 0);
        
        seg_prepare_snip_other (SNIP_REDIRECTION, (DictIdType)dict_id_FASTA_DESC, 0, redirect_snip, &redirect_snip_len);
        seg_by_did_i (vb, redirect_snip, redirect_snip_len, FASTA_LINEMETA, 0);

        SEG_EOL (FASTA_EOL);
        vb->last_line = FASTA_LINE_DESC;
    }

    // case: comment line - stored in the comment buffer
    else if (*line_start == ';' || !line_len) {
        seg_add_to_local_text ((VBlockP)vb, mtf_get_ctx (vb, (DictIdType)dict_id_FASTA_COMMENT), 
                               line_start, line_len, line_len); 

        seg_prepare_snip_other (SNIP_OTHER_LOOKUP, (DictIdType)dict_id_FASTA_COMMENT, 0, redirect_snip, &redirect_snip_len);
        seg_by_did_i (vb, redirect_snip, redirect_snip_len, FASTA_LINEMETA, 0);
        
        SEG_EOL (FASTA_EOL);
        vb->last_line = FASTA_LINE_COMMENT;
    }

    // case: sequence line
    else {
        DATA_LINE (vb->line_i)->seq_data_start = line_start - vb->txt_data.data;
        DATA_LINE (vb->line_i)->seq_len        = line_len;

        MtfContext *seq_ctx = mtf_get_ctx (vb, (DictIdType)dict_id_FASTA_SEQ);
        seq_ctx->local.len += line_len;
        seq_ctx->txt_len   += line_len;

        seg_prepare_snip_other (SNIP_OTHER_LOOKUP, (DictIdType)dict_id_FASTA_SEQ, line_len, redirect_snip, &redirect_snip_len);
        seg_by_did_i (vb, redirect_snip, redirect_snip_len, FASTA_LINEMETA, 0);

        SEG_EOL (FASTA_EOL);
        vb->last_line = FASTA_LINE_SEQ;
    }

    return next_field;
}

// callback function for compress to get data of one line (called by comp_lzma_data_in_callback)
void fast_zip_get_start_len_line_i_seq (VBlock *vb, uint32_t vb_line_i, 
                                        char **line_seq_data, uint32_t *line_seq_len,  // out 
                                        char **unused_data,  uint32_t *unused_len)
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);
    *line_seq_data = ENT (char, vb->txt_data, dl->seq_data_start);
    *line_seq_len  = dl->seq_len;
    *unused_data   = NULL;
    *unused_len    = 0;
}   

// callback function for compress to get data of one line (called by comp_compress_bzlib)
void fast_zip_get_start_len_line_i_qual (VBlock *vb, uint32_t vb_line_i, 
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

// called by I/O thread in fast_piz_read_one_vb, in case of --grep, to decompress and reconstruct the desc line, to 
// see if this vb is included. 
static bool fast_piz_test_grep (VBlockFAST *vb)
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

    MtfContext *desc_ctx = (vb->data_type == DT_FASTQ) ? &vb->mtf_ctx[FASTQ_DESC] : mtf_get_ctx (vb, (DictIdType)dict_id_FASTA_DESC);
    desc_ctx->iterator.next_b250 = FIRSTENT (uint8_t, desc_ctx->b250); 

    vb->line_i = vb->data_type == DT_FASTQ ? 4 * vb->first_line : vb->first_line;

    while (desc_ctx->iterator.next_b250 < AFTERENT (uint8_t, desc_ctx->b250)) {
        piz_reconstruct_from_ctx (vb, desc_ctx->did_i, 0);

        *AFTERENT (char, vb->txt_data) = 0; // terminate the desc string

        match = !!strstr (vb->txt_data.data, flag_grep);

        vb->txt_data.len = 0; // reset

        if (match) { // 
            found = true; // we've found a match to the grepped string
            if (vb->data_type == DT_FASTQ) break; // for FASTA, we need to go until the last line, for FASTQ, we can break here
        }

        if (vb->data_type == DT_FASTQ) vb->line_i += 4; // note: for FASTA we have no idea what txt line we're on, because we're only tracking DESC lines
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

    // reset iterators - piz_fast*_reconstruct_vb will use them again 
    mtf_init_iterator (desc_ctx);
    for (unsigned sf_i=0; sf_i < vb->desc_mapper.num_subfields; sf_i++) 
        mtf_init_iterator (&vb->mtf_ctx[vb->desc_mapper.did_i[sf_i]]);

    return found; // no match found
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

        if (flag_header_one) continue; // this is invoked by --header-only (re-written to flag_header_one in piz_read_global_area)

        piz_reconstruct_from_ctx (vb, FASTQ_SEQ,  0);
        piz_reconstruct_from_ctx (vb, FASTQ_E2L,  0);
        piz_reconstruct_from_ctx (vb, FASTQ_PLUS, 0);
        piz_reconstruct_from_ctx (vb, FASTQ_E3L,  0);
        piz_reconstruct_from_ctx (vb, FASTQ_QUAL, 0);
        piz_reconstruct_from_ctx (vb, FASTQ_E4L,  0);

        // case: we're grepping, and this line doesn't match
        *AFTERENT (char, vb->txt_data) = 0; // for strstr
        if (vb->dont_show_curr_line || (flag_grep && !strstr (ENT (char, vb->txt_data, txt_data_start_line), flag_grep)))
            vb->txt_data.len = txt_data_start_line; // rollback
    }
}

void fasta_piz_reconstruct_vb (VBlockFAST *vb)
{
    if (!flag_grep) piz_map_compound_field ((VBlockP)vb, dict_id_is_fast_desc_sf, &vb->desc_mapper); // it not already done during grep

    vb->dont_show_curr_line = flag_grep && !vb->fasta_prev_vb_last_line_was_grepped; // if we're continuing the previous VB's sequence, we obey its grep status

    for (vb->line_i=vb->first_line; vb->line_i < vb->first_line + vb->lines.len; vb->line_i++) {
        piz_reconstruct_from_ctx (vb, FASTA_LINEMETA, 0);
        piz_reconstruct_from_ctx (vb, FASTA_EOL, 0);
    }
}

bool fast_piz_read_one_vb (VBlock *vb, SectionListEntry *sl)
{ 
    // if we're grepping we we uncompress and reconstruct the DESC from the I/O thread, and terminate here if this VB is to be skipped
    if (flag_grep && !fast_piz_test_grep ((VBlockFAST *)vb)) return false; 

    return true;
}
