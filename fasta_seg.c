// ------------------------------------------------------------------
//   fasta_seg.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fast_private.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "piz.h"
#include "dict_id.h"
#include "random_access.h"
#include "strings.h"
#include "regions.h"
#include "codec.h"
#include "stats.h"

// returns true if txt_data[txt_i] is the end of a FASTA read (= next char is > or end-of-file), false if not, 
// and -1 if more data (lower first_i) is needed 
static inline int fasta_is_end_of_line (VBlock *vb, uint32_t first_i,
                                        int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    // if it is a @ then that @ is a new record EXCEPT if the previous character is + and then
    // @ is actually a quality value... (we check two previous characters, as the previous one might be \r)
    if (txt_i < vb->txt_data.len-1)
        return txt[txt_i+1] == '>';

    // if we're at the end of the line, we scan back to the previous \n and check if it NOT the >
    for (int32_t i=txt_i-1; i >= first_i; i--) 
        if (txt[i] == '\n') 
            return txt[txt_i+1] != '>'; // this row is a sequence row, not a description row

    return -1; // we need more data (lower first_i)
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i)
{
    ASSERT (*i >= 0 && *i < vb->txt_data.len, "Error in def_unconsumed: *i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    // case: reference file - we allow only one contig (or part of it) per VB - move second contig onwards to next vb
    // (note: first_i=0 when flag.make_reference)
    if (flag.make_reference) {
        bool data_found = false;
        ARRAY (char, txt, vb->txt_data);
        for (uint32_t i=0; i < vb->txt_data.len; i++) {
            // just don't allow now-obsolete ';' rather than trying to disentangle comments from descriptions
            ASSERT (txt[i] != ';', "Error: %s contains a ';' character - this is not supported for reference files. Contig descriptions must begin with a >", txt_name);
        
            // if we've encountered a new DESC line after already seeing sequence data, move this DESC line and
            // everything following to the next VB
            if (data_found && txt[i]=='>' && txt[i-1]=='\n') 
                return vb->txt_data.len - i;

            if (!data_found && (txt[i] != '\n' && txt[i] != '\r')) data_found = true; // anything, except for empty lines, is considered data
        }
    }

    // we move the final partial line to the next vb (unless we are already moving more, due to a reference file)
    for (; *i >= (int32_t)first_i; (*i)--) {

        if (vb->txt_data.data[*i] == '\n') {

            // when compressing FASTA with a reference - an "end of line" is one that the next character is >, or it is the end of the file
            // note: when compressing FASTA with a reference (eg long reads stored in a FASTA instead of a FASTQ), line cannot be too long - they 
            // must fit in a VB
            if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) 
                switch (fasta_is_end_of_line (vb, first_i, *i)) {
                    case true  : break;            // end of line
                    case false : continue;         // not end of line
                    default    : goto out_of_data; // need more data (lower first_i)
                }

            return vb->txt_data.len-1 - *i;
        }
    }

out_of_data:
    return -1; // cannot find end-of-line in the data starting first_i
}

// callback function for compress to get data of one line (called by codec_lzma_data_in_callback)
void fasta_zip_seq (VBlock *vb, uint32_t vb_line_i, 
                    char **line_seq_data, uint32_t *line_seq_len,  // out 
                    uint32_t maximum_len)
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_seq_len  = MIN (dl->seq_len, maximum_len);
    
    if (line_seq_data) // if NULL, only length was requested
        *line_seq_data = dl->seq_len ? ENT (char, vb->txt_data, dl->seq_data_start) : NULL;
}   

void fasta_seg_initialize (VBlockFAST *vb)
{
    START_TIMER;

    ASSERT (vb->vblock_i > 1 || *FIRSTENT (char, vb->txt_data) == '>' || *FIRSTENT (char, vb->txt_data) == ';',
            "Error: expecting FASTA file %s to start with a '>' or a ';'", txt_name);

    if (!flag.make_reference) {

        vb->contexts[FASTA_LINEMETA].no_stons = true; // avoid edge case where entire b250 is moved to local due to singletons, because fasta_reconstruct_vb iterates on ctx->b250
        
        codec_acgt_comp_init ((VBlockP)vb);

        if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) {
            vb->contexts[FASTA_NONREF].no_callback = true; // override callback if we are segmenting to a reference

            // in --stats, consolidate stats into SQBITMAP
            stats_set_consolidation ((VBlockP)vb, FASTA_SQBITMAP, 4, FASTA_NONREF, FASTA_NONREF_X, FASTA_GPOS, FASTA_STRAND);
        }
    }
    else { // make-reference
        vb->contexts[FASTA_CONTIG].no_vb1_sort = true; // keep contigs in the order of the reference, i.e. in the order they would appear in BAM header created with this reference
    }

    vb->contexts[FASTA_CONTIG].no_stons = true; // needs b250 node_index for reference
    vb->contexts[FASTA_TOPLEVEL].no_stons = true; // keep in b250 so it can be eliminated as all_the_same

    COPY_TIMER (seg_initialize);
}

void fasta_seg_finalize (VBlockP vb)
{
    // top level snip
    Container top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .num_items    = 2,
        .items        = { { (DictId)dict_id_fields[FASTA_LINEMETA], DID_I_NONE, ""   },
                          { (DictId)dict_id_fields[FASTA_EOL],      DID_I_NONE, ""   } }
    };

    container_seg_by_ctx (vb, &vb->contexts[FASTA_TOPLEVEL], &top_level, 0, 0, 0);
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
const char *fasta_seg_txt_line (VBlockFAST *vb, const char *line_start, uint32_t remaining_txt_len, bool *has_13) // index in vb->txt_data where this line starts
{
    // get entire line
    unsigned line_len;
    int32_t remaining_vb_txt_len = AFTERENT (char, vb->txt_data) - line_start;
    const char *next_field = seg_get_next_line (vb, line_start, &remaining_vb_txt_len, &line_len, has_13, "FASTA line");
    char special_snip[100];
    unsigned special_snip_len;

    // case: description line - we segment it to its components
    // note: we store the DESC container in its own ctx rather than just directly in LINEMETA, to make it easier to grep
    if (*line_start == '>' || (*line_start == ';' && vb->last_line == FASTA_LINE_SEQ)) {

        // we store the contig name in a dictionary only (no b250), to be used if this fasta is used as a reference
        const char *chrom_name = line_start + 1;
        unsigned chrom_name_len = strcspn (line_start + 1, " \t\r\n");

        ASSSEG0 (chrom_name_len, line_start, "contig is missing a name");

        if (!flag.make_reference) {
            // we segment using / | : . and " " as separators. 

            seg_compound_field ((VBlockP)vb, &vb->contexts[FASTA_DESC], line_start, line_len, true, 0, 0);
            
            seg_prepare_snip_other (SNIP_REDIRECTION, (DictId)dict_id_fields[FASTA_DESC], false, 0, &special_snip[2], &special_snip_len);

            special_snip[0] = SNIP_SPECIAL;
            special_snip[1] = FASTA_SPECIAL_DESC;

            seg_by_did_i (vb, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
            SEG_EOL (FASTA_EOL, true);
        }

        bool is_new;
        WordIndex chrom_node_index = ctx_evaluate_snip_seg ((VBlockP)vb, &vb->contexts[FASTA_CONTIG], chrom_name, chrom_name_len, &is_new);
        random_access_update_chrom ((VBlockP)vb, chrom_node_index, chrom_name, chrom_name_len);

        ASSERT (is_new, "Error: bad FASTA file - contig \"%.*s\" appears more than once%s", chrom_name_len, chrom_name,
                flag.bind ? " (possibly in another FASTA being bound)" : 
                (flag.reference==REF_EXTERNAL || flag.reference==REF_EXT_STORE) ? " (possibly the contig size exceeds vblock size, try enlarging with --vblock)" : "");
            
        vb->last_line = FASTA_LINE_DESC;
    }

    // case: comment line - stored in the comment buffer
    else if (*line_start == ';' || !line_len) {
        
        if (!flag.make_reference) {
            seg_add_to_local_text ((VBlockP)vb, &vb->contexts[FASTA_COMMENT], line_start, line_len, line_len); 

            seg_prepare_snip_other (SNIP_OTHER_LOOKUP, (DictId)dict_id_fields[FASTA_COMMENT], false, 0, &special_snip[2], &special_snip_len);

            special_snip[0] = SNIP_SPECIAL;
            special_snip[1] = FASTA_SPECIAL_COMMENT;

            seg_by_did_i (vb, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
            SEG_EOL (FASTA_EOL, true);
        }

        vb->last_line = FASTA_LINE_COMMENT;
    }

    // case: sequence line
    else {
        DATA_LINE (vb->line_i)->seq_data_start = line_start - vb->txt_data.data;
        DATA_LINE (vb->line_i)->seq_len        = line_len;

        Context *nonref_ctx = &vb->contexts[FASTA_NONREF];
        nonref_ctx->local.len += line_len;

        if (!flag.make_reference) {

            nonref_ctx->txt_len += line_len;

            seg_prepare_snip_other (SNIP_OTHER_LOOKUP, (DictId)dict_id_fields[FASTA_NONREF], true, (int32_t)line_len, &special_snip[3], &special_snip_len);

            special_snip[0] = SNIP_SPECIAL;
            special_snip[1] = FASTA_SPECIAL_SEQ;
            special_snip[2] = '0' + (vb->last_line != FASTA_LINE_SEQ); // first seq line in this contig
            seg_by_did_i (vb, special_snip, 3 + special_snip_len, FASTA_LINEMETA, 0);  // the payload of the special snip, is the OTHER_LOOKUP snip...

            SEG_EOL (FASTA_EOL, true); 

//            if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) {
//                aligner_seg_seq (vb, line_start, line_len, FASTA_SQBITMAP);
//                seq_ctx->local.len = 0; // we don't use FASTA_NONREF if segging to a reference
//            }
        }

        vb->last_line = FASTA_LINE_SEQ;

        // case: this sequence is continuation from the previous VB - we don't yet know the chrom - we will update it,
        // and increment the min/max_pos relative to the beginning of the seq in the vb later, in random_access_finalize_entries
        if (!vb->chrom_name && vb->line_i == 0)
            random_access_update_chrom ((VBlockP)vb, WORD_INDEX_NONE, 0, 0);

        random_access_increment_last_pos ((VBlockP)vb, line_len); 
    }

    return next_field;
}
