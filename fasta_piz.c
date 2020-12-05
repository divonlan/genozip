// ------------------------------------------------------------------
//   fasta_piz.c
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

// Called by thread I/O to initialize for a new genozip file
void fasta_piz_initialize (void)
{
    fasta_piz_initialize_contig_grepped_out (0,0,0); // initialize 
}

// returns true if section is to be skipped reading / uncompressing
bool fasta_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT sections

    if (flag.reading_reference) return false;  // doesn't apply when using FASTA as a reference

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
    if (flag.header_only_fast && 
        (dict_id.num == dict_id_fields[FASTA_NONREF] || dict_id.num == dict_id_fields[FASTA_NONREF_X] || dict_id.num == dict_id_fields[FASTA_COMMENT]))
        return true;

    // when grepping by I/O thread - skipping all sections but DESC
    if ((flag.grep || flag.regions) && (vb->grep_stages == GS_TEST) && 
        dict_id.num != dict_id_fields[FASTA_DESC] && !dict_id_is_fast_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if ((flag.grep || flag.regions) && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == dict_id_fields[FASTA_DESC] || dict_id_is_fast_desc_sf (dict_id)))
        return true;

    return false;
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat) and each toplevel item
CONTAINER_FILTER_FUNC (fasta_piz_filter)
{
    VBlockFAST *fasta_vb = (VBlockFAST *)vb;

    // if --phylip, don't reconstruct the EOL after DESC
    if (flag.out_dt == DT_PHYLIP && 
        dict_id.num == dict_id_fields[FASTA_TOPLEVEL] &&
        item == 1 /* EOL */ && 
        fasta_vb->last_line == FASTA_LINE_DESC) 
        return false; 

    return true; // show this item as normal
}

// remove trailing newline before SEQ lines in case of --sequential. Note that there might be more than one newline
// in which case the subsequent newlines were part of empty COMMENT lines
static inline void fasta_piz_remove_trailing_newlines (VBlockFAST *vb)
{
    while (*LASTENT (char, vb->txt_data) == '\n' || *LASTENT (char, vb->txt_data) == '\r')
        vb->txt_data.len--;
}

// this is used for end-of-lines of a sequence line, that are not the last line of the sequence. we skip reconstructing
// the newline if the user selected --sequential
SPECIAL_RECONSTRUCTOR (fasta_piz_special_SEQ)
{
    VBlockFAST *fasta_vb = (VBlockFAST *)vb;

    bool is_first_seq_line_in_this_contig = snip[0] - '0';

    // --sequential - if this is NOT the first seq line in the contig, we delete the previous end-of-line
    // TO DO: this doesn't yet work across vblock boundaries
    if (flag.sequential && !is_first_seq_line_in_this_contig) 
        fasta_piz_remove_trailing_newlines (fasta_vb);

    // skip showing line if this contig is grepped - but consume it anyway
    if (fasta_vb->contig_grepped_out) vb->dont_show_curr_line = true;

    // in case of not showing the SEQ in the entire file - we can skip consuming it
    if (flag.header_only_fast) // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->dont_show_curr_line = true;     
    else 
        piz_reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip+1, snip_len-1, true);    

    // case: --sequencial, and this seq line is the line in the vb, and it continues in the next vb
    if (  flag.sequential && // if we are asked for a sequential SEQ
          vb->line_i - vb->first_line == vb->lines.len-1 && // and this is the last line in this vb 
          !vb->dont_show_curr_line && 
          random_access_does_last_chrom_continue_in_next_vb (vb->vblock_i)) // and this sequence continues in the next VB 
        fasta_piz_remove_trailing_newlines (fasta_vb); // then: delete final newline if this VB ends with a 

    fasta_vb->last_line = FASTA_LINE_SEQ;

    return false; // no new value
}

SPECIAL_RECONSTRUCTOR (fasta_piz_special_COMMENT)
{
    VBlockFAST *fasta_vb = (VBlockFAST *)vb;

    // skip showing comment line in case cases - but consume it anyway:
    if (  fasta_vb->contig_grepped_out || // 1. if this contig is grepped out
          flag.out_dt == DT_PHYLIP)       // 2. if we're outputting in Phylis format
        vb->dont_show_curr_line = true;

    // in case of not showing the COMMENT in the entire file (--header-only or this is a --reference) - we can skip consuming it
    if (flag.header_only_fast)  // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->dont_show_curr_line = true;     
    else 
        piz_reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, true);    

    fasta_vb->last_line = FASTA_LINE_COMMENT;

    return false; // no new value
}

// this is called by fast_piz_test_grep - it is called sequentially for all VBs by the I/O thread
// returns true if the last contig of the previous VB was grepped-in
bool fasta_piz_initialize_contig_grepped_out (VBlockFAST *vb, bool does_vb_have_any_desc, bool last_desc_in_this_vb_matches_grep)
{
    // we pass the info from one VB to the next using this static variable
    static bool prev_vb_last_contig_grepped_out = false; 
    static uint32_t prev_vb_i = 0;

    if (!vb) { // vb=0 means initialize
        prev_vb_last_contig_grepped_out = false;
        prev_vb_i = 0;
        return 0;
    }

    // we're continuing the contig in the previous VB - until DESC is encountered
    vb->contig_grepped_out = prev_vb_last_contig_grepped_out || // last contig of previous VB had last_desc_in_this_vb_matches_grep
                             (prev_vb_i + 1 < vb->vblock_i);    // previous VB was skipped in piz_one_file due to random_access_is_vb_included
    
    // update for use of next VB, IF this VB contains any DESC line, otherwise just carry forward the current value
    if (does_vb_have_any_desc) 
        prev_vb_last_contig_grepped_out = !last_desc_in_this_vb_matches_grep; 

    prev_vb_i = vb->vblock_i;
    
    return !vb->contig_grepped_out;
}

bool fasta_piz_is_grepped_out_due_to_regions (VBlockFAST *vb, const char *line_start)
{
    const char *chrom_name = line_start + 1;
    unsigned chrom_name_len = strcspn (chrom_name, " \t\r\n");
    WordIndex chrom_index = ctx_search_for_word_index (&vb->contexts[CHROM], chrom_name, chrom_name_len);
    return !regions_is_site_included (chrom_index, 1); // we check for POS 1 to include (or not) the whole contig
}

// Phylip format mandates exact 10 space-padded characters: http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html
static inline void fasta_piz_translate_desc_to_phylip (VBlock *vb, char *desc_start)
{
    uint32_t reconstructed_len = AFTERENT (const char, vb->txt_data) - desc_start;
    *AFTERENT (char, vb->txt_data) = 0; // nul-terminate

    const char *chrom_name = desc_start + 1;
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");
    
    memcpy (desc_start, chrom_name, MIN (chrom_name_len, 10));
    if (chrom_name_len < 10) memcpy (desc_start + chrom_name_len, "          ", 10-chrom_name_len); // pad with spaces

    if (reconstructed_len > 10) vb->txt_data.len -= reconstructed_len - 10; // we do it this way to avoid signed problems
    else                        vb->txt_data.len += 10 - reconstructed_len;
}

// shorten DESC to the first white space
static inline void fasta_piz_desc_header_one (VBlock *vb, char *desc_start)
{
    uint32_t reconstructed_len = AFTERENT (const char, vb->txt_data) - desc_start;
    *AFTERENT (char, vb->txt_data) = 0; // nul-terminate
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");
    
    vb->txt_data.len -= reconstructed_len - chrom_name_len;
}

SPECIAL_RECONSTRUCTOR (fasta_piz_special_DESC)
{
    VBlockFAST *fasta_vb = (VBlockFAST *)vb;
    fasta_vb->contig_grepped_out = false;

    char *desc_start = AFTERENT (char, vb->txt_data);
    piz_reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, true);    

    // if --grep: here we decide whether to show this contig or not
    if (flag.grep) {
        *AFTERENT (char, vb->txt_data) = 0; // for strstr
        fasta_vb->contig_grepped_out = !strstr (desc_start, flag.grep);
    }

    if (flag.regions && !fasta_vb->contig_grepped_out) 
        fasta_vb->contig_grepped_out = fasta_piz_is_grepped_out_due_to_regions (fasta_vb, desc_start);

    // note: this logic allows the to grep contigs even if --no-header 
    if (fasta_vb->contig_grepped_out || flag.no_header)
        fasta_vb->dont_show_curr_line = true;     

    if (flag.out_dt == DT_PHYLIP) 
        fasta_piz_translate_desc_to_phylip (vb, desc_start);

    if (flag.header_one)
        fasta_piz_desc_header_one (vb, desc_start);

    fasta_vb->last_line = FASTA_LINE_DESC;

    return false; // no new value
}

bool fasta_piz_read_one_vb (VBlock *vb, ConstSectionListEntryP sl)
{ 
    // if we're grepping we we uncompress and reconstruct the DESC from the I/O thread, and terminate here if this VB is to be skipped
    if ((flag.grep || flag.regions) && !fast_piz_test_grep ((VBlockFAST *)vb)) return false; 

    return true;
}

// create Phylip header line
TXTHEADER_TRANSLATOR (txtheader_fa2phy)
{
    // get length of contigs and error if they are not all the same length
    uint32_t contig_len = random_access_verify_all_contigs_same_length();

    bufprintf (evb, txtheader_buf, "%"PRIu64" %u\n", z_file->contexts[FASTA_CONTIG].word_list.len, contig_len);
}
