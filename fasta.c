// ------------------------------------------------------------------
//   fasta.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fasta.h"
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
#include "reconstruct.h"
#include "vblock.h"

#define dict_id_is_fasta_desc_sf dict_id_is_type_1
#define dict_id_fasta_desc_sf dict_id_type_1

typedef struct {
    uint32_t seq_data_start, seq_len; // regular fasta and make-reference: start & length within vb->txt_data
} ZipDataLineFASTA;

// IMPORTANT: if changing fields in VBlockFASTA, also update vb_fast_release_vb 
typedef struct VBlockFASTA {    
    VBLOCK_COMMON_FIELDS

    bool contig_grepped_out;
    // note: last_line is initialized to FASTA_LINE_SEQ (=0) so that a ; line as the first line of the VB is interpreted as a description, not a comment
    enum { FASTA_LINE_SEQ, FASTA_LINE_DESC, FASTA_LINE_COMMENT } last_line; // ZIP & PIZ

    uint32_t lines_this_contig;    // ZIP

    // caching of seq line snips
    uint32_t std_line_len;         // ZIP: determined by first non-first-line seq line in VB 
    WordIndex std_line_node_index; // ZIP: node index for non-first lines, with same length as std_line_len

} VBlockFASTA;

#define DATA_LINE(i) ENT (ZipDataLineFASTA, vb->lines, (i))

unsigned fasta_vb_size (void) { return sizeof (VBlockFASTA); }
unsigned fasta_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTA); }

void fasta_vb_release_vb (VBlockFASTA *vb)
{
    memset ((char *)vb + sizeof (VBlock), 0, sizeof (VBlockFASTA) - sizeof (VBlock)); // zero all data unique to VBlockFASTA
    vb->contexts[FASTA_NONREF].local.len = 0; // len might be is used even though buffer is not allocated (in make-ref)
}

void fasta_vb_destroy_vb (VBlockFASTA *vb) {}

// used by ref_make_create_range
void fasta_get_data_line (VBlockP vb_, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len)
{
    VBlockFASTA *vb = (VBlockFASTA *)vb_;

    *seq_data_start = DATA_LINE(line_i)->seq_data_start;
    *seq_len        = DATA_LINE(line_i)->seq_len;
}

//-------------------------
// TXTFILE stuff
//-------------------------

// returns true if txt_data[txt_i] is the end of a FASTA contig (= next char is '>' or end-of-file), false if not, 
// and -1 if more data (lower first_i) is needed 
static inline int fasta_is_end_of_config (VBlock *vb, uint32_t first_i,
                                          int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
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
    ASSERTE (*i >= 0 && *i < vb->txt_data.len, "*i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    // case: reference file - we allow only one contig (or part of it) per VB - move second contig onwards to next vb
    // (note: first_i=0 when flag.make_reference)
    if (flag.make_reference) {
        bool data_found = false;
        ARRAY (char, txt, vb->txt_data);
        for (uint32_t i=0; i < vb->txt_data.len; i++) {
            // just don't allow now-obsolete ';' rather than trying to disentangle comments from descriptions
            ASSINP (txt[i] != ';', "Error: %s contains a ';' character - this is not supported for reference files. Contig descriptions must begin with a >", txt_name);
        
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

            // when compressing FASTA with a reference or --multifasta - an "end of line" means "end of contig" - 
            // i.e. one that the next character is >, or it is the end of the file
            // note: when compressing FASTA with a reference (eg long reads stored in a FASTA instead of a FASTQ), 
            // line cannot be too long - they must fit in a VB
            if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE || flag.multifasta) 
                switch (fasta_is_end_of_config (vb, first_i, *i)) {
                    case true  : break;            // end of contig
                    case false : continue;         // not end of contig
                    default    : goto out_of_data; // need more data (lower first_i)
                }

            return vb->txt_data.len-1 - *i;
        }
    }

out_of_data:
    return -1; // cannot find end-of-line in the data starting first_i
}

//-------------------------
// SEG & ZIP stuff
//-------------------------

// callback function for compress to get data of one line (called by codec_lzma_data_in_callback)
void fasta_zip_seq (VBlock *vb, uint32_t vb_line_i, 
                    char **line_seq_data, uint32_t *line_seq_len,  // out 
                    uint32_t maximum_len)
{
    ZipDataLineFASTA *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_seq_len = MIN (dl->seq_len, maximum_len);
    
    if (line_seq_data) // if NULL, only length was requested
        *line_seq_data = dl->seq_len ? ENT (char, vb->txt_data, dl->seq_data_start) : NULL;
}   

void fasta_seg_initialize (VBlockFASTA *vb)
{
    START_TIMER;

    ASSINP (vb->vblock_i > 1 || *FIRSTENT (char, vb->txt_data) == '>' || *FIRSTENT (char, vb->txt_data) == ';',
            "Error: expecting FASTA file %s to start with a '>' or a ';'", txt_name);

    if (!flag.make_reference) {

        vb->contexts[FASTA_LINEMETA].no_stons = true; // avoid edge case where entire b250 is moved to local due to singletons, because fasta_reconstruct_vb iterates on ctx->b250

        // if this FASTA contains unrelated contigs, we're better off with ACGT        
        if (!flag.multifasta)
            codec_acgt_comp_init ((VBlockP)vb);

        // if the contigs in this FASTA are related, LZMA will do a better job (but it will consume a lot of memory)
        else {
            vb->contexts[FASTA_NONREF].lcodec = CODEC_LZMA;
            vb->contexts[FASTA_NONREF].ltype  = LT_SEQUENCE;
        }
  
        if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) 
            vb->contexts[FASTA_NONREF].no_callback = true; // override callback if we are segmenting to a reference

        // in --stats, consolidate stats into FASTA_NONREF
        stats_set_consolidation ((VBlockP)vb, FASTA_NONREF, 1, FASTA_NONREF_X);
    }
    else { // make-reference
        vb->contexts[FASTA_CONTIG].no_vb1_sort = true; // keep contigs in the order of the reference, i.e. in the order they would appear in BAM header created with this reference
    }

    vb->contexts[FASTA_CONTIG  ].no_stons = true; // needs b250 node_index for reference
    vb->contexts[FASTA_TOPLEVEL].no_stons = true; // keep in b250 so it can be eliminated as all_the_same

    COPY_TIMER (seg_initialize);
}

void fasta_seg_finalize (VBlockP vb)
{
    // top level snip
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .filter_items = true,
        .nitems_lo    = 2,
        .items        = { { .dict_id = (DictId)dict_id_fields[FASTA_LINEMETA]  },
                          { .dict_id = (DictId)dict_id_fields[FASTA_EOL]       } }
    };

    container_seg_by_ctx (vb, &vb->contexts[FASTA_TOPLEVEL], (ContainerP)&top_level, 0, 0, 0);
}

bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == dict_id_fields[FASTA_TOPLEVEL] ||
           dict_id.num == dict_id_fields[FASTA_DESC]     ||
           dict_id.num == dict_id_fields[FASTA_LINEMETA] ||
           dict_id.num == dict_id_fields[FASTA_EOL];
}

// description line - we segment it to its components
// note: we store the DESC container in its own ctx rather than just directly in LINEMETA, to make it easier to grep
static void fast_seg_desc_line (VBlockFASTA *vb, const char *line_start, uint32_t line_len, bool *has_13)
{
    // we store the contig name in a dictionary only (no b250), to be used if this fasta is used as a reference
    const char *chrom_name = line_start + 1;
    unsigned chrom_name_len = strcspn (line_start + 1, " \t\r\n");

    ASSSEG0 (chrom_name_len, line_start, "contig is missing a name");

    if (!flag.make_reference) {
        SegCompoundArg arg = { .slash = true, .pipe = true, .dot = true, .colon = true, .whitespace = true };
        seg_compound_field ((VBlockP)vb, &vb->contexts[FASTA_DESC], line_start, line_len, arg, 0, 0);
        
        char special_snip[100]; unsigned special_snip_len;
        seg_prepare_snip_other (SNIP_REDIRECTION, (DictId)dict_id_fields[FASTA_DESC], false, 0, &special_snip[2], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_DESC;

        seg_by_did_i (vb, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    bool is_new;
    WordIndex chrom_node_index = ctx_evaluate_snip_seg ((VBlockP)vb, &vb->contexts[FASTA_CONTIG], chrom_name, chrom_name_len, &is_new);
    random_access_update_chrom ((VBlockP)vb, chrom_node_index, chrom_name, chrom_name_len);

    ASSINP (is_new, "Error: bad FASTA file - contig \"%.*s\" appears more than once%s", chrom_name_len, chrom_name,
            flag.bind ? " (possibly in another FASTA being bound)" : 
            (flag.reference==REF_EXTERNAL || flag.reference==REF_EXT_STORE) ? " (possibly the contig size exceeds vblock size, try enlarging with --vblock)" : "");
        
    vb->last_line = FASTA_LINE_DESC;    
}

static void fast_seg_comment_line (VBlockFASTA *vb, const char *line_start, uint32_t line_len, bool *has_13)
{
    if (!flag.make_reference) {
        seg_add_to_local_text ((VBlockP)vb, &vb->contexts[FASTA_COMMENT], line_start, line_len, line_len); 

        char special_snip[100]; unsigned special_snip_len;
        seg_prepare_snip_other (SNIP_OTHER_LOOKUP, (DictId)dict_id_fields[FASTA_COMMENT], false, 0, &special_snip[2], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_COMMENT;

        seg_by_did_i (vb, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    vb->last_line = FASTA_LINE_COMMENT;
}

static void fasta_seg_seq_line_do (VBlockFASTA *vb, uint32_t line_len, bool is_first_line_in_contig)
{
    Context *lm_ctx  = &vb->contexts[FASTA_LINEMETA];
    Context *seq_ctx = &vb->contexts[FASTA_NONREF];

    // cached node index
    if (!is_first_line_in_contig && line_len == vb->std_line_len) {
        buf_alloc_more (vb, &lm_ctx->b250, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->b250");
        NEXTENT (WordIndex, lm_ctx->b250) = vb->std_line_node_index;
    }

    else { // not cached
        char special_snip[100]; unsigned special_snip_len;
        seg_prepare_snip_other (SNIP_OTHER_LOOKUP, (DictId)dict_id_fields[FASTA_NONREF], 
                                true, (int32_t)line_len, &special_snip[3], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_SEQ;
        special_snip[2] = '0' + is_first_line_in_contig; 
        WordIndex node_index = seg_by_ctx (vb, special_snip, 3 + special_snip_len, lm_ctx, 0);  // the payload of the special snip, is the OTHER_LOOKUP snip...

        // create cache if needed
        if (!is_first_line_in_contig && !vb->std_line_len) {
            vb->std_line_len = line_len;
            vb->std_line_node_index = node_index;
        }
    }

    seq_ctx->txt_len   += line_len;
    seq_ctx->local.len += line_len;
} 

static void fasta_seg_seq_line (VBlockFASTA *vb, const char *line_start, uint32_t line_len, bool is_last_line_in_contig, bool *has_13)
{
    vb->lines_this_contig++;

    *DATA_LINE (vb->line_i) = (ZipDataLineFASTA){ .seq_data_start = line_start - vb->txt_data.data,
                                                  .seq_len        = line_len };

    if (flag.make_reference)
        vb->contexts[FASTA_NONREF].local.len += line_len;

    // after last line in contig, we know if its SIMPLE or PBWT and seg all lines in this contig into FASTA_LINEMETA
    if (!flag.make_reference && is_last_line_in_contig) {

        ZipDataLineFASTA *dl = DATA_LINE (vb->line_i - vb->lines_this_contig + 1);

        for (int32_t i=0; i < vb->lines_this_contig; i++) {
            fasta_seg_seq_line_do (vb, dl[i].seq_len, i==0);
            SEG_EOL (FASTA_EOL, true); 
        }        
        vb->lines_this_contig = 0;
    }

    vb->last_line = FASTA_LINE_SEQ;

    // case: this sequence is continuation from the previous VB - we don't yet know the chrom - we will update it,
    // and increment the min/max_pos relative to the beginning of the seq in the vb later, in random_access_finalize_entries
    if (!vb->chrom_name && vb->line_i == 0)
        random_access_update_chrom ((VBlockP)vb, WORD_INDEX_NONE, 0, 0);

    random_access_increment_last_pos ((VBlockP)vb, line_len); 
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
const char *fasta_seg_txt_line (VBlockFASTA *vb, const char *line_start, uint32_t remaining_txt_len, bool *has_13) // index in vb->txt_data where this line starts
{
    // get entire line
    unsigned line_len;
    int32_t remaining_vb_txt_len = AFTERENT (char, vb->txt_data) - line_start;
    const char *next_field = seg_get_next_line (vb, line_start, &remaining_vb_txt_len, &line_len, has_13, "FASTA line");

    // case: description line - we segment it to its components
    if (*line_start == '>' || (*line_start == ';' && vb->last_line == FASTA_LINE_SEQ)) 
        fast_seg_desc_line (vb, line_start, line_len, has_13);

    // case: comment line - stored in the comment buffer
    else if (*line_start == ';' || !line_len) 
        fast_seg_comment_line (vb, line_start, line_len, has_13);

    // case: sequence line
    else 
        fasta_seg_seq_line (vb, line_start, line_len, 
                            !remaining_vb_txt_len || *next_field == '>' || *next_field == ';' || *next_field == '\n' || *next_field == '\r', // is_last_line_in_contig
                            has_13);

    return next_field;
}

//-------------------------
// PIZ stuff
//-------------------------

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
        dict_id.num != dict_id_fields[FASTA_DESC] && !dict_id_is_fasta_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if ((flag.grep || flag.regions) && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == dict_id_fields[FASTA_DESC] || dict_id_is_fasta_desc_sf (dict_id)))
        return true;

    return false;
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat) and each toplevel item
CONTAINER_FILTER_FUNC (fasta_piz_filter)
{
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;

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
static inline void fasta_piz_remove_trailing_newlines (VBlockFASTA *vb)
{
    while (*LASTENT (char, vb->txt_data) == '\n' || *LASTENT (char, vb->txt_data) == '\r')
        vb->txt_data.len--;
}

// this is used for end-of-lines of a sequence line, that are not the last line of the sequence. we skip reconstructing
// the newline if the user selected --sequential
SPECIAL_RECONSTRUCTOR (fasta_piz_special_SEQ)
{
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;

    bool is_first_seq_line_in_this_contig = snip[0] - '0';

    // --sequential - if this is NOT the first seq line in the contig, we delete the previous end-of-line
    if (flag.sequential && !is_first_seq_line_in_this_contig) 
        fasta_piz_remove_trailing_newlines (fasta_vb);

    // skip showing line if this contig is grepped - but consume it anyway
    if (fasta_vb->contig_grepped_out) vb->dont_show_curr_line = true;

    // in case of not showing the SEQ in the entire file - we can skip consuming it
    if (flag.header_only_fast) // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->dont_show_curr_line = true;     
    else 
        reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip+1, snip_len-1, true);    

    // case: --sequential, and this seq line is the last line in the vb, and it continues in the next vb
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
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;

    // skip showing comment line in case cases - but consume it anyway:
    if (  fasta_vb->contig_grepped_out || // 1. if this contig is grepped out
          flag.out_dt == DT_PHYLIP)       // 2. if we're outputting in Phylis format
        vb->dont_show_curr_line = true;

    // in case of not showing the COMMENT in the entire file (--header-only or this is a --reference) - we can skip consuming it
    if (flag.header_only_fast)  // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->dont_show_curr_line = true;     
    else 
        reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, true);    

    fasta_vb->last_line = FASTA_LINE_COMMENT;

    return false; // no new value
}

// this is called by piz_test_grep - it is called sequentially for all VBs by the I/O thread
// returns true if the last contig of the previous VB was grepped-in
bool fasta_piz_initialize_contig_grepped_out (VBlock *vb_, bool does_vb_have_any_desc, bool last_desc_in_this_vb_matches_grep)
{
    VBlockFASTA *vb = (VBlockFASTA *)vb_;

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
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;
    fasta_vb->contig_grepped_out = false;

    char *desc_start = AFTERENT (char, vb->txt_data);
    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, true);    
    *AFTERENT (char, vb->txt_data) = 0; // for strstr and strcspn

    // if --grep: here we decide whether to show this contig or not
    if (flag.grep) 
        fasta_vb->contig_grepped_out = !strstr (desc_start, flag.grep);
    
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");
    vb->chrom_node_index = ctx_search_for_word_index (&vb->contexts[CHROM], desc_start + 1, chrom_name_len);

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
    if ((flag.grep || flag.regions) && !piz_test_grep (vb)) return false; 

    return true;
}

// create Phylip header line
TXTHEADER_TRANSLATOR (txtheader_fa2phy)
{
    // get length of contigs and error if they are not all the same length
    uint32_t contig_len = random_access_verify_all_contigs_same_length();

    bufprintf (evb, txtheader_buf, "%"PRIu64" %u\n", z_file->contexts[FASTA_CONTIG].word_list.len, contig_len);
}
