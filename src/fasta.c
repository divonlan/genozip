// ------------------------------------------------------------------
//   fasta.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fasta_private.h"
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
#include "kraken.h"
#include "reference.h"
#include "segconf.h"
#include "chrom.h"
#include "tokenizer.h"
#include "writer.h" 

#define dict_id_is_fasta_desc_sf dict_id_is_type_1
#define dict_id_fasta_desc_sf dict_id_type_1

typedef struct {
    uint32_t seq_data_start, seq_len; // regular fasta and make-reference: start & length within vb->txt_data
    bool has_13;
} ZipDataLineFASTA;
#define DATA_LINE(i) B(ZipDataLineFASTA, VB_FASTA->lines, (i))

unsigned fasta_vb_size (DataType dt) 
{ 
    return dt == DT_REF && IS_PIZ ? sizeof (VBlock) : sizeof (VBlockFASTA); 
}

unsigned fasta_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTA); }

void fasta_vb_release_vb (VBlockFASTAP vb)
{
    if (VB_DT(REF) && IS_PIZ) return; // this is actually a VBlock, not VBlockFASTA
    
    memset ((char *)vb + sizeof (VBlock), 0, sizeof (VBlockFASTA) - sizeof (VBlock)); // zero all data unique to VBlockFASTA
    CTX(FASTA_NONREF)->local.len = 0; // len might be is used even though buffer is not allocated (in make-ref)
}

void fasta_vb_destroy_vb (VBlockFASTAP vb) {}

// used by ref_make_create_range
void fasta_get_data_line (VBlockP vb, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len)
{
    *seq_data_start = DATA_LINE(line_i)->seq_data_start;
    *seq_len        = DATA_LINE(line_i)->seq_len;
}

//-------------------------
// TXTFILE stuff
//-------------------------

// returns true if txt_data[txt_i] is the end of a FASTA contig (= next char is '>' or end-of-file), false if not, 
// and -1 if more data (lower first_i) is needed 
static inline int fasta_is_end_of_contig (VBlockP vb, uint32_t first_i,
                                          int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    if (txt_i < vb->txt_data.len-1)
        return txt[txt_i+1] == '>';

    // if we're at the end of the line, we scan back to the previous \n and check if it NOT the >
    bool newline_run = true;
    for (int32_t i=txt_i-1; i >= first_i; i--) {
        if (txt[i] == '\n' && !newline_run) // ASSUMES NO COMMENT LINES (starting with ;), To do: fix this
            return txt[i+1] != '>'; // this row is a sequence row, not a description row
        
        else if (newline_run && txt[i] != '\n' && txt[i] != '\r')
            newline_run = false; // newline run ends when we encounter the first non newline
    }

    return -1; // we need more data (lower first_i)
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i, int32_t *last_i)
{
    bool is_entire_vb = (first_i == 0 && *last_i == vb->txt_data.len-1);

    ASSERT (*last_i >= 0 && *last_i < vb->txt_data.len, "*last_i=%d is out of range [0,%"PRIu64"]", *last_i, vb->txt_data.len);

    ARRAY (char, txt, vb->txt_data);

    // case: reference file - we allow only one contig (or part of it) per VB - move second contig onwards to next vb
    // (note: first_i=0 when flag.make_reference)
    if (flag.make_reference) {
        bool data_found = false;
        for (uint32_t i=0; i < vb->txt_data.len32; i++) {
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
    for (int32_t i=*last_i; i >= (int32_t)first_i; i--) {

        if (txt[i] == '\n') {

            // when compressing FASTA with a reference or --multiseq - an "end of line" means "end of contig" - 
            // i.e. one that the next character is >, or it is the end of the file
            // note: when compressing FASTA with a reference (eg long reads stored in a FASTA instead of a FASTQ), 
            // line cannot be too long - they must fit in a VB
            if ((flag.reference & REF_ZIP_LOADED) || flag.multiseq) {
                int is_end_of_contig = fasta_is_end_of_contig (vb, first_i, i);

                switch (is_end_of_contig) {
                    case true  : break;            // end of contig
                    case false : continue;         // not end of contig
                    default    : goto out_of_data; // need more data (lower first_i)
                }
            }

            // otherwise - tolerate a VB that ends part way through a SEQ
            else if (is_entire_vb && i+1 < vb->txt_data.len &&
                     txt[i+1] != ';' && txt[i+1] != '>') { // partial line isn't a Description or a Comment, hence its a Sequence
                ((VBlockFASTAP)vb)->vb_has_no_newline = true;
                return 0;                
            }
                
            *last_i = i;
            return vb->txt_data.len-1 - i;
        }
    }

out_of_data:
    // case: an entire FASTA VB without newlines (i.e. a very long sequential SEQ) - we accept a VB without newlines and deal with it in Seg
    if (is_entire_vb && !segconf.running) {
        ((VBlockFASTAP)vb)->vb_has_no_newline = true;
        return 0;
    }

    // case: a single BGZF block without newlines - test next block
    else 
        return -1; // cannot find end-of-line in the data starting first_i
}

//-------------------------
// SEG & ZIP stuff
//-------------------------

// called by main thread at the beginning of zipping this file
void fasta_zip_initialize (void)
{
}

void fasta_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeader *vb_header)
{
    if (Z_DT(GFF))
        vb_header->flags.vb_header.gff.embedded_fasta = true;
}

// callback function for compress to get data of one line 
COMPRESSOR_CALLBACK (fasta_zip_seq)
{
    ZipDataLineFASTA *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len = MIN_(dl->seq_len, maximum_size);
    
    if (line_data) // if NULL, only length was requested
        *line_data = dl->seq_len ? Bc (vb->txt_data, dl->seq_data_start) : NULL;

    if (is_rev) *is_rev = 0;
}   

void fasta_seg_initialize (VBlockP vb)
{
    START_TIMER;

    ASSINP (vb->vblock_i > 1 || *B1STtxt == '>' || *B1STtxt == ';',
            "Error: expecting FASTA file %s to start with a '>' or a ';' but seeing \"%.*s\"", txt_name, MIN_(64, vb->txt_data.len32), B1STtxt);

    CTX(FASTA_TOPLEVEL)->no_stons  = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(FASTA_CONTIG)->flags.store = STORE_INDEX; // since v12
    CTX(FASTA_CONTIG)->no_stons    = true; // needs b250 node_index for reference
    CTX(FASTA_LINEMETA)->no_stons  = true; // avoid edge case where entire b250 is moved to local due to singletons, because fasta_reconstruct_vb iterates on ctx->b250
    CTX(FASTA_COMMENT)->no_stons   = true;
    
    if (!segconf.fasta_has_contigs) 
        CTX(FASTA_DESC)->no_stons  = true;

    if (kraken_is_loaded) {
        CTX(FASTA_TAXID)->flags.store    = STORE_INT;
        CTX(FASTA_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTA_TAXID)->counts_section = true; 
    }

    // if this neocleotide FASTA of unrelated contigs, we're better off with ACGT        
    if (!flag.multiseq && segconf.seq_type == SQT_NUKE)
        codec_acgt_comp_init (VB, FASTA_NONREF);

    // if the contigs in this FASTA are related, let codec_assign_best_codec assign the bext codec 
    else 
        CTX(FASTA_NONREF)->ltype  = LT_SEQUENCE;

    if (flag.reference & REF_ZIP_LOADED) 
        CTX(FASTA_NONREF)->no_callback = true; // override callback if we are segmenting to a reference

    // in --stats, consolidate stats into FASTA_NONREF
    ctx_consolidate_stats (vb, FASTA_NONREF, FASTA_NONREF_X, DID_EOL);

    if (segconf.running)
        segconf.fasta_has_contigs = true; // initialize optimistically
        
    COPY_TIMER (seg_initialize);
}

void fasta_seg_finalize (VBlockP vb)
{
    if (!flag.make_reference) {
        // top level snip
        SmallContainer top_level = { 
            .repeats      = vb->lines.len,
            .is_toplevel  = true,
            .callback     = true,
            .nitems_lo    = 2,
            .items        = { { .dict_id = { _FASTA_LINEMETA }  },
                              { .dict_id = { _FASTA_EOL      }, .translator = FASTA2PHYLIP_EOL } }
        };

        container_seg (vb, CTX(FASTA_TOPLEVEL), (ContainerP)&top_level, 0, 0, 0);
    }

    // decide whether the sequences in this FASTA represent contigs (in which case we want a FASTA_CONTIG dictionary
    // and random access) or do they represent reads (in which case they are likely to numerous to be added to a dict)
    if (segconf.running) {
        uint64_t num_contigs_this_vb = CTX(FASTA_CONTIG)->nodes.len;
        ASSINP0 (num_contigs_this_vb, "Invalid FASTA file: no sequence description line");

        uint64_t avg_contig_size_this_vb = vb->txt_data.len / num_contigs_this_vb;
        uint64_t est_num_contigs_in_file = txtfile_get_seggable_size() / avg_contig_size_this_vb;

        // case: we've seen only characters that are both nucleotide and protein (as are A,C,G,T,N) - call it as nucleotide
        if (segconf.seq_type == SQT_NUKE_OR_AMINO) segconf.seq_type = SQT_NUKE;
        ASSINP0 (!flag.make_reference || segconf.seq_type == SQT_NUKE, "Can't use --make-reference on this file, because it contains amino acids instead of nucleotides");

        // limit the number of contigs, to avoid the FASTA_CONTIG dictionary becoming too big. note this also
        // sets a limit for fasta-to-phylip translation
        #define MAX_CONTIGS_IN_FILE 1000000 
        segconf.fasta_has_contigs &= (num_contigs_this_vb == 1 || // the entire VB is a single contig
                                      est_num_contigs_in_file <  MAX_CONTIGS_IN_FILE); 

        segconf.fasta_has_contigs &= !TXT_DT(GFF); // GFF3-embedded FASTA doesn't have contigs, because did=0 is reserved for GFF's SEQID

        ASSINP0 (!flag.make_reference || segconf.fasta_has_contigs, "Can't use --make-reference on this file, because Genozip can't find the contig names in the FASTA description lines");
    }
}

bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _FASTA_TOPLEVEL ||
           dict_id.num == _FASTA_DESC     ||
           dict_id.num == _FASTA_LINEMETA ||
           dict_id.num == _FASTA_TAXID    ||
           dict_id.num == _FASTA_EOL;
}

// description line - we segment it to its components
// note: we store the DESC container in its own ctx rather than just directly in LINEMETA, to make it easier to grep
static void fasta_seg_desc_line (VBlockFASTAP vb, rom line_start, uint32_t line_len, bool *has_13)
{
    SAFE_NUL (&line_start[line_len]);
    
    // we store the contig name in a dictionary only (no b250), to be used if this fasta is used as a reference
    rom chrom_name = line_start + 1;
    unsigned chrom_name_len = strcspn (line_start + 1, " \t\r\n");

    ASSSEG0 (chrom_name_len, line_start, "contig is missing a name");

    __atomic_add_fetch (&z_file->num_sequences, (uint64_t)1, __ATOMIC_RELAXED);

    if (!flag.make_reference) {
        if (segconf.fasta_has_contigs) 
            tokenizer_seg (VB, CTX(FASTA_DESC), line_start, line_len, sep_with_space, 0);

        // if we don't have contigs, eg this is an amino acid fasta, we're better off
        // not tokenizing the description line as its components are often correlated
        else
            seg_add_to_local_text (VB, CTX(FASTA_DESC), line_start, line_len, LOOKUP_SIMPLE, line_len);

        char special_snip[100]; unsigned special_snip_len = sizeof (special_snip);
        seg_prepare_snip_other_do (SNIP_REDIRECTION, _FASTA_DESC, false, 0, 0, &special_snip[2], &special_snip_len);
        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_DESC;

        seg_by_did (VB, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    // case make_ref: add contig metadata (the rest of the line, except for the chrom_name)
    else {
        rom md_start = chrom_name + chrom_name_len + strspn (&chrom_name[chrom_name_len], " \t");
        unsigned md_len = MIN_(strcspn (md_start, "\n\r"), REFCONTIG_MD_LEN-1);
        memcpy (vb->contig_metadata.str, md_start, md_len);
        vb->has_contig_metadata = true;        
    }

    // add contig to CONTIG dictionary (but not b250) and verify that its unique
    if (segconf.fasta_has_contigs || flag.make_reference || segconf.running) {
        
        if (segconf.running)
            seg_create_rollback_point (VB, NULL, 1, FASTA_CONTIG);
            
        bool is_new;
        chrom_seg_no_b250 (VB, STRa(chrom_name), &is_new);

        if (!is_new && segconf.running) {
            segconf.fasta_has_contigs = false; // this FASTA is not of contigs. It might be just a collection of variants of a sequence.
            seg_rollback (VB);
        }

        else {
            ASSINP (is_new || segconf.running, "Error: bad FASTA file - sequence \"%.*s\" appears more than once%s", chrom_name_len, chrom_name,
                    flag.bind ? " (possibly in another FASTA being bound)" : 
                    (flag.reference & REF_ZIP_LOADED) ? " (possibly the sequence size exceeds vblock size, try enlarging with --vblock)" : "");
         
            vb->ra_initialized = true;
        }
    }

    vb->last_line = FASTA_LINE_DESC;    
    SAFE_RESTORE;
}

static void fast_seg_comment_line (VBlockFASTAP vb, rom line_start, uint32_t line_len, bool *has_13)
{
    if (!flag.make_reference) {
        seg_add_to_local_text (VB, CTX(FASTA_COMMENT), line_start, line_len, LOOKUP_NONE, line_len); 

        char special_snip[100]; unsigned special_snip_len = sizeof (special_snip);
        seg_prepare_snip_other_do (SNIP_OTHER_LOOKUP, _FASTA_COMMENT, false, 0, 0, &special_snip[2], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_COMMENT;

        seg_by_did (VB, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    vb->last_line = FASTA_LINE_COMMENT;
}

// ZIP: main thread during segconf.running
static SeqType fasta_get_seq_type (STRp(seq))
{
    // we determine the type by characters that discriminate between protein and nucleotides, 
    // according to: https://www.bioinformatics.org/sms/iupac.html and https://en.wikipedia.org/wiki/FASTA_format
    // A,C,D,G,H,K,M,N,R,S,T,V,W,Y are can be either nucleoide or protein - in particular all of A,C,G,T,N can
    static bool uniq_amino[256]    = { ['E']=true, ['F']=true, ['I']=true, ['L']=true, ['P']=true, ['Q']=true, 
                                       ['X']=true, ['Z']=true,  // may be protein according to the FASTA_format page
                                       ['e']=true, ['f']=true, ['i']=true, ['l']=true, ['p']=true, ['q']=true, 
                                       ['x']=true, ['z']=true };
                                      
    static bool nuke_or_amino[256] = { ['A']=true, ['C']=true, ['D']=true, ['G']=true, ['H']=true, ['K']=true, ['M']=true, 
                                       ['N']=true, ['R']=true, ['S']=true, ['T']=true, ['V']=true, ['W']=true, ['Y']=true,
                                       ['U']=true, ['B']=true, // may be protein according to the FASTA_format page (in addition to standard nuke IUPACs)
                                       ['a']=true, ['c']=true, ['d']=true, ['g']=true, ['h']=true, ['k']=true, ['m']=true, 
                                       ['n']=true, ['r']=true, ['s']=true, ['t']=true, ['v']=true, ['w']=true, ['y']=true,
                                       ['u']=true, ['b']=true };
    
    bool evidence_of_amino=false, evidence_of_both=false;

    for (uint32_t i=0; i < seq_len; i++) {
        if (uniq_amino[(int)seq[i]])    evidence_of_amino = true;
        if (nuke_or_amino[(int)seq[i]]) evidence_of_both = true;

        segconf.seq_type_counter++;
    }

    if (evidence_of_amino) return SQT_AMINO;

    if (evidence_of_both) return segconf.seq_type_counter > 10000 ? SQT_NUKE : SQT_NUKE_OR_AMINO; // A,C,G,T,N are in both - call it as NUKE if we've seen enough without evidence of unique Amino characters
    
    return segconf.seq_type; // unchanged
}

static void fasta_seg_seq_line_do (VBlockFASTAP vb, uint32_t line_len, bool is_first_line_in_contig)
{
    Context *lm_ctx  = CTX(FASTA_LINEMETA);
    Context *seq_ctx = CTX(FASTA_NONREF);

    // line length is same as previous SEQ line
    if (!is_first_line_in_contig && ctx_has_value_in_line_(vb, lm_ctx) && line_len == lm_ctx->last_value.i) 
        seg_duplicate_last (VB, lm_ctx, 0);

    else { 
        char special_snip[100]; unsigned special_snip_len = sizeof (special_snip);
        seg_prepare_snip_other_do (SNIP_OTHER_LOOKUP, _FASTA_NONREF, 
                                   true, (int32_t)line_len, 0, &special_snip[3], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_SEQ;
        special_snip[2] = '0' + is_first_line_in_contig; 
        seg_by_ctx (VB, special_snip, 3 + special_snip_len, lm_ctx, 0);  // the payload of the special snip, is the OTHER_LOOKUP snip...
    }

    // note: we don't set value for first line, so that seg_duplicate_last doesn't copy it - since special_snip[2] is different 
    if (!is_first_line_in_contig) 
        ctx_set_last_value (VB, lm_ctx, (ValueType){ .i = line_len });

    seq_ctx->txt_len   += line_len;
    seq_ctx->local.len += line_len;
} 

static void fasta_seg_seq_line (VBlockFASTAP vb, STRp(line), 
                                bool is_last_line_vb_no_newline, bool is_last_line_in_contig, 
                                bool has_13)
{
    vb->lines_this_contig++;

    *DATA_LINE (vb->line_i) = (ZipDataLineFASTA){ .seq_data_start = BNUMtxt (line),
                                                  .seq_len        = line_len,
                                                  .has_13         = has_13 };

    if (flag.make_reference)
        CTX(FASTA_NONREF)->local.len += line_len;

    // after last line in contig, we know if its SIMPLE or PBWT and seg all lines in this contig into FASTA_LINEMETA
    if (!flag.make_reference && is_last_line_in_contig) {

        ZipDataLineFASTA *dl = DATA_LINE (vb->line_i - vb->lines_this_contig + 1);
    
        for (int32_t i=0; i < vb->lines_this_contig; i++) {
            fasta_seg_seq_line_do (vb, dl[i].seq_len, i==0);

            // case last Seq line of VB, and this VB ends part-way through the Seq line (no newline)
            if (is_last_line_vb_no_newline && (i == vb->lines_this_contig - 1)) 
                seg_by_did (VB, dl[i].has_13 ? "\r" : "", dl[i].has_13, FASTA_EOL, dl[i].has_13); // an EOL without \n
            else 
                seg_by_did (VB, dl[i].has_13 ? "\r\n" : "\n", 1 + dl[i].has_13, FASTA_EOL, 1 + dl[i].has_13);
        }        
        vb->lines_this_contig = 0;
    }

    if (segconf.running && (segconf.seq_type == SQT_UNKNOWN || segconf.seq_type == SQT_NUKE_OR_AMINO))
        segconf.seq_type = fasta_get_seq_type (STRa(line));

    vb->last_line = FASTA_LINE_SEQ;

    // case: this sequence is continuation from the previous VB - we don't yet know the chrom - we will update it,
    // and increment the min/max_pos relative to the beginning of the seq in the vb later, in random_access_finalize_entries
    if (segconf.fasta_has_contigs && !vb->chrom_name && !vb->ra_initialized) {
        random_access_update_chrom (VB, 0, WORD_INDEX_NONE, 0, 0);
        vb->ra_initialized = true;
        vb->chrom_node_index = WORD_INDEX_NONE; // the chrom started in a previous VB, we don't yet know its index
    }

    if (segconf.fasta_has_contigs)
        random_access_increment_last_pos (VB, 0, line_len); 
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
rom fasta_seg_txt_line (VBlockFASTAP vb, rom line_start, uint32_t remaining_txt_len, bool *has_13) // index in vb->txt_data where this line starts
{
    // get entire line
    unsigned line_len;
    int32_t remaining_vb_txt_len = BAFTtxt - line_start;
    
    rom next_field = seg_get_next_line (vb, line_start, &remaining_vb_txt_len, &line_len, !vb->vb_has_no_newline, has_13, "FASTA line");

    // case: description line - we segment it to its components
    if (*line_start == '>' || (*line_start == ';' && vb->last_line == FASTA_LINE_SEQ)) {
        fasta_seg_desc_line (vb, line_start, line_len, has_13);

        if (kraken_is_loaded) {
            unsigned qname_len = strcspn (line_start + 1, " \t\r\n"); // +1 to skip the '>' or ';'
            kraken_seg_taxid (VB, FASTA_TAXID, line_start + 1, qname_len, true);
        }
    }

    // case: comment line - stored in the comment buffer
    else if (*line_start == ';' || !line_len) 
        fast_seg_comment_line (vb, line_start, line_len, has_13);

    // case: sequence line
    else 
        fasta_seg_seq_line (vb, line_start, line_len, 
                            remaining_txt_len == line_len + *has_13, // true if this is the last line in the VB with no newline (but may or may not have \r)
                            !remaining_vb_txt_len || *next_field == '>' || *next_field == ';' || *next_field == '\n' || *next_field == '\r', // is_last_line_in_contig
                            *has_13);

    return next_field;
}

//-------------------------
// PIZ stuff
//-------------------------

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fasta_piz_is_skip_section)
{
    if (st != SEC_B250 && st != SEC_LOCAL) return false; // we only skip contexts

    if (flag.reading_reference) return false;  // doesn't apply when using FASTA as a reference

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
    if (flag.header_only_fast && 
        (dict_id.num == _FASTA_NONREF || dict_id.num == _FASTA_NONREF_X || dict_id.num == _FASTA_COMMENT))
        return true;

    // when grepping by main thread - skipping all sections but DESC
    if (purpose == SKIP_PURPOSE_PREPROC && 
        dict_id.num != _FASTA_DESC && !dict_id_is_fasta_desc_sf (dict_id))
        return true;

    // no need for the TAXID data if user didn't specify --taxid
    if (flag.kraken_taxid == TAXID_NONE && dict_id.num == _FASTA_TAXID)
        return true;

    return false;
}

// remove trailing newline before SEQ lines in case of --sequential. Note that there might be more than one newline
// in which case the subsequent newlines were part of empty COMMENT lines
static inline void fasta_piz_unreconstruct_trailing_newlines (VBlockFASTAP vb)
{
    while (*BLSTtxt == '\n' || *BLSTtxt == '\r') 
        vb->txt_data.len32--;

    // update final entries vb->lines to reflect the removal of the final newlines
    uint32_t line_index = BAFTtxt - B1STtxt;

    for (int32_t line_i = vb->line_i; line_i >= 0 && *B32(vb->lines, line_i) > line_index; line_i--)
        *B32(vb->lines, line_i) = line_index;
}

// this is used for end-of-lines of a sequence line, that are not the last line of the sequence. we skip reconstructing
// the newline if the user selected --sequential
SPECIAL_RECONSTRUCTOR_DT (fasta_piz_special_SEQ)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    bool is_first_seq_line_in_this_contig = snip[0] - '0';

    // skip showing line if this contig is grepped - but consume it anyway (possibly carried through from previous VB)
    if (vb->contig_grepped_out) vb->drop_curr_line = "grep";

    // --sequential - if this is NOT the first seq line in the contig, we delete the previous end-of-line
    else if (flag.sequential && !is_first_seq_line_in_this_contig) 
        fasta_piz_unreconstruct_trailing_newlines (vb);

    // in case of not showing the SEQ in the entire file - we can skip consuming it
    if (flag.header_only_fast) // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->drop_curr_line = "header_only_fast";     
    else 
        reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, snip+1, snip_len-1, RECON_ON);    

    // case: --sequential, and this seq line is the last line in the vb, and it continues in the next vb
    if (  flag.sequential && // if we are asked for a sequential SEQ
          vb->line_i == vb->lines.len-1 && // and this is the last line in this vb 
          !vb->drop_curr_line && 
          random_access_does_last_chrom_continue_in_next_vb (vb->vblock_i)) // and this sequence continues in the next VB 
        fasta_piz_unreconstruct_trailing_newlines (vb); // then: delete final newline if this VB ends with a newline

    vb->last_line = FASTA_LINE_SEQ;

    return NO_NEW_VALUE;
}

SPECIAL_RECONSTRUCTOR_DT (fasta_piz_special_COMMENT)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;

    // skip showing comment line in case cases - but consume it anyway:
    if (  vb->contig_grepped_out || // 1. if this contig is grepped out
          flag.out_dt == DT_PHYLIP)       // 2. if we're outputting in Phylis format
        vb->drop_curr_line = "grep";

    // in case of not showing the COMMENT in the entire file (--header-only or this is a --reference) - we can skip consuming it
    if (flag.header_only_fast)  // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->drop_curr_line = "header_only_fast";     
    else 
        reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, snip, snip_len, RECON_ON);    

    vb->last_line = FASTA_LINE_COMMENT;

    return NO_NEW_VALUE;
}

bool fasta_piz_init_vb (VBlockP vb, const SectionHeaderVbHeader *header, uint32_t *txt_data_so_far_single_0_increment)
{
    VB_FASTA->contig_grepped_out = writer_get_fasta_contig_grepped_out (vb->vblock_i);
    return true;
}

// PIZ Called sequentially for all VBs from writer. excludes VBs due to --grep, --grep-w or --regions 
bool fasta_piz_is_vb_needed (VBIType vb_i)
{
    if (!flag.grep) return true;

    TEMP_FLAG (show_containers, false); // don't show during preprocessing

    rom task_name = "fasta_filter_grep";

    VBlockP vb = vb_get_vb (POOL_MAIN, task_name, vb_i, COMP_MAIN);
    vb->preprocessing = true;

    piz_read_one_vb (vb, false); 

    // we only need room for one line for now 
    uint32_t longest_line_len = BGEN32 (B1ST (SectionHeaderVbHeader, vb->z_data)->longest_line_len);
    buf_alloc (vb, &vb->txt_data, 0, longest_line_len, char, 1.1, "txt_data");

    // uncompress & map desc field (filtered by piz_is_skip_section)
    piz_uncompress_all_ctxs (VB);

    Context *desc_ctx = CTX(FASTA_DESC);
    desc_ctx->iterator.next_b250 = B1ST8 (desc_ctx->b250); 

    uint32_t num_descs = random_access_num_chroms_start_in_this_vb (vb->vblock_i);
    ASSERT0 (vb_i > 1 || num_descs, "Expecting num_descs>0 in vb_i==1");

    // iterate on all DESCs of this VB
    bool grepped_in = false, last_contig_grepped_in = false;
    for (uint32_t desc_i=0; desc_i < num_descs; desc_i++) { 
        vb->txt_data.len = 0; // reset
        reconstruct_from_ctx (vb, FASTA_DESC, 0, true);

        bool match = flag.grep && piz_grep_match (B1STtxt, BAFTtxt);
        grepped_in |= match;
        if (match && desc_i == num_descs - 1) last_contig_grepped_in = true;
    }

    // last DESC- carry over whether its grepped to the next VB - in case next VB starts not from the description line
    // similarly, note whether the previous VB ended with a grepped sequence. If previous VB didn't have any description
    // i.e the entire VB was a sequence that started in an earlier VB - the grep status of the earlier VB is carried forward

    bool needed = (vb_i > 1 && !writer_get_fasta_contig_grepped_out(vb_i)) || // if the last contig of the previous vb was grepped in - then include this VB anyway as it *might* start with SEQ
                   grepped_in;
 
    // update the next VB - any any initial SEQ data is grepped out under the following conditions
    if (vb_i < z_file->num_vbs)  // not last VB
        if ((num_descs && !last_contig_grepped_in) || // this contig has DESCs, and the final DESC (which may continue in next VB) is grepped out
            (!num_descs && !needed))                  // this contig has no DESC, it just continues SEQ of prev VB, which is grepped out
            writer_set_fasta_contig_grepped_out (vb_i + 1); // if next VB starts with SEQ, that SEQ is part of the grepped-out contig

    vb_release_vb (&vb, task_name);

    RESTORE_FLAG (show_containers);
    return needed;
}

// PIZ: piz_process_recon callback: called by the main thread, in the order of VBs
void fasta_piz_process_recon (VBlockP vb)
{
    if (flag.reading_kraken)
        kraken_piz_handover_data (vb);
}

// Phylip format mandates exactly 10 space-padded characters: http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html
static inline void fasta_piz_translate_desc_to_phylip (VBlockFASTAP vb, char *desc_start)
{
    uint32_t recon_len = BAFTtxt - desc_start;
    *BAFTtxt = 0; // nul-terminate

    rom chrom_name = desc_start + 1;
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");

    memmove (desc_start, chrom_name, MIN_(chrom_name_len, 10));
    if (chrom_name_len < 10) memcpy (desc_start + chrom_name_len, "          ", 10-chrom_name_len); // pad with spaces

    if (recon_len > 10) vb->txt_data.len -= recon_len - 10; // we do it this way to avoid signed problems
    else                vb->txt_data.len += 10 - recon_len;
}

// shorten DESC to the first white space
static inline void fasta_piz_desc_header_one (VBlockFASTAP vb, char *desc_start)
{
    uint32_t recon_len = BAFTtxt - desc_start;
    *BAFTtxt = 0; // nul-terminate
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");
    
    vb->txt_data.len -= recon_len - chrom_name_len -1;
}

SPECIAL_RECONSTRUCTOR_DT (fasta_piz_special_DESC)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;
    vb->contig_grepped_out = false;

    char *desc_start = BAFTtxt;
    reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, snip, snip_len, RECON_ON);    
    *BAFTtxt = 0; // for strstr and strcspn

    // if --grep: here we decide whether to show this contig or not
    if (flag.grep) 
        vb->contig_grepped_out = !strstr (desc_start, flag.grep);
    
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n"); // +1 to skip the '>'
    
    if (CTX(FASTA_CONTIG)->is_loaded) { // some FASTAs may not have contigs
        // --taxid: grep out by Kraken taxid 
        if (flag.kraken_taxid != TAXID_NONE) 
            vb->contig_grepped_out |= 
                ((!kraken_is_loaded && !kraken_is_included_stored (VB, FASTA_TAXID, false)) ||
                ( kraken_is_loaded && !kraken_is_included_loaded (VB, desc_start + 1, chrom_name_len)));
        
        vb->chrom_node_index = vb->last_index(CHROM) = ctx_search_for_word_index (CTX(CHROM), desc_start + 1, chrom_name_len);

        // note: this logic allows the to grep contigs even if --no-header 
        if (vb->contig_grepped_out)
            vb->drop_curr_line = "grep";     
    }
    
    if (flag.no_header)
        vb->drop_curr_line = "no_header";     

    if (flag.out_dt == DT_PHYLIP) 
        fasta_piz_translate_desc_to_phylip (vb, desc_start);

    if (flag.header_one)
        fasta_piz_desc_header_one (vb, desc_start);

    vb->last_line = FASTA_LINE_DESC;

    return NO_NEW_VALUE;
}

//---------------------------------------
// Multifasta -> PHYLIP translation stuff
//---------------------------------------

// create Phylip header line
TXTHEADER_TRANSLATOR (txtheader_fa2phy)
{
    // get length of contigs and error if they are not all the same length
    uint32_t contig_len = random_access_verify_all_contigs_same_length();

    bufprintf (comp_vb, txtheader_buf, "%"PRIu64" %u\n", ZCTX(FASTA_CONTIG)->word_list.len, contig_len);
}

// Translating FASTA->PHYLIP: drop EOL after DESC + don't allow multiple consecutive EOL
TRANSLATOR_FUNC (fasta_piz_fa2phy_EOL)
{
    if (VB_FASTA->last_line == FASTA_LINE_DESC // line is a DESC
    ||  recon == vb->txt_data.data     // initial EOL - not allowed in Phylip)
    ||  recon[-1] == '\n')             // previous item was an EOL - remove this one then
    
        vb->txt_data.len -= recon_len;
    
    return 0;
}

CONTAINER_FILTER_FUNC (fasta_piz_filter)
{
    // currently unused, but specified in FASTA toplevel containers created up to v11, so the function needs to exist
    return true;
}
