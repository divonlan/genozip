// ------------------------------------------------------------------
//   fasta.c
//   Copyright (C) 2020-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "fasta_private.h"
#include "seg.h"
#include "piz.h"
#include "dict_id.h"
#include "random_access.h"
#include "regions.h"
#include "codec.h"
#include "stats.h"
#include "chrom.h"
#include "tokenizer.h"
#include "writer.h" 
#include "zfile.h"
#include "qname.h"

#define dict_id_is_fasta_desc_sf dict_id_is_type_1
#define dict_id_fasta_desc_sf dict_id_type_1

sSTRl(desc_redirect_snip, 24);
sSTRl(comment_redirect_snip, 24);

typedef struct {
    uint32_t seq_data_start, seq_len; // regular fasta and make-reference: start & length within vb->txt_data
    bool has_13;
} ZipDataLineFASTA;
#define DATA_LINE(i) B(ZipDataLineFASTA, VB_FASTA->lines, (i))

unsigned fasta_vb_size (DataType dt) 
{ 
    return (dt == DT_REF && IS_PIZ) ? sizeof (VBlock) : sizeof (VBlockFASTA); 
}

unsigned fasta_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTA); }

// used by ref_make_create_range
void fasta_get_data_line (VBlockP vb, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len)
{
    *seq_data_start = DATA_LINE(line_i)->seq_data_start;
    *seq_len        = DATA_LINE(line_i)->seq_len;
}

//-------------------------
// TXTFILE stuff
//-------------------------

// detect if a generic file is actually a FASTA - 
// we call based on first character being '>' or '@', and subsequent lines being equal-length nuke or amino characters
// note: we accept comment lines starting with ; but every sequence must start with >
bool is_fasta (STRp(txt_data), bool *need_more/*optional*/)
{
    // characters that may appear in a FASTA seqience
    static bool seq_char[256] = { ['E']=true, ['F']=true, ['I']=true, ['L']=true, ['P']=true, ['Q']=true, ['X']=true, ['Z']=true,
                                  ['e']=true, ['f']=true, ['i']=true, ['l']=true, ['p']=true, ['q']=true, ['x']=true, ['z']=true,                                    
                                  ['A']=true, ['C']=true, ['D']=true, ['G']=true, ['H']=true, ['K']=true, ['M']=true, 
                                  ['N']=true, ['R']=true, ['S']=true, ['T']=true, ['V']=true, ['W']=true, ['Y']=true, ['U']=true, ['B']=true,
                                  ['a']=true, ['c']=true, ['d']=true, ['g']=true, ['h']=true, ['k']=true, ['m']=true, 
                                  ['n']=true, ['r']=true, ['s']=true, ['t']=true, ['v']=true, ['w']=true, ['y']=true, ['u']=true, ['b']=true,
                                  ['-']=true, ['*']=true };

    char dc = DC ? DC : '>'; // if called from GENERIC, DC is not set yet - we test only for '>'
    
    if (dc != '>' && dc != '@' && dc != ';') return false;

    // first line must be the sequence description
    if (!txt_data_len || txt_data[0] != dc || !str_is_printable (STRa(txt_data))) 
        return false; // fail fast

    #define NUM_TEST_LINES 10
    str_split_by_lines (txt_data, txt_data_len, NUM_TEST_LINES);

    #define CANT_TELL ({ if (need_more) { \
                            *need_more = true; /* we can't tell yet - need more data */ \
                            return false;  \
                         }\
                         else return true;  /* called from Seg of vb_i=1 - VB is shorter than one whole sequence - that's ok, we can handle it */ \
                       })

    if (!n_lines) CANT_TELL;

    if (lines[0][0] != dc) 
        return false;

    // skip comment and empty lines
    int line_i=1; 
    while (line_i < n_lines && (!line_lens[line_i] || lines[line_i][0] == ';')) line_i++;

    // we need at least one sequence line
    if (n_lines - line_i < 1) CANT_TELL;
    
    // next lines must contain sequence (specifically, not DC)
    if (!line_lens[line_i] || !seq_char[(int)lines[line_i][0]]) 
        return false;
    
    int seq_line_len = line_lens[line_i];

    for (; line_i < n_lines; line_i++) {
        // we arrived at a line beyond the contig - we're done
        if (!line_lens[line_i] || lines[line_i][0] == dc) break;

        // all sequence lines, except for the last of the contig, must be equal length
        if (line_lens[line_i] != seq_line_len && // line is different length than first line of contig
            line_i != n_lines-1 && line_lens[line_i+1] && lines[line_i+1][0] != dc) // and we are sure that it is not the last line of the contig
            return false;

        // entire line must be nukes or aminos or '-' (gap) or '*' (end of sequence)
        rom after_c = lines[line_i] + line_lens[line_i]; 
        for (rom c=lines[line_i]; c < after_c; c++)
            if (!seq_char[(int)*c]) 
                return false; 
    }

    return true;
}

// returns true if txt_data[txt_i] is the end of a FASTA contig (= next char is '>' or end-of-file), false if not, 
// and -1 if more data (lower first_i) is needed 
static inline int fasta_is_end_of_contig (VBlockP vb, uint32_t first_i,
                                          int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    if (txt_i < Ltxt-1)
        return txt[txt_i+1] == DC;

    // if we're at the end of the line, we scan back to the previous \n and check if it NOT the >
    bool newline_run = true;
    for (int32_t i=txt_i-1; i >= first_i; i--) {
        if (txt[i] == '\n' && !newline_run) // ASSUMES NO COMMENT LINES (starting with ;), To do: fix this
            return txt[i+1] != DC; // this row is a sequence row, not a description row
        
        else if (newline_run && txt[i] != '\n' && txt[i] != '\r')
            newline_run = false; // newline run ends when we encounter the first non newline
    }

    return -1; // we need more data (lower first_i)
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i)
{
    ASSERTNOTZERO (Ltxt);

    int32_t last_i = Ltxt-1;
    bool is_entire_vb = (first_i == 0);

    ARRAY (char, txt, vb->txt_data);

    // case: reference file - we allow only one contig (or part of it) per VB - move second contig onwards to next vb
    // (note: first_i=0 when flag.make_reference)
    if (flag.make_reference) {
        bool data_found = false;
        for (uint32_t i=0; i < Ltxt; i++) {
            // just don't allow now-obsolete ';' rather than trying to disentangle comments from descriptions
            ASSINP (txt[i] != ';', "Error: %s contains a ';' character - this is not supported for reference files. Contig descriptions must begin with a >", txt_name);
        
            // if we've encountered a new DESC line after already seeing sequence data, move this DESC line and
            // everything following to the next VB
            if (data_found && txt[i]==DC && txt[i-1]=='\n') 
                return Ltxt - i;

            if (!data_found && (txt[i] != '\n' && txt[i] != '\r')) data_found = true; // anything, except for empty lines, is considered data
        }
    }

    // we move the final partial line to the next vb (unless we are already moving more, due to a --make-reference)
    for (int32_t i=last_i; i >= (int32_t)first_i; i--) {

        if (txt[i] == '\n') {

            // when compressing FASTA with a reference or --multiseq - an "end of line" means "end of contig" - 
            // i.e. one that the next character is >, or it is the end of the file
            // note: when compressing FASTA with a reference (eg long reads stored in a FASTA instead of a FASTQ), 
            // line cannot be too long - they must fit in a VB
            if (segconf.multiseq) {
                int is_end_of_contig = fasta_is_end_of_contig (vb, first_i, i);

                switch (is_end_of_contig) {
                    case true  : break;            // end of contig
                    case false : continue;         // not end of contig
                    default    : goto out_of_data; // need more data (lower first_i)
                }
            }

            // otherwise - tolerate a VB that ends part way through a SEQ
            else if (is_entire_vb && i+1 < Ltxt &&
                     txt[i+1] != ';' && txt[i+1] != '>' && txt[i+1] != '@') { // partial line isn't a Description or a Comment, hence its a Sequence
                ((VBlockFASTAP)vb)->vb_has_no_newline = true;
                return 0;                
            }
                
            last_i = i;
            return Ltxt-1 - i;
        }
    }

out_of_data:
    // case: an entire FASTA VB without newlines (i.e. a very long sequential SEQ) - we accept a VB without newlines and deal with it in Seg
    if (is_entire_vb && !segconf_running) {
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

extern void fastq_bamass_populate (void);

// called by main thread at the beginning of zipping this file
void fasta_zip_initialize (void)
{
    tokenizer_zip_initialize();

    DO_ONCE {
        seg_prepare_snip_other_do (SNIP_REDIRECTION, (DictId)_FASTA_DESC, false, 0, 0, &desc_redirect_snip[2], &desc_redirect_snip_len);
        desc_redirect_snip[0] = SNIP_SPECIAL;
        desc_redirect_snip[1] = FASTA_SPECIAL_DESC;
        desc_redirect_snip_len += 2;

        seg_prepare_snip_other_do (SNIP_OTHER_LOOKUP, (DictId)_FASTA_COMMENT, false, 0, 0, &comment_redirect_snip[2], &comment_redirect_snip_len);
        comment_redirect_snip[0] = SNIP_SPECIAL;
        comment_redirect_snip[1] = FASTA_SPECIAL_COMMENT;
        comment_redirect_snip_len += 2;
    }

    // with REF_EXTERNAL, user is telling as this is FAF.
    // we just copy all reference contigs. this are not needed for uncompression, just for --coverage/--idxstats
    if (IS_REF_EXTERNAL && z_file->num_txts_so_far == 1) // single file, or first of pair (and never Deep)
        ctx_populate_zf_ctx_from_contigs (FASTA_CONTIG, ref_get_ctgs()); 

    if (flag.bam_assist) {
        qname_zip_initialize(); // note: might also be called during segconf from fasta_segconf_is_qualless_fastq. no harm.
        fastq_bamass_populate();
    }
}

void fasta_zip_set_vb_header_specific (VBlockP vb, SectionHeaderVbHeaderP vb_header)
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
        *line_data = dl->seq_len ? Btxt (dl->seq_data_start) : NULL;

    if (is_rev) *is_rev = 0;
}   

// called by main thread, as VBs complete (might be out-of-order)
void fasta_zip_after_compute (VBlockP vb)
{
    z_file->num_sequences += vb->num_sequences; // for stats
}

// we consider it a "QUAL-less FASTQ" if both these conditions are met:
// 1. DESC line starts with a recognized QNAME
// 2. Testing the first bunch of lines - the consist of pairs of (DESC, SEQ) lines (i.e. all contigs are single lines)
static bool fasta_segconf_is_qualless_fastq (VBlockP vb)
{
    #define NUM_FASTQ_TEST_LINES 2000
    str_split_by_lines (vb->txt_data.data, vb->txt_data.len32, NUM_FASTQ_TEST_LINES);

    if (n_lines < 16 && !txt_file->no_more_blocks) 
        return false; // not enough lines to determine (except of segconf is the entire file) 

    n_lines = ROUNDDOWN2(n_lines) - 2; // keep whole pairs or lines, and drop last pair that might be truncated. now there are at least 4 pairs (8 lines)

    for (uint32_t i=0; i < n_lines; i += 2) 
        if (lines[i][0] != DC) 
            return false; // not all contigs are single-line of SEQ

    // test first 4 seqs to make sure they contain only characters valid FASTQ (uppercase and not amino acids)
    for (uint32_t i=1; i <= 7; i += 2) 
        if (!str_is_fastq_seq (lines[i], MIN_(256, line_lens[i]))) 
            return false; // line contains an invalid SEQ character for FASTQ  
    
    qname_zip_initialize(); 
    qname_segconf_discover_flavor (vb, QNAME1, lines[0]+1, strcspn (lines[0]+1, " \t\r\n")); 

    if (!flag.reference  && // if a reference file is given, we treat FASTA and FASTQ regardless of QNAME
        !flag.bam_assist && // --bamass, --pair or --deep: user is effectively telling us this is FAF file
        !flag.pair       && 
        !flag.deep       && 
        !segconf.qname_flavor[QNAME1]) 
        return false;

    // DECIDED! convert VB to a FASTQ VB.

    // re-initialize z_file / txt_file stuff
    txt_file->data_type = DT_FASTQ;
    
    if (!flag.deep) { // if --deep, already initialized as SAM
        z_file->data_type = DT_FASTQ;
        ctx_initialize_predefined_ctxs (z_file->contexts, DT_FASTQ, z_file->d2d_map, &z_file->num_contexts);
    }

    // re-initialize segconf VB
    vb_change_datatype_nonpool_vb (&vb, DT_FASTQ);

    buflist_destroy_private_and_context_vb_bufs (vb);
    ctx_clone (vb); // clone FASTQ contexts

    // move redundant txt to unconsumed (since FASTA allows last line to be a description line or a partial SEQ line)
    rom final_read = memrchr (B1STtxt, DC, Ltxt);
    uint32_t final_read_len = BAFTtxt - final_read;
    if (str_count_char (STRa(final_read), '\n') != 2)
        Ltxt -= final_read_len; // note: this txt_data is discards after segconf so no problem modifying it

    // call FASTQ's seg initialize instead
    FAF = true;
    fastq_seg_initialize (vb);

    return true;
}

void fasta_seg_initialize (VBlockP vb)
{
    START_TIMER;

    // NSBI FASTA download files might have "@" instead of ">"
    if (segconf_running) {
        DC = *B1STtxt; // note: common: '>', very old FASTA: ';', downloaded from NCBI: '@'
        
        ASSINP (is_fasta (STRb(vb->txt_data), NULL), "Error: %s is not a valid FASTA file. Solution: use --input generic", txt_name);
    
        // case: this FASTA looks more like a FASTQ without QUAL data - better off segged as FASTQ with reference support and QNAME handling 
        if (!flag.make_reference && !flag.no_faf && !flag.ext_indexing &&
            fasta_segconf_is_qualless_fastq (vb)) return;
    }

    CTX(FASTA_TOPLEVEL)->no_stons  = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(FASTA_CONTIG)->flags.store = STORE_INDEX; // since v12
    CTX(FASTA_CONTIG)->no_stons    = true; // needs b250 node_index for reference
    CTX(FASTA_LINEMETA)->no_stons  = true; // avoid edge case where entire b250 is moved to local due to singletons, because fasta_reconstruct_vb iterates on ctx->b250
    
    CTX(FASTA_COMMENT)->ltype = LT_STRING;
    
    if (segconf.seq_type == SQT_AMINO || Z_DT(GFF)) 
        CTX(FASTA_DESC)->ltype = LT_STRING;

    // if this neocleotide FASTA of unrelated contigs, we're better off with ACGT        
    if (!segconf_running && !segconf.multiseq && segconf.seq_type == SQT_NUKE)
        codec_acgt_seg_initialize (VB, FASTA_NONREF, true);

    else 
        CTX(FASTA_NONREF)->ltype  = LT_BLOB;

    // in --stats, consolidate stats into FASTA_NONREF
    ctx_consolidate_stats (vb, FASTA_NONREF, FASTA_NONREF_X, DID_EOL);
        
    COPY_TIMER (seg_initialize);
}

void fasta_segconf_finalize (VBlockP vb)
{
    if (segconf.seq_type == SQT_UNKNOWN) segconf.seq_type = SQT_NUKE; // sequence too short to call, and no amino-only characters
    ASSINP0 (!flag.make_reference || segconf.seq_type == SQT_NUKE, "Can't use --make-reference on this file, because it contains amino acids instead of nucleotides");

    // decide whether the sequences in this FASTA represent contigs (in which case we want a FASTA_CONTIG dictionary
    // and random access) or do they represent reads (in which case they are likely too numerous to be added to a dict)
    uint64_t num_contigs_this_vb = CTX(FASTA_CONTIG)->nodes.len;
    ASSINP0 (num_contigs_this_vb, "Invalid FASTA file: no sequence description line");

    uint64_t avg_contig_size_this_vb = Ltxt / num_contigs_this_vb;
    uint64_t est_num_contigs_in_file = txt_file->est_seggable_size / avg_contig_size_this_vb;

    ASSINP0 (!flag.make_reference || segconf.fasta_has_contigs, "Can't use --make-reference on this file, because Genozip can't find the contig names in the FASTA description lines");

    // if we have multiple shortish contigs in segconf data, this might be multiseq - test
    if (num_contigs_this_vb >= 5 && !flag.make_reference && !flag.fast) 
        segconf_test_multiseq (VB, FASTA_NONREF);

    // output FASTA_CONTIG dictionary (allowing genocat --grep and --regions) if requested
    // with --index, or if small enough
    #define MAX_CONTIGS_IN_FILE 10000
    segconf.fasta_has_contigs &= (num_contigs_this_vb == 1 || // the entire VB is a single contig
                                  flag.ext_indexing ||
                                  est_num_contigs_in_file <= MAX_CONTIGS_IN_FILE||
                                  flag.make_reference); 
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
                              { .dict_id = { _FASTA_EOL      }  } }
        };

        container_seg (vb, CTX(FASTA_TOPLEVEL), (ContainerP)&top_level, 0, 0, 0);
    }
}

bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _FASTA_TOPLEVEL ||
           dict_id.num == _FASTA_DESC     ||
           dict_id.num == _FASTA_LINEMETA ||
           dict_id.num == _FASTA_EOL;
}

bool fasta_seg_is_big (ConstVBlockP vb, DictId dict_id, DictId st_dict_id)
{
    return dict_id.num == _FASTA_CONTIG || // some FASTAs have lots of contigs. Since the FASTA data type has very few contexts, we can be generous here.
           dict_id_is_type_1 (dict_id);
}

// description line - we segment it to its components
// note: we store the DESC container in its own ctx rather than just directly in LINEMETA, to make it easier to grep
static void fasta_seg_desc_line (VBlockFASTAP vb, rom line, uint32_t line_len, bool *has_13)
{
    SAFE_NUL (&line[line_len]);
    
    // we store the contig name in a dictionary only (no b250), to be used if this fasta is used as a reference
    rom chrom_name = line + 1;
    unsigned chrom_name_len = strcspn (line + 1, " \t\r\n");

    ASSSEG0 (chrom_name_len, "contig is missing a name");

    vb->num_sequences++; // for stats

    if (!flag.make_reference) {
        if (segconf.seq_type == SQT_NUKE && !Z_DT(GFF)) 
            tokenizer_seg (VB, CTX(FASTA_DESC), STRa(line), sep_with_space, 0);

        // if this is an amino acid fasta, we're better off not tokenizing the description line as its 
        // components are often correlated
        else
            seg_add_to_local_string (VB, CTX(FASTA_DESC), STRa(line), LOOKUP_SIMPLE, line_len);

        seg_by_did (VB, STRa(desc_redirect_snip), FASTA_LINEMETA, 0);
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
    if (segconf.fasta_has_contigs || flag.make_reference || segconf_running) {
        
        if (segconf_running)
            seg_create_rollback_point (VB, NULL, 1, FASTA_CONTIG);
            
        bool is_new;
        chrom_seg_no_b250 (VB, STRa(chrom_name), &is_new);

        if (!is_new && segconf_running) {
            if (flag.ext_indexing) {
                WARN_ONCE ("FYI: Cannot index %s, because it contains duplicate contig names, for example: %.*s", 
                           txt_name, STRf(chrom_name));
                flag.ext_indexing = 0;
            }
            
            segconf.fasta_has_contigs = false; // this FASTA is not of contigs. It might be just a multiseq (collection of variants of a sequence).
            seg_rollback (VB);
        }

        else {
            ASSINP (is_new || segconf_running, "Error: bad FASTA file - sequence \"%.*s\" appears more than once", STRf(chrom_name));
         
            vb->ra_initialized = true;
        }
    }

    vb->last_line = FASTA_LINE_DESC;    
    SAFE_RESTORE;
}

static void fast_seg_comment_line (VBlockFASTAP vb, STRp (line), bool *has_13)
{
    if (!flag.make_reference) {
        seg_add_to_local_string (VB, CTX(FASTA_COMMENT), STRa(line), LOOKUP_NONE, line_len); 

        seg_by_did (VB, STRa(comment_redirect_snip), FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    vb->last_line = FASTA_LINE_COMMENT;
}

// ZIP: main thread during segconf.running
static void fasta_set_seq_type (VBlockFASTAP vb, STRp(seq))
{    
    static uint32_t n_tested;
    if (!vb->line_i) n_tested = 0;
    
    SAFE_NULT(seq);

    // see: https://www.bioinformatics.org/sms/iupac.html and https://en.wikipedia.org/wiki/FASTA_format

    // case: seq is only A,C,G,T,N - clearly nuke (most common case) (note: failing can still be nuke, bc nuke can also contain IUPAC characters)
    if (seq_len > 50 && strspn (seq, "ACGTNacgtn") == seq_len)
        segconf.seq_type = SQT_NUKE;  
    
    // case: seq contains characters unique to amino sequences (i.e. that are not IUPACs)
    else if (!!strpbrk (seq, "EFILPQXZefilpqxz")) 
        segconf.seq_type = SQT_AMINO; 

    // case: no amino-only charater in the first 4K - we assume only ACTGN and IUPAC "bases" (although we didn't test explicitly)
    else if ((n_tested += seq_len) > 4096) 
        segconf.seq_type = SQT_NUKE;  

    // we can't decide yet - wait for next seq
    else {}

    SAFE_RESTORE;
}

static void fasta_seg_seq_line_do (VBlockFASTAP vb, uint32_t line_len, bool is_first_line_in_contig)
{
    Context *lm_ctx  = CTX(FASTA_LINEMETA);
    Context *seq_ctx = CTX(FASTA_NONREF);

    // line length is same as previous SEQ line
    if (!is_first_line_in_contig && ctx_has_value_in_line_(vb, lm_ctx) && line_len == lm_ctx->last_value.i) 
        seg_duplicate_last (VB, lm_ctx, 0);

    else { 
        STRlic (special_snip, 100);
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

    seq_ctx->txt_len     += line_len;
    seq_ctx->local.len32 += line_len;
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

    if (segconf_running && segconf.seq_type == SQT_UNKNOWN)
        fasta_set_seq_type (vb, STRa(line));

    vb->last_line = FASTA_LINE_SEQ;

    // case: this sequence is continuation from the previous VB - we don't yet know the chrom - we will update it,
    // and increment the min/max_pos relative to the beginning of the seq in the vb later, in random_access_finalize_entries
    if (segconf.fasta_has_contigs && !vb->chrom_name && !vb->ra_initialized) {
        random_access_update_chrom (VB, WORD_INDEX_NONE, 0, 0);
        vb->ra_initialized = true;
        vb->chrom_node_index = WORD_INDEX_NONE; // the chrom started in a previous VB, we don't yet know its index
    }

    if (segconf.fasta_has_contigs)
        random_access_increment_last_pos (VB, line_len); 
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
rom fasta_seg_txt_line (VBlockP vb_, rom line, uint32_t remaining_txt_len, bool *has_13) // index in vb->txt_data where this line starts
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;
    
    // get entire line
    unsigned line_len;
    int32_t remaining_vb_txt_len = BAFTtxt - line;
    
    rom next_field = seg_get_next_line (VB, line, &remaining_vb_txt_len, &line_len, !vb->vb_has_no_newline, has_13, "FASTA line");

    // case: description line - we segment it to its components
    if (*line == DC && (vb->last_line == FASTA_LINE_SEQ || vb->last_line == FASTA_LINE_COMMENT))
        fasta_seg_desc_line (vb, line, line_len, has_13);

    // case: comment line - stored in the comment buffer
    else if (*line == ';' || !line_len) 
        fast_seg_comment_line (vb, STRa(line), has_13);

    // case: sequence line
    else 
        fasta_seg_seq_line (vb, STRa(line), 
                            remaining_txt_len == line_len + *has_13, // true if this is the last line in the VB with no newline (but may or may not have \r)
                            !remaining_vb_txt_len || *next_field == DC || *next_field == ';' || *next_field == '\n' || *next_field == '\r', // is_last_line_in_contig
                            *has_13);

    return next_field;
}

//-------------------------
// PIZ stuff
//-------------------------

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fasta_piz_is_skip_section)
{
    if (!ST(B250) && !ST(LOCAL)) return false; // we only skip contexts

    if (flag.reading_reference) return false;  // doesn't apply when using FASTA as a reference

    // note that flags_update_piz_one_z_file rewrites --header-only as flag.header_only_fast
    if (flag.header_only_fast && 
        (dict_id.num == _FASTA_NONREF || dict_id.num == _FASTA_NONREF_X || dict_id.num == _FASTA_COMMENT))
        return true;

    // when grepping by main thread - skipping all sections but DESC
    if (purpose == SKIP_PURPOSE_PREPROC && 
        dict_id.num != _FASTA_DESC && !dict_id_is_fasta_desc_sf (dict_id))
        return true;

    if (dict_id.num == _FASTA_TAXID) return true; // skip TAXID in old files

    return false;
}

// remove trailing newline before SEQ lines in case of --sequential. Note that there might be more than one newline
// in which case the subsequent newlines were part of empty COMMENT lines
static inline void fasta_piz_unreconstruct_trailing_newlines (VBlockFASTAP vb)
{
    while (*BLSTtxt == '\n' || *BLSTtxt == '\r') 
        Ltxt--;

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
    if (flag.header_only_fast) // note that flags_update_piz_one_z_file rewrites --header-only as flag.header_only_fast
        vb->drop_curr_line = "header_only_fast";     
    else 
        reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, snip+1, snip_len-1, RECON_ON, __FUNCLINE);    

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
    if (vb->contig_grepped_out) // 1. if this contig is grepped out
        vb->drop_curr_line = "grep";

    // in case of not showing the COMMENT in the entire file (--header-only or this is a --reference) - we can skip consuming it
    if (flag.header_only_fast)  // note that flags_update_piz_one_z_file rewrites --header-only as flag.header_only_fast
        vb->drop_curr_line = "header_only_fast";     
    else 
        reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, snip, snip_len, RECON_ON, __FUNCLINE);    

    vb->last_line = FASTA_LINE_COMMENT;

    return NO_NEW_VALUE;
}

// main thread: called for each txt file, after reading global area, before reading txt header
bool fasta_piz_initialize (CompIType comp_i)
{
    return true;
}

bool fasta_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header)
{
    VB_FASTA->contig_grepped_out = writer_get_fasta_contig_grepped_out (vb->vblock_i);
    return true;
}

// PIZ Called sequentially for all VBs from writer. excludes VBs due to --grep, --grep-w or --regions 
bool fasta_piz_is_vb_needed (VBIType vb_i)
{
    if (!flag.grep) return true;

    ASSINP0 (z_file->ra_buf.len, "--grep is not supported for this file because it was not indexed during compression (to index, compress with --index)"); // FASTA compressed without contigs

    TEMP_FLAG (show_containers, false); // don't show during preprocessing

    rom task_name = "fasta_filter_grep";

    VBlockP vb = vb_get_vb (POOL_MAIN, task_name, vb_i, COMP_MAIN);
    vb->preprocessing = true;

    piz_read_one_vb (vb, false); // preprocessing - reads only DESC data

    // we only need room for one line for now 
    uint32_t longest_line_len = BGEN32 (B1ST (SectionHeaderVbHeader, vb->z_data)->longest_line_len);
    buf_alloc (vb, &vb->txt_data, 0, longest_line_len, char, 1.1, "txt_data");

    // uncompress & map desc field (filtered by piz_is_skip_section)
    piz_uncompress_all_ctxs (VB, PUR_FASTA_WRITER_INIT);

    ContextP desc_ctx = CTX(FASTA_DESC);
    desc_ctx->iterator.next_b250 = B1ST8 (desc_ctx->b250); 

    uint32_t num_descs = random_access_num_chroms_start_in_this_vb (vb->vblock_i);
    ASSERT0 (vb_i > 1 || num_descs, "Expecting num_descs>0 in vb_i==1");

    // iterate on all DESCs of this VB
    bool grepped_in = false, last_contig_grepped_in = false;
    for (uint32_t desc_i=0; desc_i < num_descs; desc_i++) { 
        Ltxt = 0; // reset
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

// shorten DESC to the first white space
static inline void fasta_piz_desc_header_one (VBlockFASTAP vb, char *desc_start)
{
    uint32_t recon_len = BAFTtxt - desc_start;
    *BAFTtxt = 0; // nul-terminate
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");
    
    Ltxt -= recon_len - chrom_name_len -1;
}

SPECIAL_RECONSTRUCTOR_DT (fasta_piz_special_DESC)
{
    VBlockFASTAP vb = (VBlockFASTAP)vb_;
    vb->contig_grepped_out = false;

    char *desc_start = BAFTtxt;
    reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, STRa(snip), RECON_ON, __FUNCLINE);    
    *BAFTtxt = 0; // for strstr and strcspn

    // if --grep: here we decide whether to show this contig or not
    if (flag.grep) 
        vb->contig_grepped_out = !strstr (desc_start, flag.grep);
    
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n"); // +1 to skip the '>'
    
    if (CTX(FASTA_CONTIG)->is_loaded) { // some FASTAs may not have contigs
        vb->chrom_node_index = vb->last_index(CHROM) = ctx_search_for_word_index (CTX(CHROM), desc_start + 1, chrom_name_len);

        // note: this logic allows the to grep contigs even if --no-header 
        if (vb->contig_grepped_out)
            vb->drop_curr_line = "grep";     
    }
    
    if (flag.no_header)
        vb->drop_curr_line = "no_header";     

    if (flag.header_one)
        fasta_piz_desc_header_one (vb, desc_start);

    vb->last_line = FASTA_LINE_DESC;

    return NO_NEW_VALUE;
}

CONTAINER_FILTER_FUNC (fasta_piz_filter)
{
    // currently unused, but specified in FASTA toplevel containers created up to v11, so the function needs to exist
    return true;
}
