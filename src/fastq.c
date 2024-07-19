// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <dirent.h>
#include <libgen.h>
#include "fastq_private.h"
#include "codec.h"
#include "writer.h"
#include "bases_filter.h"
#include "qname.h"
#include "tip.h"
#include "deep.h"
#include "tokenizer.h"
#include "coverage.h"
#include "arch.h"
#include "filename.h"
#include "aligner.h"
#include "zfile.h"
#include "zriter.h"
#include "qname_filter.h"
#include "mgzip.h"

#define dict_id_is_fastq_qname_sf dict_id_is_type_1
#define dict_id_is_fastq_aux      dict_id_is_type_2
#define dict_id_fastq_qname_sf dict_id_type_1

sSTRl (copy_qname_snip, 30);

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

#define DT_NAME (FAF ? "FASTA" : "FASTQ")

//-----------------------
// TXTFILE stuff
//-----------------------

static inline bool is_valid_read (const rom *t,      // 4 (FASTQ) or 2 (FASTA-as-FASTQ) textual lines
                                  const uint32_t *l) // their lengths
{
    if (FAF)
        return l[0] >= 2 && t[0][0] == DC && // DESC line starts with a '> or '@' and is of length at least 2
               l[1] >= 1 && str_is_fastq_seq (t[1], MIN_(256,l[1])); // first 256 chars of SEQ consists of only valid characters
    
    else
        return l[0] >= 2 && t[0][0] == '@' && // DESC line starts with a '@' and is of length at least 2
               l[1] >= 1 && l[1] == l[3]   && // QUAL and SEQ have the same length, which is at least 1
               l[2] >= 1 && t[2][0] == '+';   // THIRD line starts with '+'
}

// detect if a GENERIC file is actually a FASTQ - we call based on first character being '@', line 3 starting with '+', and line 1 and 3 being the same length
bool is_fastq (STRp(header), bool *need_more)
{
    if (!header_len || header[0] != (DC ? DC : '@') || !str_is_printable (STRa(header))) return false; // fail fast

    #define NUM_TEST_READS 3
    str_split_by_lines (header, header_len, 4 * NUM_TEST_READS);
    n_lines = ROUNDDOWN4 (n_lines); // round to whole reads

    if (!n_lines) {
        if (need_more) *need_more = true; // we can't tell yet - need more data
        return false;
    }

    for (int i=0; i < n_lines; i += 4)
        if (!is_valid_read (&lines[i], &line_lens[i])) return false; // SEQ and QUAL lines are of equal length

    return true;
}

VBIType fastq_get_R1_vb_i (VBlockP vb)          { return VB_FASTQ->R1_vb_i; }
uint32_t fastq_get_R1_num_lines (VBlockP vb)    { return VB_FASTQ->R1_num_lines; }
rom fastq_get_R1_last_qname (VBlockP vb)        { return VB_FASTQ->R1_last_qname; }
bool is_fastq_pair_2 (VBlockP vb)               { return VB_DT(FASTQ) && VB_FASTQ->R1_vb_i > 0; }

uint32_t fastq_get_R1_txt_data_len (VBlockP vb) 
{ 
    if (!VB_FASTQ->R1_vb_i) return 0; // no corresponding R1 VB (note: we don't error here, so txtfile_read_vblock can verify that indeed R2 has no more data)

    ASSERT (z_file->R1_txt_data_lens.len32 > VB_FASTQ->R1_vb_i - z_file->R1_first_vb_i, 
            "%s: expecting z_file->R1_txt_data_lens.len=%u > R1_vb_i=%u - R1_first_vb_i=%u", 
            VB_NAME, z_file->R1_txt_data_lens.len32, VB_FASTQ->R1_vb_i, z_file->R1_first_vb_i);

    return *B32(z_file->R1_txt_data_lens, VB_FASTQ->R1_vb_i - z_file->R1_first_vb_i); 
}

// "pair assisted" is a type pairing in which R1 data is loaded to ctx->localR1/b250R1 and R2 consults with it in seg/recon.
// An R2 section that is pair-assisted will have FlagsCtx.paired set (note: in v14, due to a bug, R2.SQBITMAP.b250 should have but fails to set this flag)
#define DNUM(x) (dict_id.num == _FASTQ_##x)
bool fastq_zip_use_pair_assisted (DictId dict_id, SectionType st)
{
    if (ST(LOCAL) && (DNUM(GPOS) || DNUM(STRAND) || DNUM(DEEP)))
        return true;

    else if (ST(B250) && DNUM(SQBITMAP))
        return true;

    else
        return false;
}

// "pair identical" is a type pairing in which the b250 and/or local buffer is identical in both pairs, 
// so we store only R1 and in PIZ, load it to R2 too. in ZIP, R1 data is loaded to ctx->b250R1/localR1
// An R1 section with FlagsCtx.paired set means that it should be used for R2 too, if the corresponding R2 section is missing.
bool fastq_zip_use_pair_identical (DictId dict_id)
{
    return dict_id_is_fastq_qname_sf (dict_id) || dict_id_is_fastq_aux (dict_id) || 
           DNUM(QNAME) || DNUM(QNAME2) || DNUM(LINE3) || DNUM(EXTRA) ||
           DNUM(AUX) || DNUM(E1L) || DNUM (E2L) || DNUM(TOPLEVEL);
}

// two reads are "interleaved" if their line1 is of identical length, and differs in one character
static bool fastq_zip_is_interleaved (STRp(r1), STRp(r2))
{
    if (segconf.qname_flavor[QNAME1]) {
        r1_len = strcspn (r1, " \t\n\r"); // note: no need to nul-terminate - we always call this function when we know there are full lines, so at least a \n is there
        qname_canonize (QNAME1, r1, &r1_len);

        r2_len = strcspn (r2, " \t\n\r");
        qname_canonize (QNAME1, r2, &r2_len);
    
        return str_issame (r1, r2);
    }

    else {
       if (r1_len != r2_len) return false;
    
       uint32_t count_diff=0, diff_i[3]={};

        // we allow up to two diffs, with the same values, e.g.:
        // A00311:85:HYGWAVAXX:1:1101:3025:1000/1 1:N:0:CAACGAGAGC+GAATTGAGTG
        // A00311:85:HYGWAVAXX:1:1101:3025:1000/2 2:N:0:CAACGAGAGC+GAATTGAGTG
        for (uint32_t i=0; i < r1_len && count_diff < 3; i++)
            if (r1[i] != r2[i]) 
                diff_i[count_diff++] = i;

        if (count_diff <= 2 && !segconf.interleaved_r1) { // segconf: first pair of lines
            segconf.interleaved_r1 = r1[diff_i[0]];
            segconf.interleaved_r2 = r2[diff_i[0]];
        }

        return count_diff <= 2 && 
            (count_diff < 1 || (segconf.interleaved_r1 == r1[diff_i[0]] && segconf.interleaved_r2 == r2[diff_i[0]])) &&
            (count_diff < 2 || (segconf.interleaved_r1 == r1[diff_i[1]] && segconf.interleaved_r2 == r2[diff_i[1]]));
    }
}

static inline bool is_last_qname (VBlockFASTQP vb, STRp(qname), STRp(pair_qname))
{
    qname_canonize (QNAME1, qSTRa(qname)); // changes qname_len, but not qname

    if (qname_len != pair_qname_len) return false;

    // compare qnames in reverse - fail faster
    for (int i=qname_len-1; i >= 0; i--)
        if (qname[i] != pair_qname[i])
            return false;

    if (flag.show_bgzf)
        iprintf ("R2_SYNC_QNAME vb=%-7s R1_vb_i=%u qname=\"%.*s\"\n", VB_NAME, vb->R1_vb_i, STRf(qname));

    return true;
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fastq_unconsumed (VBlockP vb_, 
                          uint32_t first_i) // the smallest index in txt_data for which txt_data is populated (the rest might still in uncompressed MGZIP blocks) 
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;
    ASSERTNOTZERO (Ltxt);

    // if entire R2 vb->txt_data doesn't have a counterpart in R2, truncate it if we are allowed
    if (IS_R2 && !vb->R1_vb_i && flag.truncate)
        return -2;

    // initialize new R2 VB
    if (IS_R2 && !vb->R1_last_qname) {
        ASSINP (vb->R1_vb_i, NO_PAIR_FMT_PREFIX "%s doesn't have a counterpart VB in R1)%s", txt_name, VB_NAME, NO_PAIR_FMT_SUFFIX);

        ASSERT (z_file->R1_last_qname_index.len32 > VB_FASTQ->R1_vb_i - z_file->R1_first_vb_i,
                "%s: z_file->R1_last_qname_index is missing the last_qname of R1_vb_i=%u (R1_first_vb_i=%u len=%u)", 
                VB_NAME, vb->R1_vb_i, z_file->R1_first_vb_i, z_file->R1_last_qname_index.len32);

        uint64_t index = *B64 (z_file->R1_last_qname_index, vb->R1_vb_i - z_file->R1_first_vb_i); // set in fastq_zip_after_compute
        vb->R1_last_qname = Bc(z_file->R1_last_qname, index); // nul-terminated
        vb->R1_last_qname_len = strlen (vb->R1_last_qname);

        vb->R2_lowest_read = vb->R2_highest_read = -1;
    }

    // search backwards a suffient number of newlines (eg. for normal FASTQ: best case: \nD\nS\nT\nQ\n ; worst case: \nD1\nS1\nT1\nQ1\nD2\nS2\nT2\nq2 (q2 is partial Q2))
    int n=0;
    int height = (FAF ? 2 : 4);    // number of lines per read
    int min_lines = height + (segconf.is_interleaved ? height : 0); // minimum lines needed for testing
    int max_lines = min_lines + height + (segconf.is_interleaved ? height : 0); // maximum lines needed for testing (added lines in case of final partial read/interleaved-double-read that needs to be skipped)

    rom lines[9]={};          // newline pointers: nl[0] is the first from the end 
    uint32_t line_lens[9]={}; // lengths of segments excluding \n and \r
    int line_1_modulo = -1;   // value 0-3 means we know which n%4 an R2 read starts, otherwise -1.
    uint32_t highest_read_in_this_call=0;

    for (rom c = ((IS_R2 && vb->R2_lowest_read >= 0) ? Btxt(vb->R2_lowest_read-1) : BLSTtxt), one_before = Btxt (first_i)-1; 
         c >= one_before && (n <= max_lines || IS_R2); 
         c--) 
        
        if (c == one_before || *c == '\n') { // we consider character before the start to also be a "virtual newline"
            memmove (lines+1,     lines,     min_lines * sizeof(rom)); // [0] is always the current, i.e. the lowest in txt_data, line.
            memmove (line_lens+1, line_lens, min_lines * sizeof(uint32_t));
            lines[0] = c+1; // first character after \n
            line_lens[0] = n ? (lines[1] - lines[0] -1/*\n*/ - (lines[1][-2] == '\r')) : 0; // when n=0, it is the final, partial, line, so we considers its length to be 0.
            
            // case: test for valid read after reading a sufficient number of lines
            if (n >= min_lines && 
                (line_1_modulo == -1 || line_1_modulo == n % 4) && // R2: no need to call is_valid_read if we are not on a start of a read
                is_valid_read (lines, line_lens) &&
                (segconf.is_interleaved ? is_valid_read (&lines[height], &line_lens[height]) : true) &&
                (segconf.is_interleaved ? fastq_zip_is_interleaved (STRi(line, 0), STRi(line, height)) : true))
            {
                // case R2: starting with the last validated read, scan backwards until reaching (or not) the read we're looking for
                if (IS_R2) {
                    uint32_t read_bnum = BNUMtxt (lines[0]); 
                    
                    if (line_1_modulo == -1) {
                        highest_read_in_this_call = read_bnum;
                        line_1_modulo = n % 4;
                    }

                    if (vb->R2_lowest_read == -1/*uninitialized*/ || vb->R2_lowest_read > read_bnum) 
                        vb->R2_lowest_read = read_bnum; // this is the read with the lowest index in txt_data so far to be considered

                    // case: in previous calls to this function, we already tested this read and all lower reads - we need to read more data
                    if (vb->R2_highest_read == read_bnum && vb->R2_lowest_read == 0) {
                        vb->R2_highest_read = highest_read_in_this_call; // the highest read we've considered
                        return -1; // all current txt_data has been considered and matching QNAME not found, read more data from disk please
                    }

                    // case: the read does not have the same QNAME as the last read of R1 - continue searching backwards
                    if (!is_last_qname (vb, lines[0]+1/*skip @*/, strcspn (lines[0]+1, " \t\n\r"), STRa(vb->R1_last_qname)))
                        goto next_line;
                }
                
                // everything after the last full read goes to the next VB
                ASSERTNOTNULL (lines[min_lines]);
                return BAFTtxt - lines[min_lines];   // number of "unconsumed" characters remaining in txt_data after the last line of this read
            }
            
            next_line: n++;
        }

    ASSINP (n < max_lines || IS_R2, "%s: Examined %d textual lines at the end of the VB and could not find a valid read, it appears that this is not a valid %s file. Last %u lines examined:\n[0]=\"%.*s\"\n[1]=\"%.*s\"\n[2]=\"%.*s\"\n[3]=\"%.*s\"\n[4]=\"%.*s\"\n",
            VB_NAME, n, DT_NAME, MIN_(n, 5), STRfi(line,0), STRfi(line,1), STRfi(line,2), STRfi(line,3), STRfi(line,4));

    return -1; // uncompress one more mgzip block please 
}

void fastq_zip_set_txt_header_flags (struct FlagsTxtHeader *f)
{
    f->pair = flag.pair;
}

//---------------
// ZIP / SEG stuff
//---------------

// ZIP main thread: called after reading VB
void fastq_zip_init_vb (VBlockP vb)
{
    // in case we're optimizing DESC in FASTQ, we need to know the number of lines
    if (flag.zip_lines_counted_at_init_vb) {
        uint32_t num_lines = str_count_char (STRb(vb->txt_data), '\n');
        ASSERT (num_lines % 4 == 0, "expecting number of txt lines in VB to be a multiple of 4, but found %u", num_lines);

        VB_FASTQ->first_line = txt_file->num_lines + 1;
        txt_file->num_lines += num_lines / 4;     // update here instead of in zip_update_txt_counters;
    }
}

// called by main thread, as VBs complete (might be out-of-order)
void fastq_zip_after_compute (VBlockP vb)
{
    if (flag.deep)
        fastq_deep_zip_after_compute (VB_FASTQ);

    z_file->num_perfect_matches += vb->num_perfect_matches; // for stats
    z_file->num_aligned         += vb->num_aligned;

    if (!flag.deep && IS_REF_LOADED_ZIP) {
        DO_ONCE ref_verify_organism (vb);
    }

    // capture the last qname of each R1 VB, allowing the generation of the respective R2 VB in fastq_unconsumed
    if (IS_R1) {
        // note: we store qnames with indirection (index) because this function is called out-of-order
        buf_alloc_zero (evb, &z_file->R1_last_qname_index, 0, vb->vblock_i - z_file->R1_first_vb_i + 1, uint64_t, CTX_GROWTH, NULL); // pre-allocated in fastq_zip_after_segconf
        *B64(z_file->R1_last_qname_index, vb->vblock_i - z_file->R1_first_vb_i) = z_file->R1_last_qname.len;
        z_file->R1_last_qname_index.len32 = MAX_(z_file->R1_last_qname_index.len32, vb->vblock_i - z_file->R1_first_vb_i + 1);

        STRlast (qname, FASTQ_QNAME);
        qname_canonize (QNAME1, qname, &qname_len); // get canonical qname_len

        buf_add_more (evb, &z_file->R1_last_qname, qname, qname_len, NULL); // pre-allocated in fastq_zip_after_segconf
        BNXTc (z_file->R1_last_qname) = 0; // nul

        if (flag.show_bgzf)
            iprintf ("R1_LAST_QNAME vb=%-7s qname=\"%.*s\" num_lines=%u\n", VB_NAME, STRf(qname), vb->lines.len32);
    }
}

// case of --optimize-DESC: generate the prefix of the read name from the txt file name
// eg. "../../fqs/sample.1.fq.gz" -> "@sample-1."
static void fastq_get_optimized_qname_read_name (void)
{
    strncpy (segconf.optimized_qname.s, txt_file->basename, sizeof(StrText)-2); // qname prefix is truncated at 78 characters
    segconf.optimized_qname_len = 1 + file_get_raw_name_and_type (segconf.optimized_qname.s, NULL, NULL); // remove file type extension

    // replace '.' in the filename with '-' as '.' is a separator in con_genozip_opt
    str_replace_letter (segconf.optimized_qname.s, segconf.optimized_qname_len, '.', '-');

    segconf.optimized_qname.s[segconf.optimized_qname_len-1] = '.';
    memset (&segconf.optimized_qname.s[segconf.optimized_qname_len], 0, sizeof(StrText) - segconf.optimized_qname_len); // hygine
}

// test if an R2 file has been compressed after R1 without --pair - and produce tip if so
static void fastq_tip_if_should_be_pair (void)
{
    if (flag.pair || flag.deep || !txt_file->name || !flag.reference || txt_file->redirected || txt_file->is_remote) 
        return;

    char dir_name_buf[strlen (txt_file->name) + 1];
    strcpy (dir_name_buf, txt_file->name);    
    rom dir_name = dirname (dir_name_buf); // destructive - might replace the last slash with 0

    bool is_cd = dir_name[0] == '.' && dir_name[1] == 0;

    DIR *dir = opendir (dir_name);
    if (!dir) return; // can't open directory

    rom bn = txt_file->basename; // basename
    int bn_len = strlen (bn);

    // find R2 file assuming we are R1, if there is one
    struct dirent *ent=0;
    while (!segconf.r1_or_r2 && (ent = readdir(dir))) 
        segconf.r1_or_r2 = filename_is_fastq_pair (STRa(bn), ent->d_name, strlen (ent->d_name));

    if (segconf.r1_or_r2 == PAIR_R1) 
        TIP ("Using --pair to compress paired-end FASTQs can reduce the compressed file's size by 10%%. E.g.:\n"
             "%s --reference %s --pair %s %s%s%s\n",
             arch_get_argv0(), ref_get_filename (gref), txt_name, (is_cd ? "" : dir_name), (is_cd ? "" : "/"), ent->d_name);

    else if (segconf.r1_or_r2 == PAIR_R2) 
        TIP ("Using --pair to compress paired-end FASTQs can reduce the compressed file's size by 10%%. E.g.:\n"
             "%s --reference %s --pair %s%s%s %s\n",
             arch_get_argv0(), ref_get_filename (gref), (is_cd ? "" : dir_name), (is_cd ? "" : "/"), ent->d_name, txt_name);

    closedir(dir);    
}

// called by main thread at the beginning of zipping this txt file  
void fastq_zip_initialize (void)
{
    DO_ONCE {
        seg_prepare_snip_other (SNIP_COPY, _FASTQ_QNAME, 0, 0, copy_qname_snip);
    }
    
    // reset lcodec for STRAND and GPOS, as these may change between PAIR_R1 and PAIR_R2 files
    ZCTX(FASTQ_STRAND)->lcodec = CODEC_UNKNOWN;
    ZCTX(FASTQ_GPOS  )->lcodec = CODEC_UNKNOWN;

    if (IS_R1)
        z_file->R1_first_vb_i = z_file->num_vbs + 1; // 1 for --pair, >1 for --deep

    // with REF_EXTERNAL, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (IS_REF_EXTERNAL && z_file->num_txts_so_far == 1) // single file, or first of pair (and never Deep)
        ctx_populate_zf_ctx_from_contigs (gref, FASTQ_CONTIG, ref_get_ctgs (gref)); 

    qname_zip_initialize();

    if (flag.deep) 
        fastq_deep_zip_initialize();
}

// called by main thread after each txt file compressing is done
void fastq_zip_finalize (bool is_last_user_txt_file)
{
    if (is_last_user_txt_file && flag.deep)
        fastq_deep_zip_finalize();

    if (!flag.let_OS_cleanup_on_exit && flag.pair != PAIR_R1 && !flag.deep) {
        if (IS_REF_EXT_STORE)
            ref_destroy_reference (gref);
    }
}

// called by Compute thread at the beginning of this VB
void fastq_seg_initialize (VBlockP vb_)
{
    START_TIMER;

    VBlockFASTQP vb = (VBlockFASTQP)vb_;
    declare_seq_contexts;

    if (segconf.running) {
        if (!FAF) DC = '@'; // might be '>' or ';' if fasta_as_fastq (set by fasta_seg_initialize)
        
        // if no --pair, segconf to determine if file is interleaved
        segconf.is_interleaved = (flag.pair || flag.no_interleaved) ? no : unknown; 
    }

    // if this is an R2 VB that has been uncompressed in the compute thread, verify the number lines
    if (IS_R2 && TXT_IS_IN_SYNC) {
        uint32_t actual_num_lines = str_count_char (STRb(vb->txt_data), '\n') / 4;
        ASSERT (actual_num_lines == vb->lines.len32/*set in seg_all_data_lines*/, 
                "expecting n_reads=%u in %s to match n_reads=%u in corresponding R1 vb=%u. effective_codec=%s. Please report this to "EMAIL_SUPPORT". Solution: use --no-bgzf.", 
                actual_num_lines, VB_NAME, vb->lines.len32, vb->R1_vb_i, codec_name (txt_file->effective_codec)); 
    }

    vb->has_extra = segconf.has_extra; // VB-private copy

    if (!flag.deep)
        CTX(FASTQ_CONTIG)->flags.store = STORE_INDEX; // since v12

    if (flag.aligner_available) {
        strand_ctx->ltype        = LT_BITMAP;
        gpos_ctx->ltype          = LT_UINT32;
        gpos_ctx->flags.store    = STORE_INT;
        gpos_d_ctx->ltype        = LT_INT16;
        nonref_ctx->no_callback  = true; // when using aligner, nonref data is in local. without aligner, the entire SEQ is segged into nonref using a callback
        
        bitmap_ctx->ltype        = LT_BITMAP; // implies no_stons
        bitmap_ctx->local_always = true;
        
        buf_alloc (vb, &bitmap_ctx->local, 1, Ltxt / 4, uint8_t, 0, CTX_TAG_LOCAL); 

        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, 1, Ltxt / 128, char, 0, CTX_TAG_LOCAL); 

        if (segconf.is_interleaved) {
            strand_r2_ctx->ltype     = LT_BITMAP;
            gpos_r2_ctx->ltype       = LT_UINT32;            
            gpos_r2_ctx->flags.store = STORE_INT;
            
            buf_alloc (vb, &gpos_ctx->local,    1, vb->lines.len / 2, uint32_t, CTX_GROWTH, CTX_TAG_LOCAL); 
            buf_alloc (vb, &gpos_r2_ctx->local, 1, vb->lines.len / 2, uint32_t, CTX_GROWTH, CTX_TAG_LOCAL); 
            buf_alloc (vb, &gpos_d_ctx->local,  1, vb->lines.len / 2, int16_t,  CTX_GROWTH, CTX_TAG_LOCAL); 

            buf_alloc (vb, &strand_ctx->local,    0, roundup_bits2bytes64 (vb->lines.len / 2), uint8_t, 0, CTX_TAG_LOCAL); 
            buf_alloc (vb, &strand_r2_ctx->local, 0, roundup_bits2bytes64 (vb->lines.len / 2), uint8_t, 0, CTX_TAG_LOCAL); 
        }
        else {
            buf_alloc (vb, &gpos_ctx->local, 1, vb->lines.len, uint32_t, CTX_GROWTH, CTX_TAG_LOCAL); 
            buf_alloc (vb, &strand_ctx->local, 0, roundup_bits2bytes64 (vb->lines.len), uint8_t, 0, CTX_TAG_LOCAL); 

            if (vb->R1_vb_i)
                buf_alloc (vb, &gpos_d_ctx->local,  1, vb->lines.len, int16_t, CTX_GROWTH, CTX_TAG_LOCAL); 
        }
    }

    ctx_set_ltype (VB, LT_UINT8, FASTQ_SEQMIS_A, FASTQ_SEQMIS_C, FASTQ_SEQMIS_G, FASTQ_SEQMIS_T, DID_EOL);

    if (!segconf.multiseq && !segconf.running)
        codec_acgt_seg_initialize (VB, FASTQ_NONREF, true);
    else
        CTX(FASTQ_NONREF)->ltype = LT_BLOB;

    // initialize QUAL to LT_BLOB, it might be changed later to LT_CODEC (eg domq, longr)
    ctx_set_ltype (VB, LT_BLOB, FASTQ_QUAL, DID_EOL);

    if (IS_R1) 
        // cannot all_the_same with no b250 for PAIR_R1 - SQBITMAP.b250 is tested in fastq_get_pair_1_gpos_strand
        // See defect 2023-02-11. We rely on this "no_drop_b250" in fastq_piz_get_r2_is_forward 
        bitmap_ctx->no_drop_b250 = true; 
    
    else if (IS_R2) {
        ASSINP (vb->lines.len32 == vb->R1_num_lines, NO_PAIR_FMT_PREFIX "in vb=%s: lines.len=%u but R1_num_lines=%u in its corresponding R1 vb_i=%u)%s", 
                txt_name, VB_NAME, vb->lines.len32, vb->R1_num_lines, vb->R1_vb_i, NO_PAIR_FMT_SUFFIX);

        // we're pair-2, decompress all of pair-1's contexts needed for pairing
        piz_uncompress_all_ctxs (VB);

        // we've finished uncompressing the pair sections in z_data into their contexts. now, reset z_data for our compressed output coming next.
        vb->z_data.len32 = 0; 

        // always write the R2 pair-assisted LOCAL/B250 sections, as they carry flag.paired which causes PIZ to load of R1 data
        bitmap_ctx->no_drop_b250 = true;
        gpos_ctx->local_always = strand_ctx->local_always = true;
        if (flag.deep) CTX(FASTQ_DEEP)->local_always = true; 
    }

    if (!segconf.running) { // note: if segconf.running, we initialize in qname_segconf_discover_flavor
        qname_seg_initialize (VB, QNAME1, FASTQ_QNAME); 
        qname_seg_initialize (VB, QNAME2, (segconf.desc_is_l3 ? FASTQ_LINE3 : FASTQ_QNAME)); 
        qname_seg_initialize (VB, QLINE3, FASTQ_LINE3); 
        ctx_consolidate_stats (VB, (segconf.desc_is_l3 ? FASTQ_LINE3 : FASTQ_QNAME), FASTQ_AUX, FASTQ_EXTRA, DID_EOL);
    }

    if (IS_R2)
        ctx_create_node (VB, FASTQ_SQBITMAP, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_mate_lookup }, 2);

    // when pairing, we cannot have singletons, bc a singleton in R1, when appearing in R2 will not
    // be a singleton (since not the first appearance), causing both b250 and local to differ between the R's. 
    if (flag.pair)
        for (Did did_i=0; did_i < MAX_DICTS; did_i++)
            if (segconf.has[did_i]) // aux and saux contexts. note: for non-predefined contexts, their did_i will be the same in all VBs, bc the format of all lines in FASTQ is the same. Even if not, all non-predefined contexts in FASTQ are DTYPE_2, and we are setting no_stons for DTYPE_2 contexts anyway
                CTX(did_i)->no_stons = true;

    if (flag.deep)
        fastq_deep_seg_initialize (vb);

    if (segconf.seq_len_dict_id.num)
        ctx_set_store (VB, STORE_INT, ECTX(segconf.seq_len_dict_id)->did_i, DID_EOL);

    if (segconf.has_agent_trimmer) 
        agilent_seg_initialize (VB);
    
    // consolidate stats for --stats
    ctx_consolidate_stats (VB, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, 
                           FASTQ_GPOS, FASTQ_GPOS_R2, FASTQ_GPOS_DELTA, 
                           FASTQ_STRAND, FASTQ_STRAND_R2, 
                           FASTQ_SEQMIS_A, FASTQ_SEQMIS_C, FASTQ_SEQMIS_G, FASTQ_SEQMIS_T, DID_EOL);
    
    ctx_consolidate_stats (VB, FASTQ_QUAL, FASTQ_DOMQRUNS, FASTQ_QUALMPLX, FASTQ_DIVRQUAL, DID_EOL);

    COPY_TIMER (seg_initialize);
}

void fastq_segconf_finalize (VBlockP vb)
{
    segconf.longest_seq_len = vb->longest_seq_len;
    segconf.is_long_reads = segconf_is_long_reads();
    
    if (segconf.is_interleaved == unknown || segconf.is_long_reads || TECH(ULTIMA)) 
        segconf.is_interleaved = no; // no conclusive evidence of interleaving

    if (flag.deep) fastq_deep_seg_finalize_segconf (vb->lines.len32);

    // if no reference, if fastq is multiseq
    if (!flag.reference && !flag.fast) 
        segconf_test_multiseq (VB, FASTQ_NONREF);

    if (!segconf.multiseq) {
        if (!flag.reference && !txt_file->redirected && !flag.seg_only)
            TIP ("Compressing a %s file using a reference file can reduce the compressed file's size by an additional %s.\n"
                 "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n", 
                 DT_NAME, (FAF ? "60%-80%" : "20%-60%"), arch_get_argv0(), txt_file->name);

        else
            fastq_tip_if_should_be_pair();          
    }  

    // cases where aligner is available (note: called even if reference is not loaded, so that it errors in segconf_calculate)
    if (flag.best || flag.deep)
        flag.aligner_available = true;

    // cases where aligner is not available, despite setting in main_load_reference
    if (segconf.is_long_reads) 
        flag.aligner_available = false;

    if (codec_pacb_maybe_used (FASTQ_QUAL)) 
        codec_pacb_segconf_finalize (vb);

    if (codec_longr_maybe_used (vb, FASTQ_QUAL)) 
        codec_longr_segconf_calculate_bins (vb, CTX(FASTQ_QUAL + 1), fastq_zip_qual);

    // note: we calculate the smux stdv to be reported in stats, even if SMUX is not used
    codec_smux_calc_stats (vb);
    
    qname_segconf_finalize (vb);

    // set optimizations (these might tell get canceled in segconf_finalize_optimize())
    if (flag.optimize) {
        segconf.optimize[FASTQ_QNAME]  = !flag.deep; 
        segconf.optimize[FASTQ_QNAME2] = (!flag.deep && segconf.has_qname2); // note: we don't reset has_qname2 as fastq_zip_modify needs it to parse the read
        segconf.optimize[FASTQ_LINE3]  = (segconf.line3 != L3_EMPTY);

        // optimize QUAL unless already binned (8 is the number of bins in Illimina: https://sapac.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf)
        segconf.optimize[FASTQ_QUAL]   = (segconf_get_num_qual_scores(QHT_QUAL) > 8); 
    }
}

// called after segconf inc. segconf_finalize_optimize() which might remove optimizations
void fastq_zip_after_segconf (void)
{
    if (IS_R1) {
        double est_num_vbs = MAX_(1, (double)txtfile_get_seggable_size() / (double)segconf.vb_size * 1.1);

        // allocate memory to store txt_data.len32 of each R1 VB
        buf_alloc (evb, &z_file->R1_txt_data_lens, 0, est_num_vbs, uint32_t, 0, "z_file->R1_txt_data_lens");
        
        // allocate memory to store the last qname (canonized) of each R1 VB
        uint32_t canonical_len = strlen (segconf.qname_line0[QNAME1].s);
        qname_canonize (QNAME1, segconf.qname_line0[QNAME1].s, &canonical_len);
        buf_alloc (evb, &z_file->R1_last_qname, 0, est_num_vbs * (1 + canonical_len), char, 0, "z_file->R1_last_qname");
        buf_alloc (evb, &z_file->R1_last_qname_index, 0, est_num_vbs, uint64_t, 0, "z_file->R1_last_qname_index");        
    }

    if (segconf.optimize[FASTQ_LINE3]) {
         segconf.line3 = L3_EMPTY;

         if (segconf.desc_is_l3)
            segconf.seq_len_dict_id.num = 0;
    }

    if (segconf.optimize[FASTQ_QNAME2]) {
        segconf.qname_flavor[QNAME2] = NULL;
        ZCTX(FASTQ_QNAME2)->st_did_i = FASTQ_QNAME; // consolidate_stats doesn't work for QNAME2 because it is not merged if optimized 
    }

    if (segconf.optimize[FASTQ_QNAME]) {
        segconf.qname_flavor[QNAME1] = qname_get_optimize_qf();
        fastq_get_optimized_qname_read_name();
    }
}

void fastq_seg_finalize (VBlockP vb)
{
    // assign the QUAL codec
    codec_assign_best_qual_codec (vb, FASTQ_QUAL, fastq_zip_qual, segconf.deep_has_trimmed, false, NULL);
    
    if (segconf.has_agent_trimmer) 
        codec_assign_best_qual_codec (vb, OPTION_QX_Z, NULL, true, false, NULL);
    
    // top level snip
    SmallContainer top_level = { 
        .repeats        = vb->lines.len32,
        .is_toplevel    = true,
        .filter_items   = true,
        .filter_repeats = true,
        .callback       = true,
        .nitems_lo      = FASTQ_NUM_TOP_LEVEL_FIELDS, // initial before removals
        .items          = {
            { .dict_id = { _FASTQ_DEEP      }, .separator = { CI0_INVISIBLE } }, // 0:  note: will be removed if not --deep. added v15
            { .dict_id = { _FASTQ_QNAME     },                                }, // 1:
            { .dict_id = { _FASTQ_QNAME2    },                                }, // 2:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_EXTRA     },                                }, // 3:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_AUX       },                                }, // 4:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_E1L       },                                }, // 5:  note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
            { .dict_id = { _FASTQ_SQBITMAP  },                                }, // 6:  
            { .dict_id = { _FASTQ_E2L       },                                }, // 7:
            { .dict_id = { _FASTQ_LINE3     },                                }, // 8:  added in 12.0.14, before '+' was a separator of the previous E2L. v15: note: removed if not needed
            { .dict_id = { _FASTQ_QNAME2    },                                }, // 9 : note: removed if not needed. added v15
            { .dict_id = { _FASTQ_EXTRA     },                                }, // 10: note: removed if not needed. added v15
            { .dict_id = { _FASTQ_AUX       },                                }, // 11: note: removed if not needed. added v15
            { .dict_id = { _FASTQ_E2L       },                                }, // 12:
            { .dict_id = { _FASTQ_QUAL      },                                }, // 13:
            { .dict_id = { _FASTQ_E2L       },                                }, // 14:
            { segconf.seq_len_dict_id,         .separator = { CI0_INVISIBLE } }  // 15: reconstruct for seq_len in case QNAME and LINE3 are dropped with --seq-only or --qual-only. See logic in fastq_piz_initialize_item_filter
        }
    };

    // prefixes in this container were added in 12.0.14, before, '@' was part of DESC and '+' was a separator
    char prefixes[] = { CON_PX_SEP,        // 0:  initial
                        CON_PX_SEP,        // 1:  terminator for empty container-wide prefix
                        CON_PX_SEP,        // 2:  empty DEEP line prefix. note: will be removed below if not --deep
                        DC, CON_PX_SEP,    // 3:  QNAME prefix
                        ' ', CON_PX_SEP,   // 5:  QNAME2 prefix
                        ' ', CON_PX_SEP,   // 7:  EXTRA prefix
                        " \t"[segconf.saux_tab_sep], CON_PX_SEP,   // 9:  AUX or SAUX prefix
                        CON_PX_SEP,        // 11: empty E1L line prefix
                        CON_PX_SEP,        // 12: empty SQBITMAP line prefix
                        CON_PX_SEP,        // 13: empty E2L line prefix
                        '+', CON_PX_SEP,   // 14: LINE3 prefix
                        ' ', CON_PX_SEP,   // 16: QNAME2 prefix
                        ' ', CON_PX_SEP,   // 18: EXTRA prefix
                        ' ', CON_PX_SEP,   // 20: AUX prefix
                        CON_PX_SEP,        // 22: E2L prefix
                        CON_PX_SEP,        // 23: QUAL prefix
                        CON_PX_SEP,        // 24: E2L prefix
                        CON_PX_SEP,        // 25: seq_len_dict_id prefix
                        CON_PX_SEP      }; // 26: END of prefixes

    int prefixes_len = sizeof (prefixes);

    #define REMOVE(item_i,px_i,px_len) ({ \
        memmove (&top_level.items[item_i], &top_level.items[item_i+1], (top_level.nitems_lo - item_i - 1) * sizeof (ContainerItem)); \
        top_level.nitems_lo--; \
        memmove (&prefixes[px_i], &prefixes[px_i+px_len], sizeof (prefixes)-(px_i+px_len)); \
        prefixes_len -= px_len; })

    bool has_line3  = segconf.line3 != L3_EMPTY;
    
    // whether the Description (QNAME2 + EXTRA + AUX) appears on line 1 or line 3 
    // note: if on both, the line 3 is just a copy snip from line 1
    // note: if --optimize-DESC, has_line3==desc_is_l3==false
    bool desc_is_l3 = segconf.desc_is_l3; 
    
    bool has_qname2 = segconf.has_qname2;
    bool has_extra  = VB_FASTQ->has_extra;
    bool has_aux    = (segconf.has_aux || segconf.has_saux);

    // remove unneeded container and prefix items - in reverse order
    if (!segconf.seq_len_dict_id.num)      REMOVE (15, 25, 1);
    if (FAF)                               REMOVE (14, 24, 1);
    if (FAF)                               REMOVE (13, 23, 1);    
    if (FAF)                               REMOVE (12, 22, 1);    
    if (FAF || !desc_is_l3 || !has_aux)    REMOVE (11, 20, 2);
    if (FAF || !desc_is_l3 || !has_extra)  REMOVE (10, 18, 2);
    if (FAF || !desc_is_l3 || !has_qname2) REMOVE (9,  16, 2);
    if (FAF)                               REMOVE (8,  14, 2);    
    if (!has_line3)                        REMOVE (8,  15, 1); // removing CON_PX_SEP the '+' now becomes the prefix of E2L (as QNAME2, EXTRA, AUX have already been removed)
    if (desc_is_l3 || !has_aux)            REMOVE (4,  9,  2);
    if (desc_is_l3 || !has_extra)          REMOVE (3,  7,  2);
    if (desc_is_l3 || !has_qname2)         REMOVE (2,  5,  2);
    if (!flag.deep)                        REMOVE (0,  2,  1);

    container_seg (vb, CTX(FASTQ_TOPLEVEL), (ContainerP)&top_level, prefixes, prefixes_len, 0); // note: the '@', '+' and ' ' are accounted for in the QNAME, QNAME2/EXTRA/AUX and LINE3 fields respectively
}

bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _FASTQ_TOPLEVEL ||
           dict_id.num == _FASTQ_QNAME    ||
           dict_id.num == _FASTQ_E1L      ||
           dict_id.num == _FASTQ_E2L;
}

// ZIP/PIZ main thread: called ahead of zip or piz a pair 2 vb - to read data we need from the previous pair 1 file
// returns true if successful, false if there isn't a vb with vb_i in the previous file
void fastq_read_R1_data (VBlockP vb_, VBIType R1_vb_i)
{
    START_TIMER;
    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    if (flag.no_zriter) zriter_flush(); 

    Section sec = sections_vb_header (R1_vb_i);

    vb->R1_vb_i      = R1_vb_i;
    vb->R1_num_lines = sec->num_lines;

    // read into ctx->pair the data we need from our pair: QNAME,QNAME2,LINE3 and its components, GPOS and STRAND
    buf_alloc (vb, &vb->z_section_headers, MAX_DICTS * 2, 0, uint32_t, 0, "z_section_headers"); // indices into vb->z_data of section headers
        
    piz_read_all_ctxs (VB, &sec, true);
    
    if (flag.no_zriter) 
        file_seek (z_file, 0, SEEK_END, READ, HARD_FAIL); // restore

    COPY_TIMER (fastq_read_R1_data);
}

// main thread: after reading VB_HEADER and before reading local/b250 sections from z_file
void fastq_piz_before_read (VBlockP vb)
{
    if (writer_am_i_pair_2 (vb->vblock_i, &VB_FASTQ->R1_vb_i)) { // sets R1_vb_i if R2, leaves it 0 if R1
        
        // backward compatability: prior to V15, PIZ didn't rely on FlagCtx.paired, and it was not always
        // applied: SQBITMAP.b250 in v14 and GPOS.local (at least) in v14 incorrectly didn't set the flag
        if (!VER(15)) {
            CTX(FASTQ_GPOS)->pair_assist_type = CTX(FASTQ_STRAND)->pair_assist_type = SEC_LOCAL;
            if (VER(14)) CTX(FASTQ_SQBITMAP)->pair_assist_type = SEC_B250;
        }
    }
}

// main thread: called from piz_read_one_vb as DTP(piz_init_vb)
bool fastq_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header)
{
    // in case of this is a R2 of a paired fastq file, get the R1 data
    if (vb && VB_FASTQ->R1_vb_i > 0)
        fastq_read_R1_data (vb, VB_FASTQ->R1_vb_i);

    return true;
}

static rom fastq_seg_get_lines (VBlockFASTQP vb, rom line, int32_t remaining, 
                                pSTRp(qname), pSTRp(qname2), pSTRp(desc), 
                                pSTRp(seq), pSTRp(line3), pSTRp(qual), uint32_t *line1_len, bool *has_13) // out
{
    START_TIMER;

    ASSSEG0 (*line != '\n', "Invalid FASTQ file format: unexpected newline");

    ASSSEG (*line == DC, "Invalid FASTQ file format: expecting description line to start with '%c' but it starts with %c", DC, *line);     

    *qname = line+1; remaining--; // skip the '@' - its already included in the container prefix
    
    bool analyze = segconf.running && vb->line_i==0;

    // get QNAME + DESC line
    *seq = seg_get_next_line (VB, *qname, &remaining, line1_len, true, &has_13[0], "LINE1");
    
    // segconf: test for interleaved
    if (segconf.running && segconf.is_interleaved != no) {
        if (vb->line_i % 2)
            // note: a "no" is conclusive. a "yes" is evidence and remains so unless a "no" is discovered in a future line
            segconf.is_interleaved = fastq_zip_is_interleaved (STRtxt (CTX(FASTQ_QNAME)->last_line1), *qname, *line1_len); 
        else 
            CTX(FASTQ_QNAME)->last_line1 = (TxtWord){ .index = BNUMtxt(*qname), .len = *line1_len };
    }

    // analyze Line1 (segconf)
    if (analyze) {
        if (flag.deep) 
            memcpy (segconf.deep_1st_desc, *qname, MIN_(*line1_len, sizeof (segconf.deep_1st_desc)-1));

        *desc = memchr2 (*qname, ' ', '\t', *line1_len);

        // if \t was found before any space, this *might* be SAUX (or \t might belong to qname)
        if (*desc && !fastq_segconf_analyze_saux (vb, *desc + 1, *qname + *line1_len - *desc - 1) && **desc == '\t') 
            *desc = memchr (*qname, ' ', *line1_len); // \t separator but not SAUX - tab belongs to qname

        // Discover the QNAME flavor. Also discovers the original TECH, even when optimize_DESC
        qname_segconf_discover_flavor (VB, QNAME1, *qname, (*desc ? *desc - *qname : *line1_len)); 

        if (!segconf.has_saux && *desc && **desc == ' ') 
            fastq_segconf_analyze_DESC (vb, *desc + 1, *qname + *line1_len - *desc - 1);
    }

    // get desc and qname2
    if (segconf.has_desc && !segconf.desc_is_l3 && (*desc = memchr (*qname, " \t"[segconf.saux_tab_sep], *line1_len))) {
        
        (*desc)++; // skip separator 
        *desc_len = *qname + *line1_len - *desc; 
        *qname_len = *line1_len - *desc_len - 1;

        rom after;
        if (segconf.has_qname2 && // note: has_saux and has_qname2 are mutually exclusive
            (after = (segconf.has_aux || segconf.has_extra) ? memchr (*desc, ' ', *desc_len) : (*desc + *desc_len))) {
            *qname2 = *desc;
            *qname2_len = after - *qname2;
        }
        else 
            *qname2_len = 0;
    }

    // entire line will be segged in seg_qname. hopefully it only contains a qname, but if not - it will be tokenized
    else {
        *qname_len = *line1_len;
        *desc_len = *qname2_len = 0;
    }

    // get SEQ line
    rom after = seg_get_next_item (VB, *seq, &remaining, GN_SEP, GN_IGNORE, GN_IGNORE, seq_len, NULL, &has_13[1], "SEQ");

    // case FASTA as FASTQ: no Line3 and no QUAL
    if (FAF) goto done;

    *line3 = after;
    ASSSEG (remaining && (*line3)[0] == '+', "Invalid FASTQ file format (#3): expecting middle line to be a \"+\", but it starts with a '%c'", (*line3)[0]);
    (*line3)++; // skip '+';

    // get LINE3
    *qual = seg_get_next_line (VB, *line3, &remaining, line3_len, true, &has_13[2], "LINE3");

    // analyze Line3 (segconf). note: if flag.optimize, we will update in fastq_seg_finalize
    if (analyze) {
        if (*line3_len == 0) 
            segconf.line3 = L3_EMPTY;

        else if (fastq_is_line3_copy_of_line1 (STRa(*qname), STRa(*line3), *desc_len))
            segconf.line3 = L3_COPY_LINE1;

        // test for e.g. Line1: "@3/1" Line3: "+SRR2982101.3 3 length=101"
        else if (!segconf.has_desc) {
            *desc = memchr (*line3, ' ', *line3_len); // note: SAUX not supported on line 3 (never observed in the wild)

            qname_segconf_discover_flavor (VB, QLINE3, *line3, *desc - *line3);

            if (*desc) fastq_segconf_analyze_DESC (vb, *desc + 1, *line3 + *line3_len - *desc - 1);

            segconf.desc_is_l3 = true; // DESC appears on line 3 rather than line 1
            segconf.line3 = L3_NCBI;
        }

        else 
            ABORT ("Unsupported FASTQ file format (#4): expecting middle line to be a \"+\" with or without a copy of the description, or NCBI, but it is \"%.*s\"",
                   (*line3_len)+1, (*line3)-1);
    }

    if (segconf.has_desc && segconf.desc_is_l3) {
        *desc = memchr (*line3, ' ', *line3_len); // note: not \t, bc SAUX not supported on line 3 (never observed in the wild)
        if (*desc) {
            (*desc)++; // skip space
            *desc_len = *line3 + *line3_len - *desc;
            *line3_len -= *desc_len + 1; // just QLINE3 remaining
        }
        else
            *desc_len = 0;
    }
    else {
        // entire line will be segged in seg_qname. hopefully it only contains a qname, but if not - it will be tokenized
    }

    // get QUAL line
    after = seg_get_next_item (VB, *qual, &remaining, GN_SEP, GN_FORBIDEN, GN_IGNORE, qual_len, NULL, &has_13[3], "QUAL"); 

    ASSSEG (*qual_len == *seq_len, "Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ = %.*s\nQUAL= %.*s",
            *seq_len, *qual_len, STRf(*seq), STRf(*qual));

    ASSSEG (str_is_in_range (STRa(*qual), 33, 126), "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", STRf(*qual));

done:
    COPY_TIMER (fastq_seg_get_lines);

    return after;
}

// --optimize: re-write VB before digest and seg
rom fastq_zip_modify (VBlockP vb_, rom line_start, uint32_t remaining)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    // Split read to to desc, seq, line3 and qual (excluding the '@' and '+')
    STR0(qname); STR0(qname2); STR0(desc); STR0(seq); STR0(line3); STR0(qual);
    bool line_has_13[4] = {}; // initialize as expected by downstream functions
    uint32_t line1_len;
    rom after = fastq_seg_get_lines (vb, line_start, remaining, pSTRa(qname), pSTRa(qname2), pSTRa(desc), 
                                     pSTRa(seq), pSTRa(line3), pSTRa(qual), &line1_len, line_has_13);

    buf_alloc (vb, &vb->optimized_line, 0, after - line_start + segconf.optimized_qname_len + 20, char, 0, "optimized_line"); // +20 for worst case scenario that original DESC was of negligible length, and then we add a long number
    char *next = B1STc (vb->optimized_line);

    *next++ = DC;

    // case: --optimize_DESC: we replace the description with "filename.line_i" (optimized string stored in segconf.optimized_qname)
    if (segconf.optimize[FASTQ_QNAME]) {
        char *save_next = next;
        next = mempcpy (next, segconf.optimized_qname.s, segconf.optimized_qname_len);
        next += str_int (vb->first_line + vb->line_i, next);

        // preserve mate (/1 or /2) if mate is on QNAME (possibly is on QNAME2 in which case it is copied anyway)
        if (qf_is_mated(QNAME1) && qname_len > 2 && qname[qname_len-2] == '/' && (qname[qname_len-1] == '1' || qname[qname_len-1] == '2')) {
            *next++ = '/';
            *next++ = qname[qname_len-1];
        }

        CTX(FASTQ_QNAME)->txt_shrinkage += qname_len - (next - save_next);

        // case: we have QNAME2 - copy line except QNAME and QNAME2
        if (segconf.optimize[FASTQ_QNAME2] && qname2_len) {
            next = mempcpy (next, qname2 + qname2_len, line1_len - qname_len - 1/*sep*/ - qname2_len);
            
            CTX(FASTQ_QNAME2)->txt_shrinkage += qname2_len + 1/*sep*/;
        }

        // case: no QNAME2 - copy line except QNAME
        else
            next = mempcpy (next, qname + qname_len, line1_len - qname_len); // remainder of line1 is unchanged
    }

    else 
        next = mempcpy (next, qname, line1_len);

    *next++ = '\n'; // note: modify always excludes any \r

    // unmodified SEQ
    next = mempcpy (next, seq, seq_len);
    *next++ = '\n'; 

    if (!FAF) {
        *next++ = '+'; 

        // case: optimize -> empty line3
        if (segconf.optimize[FASTQ_LINE3]) 
            CTX(FASTQ_LINE3)->txt_shrinkage += line3_len + (segconf.desc_is_l3 ? (1 + desc_len) : 0);
        else {
            next = mempcpy (next, line3, line3_len);
            if (segconf.desc_is_l3 && desc_len) {
                *next++ = ' ';
                next = mempcpy (next, desc, desc_len);
            }
        }

        *next++ = '\n'; 

        // QUAL: optimize or not
        next = segconf.optimize[FASTQ_QUAL] ? optimize_phred_quality_string (STRa(qual), next, false, false) 
                                            : mempcpy (next, qual, qual_len);
        *next++ = '\n'; 
    }

    vb->optimized_line.len = BNUM (vb->optimized_line, next);
    
    // account for the shrinkage due to the \r's we got rid of
    CTX(FASTQ_E1L)->txt_shrinkage += line_has_13[0];
    CTX(FASTQ_E2L)->txt_shrinkage += line_has_13[1] + line_has_13[2] + line_has_13[3];

    return after;
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
rom fastq_seg_txt_line (VBlockP vb_, rom line_start, uint32_t remaining, bool *has_13)     // index in vb->txt_data where this line starts
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    // Split read to to desc, seq, line3 and qual (excluding the '@' and '+')
    STR0(qname); STR0(qname2); STR0(desc); STR0(seq); STR0(line3); STR0(qual);
    bool line_has_13[4] = {};
    uint32_t line1_len;
    rom after = fastq_seg_get_lines (vb, line_start, remaining, pSTRa(qname), pSTRa(qname2), pSTRa(desc), pSTRa(seq), pSTRa(line3), pSTRa(qual), 
                                     &line1_len, line_has_13);
    
    // set dl fields, consumed by fastq_zip_qual/seq
    ZipDataLineFASTQ *dl = DATA_LINE (vb->line_i);
    dl->seq  = TXTWORD(seq); 
    
    if (!FAF) 
        dl->qual = TXTWORD(qual);

    dl->monochar = str_is_monochar (STRa(seq)); // set dl->monochar

    // case --deep: compare DESC, SEQ and QUAL to the SAM/BAM data
    bool deep_qname=false, deep_seq=false;
    uint32_t uncanonical_suffix_len = 0;
    
    if (flag.deep || flag.show_deep == 2)
        fastq_seg_deep (vb, dl, STRa(qname), STRa(qname2), STRa(seq), STRa(qual), &deep_qname, &deep_seq, &dl->deep_qual, &uncanonical_suffix_len);

    fastq_seg_QNAME (vb, STRa(qname), line1_len, segconf.deep_qtype == QNAME1 && deep_qname, uncanonical_suffix_len);

    // seg SAUX
    if (segconf.has_saux)
        fastq_seg_saux (vb, STRa(desc));

    // seg DESC (i.e. QNAME2 + EXTRA + AUX), it might have come either from line1 or from line3
    else if (segconf.has_desc) 
        fastq_seg_DESC (vb, STRa(desc), segconf.deep_qtype == QNAME2 && deep_qname, uncanonical_suffix_len);

    if (!FAF)
        fastq_seg_LINE3 (vb, STRa(line3), STRa(qname), STRa(desc)); // before SEQ, in case it as a segconf_seq_len_dict_id field

    fastq_seg_SEQ (vb, dl, STRa(seq), deep_seq);

    if (!FAF)
        fastq_seg_QUAL (vb, dl, STRa(qual));

    // 4 end of lines. note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
    for (int i=0; i < (FAF ? 2 : 4); i++) {
        *has_13 = line_has_13[i]; 
        SEG_EOL (i==0 ? FASTQ_E1L : FASTQ_E2L, true);
    }

    // if seq_len_dict_id was detected in segconf line 0, it must appear in all lines
    ASSERT (!segconf.seq_len_dict_id.num || ctx_has_value_in_line (VB, segconf.seq_len_dict_id, NULL), 
            "Line missing length component (ctx=%s qname=\"%.*s\")", ECTX(segconf.seq_len_dict_id)->tag_name, STRf(qname));

    return after;
}

rom fastq_assseg_line (VBlockP vb)
{
    if (vb->line_start >= Ltxt) return "Invalid line_start";

    // set nl to the 4th newline, or end of txt_data if there are no 4 newliens
    char *nl = Btxt(vb->line_start) - 1;
    for (int i=0; i < 4; i++) {
        nl = memchr (nl+1, '\n', BAFTtxt - (nl+1));
        if (!nl) {
            nl = BAFTtxt; // possibly overwriting txt_data's overflow fence
            break;
        }
    }

    *nl = 0; // terminate string
    return Btxt(vb->line_start);
}

//-----------------
// PIZ stuff
//-----------------

// main thread: called for each txt file, after reading global area, before reading txt header
bool fastq_piz_initialize (CompIType comp_i)
{
    ASSINP0 (!(FAF && flag.qual_only), "--qual-only is only supported for FASTQ files"); // no --qual-only for FASTA (same error as in flags.c)

    return true;
}

void fastq_piz_header_init (void)
{
    if (flag.qnames_file)
        qname_filter_initialize_from_file (flag.qnames_file); // note: uses SectionHeaderTxtHeader.flav_props to canonize qnames

    if (flag.qnames_opt)
        qname_filter_initialize_from_opt (flag.qnames_opt); 
}

// PIZ: main thread: piz_process_recon callback: usually called in order of VBs, but out-of-order if --test with no writer
void fastq_piz_process_recon (VBlockP vb)
{
    if (flag.collect_coverage)    
        coverage_add_one_vb (vb);
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fastq_piz_is_skip_section)
{
    if (!ST(LOCAL) && !ST(B250) && !ST(DICT)) return false; // we only consider context data for skipping

    #define LINE3_dicts _FASTQ_LINE3,  _FASTQ_T0HIRD, _FASTQ_T1HIRD, _FASTQ_T2HIRD, _FASTQ_T3HIRD, _FASTQ_T4HIRD, _FASTQ_T5HIRD, \
                        _FASTQ_T6HIRD, _FASTQ_T7HIRD, _FASTQ_T8HIRD, _FASTQ_T9HIRD, _FASTQ_TAHIRD, _FASTQ_TBHIRD, _FASTQ_TmHIRD /* just line3, not all qnames */

    // we always keep the seq_len context, if we have one
    if (dict_id.num && dict_id.num == segconf.seq_len_dict_id.num) return false;

    if (codec_pacb_smux_is_qual (dict_id)) dict_id.num = _FASTQ_QUAL;

    // note that flags_update_piz_one_z_file rewrites --header-only as flag.header_only_fast: skip all items but DESC and E1L (except if we need them for --grep)
    if (flag.header_only_fast && 
        dict_id_is_in (dict_id, _FASTQ_DEBUG_LINES, LINE3_dicts,
                       _FASTQ_SQBITMAP, _FASTQ_NONREF, _FASTQ_NONREF_X,
                        _FASTQ_GPOS, _FASTQ_GPOS_DELTA, _FASTQ_GPOS_R2, _FASTQ_STRAND, _FASTQ_STRAND_R2,
                       _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, 
                       _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE))
        return true;

    if (flag.seq_only && 
        (   dict_id_is_in (dict_id, _FASTQ_E1L, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_EXTRA, _FASTQ_DEBUG_LINES, _FASTQ_LINE3,
                           _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE)  
         || dict_id_is_fastq_qname_sf(dict_id) || dict_id_is_fastq_aux(dict_id))) // we don't need the DESC line 
        return true;

    // note: we need SQBITMAP to extract seq_len. We even needs its local, bc otherwise SNIP_LOOKUP won't work.
    if (flag.qual_only && !FAF && 
        (   dict_id_is_in (dict_id, _FASTQ_E1L, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_EXTRA, _FASTQ_DEBUG_LINES, DICT_ID_NONE) 
         || dict_id_is_fastq_qname_sf(dict_id) || dict_id_is_fastq_aux(dict_id)
         || dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_GPOS_DELTA, _FASTQ_GPOS_R2, _FASTQ_STRAND, _FASTQ_STRAND_R2, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE))) // we don't need the SEQ line 
        return true;

    // if we're doing --sex/coverage, we only need TOPLEVEL, FASTQ_SQBITMAP and GPOS
    if (flag.collect_coverage && 
        (   dict_id_is_in (dict_id, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_EXTRA, _FASTQ_DEBUG_LINES,
                           _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE)
         || dict_id_is_fastq_qname_sf(dict_id) || dict_id_is_fastq_aux(dict_id)
         || (!flag.bases && dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE))))
        return true;

    // note: we don't SKIP for --count with an additional filter. Logic is too complicated and bug-prone.
    if (flag.count && !flag.grep && IS_DICTED_SEC (st) && !flag.bases && !flag.regions && 
        dict_id.num != _FASTQ_TOPLEVEL) 
        return true;

    if (dict_id.num == _FASTQ_TAXID) return true; // skip TAXID in older files (12.0.2 - 15.0.41)

    return false;
}

// called before reconstructing first repeat of toplevel container
static inline void fastq_piz_initialize_item_filter (VBlockFASTQP vb, ConstContainerP con)
{
    ASSERT (sizeof (vb->item_filter) >= con->nitems_lo, "top_level.nitems_lo=%u by item_filter is only %u", con->nitems_lo, (int)sizeof (vb->item_filter));

    memset (vb->item_filter, true, sizeof (vb->item_filter)); // default - reconstruct all items

    // find newline items
    int nl[4];
    int n_nl = 0;
    int expected_n_nl = (FAF ? 2 : 4);
    for (int i=0; i < con->nitems_lo && n_nl < expected_n_nl/*safety*/; i++)
        if (con->items[i].did_i_small == FASTQ_E1L || con->items[i].did_i_small == FASTQ_E2L)
            nl[n_nl++] = i;

    ASSERT (n_nl == expected_n_nl, "expecting %d newlines in a FASTQ container but found %d", expected_n_nl, n_nl);

    #define KEEP_ITEMS(first_item, last_item) \
        ({ memset (&vb->item_filter[flag.deep], false, (first_item)-flag.deep); /*always keep DEEP (item 0, if flag.deep) */\
           memset (&vb->item_filter[(last_item)+1], false, con->nitems_lo - ((last_item)+1)); })

    // --header-only: keep items up to and including first newline
    if (flag.header_only_fast) KEEP_ITEMS (flag.deep, nl[0]);

    // --seg-only: keep items after first newline, up to second newline
    if (flag.seq_only) KEEP_ITEMS(nl[0]+1, nl[1]);

    // --qual-only: keep items after third newline, up to forth newline
    if (flag.qual_only) {
        ASSINP0 (CTX(FASTQ_QUAL)->lcodec != CODEC_LONGR, 
                 "Genozip limitation: --qual-only cannot be used on this file, because base quality scores were compressed with the LONGR codec that requires sequence data to reconstruct correctly.");

        KEEP_ITEMS(nl[2]+1, nl[3]);
        vb->item_filter[nl[0]+1] = true; // keep SQBITMAP to extract seq_len, but fastq_recon_aligned_SEQ / fastq_special_unaligned_SEQ skip reconstructing if qual_only
    }

    // case: we have a seq_len item. we keep it if and only if we are dropping the first and third lines, that normally contain it
    if (nl[n_nl-1] != con->nitems_lo - 1) // can only happen since v15
        vb->item_filter[con->nitems_lo-1] = (flag.seq_only || flag.qual_only); 
}

// filtering during reconstruction: called by container_reconstruct for each fastq record (repeat) and each toplevel item
CONTAINER_FILTER_FUNC (fastq_piz_filter)
{
    if (dict_id.num == _FASTQ_TOPLEVEL) {
        
        // initialize item filter
        if (item == -1 && rep == 0) {
            fastq_piz_initialize_item_filter (VB_FASTQ, con);

            if (flag.deep) // Deep, since v15
                fastq_deep_piz_wait_for_deep_data();
        }

        // keep or drop toplevel item based on item filter
        else if (item >= 0) 
            return VB_FASTQ->item_filter[item];
    }

    return true; // reconstruct
}

// filtering during reconstruction: called by container_reconstruct for each fastq record (repeat)
CONTAINER_CALLBACK (fastq_piz_container_cb)
{
    if (is_top_level) {
        #define DROP_LINE(reason) ({ vb->drop_curr_line = (reason); goto dropped; })
        
        // --bases
        if (flag.bases && !vb->drop_curr_line &&
            !iupac_is_included_ascii (STRlst(FASTQ_SQBITMAP)))
            DROP_LINE ("bases");

        // --qnames and --qnames-file
        if (flag.qname_filter && !qname_filter_does_line_survive (STRlst (FASTQ_QNAME)))
            DROP_LINE ("qname_filter");
        
        dropped: {}
    }
}

// Used in R2: used for pair-assisted b250 reconstruction. copy parallel b250 snip from R1. 
// Since v14, used for SQBITMAP. For files up to v14, also used for all QNAME subfield contexts
SPECIAL_RECONSTRUCTOR (fastq_special_mate_lookup)
{
    ASSPIZ (ctx->b250R1.len32, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2.", ctx->tag_name);
            
    ctx_get_next_snip (vb, ctx, true, pSTRa(snip));

    if (ctx->did_i == FASTQ_SQBITMAP && VER(14))
        ctx->r1_is_aligned = (snip_len && *snip == SNIP_LOOKUP) ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;

    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE /* we can't cache pair items */, STRa(snip), reconstruct, __FUNCLINE); // might include delta etc - works because in --pair, ALL the snips in a context are FASTQ_SPECIAL_mate_lookup

    return NO_NEW_VALUE; // last_value already set (if needed) in reconstruct_one_snip
}

// used for interleaved files - demux by R (since 15.0.58)
SPECIAL_RECONSTRUCTOR (fastq_special_DEMUX_BY_R)
{
    int channel_i = vb->line_i % 2;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}

void fastq_reset_line (VBlockP vb_)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;

    if (!(segconf.is_interleaved && vb->line_i % 2)) // note: in interleaved, we an R2 lines relies on the previous, R1, line    
        CTX(FASTQ_SQBITMAP)->r1_is_aligned = VER(14) ? PAIR1_ALIGNED_UNKNOWN 
                                                     : PAIR1_ALIGNED; // up to v13, all lines had alignment data, even if unmapped

    vb->sam_seq_offset = 0;
}

//---------------------------------------------------------
// FASTQ-specific fields in genozip header
//---------------------------------------------------------

void fastq_zip_genozip_header (SectionHeaderGenozipHeader *header)
{
    header->fastq.segconf_seq_len_dict_id = segconf.seq_len_dict_id; // v14
    header->fastq.segconf_fa_as_fq        = segconf.fasta_as_fastq;  // 15.0.58
    header->fastq.segconf_is_ileaved      = segconf.is_interleaved;  // 15.0.58
}

void fastq_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header)
{
    if (VER(14)) 
        segconf.seq_len_dict_id = header->fastq.segconf_seq_len_dict_id; 
    
    if (VER(15)) {
        segconf.fasta_as_fastq  = header->fastq.segconf_fa_as_fq;
        segconf.is_interleaved  = header->fastq.segconf_is_ileaved;
    } 
}

