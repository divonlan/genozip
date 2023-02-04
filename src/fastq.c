// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <dirent.h>
#include <libgen.h>
#include "fastq_private.h"
#include "seg.h"
#include "piz.h"
#include "optimize.h"
#include "codec.h"
#include "writer.h"
#include "kraken.h"
#include "bases_filter.h"
#include "qname.h"
#include "tip.h"
#include "deep.h"
#include "tokenizer.h"
#include "coverage.h"

#define dict_id_is_fastq_desc_sf dict_id_is_type_1
#define dict_id_fastq_desc_sf dict_id_type_1

// Used by Seg
static char copy_desc_snip[30];
static unsigned copy_desc_snip_len;

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

void fastq_vb_release_vb (VBlockFASTQP vb)
{
    vb->pair_num_lines = vb->pair_vb_i = vb->optimized_desc_len = 0;
    vb->first_line = 0;
    FREE (vb->optimized_desc);
}

void fastq_vb_destroy_vb (VBlockFASTQP vb)
{
}

//-----------------------
// TXTFILE stuff
//-----------------------

// detect if a generic file is actually a FASTQ - we call based on first character being '@', line 3 starting with '+', and line 1 and 3 being the same length
bool is_fastq (STRp(header), bool *need_more)
{
    if (!header_len || header[0] != '@' || !str_is_printable (STRa(header))) return false;

    int num_newlines = str_count_char (STRa(header), '\n');
    if (num_newlines > 10000) return false; // if num_newlines is huge (for a HEADER_BLOCK=256K header), this is most likely not FASTQ, and we end here as this might cause a stack overflow

    if (num_newlines < 4) {
        *need_more = true; // we can't tell yet - need more data
        return false;
    }

    str_split (header, header_len, num_newlines+1, '\n', line, false);

    return line_lens[1] > 0 && line_lens[1] == line_lens[3] && // SEQ and QUAL lines are of equal length
           line_lens[3] > 0 && lines[3][0] == '+';
}

// returns true if txt_data[txt_i] (which is a \n) is the end of a FASTQ record (= block of 4 lines in the file); -1 if out of data
static inline int fastq_is_end_of_line (VBlockP vb, uint32_t first_i, int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we are not at the EOF and the next char is not '@', then this is for sure not the end of the FASTQ record
    if (txt_i < vb->txt_data.len-1 && txt[txt_i+1] != '@') return false; 

    // if at the end of the data or at the end of a line where the next char is '@'. Verify that the previous line
    // is the '+' line, to prevent the '@' we're seeing actually the first quality score in QUAL, or the file not ending after the quality line
 
    // move two \n's back - the char after is expected to be a '+'
    unsigned count_nl = 0;
    for (int32_t i=txt_i-1; i >= (int32_t)first_i; i--) {
        if (txt[i] == '\n') count_nl++;
        if (count_nl == 2) return txt[i+1] == '+';
    }

    return -1; 
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i /* in/out */)
{    
    ASSERT (*i >= 0 && *i < vb->txt_data.len, "*i=%d is out of range [0,%"PRIu64"]", *i, vb->txt_data.len);

    for (; *i >= (int32_t)first_i; (*i)--) {
        // in FASTQ - an "end of line" is one that the next character is @, or it is the end of the file
        if (vb->txt_data.data[*i] == '\n')
            switch (fastq_is_end_of_line (vb, first_i, *i)) {
                case true  : return vb->txt_data.len-1 - *i; // end of line
                case false : continue;                       // not end of line, continue searching
                default    : goto out_of_data;  
            }
    }

out_of_data:
    return -1; // cannot find end-of-line in the data starting first_i
}
// called by txtfile_read_vblock when reading the 2nd file in a fastq pair - counts the number of fastq "lines" (each being 4 textual lines),
// comparing to the number of lines in the first file of the pair
// returns true if we have at least as much as needed, and sets unconsumed_len to the amount of excess characters read
// returns false is we don't yet have pair_1_num_lines lines - we need to read more
bool fastq_txtfile_have_enough_lines (VBlockP vb_, uint32_t *unconsumed_len, 
                                      uint32_t *my_lines, uint32_t *her_lines) // out - only set in case of failure
{
    VBlockFASTQP vb = (VBlockFASTQP )vb_;
    ASSERTNOTZERO (vb->pair_num_lines, "");

    rom next  = B1STtxt;
    rom after = BAFTtxt;

    uint32_t line_i; for (line_i=0; line_i < vb->pair_num_lines * 4; line_i++) {
        while (*next != '\n' && next < after) next++; 
        if (next >= after) {
            *my_lines  = line_i;
            *her_lines = vb->pair_num_lines * 4;
            return false;
        }
        next++; // skip newline
    }

    vb->lines.len = line_i/4;
    *unconsumed_len = after - next;
    return true;
}

//---------------
// ZIP / SEG stuff
//---------------

void fastq_zip_init_vb (VBlockP vb)
{
    // in case we're optimizing DESC in FASTQ, we need to know the number of lines
    if (flag.optimize_DESC) {
        uint32_t num_lines = str_count_char (STRb(vb->txt_data), '\n');
        ASSERT (num_lines % 4 == 0, "expecting number of txt lines in VB to be a multiple of 4, but found %u", num_lines);

        VB_FASTQ->first_line = txt_file->num_lines + 1;
        txt_file->num_lines += num_lines / 4;     // update here instead of in zip_update_txt_counters;
    }
}

// case of --optimize-DESC: generate the prefix of the read name from the txt file name
// eg. "../../fqs/sample.1.fq.gz" -> "@sample-1."
static void fastq_get_optimized_desc_read_name (VBlockFASTQP vb)
{
    vb->optimized_desc = MALLOC (strlen (txt_file->basename) + 30); // leave room for the line number to follow
    strcpy (vb->optimized_desc, txt_file->basename);
    file_get_raw_name_and_type (vb->optimized_desc, NULL, NULL); // remove file type extension
    vb->optimized_desc_len = strlen (vb->optimized_desc) + 1; // +1 for ':'

    // replace '.' in the filename with '-' as '.' is a separator in compound_seg and would needless inflate the number of contexts
    for (unsigned i=0; i < vb->optimized_desc_len-1; i++)
        if (vb->optimized_desc[i] == '.') vb->optimized_desc[i] = '-';

    vb->optimized_desc[vb->optimized_desc_len-1] = '.';
}

// test if an R2 file has been compressed after R1 without --pair - and produce tip if so
static void fastq_tip_if_should_be_pair (void)
{
    if (flag.pair || !txt_file->name || !flag.reference || txt_file->redirected) return;

    char dir_name_buf[strlen (txt_file->name) + 1];
    strcpy (dir_name_buf, txt_file->name);    
    rom dir_name = dirname (dir_name_buf); // destructive - might replace the last slash with 0

    bool is_cd = dir_name[0] == '.' && dir_name[1] == 0;

    DIR *dir = opendir (dir_name);
    if (!dir) return; // can't open directory

    rom bn = txt_file->basename; // basename
    int bn_len = strlen (bn);

    // find R2 file assuming we are R1, if there is one
    struct dirent *ent;
    while ((ent = readdir(dir))) {
        if (strlen (ent->d_name) != bn_len) continue; 

        int mismatches = 0, mm_i=0;
        rom dn = ent->d_name;
        for (int i=0; i < bn_len && mismatches < 2; i++)
            if (dn[i] != bn[i]) {
                mismatches++;
                mm_i = i;
            }
        
        // case: we found another genozip file with the same name as txt_file, except for '1'⇄'2' switch
        if (mismatches == 1 && dn[mm_i] == '2' && bn[mm_i] == '1') {

            char other_bn[bn_len + 1];
            strcpy (other_bn, bn);
            other_bn[mm_i] = '2';

            if (!flag.deep)
                TIP ("Using --pair to compress paired-end FASTQs can reduce the compressed file's size by 10%%. E.g.:\n"
                    "genozip --reference %s --pair %s %s%s%s\n",
                    ref_get_filename (gref), txt_name, (is_cd ? "" : dir_name), (is_cd ? "" : "/"), other_bn);

            segconf.r1_or_r2 = PAIR_READ_1;
            break;
        } 

        else if (mismatches == 1 && dn[mm_i] == '1' && bn[mm_i] == '2') {
            segconf.r1_or_r2 = PAIR_READ_2;
            break;
        }
    }

    closedir(dir);    
}

// called by main thread at the beginning of zipping this file
void fastq_zip_initialize (void)
{
    if (!flag.reference && !txt_file->redirected && !flag.multiseq)
        TIP ("compressing a FASTQ file using a reference file can reduce the compressed file's size by 20%%-60%%.\nUse: \"genozip --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n", txt_file->name);

    fastq_tip_if_should_be_pair();
    
    // reset lcodec for STRAND and GPOS, as these may change between PAIR_1 and PAIR_2 files
    ZCTX(FASTQ_STRAND)->lcodec = CODEC_UNKNOWN;
    ZCTX(FASTQ_GPOS  )->lcodec = CODEC_UNKNOWN;

    seg_prepare_snip_other (SNIP_COPY, _FASTQ_DESC, 0, 0, copy_desc_snip);

    // with REF_EXTERNAL, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (IS_REF_EXTERNAL && z_file->num_txts_so_far == 1) // first file
        ctx_populate_zf_ctx_from_contigs (gref, FASTQ_CONTIG, ref_get_ctgs (gref)); 

    qname_zip_initialize (FASTQ_DESC);
}

// called by main thread after each txt file compressing is done
void fastq_zip_finalize (bool is_last_user_txt_file)
{
    if (is_last_user_txt_file && flag.deep && flag.debug_deep)
        fastq_deep_show_entries_stats();
}

// called by zfile_compress_genozip_header to set FlagsGenozipHeader.dt_specific
bool fastq_zip_dts_flag (int dts)
{
    return (dts==1) ? (flag.pair != NOT_PAIRED_END) : 0;
}

// called by Compute thread at the beginning of this VB
void fastq_seg_initialize (VBlockFASTQP vb)
{
    START_TIMER;

    CTX(FASTQ_CONTIG)->flags.store = STORE_INDEX; // since v12

    Context *gpos_ctx     = CTX(FASTQ_GPOS);
    Context *strand_ctx   = CTX(FASTQ_STRAND);
    Context *sqbitmap_ctx = CTX(FASTQ_SQBITMAP);
    Context *nonref_ctx   = CTX(FASTQ_NONREF);
    Context *seqmis_ctx   = CTX(FASTQ_SEQMIS_A);

    if (flag.aligner_available) {
        strand_ctx->ltype       = LT_BITMAP;
        gpos_ctx->ltype         = LT_UINT32;
        gpos_ctx->flags.store   = STORE_INT;
        nonref_ctx->no_callback = true; // when using aligner, nonref data is in local. without aligner, the entire SEQ is segged into nonref using a callback
        
        sqbitmap_ctx->ltype     = LT_BITMAP; // implies no_stons
        sqbitmap_ctx->local_always = true;

        buf_alloc (vb, &sqbitmap_ctx->local, 1, vb->txt_data.len / 4, uint8_t, 0, "contexts->local"); 
        buf_alloc (vb, &strand_ctx->local, 0, roundup_bits2bytes64 (vb->lines.len), uint8_t, 0, "contexts->local"); 

        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, 1, vb->txt_data.len / 128, char, 0, "contexts->local"); 

        buf_alloc (vb, &gpos_ctx->local, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->local"); 
    }

    if (!flag.multiseq)
        codec_acgt_comp_init (VB, FASTQ_NONREF);

    // initialize QUAL to LT_SEQUENCE, it might be changed later to LT_CODEC (eg domq, longr)
    ctx_set_ltype (VB, LT_SEQUENCE, FASTQ_QUAL, DID_EOL);

    if (flag.pair == PAIR_READ_2) {

        ASSERT (vb->lines.len32 == vb->pair_num_lines, "in vb=%s (PAIR_READ_2): pair_num_lines=%u but lines.len=%u",
                VB_NAME, vb->pair_num_lines, vb->lines.len32);

        gpos_ctx->pair_local = strand_ctx->pair_local = true;

        piz_uncompress_all_ctxs (VB);

        vb->z_data.len32 = 0; // we've finished reading the pair file z_data, next, we're going to write to z_data our compressed output
    }

    qname_seg_initialize (VB, FASTQ_DESC);

    if (vb->comp_i == FQ_COMP_R2)
        ctx_create_node (VB, FASTQ_SQBITMAP, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_mate_lookup }, 2); // required by ctx_convert_generated_b250_to_mate_lookup

    if (flag.optimize_DESC) 
        fastq_get_optimized_desc_read_name (vb);

    if (kraken_is_loaded) {
        CTX(FASTQ_TAXID)->flags.store    = STORE_INT;
        CTX(FASTQ_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTQ_TAXID)->counts_section = true; 
    }

    if (flag.deep)
        fastq_deep_seg_initialize (VB);

    // consolidate stats for --stats
    ctx_consolidate_stats (VB, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND, FASTQ_SEQMIS_A, FASTQ_SEQMIS_C, FASTQ_SEQMIS_G, FASTQ_SEQMIS_T, DID_EOL);
    ctx_consolidate_stats (VB, FASTQ_QUAL, FASTQ_DOMQRUNS, FASTQ_QUALMPLX, FASTQ_DIVRQUAL, DID_EOL);

    COPY_TIMER (seg_initialize);
}

void fastq_seg_finalize (VBlockP vb)
{
    if (segconf.running) {
        segconf.is_long_reads = segconf_is_long_reads();

        if (flag.deep) fastq_deep_seg_finalize_segconf (vb->lines.len32);
    }

    // assign the QUAL codec
    codec_assign_best_qual_codec (vb, FASTQ_QUAL, fastq_zip_qual, false, false);

    // top level snip
    SmallContainer top_level = { 
        .repeats        = vb->lines.len32,
        .is_toplevel    = true,
        .filter_items   = true,
        .filter_repeats = true,
        .callback       = true,
        .nitems_lo      = 8 + flag.deep,
        .items          = {
            { .dict_id = { _FASTQ_DEEP     }, }, // note: will be removed below if not --deep
            { .dict_id = { _FASTQ_DESC     }, },
            { .dict_id = { _FASTQ_E1L      }, }, // note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
            { .dict_id = { _FASTQ_SQBITMAP }, },
            { .dict_id = { _FASTQ_E2L      }, },
            { .dict_id = { _FASTQ_LINE3    }, }, // added in 12.0.14, before '+' was a separator of the previous E2L
            { .dict_id = { _FASTQ_E2L      }, },
            { .dict_id = { _FASTQ_QUAL     }, },
            { .dict_id = { _FASTQ_E2L      }, } 
        }
    };

    // prefixes in this container were added in 12.0.14, before, '@' was part of DESC and '+' was a separator
    char prefixes[] = { CON_PX_SEP,        // initial
                        CON_PX_SEP,        // terminator for empty container-wide prefix
                        CON_PX_SEP,        // empty DEEP line prefix. note: will be removed below if not --deep
                        '@', CON_PX_SEP,   // DESC prefix
                        CON_PX_SEP,        // empty E1L line prefix
                        CON_PX_SEP,        // empty SQBITMAP line prefix
                        CON_PX_SEP,        // empty E2L line prefix
                        '+', CON_PX_SEP,   // third line prefix
                        CON_PX_SEP      }; // END of prefixes

    // if not --deep, remove from container and prefixes to avoid unnecessary overhead in PIZ
    if (!flag.deep) {
        memmove (&top_level.items[0], &top_level.items[1], top_level.nitems_lo * sizeof (ContainerItem));
        memmove (&prefixes[0], &prefixes[1], sizeof (prefixes)-1);
    }

    container_seg (vb, CTX(FASTQ_TOPLEVEL), (ContainerP)&top_level, prefixes, sizeof (prefixes) - !flag.deep, 0); // note: the '@' and '+' are accounted for in the DESC and LINE3 fields respectively
}

bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _FASTQ_TOPLEVEL ||
           dict_id.num == _FASTQ_DESC     ||
           dict_id.num == _FASTQ_TAXID    ||
           dict_id.num == _FASTQ_E1L      ||
           dict_id.num == _FASTQ_E2L;
}

// ZIP/PIZ main thread: called ahead of zip or piz a pair 2 vb - to read data we need from the previous pair 1 file
// returns true if successful, false if there isn't a vb with vb_i in the previous file
bool fastq_read_pair_1_data (VBlockP vb_, uint32_t pair_vb_i, bool must_have)
{
    VBlockFASTQP vb = (VBlockFASTQP )vb_;
    uint64_t save_disk_so_far = z_file->disk_so_far;

    vb->pair_vb_i = pair_vb_i;

    Section sec = sections_vb_header (pair_vb_i, SOFT_FAIL);
    ASSERT (!must_have || sec, "file unexpectedly does not contain data for pair 1 vb_i=%u", pair_vb_i);
    if (!sec) return false;

    vb->pair_num_lines = sec->num_lines;

    // read into ctx->pair the data we need from our pair: DESC and its components, GPOS and STRAND
    buf_alloc (vb, &vb->z_section_headers, MAX_FIELDS + 10, MAX_DICTS * 2 + 50, uint32_t, 0, "z_section_headers"); // room for section headers  

    sec++;
    piz_read_all_ctxs (VB, &sec, true);

    file_seek (z_file, 0, SEEK_END, false); // restore
    z_file->disk_so_far = save_disk_so_far;

    return true;
}

// main thread: is it possible that genocat of this file will re-order lines
bool fastq_piz_maybe_reorder_lines (void)
{
    return sections_get_num_comps() == 2; // genocat always interleaves (= reorders lines) of paired FASTQs
}

// main thread: called from piz_read_one_vb as DTPZ(piz_read_one_vb)
bool fastq_piz_init_vb (VBlockP vb_, const SectionHeaderVbHeader *header, uint32_t *txt_data_so_far_single_0_increment)
{
    VBlockFASTQP vb = (VBlockFASTQP )vb_;
    uint32_t pair_vb_i=0;
    bool i_am_pair_2 = vb && writer_get_pair (vb->vblock_i, &pair_vb_i) == 2;

    // in case of this is a paired fastq file, get all the pair_1 data not already fetched for the grep above
    if (i_am_pair_2) 
        fastq_read_pair_1_data (vb_, pair_vb_i, true);

    return true;
}

static void fastq_seg_DESC (VBlockFASTQP vb, STRp(desc), unsigned qname_len, bool deep)
{
    // case: --optimize_DESC: we replace the description with "filename.line_i" (optimized string stored in vb->optimized_desc)
    unsigned optimized_len = 0; 
    if (flag.optimize_DESC) {
        optimized_len  = vb->optimized_desc_len + str_int (vb->first_line + vb->line_i, &vb->optimized_desc[vb->optimized_desc_len]);   
        vb->recon_size -= desc_len - optimized_len;

        desc = vb->optimized_desc;
        desc_len = qname_len = optimized_len;
    }

    // case: segconf && deep && DESC is made out of qname followed a whitespace and more DESC - prepare container
    if (segconf.running && flag.deep && !segconf.deep_desc_con_snip_len && desc_len != qname_len) {
        SmallContainer con = { .repeats   = 1,
                               .nitems_lo = 1,
                               .items     = { { .dict_id = { _FASTQ_Q0NAME }, .separator = { desc[qname_len] } },
                                              { .dict_id = { _FASTQ_QNAME2 }                                   } } };

        segconf.deep_desc_con_snip_len = sizeof (segconf.deep_desc_con_snip);
        container_prepare_snip ((ContainerP)&con, 0, 0, qSTRa(segconf.deep_desc_con_snip));
        segconf.deep_desc_con_snip_sep = desc[qname_len];
    }

    // case: deep: desc is the same as SAM's qname
    if (deep && desc_len == qname_len)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_copy_deep }, 2, FASTQ_DESC, qname_len + 1); // +1 for '@'

    // case: deep + qname is the same as SAM's qname, but is followed by a whitespace and more DESC - and matches deep_desc_con_snip
    else if (deep && segconf.deep_desc_con_snip_len && desc[qname_len] == segconf.deep_desc_con_snip_sep) {
        seg_by_did (VB, STRa(segconf.deep_desc_con_snip), FASTQ_DESC, 2);                                  // container: 2 = container item sep and '@'        
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_copy_deep }, 2, FASTQ_Q0NAME, qname_len);    // QNAME part
        tokenizer_seg (VB, CTX(FASTQ_QNAME2), &desc[qname_len+1], desc_len - qname_len - 1, sep_with_space, 0); // non-QNAME part
    }

    // case: normal seg: attempt to parse the desciption line as qname, and fallback to tokenizer desc mismatches the flavor
    else
        qname_seg (VB, CTX(FASTQ_DESC), STRa(desc), 1); // account for the '@' (segged as a toplevel container prefix)
}

static void fastq_seg_LINE3 (VBlockFASTQP vb, STRp(line3), STRp(desc))
{
    // first line of segconf VB - discover line3 type
    if (segconf.running && vb->line_i==0) {

        if (line3_len==0 || flag.optimize_DESC) 
            segconf.line3 = L3_EMPTY;

        else if (str_issame_(STRa(line3), STRa(desc)))
            segconf.line3 = L3_COPY_DESC;

        else if (qname_segconf_discover_fastq_line3_sra_flavor (VB, STRa(line3)))
            segconf.line3 = L3_QF;

        else 
            ASSSEG (false, line3, "Invalid FASTQ file format: expecting middle line to be a \"+\" with or without a copy of the description, but it is \"%.*s\"",
                    line3_len+1, line3-1);
    }

    if (flag.optimize_DESC) {
        seg_by_did (VB, "", 0, FASTQ_LINE3, 1); // account for the '+' (segged as a toplevel container prefix)
        vb->recon_size -= line3_len; 
    }
    
    else if (segconf.line3 == L3_EMPTY) {
        ASSSEG (!line3_len || flag.optimize_DESC, line3, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it is \"%.*s\"", line3_len+1, line3-1);
        seg_by_did (VB, "", 0, FASTQ_LINE3, 1); // account for the '+' (segged as a toplevel container prefix)
    }
    
    else if (segconf.line3 == L3_COPY_DESC) {
        ASSSEG (str_issame_(STRa(line3), STRa(desc)), line3, 
                "Invalid FASTQ file format: expecting middle line to be a \"+\" followed by a copy of the description line, but it is \"%.*s\"", line3_len+1, line3-1); 
        seg_by_did (VB, STRa(copy_desc_snip), FASTQ_LINE3, 1 + line3_len); // +1 - account for the '+' (segged as a toplevel container prefix)
    }

    else if (segconf.line3 == L3_QF) 
        ASSSEG (qname_seg_qf (VB, CTX(FASTQ_LINE3), segconf.line3_flavor, STRa(line3), false, 1), line3, // +1 - account for the '+' (segged as a toplevel container prefix)
                "Invalid FASTQ file format: expecting middle line to be a \"+\" followed by a \"%s\" flavor, but it is \"%.*s\"", qf_name(segconf.line3_flavor), line3_len+1, line3-1); 

    else 
        ASSSEG (false, line3, "Invalid FASTQ file format: expecting middle line to be a \"+\" with or without a copy of the description, but it is \"%.*s\"",
                line3_len+1, line3-1);
}

static void fastq_seg_kraken_tax_id (VBlockFASTQP vb, STRp(qname))
{
    unsigned taxid_found = kraken_seg_taxid (VB, FASTQ_TAXID, STRa(qname), false);

    // if not found tax id for this read, try again, perhaps removing /1 or /2
    if (!taxid_found) {
        if (qname_len > 2 && qname[qname_len-2] == '/' &&
            (qname[qname_len-1] == '1' || qname[qname_len-1] == '2'))
            qname_len -= 2;

        kraken_seg_taxid (VB, FASTQ_TAXID, STRa(qname), true); // this fails if missing
    }

}

static rom fastq_seg_get_lines (VBlockFASTQP vb, rom line, int32_t remaining, 
                                pSTRp(desc), pSTRp(seq), pSTRp(line3), pSTRp(qual), bool *has_13) // out
{
    ASSSEG0 (*line != '\n', line, "Invalid FASTQ file format: unexpected newline");

    ASSSEG (*line == '@', line, "Invalid FASTQ file format: expecting description line to start with @ but it starts with %c", *line);     

    *desc = line+1; remaining--; // skip the '@' - its already included in the container prefix
    
    // get DESC line
    *seq = seg_get_next_line (vb, *desc, &remaining, desc_len, true, &has_13[0], "DESC");
    
    // get SEQ line
    *line3 = seg_get_next_item (vb, *seq, &remaining, GN_SEP, GN_IGNORE, GN_IGNORE, seq_len, NULL, &has_13[1], "SEQ");
    ASSSEG (remaining && (*line3)[0] == '+', *line3, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it starts with a '%c'", (*line3)[0]);
    (*line3)++; // skip '+';

    // get LINE3
    *qual = seg_get_next_line (vb, *line3, &remaining, line3_len, true, &has_13[2], "LINE3");
    
    // get QUAL line
    rom after = seg_get_next_item (vb, *qual, &remaining, GN_SEP, GN_FORBIDEN, GN_IGNORE, qual_len, NULL, &has_13[3], "QUAL"); 

    ASSSEG (*qual_len == *seq_len, *qual, "Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ = %.*s\nQUAL= %.*s",
            *seq_len, *qual_len, STRf(*seq), STRf(*qual));

    ASSSEG (str_is_in_range (STRa(*qual), 33, 126), *qual, "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", STRf(*qual));

    return after;
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
rom fastq_seg_txt_line (VBlockFASTQP vb, rom line_start, uint32_t remaining, bool *has_13)     // index in vb->txt_data where this line starts
{
    // Split read to to desc, seq, line3 and qual (excluding the '@' and '+')
    STR(desc); STR (seq); STR(line3); STR(qual);
    bool line_has_13[4] = {};
    rom after = fastq_seg_get_lines (vb, line_start, remaining, pSTRa(desc), pSTRa(seq), pSTRa(line3), pSTRa(qual), line_has_13);
    
    // case segconf: gather data
    if (segconf.running) {
        segconf.longest_seq_len = MAX_(seq_len, segconf.longest_seq_len);
        
        segconf_update_qual (STRa(qual)); // get stats on qual scores
        
        if (vb->line_i==0) qname_segconf_discover_flavor (VB, FASTQ_DESC, STRa(desc)); // discover the original TECH, even when optimize_DESC
    }

    // set dl fields, consumed by fastq_zip_qual/seq
    ZipDataLineFASTQ *dl = DATA_LINE (vb->line_i);
    dl->seq = TXTWORD(seq); 
    dl->qual_index = BNUMtxt (qual);

    unsigned qname_len = strcspn (desc, " \t\r\n"); 

    if (kraken_is_loaded) 
        fastq_seg_kraken_tax_id (vb, desc, qname_len);

    // case --optimize_QUAL: optimize in place, must be done before fastq_seg_deep calculates a hash
    if (flag.optimize_QUAL) optimize_phred_quality_string ((char *)STRa(qual));

    dl->monobase = str_is_monochar (STRa(seq)); // set dl->monobase

    // case --deep: compare DESC, SEQ and QUAL to the SAM/BAM data
    bool deep_desc=false, deep_seq=false, deep_qual=false;
    if (flag.deep || flag.debug_deep == 2)
        fastq_seg_deep (vb, dl, desc, qname_len, STRa(seq), STRa(qual), &deep_desc, &deep_seq, &deep_qual);

    fastq_seg_DESC (vb, STRa(desc), qname_len, deep_desc);

    fastq_seg_SEQ (vb, dl, STRa(seq), deep_seq);

    fastq_seg_LINE3 (vb, STRa(line3), STRa(desc)); // note: LINE3 can be either empty or a copy of the desc

    fastq_seg_QUAL (vb, dl, qual_len, deep_qual);

    // 4 end of lines. note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
    for (int i=0; i < 4; i++) {
        *has_13 = line_has_13[i]; 
        SEG_EOL (i==0 ? FASTQ_E1L : FASTQ_E2L, true);
    }

    return after;
}

//-----------------
// PIZ stuff
//-----------------

// PIZ: piz_process_recon callback: called by the main thread, in the order of VBs
void fastq_piz_process_recon (VBlockP vb)
{
    if (flag.collect_coverage)    
        coverage_add_one_vb (vb);
    
    if (flag.reading_kraken)
        kraken_piz_handover_data (vb);
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fastq_piz_is_skip_section)
{
    if (st != SEC_LOCAL && st != SEC_B250 && st != SEC_DICT) return false; // we only skip contexts

    // case: this is pair-2 loading pair-1 sections
    if (purpose == SKIP_PURPOSE_PREPROC)  
        return !(dict_id_is_fastq_desc_sf (dict_id) || 
                 dict_id_is_in (dict_id, _FASTQ_DESC, _FASTQ_GPOS, _FASTQ_STRAND, DICT_ID_NONE) ||
                 (dict_id.num == _FASTQ_SQBITMAP && st != SEC_LOCAL));

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast: skip all items but DESC and E1L (except if we need them for --grep)
    if (flag.header_only_fast && !flag.grep &&
        dict_id_is_in (dict_id, _FASTQ_E2L, _FASTQ_TAXID, _FASTQ_DEBUG_LINES,
                       _FASTQ_SQBITMAP, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_STRAND,
                       _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, 
                       _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE))
        return true;

    if (flag.seq_only && !flag.grep && 
        (   dict_id_is_in (dict_id, _FASTQ_E1L, _FASTQ_DESC, _FASTQ_TAXID, _FASTQ_DEBUG_LINES, _FASTQ_LINE3,
                           _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE)  
         || dict_id_is_fastq_desc_sf(dict_id))) // we don't need the DESC line 
        return true;

    if (flag.qual_only && !flag.grep &&
        (   dict_id_is_in (dict_id, _FASTQ_E1L, _FASTQ_DESC, _FASTQ_TAXID, _FASTQ_DEBUG_LINES, DICT_ID_NONE) 
         || dict_id_is_fastq_desc_sf(dict_id) 
         || (!flag.bases && dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE)))) // we don't need the SEQ line 
        return true;

    // if we're doing --sex/coverage, we only need TOPLEVEL, FASTQ_SQBITMAP and GPOS
    if (flag.collect_coverage && 
        (   dict_id_is_in (dict_id, _FASTQ_DESC, _FASTQ_TAXID, _FASTQ_DEBUG_LINES, _FASTQ_LINE3,
                           _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE)
         || dict_id_is_fastq_desc_sf(dict_id) 
         || (!flag.bases && dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE))))
        return true;

    // no need for the TAXID data if user didn't specify --taxid
    if (flag.kraken_taxid==TAXID_NONE && dict_id.num == _FASTQ_TAXID)
        return true;

    // if --count, we only need TOPLEVEL and the fields needed for the available filters (--taxid, --kraken, --grep, --bases)
    if (flag.count && !flag.grep && sections_has_dict_id (st) &&
         (dict_id.num != _FASTQ_TOPLEVEL && 
         (dict_id.num != _FASTQ_TAXID    || flag.kraken_taxid == TAXID_NONE) && 
         (!flag.bases || !dict_id_is_in (dict_id, _FASTQ_SQBITMAP, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE)) &&
         (dict_id.num != _FASTQ_DESC     || !kraken_is_loaded) && 
         (!dict_id_is_fastq_desc_sf(dict_id) || !kraken_is_loaded))) 
        return true;

    return false;
}

// inspects z_file flags and if needed reads additional data, and returns true if the z_file consists of FASTQs compressed with --pair
bool fastq_piz_is_paired (void)
{
    if (z_file->data_type != DT_FASTQ || 
        (VER(15) && flag.deep)        ||          // --deep files are not paired
        z_file->num_txt_files == 1) return false; // note: %2 and not ==1 for back-comp for concatenated files up to v13

    // this is a FASTQ genozip file. Now we can check dts_paired
    if (z_file->z_flags.dts_paired) return true;  

    // dts_paired is not set. This flag was introduced in 9.0.13 - if file is compressed with genozip version 10+, then for sure the file is not paired
    if (VER(10)) return false;
    
    // for v8, and v9 up to 9.0.12 it is paired iff user is tell us explicitly that this is a paired file
    z_file->z_flags.dts_paired = flag.undocumented_dts_paired;
    return flag.undocumented_dts_paired; 
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat) and each toplevel item
CONTAINER_FILTER_FUNC (fastq_piz_filter)
{
    if (dict_id.num == _FASTQ_TOPLEVEL) {
        
        // case: --header-only: dont show items 2+. note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        if (flag.header_only_fast && !flag.grep && item >= 2) 
            return false; // don't reconstruct this item (non-header textual)

        // note: for seq_only and qual_only we show a newline from E2L, but we don't consume the other two newlines produced by this
        // read, so we are going to show a newline from another line. only a potenial but highly unlikely problem if some lines have \n and some \r\n.
        else if (flag.seq_only && !flag.grep && item != -1 && item != 2 && item != 3) 
            return false; 

        else if (flag.qual_only && !flag.grep && item != -1 && item != 2 && item != 6 && item != 7) 
            return false; 
    }

    return true; // reconstruct
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat)
CONTAINER_CALLBACK (fastq_piz_container_cb)
{
    // --taxid: filter out by Kraken taxid 
    if (flag.kraken_taxid && is_top_level) {
        
        if (!kraken_is_loaded && !kraken_is_included_stored (vb, FASTQ_TAXID, false))
            vb->drop_curr_line = "taxid";
        
        else if (kraken_is_loaded && flag.kraken_taxid != TAXID_NONE) {
            rom qname = last_txt (vb, FASTQ_DESC) + !prefixes_len; // +1 to skip the "@", if '@' is in DESC and not in prefixes (for files up to version 12.0.13)
            unsigned qname_len = strcspn (qname, " \t\n\r");

            if (!kraken_is_included_loaded (vb, STRa(qname))) 
                vb->drop_curr_line = "taxid";
        }
    }

    // --bases
    if (flag.bases && is_top_level && !vb->drop_curr_line &&
        !iupac_is_included_ascii (last_txt (vb, FASTQ_SQBITMAP), vb->last_txt_len (FASTQ_SQBITMAP)))
        vb->drop_curr_line = "bases";
}

// Used in pair-2: reconstruct the parallel b250 snip from pair-1
SPECIAL_RECONSTRUCTOR (fastq_special_mate_lookup)
{
    ASSPIZ (ctx->pair_b250, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2. %s%s", ctx->tag_name,
            !VER(10) ? "If this file was compressed with Genozip version 9.0.12 or older use the --dts_paired command line option. You can see the Genozip version used to compress it with 'genocat -w '" : "",  
            !VER(10) ? z_name : "");
            
    ctx_get_next_snip (vb, ctx, true, pSTRa(snip));

    if (ctx->did_i == FASTQ_SQBITMAP && VER(14))
        ctx->pair1_is_aligned = (snip_len && *snip == SNIP_LOOKUP) ? PAIR1_ALIGNED : PAIR1_NOT_ALIGNED;

    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE /* we can't cache pair items */, STRa(snip), reconstruct); // might include delta etc - works because in --pair, ALL the snips in a context are FASTQ_SPECIAL_mate_lookup

    return NO_NEW_VALUE; // last_value already set (if needed) in reconstruct_one_snip
}

void fastq_reset_line (VBlockP vb)
{
    CTX(FASTQ_SQBITMAP)->pair1_is_aligned = VER(14) ? PAIR1_ALIGNED_UNKNOWN 
                                                    : PAIR1_ALIGNED; // up to v13, all lines had alignment data, even if unmapped
}

//---------------------------------------------------------
// FASTQ-specific fields in genozip header
//---------------------------------------------------------

void fastq_zip_genozip_header (SectionHeaderGenozipHeader *header)
{
    header->fastq.segconf_seq_len_dict_id = segconf.qname_seq_len_dict_id; // v14
}

void fastq_piz_genozip_header (const SectionHeaderGenozipHeader *header)
{
    if (VER(14)) segconf.qname_seq_len_dict_id = header->fastq.segconf_seq_len_dict_id; 
}
