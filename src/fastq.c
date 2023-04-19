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
#include "arch.h"
#include "filename.h"
#include "aligner.h"
#include "zfile.h"

#define dict_id_is_fastq_qname_sf dict_id_is_type_1
#define dict_id_is_fastq_aux      dict_id_is_type_2
#define dict_id_fastq_qname_sf dict_id_type_1

sSTRl (copy_qname_snip, 30);

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

void fastq_vb_release_vb (VBlockFASTQP vb)
{
    vb->pair_num_lines = vb->pair_vb_i = vb->optimized_desc_len = vb->pair_txt_data_len = 0;
    vb->first_line = 0;
    vb->has_extra = 0;
    memset (vb->item_filter, 0, sizeof(vb->item_filter));
    memset (vb->deep_stats,  0, sizeof(vb->deep_stats));
    FREE (vb->optimized_desc);
}

//-----------------------
// TXTFILE stuff
//-----------------------

static inline bool is_valid_read (rom t[4],      // for textual lines
                                  uint32_t l[4]) // their lengths
{
   return l[0] >= 2 && t[0][0] == '@' &&  // DESC line starts with a '@' and is of length at least 2
          l[1] >= 1 && l[1] == l[3]   && // QUAL and SEQ have the same length, which is at least 1
          l[2] >= 1 && t[2][0] == '+';   // THIRD line starts with '+'
}

// detect if a GENERIC file is actually a FASTQ - we call based on first character being '@', line 3 starting with '+', and line 1 and 3 being the same length
bool is_fastq (STRp(header), bool *need_more)
{
    if (!header_len || header[0] != '@' || !str_is_printable (STRa(header))) return false; // fail fast

    #define NUM_TEST_READS 3
    str_split_by_lines (header, header_len, 4 * NUM_TEST_READS);
    n_lines = (n_lines / 4) * 4; // round to whole reads

    if (!n_lines) {
        *need_more = true; // we can't tell yet - need more data
        return false;
    }

    for (int i=0; i < n_lines; i += 4)
        if (!is_valid_read (&lines[i], &line_lens[i])) return false; // SEQ and QUAL lines are of equal length

    return true;
}

bool is_fastq_pair_2 (VBlockP vb)
{
    return VB_DT(FASTQ) && VB_FASTQ->pair_vb_i > 0;
}

// "pair assisted" is a type pairing in which R1 data is loaded to ctx->localR1/b250R1 and R2 consults with it in seg/recon.
// An R2 section that is pair-assisted will have FlagsCtx.paired set (note: in v14, due to a bug, R2.SQBITMAP.b250 should have but failed to set this flag)
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
           DNUM(AUX) || DNUM(AUX_LENGTH) || DNUM(E1L) || DNUM (E2L) || DNUM(TAXID) || DNUM(TOPLEVEL);
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fastq_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i_out /* in/out */)
{    
    ASSERT (*i_out >= 0 && *i_out < Ltxt, "*i=%d is out of range [0,%u]", *i_out, Ltxt);

    rom nl[8];     // newline pointers: nl[0] is the first from the end
    uint32_t l[8]; // lengths of segments excluding \n and \r: l[1] is the segment that starts at nl[1]+1 until nl[0]-1 (or nl[0]-2 if there is a \r). l[0] is not used.

    // search backwards for up to 8 newlines (best case: \nD\nS\nT\nQ\n ; worst case: \nD1\nS1\nT1\nQ1\nD2\nS2\nT2\nq2 (q2 is partial Q2))
    int n=0;
    for (rom c=Btxt (*i_out), first_c=Btxt (first_i) ; c >= first_c-1/*one beyond*/ && n < 8; c--) 
        if (c == (first_c-1) || *c == '\n') { // we consider character before the start to also be a "virtual newline"
            nl[n] = c;
            if (n) l[n] = ((nl[n-1]) - (nl[n-1][-1] == '\r')) - (nl[n] + 1);

            // 5th+ newline - test for valid read
            if (n >= 4 && is_valid_read ((rom[]){ nl[n]+1, nl[n-1]+1, nl[n-2]+1, nl[n-3]+1 }, (uint32_t[]){ l[n], l[n-1], l[n-2], l[n-3]})) {
                *i_out = BNUMtxt (nl[n-4]); // the final newline of this read (everything beyond is "unconsumed" and moved to the next VB) 
                return BLSTtxt - nl[n-4];   // number of "unconsumed" characters remaining in txt_data after the last \n of this read
            }
        
            n++;
        }

    ASSINP (n < 8, "Examined 7 textual lines and could not find a valid FASTQ read, it appears that this is not a valid FASTQ file. Data examined:\n%.*s",
            (int)(nl[0] - nl[7]), nl[7] + 1); // 7 lines and their newlines

    // case: the data provided has less than 8 newlines, and within it we didn't find a read. need more data.
    *i_out = (int32_t)first_i - 1; // next index to test - one before first_i
    return -1;
}

// called by txtfile_read_vblock when reading the 2nd file in a fastq pair - counts the number of fastq "lines" (each being 4 textual lines),
// comparing to the number of lines in the first file of the pair
// returns true if we have at least as much as needed, and sets unconsumed_len to the amount of excess characters read
// returns false is we don't yet have pair_1_num_lines lines - we need to read more
bool fastq_txtfile_have_enough_lines (VBlockP vb_, uint32_t *unconsumed_len, 
                                      uint32_t *my_lines, VBIType *pair_vb_i, uint32_t *pair_lines, uint32_t *pair_txt_data_len) // out - only set in case of failure
{
    START_TIMER;

    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
    ASSERTNOTZERO (vb->pair_num_lines, "");

    rom next  = B1STtxt;
    rom after = BAFTtxt;

    uint32_t pair_num_txt_lines = vb->pair_num_lines * 4, line_i;
    for (line_i=0; line_i < pair_num_txt_lines; line_i++) {
        if (!(next = memchr (next, '\n', after - next))) {
            *my_lines          = line_i;
            *pair_vb_i         = vb->pair_vb_i;
            *pair_lines        = pair_num_txt_lines;
            *pair_txt_data_len = vb->pair_txt_data_len;
            return false;
        }
        next++; // skip newline
    }

    vb->lines.len32 = line_i / 4;
    *unconsumed_len = after - next;

    COPY_TIMER (fastq_txtfile_have_enough_lines);
    return true;
}

void fastq_zip_set_txt_header_flags (struct FlagsTxtHeader *f)
{
    f->pair = flag.pair;
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

// called by main thread, as VBs complete (might be out-of-order)
void fastq_zip_after_compute (VBlockP vb)
{
    if (flag.deep)
        fastq_deep_zip_after_compute (VB_FASTQ);
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

    if (segconf.r1_or_r2)
        TIP ("Using --pair to compress paired-end FASTQs can reduce the compressed file's size by 10%%. E.g.:\n"
             "%s --reference %s --pair %s %s%s%s\n",
             arch_get_argv0(), ref_get_filename (gref), txt_name, (is_cd ? "" : dir_name), (is_cd ? "" : "/"), ent->d_name);

    closedir(dir);    
}

// called by main thread at the beginning of zipping this txt file  
void fastq_zip_initialize (void)
{
    if (!flag.reference && !txt_file->redirected && !flag.multiseq)
        TIP ("Compressing a FASTQ file using a reference file can reduce the compressed file's size by 20%%-60%%.\n"
             "Use: \"%s --reference <ref-file> %s\". ref-file may be a FASTA file or a .ref.genozip file.\n", 
             arch_get_argv0(), txt_file->name);

    else
        fastq_tip_if_should_be_pair();
    
    // reset lcodec for STRAND and GPOS, as these may change between PAIR_1 and PAIR_2 files
    ZCTX(FASTQ_STRAND)->lcodec = CODEC_UNKNOWN;
    ZCTX(FASTQ_GPOS  )->lcodec = CODEC_UNKNOWN;

    seg_prepare_snip_other (SNIP_COPY, _FASTQ_QNAME, 0, 0, copy_qname_snip);

    // with REF_EXTERNAL, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (IS_REF_EXTERNAL && z_file->num_txts_so_far == 1) // single file, or first of pair (and never Deep)
        ctx_populate_zf_ctx_from_contigs (gref, FASTQ_CONTIG, ref_get_ctgs (gref)); 

    qname_zip_initialize();

    if (flag.deep) fastq_deep_zip_initialize();
}

// called by main thread after each txt file compressing is done
void fastq_zip_finalize (bool is_last_user_txt_file)
{
    if (is_last_user_txt_file && flag.deep && flag.show_deep)
        fastq_zip_deep_show_entries_stats();

    // alternate pair - also works in case of multiple pairs of files
    flag.pair = (flag.pair == PAIR_R1) ? PAIR_R2
              : (flag.pair == PAIR_R2) ? PAIR_R1
              :                              NOT_PAIRED_END;
}

// called by Compute thread at the beginning of this VB
void fastq_seg_initialize (VBlockFASTQP vb)
{
    START_TIMER;
    declare_seq_contexts;

    vb->has_extra = segconf.has_extra; // VB-private copy

    if (!flag.deep)
        CTX(FASTQ_CONTIG)->flags.store = STORE_INDEX; // since v12

    if (flag.aligner_available) {
        strand_ctx->ltype       = LT_BITMAP;
        gpos_ctx->ltype         = LT_UINT32;
        gpos_ctx->flags.store   = STORE_INT;
        gpos_d_ctx->ltype       = LT_INT16;
        nonref_ctx->no_callback = true; // when using aligner, nonref data is in local. without aligner, the entire SEQ is segged into nonref using a callback
        
        bitmap_ctx->ltype     = LT_BITMAP; // implies no_stons
        bitmap_ctx->local_always = true;

        buf_alloc (vb, &bitmap_ctx->local, 1, Ltxt / 4, uint8_t, 0, "contexts->local"); 
        buf_alloc (vb, &strand_ctx->local, 0, roundup_bits2bytes64 (vb->lines.len), uint8_t, 0, "contexts->local"); 

        for (int i=0; i < 4; i++)
            buf_alloc (vb, &seqmis_ctx[i].local, 1, Ltxt / 128, char, 0, "contexts->local"); 

        buf_alloc (vb, &gpos_ctx->local, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->local"); 
    }

    if (!flag.multiseq)
        codec_acgt_comp_init (VB, FASTQ_NONREF);

    // initialize QUAL to LT_SEQUENCE, it might be changed later to LT_CODEC (eg domq, longr)
    ctx_set_ltype (VB, LT_SEQUENCE, FASTQ_QUAL, DID_EOL);

    if (flag.pair == PAIR_R1) 
        // cannot all_the_same with no b250 for PAIR_1 - SQBITMAP.b250 is tested in fastq_get_pair_1_gpos_strand
        // See defect 2023-02-11. We rely on this "no_drop_b250" in fastq_piz_get_pair2_is_forward 
        bitmap_ctx->no_drop_b250 = true; 
    
    else if (flag.pair == PAIR_R2) {
        ASSERT (vb->lines.len32 == vb->pair_num_lines, "in vb=%s (PAIR_R2): pair_num_lines=%u but lines.len=%u",
                VB_NAME, vb->pair_num_lines, vb->lines.len32);

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

    if (flag.pair == PAIR_R2)
        ctx_create_node (VB, FASTQ_SQBITMAP, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_mate_lookup }, 2);

    if (flag.optimize_DESC) 
        fastq_get_optimized_desc_read_name (vb);

    if (kraken_is_loaded) {
        CTX(FASTQ_TAXID)->flags.store    = STORE_INT;
        CTX(FASTQ_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTQ_TAXID)->counts_section = true; 
    }

    if (flag.deep)
        fastq_deep_seg_initialize (vb);

    if (segconf.seq_len_dict_id.num)
        ctx_set_store (VB, STORE_INT, ECTX(segconf.seq_len_dict_id)->did_i, DID_EOL);

    // consolidate stats for --stats
    ctx_consolidate_stats (VB, FASTQ_SQBITMAP, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_GPOS_DELTA, FASTQ_STRAND, FASTQ_SEQMIS_A, FASTQ_SEQMIS_C, FASTQ_SEQMIS_G, FASTQ_SEQMIS_T, DID_EOL);
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
        .nitems_lo      = 16, // initial before removals
        .items          = {
            { .dict_id = { _FASTQ_DEEP     }, .separator = { CI0_INVISIBLE } }, // 0:  note: will be removed if not --deep. added v15
            { .dict_id = { _FASTQ_QNAME    },                                }, // 1:
            { .dict_id = { _FASTQ_QNAME2   },                                }, // 2:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_EXTRA    },                                }, // 3:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_AUX      },                                }, // 4:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_E1L      },                                }, // 5:  note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
            { .dict_id = { _FASTQ_SQBITMAP },                                }, // 6:  
            { .dict_id = { _FASTQ_E2L      },                                }, // 7:
            { .dict_id = { _FASTQ_LINE3    },                                }, // 8:  added in 12.0.14, before '+' was a separator of the previous E2L. v15: note: removed if not needed
            { .dict_id = { _FASTQ_QNAME2   },                                }, // 9:  note: removed if not needed. added v15
            { .dict_id = { _FASTQ_EXTRA    },                                }, // 10: note: removed if not needed. added v15
            { .dict_id = { _FASTQ_AUX      },                                }, // 11: note: removed if not needed. added v15
            { .dict_id = { _FASTQ_E2L      },                                }, // 12:
            { .dict_id = { _FASTQ_QUAL     },                                }, // 13:
            { .dict_id = { _FASTQ_E2L      },                                }, // 14:
            { segconf.seq_len_dict_id,        .separator = { CI0_INVISIBLE } }  // 15: reconstruct for seq_len in case QNAME and LINE3 are dropped with --seq-only or --qual-only. See logic in fastq_piz_initialize_item_filter
        }
    };

    ASSERT0 (top_level.nitems_lo <= sizeof (VB_FASTQ->item_filter), "length of VBlockFASTQ.item_filter needs increasing");

    // prefixes in this container were added in 12.0.14, before, '@' was part of DESC and '+' was a separator
    char prefixes[] = { CON_PX_SEP,        // 0:  initial
                        CON_PX_SEP,        // 1:  terminator for empty container-wide prefix
                        CON_PX_SEP,        // 2:  empty DEEP line prefix. note: will be removed below if not --deep
                        '@', CON_PX_SEP,   // 3:  QNAME prefix
                        ' ', CON_PX_SEP,   // 5:  QNAME2 prefix
                        ' ', CON_PX_SEP,   // 7:  EXTRA prefix
                        ' ', CON_PX_SEP,   // 9:  AUX prefix
                        CON_PX_SEP,        // 11: empty E1L line prefix
                        CON_PX_SEP,        // 12: empty SQBITMAP line prefix
                        CON_PX_SEP,        // 13: empty E2L line prefix
                        '+', CON_PX_SEP,   // 14: LINE3 prefix
                        ' ', CON_PX_SEP,   // 16: QNAME2 prefix
                        ' ', CON_PX_SEP,   // 18: EXTRA prefix
                        ' ', CON_PX_SEP,   // 20: AUX prefix
                        CON_PX_SEP,        // 22: E2L prefix
                        CON_PX_SEP,        // 23: seq_len_dict_id prefix
                        CON_PX_SEP      }; // 24: END of prefixes

    int prefixes_len = sizeof (prefixes);

    #define REMOVE(item_i,px_i,px_len) ({ \
        memmove (&top_level.items[item_i], &top_level.items[item_i+1], (top_level.nitems_lo - item_i - 1) * sizeof (ContainerItem)); \
        top_level.nitems_lo--; \
        memmove (&prefixes[px_i], &prefixes[px_i+px_len], sizeof (prefixes)-(px_i+px_len)); \
        prefixes_len -= px_len; })

    bool has_line3  = segconf.line3 != L3_EMPTY;
    bool has_qname2 = segconf.has_qname2  && !flag.optimize_DESC;
    bool has_extra  = VB_FASTQ->has_extra && !flag.optimize_DESC;
    bool has_aux    = segconf.has_aux     && !flag.optimize_DESC;
    bool desc_is_l3 = segconf.desc_is_l3; // whether the Description (QNAME + DESC + AUX) appears on line 1 or line 3 (note: if on both, the line 3 is just a copy snip from line 1)

    // remove unneeded container and prefix items - in reverse order
    if (!segconf.seq_len_dict_id.num) REMOVE (15, 23, 1);
    if (!desc_is_l3 || !has_aux)      REMOVE (11, 20, 2);
    if (!desc_is_l3 || !has_extra)    REMOVE (10, 18, 2);
    if (!desc_is_l3 || !has_qname2)   REMOVE (9,  16, 2);
    if (!has_line3)                   REMOVE (8,  15, 1); // removing CON_PX_SEP the '+' now becomes the prefix of E2L (as QNAME2, EXTRA, AUX have already been removed)
    if (desc_is_l3  || !has_aux)      REMOVE (4,  9,  2);
    if (desc_is_l3  || !has_extra)    REMOVE (3,  7,  2);
    if (desc_is_l3  || !has_qname2)   REMOVE (2,  5,  2);
    if (!flag.deep)                   REMOVE (0,  2,  1);

    container_seg (vb, CTX(FASTQ_TOPLEVEL), (ContainerP)&top_level, prefixes, prefixes_len, 0); // note: the '@', '+' and ' ' are accounted for in the QNAME, QNAME2/EXTRA/AUX and LINE3 fields respectively
}

bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _FASTQ_TOPLEVEL ||
           dict_id.num == _FASTQ_QNAME    ||
           dict_id.num == _FASTQ_TAXID    ||
           dict_id.num == _FASTQ_E1L      ||
           dict_id.num == _FASTQ_E2L;
}

// ZIP/PIZ main thread: called ahead of zip or piz a pair 2 vb - to read data we need from the previous pair 1 file
// returns true if successful, false if there isn't a vb with vb_i in the previous file
void fastq_read_pair_1_data (VBlockP vb_, VBIType pair_vb_i)
{
    START_TIMER;

    VBlockFASTQP vb = (VBlockFASTQP)vb_;
    uint64_t save_disk_so_far = z_file->disk_so_far;

    vb->pair_vb_i = pair_vb_i;

    Section sec = sections_vb_header (pair_vb_i);
    vb->pair_num_lines = sec->num_lines;

    if (flag.debug)  // use --debug to access - displays in errors in txtfile_read_vblock
        vb->pair_txt_data_len = BGEN32 (zfile_read_section_header (vb, sec, SEC_VB_HEADER).vb_header.recon_size_prim);
    
    // read into ctx->pair the data we need from our pair: QNAME,QNAME2,LINE3 and its components, GPOS and STRAND
    buf_alloc (vb, &vb->z_section_headers, MAX_DICTS * 2, 0, uint32_t, 0, "z_section_headers"); // indices into vb->z_data of section headers
        
    piz_read_all_ctxs (VB, &sec, true);
    
    file_seek (z_file, 0, SEEK_END, false); // restore
    z_file->disk_so_far = save_disk_so_far;

    COPY_TIMER (fastq_read_pair_1_data);
}

// main thread: after reading VB_HEADER and before reading local/b250 sections from z_file
void fastq_piz_before_read (VBlockP vb)
{
    if (writer_am_i_pair_2 (vb->vblock_i, &VB_FASTQ->pair_vb_i)) { // sets pair_vb_i if R2, leaves it 0 if R1
        
        // backward compatability: prior to V15, PIZ didn't rely on FlagCtx.paired, and it was not always
        // applied: SQBITMAP.b250 in v14 and GPOS.local (at least) in v14 incorrectly didn't set the flag
        if (!VER(15)) {
            CTX(FASTQ_GPOS)->pair_assist_type = CTX(FASTQ_STRAND)->pair_assist_type = SEC_LOCAL;
            if (VER(14)) CTX(FASTQ_SQBITMAP)->pair_assist_type = SEC_B250;
        }
    }
}

// main thread: called from piz_read_one_vb as DTP(piz_init_vb)
bool fastq_piz_init_vb (VBlockP vb, ConstSectionHeaderVbHeaderP header, uint32_t *txt_data_so_far_single_0_increment)
{
    // in case of this is a R2 of a paired fastq file, get the R1 data
    if (vb && VB_FASTQ->pair_vb_i > 0)
        fastq_read_pair_1_data (vb, VB_FASTQ->pair_vb_i);

    return true;
}

// seg taxid, or abort
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
                                pSTRp(qname), pSTRp(desc), pSTRp(seq), pSTRp(line3), pSTRp(qual), uint32_t *line1_len, bool *has_13) // out
{
    ASSSEG0 (*line != '\n', line, "Invalid FASTQ file format: unexpected newline");

    ASSSEG (*line == '@', line, "Invalid FASTQ file format: expecting description line to start with @ but it starts with %c", *line);     

    *qname = line+1; remaining--; // skip the '@' - its already included in the container prefix
    
    bool analyze = segconf.running && vb->line_i==0;

    // get QNAME + DESC line
    *seq = seg_get_next_line (VB, *qname, &remaining, line1_len, true, &has_13[0], "LINE1");
    
    // analyze Line1 (segconf)
    if (analyze) {
        *desc = memchr (*qname, ' ', *line1_len);

        // Discover the QNAME flavor. Also discovers the original TECH, even when optimize_DESC
        qname_segconf_discover_flavor (VB, QNAME1, *qname, (*desc ? *desc - *qname : *line1_len)); 

        if (*desc) fastq_segconf_analyze_DESC (vb, *desc + 1, *qname + *line1_len - *desc - 1);
    }

    if (segconf.has_desc && !segconf.desc_is_l3 &&
        (*desc = memchr (*qname, ' ', *line1_len))) {
        
        (*desc)++; // skip space
        *desc_len = *qname + *line1_len - *desc; 
        *qname_len = *line1_len - *desc_len - 1;
    }

    // entire line will be segged in seg_qname. hopefully it only contains a qname, but if not - it will be tokenized
    else {
        *qname_len = *line1_len;
        *desc_len = 0;
    }

    // get SEQ line
    *line3 = seg_get_next_item (VB, *seq, &remaining, GN_SEP, GN_IGNORE, GN_IGNORE, seq_len, NULL, &has_13[1], "SEQ");
    ASSSEG (remaining && (*line3)[0] == '+', *line3, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it starts with a '%c'", (*line3)[0]);
    (*line3)++; // skip '+';

    // get LINE3
    *qual = seg_get_next_line (VB, *line3, &remaining, line3_len, true, &has_13[2], "LINE3");

    // analyze Line3 (segconf)
    if (analyze) {
        if (*line3_len==0 || flag.optimize_DESC) 
            segconf.line3 = L3_EMPTY;

        else if (fastq_is_line3_copy_of_line1 (STRa(*qname), STRa(*line3), *desc_len))
            segconf.line3 = L3_COPY_LINE1;

        // test for e.g. Line1: "@3/1" Line3: "+SRR2982101.3 3 length=101"
        else if (!segconf.has_desc) {
            *desc = memchr (*line3, ' ', *line3_len);

            qname_segconf_discover_flavor (VB, QLINE3, *line3, *desc - *line3);

            if (*desc) fastq_segconf_analyze_DESC (vb, *desc + 1, *line3 + *line3_len - *desc - 1);

            segconf.desc_is_l3 = true; // DESC appears on line 3 rather than line 1
            segconf.line3 = L3_NCBI;
        }

        else 
            ABORT ("Unsupported FASTQ file format: expecting middle line to be a \"+\" with or without a copy of the description, or NCBI, but it is \"%.*s\"",
                    (*line3_len)+1, (*line3)-1);
    }
    
    if (segconf.has_desc && segconf.desc_is_l3) {
        *desc = memchr (*line3, ' ', *line3_len);
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
    rom after = seg_get_next_item (VB, *qual, &remaining, GN_SEP, GN_FORBIDEN, GN_IGNORE, qual_len, NULL, &has_13[3], "QUAL"); 

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
    STR0(qname); STR0(desc); STR0(seq); STR0(line3); STR0(qual);
    bool line_has_13[4] = {};
    uint32_t line1_len;
    rom after = fastq_seg_get_lines (vb, line_start, remaining, pSTRa(qname), pSTRa(desc), pSTRa(seq), pSTRa(line3), pSTRa(qual), 
                                     &line1_len, line_has_13);
    
    // case segconf: gather data
    if (segconf.running) 
        segconf.longest_seq_len = MAX_(seq_len, segconf.longest_seq_len);
        
    // set dl fields, consumed by fastq_zip_qual/seq
    ZipDataLineFASTQ *dl = DATA_LINE (vb->line_i);
    dl->seq = TXTWORD(seq); 
    dl->qual_index = BNUMtxt (qual);

    if (kraken_is_loaded) 
        fastq_seg_kraken_tax_id (vb, STRa(qname));

    // case --optimize_QUAL: optimize in place, must be done before fastq_seg_deep calculates a hash
    if (flag.optimize_QUAL) optimize_phred_quality_string ((char *)STRa(qual));

    dl->monochar = str_is_monochar (STRa(seq)); // set dl->monochar

    // case --deep: compare DESC, SEQ and QUAL to the SAM/BAM data
    bool deep_qname=false, deep_seq=false, deep_qual=false;
    if (flag.deep || flag.show_deep == 2)
        fastq_seg_deep (vb, dl, STRa(qname), STRa(seq), STRa(qual), &deep_qname, &deep_seq, &deep_qual);

    fastq_seg_QNAME (vb, STRa(qname), line1_len, deep_qname);

    // seg DESC (i.e. QNAME2 + EXTRA + AUX), it might have come either from line1 or from line3
    if (segconf.has_desc && !flag.optimize_DESC)
        fastq_seg_DESC (vb, STRa(desc));

    fastq_seg_LINE3 (vb, STRa(line3), STRa(qname), STRa(desc)); // before SEQ, in case it as a segconf_seq_len_dict_id field

    fastq_seg_SEQ (vb, dl, STRa(seq), deep_seq);

    fastq_seg_QUAL (vb, dl, STRa(qual), deep_qual);

    // 4 end of lines. note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
    for (int i=0; i < 4; i++) {
        *has_13 = line_has_13[i]; 
        SEG_EOL (i==0 ? FASTQ_E1L : FASTQ_E2L, true);
    }

    // if seq_len_dict_id was detected in segconf line 0, it must appear in all lines
    ASSERT (!segconf.seq_len_dict_id.num || ctx_has_value_in_line (VB, segconf.seq_len_dict_id, NULL), 
            "Line missing length component (ctx=%s qname=\"%.*s\")", ECTX(segconf.seq_len_dict_id)->tag_name, STRf(qname));

    return after;
}

//-----------------
// PIZ stuff
//-----------------

// PIZ: main thread: piz_process_recon callback: usually called in order of VBs, but out-of-order if --test with no writer
void fastq_piz_process_recon (VBlockP vb)
{
    if (flag.collect_coverage)    
        coverage_add_one_vb (vb);
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fastq_piz_is_skip_section)
{
    if (!ST(LOCAL) && !ST(B250) && !ST(DICT)) return false; // we only skip contexts

    if (dict_id.num == _FASTA_TAXID && flag.kraken_taxid) return false; // TAXID normally skipped in --seq-only etc, but not if filtering for taxid

    #define LINE3_dicts _FASTQ_LINE3,  _FASTQ_T0HIRD, _FASTQ_T1HIRD, _FASTQ_T2HIRD, _FASTQ_T3HIRD, _FASTQ_T4HIRD, _FASTQ_T5HIRD, \
                        _FASTQ_T6HIRD, _FASTQ_T7HIRD, _FASTQ_T8HIRD, _FASTQ_T9HIRD, _FASTQ_TAHIRD, _FASTQ_TBHIRD, _FASTQ_TmHIRD /* just line3, not all qnames */

    // we always keep the seq_len context, if we have one
    if (dict_id.num && dict_id.num == segconf.seq_len_dict_id.num) return false;

    // note that flags_update_piz_one_z_file rewrites --header-only as flag.header_only_fast: skip all items but DESC and E1L (except if we need them for --grep)
    if (flag.header_only_fast && 
        dict_id_is_in (dict_id, _FASTQ_E2L, _FASTQ_DEBUG_LINES, LINE3_dicts,
                       _FASTQ_SQBITMAP, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_GPOS_DELTA, _FASTQ_STRAND,
                       _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, 
                       _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE))
        return true;

    if (flag.seq_only && 
        (   dict_id_is_in (dict_id, _FASTQ_E1L, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_EXTRA, _FASTQ_TAXID, _FASTQ_DEBUG_LINES, _FASTQ_LINE3,
                           _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE)  
         || dict_id_is_fastq_qname_sf(dict_id))) // we don't need the DESC line 
        return true;

    // note: we need SQBITMAP to extract seq_len. We even needs its local, bc otherwise SNIP_LOOKUP won't work.
    if (flag.qual_only && 
        (   dict_id_is_in (dict_id, _FASTQ_E1L, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_EXTRA, _FASTQ_TAXID, _FASTQ_DEBUG_LINES, DICT_ID_NONE) 
         || dict_id_is_fastq_qname_sf(dict_id) 
         || dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_GPOS_DELTA, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE))) // we don't need the SEQ line 
        return true;

    // if we're doing --sex/coverage, we only need TOPLEVEL, FASTQ_SQBITMAP and GPOS
    if (flag.collect_coverage && 
        (   dict_id_is_in (dict_id, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_EXTRA, _FASTQ_TAXID, _FASTQ_DEBUG_LINES,
                           _FASTQ_QUAL, _FASTQ_DOMQRUNS, _FASTQ_QUALMPLX, _FASTQ_DIVRQUAL, DICT_ID_NONE)
         || dict_id_is_fastq_qname_sf(dict_id) 
         || (!flag.bases && dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE))))
        return true;

    // note: we don't SKIP for --count with an additional filter. Logic is too complicated and buggy.
    if (flag.count && !flag.grep && sections_has_dict_id (st) && 
        !flag.bases && !flag.regions && !flag.kraken_taxid && !kraken_is_loaded &&
        dict_id.num != _FASTQ_TOPLEVEL) 
        return true;

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
    for (int i=0; i < con->nitems_lo && n_nl < 4/*safety*/; i++)
        if (con->items[i].did_i_small == FASTQ_E1L || con->items[i].did_i_small == FASTQ_E2L)
            nl[n_nl++] = i;

    ASSERT0 (n_nl == 4, "expecting 4 newlines in a FASTQ container");

    #define KEEP_ITEMS(first_item, last_item) \
        ({ memset (vb->item_filter, false, first_item); \
           memset (&vb->item_filter[last_item+1], false, con->nitems_lo - (last_item+1)); })

    // --header-only: keep items up to and including first newline
    if (flag.header_only_fast) KEEP_ITEMS (0, nl[0]);

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
    if (nl[3] != con->nitems_lo - 1) // can only happen since v15
        vb->item_filter[con->nitems_lo-1] = (flag.seq_only || flag.qual_only); 
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat) and each toplevel item
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

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat)
CONTAINER_CALLBACK (fastq_piz_container_cb)
{
    // --taxid: filter out by Kraken taxid 
    if (flag.kraken_taxid && is_top_level) {
        
        if (!kraken_is_loaded && !kraken_is_included_stored (vb, FASTQ_TAXID, !flag.seq_only && !flag.qual_only))
            vb->drop_curr_line = "taxid";
        
        else if (kraken_is_loaded) {
            STR(qname);

            if (VER(15)) 
                CTXlast (qname, CTX(FASTQ_QNAME));
            
            else { // up to 14 the entire line1 was segged in FASTQ_EXTRA
                qname = last_txt (vb, FASTQ_EXTRA) + !prefixes_len; // +1 to skip the "@", if '@' is in DESC and not in prefixes (for files up to version 12.0.13)
                qname_len = strcspn (qname, " \t\n\r");
            }

            if (!kraken_is_included_loaded (vb, STRa(qname))) 
                vb->drop_curr_line = "taxid";
        }
    }

    // --bases
    if (flag.bases && is_top_level && !vb->drop_curr_line &&
        !iupac_is_included_ascii (STRlst(FASTQ_SQBITMAP)))
        vb->drop_curr_line = "bases";
}

// Used in R2: used for pair-assisted b250 reconstruction. copy parallel b250 snip from R1. 
// Since v14, used for SQBITMAP. For files up to v14, also used for all QNAME subfield contexts
SPECIAL_RECONSTRUCTOR (fastq_special_mate_lookup)
{
    ASSPIZ (ctx->b250R1.len32, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2.", ctx->tag_name);
            
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
    header->fastq.segconf_seq_len_dict_id = segconf.seq_len_dict_id; // v14
}

void fastq_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header)
{
    if (VER(14)) segconf.seq_len_dict_id = header->fastq.segconf_seq_len_dict_id; 
}
