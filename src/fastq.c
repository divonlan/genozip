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

#define dict_id_is_fastq_qname_sf dict_id_is_type_1
#define dict_id_fastq_qname_sf dict_id_type_1

sSTRl (copy_qname_snip, 30);

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

void fastq_vb_release_vb (VBlockFASTQP vb)
{
    vb->pair_num_lines = vb->pair_vb_i = vb->optimized_desc_len = 0;
    vb->first_line = 0;
    vb->has_extra = 0;
    memset (vb->item_filter, 0, sizeof(vb->item_filter));
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
    if (!header_len || header[0] != '@' || !str_is_printable (STRa(header))) return false; // fail fast

    #define NUM_TEST_READS 3
    str_split_by_lines (header, header_len, 4 * NUM_TEST_READS);
    n_lines = (n_lines / 4) * 4; // round to whole reads

    if (!n_lines) {
        *need_more = true; // we can't tell yet - need more data
        return false;
    }

    for (int i=0; i < n_lines; i += 4)
        if (!line_lens[i+0] || lines[i+0][0] != '@' ||
            !line_lens[i+2] || lines[i+2][0] != '+' ||
            !line_lens[i+1] || line_lens[i+1] != line_lens[i+3]) return false; // SEQ and QUAL lines are of equal length

    return true;
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

    vb->lines.len = line_i / 4;
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

    seg_prepare_snip_other (SNIP_COPY, _FASTQ_QNAME, 0, 0, copy_qname_snip);

    // with REF_EXTERNAL, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (IS_REF_EXTERNAL && z_file->num_txts_so_far == 1 && !flag.deep) // first file
        ctx_populate_zf_ctx_from_contigs (gref, FASTQ_CONTIG, ref_get_ctgs (gref)); 

    qname_zip_initialize();
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

void fastq_initialize_ctx_for_pairing (VBlockP vb, Did did_i)
{
    CTX(did_i)->no_stons = true; // prevent singletons, so pair_1 and pair_2 are comparable based on b250 only
    
    if (flag.pair == PAIR_READ_2)
        ctx_create_node (vb, did_i, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_mate_lookup }, 2); // required by ctx_convert_generated_b250_to_mate_lookup
}

// called by Compute thread at the beginning of this VB
void fastq_seg_initialize (VBlockFASTQP vb)
{
    START_TIMER;

    vb->has_extra = segconf.has_extra; // VB-private copy

    if (!flag.deep)
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

        // cannot all_the_same with no b250 for PAIR_1 - SQBITMAP.b250 is tested in fastq_get_pair_1_gpos_strand
        // See defect 2023-02-11. We rely on this "no_drop_b250" in fastq_piz_get_pair2_is_forward 
        if (vb->comp_i == FQ_COMP_R1)
            sqbitmap_ctx->no_drop_b250 = true; 

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

    if (!segconf.running)  // note: if segconf.running, we initialize in qname_segconf_discover_flavor
        for (QType q=QNAME1; q < NUM_QTYPES; q++)
            qname_seg_initialize (VB, q, flag.pair); 

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

    if (segconf.seq_len_dict_id.num)
        ctx_set_store (VB, STORE_INT, ECTX(segconf.seq_len_dict_id)->did_i, DID_EOL);

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
bool fastq_read_pair_1_data (VBlockP vb_, uint32_t pair_vb_i, bool must_have)
{
    VBlockFASTQP vb = (VBlockFASTQP)vb_;
    uint64_t save_disk_so_far = z_file->disk_so_far;

    vb->pair_vb_i = pair_vb_i;

    Section sec = sections_vb_header (pair_vb_i, SOFT_FAIL);
    ASSERT (!must_have || sec, "Cannot find VB_HEADER of vb_i=%u section for pair 1", pair_vb_i);
    if (!sec) return false;

    vb->pair_num_lines = sec->num_lines;

    // read into ctx->pair the data we need from our pair: QNAME,QNAME2,LINE3 and its components, EXTRA, GPOS and STRAND
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

    dl->monobase = str_is_monochar (STRa(seq)); // set dl->monobase

    // case --deep: compare DESC, SEQ and QUAL to the SAM/BAM data
    bool deep_qname=false, deep_seq=false, deep_qual=false;
    if (flag.deep || flag.debug_deep == 2)
        fastq_seg_deep (vb, dl, desc, qname_len, STRa(seq), STRa(qual), &deep_qname, &deep_seq, &deep_qual);

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

// PIZ: piz_process_recon callback: called by the main thread, in the order of VBs
void fastq_piz_process_recon (VBlockP vb)
{
    if (flag.collect_coverage)    
        coverage_add_one_vb (vb);
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fastq_piz_is_skip_section)
{
    if (st != SEC_LOCAL && st != SEC_B250 && st != SEC_DICT) return false; // we only skip contexts

    if (dict_id.num == _FASTA_TAXID && flag.kraken_taxid) return false; // TAXID normally skipped in --seq-only etc, but not if filtering for taxid

    // case: this is pair-2 loading pair-1 sections
    if (purpose == SKIP_PURPOSE_PREPROC)  
        return !(dict_id_is_fastq_qname_sf (dict_id) || 
                 dict_id_is_in (dict_id, _FASTQ_QNAME, _FASTQ_QNAME2, _FASTQ_LINE3, _FASTQ_GPOS, _FASTQ_STRAND, DICT_ID_NONE) ||
                 (dict_id.num == _FASTQ_SQBITMAP && st != SEC_LOCAL));

    #define LINE3_dicts _FASTQ_LINE3,  _FASTQ_T0HIRD, _FASTQ_T1HIRD, _FASTQ_T2HIRD, _FASTQ_T3HIRD, _FASTQ_T4HIRD, _FASTQ_T5HIRD, \
                        _FASTQ_T6HIRD, _FASTQ_T7HIRD, _FASTQ_T8HIRD, _FASTQ_T9HIRD, _FASTQ_TAHIRD, _FASTQ_TBHIRD, _FASTQ_TmHIRD /* just line3, not all qnames */

    // we always keep the seq_len context, if we have one
    if (dict_id.num && dict_id.num == segconf.seq_len_dict_id.num) return false;

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast: skip all items but DESC and E1L (except if we need them for --grep)
    if (flag.header_only_fast && 
        dict_id_is_in (dict_id, _FASTQ_E2L, _FASTQ_DEBUG_LINES, LINE3_dicts,
                       _FASTQ_SQBITMAP, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_STRAND,
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
         || dict_id_is_in (dict_id, _FASTQ_NONREF, _FASTQ_NONREF_X, _FASTQ_GPOS, _FASTQ_STRAND, _FASTQ_SEQMIS_A, _FASTQ_SEQMIS_C, _FASTQ_SEQMIS_G, _FASTQ_SEQMIS_T, DICT_ID_NONE))) // we don't need the SEQ line 
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
        if (item == -1 && rep == 0)
            fastq_piz_initialize_item_filter (VB_FASTQ, con);
        
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
    header->fastq.segconf_seq_len_dict_id = segconf.seq_len_dict_id; // v14
}

void fastq_piz_genozip_header (const SectionHeaderGenozipHeader *header)
{
    if (VER(14)) segconf.seq_len_dict_id = header->fastq.segconf_seq_len_dict_id; 
}
