// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "fastq.h"
#include "vblock.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "strings.h"
#include "piz.h"
#include "optimize.h"
#include "dict_id.h"
#include "refhash.h"
#include "endianness.h"
#include "zfile.h"
#include "piz.h"
#include "buffer.h"
#include "codec.h"
#include "aligner.h"
#include "stats.h"
#include "reconstruct.h"
#include "coverage.h"
#include "writer.h"
#include "kraken.h"
#include "bases_filter.h"
#include "qname.h"

#define dict_id_is_fastq_desc_sf dict_id_is_type_1
#define dict_id_fastq_desc_sf dict_id_type_1

#define pair1_is_aligned ctx_specific // PIZ: used in SQBITMAP when reconstructing pair-2

typedef struct {
    uint32_t seq_data_start, qual_data_start; // start within vb->txt_data
    uint32_t seq_len;                         // length of SEQ and QUAL within vb->txt_data (they are identical per FASTQ spec)
} ZipDataLineFASTQ;

// Used by Seg
static char copy_desc_snip[30];
static unsigned copy_desc_snip_len;

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockFASTQ {
    VBLOCK_COMMON_FIELDS

    // pairing stuff - used if we are the 2nd file in the pair 
    uint32_t pair_vb_i;      // the equivalent vb_i in the first file, or 0 if this is the first file
    uint32_t pair_num_lines; // number of lines in the equivalent vb in the first file
    char *optimized_desc;    // base of desc in flag.optimize_DESC 
    uint32_t optimized_desc_len;
    uint64_t first_line;     // ZIP: used for optimize_DESC  
} VBlockFASTQ;

#define VB_FASTQ ((VBlockFASTQ *)vb)

#define DATA_LINE(i) B(ZipDataLineFASTQ, vb->lines, i)

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

void fastq_vb_release_vb (VBlockFASTQ *vb)
{
    vb->pair_num_lines = vb->pair_vb_i = vb->optimized_desc_len = 0;
    vb->first_line = 0;
    FREE (vb->optimized_desc);
}

void fastq_vb_destroy_vb (VBlockFASTQ *vb)
{
}

//-----------------------
// TXTFILE stuff
//-----------------------

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
    for (int32_t i=txt_i-1; i >= first_i; i--) {
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

    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;

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
static void fastq_get_optimized_desc_read_name (VBlockFASTQ *vb)
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

// called by main thread at the beginning of zipping this file
void fastq_zip_initialize (void)
{
    // reset lcodec for STRAND and GPOS, as these may change between PAIR_1 and PAIR_2 files
    ZCTX(FASTQ_STRAND)->lcodec = CODEC_UNKNOWN;
    ZCTX(FASTQ_GPOS  )->lcodec = CODEC_UNKNOWN;

    seg_prepare_snip_other (SNIP_COPY, _FASTQ_DESC, 0, 0, copy_desc_snip);

    // with REF_EXTERNAL, we don't know which chroms are seen (bc unlike REF_EXT_STORE, we don't use is_set), so
    // we just copy all reference contigs. this are not needed for decompression, just for --coverage/--sex/--idxstats
    if (flag.reference == REF_EXTERNAL && z_file->num_txts_so_far == 1) // first file
        ctx_populate_zf_ctx_from_contigs (gref, FASTQ_CONTIG, ref_get_ctgs (gref)); 

    qname_zip_initialize (FASTQ_DESC);
}

// called by zfile_compress_genozip_header to set FlagsGenozipHeader.dt_specific
bool fastq_zip_dts_flag (void)
{
    return flag.pair != NOT_PAIRED_END;
}

// called by Compute thread at the beginning of this VB
void fastq_seg_initialize (VBlockFASTQ *vb)
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
        codec_acgt_comp_init (VB);

    if (flag.pair == PAIR_READ_2) {

        ASSERT (vb->lines.len32 == vb->pair_num_lines, "in vb=%s (PAIR_READ_2): pair_num_lines=%u but lines.len=%u",
                VB_NAME, vb->pair_num_lines, (unsigned)vb->lines.len);

        gpos_ctx->pair_local = strand_ctx->pair_local = true;

        piz_uncompress_all_ctxs (VB);

        vb->z_data.len32 = 0; // we've finished reading the pair file z_data, next, we're going to write to z_data our compressed output
    }

    qname_seg_initialize (VB, FASTQ_DESC);

    if (vb->comp_i == FQ_COMP_R2)
        ctx_create_node (VB, FASTQ_SQBITMAP, (char[]){ SNIP_MATE_LOOKUP }, 1); // required by ctx_convert_generated_b250_to_mate_lookup

    if (flag.optimize_DESC) 
        fastq_get_optimized_desc_read_name (vb);

    if (kraken_is_loaded) {
        CTX(FASTQ_TAXID)->flags.store    = STORE_INT;
        CTX(FASTQ_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTQ_TAXID)->counts_section = true; 
    }

    // in --stats, consolidate stats into SQBITMAP
    stats_set_consolidation (VB, FASTQ_SQBITMAP, 8, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND, FASTQ_SEQMIS_A, FASTQ_SEQMIS_C, FASTQ_SEQMIS_G, FASTQ_SEQMIS_T);

    COPY_TIMER (seg_initialize);
}

void fastq_seg_finalize (VBlockP vb)
{
    // assign the QUAL codec
    codec_assign_best_qual_codec (vb, FASTQ_QUAL, fastq_zip_qual, false);

    // top level snip
    SmallContainer top_level = { 
        .repeats        = vb->lines.len,
        .is_toplevel    = true,
        .filter_items   = true,
        .filter_repeats = true,
        .callback       = true,
        .nitems_lo      = 8,
        .items          = { 
            { .dict_id = { _FASTQ_DESC },     },
            { .dict_id = { _FASTQ_E1L },      }, // note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
            { .dict_id = { _FASTQ_SQBITMAP }, },
            { .dict_id = { _FASTQ_E2L },      },
            { .dict_id = { _FASTQ_LINE3 },    }, // added in 12.0.14, before '+' was a separator of the previous E2L
            { .dict_id = { _FASTQ_E2L },      },
            { .dict_id = { _FASTQ_QUAL },     },
            { .dict_id = { _FASTQ_E2L },      } 
        }
    };

    // prefixes in this container were added in 12.0.14, before, '@' was part of DESC and '+' was a separator
    static char prefixes[] = { CON_PX_SEP,        // initial
                               CON_PX_SEP,        // terminator for empty container-wide prefix
                               '@', CON_PX_SEP,   // DESC prefix
                               CON_PX_SEP,        // empty E1L line prefix
                               CON_PX_SEP,        // empty SQBITMAP line prefix
                               CON_PX_SEP,        // empty E2L line prefix
                               '+', CON_PX_SEP,   // third line prefix
                               CON_PX_SEP      }; // END of prefixes

    container_seg (vb, CTX(FASTQ_TOPLEVEL), (ContainerP)&top_level, prefixes, sizeof (prefixes), 0); // note: the '@' and '+' are accounted for in the DESC and LINE3 fields respectively
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
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
    uint64_t save_disk_so_far = z_file->disk_so_far;

    vb->pair_vb_i = pair_vb_i;

    Section sec = sections_vb_header (pair_vb_i, true);
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
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
    uint32_t pair_vb_i=0;
    bool i_am_pair_2 = vb && writer_get_pair (vb->vblock_i, &pair_vb_i) == 2;

    // in case of this is a paired fastq file, get all the pair_1 data not already fetched for the grep above
    if (i_am_pair_2) 
        fastq_read_pair_1_data (vb_, pair_vb_i, true);

    return true;
}

static void fastq_seg_line3 (VBlockFASTQ *vb, STRp(line3), STRp(desc))
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
        seg_by_did_i (VB, "", 0, FASTQ_LINE3, 1); // account for the '+' (segged as a toplevel container prefix)
        vb->recon_size -= line3_len; 
    }
    
    else if (segconf.line3 == L3_EMPTY) {
        ASSSEG (!line3_len || flag.optimize_DESC, line3, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it is \"%.*s\"", line3_len+1, line3-1);
        seg_by_did_i (VB, "", 0, FASTQ_LINE3, 1); // account for the '+' (segged as a toplevel container prefix)
    }
    
    else if (segconf.line3 == L3_COPY_DESC) {
        ASSSEG (str_issame_(STRa(line3), STRa(desc)), line3, 
                "Invalid FASTQ file format: expecting middle line to be a \"+\" followed by a copy of the description line, but it is \"%.*s\"", line3_len+1, line3-1); 
        seg_by_did_i (VB, STRa(copy_desc_snip), FASTQ_LINE3, 1 + line3_len); // +1 - account for the '+' (segged as a toplevel container prefix)
    }

    else if (segconf.line3 == L3_QF) 
        ASSSEG (qname_seg_qf (VB, CTX(FASTQ_LINE3), segconf.line3_flavor, STRa(line3), false, 0), line3, 
                "Invalid FASTQ file format: expecting middle line to be a \"+\" followed by a \"%s\" flavor, but it is \"%.*s\"", qf_name(segconf.line3_flavor), line3_len+1, line3-1); 

    else 
        ASSSEG (false, line3, "Invalid FASTQ file format: expecting middle line to be a \"+\" with or without a copy of the description, but it is \"%.*s\"",
                line3_len+1, line3-1);
}

static void fastq_get_pair_1_gpos_strand (VBlockFASTQ *vb, PosType *pair_gpos, bool *pair_is_forward)
{
    ContextP bitmap_ctx = CTX(FASTQ_SQBITMAP);

    // case: we are pair-1 OR we are pair-2, but this line in pair-1 is not aligned
    if (!bitmap_ctx->pair_b250 || 
        ({ STR(snip); ctx_get_next_snip (VB, bitmap_ctx, true, pSTRa(snip)); *snip != SNIP_LOOKUP; })) { // note: SNIP_LOOKUP in case of aligned, SNIP_SPECIAL in case of unaligned

        *pair_gpos = NO_GPOS;
        *pair_is_forward = 0;
    }

    // case: we are pair-2, and the corresponding line in pair-1 is aligned: get its gpos and is_forward
    else {
        ContextP gpos_ctx   = CTX(FASTQ_GPOS);
        ContextP strand_ctx = CTX(FASTQ_STRAND);

        ASSERT (gpos_ctx->pair.next < gpos_ctx->pair.len32, "%s: not enough data GPOS.pair (len=%u)", LN_NAME, gpos_ctx->pair.len32); 

        ASSERT (gpos_ctx->pair.next < strand_ctx->pair.nbits, "%s: cannot get pair_1 STRAND bit because pair_1 strand bitarray has only %u bits",
                LN_NAME, (unsigned)strand_ctx->pair.nbits);

        *pair_gpos = (PosType)*B32 (gpos_ctx->pair, gpos_ctx->pair.next); 
        *pair_is_forward = bit_array_get ((BitArrayP)&strand_ctx->pair, gpos_ctx->pair.next); 
        gpos_ctx->pair.next++;
    }
}

static void fastq_seg_sequence (VBlockFASTQ *vb, STRp(seq))
{
    // get pair-1 gpos and is_forward, but only if we are pair-2 and the corresponding pair-1 line is aligned
    PosType pair_gpos; 
    bool pair_is_forward;
    fastq_get_pair_1_gpos_strand (vb, &pair_gpos, &pair_is_forward); // does nothing if we are pair-1
               
    // case: aligned - lookup from SQBITMAP
    MappingType aln_res;
    if (flag.aligner_available && 
        ((aln_res = aligner_seg_seq (VB, CTX(FASTQ_SQBITMAP), STRa(seq), true, (vb->comp_i == FQ_COMP_R2), pair_gpos, pair_is_forward)))) {
        
        SNIPi1 (SNIP_LOOKUP, (aln_res==MAPPING_PERFECT ? -(int64_t)seq_len : (int64_t)seq_len)); // express perfect alignment by passing a negative seq_len 
        seg_by_ctx (VB, STRa(snip), CTX(FASTQ_SQBITMAP), seq_len); 
    }

    // case: not aligned - just add data to NONREF
    else {
        if (flag.show_aligner && !segconf.running) iprintf ("%s: unaligned\n", LN_NAME);

        Context *nonref_ctx = CTX(FASTQ_NONREF);
        buf_alloc (VB, &nonref_ctx->local, seq_len + 3, vb->txt_data.len / 64, char, CTX_GROWTH, "contexts->local"); 
        buf_add (&nonref_ctx->local, seq, seq_len);

        // case: we don't need to consume pair-1 gpos (bc we are pair-1, or pair-1 was not aligned): look up from NONREF
        SNIPi2 (SNIP_SPECIAL, FASTQ_SPECIAL_unaligned_SEQ, seq_len);
        seg_by_ctx (VB, STRa(snip), CTX(FASTQ_SQBITMAP), 0); 
        CTX(FASTQ_NONREF)->txt_len += seq_len; // account for the txt data in NONREF
    }

}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
rom fastq_seg_txt_line (VBlockFASTQ *vb, rom line_start, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb->line_i);

    ASSSEG0 (*line_start != '\n', line_start, "Invalid FASTQ file format: unexpected newline");

    ASSSEG (*line_start == '@', line_start, "Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            *line_start);

    rom FASTQ_DESC_str = line_start + 1; // skip the '@' - its already included in the container prefix
    char separator;

    int32_t len = (int32_t)(BAFTtxt - FASTQ_DESC_str);
    
    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG <-- Illumina, See: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // @20A_V100002704L1C001R012000000/1 <-- BGI, see: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
    unsigned FASTQ_DESC_len;
    rom FASTQ_SEQ_str = seg_get_next_line (vb, FASTQ_DESC_str, &len, &FASTQ_DESC_len, true, has_13, "DESC");
 
    if (kraken_is_loaded) {
        unsigned qname_len = strcspn (FASTQ_DESC_str, " \t\r\n"); 

        unsigned taxid_found = kraken_seg_taxid (VB, FASTQ_TAXID, FASTQ_DESC_str, qname_len, false);

        // if not found tax id for this read, try again, perhaps removing /1 or /2
        if (!taxid_found) {
            if (qname_len > 2 && FASTQ_DESC_str[qname_len-2] == '/' &&
                (FASTQ_DESC_str[qname_len-1] == '1' || FASTQ_DESC_str[qname_len-1] == '2'))
                qname_len -= 2;

            kraken_seg_taxid (VB, FASTQ_TAXID, FASTQ_DESC_str, qname_len, true); // this fails if missing
        }
    }

    // if flag.optimize_DESC is on, we replace the description with "filename.line_i" 
    unsigned optimized_len = 0; 
    if (flag.optimize_DESC) {
        optimized_len  = vb->optimized_desc_len + str_int (vb->first_line + vb->line_i, &vb->optimized_desc[vb->optimized_desc_len]);   
        vb->recon_size -= FASTQ_DESC_len - optimized_len;
    }

    if (segconf.running && vb->line_i==0) 
        // discover the original TECH, even when optimize_DESC
        qname_segconf_discover_flavor (VB, FASTQ_DESC, STRd(FASTQ_DESC));
    
    // attempt to parse the desciption line as qname, and fallback to tokenizer no qname format matches
    qname_seg (VB, CTX(FASTQ_DESC), 
               flag.optimize_DESC ? vb->optimized_desc : FASTQ_DESC_str, 
               flag.optimize_DESC ? optimized_len      : FASTQ_DESC_len,
               1); // account for the '@' (segged as a toplevel container prefix)
    SEG_EOL (FASTQ_E1L, true);

    // SEQ - just get the whole line
    dl->seq_data_start = BNUMtxt (FASTQ_SEQ_str);
    rom FASTQ_LINE3_str = seg_get_next_item (vb, FASTQ_SEQ_str, &len, GN_SEP, GN_IGNORE, GN_IGNORE, &dl->seq_len, &separator, has_13, "SEQ");

    if (segconf.running) 
        segconf.longest_seq_len = MAX_(dl->seq_len, segconf.longest_seq_len);

    fastq_seg_sequence (vb, FASTQ_SEQ_str, dl->seq_len);

    SEG_EOL (FASTQ_E2L, true);

    // PLUS - next line is expected to be a "+", and be either empty or a copy of the DESC line (except for the '@' / '+' prefix)
    ASSSEG (*FASTQ_LINE3_str == '+', FASTQ_LINE3_str, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it starts with a '%c'", *FASTQ_LINE3_str);
    FASTQ_LINE3_str++; len--; // skip the prefix

    unsigned FASTQ_LINE3_len;
    rom FASTQ_QUAL_str = seg_get_next_line (vb, FASTQ_LINE3_str, &len, &FASTQ_LINE3_len, true, has_13, "LINE3");

    fastq_seg_line3 (vb, STRd(FASTQ_LINE3), STRd(FASTQ_DESC));

    SEG_EOL (FASTQ_E2L, true); // account for ascii-10

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = BNUMtxt (FASTQ_QUAL_str);

    unsigned FASTQ_QUAL_len;
    rom after = seg_get_next_item (vb, FASTQ_QUAL_str, &len, GN_SEP, GN_FORBIDEN, GN_IGNORE, &FASTQ_QUAL_len, &separator, has_13, "QUAL"); \

    CTX(FASTQ_QUAL)->local.len += dl->seq_len;
    CTX(FASTQ_QUAL)->txt_len   += dl->seq_len;

    // get stats on qual scores
    if (segconf.running)
        segconf_update_qual (STRd (FASTQ_QUAL));

    // End Of Line    
    SEG_EOL (FASTQ_E2L, true);

    ASSSEG (str_is_in_range (FASTQ_QUAL_str, FASTQ_QUAL_len, 33, 126), FASTQ_QUAL_str, "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", 
            FASTQ_QUAL_len, FASTQ_QUAL_str);

    ASSSEG (FASTQ_QUAL_len == dl->seq_len, FASTQ_QUAL_str, "Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ = %.*s\nQUAL= %.*s",
            dl->seq_len, FASTQ_QUAL_len, dl->seq_len, FASTQ_SEQ_str, FASTQ_QUAL_len, FASTQ_QUAL_str);
 
    return after;
}

// callback function for compress to get data of one line (called by codec_bz2_compress)
COMPRESSOR_CALLBACK (fastq_zip_qual) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len  = MIN_(dl->seq_len, maximum_size);
    
    if (!line_data) return; // only lengths were requested

    *line_data = Btxt (dl->qual_data_start);

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag.optimize_QUAL) optimize_phred_quality_string (STRa(*line_data));

    if (is_rev) *is_rev = 0;
}

COMPRESSOR_CALLBACK (fastq_zip_seq) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);
    *line_data_len = dl->seq_len;
    *line_data     = Btxt (dl->seq_data_start);
    if (is_rev) *is_rev = 0;
}

//-----------------
// PIZ stuff
//-----------------

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (fastq_piz_is_skip_section)
{
    if (st != SEC_LOCAL && st != SEC_B250 && st != SEC_DICT) return false; // we only skip contexts

    // case: this is pair-2 loading pair-1 sections
    if (purpose == SKIP_PURPOSE_PREPROC)  
        return !(dict_id_is_fastq_desc_sf (dict_id) || 
                 dict_id.num  == _FASTQ_DESC   ||
                 dict_id.num  == _FASTQ_GPOS   || 
                 dict_id.num  == _FASTQ_STRAND ||
                 (dict_id.num == _FASTQ_SQBITMAP && st != SEC_LOCAL));

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast: skip all items but DESC and E1L (except if we need them for --grep)
    if (flag.header_only_fast && !flag.grep &&
        (dict_id.num == _FASTQ_E2L    || dict_id.num == _FASTQ_SQBITMAP || 
         dict_id.num == _FASTQ_NONREF || dict_id.num == _FASTQ_NONREF_X || 
         dict_id.num == _FASTQ_GPOS   || dict_id.num == _FASTQ_STRAND   || 
         dict_id.num == _FASTQ_QUAL   || dict_id.num == _FASTQ_DOMQRUNS ))
        return true;

    if (flag.seq_only && !flag.grep &&
        (dict_id.num == _FASTQ_E1L    || dict_id.num == _FASTQ_DESC || dict_id_is_fastq_desc_sf(dict_id) || // we don't need the DESC line 
         dict_id.num == _FASTQ_QUAL   || dict_id.num == _FASTQ_DOMQRUNS )) // we don't need the QUAL line
        return true;

    if (flag.qual_only && !flag.grep &&
        (dict_id.num == _FASTQ_E1L || dict_id.num == _FASTQ_DESC || dict_id_is_fastq_desc_sf(dict_id) || // we don't need the DESC line 
         (!flag.bases && (dict_id.num == _FASTQ_NONREF || dict_id.num == _FASTQ_NONREF_X || dict_id.num == _FASTQ_GPOS || dict_id.num == _FASTQ_STRAND)))) // we don't need the SEQ line        return true;
        return true;

    // if we're doing --show-sex/coverage, we only need TOPLEVEL, FASTQ_SQBITMAP and GPOS
    if (flag.collect_coverage && 
        (dict_id.num == _FASTQ_DESC     || 
         dict_id.num == _FASTQ_QUAL     || 
         dict_id.num == _FASTQ_DOMQRUNS || 
         (dict_id.num == _FASTQ_STRAND   && !flag.bases)  ||
         (dict_id.num == _FASTQ_NONREF   && !flag.bases)  || 
         (dict_id.num == _FASTQ_NONREF_X && !flag.bases))) 
        return true;

    // no need for the TAXID data if user didn't specify --taxid
    if (flag.kraken_taxid==TAXID_NONE && dict_id.num == _FASTQ_TAXID)
        return true;

    // if --count, we only need TOPLEVEL and the fields needed for the available filters (--taxid, --kraken, --grep, --bases)
    if (flag.count && sections_has_dict_id (st) &&
         (dict_id.num != _FASTQ_TOPLEVEL && 
         (dict_id.num != _FASTQ_TAXID    || flag.kraken_taxid == TAXID_NONE) && 
         (dict_id.num != _FASTQ_DESC     || !kraken_is_loaded) && 
         (dict_id.num != _FASTQ_SQBITMAP || !flag.bases) && 
         (dict_id.num != _FASTQ_NONREF   || !flag.bases) && 
         (dict_id.num != _FASTQ_NONREF_X || !flag.bases) && 
         (dict_id.num != _FASTQ_GPOS     || !flag.bases) && 
         (dict_id.num != _FASTQ_STRAND   || !flag.bases) && 
         (!dict_id_is_fastq_desc_sf(dict_id) || !kraken_is_loaded))) 
        return true;

    return false;
}

// inspects z_file flags and if needed reads additional data, and returns true if the z_file consists of FASTQs compressed with --pair
bool fastq_piz_is_paired (void)
{
    if (z_file->data_type != DT_FASTQ || z_file->num_components % 2) return false; // quick check to avoid the need for zfile_is_paired is most cases that dts_paired is missing

    // this is a FASTQ genozip file. Now we can check dts_paired
    if (z_file->z_flags.dts_paired) return true;  

    // dts_paired is not set. This flag was introduced in 9.0.13 - if file is compressed with genozip version 10+, then for sure the file is not paired
    if (z_file->genozip_version >= 10) return false;
    
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
        else if (flag.seq_only && !flag.grep && item != 2 && item != 3) 
            return false; 

        else if (flag.qual_only && !flag.grep && item != 2 && item != 6 && item != 7) 
            return false; 
    }

    return true; // reconstruct
}

static void fastq_update_coverage_aligned (VBlockFASTQ *vb)
{
    PosType gpos;
    Context *gpos_ctx = CTX(FASTQ_GPOS);

    if (vb->comp_i == FQ_COMP_R1) 
        gpos = NEXTLOCAL (uint32_t, gpos_ctx);

    else { // pair-2
        reconstruct_from_ctx (VB, FASTQ_GPOS, 0, false); // calls fastq_special_PAIR2_GPOS
        gpos = gpos_ctx->last_value.i;
    }

    ASSPIZ0 (gpos != NO_GPOS, "expecting a GPOS, because sequence is aligned");

    WordIndex ref_index = ref_contig_get_by_gpos (gref, gpos, NULL);
    ASSPIZ0 (ref_index != WORD_INDEX_NONE, "expecting ref_index, because sequence is aligned");

    if (flag.show_coverage || flag.show_sex)
        *B64 (vb->coverage, ref_index) += vb->seq_len;

    if (flag.show_coverage || flag.idxstats)
        (*B64 (vb->read_count, ref_index))++;
}

// PIZ: aligned SEQ reconstruction - called by reconstructing FASTQ_SQBITMAP which is a LOOKUP
void fastq_recon_aligned_SEQ (VBlockP vb_, ContextP bitmap_ctx, STRp(seq_len_str), bool reconstruct)
{
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;

    // v14: perfect alignment is expressed by a negative seq_len
    bool perfect_alignment = (seq_len_str[0] == '-');
    if (perfect_alignment) { seq_len_str++; seq_len_str_len--; }

    ASSERT (str_get_int_range32 (STRa(seq_len_str), 0, 0x7fffffff, (int32_t*)&vb->seq_len), "could not parse integer \"%.*s\"", STRf(seq_len_str));

    if (vb->comp_i == FQ_COMP_R2 && z_file->genozip_version >= 14) { 
        STR(snip); 
        ctx_get_next_snip (VB, CTX(FASTQ_SQBITMAP), true, pSTRa(snip)); 
        bitmap_ctx->pair1_is_aligned = (*snip == SNIP_LOOKUP);  
    }

    // just update coverage
    if (flag.collect_coverage) 
        fastq_update_coverage_aligned (vb);

    // --qual-only: only set vb->seq_len without reconstructing
    else if (flag.qual_only) {}

    // normal reconstruction
    else 
        aligner_reconstruct_seq (vb_, bitmap_ctx, vb->seq_len, vb->comp_i == FQ_COMP_R2, perfect_alignment, reconstruct);
}

// PIZ: SEQ reconstruction - in case of unaligned sequence 
SPECIAL_RECONSTRUCTOR (fastq_special_unaligned_SEQ)
{
    // case we are pair-2: advance pair-1 SQBITMAP iterator, and if pair-1 is aligned - also its GPOS iterator
    if (vb->comp_i == FQ_COMP_R2 && 
        ({ STR(snip); ctx_get_next_snip (VB, CTX(FASTQ_SQBITMAP), true, pSTRa(snip)); *snip == SNIP_LOOKUP; }))   // note: SNIP_LOOKUP in case of aligned, SNIP_SPECIAL in case of unaligned
        CTX(FASTQ_GPOS)->pair.next++;

    // just update coverage (unaligned)
    if (flag.collect_coverage) {
        if (flag.show_coverage || flag.show_sex)
            *(BAFT64 (vb->coverage) - NUM_COVER_TYPES + CVR_UNMAPPED) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*(BAFT64 (vb->read_count) - NUM_COVER_TYPES + CVR_UNMAPPED))++;
    }

    else {
        if (flag.show_aligner) iprintf ("%s: unaligned\n", LN_NAME);

        reconstruct_from_local_sequence (vb, CTX(FASTQ_NONREF), STRa(snip), reconstruct);
    }

    return NO_NEW_VALUE;
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

//------------------------
// Pair-2 GPOS (ZIP & PIZ)
//------------------------

void fastq_seg_pair2_gpos (VBlockP vb, PosType pair1_pos, PosType pair2_gpos)
{
    #define MAX_GPOS_DELTA 1000 // paired reads are usually with a delta less than 300 - so this is more than enough

    ContextP ctx = CTX(FASTQ_GPOS);

    // case: we are pair-2 ; pair-1 is aligned ; delta is small enough : store a delta
    PosType gpos_delta = pair2_gpos - pair1_pos; 
    if (pair1_pos != NO_GPOS && pair2_gpos != NO_GPOS &&
        gpos_delta <= MAX_GPOS_DELTA && gpos_delta >= -MAX_GPOS_DELTA) {
    
        SNIPi2(SNIP_SPECIAL, FASTQ_SPECIAL_PAIR2_GPOS, gpos_delta);
        seg_by_ctx (VB, STRa(snip), ctx, 0);
    }

    // case: otherwise, store verbatim
    else {
        BNXT32 (ctx->local) = (uint32_t)pair2_gpos;
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, FASTQ_SPECIAL_PAIR2_GPOS }, 2, ctx, 0); // lookup from local and advance pair.next to consume gpos
    }
}

// note: in files up to v13, this is triggereed by v13_SNIP_FASTQ_PAIR2_GPOS
SPECIAL_RECONSTRUCTOR (fastq_special_PAIR2_GPOS)
{
    // case: no delta
    if (!snip_len) {
        if (CTX(FASTQ_SQBITMAP)->pair1_is_aligned)  // always false for V<=13
            ctx->pair.next++; // we didn't use this pair value

        new_value->i = NEXTLOCAL(uint32_t, ctx);
    }

    // case: pair-1 is aligned, and GPOS is a delta vs pair-1
    else {
        int64_t pair_value = (int64_t)(z_file->genozip_version >= 14 ? *B32 (ctx->pair, ctx->pair.next++) // starting v14, only aligned lines have GPOS
                                                                     : *B32 (ctx->pair, vb->line_i));     // up to v13, all lines segged GPOS (possibly NO_GPOS value, but not in this case, since we have a delta)
        int64_t delta = (int64_t)strtoull (snip, NULL, 10 /* base 10 */); 
        new_value->i = pair_value + delta; // just sets value, doesn't reconstruct

        ASSPIZ (pair_value != NO_GPOS, "pair_value=NO_GPOS - not expected as we have delta=%d", (int)delta);
    }
 
    return HAS_NEW_VALUE;
}

// can only be called before fastq_special_PAIR2_GPOS, because it inquires GPOS.pair.next
bool fastq_piz_get_pair2_is_forward (VBlockP vb)
{
    ContextP ctx = CTX(FASTQ_STRAND);

    if (z_file->genozip_version <= 13) {
        bool is_forward_pair_1 = bit_array_get ((BitArrayP)&ctx->pair, vb->line_i); // up to v13, all lines had strand, which was 0 if unaligned
        return NEXTLOCALBIT (ctx) ? is_forward_pair_1 : !is_forward_pair_1;
    }

    else if (CTX(FASTQ_SQBITMAP)->pair1_is_aligned) {
        bool is_forward_pair_1 = bit_array_get ((BitArrayP)&ctx->pair, CTX(FASTQ_GPOS)->pair.next); // gpos_ctx->pair.next is an iterator for both gpos and strand, and is incremented in fastq_special_PAIR2_GPOS
        return NEXTLOCALBIT (ctx) ? is_forward_pair_1 : !is_forward_pair_1;
    }

    else
        return NEXTLOCALBIT (ctx);
}
