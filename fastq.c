// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
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
    bool qual_codec_enano;         // true if we can compress qual with CODEC_ENANO

} VBlockFASTQ;

#define VB_FASTQ ((VBlockFASTQ *)vb)

#define DATA_LINE(i) ENT (ZipDataLineFASTQ, vb->lines, i)

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

void fastq_vb_release_vb (VBlockFASTQ *vb)
{
    vb->pair_num_lines = vb->pair_vb_i = vb->optimized_desc_len = 0;
    vb->qual_codec_enano = 0;

    FREE (vb->optimized_desc);
}

void fastq_vb_destroy_vb (VBlockFASTQ *vb)
{
}

//-----------------------
// TXTFILE stuff
//-----------------------

// returns true if txt_data[txt_i] (which is a \n) is the end of a FASTQ record (= block of 4 lines in the file); -1 if out of data
static inline int fastq_is_end_of_line (VBlock *vb, uint32_t first_i, int32_t txt_i) // index of a \n in txt_data
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

    const char *next  = FIRSTENT (const char, vb->txt_data);
    const char *after = AFTERENT (const char, vb->txt_data);

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

void fastq_zip_read_one_vb (VBlockP vb)
{
    // in case we're optimizing DESC in FASTQ, we need to know the number of lines
    if (flag.optimize_DESC) {
        uint32_t num_lines = str_count_char (STRb(vb->txt_data), '\n');
        ASSERT (num_lines % 4 == 0, "expecting number of txt lines in VB to be a multiple of 4, but found %u", num_lines);

        vb->first_line = txt_file->num_lines + 1; // this is normally not used in ZIP
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
    // we just copy all reference contigs. this are not need for decompression, just for --coverage/--sex/--idxstats
    if (flag.reference == REF_EXTERNAL && z_file->num_txt_components_so_far == 1) // first file
        ctx_populate_zf_ctx_from_contigs (gref, FASTQ_CONTIG, ref_get_ctgs (gref)); 

    qname_zip_initialize ((DictId)_FASTQ_DESC);
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

    CTX(FASTQ_TOPLEVEL)->no_stons  = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(FASTQ_CONTIG)->flags.store = STORE_INDEX; // since v12

    Context *gpos_ctx     = CTX(FASTQ_GPOS);
    Context *strand_ctx   = CTX(FASTQ_STRAND);
    Context *sqbitmap_ctx = CTX(FASTQ_SQBITMAP);

    if (flag.reference & REF_ZIP_LOADED) {
        strand_ctx->ltype = LT_BITMAP;
        gpos_ctx->ltype   = LT_UINT32;
        gpos_ctx->flags.store = STORE_INT;
    }

    sqbitmap_ctx->ltype = LT_BITMAP; 
    sqbitmap_ctx->local_always = true;

    codec_acgt_comp_init (VB);

     if (flag.pair == PAIR_READ_2) {

        ASSERT (vb->lines.len == vb->pair_num_lines, "in vb=%u (PAIR_READ_2): pair_num_lines=%u but lines.len=%u",
                vb->vblock_i, vb->pair_num_lines, (unsigned)vb->lines.len);

        gpos_ctx->pair_local = strand_ctx->pair_local = true;

        piz_uncompress_all_ctxs (VB, vb->pair_vb_i);

        vb->z_data.len = 0; // we've finished reading the pair file z_data, next, we're going to write to z_data our compressed output
    }

    qname_seg_initialize (VB, FASTQ_DESC);

    if (flag.optimize_DESC) 
        fastq_get_optimized_desc_read_name (vb);

    if (kraken_is_loaded) {
        CTX(FASTQ_TAXID)->flags.store    = STORE_INT;
        CTX(FASTQ_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTQ_TAXID)->counts_section = true; 
    }

    // in --stats, consolidate stats into SQBITMAP
    stats_set_consolidation (VB, FASTQ_SQBITMAP, 4, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND);

    if (segconf.running) 
        segconf.qname_flavor = 0; // unknown

//xxx need to fix qname to identify "@ERR3278978.1 82a6ce63-eb2d-4812-ad19-136092a95f3d/1"    
    if (segconf.tech == TECH_ONP)
        vb->qual_codec_enano = true; // unless we find out during Seg that the data won't allow us

    COPY_TIMER (seg_initialize);
}

void fastq_seg_finalize (VBlockP vb)
{
    // for qual data - select domqual compression if possible, or fallback 
    if (!codec_domq_comp_init (vb, FASTQ_QUAL, fastq_zip_qual)) 
        CTX(FASTQ_QUAL)->ltype  = LT_SEQUENCE; // might be overridden by codec_domq_compress

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

    container_seg (vb, CTX(FASTQ_TOPLEVEL), (ContainerP)&top_level, prefixes, sizeof (prefixes), 2 * vb->lines.len); // account for the '@' and '+' - one of each for each line

//    if (VB_FASTQ->qual_codec_enano) 
//        codec_enano_seg_init (vb, CTX(FASTQ_QUAL));
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
    uint64_t save_offset = (uint64_t)file_tell (z_file, false);
    uint64_t save_disk_so_far = z_file->disk_so_far;

    vb->pair_vb_i = pair_vb_i;

    Section sl = sections_vb_first (vb->pair_vb_i, true);
    ASSERT (!must_have || sl, "file unexpectedly does not contain data for pair 1 vb_i=%u", pair_vb_i);
    if (!sl) return false;

    // get num_lines from vb header
    SectionHeaderVbHeader vb_header = zfile_read_section_header (vb_, sl->offset, vb->pair_vb_i, SEC_VB_HEADER).vb_header;
    vb->pair_num_lines = BGEN32 (vb_header.num_lines_prim);

    // read into ctx->pair the data we need from our pair: DESC and its components, GPOS and STRAND
    buf_alloc (vb, &vb->z_section_headers, MAX_FIELDS + 10, MAX_DICTS * 2 + 50, uint32_t, 0, "z_section_headers"); // room for section headers  

    for (sl++; sl->st == SEC_B250 || sl->st == SEC_LOCAL; sl++) {
        
        if (dict_id_is_fastq_desc_sf (sl->dict_id) || sl->dict_id.num == _FASTQ_DESC ||
            sl->dict_id.num == _FASTQ_GPOS || sl->dict_id.num == _FASTQ_STRAND) { 
            
            NEXTENT (uint32_t, vb->z_section_headers) = vb->z_data.len; 
            int32_t offset = zfile_read_section (z_file, vb, vb->pair_vb_i, &vb->z_data, "z_data", sl->st, sl); // returns 0 if section is skipped

            if (offset == SECTION_SKIPPED) vb->z_section_headers.len--; // section skipped
        }
    }

    file_seek (z_file, save_offset, SEEK_SET, false); // restore
    z_file->disk_so_far = save_disk_so_far;

    return true;
}

// main thread: called from piz_read_one_vb as DTPZ(piz_read_one_vb)
bool fastq_piz_read_one_vb (VBlockP vb_, Section sl)
{
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
    uint32_t pair_vb_i=0;
    bool i_am_pair_2 = vb && writer_get_pair (vb->vblock_i, &pair_vb_i) == 2;

    // in case of this is a paired fastq file, get all the pair_1 data not already fetched for the grep above
    if (i_am_pair_2) 
        fastq_read_pair_1_data (vb_, pair_vb_i, true);

    return true;
}

uint32_t fastq_get_pair_vb_i (VBlockP vb)
{
    return VB_FASTQ->pair_vb_i;
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *fastq_seg_txt_line (VBlockFASTQ *vb, const char *line_start, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb->line_i);

    ASSSEG0 (*line_start != '\n', line_start, "Invalid FASTQ file format: unexpected newline");

    ASSSEG (*line_start == '@', line_start, "Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            *line_start);

    const char *FASTQ_DESC_str = line_start + 1; // skip the '@' - its already included in the container prefix
    char separator;

    int32_t len = (int32_t)(AFTERENT (char, vb->txt_data) - FASTQ_DESC_str);
    
    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG <-- Illumina, See: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // @20A_V100002704L1C001R012000000/1 <-- BGI, see: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
    unsigned FASTQ_DESC_len;
    const char *FASTQ_SEQ_str = seg_get_next_line (vb, FASTQ_DESC_str, &len, &FASTQ_DESC_len, true, has_13, "DESC");
 
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

    // if flag.optimize_DESC is on, we replace the description with filename:line_i 
    unsigned optimized_len = 0; 
    if (flag.optimize_DESC) {
        optimized_len  = vb->optimized_desc_len + str_int (vb->first_line + vb->line_i, &vb->optimized_desc[vb->optimized_desc_len]);   
        vb->recon_size -= FASTQ_DESC_len - optimized_len;
    }

    if (segconf.running && vb->line_i==0)
        qname_segconf_discover_flavor (VB, FASTQ_DESC, 
                                       flag.optimize_DESC ? vb->optimized_desc : FASTQ_DESC_str, 
                                       flag.optimize_DESC ? optimized_len      : FASTQ_DESC_len);

    // attempt to parse the desciption line as qname, and fallback to tokenizer no qname format matches
    qname_seg (VB, CTX(FASTQ_DESC), 
               flag.optimize_DESC ? vb->optimized_desc : FASTQ_DESC_str, 
               flag.optimize_DESC ? optimized_len      : FASTQ_DESC_len, 
               0);
    SEG_EOL (FASTQ_E1L, true);

    // SEQ - just get the whole line
    dl->seq_data_start = ENTNUM (vb->txt_data, FASTQ_SEQ_str);
    const char *FASTQ_LINE3_str = seg_get_next_item (vb, FASTQ_SEQ_str, &len, GN_SEP, GN_IGNORE, GN_IGNORE, &dl->seq_len, &separator, has_13, "SEQ");

    // case: compressing without a reference - all data goes to "nonref", and we have no bitmap
    if (flag.aligner_available) 
        aligner_seg_seq (VB, CTX(FASTQ_SQBITMAP), FASTQ_SEQ_str, dl->seq_len);

    else {
        Context *nonref_ctx = CTX(FASTQ_NONREF);
        buf_alloc (VB, &nonref_ctx->local, 0, MAX_(nonref_ctx->local.len + dl->seq_len + 3, vb->lines.len * (dl->seq_len + 5)), char, CTX_GROWTH, "contexts->local"); 
        buf_add (&nonref_ctx->local, FASTQ_SEQ_str, dl->seq_len);
    }

    // Add LOOKUP snip with seq_len
    char snip[10];
    snip[0] = SNIP_LOOKUP;
    unsigned seq_len_str_len = str_int (dl->seq_len, &snip[1]);
    seg_by_ctx (VB, snip, 1 + seq_len_str_len, CTX(FASTQ_SQBITMAP), 0); 
    CTX(FASTQ_NONREF)->txt_len += dl->seq_len; // account for the txt data in NONREF

    SEG_EOL (FASTQ_E2L, true);

    // PLUS - next line is expected to be a "+", and be either empty or a copy of the DESC line (except for the '@' / '+' prefix)
    ASSSEG (*FASTQ_LINE3_str == '+', FASTQ_LINE3_str, "Invalid FASTQ file format: expecting middle line to be a \"+\", but it starts with a '%c'", *FASTQ_LINE3_str);
    FASTQ_LINE3_str++; len--; // skip the prefix

    unsigned FASTQ_LINE3_len;
    const char *FASTQ_QUAL_str = seg_get_next_line (vb, FASTQ_LINE3_str, &len, &FASTQ_LINE3_len, true, has_13, "LINE3");

    // line3 can be either empty, or a copy of DESC.
    if (!FASTQ_LINE3_len) 
        seg_by_did_i (VB, "", 0, FASTQ_LINE3, 0);

    else if (str_issame_ (FASTQ_LINE3_str, FASTQ_LINE3_len, FASTQ_DESC_str, FASTQ_DESC_len)) {

        // if --optimize-DESC, we always produce an empty line.
        if (flag.optimize_DESC) {
            seg_by_did_i (VB, "", 0, FASTQ_LINE3, 0);
            vb->recon_size -= FASTQ_DESC_len;
        }
        else 
            seg_by_did_i (VB, copy_desc_snip, copy_desc_snip_len, FASTQ_LINE3, (flag.optimize_DESC ? optimized_len : FASTQ_LINE3_len));
    }

    else 
        ASSSEG (false, FASTQ_LINE3_str, "Invalid FASTQ file format: expecting middle line to be a \"+\" with or without a copy of the description, but it is \"%.*s\"",
                FASTQ_LINE3_len+1, FASTQ_QUAL_str-1);

    SEG_EOL (FASTQ_E2L, true); // account for ascii-10

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = ENTNUM (vb->txt_data, FASTQ_QUAL_str);

    unsigned FASTQ_QUAL_len;
    const char *after = seg_get_next_item (vb, FASTQ_QUAL_str, &len, GN_SEP, GN_FORBIDEN, GN_IGNORE, &FASTQ_QUAL_len, &separator, has_13, "QUAL"); \

    CTX(FASTQ_QUAL)->local.len += dl->seq_len;
    CTX(FASTQ_QUAL)->txt_len   += dl->seq_len;

    // End Of Line    
    SEG_EOL (FASTQ_E2L, true);

    ASSSEG (str_is_in_range (FASTQ_QUAL_str, FASTQ_QUAL_len, 33, 126), FASTQ_QUAL_str, "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", 
            FASTQ_QUAL_len, FASTQ_QUAL_str);

    ASSSEG (FASTQ_QUAL_len == dl->seq_len, FASTQ_QUAL_str, "Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ = %.*s\nQUAL= %.*s",
            dl->seq_len, FASTQ_QUAL_len, dl->seq_len, FASTQ_SEQ_str, FASTQ_QUAL_len, FASTQ_QUAL_str);
 
     if (FASTQ_QUAL_len < ENANO_MIN_READ_LEN)
        vb->qual_codec_enano = false; // we cannot compress QUAL with CODEC_ENANO as prerequisites are not met

    return after;
}

// callback function for compress to get data of one line (called by codec_bz2_compress)
COMPRESSOR_CALLBACK (fastq_zip_qual) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len  = MIN_(dl->seq_len, maximum_size);
    
    if (!line_data) return; // only lengths were requested

    *line_data = ENT (char, vb->txt_data, dl->qual_data_start);

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag.optimize_QUAL) optimize_phred_quality_string (STRa(*line_data));

    if (is_rev) *is_rev = 0;
}

COMPRESSOR_CALLBACK (fastq_zip_seq) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);
    *line_data_len = dl->seq_len;
    *line_data     = ENT (char, vb->txt_data, dl->seq_data_start);
    if (is_rev) *is_rev = 0;
}

//-----------------
// PIZ stuff
//-----------------

// returns true if section is to be skipped reading / uncompressing
bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT/SEC_COUNTS sections

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
         (!dict_id_is_fastq_desc_sf(dict_id)            || !kraken_is_loaded))) 
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
        if (item < 0)   // filter for repeat (FASTQ record)
            vb->line_i = 4 * (vb->first_line + rep); // each vb line is a fastq record which is 4 txt lines (needed for --pair)

        else { // filter for item
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
    }

    return true; // reconstruct
}

static void fastq_update_coverage (VBlockFASTQ *vb)
{
    PosType gpos;
    Context *gpos_ctx = CTX(FASTQ_GPOS);

    if (!vb->pair_vb_i) { // first file of a pair ("pair 1") or a non-pair fastq 
        if (!gpos_ctx->local.len) return;
        gpos = NEXTLOCAL (uint32_t, gpos_ctx);
    }

    // 2nd file of a pair ("pair 2")
    else {
        // gpos: reconstruct, then cancel the reconstruction and just use last_value
        reconstruct_from_ctx (VB, FASTQ_GPOS, 0, false);
        gpos = gpos_ctx->last_value.i;
    }

    WordIndex ref_index;
    if (gpos != NO_GPOS && 
        (ref_index = ref_contig_get_by_gpos (gref, gpos, NULL)) != WORD_INDEX_NONE) {

        if (flag.show_coverage || flag.show_sex)
            *ENT (uint64_t, vb->coverage, ref_index) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*ENT (uint64_t, vb->read_count, ref_index))++;
    }
    else {
        if (flag.show_coverage || flag.show_sex)
            *(AFTERENT (uint64_t, vb->coverage) - NUM_COVER_TYPES + CVR_UNMAPPED) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*(AFTERENT (uint64_t, vb->read_count) - NUM_COVER_TYPES + CVR_UNMAPPED))++;
    }
}

// PIZ: SEQ reconstruction 
void fastq_reconstruct_seq (VBlock *vb_, Context *bitmap_ctx, STRp(seq_len_str))
{
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
 
    int64_t seq_len_64;
    ASSERT (str_get_int (STRa(seq_len_str), &seq_len_64), "could not parse integer \"%.*s\"", STRf(seq_len_str));
    vb->seq_len = (uint32_t)seq_len_64;

    // just update coverage
    if (flag.collect_coverage) 
        fastq_update_coverage (vb);

    // --qual-only: only set vb->seq_len without reconstructing
    else if (flag.qual_only) {}

    // normal reconstruction
    else 
        aligner_reconstruct_seq (vb_, bitmap_ctx, vb->seq_len, (vb->pair_vb_i > 0));
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat)
CONTAINER_CALLBACK (fastq_piz_container_cb)
{
    // --taxid: filter out by Kraken taxid 
    if (flag.kraken_taxid && is_top_level) {
        
        if (!kraken_is_loaded && !kraken_is_included_stored (vb, FASTQ_TAXID, false))
            vb->drop_curr_line = "taxid";
        
        else if (kraken_is_loaded && flag.kraken_taxid >= 0) {
            const char *qname = last_txt (vb, FASTQ_DESC) + !prefixes_len; // +1 to skip the "@", if '@' is in DESC and not in prefixes (for files up to version 12.0.13)
            unsigned qname_len = strcspn (qname, " \t\n\r");

            if (!kraken_is_included_loaded (vb, qname, qname_len)) 
                vb->drop_curr_line = "taxid";
        }
    }

    // --bases
    if (flag.bases && is_top_level && !vb->drop_curr_line &&
        !iupac_is_included_ascii (last_txt (vb, FASTQ_SQBITMAP), vb->last_txt_len (FASTQ_SQBITMAP)))
        vb->drop_curr_line = "bases";
}

