// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

#define dict_id_is_fastq_desc_sf dict_id_is_type_1
#define dict_id_fastq_desc_sf dict_id_type_1

typedef struct {
    uint32_t seq_data_start, qual_data_start; // start within vb->txt_data
    uint32_t seq_len;                         // length of SEQ and QUAL within vb->txt_data (they are identical per FASTQ spec)
} ZipDataLineFASTQ;

// IMPORTANT: if changing fields in VBlockFASTQ, also update vb_fast_release_vb 
typedef struct VBlockFASTQ {
    VBLOCK_COMMON_FIELDS

    // pairing stuff - used if we are the 2nd file in the pair 
    uint32_t pair_vb_i;      // the equivalent vb_i in the first file, or 0 if this is the first file
    uint32_t pair_num_lines; // number of lines in the equivalent vb in the first file
    char *optimized_desc;    // base of desc in flag.optimize_DESC 
    uint32_t optimized_desc_len;

} VBlockFASTQ;

#define DATA_LINE(i) ENT (ZipDataLineFASTQ, vb->lines, i)

unsigned fastq_vb_size (DataType dt) { return sizeof (VBlockFASTQ); }
unsigned fastq_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTQ); }

void fastq_vb_release_vb (VBlockFASTQ *vb)
{
    vb->pair_num_lines = vb->pair_vb_i = vb->optimized_desc_len = 0;
    FREE (vb->optimized_desc);
}

void fastq_vb_destroy_vb (VBlockFASTQ *vb)
{
}

//-----------------------
// TXTFILE stuff
//-----------------------

// returns true if txt_data[txt_i] is the end of a FASTQ line (= block of 4 lines in the file); -1 if out of data
static inline int fastq_is_end_of_line (VBlock *vb, uint32_t first_i, int32_t txt_i) // index of a \n in txt_data
{
#   define IS_NL_BEFORE_QUAL_LINE(i) \
        ((i > 3) && ((txt[i-2] == '\n' && txt[i-1] == '+') || /* \n line ending case */ \
                     (txt[i-3] == '\n' && txt[i-2] == '+' && txt[i-1] == '\r'))) /* \r\n line ending case */
    
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    // if it is a @ then that @ is a new record EXCEPT if the previous character is + and then
    // @ is actually a quality value... (we check two previous characters, as the previous one might be \r)
    if (txt_i < vb->txt_data.len-1)
        return txt[txt_i+1] == '@' && !IS_NL_BEFORE_QUAL_LINE(txt_i);

    // if we're at the end of the line, we scan back to the previous \n and check if it is at the +
    for (int32_t i=txt_i-1; i >= first_i; i--) 
        if (txt[i] == '\n') 
            return IS_NL_BEFORE_QUAL_LINE(i);

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
                case false : continue;                      // not end of line, continue searching
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
        uint32_t num_lines = str_count_char (vb->txt_data.data, vb->txt_data.len, '\n');
        ASSERT (num_lines % 4 == 0, "expecting number of txt lines in VB to be a multiple of 4, but found %u", num_lines);

        vb->first_line = txt_file->num_lines + 1; // this is normally not used in ZIP
        txt_file->num_lines += num_lines / 4;     // update here instead of in zip_update_txt_counters;
    }
}

// case of --optimize-DESC: generate the prefix of the read name from the txt file name
// eg. "../../fqs/sample.1.fq.gz" -> "@sample-1:"
static void fastq_get_optimized_desc_read_name (VBlockFASTQ *vb)
{
    vb->optimized_desc = MALLOC (strlen (txt_file->basename) + 30); // leave room for the line number to follow
    vb->optimized_desc[0] = '@';
    strcpy (&vb->optimized_desc[1], txt_file->basename);
    file_get_raw_name_and_type (&vb->optimized_desc[1], NULL, NULL); // remove file type extension
    vb->optimized_desc_len = strlen (vb->optimized_desc) + 1; // +1 :
    vb->optimized_desc[vb->optimized_desc_len-1] = ':';

    // replace '.' in the filename with '-' as '.' is a separator in seg_compound_field and would needless inflate the number of contexts
    for (unsigned i=0; i < vb->optimized_desc_len; i++)
        if (vb->optimized_desc[i] == '.') vb->optimized_desc[i] = '-';
}

// called by main thread at the beginning of zipping this file
void fastq_zip_initialize (void)
{
    // reset lcodec for STRAND and GPOS, as these may change between PAIR_1 and PAIR_2 files
    ZCTX(FASTQ_STRAND)->lcodec = CODEC_UNKNOWN;
    ZCTX(FASTQ_GPOS  )->lcodec = CODEC_UNKNOWN;
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

    if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) {
        strand_ctx->ltype = LT_BITMAP;
        gpos_ctx->ltype   = LT_UINT32;
        gpos_ctx->flags.store = STORE_INT;
    }

    sqbitmap_ctx->ltype = LT_BITMAP; 
    sqbitmap_ctx->local_always = true;

    codec_acgt_comp_init ((VBlockP)vb);

     if (flag.pair == PAIR_READ_2) {

        ASSERT (vb->lines.len == vb->pair_num_lines, "in vb=%u (PAIR_READ_2): pair_num_lines=%u but lines.len=%u",
                vb->vblock_i, vb->pair_num_lines, (unsigned)vb->lines.len);

        gpos_ctx->pair_local = strand_ctx->pair_local = true;

        piz_uncompress_all_ctxs ((VBlockP)vb, vb->pair_vb_i);

        vb->z_data.len = 0; // we've finished reading the pair file z_data, next, we're going to write to z_data our compressed output
    }

    if (flag.optimize_DESC) 
        fastq_get_optimized_desc_read_name (vb);

    if (kraken_is_loaded) {
        CTX(FASTQ_TAXID)->flags.store    = STORE_INT;
        CTX(FASTQ_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTQ_TAXID)->counts_section = true; 
    }

    // in --stats, consolidate stats into SQBITMAP
    stats_set_consolidation ((VBlockP)vb, FASTQ_SQBITMAP, 4, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND);

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
        .nitems_lo      = 7,
        .items          = { 
            { .dict_id = (DictId)dict_id_fields[FASTQ_DESC],     },
            { .dict_id = (DictId)dict_id_fields[FASTQ_E1L],      }, // note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
            { .dict_id = (DictId)dict_id_fields[FASTQ_SQBITMAP], },
            { .dict_id = (DictId)dict_id_fields[FASTQ_E2L],      .seperator = "+" }, // + is the "separator" after the 2nd end-of-line
            { .dict_id = (DictId)dict_id_fields[FASTQ_E2L],      },
            { .dict_id = (DictId)dict_id_fields[FASTQ_QUAL],     },
            { .dict_id = (DictId)dict_id_fields[FASTQ_E2L],      } 
        }
    };

    container_seg_by_ctx (vb, CTX(FASTQ_TOPLEVEL), (ContainerP)&top_level, 0, 0, vb->lines.len); // account for '+' - one for each line
}

bool fastq_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == dict_id_fields[FASTQ_TOPLEVEL] ||
           dict_id.num == dict_id_fields[FASTQ_DESC]     ||
           dict_id.num == dict_id_fields[FASTQ_TAXID]    ||
           dict_id.num == dict_id_fields[FASTQ_E1L]      ||
           dict_id.num == dict_id_fields[FASTQ_E2L];
}

// ZIP/PIZ main thread: called ahead of zip or piz a pair 2 vb - to read data we need from the previous pair 1 file
// returns true if successful, false if there isn't a vb with vb_i in the previous file
bool fastq_read_pair_1_data (VBlockP vb_, uint32_t pair_vb_i, bool must_have)
{
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
    uint64_t save_offset = file_tell (z_file);
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
        
        if (dict_id_is_fastq_desc_sf (sl->dict_id) || sl->dict_id.num == dict_id_fields[FASTQ_DESC] ||
            sl->dict_id.num == dict_id_fields[FASTQ_GPOS] || sl->dict_id.num == dict_id_fields[FASTQ_STRAND]) { 
            
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

    if (flag.grep) {
        // in case of this is a paired fastq file, get just the pair_1 data that is needed to resolve the grep
        if (i_am_pair_2) {
            vb->grep_stages = GS_TEST; // tell piz_is_skip_section to skip decompressing sections not needed for determining the grep
            fastq_read_pair_1_data (vb_, pair_vb_i, true);
        }

        // if we're grepping we we uncompress and reconstruct the DESC from the main thread, and terminate here if this VB is to be skipped
        if (!piz_test_grep (vb_)) return false; // also updates vb->grep_stages
    }

    // in case of this is a paired fastq file, get all the pair_1 data not already fetched for the grep above
    if (i_am_pair_2) 
        fastq_read_pair_1_data (vb_, pair_vb_i, true);

    return true;
}

uint32_t fastq_get_pair_vb_i (VBlockP vb)
{
    return ((VBlockFASTQ *)vb)->pair_vb_i;
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *fastq_seg_txt_line (VBlockFASTQ *vb, const char *line_start, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb->line_i);

    const char *next_field, *field_start=line_start;
    unsigned field_len=0;
    char separator;

    int32_t len = (int32_t)(AFTERENT (char, vb->txt_data) - line_start);

    // the leading @ - just verify it (it will be included in D0ESC subfield)
    ASSSEG0 (*field_start != '\n', field_start, "Invalid FASTQ file format: unexpected newline");

    ASSSEG (*field_start == '@', field_start, "Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            *field_start);

    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG (<-- this is Illumina format)
    // See here for details of Illumina subfields: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    next_field = seg_get_next_line (vb, field_start, &len, &field_len, has_13, "DESC");
 
    if (kraken_is_loaded) {
        unsigned qname_len = strcspn (field_start + 1, " \t\r\n"); // +1 to skip the '@'

        unsigned taxid_found = kraken_seg_taxid ((VBlockP)vb, FASTQ_TAXID, field_start + 1, qname_len, false);

        // if not found tax id for this read, try again, perhaps removing /1 or /2
        if (!taxid_found) {
            if (qname_len > 2 && field_start[qname_len-2] == '/' &&
                (field_start[qname_len-1] == '1' ||field_start[qname_len-1] == '2'))
                qname_len -= 2;

            kraken_seg_taxid ((VBlockP)vb, FASTQ_TAXID, field_start + 1, qname_len, true); // this fail if missing
        }
    }

    // if flag.optimize_DESC is on, we replace the description with filename:line_i 
    unsigned unoptimized_len = 0; // 0 unless optimized
    if (flag.optimize_DESC) {
        unoptimized_len = field_len;
        field_start = vb->optimized_desc;
        field_len = vb->optimized_desc_len + str_int (vb->first_line + vb->line_i, &vb->optimized_desc[vb->optimized_desc_len]);   

        vb->recon_size -= unoptimized_len - field_len;
    }

    // we segment it using / | : and " " as separators. 
    static SegCompoundArg arg = { .slash = true, .pipe = true, .dot = true, .colon = true, .whitespace = true };
    seg_compound_field ((VBlockP)vb, CTX(FASTQ_DESC), field_start, field_len, arg, 0, 0);
    SEG_EOL (FASTQ_E1L, true);

    // SEQ - just get the whole line
    const char *seq_start = next_field;
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (vb, next_field, &len, GN_SEP, GN_IGNORE, GN_IGNORE, &dl->seq_len, &separator, has_13, "SEQ");

    // case: compressing without a reference - all data goes to "nonref", and we have no bitmap
    if (flag.ref_use_aligner) 
        aligner_seg_seq ((VBlockP)vb, CTX(FASTQ_SQBITMAP), seq_start, dl->seq_len);

    else {
        Context *nonref_ctx = CTX(FASTQ_NONREF);
        buf_alloc_old ((VBlockP)vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + dl->seq_len + 3, vb->lines.len * (dl->seq_len + 5)), CTX_GROWTH, "contexts->local"); 
        buf_add (&nonref_ctx->local, seq_start, dl->seq_len);
    }

    // Add LOOKUP snip with seq_len
    char snip[10];
    snip[0] = SNIP_LOOKUP;
    unsigned seq_len_str_len = str_int (dl->seq_len, &snip[1]);
    seg_by_ctx (vb, snip, 1 + seq_len_str_len, CTX(FASTQ_SQBITMAP), 0); 
    CTX(FASTQ_NONREF)->txt_len += dl->seq_len; // account for the txt data in NONREF

    SEG_EOL (FASTQ_E2L, true);

    // PLUS - next line is expected to be a "+" (note: we don't seg the +, it is recorded a separator in the top level Container)
    GET_LAST_ITEM (plus);
    ASSSEG (*field_start=='+' && field_len==1, field_start, "Invalid FASTQ file format: expecting middle line to be a \"+\" (with no spaces) but it is \"%.*s\"",
            field_len, field_start);

    SEG_EOL (FASTQ_E2L, true); // account for ascii-10

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = next_field - vb->txt_data.data;
    GET_LAST_ITEM (SAM_QUAL);
    CTX(FASTQ_QUAL)->local.len += dl->seq_len;
    CTX(FASTQ_QUAL)->txt_len   += dl->seq_len;

    // End Of Line    
    SEG_EOL (FASTQ_E2L, true);

    ASSSEG (str_is_in_range (field_start, field_len, 33, 126), field_start, "Invalid QUAL - it contains non-Phred characters: \"%.*s\"", 
            field_len, field_start);

    ASSSEG (field_len == dl->seq_len, field_start, "Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ = %.*s\nQUAL= %.*s",
            dl->seq_len, field_len, dl->seq_len, seq_start, field_len, field_start);
 
    return next_field;
}

// callback function for compress to get data of one line (called by codec_bz2_compress)
void fastq_zip_qual (VBlock *vb, uint64_t vb_line_i, 
                     char **line_qual_data, uint32_t *line_qual_len, // out
                     uint32_t maximum_len) 
{
    ZipDataLineFASTQ *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_qual_len  = MIN (dl->seq_len, maximum_len);
    
    if (!line_qual_data) return; // only lengths were requested

    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag.optimize_QUAL) optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

//-----------------
// PIZ stuff
//-----------------

// returns true if section is to be skipped reading / uncompressing
bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT/SEC_COUNTS sections

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast: skip all items but DESC and E1L
    if (flag.header_only_fast && 
        (dict_id.num == dict_id_fields[FASTQ_E2L]    || dict_id.num == dict_id_fields[FASTQ_SQBITMAP] || 
         dict_id.num == dict_id_fields[FASTQ_NONREF] || dict_id.num == dict_id_fields[FASTQ_NONREF_X] || 
         dict_id.num == dict_id_fields[FASTQ_GPOS]   || dict_id.num == dict_id_fields[FASTQ_STRAND]   || 
         dict_id.num == dict_id_fields[FASTQ_QUAL]   || dict_id.num == dict_id_fields[FASTQ_DOMQRUNS] ))
        return true;
        
    // when grepping by main thread - skipping all sections but DESC
    if (flag.grep && (vb->grep_stages == GS_TEST) && 
        dict_id.num != dict_id_fields[FASTQ_DESC] && !dict_id_is_fastq_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if (flag.grep && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == dict_id_fields[FASTQ_DESC] || dict_id_is_fastq_desc_sf (dict_id)))
        return true;

    // if we're doing --show-sex/coverage, we only need TOPLEVEL, FASTQ_SQBITMAP and GPOS
    if (flag.collect_coverage && 
        (dict_id.num == dict_id_fields[FASTQ_DESC]     || 
         dict_id.num == dict_id_fields[FASTQ_QUAL]     || 
         dict_id.num == dict_id_fields[FASTQ_DOMQRUNS] || 
         (dict_id.num == dict_id_fields[FASTQ_STRAND]   && !flag.bases)  ||
         (dict_id.num == dict_id_fields[FASTQ_NONREF]   && !flag.bases)   || 
         (dict_id.num == dict_id_fields[FASTQ_NONREF_X] && !flag.bases))) 
        return true;

    // no need for the TAXID data if user didn't specify --taxid
    if (flag.kraken_taxid==TAXID_NONE && dict_id.num == dict_id_fields[FASTQ_TAXID])
        return true;

    // if --count, we only need TOPLEVEL and the fields needed for the available filters (--taxid, --kraken, --grep, --bases)
    if (flag.count && sections_has_dict_id (st) &&
         (dict_id.num != dict_id_fields[FASTQ_TOPLEVEL] && 
         (dict_id.num != dict_id_fields[FASTQ_TAXID]    || flag.kraken_taxid == TAXID_NONE) && 
         (dict_id.num != dict_id_fields[FASTQ_DESC]     || (!kraken_is_loaded && !flag.grep)) && 
         (dict_id.num != dict_id_fields[FASTQ_SQBITMAP] || !flag.bases) && 
         (dict_id.num != dict_id_fields[FASTQ_NONREF]   || !flag.bases) && 
         (dict_id.num != dict_id_fields[FASTQ_NONREF_X] || !flag.bases) && 
         (dict_id.num != dict_id_fields[FASTQ_GPOS]     || !flag.bases) && 
         (dict_id.num != dict_id_fields[FASTQ_STRAND]   || !flag.bases) && 
         (!dict_id_is_fastq_desc_sf(dict_id)            || (!kraken_is_loaded && !flag.grep)))) 
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

    // dts_paired is not set, and this is v8 for v9. We proceed to inspect GPOS.local of the 2nd component to see if it is paired
    Section sl = NULL;
    sections_next_sec (&sl, SEC_TXT_HEADER); // first component txt header
    sections_next_sec (&sl, SEC_TXT_HEADER); // second component txt header
    sections_next_sec (&sl, SEC_VB_HEADER);  // first VB of second component txt header
    // edge cases not handled well here ^ 1. VB header not found 2. 2nd component has no VBs, and VB header actually belongs to the 3rd component
    
    // scan all B250 and Local looking for evidence of pairing
    while ((++sl)->st == SEC_B250 || sl->st == SEC_LOCAL) 
        if (sl->flags.ctx.paired) 
            return (z_file->z_flags.dts_paired = true); // assign and return
   
    return false; // no evidence of pairing
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat) and each toplevel item
CONTAINER_FILTER_FUNC (fastq_piz_filter)
{
    if (dict_id.num == dict_id_fields[FASTQ_TOPLEVEL]) {
        if (item < 0)   // filter for repeat (FASTQ record)
            vb->line_i = 4 * (vb->first_line + rep); // each vb line is a fastq record which is 4 txt lines (needed for --pair)

        else { // filter for item

            // case: --grep (note: appears before --header-only filter below, so both can be used together)
            if (flag.grep && item == 2 /* first EOL */) {
                *AFTERENT (char, vb->txt_data) = 0; // for strstr
                
                if (!vb->drop_curr_line && !strstr (ENT (char, vb->txt_data, vb->line_start), flag.grep))
                    vb->drop_curr_line = "grep";
            }

            // case: --header-only: dont show items 2+. note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
            if (flag.header_only_fast && item >= 2) 
                return false; // don't reconstruct this item (non-header textual)
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
        reconstruct_from_ctx ((VBlockP)vb, FASTQ_GPOS, 0, false);
        gpos = gpos_ctx->last_value.i;
    }

    WordIndex chrom_index;
    if (gpos != NO_GPOS && 
        (chrom_index = ref_contig_get_by_gpos (gref, gpos)) != WORD_INDEX_NONE) {

        if (flag.show_coverage || flag.show_sex)
            *ENT (uint64_t, vb->coverage, chrom_index) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*ENT (uint64_t, vb->read_count, chrom_index))++;
    }
    else {
        if (flag.show_coverage || flag.show_sex)
            *(AFTERENT (uint64_t, vb->coverage) - NUM_COVER_TYPES + CVR_UNMAPPED) += vb->seq_len;

        if (flag.show_coverage || flag.idxstats)
            (*(AFTERENT (uint64_t, vb->read_count) - NUM_COVER_TYPES + CVR_UNMAPPED))++;
    }
}

// PIZ: SEQ reconstruction 
void fastq_reconstruct_seq (VBlock *vb_, Context *bitmap_ctx, const char *seq_len_str, unsigned seq_len_str_len)
{
    VBlockFASTQ *vb = (VBlockFASTQ *)vb_;
 
    int64_t seq_len_64;
    ASSERT (str_get_int (seq_len_str, seq_len_str_len, &seq_len_64), "could not parse integer \"%.*s\"", seq_len_str_len, seq_len_str);
    vb->seq_len = (uint32_t)seq_len_64;

    // normal reconstruction
    if (!flag.collect_coverage) 
        aligner_reconstruct_seq (vb_, bitmap_ctx, vb->seq_len, (vb->pair_vb_i > 0));

    // just update coverage
    else
        fastq_update_coverage (vb);
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat)
CONTAINER_CALLBACK (fastq_piz_container_cb)
{
    // --taxid: filter out by Kraken taxid 
    if (flag.kraken_taxid && is_top_level) {
        
        if (!kraken_is_loaded && !kraken_is_included_stored (vb, FASTQ_TAXID, false))
            vb->drop_curr_line = "taxid";
        
        else if (kraken_is_loaded && flag.kraken_taxid >= 0) {
            const char *qname = last_txt (vb, FASTQ_DESC) + 1; // +1 to skip the "@"
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

