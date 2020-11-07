// ------------------------------------------------------------------
//   fast.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fast_private.h"
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

// returns true if txt_data[txt_i] is the end of a FASTQ line (= block of 4 lines in the file)
static inline bool fastq_is_end_of_line (VBlock *vb, int32_t txt_i) // index of a \n in txt_data
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
    for (int32_t i=txt_i-1; i >= 0; i--) 
        if (txt[i] == '\n') 
            return IS_NL_BEFORE_QUAL_LINE(i);

    // we can't find a complete FASTQ block in the entire vb data
    ABORT ("Error when reading %s: last FASTQ record appears truncated, or the record is bigger than vblock", txt_name);
    return 0; // squash compiler warning ; never reaches here
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
uint32_t fastq_unconsumed (VBlockP vb)
{    
    for (int32_t i=vb->txt_data.len-1; i >= 0; i--) {
        // in FASTQ - an "end of line" is one that the next character is @, or it is the end of the file
        if (vb->txt_data.data[i] == '\n' && fastq_is_end_of_line (vb, i)) 
            return vb->txt_data.len-1 - i;
    }

    ABORT ("Error in fastq_unconsumed: no new line was found in vb=%u", vb->vblock_i);
    return 0; // quiet compiler warning
}

// called by txtfile_read_vblock when reading the 2nd file in a fastq pair - counts the number of fastq "lines" (each being 4 textual lines),
// comparing to the number of lines in the first file of the pair
// returns true if we have at least as much as needed, and sets unconsumed_len to the amount of excess characters read
// returns false is we don't yet have pair_1_num_lines lines - we need to read more
bool fastq_txtfile_have_enough_lines (VBlockP vb_, uint32_t *unconsumed_len)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;

    const char *next  = FIRSTENT (const char, vb->txt_data);
    const char *after = AFTERENT (const char, vb->txt_data);

    for (uint32_t line_i=0; line_i < vb->pair_num_lines * 4; line_i++) {
        while (*next != '\n' && next < after) next++; // note: next[-1] in the first iteration an underflow byte of the buffer, so no issue
        if (next >= after) 
            return false;
        next++; // skip newline
    }

    *unconsumed_len = after - next;
    return true;
}

// used in case of flag.optimize_DESC to count the number of lines, as we need it for the description
static void fastq_txtfile_count_lines (VBlockP vb)
{
    uint32_t num_lines = 0;
    ARRAY (const char, txt, vb->txt_data);

    for (uint32_t i=0; i < vb->txt_data.len; i++)
        if (txt[i] == '\n') num_lines++;

    ASSERT (num_lines % 4 == 0, "Error in fastq_txtfile_count_lines: expecting number of txt lines in VB to be a multiple of 4, but found %u", num_lines);

    vb->first_line = txt_file->num_lines + 1; // this is normally not used in ZIP
    txt_file->num_lines += num_lines / 4;     // update here instead of in zip_update_txt_counters;
}

void fastq_zip_read_one_vb (VBlockP vb)
{
    // in case we're optimizing DESC in FASTQ, we need to know the number of lines
    if (flag.optimize_DESC)
        fastq_txtfile_count_lines (vb);
}

// case of --optimize-DESC: generate the prefix of the read name from the txt file name
// eg. "../../fqs/sample.1.fq.gz" -> "@sample-1:"
static void fastq_get_optimized_desc_read_name (VBlockFAST *vb)
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

// called by I/O thread at the beginning of zipping this file
void fastq_zip_initialize (void)
{
    // reset lcodec for STRAND and GPOS, as these may change between PAIR_1 and PAIR_2 files
    z_file->contexts[FASTQ_STRAND].lcodec = CODEC_UNKNOWN;
    z_file->contexts[FASTQ_GPOS  ].lcodec = CODEC_UNKNOWN;
}

// called by Compute thread at the beginning of this VB
void fastq_seg_initialize (VBlockFAST *vb)
{
    if (flag.reference == REF_EXTERNAL || flag.reference == REF_EXT_STORE) {
        vb->contexts[FASTQ_STRAND].ltype = LT_BITMAP;
        vb->contexts[FASTQ_GPOS  ].ltype = LT_UINT32;
        vb->contexts[FASTQ_GPOS  ].flags = CTX_FL_STORE_INT;
    }

    vb->contexts[FASTQ_SQBITMAP].ltype = LT_BITMAP; 
    vb->contexts[FASTQ_SQBITMAP].inst  = CTX_INST_LOCAL_ALWAYS;

    codec_acgt_comp_init ((VBlockP)vb);

     if (flag.pair == PAIR_READ_2) {
        vb->contexts[FASTQ_GPOS]  .inst   = CTX_INST_PAIR_LOCAL;
        vb->contexts[FASTQ_STRAND].inst   = CTX_INST_PAIR_LOCAL; 

        piz_uncompress_all_ctxs ((VBlockP)vb, vb->pair_vb_i);

        vb->z_data.len = 0; // we've finished reading the pair file z_data, next, we're going to write to z_data our compressed output
    }

    if (flag.optimize_DESC) 
        fastq_get_optimized_desc_read_name (vb);

    // in --stats, consolidate stats into SQBITMAP
    stats_set_consolidation ((VBlockP)vb, FASTQ_SQBITMAP, 4, FASTQ_NONREF, FASTQ_NONREF_X, FASTQ_GPOS, FASTQ_STRAND);
}

void fastq_seg_finalize (VBlockP vb)
{
    // for qual data - select domqual compression if possible, or fallback 
    if (!codec_domq_comp_init (vb, FASTQ_QUAL, fastq_zip_qual)) {
        vb->contexts[FASTQ_QUAL].ltype  = LT_SEQUENCE; // might be overridden by codec_domq_compress
        vb->contexts[FASTQ_QUAL].inst   = 0; // don't inherit from previous file 
    }

    // top level snip
    Container top_level = { 
        .repeats   = vb->lines.len,
        .flags     = CONTAINER_TOPLEVEL | CONTAINER_FILTER_ITEMS | CONTAINER_FILTER_REPEATS,
        .num_items = 7,
        .items     = { { (DictId)dict_id_fields[FASTQ_DESC],     DID_I_NONE, ""  },
                       { (DictId)dict_id_fields[FASTQ_E1L],      DID_I_NONE, ""  }, // note: we have 2 EOL contexts, so we can show the correct EOL if in case of --header-only
                       { (DictId)dict_id_fields[FASTQ_SQBITMAP], DID_I_NONE, ""  },
                       { (DictId)dict_id_fields[FASTQ_E2L],      DID_I_NONE, "+" }, // + is the "separator" after the 2nd end-of-line
                       { (DictId)dict_id_fields[FASTQ_E2L],      DID_I_NONE, ""  },
                       { (DictId)dict_id_fields[FASTQ_QUAL],     DID_I_NONE, ""  },
                       { (DictId)dict_id_fields[FASTQ_E2L],      DID_I_NONE, ""  } }
    };

    container_seg_by_ctx (vb, &vb->contexts[FASTQ_TOPLEVEL], &top_level, 0, 0, vb->lines.len); // account for '+' - one for each line
}


// called from I/O thread ahead of zip or piz a pair 2 vb - to read data we need from the previous pair 1 file
// returns true if successful, false if there isn't a vb with vb_i in the previous file
bool fastq_read_pair_1_data (VBlockP vb_, uint32_t first_vb_i_of_pair_1, uint32_t last_vb_i_of_pair_1)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;
    uint64_t save_offset = file_tell (z_file);
    uint64_t save_disk_so_far = z_file->disk_so_far;

    vb->pair_vb_i = first_vb_i_of_pair_1 + (vb->vblock_i - last_vb_i_of_pair_1 - 1);

    const SectionListEntry *sl = sections_vb_first (vb->pair_vb_i, true);
    if (!sl) return false;

    // get num_lines from vb header
    SectionHeaderVbHeader *vb_header = zfile_read_section_header (vb_, sl->offset, vb->pair_vb_i, SEC_VB_HEADER);
    vb->pair_num_lines = BGEN32 (vb_header->num_lines);

    buf_free (&vb_->compressed); // allocated by zfile_read_section_header

    // read into ctx->pair the data we need from our pair: DESC and its components, GPOS and STRAND
    sl++;
    buf_alloc (vb, &vb->z_section_headers, MAX ((MAX_DICTS * 2 + 50),  vb->z_section_headers.len + MAX_SUBFIELDS + 10) * sizeof(uint32_t), 0, "z_section_headers", 1); // room for section headers  

    while (sl->section_type == SEC_B250 || sl->section_type == SEC_LOCAL) {
        
        if (((dict_id_is_type_1 (sl->dict_id) || sl->dict_id.num == dict_id_fields[FASTQ_DESC]) && sl->section_type == SEC_B250) ||
            ((sl->dict_id.num == dict_id_fields[FASTQ_GPOS] || sl->dict_id.num == dict_id_fields[FASTQ_STRAND]) && sl->section_type == SEC_LOCAL)) { // these are local sections
            
            NEXTENT (uint32_t, vb->z_section_headers) = vb->z_data.len; 
            int32_t ret = zfile_read_section (z_file, vb, vb->pair_vb_i, &vb->z_data, "data", sl->section_type, sl); // returns 0 if section is skipped
            ASSERT (ret != EOF, "Error in fastq_read_pair_1_data: vb_i=%u failed to read from pair_vb=%u dict_id=%s", vb->vblock_i, vb->pair_vb_i, err_dict_id (sl->dict_id));
        }
        
        sl++;
    }

    file_seek (z_file, save_offset, SEEK_SET, false); // restore
    z_file->disk_so_far = save_disk_so_far;
    
    return true;
}

// I/O thread: called from piz_read_one_vb as DTPZ(piz_read_one_vb)
bool fastq_piz_read_one_vb (VBlockP vb, ConstSectionListEntryP sl)
{
    // if we're grepping we we uncompress and reconstruct the DESC from the I/O thread, and terminate here if this VB is to be skipped
    if (flag.grep && !fast_piz_test_grep ((VBlockFAST *)vb)) return false; 

    // in case of this is a paired fastq file, and this is the 2nd file of a pair - 
    // we also read the equivalent sections from the first (bound) file
    ARRAY (const unsigned, section_index, vb->z_section_headers);
    for (uint32_t sec_i=1; sec_i < vb->z_section_headers.len; sec_i++) {
        SectionHeaderCtx *header = (SectionHeaderCtx *)ENT (char, vb->z_data, section_index[sec_i]);
        if (header->h.flags & CTX_FL_PAIRED) {
            uint32_t prev_file_first_vb_i, prev_file_last_vb_i;
            sections_get_prev_file_vb_i (sl, &prev_file_first_vb_i, &prev_file_last_vb_i);

            fastq_read_pair_1_data (vb, prev_file_first_vb_i, prev_file_last_vb_i);
            break;
        }
    }

    return true;
}

uint32_t fastq_get_pair_vb_i (VBlockP vb)
{
    return ((VBlockFAST *)vb)->pair_vb_i;
}

// concept: we treat every 4 lines as a "line". the Description/ID is stored in DESC dictionary and segmented to subfields D?ESC.
// The sequence is stored in SEQ data. In addition, we utilize the TEMPLATE dictionary for metadata on the line, namely
// the length of the sequence and whether each line has a \r.
const char *fastq_seg_txt_line (VBlockFAST *vb, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    ZipDataLineFAST *dl = DATA_LINE (vb->line_i);

    const char *next_field, *field_start=field_start_line;
    unsigned field_len=0;
    char separator;

    int32_t len = (int32_t)(AFTERENT (char, vb->txt_data) - field_start_line);

    // the leading @ - just verify it (it will be included in D0ESC subfield)
    ASSSEG (*field_start != '\n', field_start, "%s: Invalid FASTQ file format: unexpected newline", global_cmd);

    ASSSEG (*field_start == '@', field_start, "%s: Invalid FASTQ file format: expecting description line to start with @ but it starts with %c",
            global_cmd, *field_start);

    // DESC - the description/id line is vendor-specific. example:
    // @A00910:85:HYGWJDSXX:1:1101:3025:1000 1:N:0:CAACGAGAGC+GAATTGAGTG (<-- this is Illumina format)
    // See here for details of Illumina subfields: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    next_field = seg_get_next_line (vb, field_start, &len, &field_len, has_13, "DESC");
 
    // if flag.optimize_DESC is on, we replace the description with filename:line_i 
    unsigned unoptimized_len = 0; // 0 unless optimized
    if (flag.optimize_DESC) {
        unoptimized_len = field_len;
        field_start = vb->optimized_desc;
        field_len = vb->optimized_desc_len + str_int (vb->first_line + vb->line_i, &vb->optimized_desc[vb->optimized_desc_len]);   

        vb->vb_data_size -= unoptimized_len - field_len;
    }

    // we segment it using / | : and " " as separators. 
    seg_compound_field ((VBlockP)vb, &vb->contexts[FASTQ_DESC], field_start, field_len, true, unoptimized_len, 0);
    SEG_EOL (FASTQ_E1L, true);

    // SEQ - just get the whole line
    const char *seq_start = next_field;
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (vb, next_field, &len, true, false, false, &dl->seq_len, &separator, has_13, "SEQ");

    // case: compressing without a reference - all data goes to "nonref", and we have no bitmap
    if (flag.ref_use_aligner) 
        aligner_seg_seq ((VBlockP)vb, &vb->contexts[FASTQ_SQBITMAP], seq_start, dl->seq_len);

    else {
        Context *nonref_ctx = &vb->contexts[FASTQ_NONREF];
        buf_alloc ((VBlockP)vb, &nonref_ctx->local, MAX (nonref_ctx->local.len + dl->seq_len + 3, vb->lines.len * (dl->seq_len + 5)), CTX_GROWTH, "context->local", FASTQ_NONREF); 
        buf_add (&nonref_ctx->local, seq_start, dl->seq_len);
    }

    // Add LOOKUP snip with seq_len
    char snip[10];
    snip[0] = SNIP_LOOKUP;
    unsigned seq_len_str_len = str_int (dl->seq_len, &snip[1]);
    seg_by_ctx ((VBlockP)vb, snip, 1 + seq_len_str_len, &vb->contexts[FASTQ_SQBITMAP], 0, 0); 
    vb->contexts[FASTQ_NONREF].txt_len += dl->seq_len; // account for the txt data in NONREF

    SEG_EOL (FASTQ_E2L, true);

    // PLUS - next line is expected to be a "+" (note: we don't seg the +, it is recorded a separator in the top level Container)
    GET_LAST_ITEM ("+");
    ASSSEG (*field_start=='+' && field_len==1, field_start, "%s: Invalid FASTQ file format: expecting middle line to be a \"+\" (with no spaces) but it is \"%.*s\"",
            global_cmd, field_len, field_start);

    SEG_EOL (FASTQ_E2L, true); // account for ascii-10

    // QUAL - just get the whole line and make sure its length is the same as SEQ
    dl->qual_data_start = next_field - vb->txt_data.data;
    GET_LAST_ITEM ("QUAL");
    vb->contexts[FASTQ_QUAL].local.len += dl->seq_len;
    vb->contexts[FASTQ_QUAL].txt_len   += dl->seq_len;

    // End Of Line    
    SEG_EOL (FASTQ_E2L, true);

    ASSSEG (field_len == dl->seq_len, field_start, "%s: Invalid FASTQ file format: sequence_len=%u and quality_len=%u. Expecting them to be the same.\nSEQ=%.*s\nQUAL==%.*s",
            global_cmd, dl->seq_len, field_len, dl->seq_len, seq_start, field_len, field_start);
 
    return next_field;
}

// callback function for compress to get data of one line (called by codec_bz2_compress)
void fastq_zip_qual (VBlock *vb, uint32_t vb_line_i, 
                                          char **line_qual_data, uint32_t *line_qual_len, // out
                                          uint32_t maximum_len) 
{
    ZipDataLineFAST *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_qual_len  = MIN (dl->seq_len, maximum_len);
    
    if (!line_qual_data) return; // only lengths were requested

    *line_qual_data = ENT (char, vb->txt_data, dl->qual_data_start);

    // note - we optimize just before compression - likely the string will remain in CPU cache
    // removing the need for a separate load from RAM
    if (flag.optimize_QUAL) optimize_phred_quality_string (*line_qual_data, *line_qual_len);
}

// returns true if section is to be skipped reading / uncompressing
bool fastq_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT sections

    // note that piz_read_global_area rewrites --header-only as flag.header_one: skip all items but DESC and E1L
    if (flag.header_one && 
        (dict_id.num == dict_id_fields[FASTQ_E2L]      || dict_id.num == dict_id_fields[FASTQ_SQBITMAP] || 
         dict_id.num == dict_id_fields[FASTQ_NONREF]   || dict_id.num == dict_id_fields[FASTQ_NONREF_X] || 
         dict_id.num == dict_id_fields[FASTQ_GPOS]     || dict_id.num == dict_id_fields[FASTQ_STRAND]   || 
         dict_id.num == dict_id_fields[FASTQ_QUAL]     || dict_id.num == dict_id_fields[FASTQ_DOMQRUNS] ))
        return true;
        
    // when grepping by I/O thread - skipping all sections but DESC
    if (flag.grep && (vb->grep_stages == GS_TEST) && 
        dict_id.num != dict_id_fields[FASTQ_DESC] && !dict_id_is_fast_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if (flag.grep && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == dict_id_fields[FASTQ_DESC] || dict_id_is_fast_desc_sf (dict_id)))
        return true;

    return false;
}

// filtering during reconstruction: called by container_reconstruct_do for each fastq record (repeat) and each toplevel item
CONTAINER_FILTER_FUNC (fastq_piz_filter)
{
    if (dict_id.num == dict_id_fields[FASTQ_TOPLEVEL]) {
        if (item < 0)   // filter for repeat (FASTQ record)
            vb->line_i = 4 * (vb->first_line + rep); // each vb line is a fastq record which is 4 txt lines (needed for --pair)

        else { // filter for item

            // case: --grep (note: appears before --header-one filter below, so both can be used together)
            if (flag.grep && item == 2 /* first EOL */) {
                *AFTERENT (char, vb->txt_data) = 0; // for strstr
                if (!strstr (ENT (char, vb->txt_data, vb->line_start), flag.grep))
                    vb->dont_show_curr_line = true; // container_reconstruct_do will rollback the line
            }

            // case: --header-one or --header-only: dont show items 2+. note that piz_read_global_area rewrites --header-only as flag.header_one
            if (flag.header_one && item >= 2) return false; // skip this item
        }
    }

    return true; // show this item as normal
}

// PIZ: SEQ reconstruction 
void fastq_piz_reconstruct_seq (VBlock *vb_, Context *bitmap_ctx, const char *seq_len_str, unsigned seq_len_str_len)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;
 
    int64_t seq_len_64;
    ASSERT (str_get_int (seq_len_str, seq_len_str_len, &seq_len_64), "Error in fastq_piz_reconstruct_seq: could not parse integer \"%.*s\"", seq_len_str_len, seq_len_str);
    vb->seq_len = (uint32_t)seq_len_64;

    aligner_reconstruct_seq (vb_, bitmap_ctx, vb->seq_len, (vb->pair_vb_i > 0));
}
