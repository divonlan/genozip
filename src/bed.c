// ------------------------------------------------------------------
//   bed.c
//   Copyright (C) 2023-2023 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "random_access.h"
#include "file.h"
#include "strings.h"
#include "piz.h"
#include "bed.h"
#include "dict_id.h"
#include "codec.h"
#include "dict_id_gen.h"
#include "chrom.h"
#include "stats.h"
#include "generic.h"

sSTRl(copy_TSTART_snip, 30);
sSTRl(copy_TEND_snip, 30);

// detect if a generic file is actually a BED, based on "browser" optional keyword
bool is_bed (STRp(header), bool *need_more)
{
    bool is_bed = (header_len >= 9 && !memcmp (header, "browser", 7));

    return is_bed;
}

// ZIP: called from txtfile_read_header
// returns header length if header read is complete, -1 not complete yet 
int32_t bed_is_header_done (bool is_eof)
{
    if (!evb->txt_data.len) return HEADER_NEED_MORE;

    ASSERT0 (!is_eof || *BLSTc(evb->txt_data) == '\n', "Missing newline at end of file");
    
    str_split (B1STc(evb->txt_data), evb->txt_data.len32, 0, '\n', line, true);

    n_lines--; // ignore empty last entry (if txt_data ends with \n) or partial last line (if it doesn't)

    uint32_t header_len = 0;
    int i; for (i=0; i < n_lines; i++) {
        if (!(line_lens[i] >= 7 && !memcmp (lines[i], "browser", 7)) &&
            !(line_lens[i] >= 5 && !memcmp (lines[i], "track", 5)) &&
            !(line_lens[i] >= 1 && lines[i][0] == '#') &&
            line_lens[i] != 0) break; // not a header line

        header_len += line_lens[i] + 1; // +1 for newline
    }

    return (i == n_lines && !is_eof) ? HEADER_NEED_MORE : header_len; // perhaps header is longer if all lines are header lines
}

void bed_zip_initialize (void)
{
    seg_prepare_snip_other (SNIP_COPY, _BED_TSTART, false, 0, copy_TSTART_snip);
    seg_prepare_snip_other (SNIP_COPY, _BED_TEND,   false, 0, copy_TEND_snip);
}

// called from seg_all_data_lines
void bed_seg_initialize (VBlockP vb)
{
    ctx_set_store (VB, STORE_INT, BED_START, BED_END, DID_EOL);

    ctx_set_no_stons (vb, 
                      BED_CHROM, // needs b250 node_index for random access
                      BED_NAME,  // required by seg_add_to_local_text
                      BED_START, BED_END, BED_TSTART, BED_TEND, DID_EOL); // as requied by seg_pos_field

    if (!segconf.is_sorted || segconf.running) 
        ctx_set_ltype (vb, LT_UINT32, BED_START, DID_EOL);

    ctx_set_ltype (vb, LT_DYN_INT, BED_SCORE, DID_EOL);
}

void bed_seg_finalize (VBlockP vb)
{
    if (segconf.running)
        segconf_finalize_is_sorted();

    // top level snip
    SmallContainer top_level = { 
        .repeats      = vb->lines.len32,
        .is_toplevel  = true,
        .nitems_lo    = 13,
        .items        = { { .dict_id = { _BED_CHROM   }, .separator = "\t" },
                          { .dict_id = { _BED_START   }, .separator = "\t" },
                          { .dict_id = { _BED_END     }, .separator = "\t" },
                          { .dict_id = { _BED_NAME    }, .separator = "\t" },
                          { .dict_id = { _BED_SCORE   }, .separator = "\t" },
                          { .dict_id = { _BED_STRAND  }, .separator = "\t" },
                          { .dict_id = { _BED_TSTART  }, .separator = "\t" },
                          { .dict_id = { _BED_TEND    }, .separator = "\t" },
                          { .dict_id = { _BED_RGB     }, .separator = "\t" },
                          { .dict_id = { _BED_BCOUNT  }, .separator = "\t" },
                          { .dict_id = { _BED_BSIZES  }, .separator = "\t" },
                          { .dict_id = { _BED_BSTARTS },                   },
                          { .dict_id = { _BED_EOL     },                   } }
    };

    // BED can have between 3 and 12 columns - remove redundant columns
    if (segconf.bed_num_columns < 12) {
        memmove (&top_level.items[segconf.bed_num_columns], &top_level.items[12], sizeof (ContainerItem));
        top_level.nitems_lo = segconf.bed_num_columns + 1;
        top_level.items[segconf.bed_num_columns-1].separator[0] = 0;    
    }

    container_seg (vb, CTX(BED_TOPLEVEL), (Container *)&top_level, 0, 0, 0);
}

bool bed_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _BED_STRAND ||
           dict_id.num == _BED_RGB;
}

static void bed_seg_START (VBlockP vb, STRp(start), WordIndex prev_line_chrom)
{
    PosType32 start_val;
    ContextP ctx = CTX(BED_START);

    if (segconf.is_sorted && !segconf.running)
        start_val = seg_pos_field (VB, BED_START, BED_END, 0, 0, STRa(start), 0, start_len+1);
    
    else {
        ASSSEG (str_get_int_range32 (STRa(start), 0, MAX_POS32, &start_val), start, "Invalid START=\"%.*s\"", STRf(start));

        seg_integer_fixed (vb, ctx, &start_val, false, start_len+1);
    }

    if (segconf.running) 
        segconf_test_sorted (vb, prev_line_chrom, start_val, ctx->last_value.i);

    ctx_set_last_value (vb, ctx, (int64_t)start_val);

    random_access_update_pos (vb, 0, BED_START); // after setting last_value
}

static void bed_seg_TSTART_TEND (VBlockP vb, Did did_i, Did base_did_i, STRp(value), STRp(copy_snip))
{
    if (str_issame_(STRa(value), STRtxtw(CTX(did_i)->last_txt)))
        seg_by_did (vb, STRa(copy_snip), did_i, value_len + 1);

    else 
        seg_pos_field (vb, did_i, base_did_i, 0, 0, STRa(value), 0, value_len + 1);

    set_last_txt(did_i, value);
}

SPECIAL_RECONSTRUCTOR (bed_piz_special_BCOUNT)
{
    new_value->i = container_peek_repeats (vb, CTX(BED_BSIZES), ',');
    
    if (reconstruct) 
        RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

bool bed_segconf_set_n_columns (STRp(line))
{
    str_split (line, line_len, 13, '\t', column, false); // 13 and not 12, to detect lines that have more than 12 columns

    static int valid_num_columns[] = { 3, 4, 5, 6, 8, 9, 12 };
    for (int i=0; i < ARRAY_LEN(valid_num_columns); i++)
        if (n_columns == valid_num_columns[i]) {
            segconf.bed_num_columns = n_columns;
            return true; // judging by the number of columns, so far this is a valid standard BED file
        }

    return false; // invalid number of columns for a standard BED file
}

rom bed_seg_txt_line (VBlockP vb, rom line, uint32_t remaining_txt_len, bool *has_13)     // index in vb->txt_data where this line starts
{
    typedef enum { BED_CHROM, START, END, NAME, SCORE, STRAND, TSTART, TEND, RGB, BCOUNT, BSIZES, BSTARTS } Columns __attribute__((unused));
    #define MAYBE_END_HERE(f) if (segconf.bed_num_columns == (f)+1) goto eol;

    int line_len = (rom)memchr (line, '\n', remaining_txt_len) - line;
    if (line[line_len-1] == '\r') {
        *has_13 = true;
        line_len--;
    }

    if (segconf.running && vb->line_i == 0 && !bed_segconf_set_n_columns (STRa(line)))
        return fallback_to_generic (vb);

    str_split (line, line_len, segconf.bed_num_columns, '\t', column, false);
    ASSSEG (n_columns == segconf.bed_num_columns, line, "Invalid BED sequence: expecting %u tab-separated columns, but found %u", 
            segconf.bed_num_columns, n_columns);

    // CHROM
    WordIndex prev_line_chrom = vb->chrom_node_index;
    chrom_seg (vb, STRi(column, BED_CHROM));

    // START
    bed_seg_START (vb, STRi(column, START), prev_line_chrom);

    // END
    seg_pos_field (vb, BED_END, BED_START, 0, 0, STRi(column, END), 0, column_lens[END] + 1);
    MAYBE_END_HERE(END);

    // NAME
    seg_add_to_local_text (vb, CTX(BED_NAME), STRi(column, NAME), LOOKUP_NONE, column_lens[NAME] + 1);
    MAYBE_END_HERE(NAME);

    // SCORE
    int64_t score;
    ASSSEG (str_get_int (STRi(column, SCORE), &score), columns[SCORE], "Invalid SCORE=\"%.*s\"", STRfi(column,SCORE));
    seg_add_to_local_resizable (vb, CTX(BED_SCORE), score, column_lens[SCORE]+1);
    MAYBE_END_HERE(SCORE);

    // STRAND
    seg_by_did (vb, STRi(column,STRAND), BED_STRAND, column_lens[STRAND] + 1);
    MAYBE_END_HERE(STRAND);

    // TSTART + TEND
    bed_seg_TSTART_TEND (vb, BED_TSTART, BED_START, STRi(column, TSTART), STRa(copy_TSTART_snip));
    bed_seg_TSTART_TEND (vb, BED_TEND, BED_END, STRi(column, TEND), STRa(copy_TEND_snip));
    MAYBE_END_HERE(TEND);

    // RGB
    seg_by_did (vb, STRi(column,RGB), BED_RGB, column_lens[RGB] + 1);
    MAYBE_END_HERE(RGB);

    // BCOUNT
    int64_t predicted_bcount = 1 + str_count_char (STRi(column,BSIZES)-1, ','); // without last char, to avoid counting a terminating ,
    int64_t bcount;

    if (!str_get_int(STRi(column,BCOUNT), &bcount) || bcount != predicted_bcount) {
        if (segconf.running) return fallback_to_generic (vb);
        ASSSEG (false, columns[BCOUNT], "Encountered BCOUNT=\"%.*s\" - but expecting it to be the number of items in BSIZES which is %u", STRfi(column,BCOUNT), (unsigned)predicted_bcount);
    }

    seg_by_did (vb, (char[]){ SNIP_SPECIAL, BED_SPECIAL_BCOUNT }, 2, BED_BCOUNT, column_lens[BED_BCOUNT] + 1);

    // BSIZES + BSTARTS
    seg_array (vb, CTX(BED_BSIZES),  BED_BSIZES,  STRi(column,BSIZES),  ',', 0, false, true, DICT_ID_NONE, column_lens[BSIZES]  + 1);
    seg_array (vb, CTX(BED_BSTARTS), BED_BSTARTS, STRi(column,BSTARTS), ',', 0, true,  true, DICT_ID_NONE, column_lens[BSTARTS] + 1);

eol:             
    SEG_EOL (BED_EOL, false);

    return &line[line_len] + 1;
}
