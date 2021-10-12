// ------------------------------------------------------------------
//   compound.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "compound.h"
#include "buffer.h"
#include "vblock.h"
#include "segconf.h"
#include "strings.h"
#include "seg.h"
#include "sam.h"
#include "stats.h"

// We break down the field (eg QNAME in SAM/BAM/Kraken or Description in FASTA/FASTQ) into subfields separated by / and/or : - and/or whitespace 
// these are vendor-defined strings.
// Up to MAX_ARRAY_ITEMS subfields are permitted - if there are more, then all the trailing part is just
// consider part of the last component.
// each subfield is stored in its own dictionary- the second character of the dict_id  the subfield number starting
// from 0 (0->9,a->z)
const char sep_with_space[256]    = { [':']=true, [';']=true, ['/']=true, ['|']=true, ['.']=true, ['_']=true, ['#']=true, [' ']=true, ['\t']=true, [1]=true };
const char sep_without_space[256] = { [':']=true, [';']=true, ['/']=true, ['|']=true, ['.']=true, ['_']=true, ['#']=true, };

static char illumina_7_snip[200], bgi_snip[200], copy_qname[50];
static uint32_t illumina_7_snip_len, bgi_snip_len, copy_qname_len;

// note: we need to re-initialize in each file, lest the data type has changed
void compound_zip_initialize (DictId qname_dict_id)
{
    // Read names look like: "A00488:61:HMLGNDSXX:4:1101:4345:1000". First five go into item[0]
    SmallContainer illumina_7_con = {
        .repeats             = 1,
        .nitems_lo           = 3,
        .items               = { { .dict_id = { _SAM_Q0NAME }, .separator = ":"  }, // note: _FASTQ_Q0NAME is the same dict_id "Q0NAME"
                                 { .dict_id = { _SAM_Q1NAME }, .separator = ":"  },
                                 { .dict_id = { _SAM_Q2NAME },                   } } };
    illumina_7_snip_len = sizeof illumina_7_snip;
    container_prepare_snip ((Container*)&illumina_7_con, 0, 0, illumina_7_snip, &illumina_7_snip_len);

    // Read names look like: Read names have a fixed-length format that looks like "E100020409L1C001R0030000801".
    SmallContainer bgi_con = {
        .repeats             = 1,
        .nitems_lo           = 5,
        .items               = { { .dict_id = { _SAM_Q0NAME } }, 
                                 { .dict_id = { _SAM_Q1NAME } },
                                 { .dict_id = { _SAM_Q2NAME } },
                                 { .dict_id = { _SAM_Q3NAME } },
                                 { .dict_id = { _SAM_Q4NAME } } } };

    bgi_snip_len = sizeof bgi_snip;
    container_prepare_snip ((Container*)&bgi_con, "\4\4E\4L\4C\4R\4", 10, bgi_snip, &bgi_snip_len);

    seg_prepare_snip_other (SNIP_COPY, qname_dict_id, false, 0, copy_qname);
}

void compound_seg_initialize (VBlockP vb, DidIType qname_did_i)
{
    if (segconf.running) return;

    if (segconf.is_illumina_7)
        stats_set_consolidation (vb, qname_did_i, 3, qname_did_i+1, qname_did_i+2, qname_did_i+3);

    if (segconf.is_bgi_E9L1C3R10)
        stats_set_consolidation (vb, qname_did_i, 5, qname_did_i+1, qname_did_i+2, qname_did_i+3, qname_did_i+4, qname_did_i+5);
}

//--------------------------------------------
// Default compound segger
//--------------------------------------------

typedef struct { 
    const char *item;       // pointer into field
    unsigned item_len;
    char sep;               // separator after item - 0 if none
    char sep2;              // double space separator - common in FASTA headers
    bool is_int;
    int64_t value;          // in case of an int
    unsigned leading_zeros; // in case of an int - these will go into the container prefix
} CompoundItem;

static void compound_field_split (const char *field, unsigned field_len, const char *is_sep,
                                  CompoundItem *items, uint8_t *n_items) // out: array of MAX_ARRAY_ITEMS
{
    bool in_digits = field_len ? IS_DIGIT (field[0]) : false;
    bool is_final_item = false;
    const char *snip = field;
    unsigned snip_len = 0;
    *n_items = 0;

    for (unsigned i=0; i < field_len; i++) { // one more than field_len - to finalize the last subfield
    
        int next_c = (i==field_len-1) ? 0 : field[i+1];
        bool next_is_digit = IS_DIGIT (next_c);

        bool is_last_char_of_snip = (i==field_len-1)           ? true   // last character of data - definitely also last of this snip
                                  : is_final_item              ? false  // final item - we don't stop until the end of the data
                                  : (segconf.is_bgi_E9L1C3R10 && in_digits != next_is_digit) ? true   // digits / non-digits boundary is always terminating a snip
                                  : is_sep[next_c]             ? true   // a valid separator terminates a snip
                                  :                              false; // no reason to terminate this snip now;
        snip_len++;

        if (!is_last_char_of_snip) continue;

        CompoundItem *ci = &items[(*n_items)++];
        ci->item        = snip;
        ci->item_len    = snip_len; 
        ci->sep         = is_sep[next_c] ? (char)next_c : 0;

        // next, find out if item is an integer
        if (in_digits) {
            unsigned leading_zeros=0;
            for (; leading_zeros < snip_len-1 && snip[leading_zeros]=='0'; leading_zeros++) {} // last digit is never a leading 0, if it is 0, it is the integer

            // now set is_int, value and leading_zero. note that if is_int is false, the other two are meaningless
            ci->is_int = str_get_int (&snip[leading_zeros], snip_len - leading_zeros, &ci->value);
            ci->leading_zeros = leading_zeros;
        }
        else
            ci->is_int = false;

        if (*n_items == MAX_ARRAY_ITEMS-1) is_final_item = true;

        if (ci->sep) i++; // skip separator

        // check for double space separator
        if (ci->sep == ' ' && i < field_len-1 && field[i+1] == ' ') {
            ci->sep2 = ' ';
            i++;
        }
        else
            ci->sep2 = 0;

        if (i < field_len-1) in_digits = IS_DIGIT (field[i+1]);

        snip = &field[i+1];
        snip_len = 0;
    }
}                                      

static void compound_seg_default (VBlock *vb, 
                                  Context *field_ctx, const char *field, unsigned field_len, 
                                  const char *is_sep,   
                                  unsigned nonoptimized_len,      // if non-zero, we don't account for the string given, instead, only for this amount (+add_for_eol)
                                  unsigned add_additional_bytes)  // account for characters in addition to the field
{
    // we use nodes.param in D?ESC contexts to track whether all snips in in this VB are the same
    // defaults to 0 (the same) and we set it to 1 if we encounter a different one
    #define not_all_the_same nodes.param 

    Container con = seg_initialize_container_array (field_ctx->dict_id, true, false); 
    CompoundItem items[MAX_ARRAY_ITEMS];

    compound_field_split (field, field_len, is_sep, items, &con.nitems_lo);
    
    char prefixes[field_len + con.nitems_lo + 2];
    prefixes[0] = prefixes[1] = CON_PREFIX_SEP;
    unsigned prefixes_len = 2;
    unsigned num_seps = 0;

    for (unsigned i=0; i < con.nitems_lo; i++) {
        CompoundItem  *ci = &items[i];
        ContainerItem *CI = &con.items[i];

        // process the subfield that just ended
        Context *item_ctx = ctx_get_ctx (vb, CI->dict_id);
        ASSERT (item_ctx, "item_ctx for %s is NULL", dis_dict_id (CI->dict_id).s);

        item_ctx->st_did_i = field_ctx->did_i;

        // allocate memory if needed
        buf_alloc (vb, &item_ctx->b250, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->b250");

        // if snip is an integer, we store a delta
        char delta_snip[30];
        unsigned original_item_len = ci->item_len;

        if (ci->is_int) {
            delta_snip[0] = SNIP_SELF_DELTA;

            PosType delta = ci->value - item_ctx->last_value.i;

            // note: if all the snips so far in this VB are the same - store just the snip, so that if the 
            // entire b250 is the same, it can be removed
            if (!segconf.sam_is_sorted && (delta || item_ctx->not_all_the_same)) {
                ci->item     = delta_snip;
                ci->item_len = 1 + str_int (delta, &delta_snip[1]);

                item_ctx->flags.store = STORE_INT;
                item_ctx->not_all_the_same = true;
            }
            else
                ci->leading_zeros = 0; // we store the snip as-is

            ctx_set_last_value (vb, item_ctx, ci->value);
        }
        
        if (flag.pair == PAIR_READ_1)
            item_ctx->no_stons = true; // prevent singletons, so pair_2 can compare to us
        
        // we are evaluating but might throw away this snip and use SNIP_PAIR_LOOKUP instead - however, we throw away if its in the pair file,
        // i.e. its already in the dictionary and hash table - so no resources wasted
        WordIndex node_index = ctx_evaluate_snip_seg (vb, item_ctx, ci->item, ci->item_len, NULL);

        // case we are compressing fastq pairs - read 1 is the basis and thus must have a b250 node index,
        // and read 2 might have SNIP_PAIR_LOOKUP
        if (flag.pair == PAIR_READ_2) {

            // if the number of components in the compound is not exactly the same for every line of
            // pair 1 and pair 2 for this vb, the readings from the b250 will be incorrect, causing missed opportunities 
            // for SNIP_PAIR_LOOKUP and hence worse compression. This conditions makes sure this situation
            // doesn't result in an error (TO DO: overcome this, bug 159)
            // note: this can also happen if the there is no item_ctx->pair do it being fully singletons and moved to local 
            if (item_ctx->pair_b250_iter.next_b250 < AFTERENT (uint8_t, item_ctx->pair) &&
                // for pairing to word with SNIP_DELTA, if we have SNIP_PAIR_LOOKUP then all previous lines
                // this VB must have been SNIP_PAIR_LOOKUP as well. Therefore, the first time we encounter an
                // inequality - we stop the pairing going forward till the end of this VB
                !item_ctx->stop_pairing) {
                
                WordIndex pair_word_index = base250_decode (&item_ctx->pair_b250_iter.next_b250, !item_ctx->pair_flags.all_the_same, item_ctx->tag_name);  
                
                if (pair_word_index == WORD_INDEX_ONE_UP) 
                    pair_word_index = item_ctx->pair_b250_iter.prev_word_index + 1;
                
                item_ctx->pair_b250_iter.prev_word_index = pair_word_index;
                
                // note: if the pair word is a singleton in pair_1 file, then pair_word_index will be the index of {SNIP_LOOKUP}
                // rather than the snip (as replaced in ctx_evaluate_snip_merge), therefore this condition will fail. This is quite
                // rare, so not worth handling this case
                if (node_index == pair_word_index) {

                    // discard node_index - decrement its count
                    ctx_decrement_count (vb, item_ctx, node_index);

                    // get new node_index instead
                    item_ctx->pair_b250 = true;
                    static const char lookup_pair_snip[1] = { SNIP_PAIR_LOOKUP };
                    node_index = ctx_evaluate_snip_seg (vb, item_ctx, lookup_pair_snip, 1, NULL);
                } 
                else
                    // To improve: currently, pairing stops at the first non-match
                    item_ctx->stop_pairing = true;
            }
            else
                item_ctx->stop_pairing = true;
        }

        NEXTENT (uint32_t, item_ctx->b250) = node_index;

        item_ctx->txt_len += nonoptimized_len ? 0 : original_item_len;

        // set separators
        CI->separator[0] = ci->sep;
        CI->separator[1] = ci->sep2;
        num_seps += (ci->sep != 0) + (ci->sep2 != 0);

        // add the leading zeros to the container prefixes
        if (ci->is_int && ci->leading_zeros) {
            memset (&prefixes[prefixes_len], '0', ci->leading_zeros);
            prefixes_len += ci->leading_zeros + 1;
            prefixes[prefixes_len-1] = CON_PREFIX_SEP;
        }
        else
            prefixes[prefixes_len++] = CON_PREFIX_SEP;
    }

    container_seg (vb, field_ctx, &con, prefixes, prefixes_len, (nonoptimized_len ? nonoptimized_len : num_seps) + add_additional_bytes);
}

//-------------------------------------------------------------------
// Seggers for specific formats, that perform better than the default
//-------------------------------------------------------------------

static bool __attribute__((unused))compound_is_illumina_7 (STRp(qname))
{
    return str_count_char (STRa(qname), ':') == 6;
}

// Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
static bool __attribute__((unused))compound_seg_illumnina_7 (VBlock *vb, Context *ctx, STRp(qname),
                                      unsigned add_additional_bytes)  // account for characters in addition to the field
{
    str_split (qname, qname_len, 7, ':', item, true);
    if (!n_items) return false;

    uint32_t five_items_len = qname_len - (1 + item_lens[5] + 1 + item_lens[6]);
    seg_by_ctx (vb, STRa(illumina_7_snip), ctx,   2 + add_additional_bytes); // 2= two ':'
    seg_by_ctx (vb, qname, five_items_len, ctx+1, five_items_len); // first 5 items concatenated (as they are correlated)
    seg_by_ctx (vb, STRi(item, 5),         ctx+2, item_lens[5]); 
    seg_by_ctx (vb, STRi(item, 6),         ctx+3, item_lens[6]); 
    
    return true;
}

static bool __attribute__((unused))compound_is_bgi_E9L1C3R10 (STRp(qname))
{
    static const unsigned digit_i[23] = { /*E*/1,2,3,4,5,6,7,8,9, /*L*/11, /*C*/13,14,15, /*R*/17,18,19,20,21,22,23,24,25,26 };

    if (qname_len != 27 || qname[0] != 'E' || qname[10] != 'L' || qname[12] != 'C' || qname[16] != 'R') return false;

    for (unsigned i=0; i < 23; i++)
        if (!IS_DIGIT(qname[digit_i[i]])) return false;

    return true;
}

// BGI: E100020409L1C001R0030000234 (E100020409=Flow cell serial number, L1=Lane 1, C001R003=column 1 row 3, 00q00234=Tile) Also see: https://github.com/IMB-Computational-Genomics-Lab/BGIvsIllumina_scRNASeq
static bool __attribute__((unused)) compound_seg_bgi_E9L1C3R10 (VBlockP vb, ContextP ctx, STRp(qname), 
                                        unsigned add_additional_bytes)
{
    if (!compound_is_bgi_E9L1C3R10 (STRa(qname))) return false;

    seg_by_ctx (vb, STRa(bgi_snip), ctx,   4 + add_additional_bytes); // 4=ELCR
    seg_by_ctx (vb, qname+1,  9,    ctx+1, 9);
    seg_by_ctx (vb, qname+11, 1,    ctx+2, 1);
    seg_by_ctx (vb, qname+13, 3,    ctx+3, 3);
    seg_by_ctx (vb, qname+17, 3,    ctx+4, 3);
    seg_by_ctx (vb, qname+20, 7,    ctx+5, 7);

    return true;
}

void compound_seg (VBlock *vb, 
                   Context *ctx, STRp (qname),
                   const char *default_is_sep,        
                   unsigned nonoptimized_len,      // if non-zero, we don't account for the string given, instead, only for this amount (+add_for_eol)
                   unsigned add_additional_bytes)  // account for characters in addition to the field
{
    START_TIMER;

    bool success = false;
    
    // in sorted SAMs, there is a tiny advantage of copying consecutive identical QNAMEs. In collated mode,
    // we're better off not copying as the entroy added QNAME is more than is saved in Q?NAME 
    if (segconf.sam_is_sorted && vb->line_i && is_same_last_txt (vb, ctx, STRa(qname))) {
        seg_by_ctx (vb, STRa(copy_qname), ctx, (nonoptimized_len ? nonoptimized_len : qname_len) + add_additional_bytes);
        return;
    }

    /* not active yet - specific seggers don't yet perform better than the default (need to add diffs etc)
    if (!nonoptimized_len) {
        if (segconf.is_illumina_7)         success = compound_seg_illumnina_7   (VB, ctx, STRa(qname), add_additional_bytes);
        else if (segconf.is_bgi_E9L1C3R10) success = compound_seg_bgi_E9L1C3R10 (VB, ctx, STRa(qname), add_additional_bytes);
    }*/

    if (!success) compound_seg_default (VB, ctx, STRa(qname), default_is_sep, nonoptimized_len, add_additional_bytes);

    COPY_TIMER (compound_seg);
}

void compound_segconf_test (STRp(qname))
{
    segconf.is_bgi_E9L1C3R10 = segconf.is_bgi_E9L1C3R10 && compound_is_bgi_E9L1C3R10 (STRa(qname));
    segconf.is_illumina_7 = segconf.is_illumina_7 && compound_is_illumina_7 (STRa(qname));
}
