// ------------------------------------------------------------------
//   tokenizer.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "tokenizer.h"
#include "buffer.h"
#include "vblock.h"
#include "strings.h"
#include "container.h"
#include "seg.h"

#define MAX_TOKENS 32

// We break down the field (eg QNAME in SAM/BAM/Kraken or Description in FASTA/FASTQ) into subfields separated by / and/or : - and/or whitespace 
// these are vendor-defined strings.
// Up to MAX_TOKENS subfields are permitted - if there are more, then all the trailing part is just
// consider part of the last component.
// each subfield is stored in its own dictionary
const char sep_with_space[256]    = { [':']=true, [';']=true, ['/']=true, ['|']=true, ['_']=true, ['#']=true, [' ']=true, ['\t']=true, [1]=true };
const char sep_without_space[256] = { [':']=true, [';']=true, ['/']=true, ['|']=true, ['_']=true, ['#']=true, };

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
} Token;

static void tokenizer_split (STRp(field), const char *is_sep, bool split_on_digit_boundary,
                             Token *items, uint8_t *n_items) // out: array of MAX_TOKENS
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
                                  : (split_on_digit_boundary && in_digits != next_is_digit) ? true   // digits / non-digits boundary terminates a snip
                                  : is_sep[next_c]             ? true   // a valid separator terminates a snip
                                  :                              false; // no reason to terminate this snip now;
        snip_len++;

        if (!is_last_char_of_snip) continue;

        Token *ci = &items[(*n_items)++];
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

        if (*n_items == MAX_TOKENS-1) is_final_item = true;

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

Container tokenizer_initialize_container_array (DictId dict_id)
{
    Container con = (Container){ .repeats = 1 };

    for (unsigned i=0; i < MAX_TOKENS; i++) {
        const uint8_t *id = dict_id.id;
        
        char dict_id_str[8] = { id[0], base32(i), id[1], id[2], id[3], id[4], id[5], id[6] };
        
        con.items[i].dict_id = dict_id_make (dict_id_str, 8, DTYPE_1);
    }

    return con;
}

void tokenizer_seg (VBlockP vb, ContextP field_ctx, STRp(field), 
                    const char *is_sep,   
                    unsigned add_additional_bytes)  // account for characters in addition to the field
{
    // we use nodes.param in D?ESC contexts to track whether all snips in in this VB are the same
    // defaults to 0 (the same) and we set it to 1 if we encounter a different one
    #define not_all_the_same nodes.param 

    // when generate the items from "TOKEN" - from from field_ctx's dict_id to avoid clashing
    // with qname in case tokenizer is used as a fallback
    DictId generator_dict_id = dict_id_make ("TOKEN", 5, DTYPE_1); 
    Container con = tokenizer_initialize_container_array (generator_dict_id); 
    Token items[MAX_TOKENS];

    tokenizer_split (field, field_len, is_sep, false, items, &con.nitems_lo);
    
    char prefixes[field_len + con.nitems_lo + 2];
    prefixes[0] = prefixes[1] = CON_PX_SEP;
    unsigned prefixes_len = 2;
    unsigned num_seps = 0;

    for (unsigned i=0; i < con.nitems_lo; i++) {
        Token  *ci = &items[i];
        ContainerItem *CI = &con.items[i];

        // process the subfield that just ended
        Context *item_ctx = ctx_get_ctx (vb, CI->dict_id);
        ASSERT (item_ctx, "item_ctx for %s is NULL", dis_dict_id (CI->dict_id).s);

        item_ctx->st_did_i = field_ctx->did_i;

        // if snip is an integer, we store a delta
        char delta_snip[30];
        unsigned original_item_len = ci->item_len;
        PosType delta;

        #define MAX_TOKENIZER_DETLA 16384 // arbitrary (Illumina ~= 100-800)

        if (ci->is_int && ctx_has_value_in_prev_line_(vb, item_ctx) && ABS((delta = ci->value - item_ctx->last_value.i)) < MAX_TOKENIZER_DETLA) {
            delta_snip[0] = SNIP_SELF_DELTA;

            // note: if all the snips so far in this VB are the same - store just the snip, so that if the 
            // entire b250 is the same, it can be removed
            if (delta || item_ctx->not_all_the_same) {
                ci->item     = delta_snip;
                ci->item_len = 1 + str_int (delta, &delta_snip[1]);

                item_ctx->flags.store = STORE_INT;
                item_ctx->not_all_the_same = true;
            }
            else
                ci->leading_zeros = 0; // we store the snip as-is
        }

        if (ci->is_int) ctx_set_last_value (vb, item_ctx, ci->value);
        
        if (flag.pair == PAIR_READ_1)
            item_ctx->no_stons = true; // prevent singletons, so pair_2 can compare to us
        
        // we are evaluating but might throw away this snip and use SNIP_PAIR_LOOKUP instead - however, we throw away if its in the pair file,
        // i.e. its already in the dictionary and hash table - so no resources wasted
        WordIndex node_index = ctx_create_node_do (vb, item_ctx, ci->item, ci->item_len, NULL);

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
                    node_index = ctx_create_node_do (vb, item_ctx, lookup_pair_snip, 1, NULL);
                } 
                else
                    // To improve: currently, pairing stops at the first non-match
                    item_ctx->stop_pairing = true;
            }
            else
                item_ctx->stop_pairing = true;
        }

        ctx_append_b250 (vb, item_ctx, node_index);

        item_ctx->txt_len += original_item_len;

        // set separators
        CI->separator[0] = ci->sep;
        CI->separator[1] = ci->sep2;
        num_seps += (ci->sep != 0) + (ci->sep2 != 0);

        // add the leading zeros to the container prefixes
        if (ci->is_int && ci->leading_zeros) {
            memset (&prefixes[prefixes_len], '0', ci->leading_zeros);
            prefixes_len += ci->leading_zeros + 1;
            prefixes[prefixes_len-1] = CON_PX_SEP;
        }
        else
            prefixes[prefixes_len++] = CON_PX_SEP;
    }

    container_seg (vb, field_ctx, &con, prefixes, prefixes_len, num_seps + add_additional_bytes);
}
