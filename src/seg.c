// ------------------------------------------------------------------
//   seg.c
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <stdarg.h>
#include "genozip.h"
#include "profiler.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "endianness.h"
#include "file.h"
#include "strings.h"
#include "optimize.h"
#include "random_access.h"
#include "dict_id.h"
#include "base64.h"
#include "piz.h"
#include "zfile.h"
#include "data_types.h"
#include "container.h"
#include "codec.h"
#include "reference.h"
#include "zip.h"
#include "segconf.h"
#include "website.h"
#include "stats.h"
#include "bgzf.h"
#include "dispatcher.h"
#include "buf_list.h"
#include "b250.h"
#include "zip_dyn_int.h"
#include "libdeflate/libdeflate.h"

WordIndex seg_by_ctx_ex (VBlockP vb, STRp(snip), ContextP ctx, uint32_t add_bytes,
                         bool *is_new) // optional out
{
    ASSERTNOTNULL (ctx);

    WordIndex node_index = ctx_create_node_do (VB, ctx, STRa(snip), is_new);

    ASSERT (node_index < ctx->nodes.len32 + ctx->ol_nodes.len32 || node_index == WORD_INDEX_EMPTY || node_index == WORD_INDEX_MISSING, 
            "out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u",  
            ctx->tag_name, node_index, ctx->nodes.len32, ctx->ol_nodes.len32);
    
    b250_seg_append (vb, ctx, node_index);
    ctx->txt_len += add_bytes;
    
    return node_index;
} 

// segs the same node as previous seg
WordIndex seg_known_node_index (VBlockP vb, ContextP ctx, WordIndex node_index, unsigned add_bytes) 
{ 
    b250_seg_append (vb, ctx, node_index);
    ctx->txt_len += add_bytes;
    
    ctx_increment_count (vb, ctx, node_index);

    return node_index;
}

WordIndex seg_duplicate_last (VBlockP vb, ContextP ctx, unsigned add_bytes) 
{ 
    ASSERTISALLOCED (ctx->b250);
    return seg_known_node_index (vb, ctx, b250_seg_get_last (ctx), add_bytes); 
}

#define MAX_ROLLBACK_CTXS ARRAY_LEN(vb->rollback_dids)

// NOTE: this does not save recon_size - if the processing may change recon_size it must be saved and rolled back separately
void seg_create_rollback_point (VBlockP vb, 
                                ContainerP con,         // option 1 - all the items in a container
                                unsigned num_ctxs, ...) // option 2 - explict list of Dids
{
    if (con) num_ctxs = con_nitems (*con);

    ASSERT (num_ctxs <= MAX_ROLLBACK_CTXS, "num_ctxs=%u > MAX_ROLLBACK_CTXS=%u", num_ctxs, MAX_ROLLBACK_CTXS);

    va_list args;
    va_start (args, num_ctxs);

    vb->rback_id++; // new rollback point
    vb->num_rollback_ctxs = 0;

    for (uint32_t i=0; i < num_ctxs; i++) {
        if (con && !con->items[i].dict_id.num) continue;

        vb->rollback_dids[vb->num_rollback_ctxs] = con ? ctx_get_ctx (vb, con->items[i].dict_id)->did_i : va_arg (args, int); // note: 'Did' {aka 'short unsigned int'} is promoted to 'int' when passed through '...'

        if (ctx_set_rollback (vb, CTX(vb->rollback_dids[vb->num_rollback_ctxs]), false)) // added
            vb->num_rollback_ctxs++;
    }
}

// adds a ctx to the current rollback point, if its not already there
void seg_add_ctx_to_rollback_point (VBlockP vb, ContextP ctx)
{
    ASSERT (vb->num_rollback_ctxs < MAX_ROLLBACK_CTXS, "num_ctxs=%u > MAX_ROLLBACK_CTXS=%u", vb->num_rollback_ctxs+1, MAX_ROLLBACK_CTXS);

    vb->rollback_dids[vb->num_rollback_ctxs] = ctx->did_i;
    
    if (ctx_set_rollback (vb, ctx, false)) // added
        vb->num_rollback_ctxs++;
}

void seg_rollback (VBlockP vb)
{
    for (uint32_t i=0; i < vb->num_rollback_ctxs; i++) 
        ctx_rollback (vb, CTX(vb->rollback_dids[i]), false);
}

void seg_simple_lookup (VBlockP vb, ContextP ctx, unsigned add_bytes)
{
    seg_by_ctx (VB, (char[]){ SNIP_LOOKUP }, 1, ctx, add_bytes);
}

void seg_lookup_with_length (VBlockP vb, ContextP ctx, int32_t length/*can be negative*/, unsigned add_bytes)
{
    SNIPi1 (SNIP_LOOKUP, length);
    seg_by_ctx (VB, STRa(snip), ctx, add_bytes);
}

rom seg_get_next_item (VBlockP vb, rom str, int *str_len, 
                       GetNextAllow newline, GetNextAllow tab, GetNextAllow space,
                       unsigned *len, 
                       char *separator,  // optional out
                       bool *has_13,     // out - only modified if '\r' detected ; only needed if newline=GN_SEP
                       rom item_name)
{
    START_TIMER;

    rom c, after;
    for (c = str, after = str + *str_len; c < after; c++) 
        switch (*c) {
            case ' '  : if      (space   == GN_SEP     ) goto sep_found; 
                        else if (space   == GN_FORBIDEN) goto not_found;
                        break;
                       
            case '\t' : if      (tab     == GN_SEP     ) goto sep_found; 
                        else if (tab     == GN_FORBIDEN) goto not_found;
                        break;

            case '\n' : if      (newline == GN_SEP     ) goto sep_found; 
                        else if (newline == GN_FORBIDEN) goto not_found;
                        break;

            default   : break;
        }

not_found: // no sep found in entire string, or forbidden separator encountered
    ABOSEG ("while segmenting %s: expecting a %s%s%s in \"%.*s\"", 
            item_name,
            newline==GN_SEP ? "NEWLINE " : "", tab==GN_SEP ? "TAB " : "", space==GN_SEP ? "\" \" " : "", 
            (int)MIN_(after - str, 1000), str);
    
sep_found:
    *len = c - str;
    if (separator) *separator = *c;
    *str_len -= *len + 1;

    // check for Windows-style '\r\n' end of line 
    if (c != str && *c == '\n' && *(c-1) == '\r') {
        ASSERTNOTNULL (has_13);
        (*len)--;
        *has_13 = true;
    }

    COPY_TIMER (seg_get_next_item);
    return str + *len + 1 + (*c == '\n' && has_13 && *has_13); // beyond the separator
}

// returns first character after current line
rom seg_get_next_line (VBlockP vb, rom str, 
                       int *remaining_len, // in/out
                       unsigned *len,      // out - length of line exluding \r and \n
                       bool must_have_newline, bool *has_13 /* out */, rom item_name)
{
    START_TIMER;

    rom after = str + *remaining_len;
    for (rom s=str; s < after; s++)
        if (*s == '\n') {
                *len = s - str;
                *remaining_len -= *len + 1;

                // check for Windows-style '\r\n' end of line 
                if (s > str && s[-1] == '\r') {
                    (*len)--;
                    *has_13 = true;
                }

                COPY_TIMER (seg_get_next_line);
                return str + *len + *has_13 + 1; // beyond the separator
        }
    
    ASSSEG (*remaining_len, "missing %s field", item_name);

    // if we reached here - line doesn't end with a newline
    ASSSEG (!must_have_newline, "while segmenting %s: expecting a NEWLINE after (showing at most 1000 characters): \"%.*s\"", 
            item_name, MIN_(*remaining_len, 1000), str);

    // we have no newline, but check if last character is a \r
    *has_13 = (after > str) && (after[-1] == '\r');
    *len = *remaining_len - *has_13;
    *remaining_len = 0;

    COPY_TIMER (seg_get_next_line);
    return after;
}

// returns true is value is of type store_type and stored in last_value
bool seg_set_last_txt_store_value (VBlockP vb, ContextP ctx, STRp(value), StoreType store_type)
{
    bool is_value_in_txt_data = value >= B1STtxt &&
                                value <= BLST  (char, vb->txt_data);

    ctx->last_txt = (TxtWord) { .index = is_value_in_txt_data ? BNUMtxt (value) : INVALID_LAST_TXT_INDEX,
                                .len   = value_len };

    bool stored = false;

    if (store_type == STORE_INT) 
        stored = str_get_int (value, value_len, &ctx->last_value.i);
    
    else if (store_type == STORE_FLOAT) {
        char *after; 
        double f = strtod (value, &after); 
        if ((stored = (after == value + value_len))) 
            ctx->last_value.f = f;
    }

    if (stored) {
        ctx->last_line_i   = vb->line_i;    
        ctx->last_sample_i = vb->sample_i;
    }
    else 
        ctx_set_encountered (vb, ctx);

    return stored;
}

WordIndex seg_integer_as_snip_do (VBlockP vb, ContextP ctx, int64_t n, unsigned add_bytes)
{
    char snip[24];
    unsigned snip_len = str_int (n, snip);
    return seg_by_ctx (vb, STRa(snip), ctx, add_bytes);
}

// prepare snips that contain code + dict_id + optional parameter (SNIP_OTHER_LOOKUP, SNIP_OTHER_DELTA, SNIP_REDIRECTION...)
void seg_prepare_snip_other_do (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int64_t int_param, char char_param, 
                                qSTRp(snip)) //  in / out
{
    // make sure we have enough memory
    unsigned required_len = 1/*snip code*/ + 1/*\0*/ + base64_size (DICT_ID_LEN) + (has_parameter ?  11 : 0)/* max length of a int32 -1000000000 */;
    ASSERT (*snip_len >= required_len, "*snip_len=%u, but it needs to be at least %u", *snip_len, required_len);

    snip[0] = snip_code;
    *snip_len = 1 + base64_encode (other_dict_id.id, DICT_ID_LEN, &snip[1]);

    if (has_parameter) {
        if (char_param) 
            snip[(*snip_len)++] = char_param;
        else
            *snip_len += str_int (int_param, &snip[*snip_len]);
    }
}

void seg_prepare_multi_dict_id_special_snip (uint8_t special_code, unsigned num_dict_ids, DictId *dict_ids,
                                             qSTRp (out_snip)) // in/out - allocated by caller
{
    // prepare snip - a string consisting of num_dict_ids x { VCF_SPECIAL_MUX_BY_DOSAGE , base64(dict_id.id) }   
    char *snip = out_snip;
    snip[0] = SNIP_SPECIAL;
    unsigned snip_len = 1;

    for (int i=0; i < num_dict_ids; i++) {
        snip += snip_len;
        snip_len = out_snip + *out_snip_len - snip;
        seg_prepare_snip_other_do (i ? '\t' : special_code, dict_ids[i], 0, 0, 0, snip, &snip_len);
    }

    *out_snip_len = (snip + snip_len) - out_snip;
}

// prepare snip of A - B
void seg_prepare_array_dict_id_special_snip (int num_dict_ids, DictId *dict_ids, uint8_t special_code, qSTRp(snip))
{
    // make sure we have enough memory
    unsigned required_len = 2 + base64_size (num_dict_ids * DICT_ID_LEN); 
    ASSERT (*snip_len >= required_len, "*snip_len=%u, but it needs to be at least %u", *snip_len, required_len);

    snip[0] = SNIP_SPECIAL;
    snip[1] = special_code;
    
    *snip_len = 2 + base64_encode ((uint8_t *)dict_ids, num_dict_ids * DICT_ID_LEN, &snip[2]);
}

void seg_mux_display (MultiplexerP mux)
{
    iprintf ("ctx=%s num_channels=%u snip_len=%u snip=%.100s", 
             mux->ctx->tag_name, mux->num_channels, mux->snip_len, MUX_SNIP(mux));

    for (unsigned i=0; i < MIN_(mux->num_channels, 1024); i++)
        iprintf ("[%u] dict_id=%s channel_ctx=%p\n", i, dis_dict_id (mux->dict_ids[i]).s, MUX_CHANNEL_CTX(mux)[i]);
}

void seg_mux_init (VBlockP vb, ContextP ctx, unsigned num_channels, uint8_t special_code, 
                   bool no_stons, MultiplexerP mux) // optional - a string with num_channels unique characters
{
    // note: if we ever need more than 256, we need to update the dict_id generation formula below
    ASSERT (num_channels <= 256, "num_channels=%u exceeds maximum of 256", num_channels);

    bytes id = ctx->dict_id.id;
    DictId dict_id_template = (DictId){ .id = { id[0], 0, id[1], id[2], id[3], id[4], id[5], id[6] } };
    unsigned id_len = 1 + strnlen ((rom)ctx->dict_id.id, 7);

    mux->ctx          = ctx;
    mux->num_channels = num_channels;
    mux->no_stons     = no_stons;
    mux->special_code = special_code;
    
    // caller gave st_did_i, but st_did_i might itself consolidate 
    if (ctx->st_did_i == DID_NONE) ctx->is_stats_parent = true;

    // calculate dict_ids, eg: PL -> PlL, PmL, PnL, PoL
    for (int i=0; i < num_channels; i++) {
        mux->dict_ids[i] = dict_id_template;
        mux->dict_ids[i].id[1] = ((uint16_t)'0' + (uint16_t)i) & 0xff;
        if (id_len <= 5) str_int (i, (char*)&mux->dict_ids[i].id[id_len]); // add numerical digits for ease of debugging if we have room
    }

    // prepare snip - a string consisting of num_channels x { special_code or \t , base64(dict_id.id) }
    mux->snip_len = MUX_SIZEOF_SNIP(mux);
    seg_prepare_multi_dict_id_special_snip (special_code, num_channels, mux->dict_ids, MUX_SNIP(mux), &mux->snip_len);

    // note: snips prepared with seg_prepare_multi_dict_id_special_snip can add 
    // an optional 3 chars + SNIP_SPECIAL in order to display the dicts in --show-snip 
    memcpy (&MUX_SNIP(mux)[mux->snip_len], ((char[]){'\t', 'M','U','X', '0', SNIP_SPECIAL }), 6); // for --show-snip
    mux->snip_len += 6;

    // seg_mux_display (mux);
}

ContextP seg_mux_get_channel_ctx (VBlockP vb, Did did_i, MultiplexerP mux, uint32_t channel_i)
{
    ASSERT (mux->ctx, "%s: mux of %s is not initialized", LN_NAME, CTX(did_i)->tag_name);
    
    ASSERT (channel_i < mux->num_channels, "channel_i=%u out of range, multiplexer of \"%s\" has only %u channels",
            channel_i, mux->ctx->tag_name, mux->num_channels);

    ContextP channel_ctx = ctx_get_ctx (vb, mux->dict_ids[channel_i]);

    if (!channel_ctx->is_initialized) { // note: context may exist (cloned from z_ctx) but not yet initialized

        MUX_CHANNEL_CTX(mux)[channel_i] = channel_ctx;

        channel_ctx->is_initialized = true;
        channel_ctx->flags.store    = mux->ctx->flags.store;
        channel_ctx->no_stons       = mux->no_stons; // not inherited from parent
        channel_ctx->st_did_i       = (mux->ctx->st_did_i == DID_NONE) ? mux->ctx->did_i : mux->ctx->st_did_i;
        channel_ctx->flags.ctx_specific_flag = mux->ctx->flags.ctx_specific_flag; // inherit - needed for "same_line"        

        if (channel_ctx->ltype == LT_SINGLETON) { // individual channel might have been already set in seg_initialize to something else 
            if (IS_LT_DYN (mux->ctx->ltype))
                dyn_int_init_ctx (vb, channel_ctx, 0);
            else
                channel_ctx->ltype = mux->ctx->ltype; 
        }
    }

    return channel_ctx;
}

// scans a pos field - in case of non-digit or not in the range [0,MAX_POS], either returns -1
// (if allow_nonsense) or errors
static PosType64 seg_scan_pos_snip (VBlockP vb, rom snip, unsigned snip_len, bool zero_is_bad, 
                                  SegError *err) // out - if NULL, exception if error
{
    PosType64 value=0;
    bool is_int = str_get_int (snip, snip_len, &value); // value unchanged if not integer

    if (is_int && value >= !!zero_is_bad && value <= MAX_POS) {
        *err = ERR_SEG_NO_ERROR;
        return value; // all good
    }

    ASSSEG (err, "position field must be an integer number between %u and %"PRId64", seeing: %.*s", 
            !!zero_is_bad, MAX_POS, STRf(snip));

    *err = is_int ? ERR_SEG_OUT_OF_RANGE : ERR_SEG_NOT_INTEGER;
    return value; // in- or out-of- range integer, or 0 if value is not integer
}

// returns POS value if a valid pos, or 0 if not
PosType64 seg_pos_field (VBlockP vb, 
                         Did snip_did_i,    // mandatory: the ctx the snip belongs to
                         Did base_did_i,    // mandatory: base for delta
                         unsigned opt,      // a combination of SPF_* options
                         char missing,      // a character allowed (meaning "missing value"), segged as SNIP_DONT_STORE
                         STRp(pos_str),     // option 1
                         PosType64 this_pos,// option 2
                         unsigned add_bytes)     
{
    ContextP snip_ctx = CTX(snip_did_i);
    ContextP base_ctx = CTX(base_did_i);

    // note: *_seg_initialize is required to set both snip_ctx and base_ctx: no_stons=true, store=STORE_INT

    SegError err = ERR_SEG_NO_ERROR;
    if (pos_str) { // option 1
        if (str_is_1char (pos_str, missing)) 
            err = ERR_SEG_NOT_INTEGER; // not an error, just so that we seg this as SNIP_DONT_STORE
        
        else {
            this_pos = seg_scan_pos_snip (vb, STRa(pos_str), IS_FLAG (opt, SPF_ZERO_IS_BAD), &err);
            ASSERT (IS_FLAG (opt, SPF_BAD_SNIPS_TOO) || !err, "%s: invalid value %.*s in %s", 
                    LN_NAME, STRf(pos_str), snip_ctx->tag_name);
        }

        // we accept out-of-range integer values for non-self-delta
        if (snip_did_i != base_did_i && err == ERR_SEG_OUT_OF_RANGE) err = ERR_SEG_NO_ERROR;

        if (err) {
            SAFE_ASSIGN (pos_str-1, SNIP_DONT_STORE);
            seg_by_ctx (VB, pos_str-1, pos_str_len+1, snip_ctx, add_bytes); 
            SAFE_RESTORE;
            snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
            return 0; // invalid pos
        }
    }

    else { // option 2 
        // check out-of-range for self-delta
        if (snip_did_i == base_did_i && (this_pos < 0 || this_pos > MAX_POS)) {
            err = ERR_SEG_OUT_OF_RANGE;
            SNIPi1 (SNIP_DONT_STORE, this_pos);
            seg_by_ctx (VB, STRa(snip), snip_ctx, add_bytes); 
            snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
            return 0; // invalid pos
        }
    }

    PosType64 pos_delta = this_pos - base_ctx->last_value.i;
    
    ctx_set_last_value (vb, snip_ctx, this_pos);

    // store the delta in last_delta only if we're also putting in the b250
    snip_ctx->last_delta = pos_delta;
    
    // if the delta is the negative of the previous delta (as happens in collated BAM files with the second line in
    // each pair of lines) - we just store an empty snippet
    bool is_negated_last = (snip_ctx->last_delta && snip_ctx->last_delta == -pos_delta);

    // case: add a delta
    if (!is_negated_last) {
        
        char pos_delta_str[16]; // enough for the base64
        unsigned total_len = sizeof (pos_delta_str);

        if (base_ctx == snip_ctx) {
            pos_delta_str[0] = SNIP_SELF_DELTA;
            total_len = 1;
        }
        else 
            seg_prepare_snip_other_do (SNIP_OTHER_DELTA, base_ctx->dict_id, false, 0, 0, pos_delta_str, &total_len);

        // backcomp note: starting 15.0.27 we store the delta in local. up to 15.0.26, local 
        // was used to store full (not delta) values, if the delta went beyond MAX_DELTA, while
        // the delta was stored in the snip as a textual integer
        pos_delta_str[total_len++] = '$'; 

        dyn_int_append (vb, snip_ctx, pos_delta, 0); 

        seg_by_ctx (VB, pos_delta_str, total_len, snip_ctx, add_bytes);
    }

    // case: the delta is the negative of the previous delta - add a SNIP_SELF_DELTA with no payload - meaning negated delta
    else {
        char negated_delta = SNIP_SELF_DELTA; // no delta means negate the previous delta
        seg_by_ctx (VB, &negated_delta, 1, snip_ctx, add_bytes);
        snip_ctx->last_delta = 0; // no negated delta next time
    }

    return this_pos;
}

bool seg_pos_field_cb (VBlockP vb, ContextP ctx, STRp(pos_str), uint32_t repeat)
{
    seg_pos_field (vb, ctx->did_i, ctx->did_i, 0, 0, STRa(pos_str), 0, pos_str_len);
    return true; // segged successfully
}

// note: caller must set ctx->ltype=LT_DYN_INT*
void seg_integer (VBlockP vb, ContextP ctx, int64_t n, bool with_lookup, unsigned add_bytes)
{
    ctx_set_last_value (vb, ctx, n);

    dyn_int_append (vb, ctx, n, add_bytes);

    if (with_lookup)     
        seg_by_ctx (vb, (char[]){ SNIP_LOOKUP }, 1, ctx, 0);
}

// returns true if it was an integer
bool seg_integer_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned add_bytes)
{
    if (ctx->nothing_char && value_len == 1 && *value == ctx->nothing_char) { 
        dyn_int_append_nothing_char (vb, ctx, add_bytes);
        seg_by_ctx (vb, (char[]){ SNIP_LOOKUP }, 1, ctx, 0);

        return false; // not an int
    }

    int64_t n;
    bool is_int = (ctx->ltype == LT_DYN_INT_h) ? (value[0] != '0' && str_get_int_hex (STRa(value), true, false, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                : (ctx->ltype == LT_DYN_INT_H) ? (value[0] != '0' && str_get_int_hex (STRa(value), false, true, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                :                                str_get_int (STRa(value), &n);

    // case: its an integer
    if (is_int && (IS_LT_DYN(ctx->ltype) || ctx->ltype == LT_SINGLETON)) {         
        seg_integer (vb, ctx, n, true, add_bytes);
        return true;
    }

    // case: non-numeric snip 
    else { 
        seg_by_ctx (VB, STRa(value), ctx, add_bytes);
        return false;
    }
}

// segs a fixed width, 0-padded, non-negative integer - decimal, hex or HEX
void seg_numeric_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || IS_LT_DYN (ctx->ltype), "Expecting %s.ltype to be DYN_INT* but found %s", ctx->tag_name, lt_name(ctx->ltype));
#endif

    if (segconf.running) goto fallback;

    bool has_zero_x = (ctx->ltype == LT_DYN_INT_h || ctx->ltype == LT_DYN_INT_H) && value_len > 2 && value[0] == '0' && value[1] == 'x';
    if (has_zero_x) STRinc (value, 2); // skip initial "0x"

    int64_t n;
    bool is_numeric = (ctx->ltype == LT_DYN_INT_h) ? (str_get_int_hex (STRa(value), true, false, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                    : (ctx->ltype == LT_DYN_INT_H) ? (str_get_int_hex (STRa(value), false, true, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                    :                                 str_get_int_dec (STRa(value), (uint64_t*)&n);

    // case: its an integer
    if (is_numeric) { 
        seg_integer (vb, ctx, n, false, add_bytes);
        seg_by_ctx (vb, (char[]){ SNIP_NUMERIC, '0'+ (ctx->ltype - LT_DYN_INT), '0' + value_len, 'x' }, 3 + has_zero_x, ctx, 0);
    }

    // case: non-numeric snip
    else fallback: 
        seg_by_ctx (VB, STRa(value), ctx, add_bytes);
}

bool seg_integer_or_not_cb (VBlockP vb, ContextP ctx, STRp(int_str), uint32_t repeat)
{
    seg_integer_or_not (vb, ctx, STRa(int_str), int_str_len);
    return true; // segged successfully
}

// seg a textual snip that is expected to be a float, so that reconstruction reconstructs
// the exact textual representation without floating point issues
bool seg_float_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned add_bytes)
{
    // TO DO: implement reconstruction in reconstruct_one_snip-SNIP_LOOKUP
    char snip[1 + FLOAT_FORMAT_LEN];
    unsigned format_len;

    // case: its an float
    if ((ctx->ltype == LT_FLOAT32 || ctx->ltype == LT_SINGLETON) &&
        str_get_float (STRa(value), &ctx->last_value.f, &snip[1], &format_len)) {

        // verify that we can reconstruct the number precisely (not always the case, 
        // as the textual float may exceed float precision)
        float f = (float)ctx->last_value.f; // 64bit -> 32bit
        char recon[value_len+1];
        sprintf (recon, &snip[1], f);

        if (memcmp (value, recon, value_len)) goto fallback;

        if (ctx->ltype == LT_SINGLETON)
            ctx->ltype = LT_FLOAT32; // set only upon storing the first number - if there are no numbers, leave it as LT_SINGLETON so it can be used for singletons

        // add to local
        buf_alloc (vb, &ctx->local, 1, vb->lines.len, float, CTX_GROWTH, CTX_TAG_LOCAL);
        BNXT (float, ctx->local) = f;
        
        // add to b250
        snip[0] = SNIP_LOOKUP;
        seg_by_ctx (VB, snip, 1 + format_len, ctx, add_bytes);

        return true;
    }

    // case: non-float snip
    else fallback: { 
        seg_by_ctx (VB, STRa(value), ctx, add_bytes);
        return false;
    }
}

void seg_delta_vs_other_do (VBlockP vb, ContextP ctx, ContextP other_ctx, 
                            STRp(value), // if value==NULL, we use value_n
                            int64_t value_n,
                            int64_t max_delta /* max abs value of delta - beyond that, seg as is, ignored if < 0 */,
                            unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || 
            (ctx->flags.store == STORE_INT && (other_ctx->flags.store == STORE_INT || other_ctx->flags.store == STORE_FLOAT)), 
            "expecting ctx=%s and other_ctx=%s to have STORE_INT", ctx->tag_name, other_ctx->tag_name);
#endif

    ASSERTNOTNULL (other_ctx);

    if (value && !str_get_int (STRa(value), &value_n)) {
        seg_by_ctx (VB, STRa(value), ctx, add_bytes);
        return;
    }

    int64_t other_value_n = (other_ctx->flags.store == STORE_FLOAT) 
        ? (int64_t)other_ctx->last_value.f // simple cast to int (not rounding) as in reconstruct_from_delta
        : other_ctx->last_value.i;
    
    int64_t delta = value_n - other_value_n; 
    if (max_delta >= 0 && ABS(delta) > max_delta) 
        seg_integer (vb, ctx, value_n, true, add_bytes);

    else {
        SNIP(32);
        seg_prepare_snip_other (SNIP_OTHER_DELTA, other_ctx->dict_id, false, 0, snip);
        snip[snip_len++] = '$'; // lookup delta from local
        seg_by_ctx (VB, STRa(snip), ctx, add_bytes);

        dyn_int_append (vb, ctx, delta, 0); 

        ctx->last_delta = delta;
    }

    ctx_set_last_value (vb, ctx, value_n);
}

// note: seg_initialize should set STORE_INT and LT_DYN* for this ctx
WordIndex seg_self_delta (VBlockP vb, ContextP ctx, int64_t value, 
                          char format,        // 'd', 'x' or 'X'. 0 is normal, unformatted, number 
                          unsigned fixed_len, // >0 means zero-padded fixed len (legal values: 0-190) 
                          uint32_t add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || ctx->flags.store == STORE_INT, "expecting %s to have store=STORE_INT", ctx->tag_name);
    ASSERT (segconf.running || IS_LT_DYN (ctx->ltype), 
            "%s: Expecting %s.ltype=LT_DYN_INT* but found %s", LN_NAME, ctx->tag_name, lt_name (ctx->ltype));
    ASSERT (fixed_len <= 190, "fixed_len=%u is larger than 190", fixed_len);
#endif

    char delta_snip[8];
    delta_snip[0] = SNIP_SELF_DELTA;
    if (format) {
        delta_snip[1] = format; // currently only 'x' (lower case hex) is supported
        if (fixed_len) delta_snip[2] = 'A' + fixed_len; // fixed_len=190 -> snip=255
    }

    unsigned delta_snip_len = 1 + !!format + !!fixed_len;
    delta_snip[delta_snip_len++] = '$'; // lookup delta from local

    ctx->last_delta = value - ctx->last_value.i;

    dyn_int_append (vb, ctx, ctx->last_delta, 0); 

    ctx_set_last_value (vb, ctx, value);

    return seg_by_ctx (VB, STRa(delta_snip), ctx, add_bytes);
}

// an array or array of arrays
// note: if container_ctx->flags.store=STORE_INT, container_ctx->last_value.i will be set to sum of integer
// elements including recursively from sub-arrays (non-integer elements will be ignored)
WordIndex seg_array (VBlockP vb, ContextP container_ctx, Did stats_conslidation_did_i, 
                     rom value, int32_t value_len, // must be signed
                     char sep, 
                     char subarray_sep,         // if non-zero, will attempt to find internal arrays
                     bool use_integer_delta,    // first item stored as is, subsequent items stored as deltas
                     StoreType store_in_local,  // STORE_INT/FLOAT - if its an integer/float, store in local (with LOOKUP) instead of a snip
                     DictId arr_dict_id,        // item dict_id, DICT_ID_NONE if we don't care what it is 
                     int add_bytes)             // account for this much
{
    MiniContainerP con;
    ContextP arr_ctx;
    int additional_bytes = add_bytes - value_len;

    // first use in this VB - prepare context where array elements will go in
    if (!container_ctx->con_cache.len32) {
        if (!arr_dict_id.num)        
            arr_dict_id = sub_dict_id (container_ctx->dict_id, '0');
        
        buf_alloc (vb, &container_ctx->con_cache, 0, 1, MiniContainer, 1, CTX_TAG_CON_CACHE);

        con = B1ST (MiniContainer, container_ctx->con_cache);
        *con = (MiniContainer){ .nitems_lo = 1, 
                                .drop_final_repsep = true,
                                .repsep    = { sep },
                                .items[0].dict_id = arr_dict_id }; // only one item

        arr_ctx = ctx_get_ctx (vb, arr_dict_id);
    
        if (arr_ctx->st_did_i == DID_NONE) // first time
            arr_ctx->st_did_i = container_ctx->st_did_i = stats_conslidation_did_i;

        if (use_integer_delta || container_ctx->flags.store == STORE_INT) 
            arr_ctx->flags.store = STORE_INT;
    }
    else { 
        con         = B1ST (MiniContainer, container_ctx->con_cache);
        arr_dict_id = con->items[0].dict_id;
        arr_ctx     = ctx_get_ctx (vb, arr_dict_id);
    }

    // case: value ends with a separator
    if (value[value_len-1] == sep) {
        con->drop_final_repsep = false;
        value_len--;
    }
    else
        con->drop_final_repsep = true;

    // count repeats (1 + number of seperators)
    con->repeats = 1 + str_count_char (STRa(value), sep);

    if (container_ctx->flags.store == STORE_INT) 
        ctx_set_last_value (vb, container_ctx, (int64_t)0);

    for (uint32_t i=0; i < con->repeats; i++) { // value_len will be -1 after last number

        rom this_item = value;
        int64_t this_item_value=0;

        bool is_subarray = false;
        for (; value_len && *value != sep; value++, value_len--) 
            if (*value == subarray_sep) 
                is_subarray = true;

        unsigned this_item_len  = (unsigned)(value - this_item);

        // case: it has a sub-array
        if (is_subarray) {
            seg_array (vb, arr_ctx, stats_conslidation_did_i, STRa(this_item), subarray_sep, 0, use_integer_delta, store_in_local, DICT_ID_NONE, this_item_len);
            this_item_value = arr_ctx->last_value.i;
        }

        // case: its an scalar (we don't delta arrays that have sub-arrays and we don't delta the first item)
        else if (!use_integer_delta || subarray_sep || i==0) {
            
            if (store_in_local == STORE_INT) {
                bool is_int = seg_integer_or_not (vb, arr_ctx, STRa(this_item), this_item_len);
                if (is_int) this_item_value = arr_ctx->last_value.i;
            }

            else if (store_in_local == STORE_FLOAT) 
                seg_float_or_not (vb, arr_ctx, STRa(this_item), this_item_len);
            
            else {
                seg_by_ctx (VB, STRa(this_item), arr_ctx, this_item_len);
                str_get_int (STRa(this_item), &arr_ctx->last_value.i); // sets last_value only if it is indeed an integer
            }
        }

        // case: delta of 2nd+ item
        else if (str_get_int (STRa(this_item), &this_item_value)) 
            seg_self_delta (vb, arr_ctx, this_item_value, 0, 0, this_item_len);

        // non-integer that cannot be delta'd - store as-is
        else 
            seg_by_ctx (VB, STRa(this_item), arr_ctx, 0);

        if (container_ctx->flags.store == STORE_INT)
            container_ctx->last_value.i += this_item_value;

        value_len--; // skip separator
        value++;
    }

    return container_seg (vb, container_ctx, (ContainerP)con, 0, 0, con->repeats - con->drop_final_repsep + additional_bytes); // acount for separators and additional bytes
}


void seg_array_by_callback (VBlockP vb, ContextP container_ctx, STRp(arr), char sep, SegCallback item_seg, unsigned add_bytes)
{
    MiniContainerP con;
    ContextP arr_ctx;

    if (!arr_len) {
        seg_by_ctx (vb, "", 0, container_ctx, 0);
        return;
    }

    // first use in this VB - prepare context where array elements will go in
    if (!container_ctx->con_cache.len32) {
        DictId arr_dict_id = sub_dict_id (container_ctx->dict_id, '0');
        
        buf_alloc (vb, &container_ctx->con_cache, 0, 1, MiniContainer, 1, CTX_TAG_CON_CACHE);

        con = B1ST (MiniContainer, container_ctx->con_cache);
        *con = (MiniContainer){ .nitems_lo = 1, 
                                .repsep[0] = sep,
                                .items[0].dict_id = arr_dict_id }; // only one item

        ctx_consolidate_stats_(vb, container_ctx, (ContainerP)con);
    }
    else 
        con = B1ST (MiniContainer, container_ctx->con_cache);

    arr_ctx = ctx_get_ctx (vb, con->items[0].dict_id);

    con->drop_final_repsep = (arr_len && arr[arr_len-1] != sep);

    str_split (arr, arr_len - !con->drop_final_repsep, 0, sep, item, false);

    unsigned sum_lens = 0;
    for (uint32_t i=0; i < n_items; i++) {  // value_len will be -1 after last number
        item_seg (VB, arr_ctx, STRi(item,i), i); // return value ignored - must be true
    
        sum_lens += item_lens[i];
    }

    con->repeats = n_items;
    container_seg (vb, container_ctx, (ContainerP)con, 0, 0, add_bytes - sum_lens); // acount for separators and additional bytes
}

// segs an 2-dimentional array (e.g. "0|0|0,1|2|3") of uint32_t so that it is stored in local transposed. 
// The number of columns (3 in this example) must be identical for the entire VB or its an seg error.
// its also a seg error if any item is not a uint32 integer. The number of lines may vary.
// note: dyn_int_transpose might later reduce it to uint16 or uint8.
void seg_uint32_matrix (VBlockP vb, ContextP container_ctx, Did stats_conslidation_did_i, 
                        STRp(value), char row_sep/*255 if always single row*/, char col_sep, bool is_transposed, int add_bytes)
{
    bytes id = container_ctx->dict_id.id;
    DictId arr_dict_id = { .id = { id[0], 
                                   ((id[1]+1) % 256) | 128, // different IDs for top array, subarray and items
                                   id[2], id[3], id[4], id[5], id[6], id[7] } };

    ContextP arr_ctx = ctx_get_ctx (vb, arr_dict_id);

    // case: first call in this VB - initialize
    if (!arr_ctx->local.n_cols) {
        ctx_consolidate_stats (vb, stats_conslidation_did_i, container_ctx->did_i, arr_ctx->did_i, DID_EOL);

        dyn_int_init_ctx (vb, arr_ctx, 0);
        arr_ctx->flags.store = STORE_INT; // also need for BAM translation
        arr_ctx->dyn_transposed = is_transposed;

        rom first_row_sep = memchr (value, row_sep, value_len);
        
        arr_ctx->local.n_cols = 1 + str_count_char (value, first_row_sep ? (first_row_sep - value) : value_len, col_sep); // n_cols

        ASSERT (!is_transposed || arr_ctx->local.n_cols <= 255, "%s: n_cols=%d greater than max allowed 255", VB_NAME, (int)arr_ctx->local.n_cols);
    }

    seg_create_rollback_point (vb, NULL, 2, container_ctx->did_i, arr_ctx->did_i);

    // seg all integers into matrix
    str_split (value, value_len, 0, row_sep, row, false);
    if (!n_rows) goto badly_formatted;
    uint32_t n_cols_ = arr_ctx->local.n_cols;

    for (int r=0; r < n_rows; r++) {
        str_split (rows[r], row_lens[r], n_cols_, col_sep, col, true);
        if (!n_cols) goto badly_formatted;

        for (int c=0; c < n_cols; c++) {
            int64_t item;
            if (!str_get_int (cols[c], col_lens[c], &item)) goto badly_formatted;
            
            dyn_int_append (vb, arr_ctx, item, 0);
        }
    }

    unsigned account_in_arr_ctx = value_len - (n_rows * n_cols_ - 1);
    arr_ctx->txt_len += account_in_arr_ctx; // entire string except separators

    // finally, the Container snip itself - we attempt to use the known node_index if it is cached in con_index

    // in our specific case, element i of con_index contains the node_index of the snip of the container with i repeats, or WORD_INDEX_NONE.
    if (container_ctx->con_index.len32 < n_rows+1) {
        buf_alloc_255 (vb, &container_ctx->con_index, 0, n_rows+1, WordIndex, 0, CTX_TAG_CON_INDEX);
        container_ctx->con_index.len32 = container_ctx->con_index.len32 / sizeof (WordIndex);
    }

    WordIndex node_index = *B(WordIndex, container_ctx->con_index, n_rows);
    unsigned account_in_container_ctx = add_bytes - account_in_arr_ctx;

    // case: first container with this many repeats - seg and add to cache
    if (node_index == WORD_INDEX_NONE) {
        ASSSEG (n_cols_ <= MEDIUM_CON_NITEMS, "expecting n_cols=%u <= MEDIUM_CON_NITEMS=%u", n_cols_, MEDIUM_CON_NITEMS);

        MediumContainer con = { .nitems_lo           = n_cols_,
                                .repeats             = n_rows,
                                .repsep[0]           = row_sep,
                                .drop_final_repsep   = true,
                                .drop_final_item_sep = true     };

        ContainerItem con_item = { .dict_id = arr_ctx->dict_id, .separator[0] = col_sep };
        for (uint32_t c=0; c < n_cols_; c++)
            con.items[c] = con_item;

        *B(WordIndex, container_ctx->con_index, n_rows) = 
            container_seg (vb, container_ctx, (ContainerP)&con, NULL, 0, account_in_container_ctx);
    }

    // case: we already know the node index of the container with this many repeats
    else 
        seg_known_node_index (vb, container_ctx, node_index, account_in_container_ctx);

    return;

badly_formatted:
    // roll back all the changed data
    seg_rollback (vb);

    // now just seg the value
    seg_by_ctx (VB, STRa(value), container_ctx, add_bytes);
}

bool seg_do_nothing_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    return true; // "segged successfully" - do nothing
}

bool seg_struct (VBlockP vb, ContextP ctx, MediumContainer con, STRp(snip), 
                 const SegCallback *callbacks, // optional - either NULL, or contains a seg callback for each item (any callback may be NULL)
                 unsigned add_bytes,
                 bool account_in_subfields)    // true if to account in subfields, false if in parent
{
    ASSERT0 (con.repeats==1 && !con.repsep[0], "expecting con.repeats==1 and no repsep");

    if (!ctx->is_stats_parent) 
        ctx_consolidate_stats_(vb, ctx, (ContainerP)&con);

    seg_create_rollback_point (vb, (ContainerP)&con, 0);

    str_split_by_container (snip, snip_len, &con, NULL, 0, item, NULL);
    
    if (n_items != con.nitems_lo) 
        goto badly_formatted;

    for (uint32_t i=0; i < n_items; i++) {
        ContextP item_ctx = ctx_get_ctx (vb, con.items[i].dict_id);

        if (callbacks && callbacks[i]) {
            if (!callbacks[i] (vb, item_ctx, STRi(item,i), 0))
                goto badly_formatted;
        }
        
        else if (IS_LT_DYN (item_ctx->ltype))
            seg_integer_or_not (VB, item_ctx, STRi(item,i), account_in_subfields ? item_lens[i] : 0);
        
        else
            seg_by_ctx (VB, STRi(item,i), item_ctx, account_in_subfields ? item_lens[i] : 0);
    }

    // finally, the Container snip itself - we attempt to use the known node_index if it is cached in con_index

    buf_alloc_exact_255 (vb, ctx->con_index, 1, WordIndex, CTX_TAG_CON_INDEX);

    // count printable item separators
    unsigned num_printable_separators=0;
    for (unsigned i=0; i < con.nitems_lo; i++) 
        num_printable_separators += is_printable[con.items[i].separator[0]] + is_printable[con.items[i].separator[1]];

    WordIndex node_index = *B1ST(WordIndex, ctx->con_index);
    unsigned account_for = account_in_subfields ? (num_printable_separators + ((int)add_bytes - snip_len)) 
                                                : add_bytes;

    // case: first time - seg and add to cache
    if (node_index == WORD_INDEX_NONE) 
        *B1ST(WordIndex, ctx->con_index) = container_seg (vb, ctx, (ContainerP)&con, NULL, 0, account_for);
    
    // case: we already know the node index of the container with this many repeats
    else 
        seg_known_node_index (vb, ctx, node_index, account_for);

    return true;

badly_formatted:
    // roll back all the changed data
    seg_rollback (vb);

    // now just seg the entire snip
    seg_by_ctx (VB, STRa(snip), ctx, add_bytes); 

    return false; // not segged as a container
}                           

// a field that looks like: "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
// we have an array (2 entires in this example) of items (4 in this examples) - the entries are separated by comma and the items by space
// observed in Ensembel generated GVF: Variant_effect, sift_prediction, polyphen_prediction, variant_peptide
// The last item is treated as an ENST_ID (format: ENST00000399012) while the other items are regular dictionaries
// the names of the dictionaries are the same as the ctx, with the 2nd character replaced by 1,2,3...
// the field itself will contain the number of entries
int32_t seg_array_of_struct (VBlockP vb, ContextP ctx, MediumContainer con, STRp(snip), 
                             const SegCallback *callbacks, // optional - either NULL, or contains a seg callback for each item (any callback may be NULL)
                             SplitCorrectionCallback split_correction_callback, 
                             unsigned add_bytes)
{
    if (!ctx->is_stats_parent) 
        ctx_consolidate_stats_(vb, ctx, (ContainerP)&con);

    seg_create_rollback_point (vb, (ContainerP)&con, 0);

    // get repeats
    str_split (snip, snip_len, 0, con.repsep[0], repeat, false);

    // if we don't have con.drop_final_repsep, the last "repeat" should be zero length and removed
    if (!con.drop_final_repsep) {
        if (repeat_lens[n_repeats-1]) goto badly_formatted; 
        n_repeats--;
    }
    
    if (split_correction_callback)
        split_correction_callback (&n_repeats, repeats, repeat_lens);

    con.repeats = n_repeats;

    ASSSEG (n_repeats <= CONTAINER_MAX_REPEATS, "exceeded maximum repeats allowed (%u) while parsing %s",
            CONTAINER_MAX_REPEATS, ctx->tag_name);

    for (uint32_t r=0; r < n_repeats; r++) {

        // get items in each repeat 
        str_split_by_container (repeats[r], repeat_lens[r], &con, NULL, 0, item, NULL);
        
        if (n_items != con.nitems_lo) 
            goto badly_formatted;

        for (uint32_t i=0; i < n_items; i++) {
            ContextP item_ctx = ctx_get_ctx (vb, con.items[i].dict_id);
            
            if (callbacks && callbacks[i]) {
                if (!callbacks[i] (vb, item_ctx, STRi(item,i), r))
                    goto badly_formatted;
            }
            else
                seg_by_ctx (VB, STRi(item,i), item_ctx, item_lens[i]);
        }
    }

    // finally, the Container snip itself - we attempt to use the known node_index if it is cached in con_index

    // in our specific case, element i of con_index contains the node_index of the snip of the container with i repeats, or WORD_INDEX_NONE.
    if (ctx->con_index.len32 < n_repeats+1) {
        buf_alloc_255 (vb, &ctx->con_index, 0, n_repeats+1, WordIndex, 0, CTX_TAG_CON_INDEX);
        ctx->con_index.len32 = ctx->con_index.len32 / sizeof (WordIndex);
    }

    // count printable item separators
    unsigned num_printable_separators=0;
    for (unsigned i=0; i < con.nitems_lo; i++) 
        num_printable_separators += is_printable[con.items[i].separator[0]] + is_printable[con.items[i].separator[1]];

    WordIndex node_index = *B(WordIndex, ctx->con_index, n_repeats);
    unsigned account_for = con.repeats * num_printable_separators /* item seperators */ + con.repeats - con.drop_final_repsep /* repeat separators */ + ((int)add_bytes - snip_len);

    // case: first container with this many repeats - seg and add to cache
    if (node_index == WORD_INDEX_NONE) 
        *B(WordIndex, ctx->con_index, n_repeats) = container_seg (vb, ctx, (ContainerP)&con, NULL, 0, account_for);
    
    // case: we already know the node index of the container with this many repeats
    else 
        seg_known_node_index (vb, ctx, node_index, account_for);

    return n_repeats;

badly_formatted:
    // roll back all the changed data
    seg_rollback (vb);

    // now just seg the entire snip
    seg_by_ctx (VB, STRa(snip), ctx, add_bytes); 

    return -1; // not segged as a container
}                           

void seg_array_of_array_of_struct (VBlockP vb, ContextP ctx, 
                                   char outer_sep,
                                   MediumContainer inner_con, // container of array of struct
                                   STRp(snip), const SegCallback *callbacks) // optional - either NULL, or contains a seg callback for each item (any callback may be NULL)
{
    str_split (snip, snip_len, 0, outer_sep, inner, false);

    bytes d = ctx->dict_id.id;
    DictId inner_dict_id = { .id = { d[0], d[1], '-', d[2], d[3], d[4], d[5], d[6] } };

    MiniContainer outer_con = { .drop_final_repsep = true,  
                                .repeats           = n_inners, 
                                .nitems_lo         = 1,
                                .repsep            = { outer_sep },
                                .items             = { { .dict_id = inner_dict_id } }
                              };

    container_seg (vb, ctx, (ContainerP)&outer_con, 0, 0, n_inners-1); // account for outer_sep

    ContextP inner_ctx = ctx_get_ctx (vb, inner_dict_id);

    ctx_consolidate_stats (vb, ctx->did_i, inner_ctx->did_i, DID_EOL);

    for (int i=0; i < n_inners; i++)
        seg_array_of_struct (vb, inner_ctx, inner_con, STRi(inner,i), callbacks, NULL, inner_lens[i]);
}                                   

// returns true if successful
bool seg_by_container (VBlockP vb, ContextP ctx, ContainerP con, STRp(value), 
                       STRp(container_snip), // optional
                       SegCallback item_seg, // optional
                       bool normal_seg_if_fail,
                       unsigned add_bytes)
{
    ASSERT (con->repeats == 1, "repeats=%u, but currently only supports repeats=1", con->repeats);

    str_split_by_container (value, value_len, con, 0, 0, item, NULL);

    // case: value doesn't match the container
    if (!n_items) {
        if (normal_seg_if_fail) {
            seg_by_ctx (vb, STRa(value), ctx, add_bytes);
            return true;  // segged
        }
        else
            return false; // not segged
    }

    int accounted_for = 0;
    for (int i=0; i < n_items; i++) {
        ContextP item_ctx = ctx_get_ctx (vb, con->items[i].dict_id);
        
        if (item_seg)
            item_seg (vb, item_ctx, STRi(item, i), 0);
        
        else if (IS_LT_DYN (item_ctx->ltype))
            seg_integer_or_not (vb, item_ctx, STRi(item, i), item_lens[i]);
        
        else
            seg_by_ctx (vb, STRi(item, i), item_ctx, item_lens[i]);
        
        accounted_for += item_lens[i];
    }

    if (container_snip)
        seg_by_ctx (vb, STRa(container_snip), ctx, 0);
    else
        container_seg (vb, ctx, con, 0, 0, 0);

    ctx->txt_len = (int)ctx->txt_len + (int)add_bytes - accounted_for; // careful signed arithmetic as addition might be negative

    return true;
}


void seg_add_to_local_fixed_do (VBlockP vb, ContextP ctx, const void *const data, uint32_t data_len, bool add_nul, Lookup lookup_type, bool is_singleton, unsigned add_bytes) 
{
#ifdef DEBUG
    ASSERT (is_singleton || ctx->ltype != LT_SINGLETON || segconf.running, 
            "%s: ctx %s requires ltype!=LT_SINGLETON", LN_NAME, ctx->tag_name);
#endif

    buf_alloc (vb, &ctx->local, data_len + add_nul, 100000, char, CTX_GROWTH, CTX_TAG_LOCAL);
    if (data_len) buf_add (&ctx->local, data, data_len); 
    if (add_nul) BNXTc (ctx->local) = 0;

    if (add_bytes) ctx->txt_len += add_bytes;
    ctx->local_num_words++;

    if (lookup_type == LOOKUP_SIMPLE) 
        seg_simple_lookup (vb, ctx, 0);

    else if (lookup_type == LOOKUP_WITH_LENGTH) 
        seg_lookup_with_length (vb, ctx, data_len, 0);
}


void seg_integer_fixed (VBlockP vb, ContextP ctx, void *number, bool with_lookup, unsigned add_bytes) 
{
    buf_alloc (vb, &ctx->local, 0, MAX_(ctx->local.len+1, 32768) * lt_width(ctx), char, CTX_GROWTH, CTX_TAG_LOCAL);

    ASSERT (ctx->ltype >= LT_INT8 && ctx->ltype <= LT_FLOAT64, "%s: %s.ltype=%s can be segged as fixed", 
            VB_NAME, ctx->tag_name, lt_name(ctx->ltype));

    switch (lt_width(ctx)) {
        case 1  : BNXT8  (ctx->local) = *(uint8_t  *)number; break;
        case 2  : BNXT16 (ctx->local) = *(uint16_t *)number; break;
        case 4  : BNXT32 (ctx->local) = *(uint32_t *)number; break;
        case 8  : BNXT64 (ctx->local) = *(uint64_t *)number; break;
        default : ABORT ("Unexpected ltype=%s", lt_name(ctx->ltype));
    }

    if (add_bytes) ctx->txt_len += add_bytes;
    ctx->local_num_words++;

    if (with_lookup)     
        seg_by_ctx (vb, (char[]){ SNIP_LOOKUP }, 1, ctx, 0);
}

void seg_diff (VBlockP vb, ContextP ctx, 
               ContextP base_ctx,  // optional for xor against other
               STRp(value), 
               bool entire_snip_if_same, 
               unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (ctx->ltype == LT_UINT8, "expecting %s to have ltype=LT_UINT8 but found %s", ctx->tag_name, lt_name(ctx->ltype));

    // entire_snip_if_same is for self-diffing - we might be better off just segging the snip again without a diff than
    // an empty diff (lower b250 entropy)
    ASSERT (ctx == base_ctx || !entire_snip_if_same, "unexpected entire_snip_if_same in %s", ctx->tag_name);
#endif

    if (base_ctx->last_txt.len < value_len || !value_len) goto fallback;
    
    bytes last = (bytes)last_txtx (vb, base_ctx); // caller's responsibility to make sure there's a correct value in here, consistent with flags.same_line
    
    bool exact = !memcmp (value, last, value_len); 
    if (exact) { 
        if (entire_snip_if_same) goto fallback;  // seg the full snip - good for some fields in collated files - we're better off having a series of identical b250s than a intermitted "copy" snip
        else                     goto skip_diff; // seg a "copy" snip with no diff in local
    }

    buf_alloc_zero (vb, &ctx->local, value_len, 100 * value_len, uint8_t, CTX_GROWTH, CTX_TAG_LOCAL);

    uint8_t *diff = BAFT8 (ctx->local); 

    // diff against the beginning of last (last is permitted to be longer that value)
    for (uint32_t i=0; i < value_len; i++) 
        // up to v13:  diff[i] = last[i] ^ ((uint8_t *)value)[i];
        if (last[i] != ((uint8_t *)value)[i]) 
            diff[i] = ((uint8_t *)value)[i];

    ctx->local.len32 += value_len;

skip_diff:
    if (ctx == base_ctx) {
        SNIPi1 (SNIP_DIFF, exact ? -(int32_t)value_len : (int32_t)value_len); // negative marks "exact"
        seg_by_ctx (vb, STRa(snip), ctx, add_bytes);
    }
    else {
        STRli(snip, 48);
        seg_prepare_snip_other (SNIP_DIFF, base_ctx->dict_id, true, exact ? -(int32_t)value_len : (int32_t)value_len, snip);
        seg_by_ctx (vb, STRa(snip), ctx, add_bytes);
    }

    return;
    
fallback:
    seg_by_ctx (vb, STRa(value), ctx, add_bytes);
}

static void seg_set_hash_hints (VBlockP vb, int third_num)
{
    if (third_num == 1) vb->num_lines_at_1_3 = vb->line_i + 1;
    else                vb->num_lines_at_2_3 = vb->line_i + 1;

    for (Did did_i=0; did_i < vb->num_contexts; did_i++) {

        decl_ctx (did_i);
        if (!ctx->nodes.len32 || ctx->global_hash.len32) continue; // our service is not needed - global_hash for this dict already exists or no nodes

        if (third_num == 1) ctx->nodes_len_at_1_3 = ctx->nodes.len32;
        else                ctx->nodes_len_at_2_3 = ctx->nodes.len32;
    }
}

// double the number of lines if we've run out of lines
static void seg_more_lines (VBlockP vb, unsigned sizeof_line)
{
    // note: sadly, we cannot use the normal Buffer macros here because each data_type has its own line type
    buf_alloc_zero (vb, &vb->lines, 0, (vb->lines.len + 1) * sizeof_line, char, 2, "lines");    
    vb->lines.len = vb->lines.size / sizeof_line;
}

static void seg_verify_file_size (VBlockP vb)
{
    uint32_t recon_size = 0; // reconstructed size, as viewed in reconstruction

    // sanity checks
    ASSERT (vb->recon_size >= 0, "%s: recon_size=%d is negative", VB_NAME, vb->recon_size);

    for (Did sf_i=0; sf_i < vb->num_contexts; sf_i++) 
        recon_size += CTX(sf_i)->txt_len;
    
    if ((vb->recon_size != recon_size || flag.debug_recon_size) && !flag.show_bam) { 

        fprintf (stderr, "context.txt_len for vb=%s:\n", VB_NAME);
        for_ctx_that (ctx->nodes.len || ctx->local.len || ctx->txt_len)
            fprintf (stderr, "%s: %u\n", ctx_tag_name_ex (ctx).s, (uint32_t)ctx->txt_len);

        fprintf (stderr, "%s: reconstructed_vb_size=%s (calculated by adding up ctx.txt_len after segging) but vb->recon_size=%s (initialized when reading the file and adjusted for modifications) (diff=%d) (vblock_memory=%s)\n",
                 VB_NAME, str_int_commas (recon_size).s, str_int_commas (vb->recon_size).s, 
                 (int32_t)recon_size - (int32_t)vb->recon_size, str_size (segconf.vb_size).s);

        ASSERT (vb->recon_size == recon_size, "Error while verifying reconstructed size.\n"
                "To get vblock: genozip --biopsy %u %s", vb->vblock_i, txt_name);
    }
}

// increment progress indicator
static void seg_increment_progress (VBlockP vb, uint32_t bytes_so_far_this_vb, int64_t *prev_increment)
{
    // by bytes or by lines
    int64_t this_increment = (txt_file->est_num_lines ? vb->line_i : bytes_so_far_this_vb);
    
    int64_t increment = this_increment - *prev_increment;
    
    if (increment > (txt_file->est_num_lines ? 2500/*lines*/ : 1000000/*bytes*/)) {
        dispatcher_increment_progress ("seg", increment);
        *prev_increment = this_increment;
    }
}

// split each lines in this VB to its components
void seg_all_data_lines (VBlockP vb)
{
    START_TIMER;

    // sanity (leave 64b to detect bugs)
    ASSERT (vb->lines.len <= vb->txt_data.len, "%s: Expecting lines.len=%"PRIu64" < txt_data.len=%"PRIu64, 
            VB_NAME, vb->lines.len, vb->txt_data.len); // 64 bit test in case of memory corruption

    // note: empty VB is possible, for example empty SAM generated component
    // note: if re-reading, data is not loaded yet (it will be in *_seg_initialize)
    ASSERT (!Ltxt || vb->reread_prescription.len || *BLSTtxt == '\n' || !DTP(vb_end_nl), "%s: %s txt_data unexpectedly doesn't end with a newline. Last 10 chars: \"%10s\"", 
            VB_NAME, dt_name(vb->data_type), Btxt (Ltxt - MIN_(10,Ltxt)));
    
    // set estimated number of lines
    vb->lines.len32 = vb->lines.len32  ? vb->lines.len32 // already set? don't change (eg 2nd pair FASTQ, bcl_unconsumed)
                    : segconf.running  ? 10              // low number of avoid memory overallocation for PacBio arrays etc 
                    : segconf.line_len ? MAX_(1, Ltxt / segconf.line_len)
                    :                    1;              // eg DT_GENERIC

    vb->scratch.name = "scratch"; // initialize so we don't need to worry about it later
    
    ContextP debug_lines_ctx = NULL;
    if (flag.debug_lines && !segconf.running) { 
        debug_lines_ctx = ECTX (_SAM_DEBUG_LINES); // same dict_id for all data_types (ie _SAM_DEBUG_LINES==_VCF_DEBUG_LINES etc)
        if (debug_lines_ctx) debug_lines_ctx->ltype = LT_UINT32;
    }
    
    DT_FUNC (vb, seg_initialize)(vb);  // data-type specific initialization

    // if local is going to be compressed using a callback, we can't have singletons go to local
    // note: called after seg_initialize, to allow for setting of ctx->no_callback where needed 
    zip_set_no_stons_if_callback (vb);
    
    // allocate lines
    uint32_t sizeof_line = DT_FUNC_OPTIONAL (vb, sizeof_zip_dataline, 1)(); // 1 - we waste a little bit of memory to avoid making exceptions throughout the code logic
    buf_alloc_zero (vb, &vb->lines, 0, vb->lines.len * sizeof_line, char, 1, "lines");

    rom field_start = B1STtxt;
    bool hash_hints_set_1_3 = false, hash_hints_set_2_3 = false;
    int64_t progress = 0;
    int64_t n_lines_processed=0; // number of lines sent to segging. some of them might have been sent to gencomp. 

    for (vb->line_i=0; vb->line_i < vb->lines.len32; vb->line_i++, n_lines_processed++) {

        if (!segconf.running) 
            seg_increment_progress (vb, BNUMtxt(field_start), &progress);
        
        uint32_t remaining_txt_len = BREMAINS (vb->txt_data, field_start);
        
        if (!remaining_txt_len) { // we're done
            vb->lines.len32 = vb->line_i; // update to actual number of lines
            break;
        }

        bool has_13 = false;
        vb->line_start = BNUMtxt (field_start);

        if (flag.show_lines)
            iprintf ("%s byte-in-vb=%u\n", LN_NAME, vb->line_start);

        // Call the segmenter of the data type to segment one line
        rom next_field = DTP(seg_txt_line) (vb, field_start, remaining_txt_len, &has_13);

        uint32_t line_len = next_field - field_start;

        if (debug_lines_ctx) {
            if (!vb->debug_line_hash_skip) {
                uint32_t hash = DTP(seg_modifies) ? vb->debug_line_hash // if seg_modifies, Seg calculates hash before modifications
                                                  : adler32 (2, field_start, line_len); // same as in container_verify_line_integrity
                seg_integer_fixed (vb, debug_lines_ctx, &hash, false, 0);
            }
            else 
                vb->debug_line_hash_skip = false; // reset
        }

        if (flag.biopsy_line.line_i == vb->line_i && flag.biopsy_line.vb_i == vb->vblock_i && !DTP(seg_modifies)) {
            file_put_line (vb, field_start, line_len, "Line biopsy:");
            exit_ok;
        }

        if (line_len > vb->longest_line_len) vb->longest_line_len = line_len;

        field_start = next_field;

        // update line_bgzf_uoffset to next line
        if (txt_file->codec == CODEC_BGZF && vb->comp_i == COMP_MAIN) 
            bgzf_zip_advance_index (vb, line_len);
        
        // if our estimate number of lines was too small, increase it
        if (vb->line_i == vb->lines.len32-1 && field_start - vb->txt_data.data != vb->txt_data.len)         
            seg_more_lines (vb, sizeof_line);
        
        // collect stats at the approximate 1/3 or 2/3s marks of the file, to help hash_alloc_global create a hash
        // table. note: we do this for every vb, not just 1, because hash_alloc_global runs in the first
        // vb a new field/subfield is introduced
        if (!hash_hints_set_1_3 && BNUMtxt (field_start) > Ltxt / 3) {
            seg_set_hash_hints (vb, 1);
            hash_hints_set_1_3 = true;
        }
        else if (!hash_hints_set_2_3 && BNUMtxt (field_start) > 2 * Ltxt / 3) {
            seg_set_hash_hints (vb, 2);
            hash_hints_set_2_3 = true;
        }

        // --head advanced option in ZIP cuts a certain number of first lines from vb=1
        if (vb->vblock_i==1 && (flag.lines_last != NO_LINE) && flag.lines_last == vb->line_i) {
            vb->recon_size -= BREMAINS (vb->txt_data, field_start);
            vb->lines.len32 = vb->line_i + 1;
            break;
        }
    }

    ASSINP (vb->lines.len32 <= CON_MAX_REPEATS, // because top_level.repeats = vb->lines.len
            "Genozip works by dividing the file to \"VBlocks\". Unfortuantely, the VBlocks for this file are too big\n"
            "and have have too many %ss (= over the Genozip's maximum of %u). Current VBlock size is %s.\n"
            "Solution: use --vblock to set a lower value (value is in MB)",
            DTP(line_name), CON_MAX_REPEATS, str_size (segconf.vb_size).s);

    DT_FUNC (vb, seg_finalize)(vb); // data-type specific finalization

     if (!flag.make_reference && !segconf.running && !flag.biopsy && flag.biopsy_line.line_i == NO_LINE) 
        seg_verify_file_size (vb);

    // txt_size and lines.len exclude lines moved to gencomp. increment might be negative
    dispatcher_increment_progress ("seg_final", 
                                   txt_file->est_num_lines ? (n_lines_processed - (int64_t)vb->lines.len) 
                                                           : ((int64_t)vb->txt_size - progress));

    if (flag.debug_or_test) buflist_test_overflows(vb, __FUNCTION__); 

    COPY_TIMER (seg_all_data_lines);
}
