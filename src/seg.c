// ------------------------------------------------------------------
//   seg.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
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
#include "libdeflate/libdeflate.h"

WordIndex seg_by_ctx_ex (VBlockP vb, STRp(snip), ContextP ctx, uint32_t add_bytes,
                         bool *is_new) // optional out
{
    ASSERTNOTNULL (ctx);

    WordIndex node_index = ctx_create_node_do (VB, ctx, STRa(snip), is_new);

    ASSERT (node_index < ctx->nodes.len32 + ctx->ol_nodes.len32 || node_index == WORD_INDEX_EMPTY || node_index == WORD_INDEX_MISSING, 
            "out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u",  
            ctx->tag_name, node_index, ctx->nodes.len32, ctx->ol_nodes.len32);
    
    ctx_append_b250 (vb, ctx, node_index);
    ctx->txt_len += add_bytes;
    
    return node_index;
} 

// segs the same node as previous seg
WordIndex seg_known_node_index (VBlockP vb, ContextP ctx, WordIndex node_index, unsigned add_bytes) 
{ 
    ctx_append_b250 (vb, ctx, node_index);
    ctx->txt_len += add_bytes;
    
    (*B32 (ctx->counts, node_index))++;

    return node_index;
}

WordIndex seg_duplicate_last(VBlockP vb, ContextP ctx, unsigned add_bytes) 
{ 
    ASSERTISALLOCED (ctx->b250);
    return seg_known_node_index (vb, ctx, LASTb250(ctx), add_bytes); 
}

#define MAX_ROLLBACK_CTXS ARRAY_LEN(vb->rollback_ctxs)

// NOTE: this does not save recon_size - if the processing may change recon_size it must be saved and rolled back separately
void seg_create_rollback_point (VBlockP vb, 
                                ContextP *ctxs,    // option 1
                                unsigned num_ctxs, // applied to both options
                                ...)               // option 2
{
    ASSERT (num_ctxs <= MAX_ROLLBACK_CTXS, "num_ctxs=%u > MAX_ROLLBACK_CTXS=%u", num_ctxs, MAX_ROLLBACK_CTXS);

    va_list args;
    va_start (args, num_ctxs);

    vb->rback_id++; // new rollback point
    vb->num_rollback_ctxs = 0;

    for (unsigned i=0; i < num_ctxs; i++) {
        vb->rollback_ctxs[vb->num_rollback_ctxs] = ctxs ? ctxs[i] : CTX(va_arg (args, int)); // note: 'Did' {aka 'short unsigned int'} is promoted to 'int' when passed through '...'

        if (ctx_set_rollback (vb, vb->rollback_ctxs[vb->num_rollback_ctxs], false)) // added
            vb->num_rollback_ctxs++;
    }
}

// adds a ctx to the current rollback point, if its not already there
void seg_add_ctx_to_rollback_point (VBlockP vb, ContextP ctx)
{
    ASSERT (vb->num_rollback_ctxs < MAX_ROLLBACK_CTXS, "num_ctxs=%u > MAX_ROLLBACK_CTXS=%u", vb->num_rollback_ctxs+1, MAX_ROLLBACK_CTXS);

    vb->rollback_ctxs[vb->num_rollback_ctxs] = ctx;
    
    if (ctx_set_rollback (vb, ctx, false)) // added
        vb->num_rollback_ctxs++;
}

void seg_rollback (VBlockP vb)
{
    for (uint32_t i=0; i < vb->num_rollback_ctxs; i++) 
        ctx_rollback (vb, vb->rollback_ctxs[i], false);
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

rom seg_get_next_item (void *vb_, rom str, int *str_len, 
                               GetNextAllow newline, GetNextAllow tab, GetNextAllow space,
                               unsigned *len, char *separator, 
                               bool *has_13, // out - only modified if '\r' detected ; only needed if newline=GN_SEP
                               rom item_name)
{
    VBlockP vb = (VBlockP)vb_;

    unsigned i=0; for (; i < *str_len; i++) {
        char c = str[i];
        if ((tab     == GN_SEP && c == '\t') ||
            (newline == GN_SEP && c == '\n') ||
            (space   == GN_SEP && c == ' ')) {
                *len = i;
                *separator = c;
                *str_len -= i+1;

                // check for Windows-style '\r\n' end of line 
                if (i && c == '\n' && str[i-1] == '\r') {
                    (*len)--;
                    ASSERT0 (has_13, "has_13=NULL but expecting it because newline=GN_SEP");
                    *has_13 = true;
                }

                return str + i+1; // beyond the separator
        }
        else if ((tab     == GN_FORBIDEN && c == '\t') ||
                 (newline == GN_FORBIDEN && c == '\n') ||
                 (space   == GN_FORBIDEN && c == ' ' ) ) break;
    }

    ASSSEG (*str_len, str, "missing %s field", item_name);

    ABOSEG (str, "while segmenting %s: expecting a %s%s%safter \"%.*s\"", 
            item_name,
            newline==GN_SEP ? "NEWLINE " : "", tab==GN_SEP ? "TAB " : "", space==GN_SEP ? "\" \" " : "", 
            MIN_(i, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

// returns first character after current line
rom seg_get_next_line (void *vb_, rom str, 
                               int *remaining_len, // in/out
                               unsigned *len,      // out - length of line exluding \r and \n
                               bool must_have_newline, bool *has_13 /* out */, rom item_name)
{
    VBlockP vb = (VBlockP)vb_;

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

                return str + *len + *has_13 + 1; // beyond the separator
        }
    
    ASSSEG (*remaining_len, str, "missing %s field", item_name);

    // if we reached here - line doesn't end with a newline
    ASSSEG (!must_have_newline, str, "while segmenting %s: expecting a NEWLINE after (showing at most 1000 characters): \"%.*s\"", 
            item_name, MIN_(*remaining_len, 1000), str);

    // we have no newline, but check if last character is a \r
    *has_13 = (after > str) && (after[-1] == '\r');
    *len = *remaining_len - *has_13;
    *remaining_len = 0;

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

WordIndex seg_integer_as_text_do (VBlockP vb, ContextP ctx, int64_t n, unsigned add_bytes)
{
    char snip[24];
    unsigned snip_len = str_int (n, snip);
    return seg_by_ctx (vb, STRa(snip), ctx, add_bytes);
}

// prepare snips that contain code + dict_id + optional parameter (SNIP_OTHER_LOOKUP, SNIP_OTHER_DELTA, SNIP_REDIRECTION...)
void seg_prepare_snip_other_do (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int64_t int_param, char char_param, 
                                char *snip, unsigned *snip_len) //  in / out
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
                                             char *out_snip, unsigned *out_snip_len) // in/out - allocated by caller
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
void seg_prepare_minus_snip_do (DictId dict_id_a, DictId dict_id_b, uint8_t special_code, char *snip, unsigned *snip_len)
{
    snip[0] = SNIP_SPECIAL;
    snip[1] = special_code;
    
    DictId two_dicts[2] = { dict_id_a, dict_id_b };
    *snip_len = 2 + base64_encode ((uint8_t *)two_dicts, sizeof (two_dicts), &snip[2]);
}

void seg_mux_display (MultiplexerP mux)
{
    iprintf ("ctx=%s num_channels=%u snip_len=%u snip=%.100s", 
             mux->ctx->tag_name, mux->num_channels, MUX_SNIP_LEN(mux), MUX_SNIP(mux));

    for (unsigned i=0; i < MIN_(mux->num_channels, 1024); i++)
        iprintf ("[%u] dict_id=%s channel_ctx=%p\n", i, dis_dict_id (mux->dict_ids[i]).s, MUX_CHANNEL_CTX(mux)[i]);
}

void seg_mux_init (VBlockP vb, ContextP ctx, unsigned num_channels, uint8_t special_code, 
                   bool no_stons, MultiplexerP mux, rom channel_letters) // optional - a string with num_channels unique characters
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
        mux->dict_ids[i].id[1] = channel_letters ? channel_letters[i] : (((uint16_t)id[1] + (uint16_t)i) & 0xff);
        if (id_len <= 5) str_int (i, (char*)&mux->dict_ids[i].id[id_len]); // add numerical digits for ease of debugging if we have room
    }

    // prepare snip - a string consisting of num_channels x { special_code or \t , base64(dict_id.id) }
    MUX_SNIP_LEN(mux) = BASE64_DICT_ID_LEN * (unsigned)num_channels;   
    seg_prepare_multi_dict_id_special_snip (special_code, num_channels, mux->dict_ids, MUX_SNIP(mux), &MUX_SNIP_LEN(mux));

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
        channel_ctx->flags.store  = mux->ctx->flags.store;
        channel_ctx->ltype        = mux->ctx->ltype; 
        channel_ctx->no_stons     = mux->no_stons; // not inherited from parent
        channel_ctx->st_did_i     = (mux->ctx->st_did_i == DID_NONE) ? mux->ctx->did_i : mux->ctx->st_did_i;
        channel_ctx->flags.ctx_specific_flag = mux->ctx->flags.ctx_specific_flag; // inherit - needed for "same_line"        
    }

    return channel_ctx;
}

// scans a pos field - in case of non-digit or not in the range [0,MAX_POS], either returns -1
// (if allow_nonsense) or errors
static PosType seg_scan_pos_snip (VBlockP vb, rom snip, unsigned snip_len, bool zero_is_bad, 
                                  SegError *err) // out - if NULL, exception if error
{
    PosType value=0;
    bool is_int = str_get_int (snip, snip_len, &value); // value unchanged if not integer

    if (is_int && value >= !!zero_is_bad && value <= MAX_POS) {
        *err = ERR_SEG_NO_ERROR;
        return value; // all good
    }

    ASSSEG (err, snip, "position field must be an integer number between %u and %"PRId64", seeing: %.*s", 
            !!zero_is_bad, MAX_POS, STRf(snip));

    *err = is_int ? ERR_SEG_OUT_OF_RANGE : ERR_SEG_NOT_INTEGER;
    return value; // in- or out-of- range integer, or 0 if value is not integer
}

// returns POS value if a valid pos, or 0 if not
PosType seg_pos_field (VBlockP vb, 
                       Did snip_did_i,    // mandatory: the ctx the snip belongs to
                       Did base_did_i,    // mandatory: base for delta
                       unsigned opt,      // a combination of SPF_* options
                       char missing,      // a character allowed (meaning "missing value"), segged as SNIP_DONT_STORE
                       STRp(pos_str),     // option 1
                       PosType this_pos,  // option 2
                       unsigned add_bytes)     
{
    ContextP snip_ctx = CTX(snip_did_i);
    ContextP base_ctx = CTX(base_did_i);

    // note: *_seg_initialize is required to set both snip_ctx and base_ctx: no_stons=true, store=STORE_INT

    SegError err = ERR_SEG_NO_ERROR;
    if (pos_str) { // option 1
        if (pos_str_len == 1 && *pos_str == missing) 
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

    PosType pos_delta = this_pos - base_ctx->last_value.i;
    
    ctx_set_last_value (vb, snip_ctx, this_pos);

    // if the delta is too big, add this_pos (not delta) to local and put SNIP_LOOKUP in the b250
    // EXCEPT if it is the first vb (ie last_pos==0) because we want to avoid creating a whole RANDOM_POS
    // section in every VB just for a single entry in case of a nicely sorted file
    if ((!IS_FLAG (opt, SPF_UNLIMITED_DELTA) && (ABS(pos_delta) > MAX_POS_DELTA) && base_ctx->last_value.i) ||
        IS_FLAG (opt, SPF_NO_DELTA)) {

        snip_ctx->ltype = LT_UINT32;
        snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
        seg_integer_fixed (vb, snip_ctx, &this_pos, true, add_bytes);

        return this_pos;
    }

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

        unsigned delta_len = str_int (pos_delta, &pos_delta_str[total_len]);
        total_len += delta_len;

        seg_by_ctx (VB, pos_delta_str, total_len, snip_ctx, add_bytes);
    }
    // case: the delta is the negative of the previous delta - add a SNIP_SELF_DELTA with no payload - meaning negated delta
    else {
        char negated_delta = SNIP_SELF_DELTA; // no delta means negate the previous delta
        seg_by_did (VB, &negated_delta, 1, snip_did_i, add_bytes);
        snip_ctx->last_delta = 0; // no negated delta next time
    }

    return this_pos;
}

bool seg_pos_field_cb (VBlockP vb, ContextP ctx, STRp(pos_str), uint32_t repeat)
{
    seg_pos_field (vb, ctx->did_i, ctx->did_i, 0, 0, STRa(pos_str), 0, pos_str_len);
    return true; // segged successfully
}

// must be called in seg_initialize 
void seg_id_field_init (ContextP ctx) 
{ 
    ctx->no_stons = true; 
    ctx->ltype = LT_DYN_INT; 
    ctx->is_initialized = true; 
}

// Commonly (but not always), IDs are SNPid identifiers like "rs17030902". We store the ID divided to 2:
// - We store the final digits, if any exist, and up to 9 digits, as an integer in SEC_NUMERIC_ID_DATA
// - In the dictionary we store the prefix up to this number, and \1 if there is a number and a \2 
//   if the caller (as in seg_me23) wants us to store an extra bit.
// example: rs17030902 : in the dictionary we store "rs\1" or "rs\1\2" and in SEC_NUMERIC_ID_DATA we store 17030902.
//          1423       : in the dictionary we store "\1" and 1423 SEC_NUMERIC_ID_DATA
//          abcd       : in the dictionary we store "abcd" and nothing is stored SEC_NUMERIC_ID_DATA

// this version compresses better if the numeric part is expected to be variable-width without leading zeros
void seg_id_field_do (VBlockP vb, ContextP ctx, STRp(id))
{
    ASSERT (ctx->is_initialized, "%s: seg_id_field_init not called for ctx=%s", LN_NAME, ctx->tag_name);

    int i=id_len-1; for (; i >= 0; i--) 
        if (!IS_DIGIT (id[i])) break;
    
    unsigned num_digits = MIN_(id_len - (i+1), 18); // up to 0x0DE0,B6B3,A763,FFFF - fits in int64_t

    // leading zeros will be part of the dictionary data, not the number
    for (unsigned i = id_len - num_digits; i < id_len; i++) 
        if (id[i] == '0') 
            num_digits--;
        else 
            break;

    // added to local if we have a trailing number
    if (num_digits) {
        int64_t id_num;
        ASSERT (str_get_int (&id[id_len - num_digits], num_digits, &id_num), 
                "Failed str_get_int ctx=%s vb=%u line=%d", ctx->tag_name, vb->vblock_i, vb->line_i);
        seg_add_to_local_resizable (vb, ctx, id_num, 0);

        if (ctx->flags.store == STORE_INT) ctx_set_last_value (vb, ctx, id_num);
    }

    // prefix the textual part with SNIP_LOOKUP_UINT32 if needed (we temporarily overwrite the previous separator or the buffer underflow area)
    unsigned new_len = id_len - num_digits;
    SAFE_ASSIGN (&id[-1], SNIP_LOOKUP); // we assign it anyway bc of the macro convenience, but we included it only if num_digits>0
    seg_by_ctx (VB, id-(num_digits > 0), new_len + (num_digits > 0), ctx, id_len); // account for the entire length, and sometimes with \t
    SAFE_RESTORE;
}

// this version compresses better if the numeric part is expected to be fixed-width with possible leading zeros
void seg_id_field2 (VBlockP vb, ContextP ctx, STRp(id), unsigned add_bytes)
{
    ASSERT (ctx->is_initialized, "%s: seg_id_field_init not called for ctx=%s", LN_NAME, ctx->tag_name);

    int i=id_len-1; for (; i >= 0; i--) 
        if (!IS_DIGIT (id[i])) break;
    
    unsigned num_digits = MIN_(id_len - (i+1), 18); // up to 0x0DE0,B6B3,A763,FFFF - fits in int64_t (excess digits will go in the prefix)

    if (num_digits) {
        int64_t n;
        str_get_int_dec (&id[id_len - num_digits], num_digits, (uint64_t*)&n);

        seg_integer (vb, ctx, n, false, add_bytes); // integer into local

        SAFE_ASSIGNx(&id[-3], SNIP_NUMERIC,      1);
        SAFE_ASSIGNx(&id[-2], '0'/*LT_DYN_INT*/, 2);
        SAFE_ASSIGNx(&id[-1], '0' + num_digits,  3);
        
        seg_by_ctx (vb, id-3, 3 + (id_len - num_digits), ctx, 0); // SNIP_NUMERIC with prefix

        SAFE_RESTOREx(1); SAFE_RESTOREx(2); SAFE_RESTOREx(3);

        if (ctx->flags.store == STORE_INT) ctx_set_last_value (vb, ctx, n);
    }

    else
        seg_by_ctx (vb, STRa(id), ctx, add_bytes);
}

bool seg_id_field_cb (VBlockP vb, ContextP ctx, STRp(id), uint32_t repeat)
{
    seg_id_field_do (vb, ctx, STRa(id));
    return true; // segged successfully
}

bool seg_id_field2_cb (VBlockP vb, ContextP ctx, STRp(id), uint32_t repeat)
{
    seg_id_field2 (vb, ctx, STRa(id), id_len);
    return true; // segged successfully
}

// note: caller must set ctx->ltype=LT_DYN_INT*
void seg_integer (VBlockP vb, ContextP ctx, int64_t n, bool with_lookup, unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || ctx->ltype == LT_DYN_INT || ctx->ltype == LT_DYN_INT_h || ctx->ltype == LT_DYN_INT_H,
            "ctx=%s must have a LT_DYN_INT* ltype", ctx->tag_name);

ASSERT (segconf.running || !(VB_DT(SAM) || VB_DT(BAM)) || !dict_id_is_aux_sf(ctx->dict_id) || ctx->flags.store == STORE_INT,
            "ctx=%s must have a STORE_INT ltype", ctx->tag_name); // needed to allow recon to translate to BAM
#endif

    ctx_set_last_value (vb, ctx, (ValueType){.i = n});

    seg_add_to_local_resizable (vb, ctx, n, add_bytes);

    if (with_lookup)     
        seg_by_ctx (vb, (char[]){ SNIP_LOOKUP }, 1, ctx, 0);
}

// returns true if it was an integer
bool seg_integer_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned add_bytes)
{
    int64_t n;
    bool is_int = (ctx->ltype == LT_DYN_INT_h) ? (value[0] != '0' && str_get_int_hex (STRa(value), true, false, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                : (ctx->ltype == LT_DYN_INT_H) ? (value[0] != '0' && str_get_int_hex (STRa(value), false, true, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                :                                str_get_int (STRa(value), &n);

    // case: its an integer
    if (!ctx->no_stons && is_int) { // we interpret no_stons as means also no moving ints to local (one of the reasons is that an int might actually be a float)
        if (ctx->ltype < LT_DYN_INT) ctx->ltype = LT_DYN_INT; // note: the LT_DYN* types are the last in LocalType
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
void seg_numeric_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned field_width, unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || ctx->ltype >= LT_DYN_INT, "Expecting %s.ltype to be DYN_INT*", ctx->tag_name);
#endif

    if (segconf.running) goto fallback;

    int64_t n;
    bool is_numeric = (ctx->ltype == LT_DYN_INT_h) ? (str_get_int_hex (STRa(value), true, false, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                    : (ctx->ltype == LT_DYN_INT_H) ? (str_get_int_hex (STRa(value), false, true, (uint64_t*)&n)) // number with leading 0 is segged as a snip
                    :                                 str_get_int_dec (STRa(value), (uint64_t*)&n);

    // case: its an integer
    if (is_numeric) { // we interpret no_stons as means also no moving ints to local (one of the reasons is that an int might actually be a float)
        seg_integer (vb, ctx, n, false, add_bytes);
        seg_by_ctx (vb, (char[]){ SNIP_NUMERIC, '0'+ (ctx->ltype - LT_DYN_INT), '0' + field_width }, 3, ctx, 0);
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

// if its a float, stores the float in local, and a LOOKUP in b250, and returns true. if not - normal seg, and returns false.
bool seg_float_or_not (VBlockP vb, ContextP ctx, STRp(value), unsigned add_bytes)
{
    // TO DO: implement reconstruction in reconstruct_one_snip-SNIP_LOOKUP
    char snip[1 + FLOAT_FORMAT_LEN];
    unsigned format_len;

    // case: its an float
    if (!ctx->no_stons && // we interpret no_stons as means also no moving ints to local (one of the reasons is that an int might actually be a float)
        str_get_float (STRa(value), &ctx->last_value.f, &snip[1], &format_len)) {

        // verify that we can reconstruct the number precisely (not always the case, 
        // as the textual float may exceed float precision)
        float f = (float)ctx->last_value.f; // 64bit -> 32bit
        char recon[value_len+1];
        sprintf (recon, &snip[1], f);
        if (memcmp (value, recon, value_len)) goto fallback;

        ctx->ltype = LT_FLOAT32; // set only upon storing the first number - if there are no numbers, leave it as LT_TEXT so it can be used for singletons

        // add to local
        buf_alloc (vb, &ctx->local, 1, vb->lines.len, float, CTX_GROWTH, "contexts->local");
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

WordIndex seg_delta_vs_other_do (VBlockP vb, ContextP ctx, ContextP other_ctx, 
                                 STRp(value), // if value==NULL, we use value_n
                                 int64_t value_n,

                                 int64_t max_delta /* max abs value of delta - beyond that, seg as is, ignored if < 0 */,
                                 unsigned add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || (ctx->flags.store == STORE_INT && other_ctx->flags.store == STORE_INT), 
            "expecting ctx=%s and other_ctx=%s to have STORE_INT", ctx->tag_name, other_ctx->tag_name);
#endif

    ASSERTNOTNULL (other_ctx);

    if (value && !str_get_int (STRa(value), &value_n)) goto fallback;

    ctx_set_last_value (vb, ctx, value_n);

    int64_t delta = value_n - other_ctx->last_value.i; 
    if (max_delta >= 0 && ABS(delta) > max_delta) goto fallback;

    SNIP(32);
    seg_prepare_snip_other (SNIP_OTHER_DELTA, other_ctx->dict_id, true, (int32_t)delta, snip);

    return seg_by_ctx (VB, STRa(snip), ctx, add_bytes);

fallback:
    return value ? seg_by_ctx (VB, STRa(value), ctx, add_bytes) 
                 : seg_integer_as_text_do (vb, ctx, value_n, add_bytes);
}

// note: seg_initialize should set STORE_INT for this ctx
WordIndex seg_self_delta (VBlockP vb, ContextP ctx, int64_t value, char format/* 0 for normal decimal*/, uint32_t add_bytes)
{
#ifdef DEBUG
    ASSERT (segconf.running || ctx->flags.store == STORE_INT, "expecting %s to have store=STORE_INT", ctx->tag_name);
#endif
    char delta_snip[30];
    delta_snip[0] = SNIP_SELF_DELTA;
    if (format) delta_snip[1] = format; // currently only 'x' (lower case hex) is supported
    unsigned delta_snip_len = 1 + !!format + str_int (value - ctx->last_value.i, &delta_snip[1 + !!format]);

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
    MiniContainer *con;
    ContextP arr_ctx;
    int additional_bytes = add_bytes - value_len;

    // first use in this VB - prepare context where array elements will go in
    if (!container_ctx->con_cache.len32) {
        bytes id = container_ctx->dict_id.id;

        if (!arr_dict_id.num)        
            arr_dict_id = (DictId){ .id = { id[0], 
                                            ((id[1]+1) % 256) | 128, // different IDs for top array, subarray and items
                                            id[2], id[3], id[4], id[5], id[6], id[7] } };
        
        buf_alloc (vb, &container_ctx->con_cache, 0, 1, MiniContainer, 1, "contexts->con_cache");

        con = B1ST (MiniContainer, container_ctx->con_cache);
        *con = (MiniContainer){ .nitems_lo = 1, 
                                .drop_final_repsep = true,
                                .repsep    = { sep },
                                .items     = { { .dict_id = arr_dict_id } } }; // only one item

        arr_ctx = ctx_get_ctx (vb, arr_dict_id);
        arr_ctx->st_did_i = container_ctx->st_did_i = stats_conslidation_did_i;
    }
    else { 
        con         = B1ST (MiniContainer, container_ctx->con_cache);
        arr_dict_id = con->items[0].dict_id;
        arr_ctx     = ctx_get_ctx (vb, arr_dict_id);
    }

    if (use_integer_delta || container_ctx->flags.store == STORE_INT) 
        arr_ctx->flags.store = STORE_INT;

    // count repeats (1 + number of seperators)
    con->repeats=1;
    for (int32_t i=0; i < value_len; i++) 
        if (value[i] == sep) con->repeats++;

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
            seg_self_delta (vb, arr_ctx, this_item_value, 0, this_item_len);

        // non-integer that cannot be delta'd - store as-is
        else 
            seg_by_ctx (VB, STRa(this_item), arr_ctx, 0);

        if (container_ctx->flags.store == STORE_INT)
            container_ctx->last_value.i += this_item_value;

        value_len--; // skip separator
        value++;
    }

    return container_seg (vb, container_ctx, (ContainerP)con, 0, 0, con->repeats-1 + additional_bytes); // acount for separators and additional bytes
}

bool seg_do_nothing_cb (VBlockP vb, ContextP ctx, STRp(field), uint32_t rep)
{
    return true; // "segged successfully" - do nothing
}

// a field that looks like: "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
// we have an array (2 entires in this example) of items (4 in this examples) - the entries are separated by comma and the items by space
// observed in Ensembel generated GVF: Variant_effect, sift_prediction, polyphen_prediction, variant_peptide
// The last item is treated as an ENST_ID (format: ENST00000399012) while the other items are regular dictionaries
// the names of the dictionaries are the same as the ctx, with the 2nd character replaced by 1,2,3...
// the field itself will contain the number of entries
int32_t seg_array_of_struct (VBlockP vb, ContextP ctx, MediumContainer con, STRp(snip), 
                             const SegCallback *callbacks, // optional - either NULL, or contains a seg callback for each item (any callback may be NULL)
                             unsigned add_bytes)
{
    ContextP ctxs[con.nitems_lo]; 

    for (unsigned i=0; i < con.nitems_lo; i++) 
        ctxs[i] = ctx_get_ctx (vb, con.items[i].dict_id); 

    if (!ctx->is_stats_parent && ctx->st_did_i == DID_NONE) 
        ctx_consolidate_stats_(vb, ctx->did_i, con.nitems_lo, ctxs);

    seg_create_rollback_point (vb, ctxs, con.nitems_lo);

    // get repeats
    str_split (snip, snip_len, 0, con.repsep[0], repeat, false);

    // if we don't have con.drop_final_repsep, the last "repeat" should be zero length and removed
    if (!con.drop_final_repsep) {
        if (repeat_lens[n_repeats-1]) goto badly_formatted; 
        n_repeats--;
    }
    con.repeats = n_repeats;

    ASSSEG (n_repeats <= CONTAINER_MAX_REPEATS, snip, "exceeded maximum repeats allowed (%u) while parsing %s",
            CONTAINER_MAX_REPEATS, ctx->tag_name);

    for (unsigned r=0; r < n_repeats; r++) {

        // get items in each repeat
        str_split_by_container (repeats[r], repeat_lens[r], &con, NULL, 0, item);
        
        if (n_items != con.nitems_lo) 
            goto badly_formatted;

        for (unsigned i=0; i < n_items; i++)
            if (callbacks && callbacks[i]) {
                if (!callbacks[i] (vb, ctxs[i], STRi(item,i), r))
                    goto badly_formatted;
            }
            else
                seg_by_ctx (VB, STRi(item,i), ctxs[i], item_lens[i]);
    }

    // finally, the Container snip itself - we attempt to use the known node_index if it is cached in con_index

    // in our specific case, element i of con_index contains the node_index of the snip of the container with i repeats, or WORD_INDEX_NONE.
    if (ctx->con_index.len32 < n_repeats+1) {
        buf_alloc_255 (vb, &ctx->con_index, 0, n_repeats+1, WordIndex, 0, "con_index");
        ctx->con_index.len32 = ctx->con_index.len32 / sizeof (WordIndex);
    }

    // count printable item separators
    unsigned num_printable_separators=0;
    for (unsigned i=0; i < con.nitems_lo; i++) 
        num_printable_separators += is_printable[con.items[i].separator[0]] + is_printable[con.items[i].separator[1]];

    WordIndex node_index = *B(WordIndex, ctx->con_index, n_repeats);
    unsigned account_for = con.repeats * num_printable_separators /* item seperators */ + con.repeats - con.drop_final_repsep /* repeat separators */ + ((int)add_bytes - snip_len);

    // case: first container with the many repeats - seg and add to cache
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

// returns true if successful
bool seg_by_container (VBlockP vb, ContextP ctx, ContainerP con, STRp(value), 
                       STRp(container_snip), // optional
                       bool normal_seg_if_fail,
                       unsigned add_bytes)
{
    ASSERT (con->repeats == 1, "repeats=%u, but currently only supports repeats=1", con->repeats);

    str_split_by_container (value, value_len, con, 0, 0, item);

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

        if (item_ctx->ltype >= LT_DYN_INT)
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


void seg_add_to_local_fixed_do (VBlockP vb, ContextP ctx, STRp(data), bool add_nul, Lookup lookup_type, bool is_singleton, unsigned add_bytes) 
{
#ifdef DEBUG
    ASSERT (is_singleton || ctx->no_stons || ctx->ltype != LT_TEXT, "ctx %s requires no_stons or should have an ltype other than LT_TEXT", ctx->tag_name);
#endif

    buf_alloc (vb, &ctx->local, data_len + add_nul, 100000, char, CTX_GROWTH, "contexts->local");
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
    buf_alloc (vb, &ctx->local, 0, MAX_(ctx->local.len+1, 32768) * lt_desc[ctx->ltype].width, char, CTX_GROWTH, "contexts->local");

    switch (lt_desc[ctx->ltype].width) {
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
    ASSERT (ctx->ltype == LT_UINT8, "expecting %s to have ltype=LT_UINT8", ctx->tag_name);

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

    buf_alloc_zero (vb, &ctx->local, value_len, 100 * value_len, uint8_t, CTX_GROWTH, "contexts->local");

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
        STRl(snip, 48) = 48;
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

        ContextP ctx = CTX(did_i);
        if (!ctx->nodes.len || ctx->global_hash_prime) continue; // our service is not needed - global_cache for this dict already exists or no nodes

        if (third_num == 1) ctx->nodes_len_at_1_3 = ctx->nodes.len;
        else                ctx->nodes_len_at_2_3 = ctx->nodes.len;
    }
}

// double the number of lines if we've run out of lines
static void seg_more_lines (VBlockP vb, unsigned sizeof_line)
{
    // note: sadly, we cannot use the normal Buffer macros here because each data_type has its own line type
    buf_alloc_zero (vb, &vb->lines, 0, (vb->lines.len + 1) * sizeof_line, char, 2, "lines");    
    vb->lines.len = vb->lines.size / sizeof_line;

    // allocate more to the b250 buffer of the fields
    for (Did did_i=0; did_i < DTF(num_fields); did_i++) 
        if (segconf.b250_per_line[did_i])
            buf_alloc (vb, &CTX(did_i)->b250, 0, AT_LEAST(did_i), WordIndex, 1, "contexts->b250");
}

static void seg_verify_file_size (VBlockP vb)
{
    uint32_t recon_size = 0; // reconstructed size, as viewed in reconstruction

    // sanity checks
    ASSERT (vb->recon_size >= 0, "%s: recon_size=%d is negative", VB_NAME, vb->recon_size);

    for (Did sf_i=0; sf_i < vb->num_contexts; sf_i++) 
        recon_size += CTX(sf_i)->txt_len;
    
    uint32_t vb_recon_size = DTP(seg_get_vb_recon_size) ? DTP(seg_get_vb_recon_size)(vb) : vb->recon_size; // normally vb->recon_size, but callback can override

    if ((vb_recon_size != recon_size || flag.debug_recon_size) && !flag.show_bam) { 

        fprintf (stderr, "context.txt_len for vb=%s:\n", VB_NAME);
        for (Did sf_i=0; sf_i < vb->num_contexts; sf_i++) {
            ContextP ctx = CTX(sf_i);
            if (ctx->nodes.len || ctx->local.len || ctx->txt_len)
                fprintf (stderr, "%s: %u\n", ctx_tag_name_ex (ctx).s, (uint32_t)ctx->txt_len);
        }

        fprintf (stderr, "%s: reconstructed_vb_size=%s (calculated by adding up ctx.txt_len after segging) but vb->recon_size%s=%s (initialized when reading the file and adjusted for modifications) (diff=%d) (vblock_memory=%s)\n",
                 VB_NAME, str_int_commas (recon_size).s, (vb->data_type == DT_VCF && vcf_vb_is_luft(vb)) ? "_luft" : "", str_int_commas (vb_recon_size).s, 
                 (int32_t)recon_size - (int32_t)vb_recon_size, str_size (segconf.vb_size).s);

        ASSERT (vb_recon_size == recon_size, "Error while verifying reconstructed size.\n"
                "To get vblock: genozip --biopsy %u %s", vb->vblock_i, txt_name);
    }
}

// split each lines in this VB to its components
void seg_all_data_lines (VBlockP vb)
{
    START_TIMER;

    // sanity
    ASSERT (vb->lines.len <= vb->txt_data.len, "%s: Expecting lines.len=%"PRIu64" < txt_data.len=%"PRIu64, 
            VB_NAME, vb->lines.len, vb->txt_data.len);

    // note: empty VB is possible, for example empty SAM generated component
    // note: if re-reading, data is not loaded yet (it will be in *_seg_initialize)
    ASSERT (!vb->txt_data.len || vb->reread_prescription.len || *BLSTtxt == '\n' || !DTP(vb_end_nl), "%s: %s txt_data unexpectedly doesn't end with a newline. Last 10 chars: \"%10s\"", 
            VB_NAME, dt_name(vb->data_type), Bc (vb->txt_data, vb->txt_data.len - MIN_(10,vb->txt_data.len)));

    ctx_initialize_predefined_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts); // Create ctx for the fields in the correct order 
 
    // allocate the b250 for the fields which each have num_lines entries
    for (Did did_i=0; did_i < DTF(num_fields); did_i++)
        if (segconf.b250_per_line[did_i]) 
            buf_alloc (vb, &CTX(did_i)->b250, 0, AT_LEAST(did_i), WordIndex, 1, "contexts->b250");
    
    // set estimated number of lines
    vb->lines.len = vb->lines.len    ? vb->lines.len // already set? don't change (eg 2nd pair FASTQ, bcl_unconsumed)
                  : segconf.running  ? 10 // low number of avoid memory overallocation for PacBio arrays etc 
                  : segconf.line_len ? MAX_(1, vb->txt_data.len / segconf.line_len)
                  :                    1;            // eg DT_GENERIC

    vb->scratch.name = "scratch"; // initialize so we don't need to worry about it later
    
    ContextP debug_lines_ctx = NULL;
    if (flag.debug_lines && !segconf.running) { 
        debug_lines_ctx = ECTX (_SAM_DEBUG_LINES); // same dict_id for all data_types (ie _SAM_DEBUG_LINES==_VCF_DEBUG_LINES etc)
        if (debug_lines_ctx) debug_lines_ctx->ltype = LT_UINT32;
    }
    
    DT_FUNC (vb, seg_initialize)(vb);  // data-type specific initialization

    // allocate lines
    uint32_t sizeof_line = DT_FUNC_OPTIONAL (vb, sizeof_zip_dataline, 1)(); // 1 - we waste a little bit of memory to avoid making exceptions throughout the code logic
    buf_alloc_zero (vb, &vb->lines, 0, vb->lines.len * sizeof_line, char, 1, "lines");

    rom field_start = B1STtxt;
    bool hash_hints_set_1_3 = false, hash_hints_set_2_3 = false;
    int64_t prev_increment = 0;

    for (vb->line_i=0; vb->line_i < vb->lines.len; vb->line_i++) {

        // increment progress indicator
        int64_t increment = BNUMtxt(field_start) - prev_increment;
        if (increment > 1000000 && !segconf.running) {
            dispatcher_increment_progress ("seg", increment);
            prev_increment = BNUMtxt(field_start);
        }

        uint32_t remaining_txt_len = BREMAINS (vb->txt_data, field_start);
        
        if (!remaining_txt_len) { // we're done
            vb->lines.len = vb->line_i; // update to actual number of lines
            break;
        }

        bool has_13 = false;
        vb->line_start = BNUMtxt (field_start);

        if (flag.show_lines)
            iprintf ("%s byte-in-vb=%u\n", LN_NAME, vb->line_start);

        // Call the segmenter of the data type to segment one line
        rom next_field = DT_FUNC (vb, seg_txt_line) (vb, field_start, remaining_txt_len, &has_13);
        if (!next_field) next_field = field_start + remaining_txt_len; // DT_GENERIC has no segmenter

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
            exit_ok();
        }

        if (line_len > vb->longest_line_len) vb->longest_line_len = line_len;

        field_start = next_field;

        // update line_bgzf_uoffset to next line
        if (txt_file->codec == CODEC_BGZF && vb->comp_i == COMP_MAIN) 
            bgzf_zip_advance_index (vb, line_len);
        
        // if our estimate number of lines was too small, increase it
        if (vb->line_i == vb->lines.len-1 && field_start - vb->txt_data.data != vb->txt_data.len)         
            seg_more_lines (vb, sizeof_line);
        
        // collect stats at the approximate 1/3 or 2/3s marks of the file, to help hash_alloc_global create a hash
        // table. note: we do this for every vb, not just 1, because hash_alloc_global runs in the first
        // vb a new field/subfield is introduced
        if (!hash_hints_set_1_3 && (field_start - vb->txt_data.data) > vb->txt_data.len / 3) {
            seg_set_hash_hints (vb, 1);
            hash_hints_set_1_3 = true;
        }
        else if (!hash_hints_set_2_3 && (field_start - vb->txt_data.data) > 2 * vb->txt_data.len / 3) {
            seg_set_hash_hints (vb, 2);
            hash_hints_set_2_3 = true;
        }

        // --head advanced option in ZIP cuts a certain number of first lines from vb=1
        if (vb->vblock_i==1 && (flag.lines_last != NO_LINE) && flag.lines_last == vb->line_i) {
            vb->recon_size -= BREMAINS (vb->txt_data, field_start);
            vb->lines.len = vb->line_i + 1;
            break;
        }
    }

    if (segconf.running) {
        segconf.line_len = (vb->lines.len ? (vb->txt_data.len / vb->lines.len) : 500) + 1; // get average line length (rounded up ; arbitrary 500 if the segconf data ended up not having any lines (example: all lines were non-matching lines dropped by --match in a chain file))

        // limitations: only pre-defined field, not local
        for (Did did_i=0; did_i < DTF(num_fields); did_i++)
            if (CTX(did_i)->b250.len) 
                segconf.b250_per_line[did_i] = (float)CTX(did_i)->b250.len / (float)vb->lines.len;
    }

    DT_FUNC (vb, seg_finalize)(vb); // data-type specific finalization

    if (!flag.make_reference && !segconf.running && flag.biopsy_line.line_i == NO_LINE) 
        seg_verify_file_size (vb);

    dispatcher_increment_progress ("seg_final", (int64_t)vb->txt_size - prev_increment); // txt_size excludes lines moved to gencomp. increment might be negative

    if (flag.debug) buf_test_overflows(vb, __FUNCTION__); 

    COPY_TIMER (seg_all_data_lines);
}