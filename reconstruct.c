// ------------------------------------------------------------------
//   reconstruct.c
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "reconstruct.h"
#include "vblock.h"
#include "context.h"
#include "file.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h"
#include "container.h"
#include "flags.h"
#include "piz.h"
#include "base64.h"
#include "regions.h"
#include "lookback.h"

// Compute threads: decode the delta-encoded value of the POS field, and returns the new lacon_pos
// Special values:
// "-" - negated previous value
// ""  - negated previous delta
static int64_t reconstruct_from_delta (VBlockP vb, 
                                       ContextP my_ctx,   // use and store last_delta
                                       ContextP base_ctx, // get last_value
                                       STRp(delta_snip),
                                       bool reconstruct) 
{
    ASSPIZ0 (delta_snip, "delta_snip is NULL");
    ASSPIZ (base_ctx->flags.store == STORE_INT, "reconstructing %s - calculating delta \"%.*s\" from a base of %s, but %s, doesn't have STORE_INT",
            my_ctx->tag_name, delta_snip_len, delta_snip, base_ctx->tag_name, base_ctx->tag_name);

    int64_t base_value = (my_ctx->flags.delta_peek && my_ctx != base_ctx)
        ? reconstruct_peek (vb, base_ctx, 0, 0).i // value of this line/sample - whether already encountered or peek a future value
        : base_ctx->last_value.i; // use last_value even if base_ctx not encountered yet

    bool hex = delta_snip_len && delta_snip[0] == 'x';
    if (hex) {
        delta_snip_len--;
        delta_snip++;
    }

    if (delta_snip_len == 1 && delta_snip[0] == '-')
        my_ctx->last_delta = -2 * base_value; // negated previous value

    else if (!delta_snip_len)
        my_ctx->last_delta = -my_ctx->last_delta; // negated previous delta

    else 
        my_ctx->last_delta = (int64_t)strtoull (delta_snip, NULL, 10 /* base 10 */); // strtoull can handle negative numbers, despite its name

    int64_t new_value = base_value + my_ctx->last_delta;  
    if (reconstruct && !hex)
        RECONSTRUCT_INT (new_value);
    else if (reconstruct && hex)
        RECONSTRUCT_HEX (new_value, false);

    return new_value;
}

#define ASSERT_IN_BOUNDS \
    ASSPIZ (ctx->next_local < ctx->local.len, \
            "unexpected end of ctx->local data in %s (len=%u next_local=%u ltype=%s lcodec=%s did_i=%u)", \
            ctx->tag_name, (uint32_t)ctx->local.len, ctx->next_local, lt_name (ctx->ltype), codec_name (ctx->lcodec), ctx->did_i)

#define ASSERT_IN_BOUNDS_BEFORE(recon_len) \
    ASSPIZ (ctx->next_local + (recon_len) <= ctx->local.len, \
            "unexpected end of ctx->local data in %s (len=%u next_local=%u ltype=%s lcodec=%s did_i=%u)", \
            ctx->tag_name, (uint32_t)ctx->local.len, ctx->next_local, lt_name (ctx->ltype), codec_name (ctx->lcodec), ctx->did_i)

static uint32_t reconstruct_from_local_text (VBlockP vb, ContextP ctx, bool reconstruct)
{
    uint32_t start = ctx->next_local; 
    ARRAY (char, data, ctx->local);

    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != 0) ctx->next_local++;
    ASSERT_IN_BOUNDS;

    char *snip = &data[start];
    uint32_t snip_len = ctx->next_local - start; 
    ctx->next_local++; /* skip the separator */ 

    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, STRa(snip), reconstruct);

    return snip_len;
}

static void reconstruct_from_xor_diff (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct)
{
    int64_t str_len;
    ASSPIZ (str_get_int (STRa(snip), &str_len), "In ctx=%s: Invalid XOR_DIFF snip: \"%.*s", ctx->tag_name, snip_len, snip);

    ASSERT_IN_BOUNDS_BEFORE(str_len);

    bytes last = (uint8_t *)last_txtx (vb, ctx);
    bytes this = B8 (ctx->local, ctx->next_local); 
    uint8_t *recon = BAFT8 (vb->txt_data);

    for (int64_t i=0; i < str_len; i++)
        recon[i] = last[i] ^ this[i];

    vb->txt_data.len += str_len;
    ctx->next_local  += str_len;    
}

int64_t reconstruct_from_local_int (VBlockP vb, ContextP ctx, char separator /* 0 if none */, bool reconstruct)
{
    ASSERT_IN_BOUNDS;

    int64_t num=0;

    switch (ctx->ltype) {
        case LT_UINT8:  num = NEXTLOCAL(uint8_t,  ctx); break;
        case LT_UINT32: num = NEXTLOCAL(uint32_t, ctx); break;
        case LT_INT8:   num = NEXTLOCAL(int8_t,   ctx); break;
        case LT_INT32:  num = NEXTLOCAL(int32_t,  ctx); break;
        case LT_UINT16: num = NEXTLOCAL(uint16_t, ctx); break;
        case LT_INT16:  num = NEXTLOCAL(int16_t,  ctx); break;
        case LT_UINT64: num = NEXTLOCAL(uint64_t, ctx); break;
        case LT_INT64:  num = NEXTLOCAL(int64_t,  ctx); break;
        default: 
            ASSPIZ (false, "Unexpected ltype=%s(%u) for ctx=\"%s\"", lt_name(ctx->ltype), ctx->ltype, ctx->tag_name); 
    }

    if (reconstruct) { 
        if (VB_DT(DT_VCF) && num==lt_desc[ctx->ltype].max_int && dict_id_is_vcf_format_sf (ctx->dict_id)
            && !lt_desc[ctx->ltype].is_signed) {
            RECONSTRUCT1 ('.');
            num = 0; // we consider FORMAT fields that are . to be 0.
        }
        else 
            RECONSTRUCT_INT (num);
        
        if (separator) RECONSTRUCT1 (separator);
    }

    return num;
}

int64_t reconstruct_peek_local_int (VBlockP vb, ContextP ctx, int offset /*0=next_local*/)
{
    ASSERT_IN_BOUNDS;

    int64_t num=0;

    switch (ctx->ltype) {
        case LT_UINT8:  num = PEEKNEXTLOCAL(uint8_t,  ctx, offset); break;
        case LT_UINT32: num = PEEKNEXTLOCAL(uint32_t, ctx, offset); break;
        case LT_INT8:   num = PEEKNEXTLOCAL(int8_t,   ctx, offset); break;
        case LT_INT32:  num = PEEKNEXTLOCAL(int32_t,  ctx, offset); break;
        case LT_UINT16: num = PEEKNEXTLOCAL(uint16_t, ctx, offset); break;
        case LT_INT16:  num = PEEKNEXTLOCAL(int16_t,  ctx, offset); break;
        case LT_UINT64: num = PEEKNEXTLOCAL(uint64_t, ctx, offset); break;
        case LT_INT64:  num = PEEKNEXTLOCAL(int64_t,  ctx, offset); break;
        default: 
            ASSPIZ (false, "Unexpected ltype=%s(%u)", lt_name(ctx->ltype), ctx->ltype); 
    }

    if (VB_DT(DT_VCF) && num==lt_desc[ctx->ltype].max_int && dict_id_is_vcf_format_sf (ctx->dict_id)
        && !lt_desc[ctx->ltype].is_signed)
        return 0; // returns 0 if '.'

    return num;
}

static double reconstruct_from_local_float (VBlockP vb, ContextP ctx, 
                                            STRp(format), // required unless not reconstructing (eg a binary field - CI0_TRANS_NOR in the container)
                                            char separator /* 0 if none */, bool reconstruct)
{   
    ASSERT_IN_BOUNDS;

    float num=0;

    switch (ctx->ltype) {
        case LT_FLOAT32: num = NEXTLOCAL(float,  ctx); break;
        case LT_FLOAT64: num = NEXTLOCAL(double, ctx); break;
        default: 
            ASSPIZ (false, "Unexpected ltype=%s(%u) in ctx=%s", lt_name(ctx->ltype), ctx->ltype, ctx->tag_name); 
    }

    if (reconstruct) { 
        ASSPIZ (format_len, "Failed to reconstruct a float in ctx=%s, because format was not given", ctx->tag_name);

        // format is as generated by str_get_float - %8.3f. each of the two numbers can be 1 or 2 digits.
        SAFE_NULT (format);
        sprintf (BAFTc (vb->txt_data), format, num);      
        vb->txt_data.len += (format[2] == '.' ? (format[1]-'0') : ((format[1]-'0') * 10 + format[2]-'0'));
        SAFE_RESTORE;

        if (separator) RECONSTRUCT1 (separator);
    }

    return num;
}

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from vb->seq_len.
// NOTE: this serves nucleotide sequences AND qual. Bad design. Should have been two separate things.
void reconstruct_from_local_sequence (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct)
{
    ASSERTNOTNULL (ctx);

    // bool reconstruct = !piz_is_skip_section (vb, SEC_LOCAL, ctx->dict_id);
    if (!ctx->is_loaded) return;
    uint32_t len;

    // if we have length in the snip, update vb->seq_len (for example in FASTQ, we will a snip for seq but qual will use seq_len)
    if (snip_len) vb->seq_len = atoi(snip);

    // case: handle SAM missing quality (may be expressed as a ' ' or ASCII 127)
    char c = ctx->local.data[ctx->next_local];
    if (c == ' ' || c == 127) {
        len = 1;
        sam_reconstruct_missing_quality (vb, c, reconstruct);
    }

    else {
        len = vb->seq_len;
        ASSPIZ (ctx->next_local + len <= ctx->local.len, "unexpected end of %s data: expecting ctx->next_local=%u + seq_len=%u <= local.len=%u", 
                ctx->tag_name, ctx->next_local, len, ctx->local.len32);

        if (reconstruct) RECONSTRUCT (Bc(ctx->local, ctx->next_local), len);
    }

    ctx->last_value.i = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
    ctx->next_local += len;
}

ContextP reconstruct_get_other_ctx_from_snip (VBlockP vb, ContextP ctx, pSTRp (snip))
{
    unsigned b64_len = base64_sizeof (DictId);
    char err[*snip_len+20];
    ASSPIZ (b64_len + 1 <= *snip_len, "ctx=%s snip=\"%s\" snip_len=%u but expecting it to be >= %u", 
            ctx->tag_name, str_print_snip(*snip, *snip_len, err), *snip_len, b64_len + 1);

    DictId dict_id;
    base64_decode ((*snip)+1, &b64_len, dict_id.id);

    ContextP other_ctx = ECTX (dict_id);
    ASSPIZ (other_ctx, "Failed to get other context: ctx=%s snip=%.*s other_dict_id=%s", 
            ctx->tag_name, STRf(*snip), dis_dict_id(dict_id).s);
  
    *snip     += b64_len + 1;
    *snip_len -= b64_len + 1;

    ctx->other_did_i = other_ctx->did_i;

    return other_ctx;
}

// get ctx from a multi-dict_id special snip. note that we're careful to only ECTX the ctx_i requested, and not all,
// so that we don't do a full search of vb->contexts[] for a channel that was not segged and hence has no context
ContextP recon_multi_dict_id_get_ctx_first_time (VBlockP vb, ContextP ctx, STRp(snip), unsigned ctx_i)
{
    if (!ctx->con_cache.len) {
        ctx->con_cache.len = str_count_char (STRa(snip), '\t') + 1;
        buf_alloc_zero (vb, &ctx->con_cache, 0, ctx->con_cache.len, ContextP, 1, "con_cache");
    }

    // note: we get past this point only once per VB, per ctx_i
    str_split (snip, snip_len, ctx->con_cache.len, '\t', item, true);
    ASSPIZ (n_items, "Unable to decoded multi-dict-id snip for %s. snip=\"%.*s\"", ctx->tag_name, snip_len, snip);

    DictId item_dict_id;
    base64_decode (items[ctx_i], &item_lens[ctx_i], item_dict_id.id);

    return (*B(ContextP, ctx->con_cache, ctx_i) = ECTX (item_dict_id)); // NULL if no data was segged to this channel    
}

void reconstruct_set_buddy (VBlockP vb)
{
    ASSPIZ0 (vb->buddy_line_i == NO_BUDDY, "Buddy line already set for the current line");

    ContextP buddy_ctx = ECTX (_SAM_BUDDY); // all data types using buddy are expected to have a dict_id identical to _SAM_BUDDY (but different did_i)

    int32_t num_lines_back = reconstruct_from_local_int (vb, buddy_ctx, 0, false);

    // a bug that existed 12.0.41-13.0.1 (see bug 367): we stored buddy in machine endianty instead of BGEN32.
    // When we reach this point pizzing in a buggy file, num_lines_back will be in BGEN since piz_adjust_one_local set it. 
    // We detect buggy files by local.param=0 (see sam_seg_initialize) and convert it back to machine endianity.
    if (!buddy_ctx->local.param)    
        num_lines_back = BGEN32 ((uint32_t)num_lines_back);

    vb->buddy_line_i = vb->line_i - num_lines_back; // convert value passed (distance in lines to buddy) to 0-based buddy_line_i

    ASSPIZ (vb->buddy_line_i >= 0, "Expecting vb->buddy_line_i=%d to be non-negative. num_lines_back=%d buddy_ctx->local.param=%d", vb->buddy_line_i, num_lines_back, (int)buddy_ctx->local.param);
}

void reconstruct_from_buddy_get_textual_snip (VBlockP vb, ContextP ctx, pSTRp(snip))
{
    ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line");

    HistoryWord word = *B(HistoryWord, ctx->history, vb->buddy_line_i);
    Buffer *buf=NULL;
    switch (word.lookup) {
        case LookupTxtData : buf = &vb->txt_data  ; break;
        case LookupDict    : buf = &ctx->dict     ; break;
        case LookupLocal   : buf = &ctx->local    ; break;
        case LookupPerLine : buf = &ctx->per_line ; break;
        default : ASSPIZ (false, "Invalid value word.lookup=%d", word.lookup);
    }

    ASSPIZ (word.char_index < buf->len, "buddy word ctx=%s buddy_line_i=%d char_index=%"PRIu64" is out of range of buffer %s len=%"PRIu64, 
            ctx->tag_name, vb->buddy_line_i, word.char_index, buf->name, buf->len);

    *snip = Bc (*buf, word.char_index);
    *snip_len = word.snip_len;
}

// Copy from buddy: buddy is data that appears on a specific "buddy line", in this context or another one. Not all lines need
// Note of difference vs. lookback: with buddy, not all lines need to have the data (eg MC:Z), so the line number is constant,
// but if we had have used lookback, the lookback value would have been different between different fields.  
bool reconstruct_from_buddy (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct, ValueType *new_value)
{
    ContextP base_ctx = ctx;

    // set buddy if needed, and not already set (in BAM it is already set in sam_piz_filter).
    if (snip_len==1 && *snip == SNIP_COPY_BUDDY && vb->buddy_line_i == NO_BUDDY) 
        reconstruct_set_buddy (vb); 

    // optional: base context is different than ctx
    else if (snip_len > 1) {
        snip--; snip_len++; // reconstruct_get_other_ctx_from_snip skips the first char
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, ctx, pSTRa(snip));
    }

    // case: numeric value 
    if (ctx->flags.store == STORE_INT) {
        ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line"); // for textual, we test in reconstruct_from_buddy_get_textual_snip
        ASSPIZ (base_ctx->history.len, "history not set for %s, perhaps seg forgot to set store_per_line?", base_ctx->tag_name);

        new_value->i = *B(int64_t, base_ctx->history, vb->buddy_line_i);
        if (reconstruct) RECONSTRUCT_INT (new_value->i);
        return true; // has new value
    }

    // case: textual value
    else {
        if (reconstruct) {
            reconstruct_from_buddy_get_textual_snip (vb, base_ctx, pSTRa(snip));
            RECONSTRUCT (snip, snip_len);
        }

        return false; // no new value 
    }
}

static ValueType reconstruct_from_lookback (VBlockP vb, ContextP ctx, STRp(snip), bool reconstruct)
{   
    ContextP lb_ctx = SCTX(snip);
    int64_t lookback = lb_ctx->last_value.i;
    ValueType value = {};

    // a lookback by word_index
    if (!snip_len) { 
        value.i = lookback_get_index (vb, lb_ctx, ctx, lookback);
        
        STR(back_snip);
        ctx_get_snip_by_word_index (ctx, value.i, back_snip);

        if (reconstruct) RECONSTRUCT (back_snip, back_snip_len);
    }

    // a lookback by txt
    else if (snip_len == 1 && (*snip >= 'T' && *snip <= 'z')) { // maximum supported - 122(=z)-84(=T)+1 = 39
        ValueType back_value = lookback_get_value (vb, lb_ctx, ctx, lookback * (*snip - 'T' + 1));

        if (reconstruct) 
            RECONSTRUCT (Bc (vb->txt_data, back_value.txt.index), back_value.txt.len);
    }
    
    // a lookback by delta vs integer
    else { 
        ValueType back_value = lookback_get_value (vb, lb_ctx, ctx, lookback);

        PosType delta;
        ASSPIZ  (str_get_int (STRa(snip), &delta), "Invalid delta snip \"%.*s\"", STRf(snip));

        value.i = back_value.i + delta;
        
        if (reconstruct) RECONSTRUCT_INT (value.i);
    }

    return value; 
}

void reconstruct_one_snip (VBlockP vb, ContextP snip_ctx, 
                           WordIndex word_index, // WORD_INDEX_NONE if not used.
                           STRp(snip), bool reconstruct) // if false, calculates last_value but doesn't output to vb->txt_data)
{
    ValueType new_value = {0};
    bool have_new_value = false;
    int64_t prev_value = snip_ctx->last_value.i;
    ContextP base_ctx = snip_ctx; // this will change if the snip refers us to another data source
    StoreType store_type = snip_ctx->flags.store;
    bool store_delta = z_file->genozip_version >= 12 && snip_ctx->flags.store_delta; // note: the flag was used for something else in v8

    // case: empty snip
    if (!snip_len) {
        if (store_type == STORE_INDEX && word_index != WORD_INDEX_NONE) {
            new_value.i = word_index;
            have_new_value = true;
        }
        goto done;
    }

    switch (snip[0]) {

    // display the rest of the snip first, and then the lookup up text.
    case SNIP_LOOKUP:
    case SNIP_OTHER_LOOKUP: {
        if (snip[0] == SNIP_LOOKUP) 
            { snip++; snip_len--; }
        else 
            // we are request to reconstruct from another ctx
            base_ctx = reconstruct_get_other_ctx_from_snip (vb, snip_ctx, pSTRa(snip)); // also updates snip and snip_len

        switch (base_ctx->ltype) {
            case LT_TEXT:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                reconstruct_from_local_text (vb, base_ctx, reconstruct); // this will call us back recursively with the snip retrieved
                break;
                
            case LT_CODEC:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                codec_args[base_ctx->lcodec].reconstruct (vb, base_ctx->lcodec, base_ctx); break;
                break;

            case LT_INT8 ... LT_UINT64:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                new_value.i = reconstruct_from_local_int (vb, base_ctx, 0, reconstruct);
                have_new_value = true;
                break;

            case LT_FLOAT32 ... LT_FLOAT64:
                new_value.f = reconstruct_from_local_float (vb, base_ctx, STRa(snip), 0, reconstruct);
                have_new_value = true;
                break;
            // case: the snip is taken to be the length of the sequence (or if missing, the length will be taken from vb->seq_len)
            case LT_SEQUENCE: 
                reconstruct_from_local_sequence (vb, base_ctx, STRa(snip), reconstruct);
                break;
                
            case LT_BITMAP:
                ASSERT_DT_FUNC (vb, reconstruct_seq);
                DT_FUNC (vb, reconstruct_seq) (vb, base_ctx, STRa(snip));
                break;

            default: ABORT ("while reconstructing %s in vb_i=%u: Unsupported lt_type=%s (%u) for SNIP_LOOKUP or SNIP_OTHER_LOOKUP. Please upgrade to the latest version of Genozip.", 
                            base_ctx->tag_name, vb->vblock_i, lt_name(base_ctx->ltype), base_ctx->ltype);
        }

        break;
    }

    case SNIP_CONTAINER: {
        STR(prefixes);
        ContainerP con_p = container_retrieve (vb, snip_ctx, word_index, snip+1, snip_len-1, pSTRa(prefixes));
        new_value = container_reconstruct (vb, snip_ctx, con_p, prefixes, prefixes_len); 
        have_new_value = true;
        break;
    }

    case SNIP_SELF_DELTA:
        new_value.i = reconstruct_from_delta (vb, snip_ctx, base_ctx, snip+1, snip_len-1, reconstruct);
        have_new_value = true;
        break;

    case SNIP_MATE_LOOKUP: {
        ASSPIZ (z_file->data_type == DT_FASTQ, "SNIP_MATE_LOOKUP is not expected in ctx=%s since this isn't a FASTQ", snip_ctx->tag_name);

        ASSPIZ (snip_ctx->pair_b250, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2. "
                "If this file was compressed with Genozip version 9.0.12 or older use the --dts_paired command line option. "
                "You can see the Genozip version used to compress it with 'genocat -w %s'", snip_ctx->tag_name, z_name);
                
        ctx_get_next_snip (vb, snip_ctx, snip_ctx->pair_flags.all_the_same, true, &snip, &snip_len);
        reconstruct_one_snip (vb, snip_ctx, WORD_INDEX_NONE /* we can't cache pair items */, snip, snip_len, reconstruct); // might include delta etc - works because in --pair, ALL the snips in a context are PAIR_LOOKUP
        break;
    }
    case SNIP_OTHER_DELTA: 
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, snip_ctx, pSTRa(snip)); // also updates snip and snip_len
        new_value.i = reconstruct_from_delta (vb, snip_ctx, base_ctx, STRa(snip), reconstruct); 
        have_new_value = true;
        break;

    case SNIP_PAIR_DELTA: { // used for FASTQ_GPOS - uint32_t stored in originating in the pair's local
        ASSPIZ (snip_ctx->pair_local, "no pair_1 local data for ctx=%s, while reconstructing pair_2 vb=%u", snip_ctx->tag_name, vb->vblock_i);
        int32_t fastq_line_i = vb->line_i / 4; // see fastq_piz_filter for calculation
        int64_t pair_value = (int64_t) *B32 (snip_ctx->pair, fastq_line_i);  
        int64_t delta = (int64_t)strtoull (snip+1, NULL, 10 /* base 10 */); 
        new_value.i = pair_value + delta;
        if (reconstruct) RECONSTRUCT_INT (new_value.i); 
        have_new_value = true;
        break;
    }

    case SNIP_COPY: 
        base_ctx = (snip_len==1) ? snip_ctx : reconstruct_get_other_ctx_from_snip (vb, snip_ctx, pSTRa(snip)); 
        RECONSTRUCT (last_txtx (vb, base_ctx), base_ctx->last_txt_len);
        new_value = base_ctx->last_value; 
        have_new_value = true;
        break;

    case SNIP_SPECIAL:
        ASSPIZ (snip_len >= 2, "SNIP_SPECIAL expects snip_len=%u >= 2. ctx=\"%s\"", snip_len, snip_ctx->tag_name);
                
        uint8_t special = snip[1] - 32; // +32 was added by SPECIAL macro

        ASSPIZ (special < DTP (num_special), "file requires special handler %u which doesn't exist in this version of genozip - please upgrade to the latest version", special);
        ASSERT_DT_FUNC (vb, special);

        have_new_value = DT_FUNC(vb, special)[special](vb, snip_ctx, snip+2, snip_len-2, &new_value, reconstruct);  
        break;

    case SNIP_XOR_DIFF:
        reconstruct_from_xor_diff (vb, snip_ctx, snip+1, snip_len-1, reconstruct);
        break;

    case SNIP_REDIRECTION: 
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, snip_ctx, pSTRa(snip)); // also updates snip and snip_len
        reconstruct_from_ctx (vb, base_ctx->did_i, 0, reconstruct);
        break;
    
    case SNIP_DUAL: {
        str_split (&snip[1], snip_len-1, 2, SNIP_DUAL, subsnip, true);
        ASSPIZ (n_subsnips==2, "Invalid SNIP_DUAL snip in ctx=%s", snip_ctx->tag_name);

        if (vcf_vb_is_primary(vb)) // recursively call for each side 
            reconstruct_one_snip (vb, snip_ctx, word_index, STRi(subsnip,0), reconstruct);
        else
            reconstruct_one_snip (vb, snip_ctx, word_index, STRi(subsnip,1), reconstruct);
        return;
    }

    case SNIP_LOOKBACK: 
        new_value = reconstruct_from_lookback (vb, base_ctx, STRa(snip), reconstruct);
        have_new_value = true;
        break;

    case SNIP_COPY_BUDDY:
        have_new_value = reconstruct_from_buddy (vb, base_ctx, snip+1, snip_len-1, reconstruct, &new_value);
        break;

    case NUM_SNIP_CODES ... 31:
        ABORT ("File %s requires a SNIP code=%u for dict_id=%s. Please upgrade to the latest version of genozip",
               z_name, snip[0], dis_dict_id (base_ctx->dict_id).s);

    case SNIP_DONT_STORE:
        store_type  = STORE_NONE; // override store and fall through
        store_delta = false;
        snip++; snip_len--;
        // fall through
        
    default: {
        if (reconstruct) RECONSTRUCT (snip, snip_len); // simple reconstruction

        switch (store_type) {
            case STORE_INT: 
                // store the value only if the snip in its entirety is a reconstructable integer (eg NOT "21A", "-0", "012" etc)
                have_new_value = str_get_int (snip, snip_len, &new_value.i);
                break;

            case STORE_FLOAT: {
                char *after;
                new_value.f = strtod (snip, &after); // allows negative values

                // if the snip in its entirety is not a valid number, don't store the value.
                // this can happen for example when seg_pos_field stores a "nonsense" snip.
                have_new_value = (after == snip + snip_len);
                break;
            }
            case STORE_INDEX:
                new_value.i = word_index;
                have_new_value = (word_index != WORD_INDEX_NONE);
                break;

            default: {} // do nothing
        }

        snip_ctx->last_delta = 0; // delta is 0 since we didn't calculate delta
    }
    }

done:
    // update last_value if needed
    if (have_new_value && store_type) // note: we store in our own context, NOT base (a context, eg FORMAT/DP, sometimes serves as a base_ctx of MIN_DP and sometimes as the snip_ctx for INFO_DP)
        ctx_set_last_value (vb, snip_ctx, new_value); // if marely encountered it is set in is set in reconstruct_from_ctx_do

    // note: if store_delta, we do a self-delta. this overrides last_delta set by the delta snip which could be against a different
    // base_ctx. note: when Seg sets last_delta, it must also set store=STORE_INT
    if (store_delta) 
        snip_ctx->last_delta = new_value.i - prev_value;
}

// store last_value in context history - for copying by buddy line
static inline void reconstruct_store_history (VBlockP vb, ContextP ctx, uint32_t last_txt_index)
{
    // case: store last integer value
    if (ctx->flags.store == STORE_INT) 
        *B(int64_t, ctx->history, vb->line_i) = ctx->last_value.i;
    
    // case: textual value will remain in txt_data - just point to it
    else if (!vb->maybe_lines_dropped) 
        *B(HistoryWord, ctx->history, vb->line_i) = 
            (HistoryWord){ .lookup = LookupTxtData, .char_index = last_txt_index, .snip_len = ctx->last_txt_len };
    
    // case: textual value might be removed from txt_data - copy it
    else {
        *B(HistoryWord, ctx->history, vb->line_i) = 
            (HistoryWord){ .lookup = LookupPerLine, .char_index = ctx->per_line.len, .snip_len = ctx->last_txt_len };

        if (ctx->last_txt_len) {
            buf_add_more (vb, &ctx->per_line, Bc (vb->txt_data, last_txt_index), ctx->last_txt_len, "per_line");
            BNXTc (ctx->per_line) = 0; // nul-terminate
        }
    }
}

static void reconstruct_peek_add_ctx_to_frozen_state (VBlockP vb, ContextP ctx); // forward declaration

// returns reconstructed length or -1 if snip is missing and previous separator should be deleted
int32_t reconstruct_from_ctx_do (VBlockP vb, DidIType did_i, 
                                 char sep, // if non-zero, outputs after the reconstruction
                                 bool reconstruct, // if false, calculates last_value but doesn't output to vb->txt_data
                                 rom func)
{
    ASSPIZ (did_i < vb->num_contexts, "called from: %s: did_i=%u out of range: vb->num_contexts=%u for vb_i=%u", 
            func, did_i, vb->num_contexts, vb->vblock_i);

    ContextP ctx = CTX(did_i);

    // if we're peeking, freeze the context
    if (vb->frozen_state.param) 
        reconstruct_peek_add_ctx_to_frozen_state (vb, ctx);

    ASSPIZ0 (ctx->dict_id.num || ctx->did_i != DID_I_NONE, "ctx not initialized (dict_id=0)");

    // update ctx, if its an alias (only for primary field aliases as they have contexts, other alias don't have ctx)
    if (!ctx->dict_id.num) 
        ctx = CTX(ctx->did_i); // ctx->did_i is different than did_i if its an alias

    uint32_t last_txt_index = (uint32_t)vb->txt_data.len;

    // case: we have b250 data
    if (ctx->b250.len ||
        (!ctx->b250.len && !ctx->local.len && ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        STR0(snip);
        WordIndex word_index = LOAD_SNIP(ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (!snip) goto missing;

        reconstruct_one_snip (vb, ctx, word_index, STRa(snip), reconstruct);        

        // if SPECIAL function set value_is_missing (eg vcf_piz_special_PS_by_PID) - this treated as a WORD_INDEX_MISSING 
        if (ctx->value_is_missing) {
            ctx->value_is_missing = false;
            goto missing;
        }

        // for backward compatability with v8-11 that didn't yet have flags.store = STORE_INDEX for CHROM
        if (did_i == DTF(prim_chrom)) { // NOTE: CHROM cannot have aliases, because looking up the did_i by dict_id will lead to CHROM, and this code will be executed for a non-CHROM field
            if (!ctx_has_value_in_line_(vb, CTX(did_i))) 
                vb->last_index (did_i) = word_index;
            
            vb->chrom_node_index = vb->last_index (did_i); 
            if (word_index != vb->chrom_node_index) // eg if word_index is a SPECIAL and the special reconstructor returned the ultimate word index
                ctx_get_snip_by_word_index (CTX(did_i), vb->chrom_node_index, vb->chrom_name);
            else
                STRset (vb->chrom_name, snip);
        }
    }
    
    // case: all data is only in local
    else if (ctx->local.len) {
        switch (ctx->ltype) {
        case LT_INT8 ... LT_UINT64 : {
            int64_t value = reconstruct_from_local_int (vb, ctx, 0, reconstruct); 

            if (ctx->flags.store == STORE_INT) 
                ctx->last_value.i = value;

            break;
        }
        case LT_CODEC:
            codec_args[ctx->lcodec].reconstruct (vb, ctx->lcodec, ctx); break;

        case LT_SEQUENCE: 
            reconstruct_from_local_sequence (vb, ctx, NULL, 0, reconstruct); break;
                
        case LT_BITMAP:
            ASSERT_DT_FUNC (vb, reconstruct_seq);
            DT_FUNC (vb, reconstruct_seq) (vb, ctx, NULL, 0);
            break;
        
        case LT_TEXT:
            reconstruct_from_local_text (vb, ctx, reconstruct); break;

        default:
            ASSPIZ (false, "Invalid ltype=%u in ctx=%s", ctx->ltype, ctx->tag_name);
        }
    }

    // in case of LT_BITMAP, it is it is ok if the bitmap is empty and all the data is in NONREF (e.g. unaligned SAM)
    else if (ctx->ltype == LT_BITMAP && (ctx+1)->local.len) {
        ASSERT_DT_FUNC (vb, reconstruct_seq);
        DT_FUNC (vb, reconstruct_seq) (vb, ctx, NULL, 0);
    }

    // case: the entire VB was just \n - so seg dropped the ctx
    // note: for backward compatability with 8.0. for files compressed by 8.1+, it will be handled via the all_the_same mechanism
    else if (ctx->dict_id.num == DTF(eol).num) {
        if (reconstruct) { RECONSTRUCT1('\n'); }
    }

    else ASSPIZ (flag.missing_contexts_allowed,
                 "ctx %s/%s has no data (dict, b250 or local) in vb_i=%u line_i=%d did_i=%u ctx->did=%u ctx->dict_id=%s ctx->is_loaded=%s", 
                 dtype_name_z (ctx->dict_id), ctx->tag_name, vb->vblock_i, vb->line_i, did_i, ctx->did_i, dis_dict_id (ctx->dict_id).s, TF(ctx->is_loaded));
        
    if (sep && reconstruct) RECONSTRUCT1 (sep); 

    ctx->last_txt_index = last_txt_index;
    ctx->last_txt_len   = (uint32_t)vb->txt_data.len - last_txt_index;

    ctx_set_encountered (vb, ctx); // this is the ONLY place in PIZ where we set encountered
    ctx->last_encounter_was_reconstructed = reconstruct;

    // in "store per line" mode, we save one entry per line (possibly a line has no entries if it is an optional field)
    if (ctx->flags.store_per_line) 
        reconstruct_store_history (vb, ctx, last_txt_index);

    return (int32_t)ctx->last_txt_len;

missing:
    ctx->last_txt_len = 0;

    if (ctx->flags.store == STORE_INDEX) 
        ctx_set_last_value (vb, ctx, (int64_t)WORD_INDEX_MISSING);

    return reconstruct ? -1 : 0; // -1 if WORD_INDEX_MISSING - remove preceding separator
} 

static uint64_t reconstruct_state_size=0, freeze_rec_size=0;
void reconstruct_initialize (void)
{
    reconstruct_state_size = reconstruct_state_size_formula;
    freeze_rec_size = sizeof(ContextP) + reconstruct_state_size; 
}

// add a context to the state freeze, unless the context is already frozen
static void reconstruct_peek_add_ctx_to_frozen_state (VBlockP vb, ContextP ctx)
{
    ASSPIZ (!ctx->is_frozen, "Context %s is already frozen - reconstruct_peek is accessing it recursively which is not supported", ctx->tag_name);

    // add this context to the frozen state
    buf_alloc (vb, &vb->frozen_state, freeze_rec_size, 10*freeze_rec_size, char, 2, "frozen_state");

    char *reconstruct_state = BAFTc (vb->frozen_state);
    memcpy (reconstruct_state, &ctx, sizeof (ContextP)); // copy pointer to context
    memcpy (reconstruct_state + sizeof (ContextP), reconstruct_state_start(ctx), reconstruct_state_size); // copy state itself

    vb->frozen_state.len += freeze_rec_size;
    ctx->is_frozen = true;
}

static void reconstruct_peek_unfreeze_state (VBlockP vb)
{
    // unfreeze state of all involved contexts
    for (char *state=B1STc (vb->frozen_state); state < BAFTc (vb->frozen_state); state += freeze_rec_size) {
        ContextP ctx;
        memcpy (&ctx, state, sizeof (ContextP));
        memcpy (reconstruct_state_start(ctx), state + sizeof (ContextP), reconstruct_state_size);
        ctx->is_frozen = false;
    }

    // de-activate peek
    vb->frozen_state.param = 0; 
    vb->frozen_state.len   = 0; 
}

// get reconstructed text without advancing the iterator or changing last_*. context may be already reconstructed or not.
// Note: txt points into txt_data (past reconstructed or BAFT) - caller should copy it elsewhere
// LIMITATION: cannot be used with contexts that might reconstruct other contexts as that would require
// chained peeking (i.e. the secondary contexts should be peeked too, not reconstructed)
ValueType reconstruct_peek (VBlockP vb, ContextP ctx, 
                            pSTRp(txt)) // optional in / out
{
    // ASSPIZ (!vb->frozen_state.param, "Requested to peek %s, however already peeking %s. Recursive peeking not supported", 
    //         ctx->tag_name, (*B1ST(ContextP, vb->frozen_state))->tag_name);

    // case: already reconstructed in this line (or sample in the case of VCF/FORMAT)
    if (ctx_encountered (vb, ctx->did_i) && (ctx->last_encounter_was_reconstructed || (!txt && !txt_len))) {
        if (txt) *txt = last_txtx (vb, ctx);
        if (txt_len) *txt_len = ctx->last_txt_len;
        return ctx->last_value;
    }

    bool already_peeking = (bool)vb->frozen_state.param;

    // we "freeze" the reconstruction state of all contexts involved in this peek reconstruction (there could be more than
    // one context - eg items of a container, muxed contexts, "other" contexts (SNIP_OTHER etc). We freeze them by saving their state
    // in vb->frozen_state when reconstruct_from_ctx is called for this ctx.
    vb->frozen_state.param = 1; // peeking is active

    uint64_t save_txt_data_len = vb->txt_data.len;

    reconstruct_from_ctx (vb, ctx->did_i, 0, true);

    // since we are reconstructing unaccounted for data, make sure we didn't go beyond the end of txt_data (this can happen if we are close to the end
    // of the VB, and reconstructed more than OVERFLOW_SIZE allocated in piz_reconstruct_one_vb)
    ASSPIZ (vb->txt_data.len <= vb->txt_data.size, "txt_data overflow while peeking %s in vb_i=%u: len=%"PRIu64" size=%"PRIu64" last_txt_len=%u", 
            ctx->tag_name, vb->vblock_i, vb->txt_data.len, vb->txt_data.size, ctx->last_txt_len);

    if (txt) *txt = last_txtx (vb, ctx);
    if (txt_len) *txt_len = ctx->last_txt_len;

    ValueType last_value = ctx->last_value;

    // unfreeze state of all involved contexts (only in first peek in case of recursive calls to peek)
    if (!already_peeking)
        reconstruct_peek_unfreeze_state (vb);
    
    // delete data reconstructed for peeking
    vb->txt_data.len = save_txt_data_len; 

    return last_value;
}

ValueType reconstruct_peek_do (VBlockP vb, DictId dict_id, pSTRp(txt)) 
{
    ContextP ctx = ECTX (dict_id); 
    ASSPIZ (ctx, "context doesn't exist for dict_id=%s", dis_dict_id (dict_id).s);

    return reconstruct_peek (vb, ctx, txt, txt_len);
}
