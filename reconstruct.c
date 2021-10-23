// ------------------------------------------------------------------
//   reconstruct.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
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
static int64_t reconstruct_from_delta (VBlock *vb, 
                                       Context *my_ctx,   // use and store last_delta
                                       Context *base_ctx, // get last_value
                                       const char *delta_snip, unsigned delta_snip_len,
                                       bool reconstruct) 
{
    ASSERT (delta_snip, "delta_snip is NULL. vb_i=%u", vb->vblock_i);
    ASSERT (base_ctx->flags.store == STORE_INT, "reconstructing %s - calculating delta \"%.*s\" from a base of %s, but %s, doesn't have STORE_INT. vb_i=%u line_in_vb=%"PRIu64,
            my_ctx->tag_name, delta_snip_len, delta_snip, base_ctx->tag_name, base_ctx->tag_name, vb->vblock_i, vb->line_i - vb->first_line);

    if (delta_snip_len == 1 && delta_snip[0] == '-')
        my_ctx->last_delta = -2 * base_ctx->last_value.i; // negated previous value

    else if (!delta_snip_len)
        my_ctx->last_delta = -my_ctx->last_delta; // negated previous delta

    else 
        my_ctx->last_delta = (int64_t)strtoull (delta_snip, NULL, 10 /* base 10 */); // strtoull can handle negative numbers, despite its name

    int64_t new_value = base_ctx->last_value.i + my_ctx->last_delta;  
    if (reconstruct) RECONSTRUCT_INT (new_value);

    return new_value;
}

#define ASSERT_IN_BOUNDS \
    ASSERT (ctx->next_local < ctx->local.len, \
            "reconstructing txt_line=%"PRIu64" vb_i=%u: unexpected end of ctx->local data in %s (len=%u ltype=%s lcodec=%s did_i=%u)", \
            vb->line_i, vb->vblock_i, ctx->tag_name, (uint32_t)ctx->local.len, lt_name (ctx->ltype), codec_name (ctx->lcodec), ctx->did_i)

static uint32_t reconstruct_from_local_text (VBlock *vb, Context *ctx, bool reconstruct)
{
    uint32_t start = ctx->next_local; 
    ARRAY (char, data, ctx->local);

    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != 0) ctx->next_local++;
    ASSERT_IN_BOUNDS;

    char *snip = &data[start];
    uint32_t snip_len = ctx->next_local - start; 
    ctx->next_local++; /* skip the separator */ 

    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, reconstruct);

    return snip_len;
}

static int64_t reconstruct_from_local_int (VBlock *vb, Context *ctx, char separator /* 0 if none */, bool reconstruct)
{
#   define GETNUMBER(signedtype) { \
        u ## signedtype unum = NEXTLOCAL (u ## signedtype, ctx); \
        num = (int64_t)(lt_desc[ctx->ltype].is_signed ? (signedtype)unum : unum); \
        is_minus_1 = (unum == (u ## signedtype)-1); \
    }

    ASSERT_IN_BOUNDS;

    int64_t num=0;
    bool is_minus_1=false;
    switch (lt_desc[ctx->ltype].width) {
        case 4: GETNUMBER (int32_t); break;
        case 2: GETNUMBER (int16_t); break;
        case 1: GETNUMBER (int8_t ); break;
        case 8: GETNUMBER (int64_t); break;
        default: break; // never reached
    }

    // TO DO: RECONSTRUCT_INT won't reconstruct large uint64_t correctly
    if (reconstruct) { 
        if (is_minus_1 && VB_DT(DT_VCF) && dict_id_is_vcf_format_sf (ctx->dict_id)) 
            RECONSTRUCT1 ('.');
        else 
            RECONSTRUCT_INT (num);
        
        if (separator) RECONSTRUCT1 (separator);
    }

    return num;
}

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
void reconstruct_from_local_sequence (VBlock *vb, Context *ctx, STRp(snip))
{
    ASSERTNOTNULL (ctx);

    bool reconstruct = !piz_is_skip_section (vb, SEC_LOCAL, ctx->dict_id);
    uint32_t len;

    // if we have length in the snip, update vb->seq_len (for example in FASTQ, we will a snip for seq but qual will use seq_len)
    if (snip_len) vb->seq_len = atoi(snip);

    // special case: in SAM, sam_zip_qual re-wrote a '*' marking 'unavailable' as ' ' to avoid confusing with '*' as a valid quality score
    if (ctx->local.data[ctx->next_local] == ' ') {
        len = 1;
        if (reconstruct) RECONSTRUCT1 ('*');
    }
    else {
        len = vb->seq_len;
        ASSERT (ctx->next_local + len <= ctx->local.len, "reading txt_line=%"PRIu64" vb_i=%u: unexpected end of %s data", 
                vb->line_i, vb->vblock_i, ctx->tag_name);

        if (reconstruct) RECONSTRUCT (&ctx->local.data[ctx->next_local], len);
    }

    ctx->last_value.i = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
    ctx->next_local += len;
}

Context *reconstruct_get_other_ctx_from_snip (VBlockP vb, pSTRp (snip))
{
    unsigned b64_len = base64_sizeof (DictId);
    ASSERT (b64_len + 1 <= *snip_len, "snip_len=%u but expecting it to be >= %u", *snip_len, b64_len + 1);

    DictId dict_id;
    base64_decode ((*snip)+1, &b64_len, dict_id.id);

    Context *other_ctx = ECTX (dict_id);
    ASSERT (other_ctx, "Failed to get other context: snip=%.*s dict_id=%s", STRf(*snip), dis_dict_id(dict_id).s);
  
    *snip     += b64_len + 1;
    *snip_len -= b64_len + 1;
    
    return other_ctx;
}

void reconstruct_set_buddy (VBlockP vb)
{
    ASSPIZ0 (vb->buddy_line_i == NO_BUDDY, "Buddy line already set for the current line");

    ContextP buddy_ctx = ECTX (_SAM_BUDDY); // all data types using buddy are expected to have a dict_id identical to _SAM_BUDDY (but different did_i)
    int32_t num_lines_back = BGEN32 (NEXTLOCAL (uint32_t, buddy_ctx));
    vb->buddy_line_i = (vb->line_i - vb->first_line) - num_lines_back; // convert value passed (distance in lines to buddy) to 0-based buddy_line_i
    ASSPIZ (vb->buddy_line_i >= 0, "Expecting vb->buddy_line_i=%d to be non-negative", vb->buddy_line_i);
}

void reconstruct_from_buddy_get_textual_snip (VBlockP vb, ContextP ctx, pSTRp(snip))
{
    ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line");

    HistoryWord word = *ENT (HistoryWord, ctx->history, vb->buddy_line_i);
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

    *snip = ENT (char, *buf, word.char_index);
    *snip_len = word.snip_len;
}

// Copy from buddy: buddy is data that appears on a specific "buddy line", in this context or another one. Not all lines need
// Note of difference vs. lookback: with buddy, not all lines need to have the data (eg MC:Z), so the line number is constant,
// but if we had have used lookback, the lookback value would have been different between different fields.  
bool reconstruct_from_buddy (VBlock *vb, Context *ctx, STRp(snip), bool reconstruct, LastValueType *new_value)
{
    ContextP base_ctx = ctx;

    // set buddy if needed, and not already set (in BAM it is already set in sam_piz_filter).
    if (snip_len==1 && *snip == SNIP_COPY_BUDDY && vb->buddy_line_i == NO_BUDDY) 
        reconstruct_set_buddy (vb); 

    // optional: base context is different than ctx
    else if (snip_len > 1) {
        snip--; snip_len++; // reconstruct_get_other_ctx_from_snip skips the first char
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, &snip, &snip_len);
    }

    // case: numeric value 
    if (ctx->flags.store == STORE_INT) {
        ASSPIZ0 (vb->buddy_line_i >= 0, "No buddy line is set for the current line"); // for textual, we test in reconstruct_from_buddy_get_textual_snip

        new_value->i = *ENT (int64_t, base_ctx->history, vb->buddy_line_i);
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

static LastValueType reconstruct_from_lookback (VBlock *vb, Context *ctx, STRp(snip), bool reconstruct)
{   
    int64_t lookback = reconstruct_get_other_ctx_from_snip (vb, pSTRa(snip))->last_value.i;
    LastValueType value;

    if (!snip_len) { // a lookback by index
        value.i = lookback_get_index (vb, ctx, lookback);
        
        STR(back_snip);
        ctx_get_snip_by_word_index (ctx, value.i, pSTRa(back_snip));

        if (reconstruct) RECONSTRUCT (back_snip, back_snip_len);
    }
    else { // a lookback by delta
        PosType back_value = lookback_get_value (vb, ctx, lookback);
        PosType delta;
        ASSPIZ  (str_get_int (STRa(snip), &delta), "Invalid delta snip \"%.*s\"", STRf(snip));

        value.i = back_value + delta;
        
        if (reconstruct) RECONSTRUCT_INT (value.i);
    }
    
//printf ("xxx pos lookback=%u  pos=%u prev_pos=%u delta=%d\n", lookback, prev_pos+delta, prev_pos, delta);
    return value; 
}

void reconstruct_one_snip (VBlock *vb, Context *snip_ctx, 
                           WordIndex word_index, // WORD_INDEX_NONE if not used.
                           STRp(snip), bool reconstruct) // if false, calculates last_value but doesn't output to vb->txt_data)
{
    LastValueType new_value = {0};
    bool have_new_value = false;
    int64_t prev_value = snip_ctx->last_value.i;
    Context *base_ctx = snip_ctx; // this will change if the snip refers us to another data source
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
            base_ctx = reconstruct_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len

        switch (base_ctx->ltype) {
            case LT_TEXT:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                reconstruct_from_local_text (vb, base_ctx, reconstruct); // this will call us back recursively with the snip retrieved
                break;
                
            case LT_INT8 ... LT_UINT64:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                new_value.i = reconstruct_from_local_int (vb, base_ctx, 0, reconstruct);
                have_new_value = true;
                break;

            // case: the snip is taken to be the length of the sequence (or if missing, the length will be taken from vb->seq_len)
            case LT_SEQUENCE: 
                reconstruct_from_local_sequence (vb, base_ctx, STRa(snip));
                break;
                
            case LT_BITMAP:
                ASSERT_DT_FUNC (vb, reconstruct_seq);
                DT_FUNC (vb, reconstruct_seq) (vb, base_ctx, STRa(snip));
                break;

            case LT_FLOAT32:
                // TO DO - not implemented yet - see seg_float_or_not_do
                //new_value.f = reconstruct_from_local_float (vb, base_ctx, snip, snip_len, reeconstruct); // snip contains printf format for this float
                have_new_value = true;
                break;

            default: ABORT ("Error in reconstruct_one_snip of %s in vb_i=%u: Unsupported lt_type=%s (%u) for SNIP_LOOKUP or SNIP_OTHER_LOOKUP. Please upgrade to the latest version of Genozip.", 
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

    case SNIP_PAIR_LOOKUP: {
        ASSERT (snip_ctx->pair_b250, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2 vb=%u. "
                "If this file was compressed with Genozip version 9.0.12 or older use the --dts_paired flag. "
                "You can see the Genozip version used to compress it with 'genocat -w %s'", snip_ctx->tag_name, vb->vblock_i, z_name);
                
        ctx_get_next_snip (vb, snip_ctx, snip_ctx->pair_flags.all_the_same, true, &snip, &snip_len);
        reconstruct_one_snip (vb, snip_ctx, WORD_INDEX_NONE /* we can't cache pair items */, snip, snip_len, reconstruct); // might include delta etc - works because in --pair, ALL the snips in a context are PAIR_LOOKUP
        break;
    }
    case SNIP_OTHER_DELTA: 
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, pSTRa(snip)); // also updates snip and snip_len
        new_value.i = reconstruct_from_delta (vb, snip_ctx, base_ctx, STRa(snip), reconstruct); 
        have_new_value = true;
        break;

    case SNIP_PAIR_DELTA: { // used for FASTQ_GPOS - uint32_t stored in originating in the pair's local
        ASSERT (snip_ctx->pair_local, "no pair_1 local data for ctx=%s, while reconstructing pair_2 vb=%u", snip_ctx->tag_name, vb->vblock_i);
        uint32_t fastq_line_i = vb->line_i / 4 - vb->first_line; // see fastq_piz_filter for calculation
        int64_t pair_value = (int64_t) *ENT (uint32_t, snip_ctx->pair, fastq_line_i);  
        int64_t delta = (int64_t)strtoull (snip+1, NULL, 10 /* base 10 */); 
        new_value.i = pair_value + delta;
        if (reconstruct) RECONSTRUCT_INT (new_value.i); 
        have_new_value = true;
        break;
    }

    case SNIP_COPY: 
        base_ctx = (snip_len==1) ? snip_ctx : reconstruct_get_other_ctx_from_snip (vb, pSTRa(snip)); 
        RECONSTRUCT (last_txtx (vb, base_ctx), base_ctx->last_txt_len);
        new_value = base_ctx->last_value; 
        have_new_value = true;
        break;

    case SNIP_SPECIAL:
        ASSERT (snip_len >= 2, "SNIP_SPECIAL expects snip_len=%u >= 2. ctx=%s vb_i=%u line_i=%"PRIu64, 
                snip_len, snip_ctx->tag_name, vb->vblock_i, vb->line_i);
                
        uint8_t special = snip[1] - 32; // +32 was added by SPECIAL macro

        ASSERT (special < DTP (num_special), "file requires special handler %u which doesn't exist in this version of genozip - please upgrade to the latest version", special);
        ASSERT_DT_FUNC (vb, special);

        have_new_value = DT_FUNC(vb, special)[special](vb, snip_ctx, snip+2, snip_len-2, &new_value, reconstruct);  
        break;

    case SNIP_REDIRECTION: 
        base_ctx = reconstruct_get_other_ctx_from_snip (vb, pSTRa(snip)); // also updates snip and snip_len
        reconstruct_from_ctx (vb, base_ctx->did_i, 0, reconstruct);
        break;
    
    case SNIP_DUAL: {
        str_split (&snip[1], snip_len-1, 2, SNIP_DUAL, subsnip, true);
        ASSPIZ (n_subsnips==2, "Invalid SNIP_DUAL snip in ctx=%s", snip_ctx->tag_name);

        if (vb->vb_coords == DC_PRIMARY) // recursively call for each side 
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
        ctx_set_last_value (vb, snip_ctx, new_value);
    else
        ctx_set_encountered_in_line (snip_ctx);

    // note: if store_delta, we do a self-delta. this overrides last_delta set by the delta snip which could be against a different
    // base_ctx. note: when Seg sets last_delta, it must also set store=STORE_INT
    if (store_delta) 
        snip_ctx->last_delta = new_value.i - prev_value;
}

// returns reconstructed length or -1 if snip is missing and previous separator should be deleted
int32_t reconstruct_from_ctx_do (VBlock *vb, DidIType did_i, 
                                 char sep, // if non-zero, outputs after the reconstruction
                                 bool reconstruct, // if false, calculates last_value but doesn't output to vb->txt_data
                                 const char *func)
{
    ASSERT (did_i < vb->num_contexts, "called from: %s: did_i=%u out of range: vb->num_contexts=%u for vb_i=%u", 
            func, did_i, vb->num_contexts, vb->vblock_i);

    Context *ctx = CTX(did_i);

    ASSERT0 (ctx->dict_id.num || ctx->did_i != DID_I_NONE, "ctx not initialized (dict_id=0)");

    // update ctx, if its an alias (only for primary field aliases as they have contexts, other alias don't have ctx)
    if (!ctx->dict_id.num) 
        ctx = CTX(ctx->did_i); // ctx->did_i is different than did_i if its an alias

    uint32_t last_txt_index = (uint32_t)vb->txt_data.len;

    // case: we have b250 data
    if (ctx->b250.len ||
        (!ctx->b250.len && !ctx->local.len && ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        STR0(snip);
        WordIndex word_index = LOAD_SNIP(ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (!snip) {
            ctx->last_txt_len = 0;

            if (ctx->flags.store == STORE_INDEX) 
                ctx_set_last_value (vb, ctx, (int64_t)WORD_INDEX_MISSING);

            return reconstruct ? -1 : 0; // -1 if WORD_INDEX_MISSING - remove preceding separator
        }

        reconstruct_one_snip (vb, ctx, word_index, STRa(snip), reconstruct);        

        // for backward compatability with v8-11 that didn't yet have flags.store = STORE_INDEX for CHROM
        if (did_i == CHROM) { // NOTE: CHROM cannot have aliases, because looking up the did_i by dict_id will lead to CHROM, and this code will be executed for a non-CHROM field
            vb->chrom_node_index = vb->last_index (CHROM) = word_index;
            vb->chrom_name       = snip; // used for reconstruction from external reference
            vb->chrom_name_len   = snip_len;
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
            reconstruct_from_local_sequence (vb, ctx, NULL, 0); break;
                
        case LT_BITMAP:
            ASSERT_DT_FUNC (vb, reconstruct_seq);
            DT_FUNC (vb, reconstruct_seq) (vb, ctx, NULL, 0);
            break;
        
        case LT_TEXT:
            reconstruct_from_local_text (vb, ctx, reconstruct); break;

        default:
            ABORT ("Invalid ltype=%u in ctx=%s of vb_i=%u line_i=%"PRIu64, ctx->ltype, ctx->tag_name, vb->vblock_i, vb->line_i);
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

    else ASSERT (flag.missing_contexts_allowed,
                 "Error in reconstruct_from_ctx_do: ctx %s/%s has no data (dict, b250 or local) in vb_i=%u line_i=%"PRIu64" did_i=%u ctx->did=%u ctx->dict_id=%s", 
                 dtype_name_z (ctx->dict_id), ctx->tag_name, vb->vblock_i, vb->line_i, did_i, ctx->did_i, dis_dict_id (ctx->dict_id).s);

    if (sep && reconstruct) RECONSTRUCT1 (sep); 

    ctx->last_txt_index = last_txt_index;
    ctx->last_txt_len   = (uint32_t)vb->txt_data.len - ctx->last_txt_index;
    ctx->last_line_i    = vb->line_i; // reconstructed on this line

    // in "store per line" mode, we save one entry per line (possibly a line has no entries if it is an optional field)
    if (ctx->flags.store_per_line) {
        // case: store last integer value
        if (ctx->flags.store == STORE_INT) 
            *ENT (int64_t, ctx->history, vb->line_i - vb->first_line) = ctx->last_value.i;
        
        // case: textual value will remain in txt_data - just point to it
        else if (!flag.maybe_lines_dropped_by_reconstructor) 
            *ENT (HistoryWord, ctx->history, vb->line_i - vb->first_line) = 
                (HistoryWord){ .lookup = LookupTxtData, .char_index = last_txt_index, .snip_len = ctx->last_txt_len };
        
        // case: textual value might be removed from txt_data - copy it
        else {
            *ENT (HistoryWord, ctx->history, vb->line_i - vb->first_line) = 
                (HistoryWord){ .lookup = LookupPerLine, .char_index = ctx->per_line.len, .snip_len = ctx->last_txt_len };

            if (ctx->last_txt_len) {
                buf_add_more (vb, &ctx->per_line, ENT (char, vb->txt_data, last_txt_index), ctx->last_txt_len, "per_line");
                NEXTENT (char, ctx->per_line) = 0; // nul-terminate
            }
        }
    }

    return (int32_t)ctx->last_txt_len;
} 

// get reconstructed text without advancing the iterator or changing last_*. context may be already reconstructed or not.
// Note: txt points into txt_data (past reconstructed or AFTERENT) - caller should copy it elsewhere
LastValueType reconstruct_peek (VBlock *vb, Context *ctx, 
                                pSTRp(txt)) // optional in / out
{
    // case: already reconstructed 
    if (ctx->last_line_i == vb->line_i) {
        if (txt) *txt = last_txtx (vb, ctx);
        if (txt_len) *txt_len = ctx->last_txt_len;
        return ctx->last_value;
    }

    // case: ctx is not reconstructed yet in this last - we reconstruct, but then roll back state
    char reconstruct_state[reconstruct_state_size(ctx)];
    memcpy (reconstruct_state, reconstruct_state_start(ctx), reconstruct_state_size(ctx)); // save
    uint64_t save_txt_data_len = vb->txt_data.len;

    reconstruct_from_ctx (vb, ctx->did_i, 0, true);

    // since we are reconstructing unaccounted for data, make sure we didn't go beyond the end of txt_data (this can happen if we are close to the end
    // of the VB, and reconstructed more than OVERFLOW_SIZE allocated in piz_reconstruct_one_vb)
    ASSERT (vb->txt_data.len <= vb->txt_data.size, "txt_data overflow while peeking %s in vb_i=%u: len=%"PRIu64" size=%"PRIu64" last_txt_len=%u", 
            ctx->tag_name, vb->vblock_i, vb->txt_data.len, vb->txt_data.size, ctx->last_txt_len);

    if (txt) *txt = last_txtx (vb, ctx);
    if (txt_len) *txt_len = ctx->last_txt_len;

    LastValueType last_value = ctx->last_value;

    memcpy (reconstruct_state_start(ctx), reconstruct_state, reconstruct_state_size(ctx)); // restore
    vb->txt_data.len = save_txt_data_len;

    return last_value;
}

LastValueType reconstruct_peek_do (VBlockP vb, DictId dict_id, pSTRp(txt)) 
{
    Context *ctx = ECTX (dict_id); 
    ASSPIZ (ctx, "context doesn't exist for dict_id=%s", dis_dict_id (dict_id).s);

    return reconstruct_peek (vb, ctx, txt, txt_len);
}
