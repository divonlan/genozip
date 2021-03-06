// ------------------------------------------------------------------
//   reconstruct.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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
            my_ctx->name, delta_snip_len, delta_snip, base_ctx->name, base_ctx->name, vb->vblock_i, vb->line_i - vb->first_line);

    if (delta_snip_len == 1 && delta_snip[0] == '-')
        my_ctx->last_delta = -2 * base_ctx->last_value.i; // negated previous value

    else if (!delta_snip_len)
        my_ctx->last_delta = -my_ctx->last_delta; // negated previous delta

    else 
        my_ctx->last_delta = (int64_t)strtoull (delta_snip, NULL, 10 /* base 10 */); // strtoull can handle negative numbers, despite its name

    int64_t new_value = base_ctx->last_value.i + my_ctx->last_delta;  
    if (reconstruct) { RECONSTRUCT_INT (new_value) };

    return new_value;
}

#define ASSERT_IN_BOUNDS \
    ASSERT (ctx->next_local < ctx->local.len, \
            "reconstructing txt_line=%"PRIu64" vb_i=%u: unexpected end of ctx->local data in %s (len=%u ltype=%s lcodec=%s)", \
            vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len, lt_name (ctx->ltype), codec_name (ctx->lcodec))

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

static int64_t reconstruct_from_local_int (VBlock *vb, Context *ctx, char seperator /* 0 if none */, bool reconstruct)
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
        if (is_minus_1 && vb->data_type == DT_VCF && dict_id_is_vcf_format_sf (ctx->dict_id)) {
            RECONSTRUCT1 ('.');
        } else {
            RECONSTRUCT_INT (num);
        }
        if (seperator) RECONSTRUCT1 (seperator);
    }

    return num;
}

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
static void reconstruct_from_local_sequence (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
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
                vb->line_i, vb->vblock_i, ctx->name);

        if (reconstruct) RECONSTRUCT (&ctx->local.data[ctx->next_local], len);
    }

    ctx->last_value.i = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
    ctx->next_local += len;
}

static Context *piz_get_other_ctx_from_snip (VBlockP vb, const char **snip, unsigned *snip_len)
{
    unsigned b64_len = base64_sizeof (DictId);
    ASSERT (b64_len + 1 <= *snip_len, "snip_len=%u but expecting it to be >= %u", *snip_len, b64_len + 1);

    DictId dict_id;
    base64_decode ((*snip)+1, &b64_len, dict_id.id);

    Context *other_ctx = ctx_get_existing_ctx (vb, dict_id);

    *snip     += b64_len + 1;
    *snip_len -= b64_len + 1;
    
    return other_ctx;
}

void reconstruct_one_snip (VBlock *vb, Context *snip_ctx, 
                           WordIndex word_index, // WORD_INDEX_NONE if not used.
                           const char *snip, unsigned snip_len,
                           bool reconstruct) // if false, calculates last_value but doesn't output to vb->txt_data)
{
    if (!snip_len) return; // nothing to do
    
    LastValueType new_value = {0};
    bool have_new_value = false;
    Context *base_ctx = snip_ctx; // this will change if the snip refers us to another data source
    StoreType store_type = snip_ctx->flags.store;

    switch (snip[0]) {

    // display the rest of the snip first, and then the lookup up text.
    case SNIP_LOOKUP:
    case SNIP_OTHER_LOOKUP: {
        if (snip[0] == SNIP_LOOKUP) 
            { snip++; snip_len--; }
        else 
            // we are request to reconstruct from another ctx
            base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len

        switch (base_ctx->ltype) {
            case LT_TEXT:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                reconstruct_from_local_text (vb, base_ctx, reconstruct); // this will call us back recursively with the snip retrieved
                break;
                
            case LT_INT8 ...LT_UINT64:
                if (reconstruct && snip_len) RECONSTRUCT (snip, snip_len); // reconstruct this snip before adding the looked up data
                new_value.i = reconstruct_from_local_int (vb, base_ctx, 0, reconstruct);
                have_new_value = true;
                break;

            // case: the snip is taken to be the length of the sequence (or if missing, the length will be taken from vb->seq_len)
            case LT_SEQUENCE: 
                reconstruct_from_local_sequence (vb, base_ctx, snip, snip_len);
                break;
                
            case LT_BITMAP:
                ASSERT_DT_FUNC (vb, reconstruct_seq);
                DT_FUNC (vb, reconstruct_seq) (vb, base_ctx, snip, snip_len);
                break;

            case LT_FLOAT32:
                // TO DO - not implemented yet - see seg_float_or_not_do
                //new_value.f = reconstruct_from_local_float (vb, base_ctx, snip, snip_len, reeconstruct); // snip contains printf format for this float
                have_new_value = true;
                break;

            default: ABORT ("Error in reconstruct_one_snip: Unsupported lt_type=%u for SNIP_LOOKUP or SNIP_OTHER_LOOKUP", base_ctx->ltype);
        }

        break;
    }
    case SNIP_PAIR_LOOKUP: {
        ASSERT (snip_ctx->pair_b250, "no pair_1 b250 data for ctx=%s, while reconstructing pair_2 vb=%u", snip_ctx->name, vb->vblock_i);
        ctx_get_next_snip (vb, snip_ctx, snip_ctx->pair_flags.all_the_same, true, &snip, &snip_len);
        reconstruct_one_snip (vb, snip_ctx, WORD_INDEX_NONE /* we can't cache pair items */, snip, snip_len, reconstruct); // might include delta etc - works because in --pair, ALL the snips in a context are PAIR_LOOKUP
        break;
    }
    case SNIP_SELF_DELTA:
        new_value.i = reconstruct_from_delta (vb, snip_ctx, base_ctx, snip+1, snip_len-1, reconstruct);
        have_new_value = true;
        break;

    case SNIP_OTHER_DELTA: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        new_value.i = reconstruct_from_delta (vb, snip_ctx, base_ctx, snip, snip_len, reconstruct); 
        have_new_value = true;
        break;

    case SNIP_PAIR_DELTA: { // used for FASTQ_GPOS - uint32_t stored in originating in the pair's local
        ASSERT (snip_ctx->pair_local, "no pair_1 local data for ctx=%s, while reconstructing pair_2 vb=%u", snip_ctx->name, vb->vblock_i);
        uint32_t fastq_line_i = vb->line_i / 4 - vb->first_line; // see fastq_piz_filter for calculation
        int64_t pair_value = (int64_t) *ENT (uint32_t, snip_ctx->pair, fastq_line_i);  
        int64_t delta = (int64_t)strtoull (snip+1, NULL, 10 /* base 10 */); 
        new_value.i = pair_value + delta;
        if (reconstruct) { RECONSTRUCT_INT (new_value.i); }
        have_new_value = true;
        break;
    }

    case SNIP_OTHER_COPY: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        RECONSTRUCT (last_txtx (vb, base_ctx), base_ctx->last_txt_len);
        new_value = base_ctx->last_value; 
        have_new_value = true;
        break;

    case SNIP_CONTAINER:
        new_value = container_reconstruct (vb, snip_ctx, word_index, snip+1, snip_len-1);
        have_new_value = true;
        break;

    case SNIP_SPECIAL:
        ASSERT (snip_len >= 2, "SNIP_SPECIAL expects snip_len=%u >= 2. ctx=%s vb_i=%u line_i=%"PRIu64, 
                snip_len, snip_ctx->name, vb->vblock_i, vb->line_i);
                
        uint8_t special = snip[1] - 32; // +32 was added by SPECIAL macro

        ASSERT (special < DTP (num_special), "file requires special handler %u which doesn't exist in this version of genozip - please upgrade to the latest version", special);
        ASSERT_DT_FUNC (vb, special);

        have_new_value = DT_FUNC(vb, special)[special](vb, snip_ctx, snip+2, snip_len-2, &new_value, reconstruct);  
        break;

    case SNIP_REDIRECTION: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        reconstruct_from_ctx (vb, base_ctx->did_i, 0, reconstruct);
        break;
    
    case SNIP_DONT_STORE:
        store_type = STORE_NONE; // override store and fall through
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

    // update last_value if needed
    if (have_new_value && store_type) // note: we store in our own context, NOT base (a context, eg FORMAT/DP, sometimes serves as a base_ctx of MIN_DP and sometimes as the snip_ctx for INFO_DP)
        ctx_set_last_value (vb, snip_ctx, new_value);
    else
        ctx_set_encountered_in_line (snip_ctx);
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

    ctx->last_txt_index = (uint32_t)vb->txt_data.len;

    // case: we have b250 data
    if (ctx->b250.len ||
        (!ctx->b250.len && !ctx->local.len && ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        DECLARE_SNIP;
        WordIndex word_index = LOAD_SNIP(ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (!snip) {
            ctx->last_txt_len = 0;
            return reconstruct ? -1 : 0; // -1 if WORD_INDEX_MISSING_SF - remove preceding separator
        }

        reconstruct_one_snip (vb, ctx, word_index, snip, snip_len, reconstruct);        

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
            ABORT ("Invalid ltype=%u in ctx=%s of vb_i=%u line_i=%"PRIu64, ctx->ltype, ctx->name, vb->vblock_i, vb->line_i);
        }
    }

    // in case of LT_BITMAP, it is it is ok if the bitmap is empty and all the data is in NONREF (e.g. unaligned SAM)
    else if (ctx->ltype == LT_BITMAP && (ctx+1)->local.len) {
        ASSERT_DT_FUNC (vb, reconstruct_seq);
        DT_FUNC (vb, reconstruct_seq) (vb, ctx, NULL, 0);
    }

    // case: the entire VB was just \n - so seg dropped the ctx
    // note: for backward compatability with 8.0. for files compressed by 8.1+, it will be handled via the all_the_same mechanism
    else if (ctx->did_i == DTF(eol)) {
        if (reconstruct) { RECONSTRUCT1('\n'); }
    }

    else ASSERT (flag.missing_contexts_allowed,
                 "Error in reconstruct_from_ctx_do: ctx %s has no data (dict, b250 or local) in vb_i=%u line_i=%"PRIu64" did_i=%u ctx->did=%u ctx->dict_id=%s", 
                 ctx->name, vb->vblock_i, vb->line_i, did_i, ctx->did_i, dis_dict_id (ctx->dict_id).s);

    if (sep && reconstruct) RECONSTRUCT1 (sep); 

    ctx->last_txt_len = (uint32_t)vb->txt_data.len - ctx->last_txt_index;
    ctx->last_line_i  = vb->line_i; // reconstructed on this line

    return (int32_t)ctx->last_txt_len;
} 

// get reconstructed text without advancing the iterator or changing last_*. context may be already reconstructed or not.
// Note: txt points into txt_data (past reconstructed or AFTERENT) - caller should copy it elsewhere
LastValueType reconstruct_peek (VBlock *vb, Context *ctx, 
                                const char **txt, unsigned *txt_len) // optional in / out
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
            ctx->name, vb->vblock_i, vb->txt_data.len, vb->txt_data.size, ctx->last_txt_len);

    if (txt) *txt = last_txtx (vb, ctx);
    if (txt_len) *txt_len = ctx->last_txt_len;

    LastValueType last_value = ctx->last_value;

    memcpy (reconstruct_state_start(ctx), reconstruct_state, reconstruct_state_size(ctx)); // restore
    vb->txt_data.len = save_txt_data_len;

    return last_value;
}

LastValueType reconstruct_peek__do (VBlockP vb, DictId dict_id, const char **txt, unsigned *txt_len) 
{
    Context *ctx = ctx_get_existing_ctx (vb, dict_id); 
    ASSPIZ (ctx, "context doesn't exist for dict_id=%s", dis_dict_id (dict_id).s);

    return reconstruct_peek (vb, ctx, txt, txt_len);
}
