// ------------------------------------------------------------------
//   reconstruct.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
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
    ASSERTE (delta_snip, "delta_snip is NULL. vb_i=%u", vb->vblock_i);
    ASSERTE (base_ctx->flags.store == STORE_INT, "attempting calculate delta from a base of \"%s\", but this context doesn't have STORE_INT",
             base_ctx->name);

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
    ASSERTE (ctx->next_local < ctx->local.len, \
             "reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u ltype=%s lcodec=%s)", \
             vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len, lt_name (ctx->ltype), codec_name (ctx->lcodec))

static uint32_t reconstruct_from_local_text (VBlock *vb, Context *ctx, bool reconstruct)
{
    uint32_t start = ctx->next_local; 
    ARRAY (char, data, ctx->local);

    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != SNIP_SEP) ctx->next_local++;
    ASSERT_IN_BOUNDS;

    char *snip = &data[start];
    uint32_t snip_len = ctx->next_local - start; 
    ctx->next_local++; /* skip the tab */ 

    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, reconstruct);

    return snip_len;
}

int64_t reconstruct_from_local_int (VBlock *vb, Context *ctx, char seperator /* 0 if none */, bool reconstruct)
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
    ASSERTE0 (ctx, "ctx is NULL");

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
        ASSERTE (ctx->next_local + len <= ctx->local.len, "reading txt_line=%u vb_i=%u: unexpected end of %s data", 
                 vb->line_i, vb->vblock_i, ctx->name);

        if (reconstruct) RECONSTRUCT (&ctx->local.data[ctx->next_local], len);
    }

    ctx->last_value.i = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
    ctx->next_local += len;
}

static Context *piz_get_other_ctx_from_snip (VBlockP vb, const char **snip, unsigned *snip_len)
{
    unsigned b64_len = base64_sizeof (DictId);
    ASSERTE (b64_len + 1 <= *snip_len, "snip_len=%u but expecting it to be >= %u", *snip_len, b64_len + 1);

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

        // case: LOCAL is not LT_SEQUENCE/LT_BITMAP - we reconstruct this snip before adding the looked up data
        if (snip_len && base_ctx->ltype != LT_SEQUENCE && base_ctx->ltype != LT_BITMAP) 
            if (reconstruct) RECONSTRUCT (snip, snip_len);
        
        switch (base_ctx->ltype) {
            case LT_TEXT:
                reconstruct_from_local_text (vb, base_ctx, reconstruct); // this will call us back recursively with the snip retrieved
                break;
                
            case LT_INT8 ...LT_UINT64:
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

            default: ABORT ("Error in reconstruct_one_snip: Unsupported lt_type=%u for SNIP_LOOKUP or SNIP_OTHER_LOOKUP", base_ctx->ltype);
        }

        break;
    }
    case SNIP_PAIR_LOOKUP: {
        ASSERTE (snip_ctx->pair_b250_iter.next_b250, "no pair_1 data for ctx=%s, while reconstructing pair_2 vb=%u", 
                 snip_ctx->name, vb->vblock_i);

        ctx_get_next_snip (vb, snip_ctx, snip_ctx->pair_flags.all_the_same, &snip_ctx->pair_b250_iter, &snip, &snip_len);
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
        uint32_t fastq_line_i = vb->line_i / 4 - vb->first_line; // see fastq_piz_filter for calculation
        int64_t pair_value = (int64_t) *ENT (uint32_t, snip_ctx->pair, fastq_line_i);  
        int64_t delta = (int64_t)strtoull (snip+1, NULL, 10 /* base 10 */); 
        new_value.i = pair_value + delta;
        if (reconstruct) { RECONSTRUCT_INT (new_value.i); }
        have_new_value = true;
        break;
    }

    case SNIP_CONTAINER:
        new_value = container_reconstruct (vb, snip_ctx, word_index, snip+1, snip_len-1);
        have_new_value = true;
        break;

    case SNIP_SPECIAL:
        ASSERTE (snip_len >= 2, "SNIP_SPECIAL expects snip_len >= 2. ctx=%s", snip_ctx->name);
        uint8_t special = snip[1] - 32; // +32 was added by SPECIAL macro

        ASSERTE (special < DTP (num_special), "file requires special handler %u which doesn't exist in this version of genounzip - please upgrade to the latest version", special);
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
        snip_ctx->last_value = new_value;

    snip_ctx->last_line_i = vb->line_i;
}

// returns reconstructed length or -1 if snip is missing and item's operator should not be emitted
int32_t reconstruct_from_ctx_do (VBlock *vb, DidIType did_i, 
                                 char sep, // if non-zero, outputs after the reconstruction
                                 bool reconstruct) // if false, calculates last_value but doesn't output to vb->txt_data
{
    ASSERTE (did_i < vb->num_contexts, "did_i=%u out of range: vb->num_contexts=%u", did_i, vb->num_contexts);

    Context *ctx = &vb->contexts[did_i];

    ASSERTE0 (ctx->dict_id.num || ctx->did_i != DID_I_NONE, "ctx not initialized (dict_id=0)");

    // update ctx, if its an alias (only for primary field aliases as they have contexts, other alias don't have ctx)
    if (!ctx->dict_id.num) 
        ctx = &vb->contexts[ctx->did_i]; // ctx->did_i is different than did_i if its an alias

    ctx->last_txt = (uint32_t)vb->txt_data.len;

    // case: we have b250 data
    if (ctx->b250.len ||
        (!ctx->b250.len && !ctx->local.len && ctx->dict.len)) {          
        DECLARE_SNIP;
        uint32_t word_index = LOAD_SNIP(ctx->did_i); // if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (!snip) {
            ctx->last_txt_len = 0;
            return -1; // WORD_INDEX_MISSING_SF - remove preceding separator
        }

        reconstruct_one_snip (vb, ctx, word_index, snip, snip_len, reconstruct);        

        // handle chrom and pos to determine whether this line should be grepped-out in case of --regions
        if (did_i == CHROM) { // NOTE: CHROM cannot have aliases, because looking up the did_i by dict_id will lead to CHROM, and this code will be executed for a non-CHROM field
            vb->chrom_node_index = word_index;
            vb->chrom_name       = snip; // used for reconstruction from external reference
            vb->chrom_name_len   = snip_len;
        }

        if (flag.regions && did_i == DTF(pos) && !regions_is_site_included (vb->chrom_node_index, ctx->last_value.i)) 
            vb->dont_show_curr_line = true;
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
            ABORT ("Invalid ltype=%u in ctx=%s of vb_i=%u line_i=%u", ctx->ltype, ctx->name, vb->vblock_i, vb->line_i);
        }
    }

    // in case of LT_BITMAP, it is it is ok if the bitmap is empty and all the data is in NONREF (e.g. unaligned SAM)
    else if (ctx->ltype == LT_BITMAP && (ctx+1)->local.len) {
        ASSERT_DT_FUNC (vb, reconstruct_seq);
        DT_FUNC (vb, reconstruct_seq) (vb, ctx, NULL, 0);
    }

    // case: the entire VB was just \n - so seg dropped the ctx
    // note: for backward compatability with 8.0. for files compressed by 8.1+, it will be handled via a dictionary but no b250
    else if (ctx->did_i == DTF(eol)) {
        if (reconstruct) { RECONSTRUCT1('\n'); }
    }

    else ASSERTE (flag.show_sex || flag.show_coverage, // in --show-sex/coverage, we filtered out most contexts in sam_piz_is_skip_section, so this is expected
                  "Error in reconstruct_from_ctx_do: ctx %s has no data (dict, b250 or local) in vb_i=%u line_i=%u did_i=%u ctx->did=%u ctx->dict_id=%s", 
                  ctx->name, vb->vblock_i, vb->line_i, did_i, ctx->did_i, dis_dict_id (ctx->dict_id).s);

    if (sep && reconstruct) RECONSTRUCT1 (sep); 

    ctx->last_txt_len = (uint32_t)vb->txt_data.len - ctx->last_txt;

    return (int32_t)ctx->last_txt_len;
} 
