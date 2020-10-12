
// ------------------------------------------------------------------
//   piz.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "zfile.h"
#include "txtfile.h"
#include "vblock.h"
#include "base250.h"
#include "base64.h"
#include "dispatcher.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "piz.h"
#include "sections.h"
#include "random_access.h"
#include "regions.h"
#include "strings.h"
#include "seg.h"
#include "dict_id.h"
#include "reference.h"
#include "refhash.h"
#include "progress.h"
#include "domqual.h"
#include "profiler.h"
#include "stats.h"

// Compute threads: decode the delta-encoded value of the POS field, and returns the new last_pos
// Special values:
// "-" - negated previous value
// ""  - negated previous delta
static int64_t piz_reconstruct_from_delta (VBlock *vb, 
                                           Context *my_ctx,   // use and store last_delta
                                           Context *base_ctx, // get last_value
                                           const char *delta_snip, unsigned delta_snip_len) 
{
    ASSERT (delta_snip, "Error in piz_reconstruct_from_delta: delta_snip is NULL. vb_i=%u", vb->vblock_i);
    ASSERT (ctx_is_store (base_ctx, CTX_FL_STORE_INT), "Error in piz_reconstruct_from_delta: attempting calculate delta from a base of \"%s\", but this context doesn't have CTX_FL_STORE_INT",
            base_ctx->name);

    if (delta_snip_len == 1 && delta_snip[0] == '-')
        my_ctx->last_delta = -2 * base_ctx->last_value.i; // negated previous value

    else if (!delta_snip_len)
        my_ctx->last_delta = -my_ctx->last_delta; // negated previous delta

    else 
        my_ctx->last_delta = (int64_t)strtoull (delta_snip, NULL, 10 /* base 10 */); // strtoull can handle negative numbers, despite its name

    int64_t new_value = base_ctx->last_value.i + my_ctx->last_delta;  
    RECONSTRUCT_INT (new_value);

    return new_value;
}

#define ASSERT_IN_BOUNDS \
    ASSERT (ctx->next_local < ctx->local.len, \
            "Error reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u)", \
            vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len);

static uint32_t piz_reconstruct_from_local_text (VBlock *vb, Context *ctx)
{
    uint32_t start = ctx->next_local; 
    ARRAY (char, data, ctx->local);

    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != SNIP_SEP) ctx->next_local++;
    ASSERT_IN_BOUNDS;

    char *snip = &data[start];
    uint32_t snip_len = ctx->next_local - start; 
    ctx->next_local++; /* skip the tab */ 

    piz_reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len);

    return snip_len;
}

static int64_t piz_reconstruct_from_local_int (VBlock *vb, Context *ctx, char seperator /* 0 if none */)
{
#   define GETNUMBER(signedtype) { \
        u ## signedtype unum = NEXTLOCAL (u ## signedtype, ctx); \
        num = (int64_t)(lt_desc[ctx->ltype].is_signed ? (signedtype)unum : unum); \
    }

    ASSERT_IN_BOUNDS;

    int64_t num=0;
    switch (lt_desc[ctx->ltype].width) {
        case 4: GETNUMBER (int32_t); break;
        case 2: GETNUMBER (int16_t); break;
        case 1: GETNUMBER (int8_t ); break;
        case 8: GETNUMBER (int64_t); break;
        default: break; // never reached
    }

    // TO DO: RECONSTRUCT_INT won't reconstruct large uint64_t correctly
    RECONSTRUCT_INT (num);
    if (seperator) RECONSTRUCT1 (seperator);

    return num;
}

// reconstruct a allele value from haplotype matrix (used in VCF)
static void piz_reconstruct_from_local_ht (VBlock *vb, Context *ctx)
{
    if (vb->dont_show_curr_line) return;

    // get one row of the haplotype matrix for this line into vb->ht_one_array if we don't have it already
    if (vb->ht_one_array_line_i != vb->line_i) {
        comp_ht_piz_get_one_line (vb);
        vb->ht_one_array.len = 0; // length of data consumed
        vb->ht_one_array_line_i = vb->line_i;
    }

    // find next allele - skipping unused spots ('*')
    uint8_t ht = '*';
    while (ht == '*' && vb->ht_one_array.len < vb->num_haplotypes_per_line)
        ht = *ENT(uint8_t, vb->ht_one_array, vb->ht_one_array.len++);

    if (ht == '.' || IS_DIGIT(ht)) 
        RECONSTRUCT1 (ht);
    
    else if (ht == '*') 
        ABORT ("Error in piz_reconstruct_from_local_ht: reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u)", 
               vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len)
    
    else { // allele 10 to 99 (ascii 58 to 147)
        RECONSTRUCT_INT (ht - '0');
    }
}

// two options: 1. the length maybe given (textually) in snip/snip_len. in that case, it is used and vb->seq_len is updated.
// if snip_len==0, then the length is taken from seq_len.
static void piz_reconstruct_from_local_sequence (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    ASSERT0 (ctx, "Error in piz_reconstruct_from_local_sequence: ctx is NULL");

    bool reconstruct = !piz_is_skip_section (vb, SEC_LOCAL, ctx->dict_id);
    uint32_t len;

    // if we have length in the snip, update vb->seq_len (for example in FASTQ, we will a snip for seq but qual will use seq_len)
    if (snip_len) vb->seq_len = atoi(snip);

    // special case: it is "*" that was written to " " - we reconstruct it
    if (ctx->local.data[ctx->next_local] == ' ') {
        len = 1;
        if (reconstruct) RECONSTRUCT1 ('*');
    }
    else {
        len = vb->seq_len;
        ASSERT (ctx->next_local + len <= ctx->local.len, "Error reading txt_line=%u vb_i=%u: unexpected end of %s data", 
                vb->line_i, vb->vblock_i, ctx->name);

        if (reconstruct) RECONSTRUCT (&ctx->local.data[ctx->next_local], len);
    }

    ctx->last_value.i = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
    ctx->next_local += len;
}

static inline void piz_reconstruct_structured_prefix (VBlockP vb, const char **prefixes, uint32_t *prefixes_len)
{
    if (! (*prefixes_len)) return; // nothing to do
    
    const char *start = *prefixes;
    while (**prefixes != SNIP_STRUCTURED) (*prefixes)++; // prefixes are terminated by SNIP_STRUCTURED
    uint32_t len = (unsigned)((*prefixes) - start);

    RECONSTRUCT (start, len);

    (*prefixes)++; // skip SNIP_STRUCTURED seperator
    (*prefixes_len) -= len + 1;
}

void piz_reconstruct_structured_do (VBlock *vb, DictId dict_id, const Structured *st, const char *prefixes, uint32_t prefixes_len)
{
    TimeSpecType profiler_timer={0}; 
    if (flag_show_time && (st->flags & STRUCTURED_TOPLEVEL)) 
        clock_gettime (CLOCK_REALTIME, &profiler_timer);

    // structured wide prefix - it will be missing if Structured has no prefixes, or empty if it has only items prefixes
    piz_reconstruct_structured_prefix (vb, &prefixes, &prefixes_len); // item prefix (we will have one per item or none at all)

    ASSERT (DTP (structured_filter) || !(st->flags & STRUCTURED_FILTER_REPEATS) || !(st->flags & STRUCTURED_FILTER_ITEMS), 
            "Error: data_type=%s doesn't support structured_filter", dt_name (vb->data_type));

    for (uint32_t rep_i=0; rep_i < st->repeats; rep_i++) {

        // case this is the top-level snip
        if (st->flags & STRUCTURED_TOPLEVEL) {
            vb->line_i = vb->first_line + rep_i;
            vb->line_start = vb->txt_data.len;
            vb->dont_show_curr_line = false; 
        }
    
        if ((st->flags & STRUCTURED_FILTER_REPEATS) && !DTP (structured_filter) (vb, dict_id, st, rep_i, -1)) continue; // repeat is filtered out

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t item_prefixes_len = prefixes_len;

        for (unsigned i=0; i < st->num_items; i++) {

            if ((st->flags & STRUCTURED_FILTER_ITEMS) && !DTP (structured_filter) (vb, dict_id, st, rep_i, i)) continue; // item is filtered out

            piz_reconstruct_structured_prefix (vb, &item_prefixes, &item_prefixes_len); // item prefix

            const StructuredItem *item = &st->items[i];
            int32_t reconstructed_len=0;
            if (item->dict_id.num)  // not a prefix-only item
                reconstructed_len = piz_reconstruct_from_ctx (vb, item->did_i, 0);
            
            if (reconstructed_len == -1 && i > 0)  // not WORD_INDEX_MISSING_SF - delete previous item's separator
                vb->txt_data.len -= ((item-1)->seperator[0] != 0) + ((item-1)->seperator[1] != 0);

            // emit this item's separator even if this item is missing - next item will delete it if also missing (last item in Samples doesn't have a seperator)
            if (item->seperator[0]) RECONSTRUCT1 (item->seperator[0]);
            if (item->seperator[1]) RECONSTRUCT1 (item->seperator[1]);
        }

        if (rep_i+1 < st->repeats || !(st->flags & STRUCTURED_DROP_FINAL_REPEAT_SEP)) {
            if (st->repsep[0]) RECONSTRUCT1 (st->repsep[0]);
            if (st->repsep[1]) RECONSTRUCT1 (st->repsep[1]);
        }

        // in top level: after consuming the line's data, if it is not to be outputted - trim txt_data back to start of line
        if ((st->flags & STRUCTURED_TOPLEVEL) && vb->dont_show_curr_line) 
            vb->txt_data.len = vb->line_start; 
    }

    if (st->flags & STRUCTURED_DROP_FINAL_ITEM_SEP)
        vb->txt_data.len -= (st->items[st->num_items-1].seperator[0] != 0) + 
                            (st->items[st->num_items-1].seperator[1] != 0);

    if (st->flags & STRUCTURED_TOPLEVEL)   
        COPY_TIMER (vb->profile.piz_reconstruct_vb);
}

static void piz_reconstruct_structured (VBlock *vb, Context *ctx, WordIndex word_index, const char *snip, unsigned snip_len)
{
    ASSERT (snip_len <= base64_sizeof(Structured), "Error in piz_reconstruct_structured: snip_len=%u exceed base64_sizeof(Structured)=%u",
            snip_len, base64_sizeof(Structured));

    Structured st, *st_p=NULL;
    const char *prefixes;
    uint32_t prefixes_len;

    bool cache_exists = buf_is_allocated (&ctx->struct_cache);
    uint16_t cache_item_len;

    // if this structured exists in the cache - use the cached one
    if (cache_exists && word_index != WORD_INDEX_NONE && ((cache_item_len = *ENT (uint16_t, ctx->struct_len, word_index)))) {
        st_p = (Structured *)ENT (char, ctx->struct_cache, *ENT (uint32_t, ctx->struct_index, word_index));
        
        unsigned st_size = sizeof_structured (*st_p);
        prefixes = (char *)st_p + st_size; // prefixes are stored after the Structured
        prefixes_len = cache_item_len - st_size;
    }

    // case: not cached - decode it and optionally cache it
    if (!st_p) {
        // decode
        unsigned b64_len = snip_len;
        base64_decode (snip, &b64_len, (uint8_t*)&st);
        st.repeats = BGEN32 (st.repeats);

        // get the did_i for each dict_id
        for (uint8_t item_i=0; item_i < st.num_items; item_i++)
            if (st.items[item_i].dict_id.num)  // not a prefix-only item
                st.items[item_i].did_i = mtf_get_existing_did_i (vb, st.items[item_i].dict_id);

        // get prefixes
        unsigned st_size = sizeof_structured (st);
        st_p         = &st;
        prefixes     = (b64_len < snip_len) ? &snip[b64_len+1]       : NULL;
        prefixes_len = (b64_len < snip_len) ? snip_len - (b64_len+1) : 0;

        ASSERT (prefixes_len <= STRUCTURED_MAX_PREFIXES_LEN, "Error in piz_reconstruct_structured: prefixes_len=%u longer than STRUCTURED_MAX_PREFIXES_LEN=%u", 
                prefixes_len, STRUCTURED_MAX_PREFIXES_LEN);

        ASSERT (!prefixes_len || prefixes[prefixes_len-1] == SNIP_STRUCTURED, "Error in piz_reconstruct_structured: prefixes array does end with a SNIP_STRUCTURED: %.*s",
                prefixes_len, prefixes);

        // cache it if it is cacheable 
        if (word_index != WORD_INDEX_NONE) {

            ASSERT (st_size + prefixes_len <= 65535, "st_size=%u + prefixes_len=%u too large", st_size, prefixes_len);

            // first encounter with Structured for this context - allocate the cache
            if (!cache_exists) {
                buf_alloc (vb, &ctx->struct_index, ctx->word_list.len * sizeof (uint32_t), 1, "context->struct_index", ctx->did_i);
                buf_alloc (vb, &ctx->struct_len,   ctx->word_list.len * sizeof (uint16_t),  1, "context->struct_len",   ctx->did_i);
                buf_zero (&ctx->struct_len);
            }

            // place Structured followed by prefix in the cache
            *ENT (uint32_t, ctx->struct_index, word_index) = (uint32_t)ctx->struct_cache.len;
            *ENT (uint16_t, ctx->struct_len,   word_index) = (uint16_t)(st_size + prefixes_len);

            buf_alloc (vb, &ctx->struct_cache, ctx->struct_cache.len + st_size + prefixes_len, 2, "context->struct_cache", ctx->did_i);
            buf_add (&ctx->struct_cache, st_p, st_size);
            if (prefixes_len) buf_add (&ctx->struct_cache, prefixes, prefixes_len);
        }
    }

    piz_reconstruct_structured_do (vb, ctx->dict_id, st_p, prefixes, prefixes_len); 
}

static Context *piz_get_other_ctx_from_snip (VBlockP vb, const char **snip, unsigned *snip_len)
{
    unsigned b64_len = base64_sizeof (DictId);
    ASSERT (b64_len + 1 <= *snip_len, "Error in piz_get_other_ctx_from_snip: snip_len=%u but expecting it to be >= %u",
            *snip_len, b64_len + 1);

    DictId dict_id;
    base64_decode ((*snip)+1, &b64_len, dict_id.id);

    Context *other_ctx = mtf_get_existing_ctx (vb, dict_id);

    *snip     += b64_len + 1;
    *snip_len -= b64_len + 1;
    
    return other_ctx;
}

void piz_reconstruct_one_snip (VBlock *vb, Context *snip_ctx, 
                               WordIndex word_index, // WORD_INDEX_NONE if not used. Needed only if this snip might be a Structured that needs to be cached
                               const char *snip, unsigned snip_len)
{
    if (!snip_len) return; // nothing to do
    
    LastValueType new_value = {0};
    bool have_new_value = false;
    Context *base_ctx = snip_ctx; // this will change if the snip refers us to another data source
    bool store_uint  = ctx_is_store (snip_ctx, CTX_FL_STORE_INT);
    bool store_float = ctx_is_store (snip_ctx, CTX_FL_STORE_FLOAT);

    switch (snip[0]) {

    // display the rest of the snip first, and then the lookup up text.
    case SNIP_LOOKUP:
    case SNIP_OTHER_LOOKUP: {

        if (snip[0] == SNIP_LOOKUP) 
            { snip++; snip_len--; }
        else 
            // we are request to reconstruct from another ctx
            base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len

        // case 1: LOCAL is not LT_SEQ* - we reconstruct this snip before adding the looked up data
        if (snip_len && base_ctx->ltype != LT_SEQUENCE && base_ctx->ltype != LT_BITMAP) 
            RECONSTRUCT (snip, snip_len);
        
        if (base_ctx->ltype >= LT_INT8 && base_ctx->ltype <= LT_UINT64) {
            new_value.i = piz_reconstruct_from_local_int (vb, base_ctx, 0);
            have_new_value = true;
        }

        // case 2: LT_SEQUENCE  - the snip is taken to be the length of the sequence (or if missing, the length will be taken from vb->seq_len)
        else if (base_ctx->ltype == LT_SEQUENCE) 
            piz_reconstruct_from_local_sequence (vb, base_ctx, snip, snip_len);

        else if (base_ctx->ltype == LT_BITMAP) {
            ASSERT (DTP (reconstruct_seq), "Error: data_type=%s doesn't support reconstruct_seq", dt_name (vb->data_type));
            DTP (reconstruct_seq) (vb, base_ctx, snip, snip_len);
        }

        else piz_reconstruct_from_local_text (vb, base_ctx); // this will call us back recursively with the snip retrieved
                
        break;
    }
    case SNIP_PAIR_LOOKUP:
        mtf_get_next_snip (vb, snip_ctx, &snip_ctx->pair_b250_iter, &snip, &snip_len);
        piz_reconstruct_one_snip (vb, snip_ctx, WORD_INDEX_NONE /* we can't cache pair items */, snip, snip_len); // might include delta etc - works because in --pair, ALL the snips in a context are PAIR_LOOKUP
        break;

    case SNIP_SELF_DELTA:
        new_value.i = piz_reconstruct_from_delta (vb, snip_ctx, base_ctx, snip+1, snip_len-1);
        have_new_value = true;
        break;

    case SNIP_OTHER_DELTA: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        new_value.i = piz_reconstruct_from_delta (vb, snip_ctx, base_ctx, snip, snip_len); 
        have_new_value = true;
        break;

    case SNIP_PAIR_DELTA: { // used for FASTQ_GPOS 
        uint32_t fastq_line_i = vb->line_i / 4 - vb->first_line; 
        int64_t pair_value = (int64_t) *ENT (uint32_t, snip_ctx->pair, fastq_line_i);  
        int64_t delta = (int64_t)strtoull (snip+1, NULL, 10 /* base 10 */); 
        new_value.i = pair_value + delta;
        RECONSTRUCT_INT (new_value.i);
        have_new_value = true;
        break;
    }

    case SNIP_STRUCTURED:
        piz_reconstruct_structured (vb, snip_ctx, word_index, snip+1, snip_len-1);
        break;

    case SNIP_SPECIAL:
        ASSERT (snip_len >= 2, "Error: SNIP_SPECIAL expects snip_len >= 2. ctx=%s", snip_ctx->name);
        uint8_t special = snip[1] - 32; // +32 was added by SPECIAL macro
        ASSERT (special < DTP (num_special), "Error: file requires special handler %u which doesn't exist in this version of genounzip - please upgrade to the latest version", special);
        DTP(special)[special](vb, snip_ctx, snip+2, snip_len-2);  
        break;

    case SNIP_REDIRECTION: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        piz_reconstruct_from_ctx (vb, base_ctx->did_i, 0);
        break;
    
    case SNIP_DONT_STORE:
        store_uint = store_float = false; // override CTX_FL_STORE_* and fall through
        snip++; snip_len--;
        
    default: {
        RECONSTRUCT (snip, snip_len); // simple reconstruction

        if (store_uint) 
            // store the value only if the snip in its entirety is a reconstructable integer (eg NOT "21A", "-0", "012" etc)
            have_new_value = str_get_int (snip, snip_len, &new_value.i);

        else if (store_float) {
            char *after;
            new_value.d = strtod (snip, &after); // allows negative values

            // if the snip in its entirety is not a valid integer, don't store the value.
            // this can happen for example when seg_pos_field stores a "nonsense" snip.
            have_new_value = (after == snip + snip_len);
            
        }

        snip_ctx->last_delta = 0; // delta is 0 since we didn't calculate delta
    }
    }

    // update last_value if needed
    if (have_new_value && (store_uint || store_float)) // note: we store in our own context, NOT base (a context, eg FORMAT/DP, sometimes serves as a base_ctx of MIN_DP and sometimes as the snip_ctx for INFO_DP)
        snip_ctx->last_value = new_value;

    snip_ctx->last_line_i = vb->line_i;
}

// returns reconstructed length or -1 if snip is missing and item's operator should not be emitted
int32_t piz_reconstruct_from_ctx_do (VBlock *vb, DidIType did_i, char sep)
{
    ASSERT (did_i < vb->num_contexts, "Error in piz_reconstruct_from_ctx_do: did_i=%u out of range: vb->num_contexts=%u", did_i, vb->num_contexts);

    Context *ctx = &vb->contexts[did_i];

    ASSERT0 (ctx->dict_id.num || ctx->did_i != DID_I_NONE, "Error in piz_reconstruct_from_ctx: ctx not initialized (dict_id=0)");

    // update ctx, if its an alias (only for primary field aliases as they have contexts, other alias don't have ctx)
    if (!ctx->dict_id.num) 
        ctx = &vb->contexts[ctx->did_i]; // ctx->did_i is different than did_i if its an alias

    uint64_t start = vb->txt_data.len;

    // case: we have dictionary data
    if (ctx->b250.len) {         
        DECLARE_SNIP;
        uint32_t word_index = LOAD_SNIP(ctx->did_i); 

        if (!snip) return -1; // WORD_INDEX_MISSING_SF - remove preceding separator
        
        piz_reconstruct_one_snip (vb, ctx, word_index, snip, snip_len);        

        // handle chrom and pos to determine whether this line should be grepped-out in case of --regions
        if (did_i == CHROM) { // NOTE: CHROM cannot have aliases, because looking up the did_i by dict_id will lead to CHROM, and this code will be executed for a non-CHROM field
            vb->chrom_node_index = word_index;
            vb->chrom_name       = snip; // used for reconstruction from external reference
            vb->chrom_name_len   = snip_len;
        }

        if (flag_regions && did_i == DTF(pos) && !regions_is_site_included (vb->chrom_node_index, ctx->last_value.i)) 
            vb->dont_show_curr_line = true;
    }
    
    // case: all data is only in local
    else if (ctx->local.len) {
        if (ctx->ltype >= LT_INT8 && ctx->ltype <= LT_UINT64)
            piz_reconstruct_from_local_int(vb, ctx, 0);
        
        else if (ctx->ltype == LT_HT)
            piz_reconstruct_from_local_ht (vb, ctx);

        else if (ctx->ltype == LT_SEQUENCE) 
            piz_reconstruct_from_local_sequence (vb, ctx, NULL, 0);
        
        else if (ctx->ltype == LT_BITMAP) {
            ASSERT (DTP (reconstruct_seq), "Error: data_type=%s doesn't support reconstruct_seq", dt_name (vb->data_type));
            DTP (reconstruct_seq) (vb, ctx, NULL, 0);
        }
        
        else if (ctx->ltype == LT_TEXT)
            piz_reconstruct_from_local_text (vb, ctx);

        else if (ctx->ltype == LT_DOMQUAL)
            domqual_reconstruct (vb, ctx);

        else ABORT ("Invalid ltype=%u in ctx=%s of vb_i=%u line_i=%u", ctx->ltype, ctx->name, vb->vblock_i, vb->line_i);
    }

    // in case of LT_BITMAP, it is it is ok if the bitmap is empty and all the data is in SEQ_NOREF (e.g. unaligned SAM)
    else if (ctx->ltype == LT_BITMAP && (ctx+1)->local.len) {
        ASSERT (DTP (reconstruct_seq), "Error: data_type=%s doesn't support reconstruct_seq", dt_name (vb->data_type));
        DTP (reconstruct_seq) (vb, ctx, NULL, 0);
    }

    // case: the entire VB was just \n - so seg dropped the ctx
    else if (ctx->did_i == DTF(eol))
        RECONSTRUCT1('\n');

    else ABORT("Error in piz_reconstruct_from_ctx_do: ctx %s has no data (b250 or local) in vb_i=%u line_i=%u did_i=%u ctx->did=%u ctx->dict_id=%s", 
                ctx->name, vb->vblock_i, vb->line_i, did_i, ctx->did_i, err_dict_id (ctx->dict_id));

    if (sep) RECONSTRUCT1 (sep); 

    return (int32_t)(vb->txt_data.len - start);
}

uint32_t piz_uncompress_all_ctxs (VBlock *vb, 
                                  uint32_t pair_vb_i) // used in ZIP when uncompressing previous file's paired sections
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);
    
    uint32_t section_i = pair_vb_i ? 0 : 1; // normally, we skip the VB header, but when uncompressing paired sections there is no VB header
    while (section_i < vb->z_section_headers.len) {

        if (section_i == vb->z_section_headers.len) break; // no more sections left

        SectionHeaderCtx *header = (SectionHeaderCtx *)ENT (char, vb->z_data, section_index[section_i]);

        bool is_local = header->h.section_type == SEC_LOCAL;
        if (header->h.section_type == SEC_B250 || is_local) {
            Context *ctx = mtf_get_ctx (vb, header->dict_id); // creates the context
            ctx->flags |= header->h.flags;
            ctx->ltype = header->ltype;

            // case: in PIZ: CTX_FL_PAIRED appears on the sections the "pair 2" VB (that come first in section_index)
            if ((ctx->flags & CTX_FL_PAIRED) && !pair_vb_i) 
                pair_vb_i = fastq_get_pair_vb_i (vb);

            bool is_pair_section = (BGEN32 (header->h.vblock_i) == pair_vb_i); // is this a section of "pair 1" 

            Buffer *target_buf = is_local ? &ctx->local : &ctx->b250;

            zfile_uncompress_section (vb, header, 
                                      is_pair_section ? &ctx->pair      : target_buf, 
                                      is_pair_section ? "contexts.pair" : is_local ? "contexts.local" : "contexts.b250", 
                                      is_pair_section ? pair_vb_i : vb->vblock_i,
                                      header->h.section_type); 

            if (!is_pair_section && is_local && dict_id_printable (ctx->dict_id).num == dump_one_local_dict_id.num) 
                mtf_dump_local (ctx, true);

            if (!is_pair_section && !is_local && dict_id_printable (ctx->dict_id).num == dump_one_b250_dict_id.num) 
                mtf_dump_local (ctx, false);

#           define adjust_lens(buf) { \
                buf.len /= lt_desc[ctx->ltype].width; \
                if (ctx->ltype == LT_BITMAP) { \
                    buf.param = buf.len * 64; /* number of bits. note: this might be higher than the number of bits on the ZIP side, since we are rounding up the word boundary */ \
                    LTEN_bit_array (buf_get_bitarray (&buf), true); \
                } \
                else if (ctx->ltype >= LT_INT8 && ctx->ltype <= LT_UINT64)    \
                    lt_desc[ctx->ltype].file_to_native (&buf); \
            }

            if      (is_pair_section) adjust_lens (ctx->pair)
            else if (is_local)        adjust_lens (ctx->local);

            if (header->h.flags & CTX_FL_COPY_PARAM)
                target_buf->param = header->param;

            section_i++;
        }    

        else break;
    }

    // if all we wanted is to dump some data, we're done
    if (exe_type == EXE_GENOCAT && (dump_one_b250_dict_id.num || dump_one_local_dict_id.num)) exit_ok;

    // initialize pair iterators (pairs only exist in fastq)
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        Context *ctx = &vb->contexts[did_i];
        if (buf_is_allocated (&ctx->pair))
            ctx->pair_b250_iter = (SnipIterator){ .next_b250 = FIRSTENT (uint8_t, ctx->pair),
                                                  .prev_word_index = -1 };
    }

    mtf_map_aliases (vb);

    return section_i;
}

static void piz_uncompress_one_vb (VBlock *vb)
{
    START_TIMER;

    ASSERT0 (!flag_reference || genome.ref.num_of_bits, "Error in piz_uncompress_one_vb: reference is not loaded correctly");

    // we read the header and ctxs for all data_types
    ARRAY (const uint32_t, section_index, vb->z_section_headers); 

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);      
    vb->lines.len        = BGEN32 (header->num_lines);       
    vb->vb_data_size     = BGEN32 (header->vb_data_size);    
    vb->longest_line_len = BGEN32 (header->longest_line_len);
    vb->md5_hash_so_far  = header->md5_hash_so_far;

    // in case of --unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every vcf component
    if (flag_unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);

    if (flag_show_vblocks) 
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u txt_size=%u genozip_size=%u longest_line_len=%u\n",
                    vb->vblock_i, vb->first_line, (uint32_t)vb->lines.len, vb->vb_data_size, BGEN32 (header->z_data_bytes), vb->longest_line_len);

    buf_alloc (vb, &vb->txt_data, vb->vb_data_size + 10000, 1.1, "txt_data", vb->vblock_i); // +10000 as sometimes we pre-read control data (eg structured templates) and then roll back

    piz_uncompress_all_ctxs (vb, 0);

    // genocat --show-b250 only shows the b250 info and not the file data (shown when uncompressing via zfile_uncompress_section)
    if (flag_show_b250 && exe_type == EXE_GENOCAT) goto done;

    // reconstruct from top level snip
    DidIType toplevel_did_i = mtf_get_existing_did_i (vb, dict_id_field (dict_id_make (TOPLEVEL, strlen(TOPLEVEL))));
    piz_reconstruct_from_ctx (vb, toplevel_did_i, 0);

done:
    vb->is_processed = true; /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (vb->profile.compute);
}

static void piz_read_all_ctxs (VBlock *vb, SectionListEntry **next_sl)
{
    // ctxs that have dictionaries are already initialized, but others (eg local data only) are not
    mtf_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts);

    while ((*next_sl)->section_type == SEC_B250 || (*next_sl)->section_type == SEC_LOCAL) {
        uint32_t section_start = vb->z_data.len;
        *ENT (uint32_t, vb->z_section_headers, vb->z_section_headers.len) = section_start; 

        int32_t ret = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", sizeof(SectionHeaderCtx), 
                                          (*next_sl)->section_type, *next_sl); // returns 0 if section is skipped

        if (ret) vb->z_section_headers.len++;
        (*next_sl)++;                             
    }
}

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static DataType piz_read_global_area (Md5Hash *original_file_digest) // out
{
    DataType data_type = zfile_read_genozip_header (original_file_digest);
    
    if (flag_show_stats) stats_read_and_display();

    if (data_type == DT_NONE) return DT_NONE;

    dict_id_initialize (data_type); // must run after zfile_read_genozip_header that sets z_file->data_type
    
    // for FASTA and FASTQ we convert a "header_only" flag to "header_one" as flag_header_only has some additional logic
    // that doesn't work for FASTA / FASTQ
    if (flag_header_only && (data_type == DT_FASTA || data_type == DT_FASTQ)) {
        flag_header_only = false;
        flag_header_one  = true;
    }

    // check if the genozip file includes a reference
    bool has_ref_sections = sections_has_reference();

    ASSERT (!has_ref_sections || flag_reference != REF_EXTERNAL || flag_reading_reference, 
            "Error: cannot use --reference with %s because it was not compressed with --reference", z_name);

    if (!flag_reading_reference && has_ref_sections) 
        flag_reference = REF_STORED; // possibly override REF_EXTERNAL (it will be restored for the next file in )

    // if the user wants to see only the header, we can skip the dictionaries, regions and random access
    if (!flag_header_only) {
        
        zfile_read_all_dictionaries (0, DICTREAD_ALL); // read all CHROM/RNAME dictionaries - needed for regions_make_chregs()

        // update chrom node indices using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
        regions_make_chregs();

        // if the regions are negative, transform them to the positive complement instead
        regions_transform_negative_to_positive_complement();

        // if this is a stored reference we load the reference random access that will determined which reference sections
        // should be read & uncompressed in case of --regions.
        // note: in case of a data file with stored reference - SEC_REF_RAND_ACC will contain the random access of the reference
        // and SEC_RANDOM_ACCESS will contain the random access of the data. In case of a .ref.genozip file, both sections exist 
        // and are identical. It made the coding easier and their size is negligible.
        if (sections_has_random_access())
            random_access_load_ra_section (SEC_RANDOM_ACCESS, &z_file->ra_buf, "z_file->ra_buf", 
                                            flag_show_index ? "Random-access index contents (result of --show-index)" : NULL);

        if (has_ref_sections) 
            random_access_load_ra_section (SEC_REF_RAND_ACC, &ref_stored_ra, "ref_stored_ra", 
                                            flag_show_ref_index && !flag_reading_reference ? "Reference random-access index contents (result of --show-index)" : NULL);

        if ((flag_reference == REF_STORED || flag_reference == REF_EXTERNAL) && !flag_reading_reference)
            ref_contigs_sort_chroms(); // create alphabetically sorted index for user file chrom word list

        if (has_ref_sections) // note: in case of REF_EXTERNAL, reference is already pre-loaded
            ref_contigs_load_contigs();

        // mapping of the file's chroms to the reference chroms (for files originally compressed with REF_EXTERNAL/EXT_STORE and have alternative chroms)
        ref_alt_chroms_load();

        if (has_ref_sections) { // note: in case of REF_EXTERNAL, reference is already pre-loaded
            ref_load_stored_reference();

            // load the refhash, if we are compressing FASTA or FASTQ
            if (flag_reading_reference && primary_command == ZIP && flag_ref_use_aligner) 
                refhash_load();

            // exit now if all we wanted was just to see the reference (we've already shown it)
            if ((flag_show_reference || flag_show_is_set || flag_show_ref_hash) && exe_type == EXE_GENOCAT) exit_ok;

            if (flag_reading_reference) 
                progress_finalize_component ((flag_test && !flag_reading_reference) ? "Success" : "Done");
        }

        // read dict_id aliases, if there are any
        dict_id_read_aliases();
    }
    
    file_seek (z_file, 0, SEEK_SET, false);

    return data_type;
}

static bool piz_read_one_vb (VBlock *vb)
{
    START_TIMER; 

    SectionListEntry *sl = sections_vb_first (vb->vblock_i, false); 

    int32_t vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", sizeof (SectionHeaderVbHeader), 
                                                   SEC_VB_HEADER, sl++); 

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ 

    buf_alloc (vb, &vb->z_section_headers, (MAX_DICTS * 2 + 50) * sizeof(uint32_t), 0, "z_section_headers", 1); // room for section headers  
    NEXTENT (uint32_t, vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    // read all b250 and local of all fields and subfields
    piz_read_all_ctxs (vb, &sl);

    // read additional sections and other logic specific to this data type
    bool ok_to_compute = DTPZ(piz_read_one_vb) ? DTPZ(piz_read_one_vb)(vb, sl) : true; // true if we should go forward with computing this VB (otherwise skip it)

    COPY_TIMER (vb->profile.piz_read_one_vb); 

    return ok_to_compute;
}

// returns true is successfully outputted a txt file
bool piz_dispatcher (bool is_first_component, bool is_last_file)
{
    // static dispatcher - with flag_unbind, we use the same dispatcher when unzipping components
    static Dispatcher dispatcher = NULL;
    SectionListEntry *sl_ent = NULL;
    bool piz_successful = true;

    if (flag_unbind && !sections_has_more_components()) return false; // no more components

    // read genozip header
    Md5Hash original_file_digest;

    // read genozip header, dictionaries etc and set the data type when reading the first component of in case of --unbind, 
    static DataType data_type = DT_NONE; 

    if (is_first_component) {

        data_type = piz_read_global_area (&original_file_digest);
        if (data_type == DT_NONE || flag_reading_reference) goto finish; // reference file has no VBs

        ASSERT (sections_get_next_header_type(&sl_ent, NULL, NULL) == SEC_TXT_HEADER, "Error: unable to find TXT Header data in %s", z_name);

        ASSERT (!flag_test || !md5_is_zero (original_file_digest), 
                "Error testing %s: --test cannot be used with this file, as it was not compressed with --md5 or --test", z_name);
    }

    if (!dispatcher) 
        dispatcher = dispatcher_init (global_max_threads, 0, flag_test, is_last_file, z_file->basename, PROGRESS_PERCENT, 0);

    // case: we couldn't open the file because we didn't know what type it is - open it now
    if (!flag_unbind && !flag_reading_reference && !txt_file->file) file_open_txt (txt_file);

    // read and write txt header. in unbind mode this also opens txt_file
    piz_successful = txtfile_genozip_to_txt_header (NULL, &original_file_digest);
    
    ASSERT (piz_successful || !is_first_component, "Error: failed to read %s header in %s", dt_name (z_file->data_type), z_name);

    if (!piz_successful || flag_header_only) goto finish;

    if (flag_unbind) dispatcher_resume (dispatcher); // accept more input 

    if (DTPZ(piz_initialize)) DTPZ(piz_initialize)();

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. Wait for the first thread (by sequential order) to complete and write data

    bool header_only_file = true; // initialize
    do {
        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a variant block and compute it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher)) {

            bool still_more_data = false, grepped_out = false;
            bool skipped_vb;
            static Buffer region_ra_intersection_matrix = EMPTY_BUFFER; // we will move the data to the VB when we get it
            SectionType header_type = sections_get_next_header_type (&sl_ent, &skipped_vb, &region_ra_intersection_matrix);
            switch (header_type) {
                case SEC_VB_HEADER: {

                    // if we skipped VBs or we skipped the sample sections in the last vb (flag_drop_genotypes), we need to seek forward 
                    if (skipped_vb || flag_drop_genotypes) file_seek (z_file, sl_ent->offset, SEEK_SET, false); 

                    VBlock *next_vb = dispatcher_generate_next_vb (dispatcher, sl_ent->vblock_i);
                    
                    if (region_ra_intersection_matrix.data) {
                        buf_copy (next_vb, &next_vb->region_ra_intersection_matrix, &region_ra_intersection_matrix, 0,0,0, "region_ra_intersection_matrix", next_vb->vblock_i);
                        buf_free (&region_ra_intersection_matrix); // note: copy & free rather than move - so memory blocks are preserved for VB re-use
                    }
                    
                    // read one VB's genozip data
                    grepped_out = !piz_read_one_vb (next_vb);

                    if (grepped_out) dispatcher_abandon_next_vb (dispatcher); 

                    still_more_data = true; // not eof yet

                    break;
                }

                case SEC_TXT_HEADER: // 2nd+ txt header of a bound file
                    if (!flag_unbind) {
                        txtfile_genozip_to_txt_header (sl_ent, NULL); // skip 2nd+ txt header if unbinding
                        continue;
                    }
                    break; // eof if splitting
                
                case SEC_NONE: 
                    break; 
                
                default: ABORT ("Error in piz_dispatcher: unexpected section_type=%s", st_name (header_type));
            }
            
            if (still_more_data) {
                if (!grepped_out) 
                    dispatcher_compute (dispatcher, piz_uncompress_one_vb);
                    
                header_only_file = false;                
            }
            else { // eof
                dispatcher_input_exhausted (dispatcher);

                if (header_only_file)
                    dispatcher_finalize_one_vb (dispatcher);
            }
        }

        // PRIORITY 2: Wait for the first thread (by sequential order) to complete and write data
        else { // if (dispatcher_has_processed_vb (dispatcher, NULL)) {
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); 

            // read of a normal file - output uncompressed block (unless we're reading a reference - we don't need to output it)
            if (!flag_reading_reference) txtfile_write_one_vblock (processed_vb);

            z_file->num_vbs++;
            
            z_file->txt_data_so_far_single += processed_vb->vb_data_size; 

            dispatcher_finalize_one_vb (dispatcher);
        }

    } while (!dispatcher_is_done (dispatcher));

    // verify file integrity, if the genounzip compress was run with --md5 or --test
    Md5Hash decompressed_file_digest = MD5HASH_NONE;
    if (flag_md5) {
        decompressed_file_digest = md5_finalize (&txt_file->md5_ctx_bound); // z_file might be a bound file - this is the MD5 of the entire bounding file
        char s[200]; 

        if (md5_is_zero (original_file_digest)) { 
            sprintf (s, "MD5 = %s", md5_display (decompressed_file_digest));
            progress_finalize_component (s); 
        }

        else if (md5_is_equal (decompressed_file_digest, original_file_digest)) {

            if (flag_md5) { 
                sprintf (s, "MD5 = %s verified as identical to the original %s", 
                         md5_display (decompressed_file_digest), dt_name (z_file->data_type));
                progress_finalize_component (s); 
            }

            //else if (flag_test) progress_finalize_component ("Success");
        }
        else if (flag_test) {
            progress_finalize_component ("FAILED!!!");
            ABORT ("Error: MD5 of original file=%s is different than decompressed file=%s\nPlease contact bugs@genozip.com to help fix this bug in genozip\n",
                   md5_display (original_file_digest), md5_display (decompressed_file_digest));
        }

        else ASSERT (md5_is_zero (original_file_digest), // its ok if we decompressed only a partial file
                     "File integrity error: MD5 of decompressed file %s is %s, but the original %s file's was %s", 
                     txt_file->name, md5_display (decompressed_file_digest), dt_name (txt_file->data_type), 
                     md5_display (original_file_digest));
    }

    if (flag_unbind) file_close (&txt_file, true); // close this component file

    if (!flag_test) 
        progress_finalize_component_time ("Done", decompressed_file_digest);

finish:
    // in unbind mode - we continue with the same dispatcher in the next component. otherwise, we finish with it here
    if (!flag_unbind || !piz_successful) {
        if (dispatcher) dispatcher_finish (&dispatcher, NULL);
    } else
        dispatcher_pause (dispatcher);

    return piz_successful;
}
