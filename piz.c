
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
#include "profiler.h"
#include "stats.h"
#include "codec.h"
#include "container.h"
#include "bgzf.h"
#include "flags.h"

// Compute threads: decode the delta-encoded value of the POS field, and returns the new lacon_pos
// Special values:
// "-" - negated previous value
// ""  - negated previous delta
static int64_t piz_reconstruct_from_delta (VBlock *vb, 
                                           Context *my_ctx,   // use and store last_delta
                                           Context *base_ctx, // get last_value
                                           const char *delta_snip, unsigned delta_snip_len,
                                           bool reconstruct) 
{
    ASSERT (delta_snip, "Error in piz_reconstruct_from_delta: delta_snip is NULL. vb_i=%u", vb->vblock_i);
    ASSERT (base_ctx->flags.store == STORE_INT, "Error in piz_reconstruct_from_delta: attempting calculate delta from a base of \"%s\", but this context doesn't have STORE_INT",
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
    ASSERT (ctx->next_local < ctx->local.len, \
            "Error in %s:%u reconstructing txt_line=%u vb_i=%u: unexpected end of ctx->local data in %s (len=%u ltype=%s lcodec=%s)", \
            __FUNCTION__, __LINE__, vb->line_i, vb->vblock_i, ctx->name, (uint32_t)ctx->local.len, lt_desc[ctx->ltype].name, codec_name (ctx->lcodec));

static uint32_t piz_reconstruct_from_local_text (VBlock *vb, Context *ctx, bool reconstruct)
{
    uint32_t start = ctx->next_local; 
    ARRAY (char, data, ctx->local);

    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != SNIP_SEP) ctx->next_local++;
    ASSERT_IN_BOUNDS;

    char *snip = &data[start];
    uint32_t snip_len = ctx->next_local - start; 
    ctx->next_local++; /* skip the tab */ 

    piz_reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, reconstruct);

    return snip_len;
}

static int64_t piz_reconstruct_from_local_int (VBlock *vb, Context *ctx, char seperator /* 0 if none */, bool reconstruct)
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
    if (reconstruct) { 
        RECONSTRUCT_INT (num);
        if (seperator) RECONSTRUCT1 (seperator);
    }

    return num;
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

    // special case: in SAM, sam_zip_qual re-wrote a '*' marking 'unavailable' as ' ' to avoid confusing with '*' as a valid quality score
    if (ctx->local.data[ctx->next_local] == ' ') {
        len = 1;
        if (reconstruct) RECONSTRUCT1 ('*');
    }
    else {
        len = vb->seq_len;
        ASSERT (ctx->next_local + len <= ctx->local.len, "Error in piz_reconstruct_from_local_sequence: reading txt_line=%u vb_i=%u: unexpected end of %s data", 
                vb->line_i, vb->vblock_i, ctx->name);

        if (reconstruct) RECONSTRUCT (&ctx->local.data[ctx->next_local], len);
    }

    ctx->last_value.i = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
    ctx->next_local += len;
}

static Context *piz_get_other_ctx_from_snip (VBlockP vb, const char **snip, unsigned *snip_len)
{
    unsigned b64_len = base64_sizeof (DictId);
    ASSERT (b64_len + 1 <= *snip_len, "Error in piz_get_other_ctx_from_snip: snip_len=%u but expecting it to be >= %u",
            *snip_len, b64_len + 1);

    DictId dict_id;
    base64_decode ((*snip)+1, &b64_len, dict_id.id);

    Context *other_ctx = ctx_get_existing_ctx (vb, dict_id);

    *snip     += b64_len + 1;
    *snip_len -= b64_len + 1;
    
    return other_ctx;
}

void piz_reconstruct_one_snip (VBlock *vb, Context *snip_ctx, 
                               WordIndex word_index, // WORD_INDEX_NONE if not used.
                               const char *snip, unsigned snip_len,
                               bool reconstruct) // if false, calculates last_value but doesn't output to vb->txt_data)
{
    if (!snip_len) return; // nothing to do
    
    LastValueType new_value = {0};
    bool have_new_value = false;
    Context *base_ctx = snip_ctx; // this will change if the snip refers us to another data source
    enum StoreType store_type = snip_ctx->flags.store;

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
        
        if (base_ctx->ltype >= LT_INT8 && base_ctx->ltype <= LT_UINT64) {
            new_value.i = piz_reconstruct_from_local_int (vb, base_ctx, 0, reconstruct);
            have_new_value = true;
        }

        // case: the snip is taken to be the length of the sequence (or if missing, the length will be taken from vb->seq_len)
        else if (base_ctx->ltype == LT_SEQUENCE) 
            piz_reconstruct_from_local_sequence (vb, base_ctx, snip, snip_len);

        else if (base_ctx->ltype == LT_BITMAP) {
            ASSERT_DT_FUNC (vb, reconstruct_seq);
            DT_FUNC (vb, reconstruct_seq) (vb, base_ctx, snip, snip_len);
        }
        else piz_reconstruct_from_local_text (vb, base_ctx, reconstruct); // this will call us back recursively with the snip retrieved
                
        break;
    }
    case SNIP_PAIR_LOOKUP:
        ctx_get_next_snip (vb, snip_ctx, snip_ctx->pair_flags.all_the_same, &snip_ctx->pair_b250_iter, &snip, &snip_len);
        piz_reconstruct_one_snip (vb, snip_ctx, WORD_INDEX_NONE /* we can't cache pair items */, snip, snip_len, reconstruct); // might include delta etc - works because in --pair, ALL the snips in a context are PAIR_LOOKUP
        break;

    case SNIP_SELF_DELTA:
        new_value.i = piz_reconstruct_from_delta (vb, snip_ctx, base_ctx, snip+1, snip_len-1, reconstruct);
        have_new_value = true;
        break;

    case SNIP_OTHER_DELTA: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        new_value.i = piz_reconstruct_from_delta (vb, snip_ctx, base_ctx, snip, snip_len, reconstruct); 
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
        container_reconstruct (vb, snip_ctx, word_index, snip+1, snip_len-1);
        break;

    case SNIP_SPECIAL:
        ASSERT (snip_len >= 2, "Error: SNIP_SPECIAL expects snip_len >= 2. ctx=%s", snip_ctx->name);
        uint8_t special = snip[1] - 32; // +32 was added by SPECIAL macro

        ASSERT (special < DTP (num_special), "Error: file requires special handler %u which doesn't exist in this version of genounzip - please upgrade to the latest version", special);
        ASSERT_DT_FUNC (vb, special);

        have_new_value = DT_FUNC(vb, special)[special](vb, snip_ctx, snip+2, snip_len-2, &new_value, reconstruct);  
        break;

    case SNIP_REDIRECTION: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        piz_reconstruct_from_ctx (vb, base_ctx->did_i, 0, reconstruct);
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
                new_value.d = strtod (snip, &after); // allows negative values

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
int32_t piz_reconstruct_from_ctx_do (VBlock *vb, DidIType did_i, 
                                     char sep, // if non-zero, outputs after the reconstructino
                                     bool reconstruct) // if false, calculates last_value but doesn't output to vb->txt_data
{
    ASSERT (did_i < vb->num_contexts, "Error in piz_reconstruct_from_ctx_do: did_i=%u out of range: vb->num_contexts=%u", did_i, vb->num_contexts);

    Context *ctx = &vb->contexts[did_i];

    ASSERT0 (ctx->dict_id.num || ctx->did_i != DID_I_NONE, "Error in piz_reconstruct_from_ctx: ctx not initialized (dict_id=0)");

    // update ctx, if its an alias (only for primary field aliases as they have contexts, other alias don't have ctx)
    if (!ctx->dict_id.num) 
        ctx = &vb->contexts[ctx->did_i]; // ctx->did_i is different than did_i if its an alias

    uint64_t start = vb->txt_data.len;

    // case: we have b250 data
    if (ctx->b250.len ||
        (!ctx->b250.len && !ctx->local.len && ctx->dict.len)) {          
        DECLARE_SNIP;
        uint32_t word_index = LOAD_SNIP(ctx->did_i); // if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (!snip) return -1; // WORD_INDEX_MISSING_SF - remove preceding separator
        
        piz_reconstruct_one_snip (vb, ctx, word_index, snip, snip_len, reconstruct);        

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
        case LT_INT8 ... LT_UINT64 :
            piz_reconstruct_from_local_int(vb, ctx, 0, reconstruct); break;
        
        case LT_CODEC:
            codec_args[ctx->lcodec].reconstruct (vb, ctx->lcodec, ctx); break;

        case LT_SEQUENCE: 
            piz_reconstruct_from_local_sequence (vb, ctx, NULL, 0); break;
                
        case LT_BITMAP:
            ASSERT_DT_FUNC (vb, reconstruct_seq);
            DT_FUNC (vb, reconstruct_seq) (vb, ctx, NULL, 0);
            break;
        
        case LT_TEXT:
            piz_reconstruct_from_local_text (vb, ctx, reconstruct); break;

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

    else ABORT("Error in piz_reconstruct_from_ctx_do: ctx %s has no data (dict, b250 or local) in vb_i=%u line_i=%u did_i=%u ctx->did=%u ctx->dict_id=%s", 
                ctx->name, vb->vblock_i, vb->line_i, did_i, ctx->did_i, dis_dict_id (ctx->dict_id).s);

    if (sep && reconstruct) RECONSTRUCT1 (sep); 

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
            Context *ctx = ctx_get_ctx (vb, header->dict_id); // creates the context

            // case: in ZIP & we are PAIR_2, reading PAIR_1 info - save flags, ltype, lcodec to restore later
            struct FlagsCtx save_flags={}; LocalType save_ltype=0; Codec save_lcodec=0;
            if (pair_vb_i && command == ZIP) { save_flags=ctx->flags; save_ltype=ctx->ltype; save_lcodec=ctx->lcodec; }
            
            *(uint8_t*)&ctx->flags |= header->h.flags.flags; // ??? is |= really needed? or is assignment sufficient?

            ctx->ltype  = header->ltype;
            if (is_local) ctx->lcodec = header->h.codec;

            // case: in PIZ: flags.paired appears on the sections the "pair 2" VB (that come first in section_index)
            if (ctx->flags.paired && !pair_vb_i) 
                pair_vb_i = fastq_get_pair_vb_i (vb);

            bool is_pair_section = (BGEN32 (header->h.vblock_i) == pair_vb_i); // is this a section of "pair 1" 
            
            if (is_pair_section) ctx->pair_flags = header->h.flags.ctx;            

            Buffer *target_buf = is_local ? &ctx->local : &ctx->b250;

            zfile_uncompress_section (vb, header, 
                                      is_pair_section ? &ctx->pair      : target_buf, 
                                      is_pair_section ? "contexts.pair" : is_local ? "contexts.local" : "contexts.b250", 
                                      is_pair_section ? pair_vb_i : vb->vblock_i,
                                      header->h.section_type); 

            if (!is_pair_section && is_local && dict_id_printable (ctx->dict_id).num == flag.dump_one_local_dict_id.num) 
                ctx_dump_binary (vb, ctx, true);

            if (!is_pair_section && !is_local && dict_id_printable (ctx->dict_id).num == flag.dump_one_b250_dict_id.num) 
                ctx_dump_binary (vb, ctx, false);

#           define adjust_lens(buf) { \
                buf.len /= lt_desc[ctx->ltype].width; \
                if (ctx->ltype == LT_BITMAP) { \
                    buf.param = buf.len * 64; /* number of bits. note: this might be higher than the number of bits on the ZIP side, since we are rounding up the word boundary */ \
                    LTEN_bit_array (buf_get_bitarray (&buf)); \
                } \
                else if (ctx->ltype >= LT_INT8 && ctx->ltype <= LT_UINT64)    \
                    lt_desc[ctx->ltype].file_to_native (&buf); \
            }

            if      (is_pair_section) adjust_lens (ctx->pair)
            else if (is_local)        adjust_lens (ctx->local);

            if (header->h.flags.ctx.copy_param)
                target_buf->param = header->param;

            // restore
            if (pair_vb_i && command == ZIP) { ctx->flags=save_flags; ctx->ltype=save_ltype; ctx->lcodec=save_lcodec; }

            section_i++;
        }    

        else break;
    }

    // initialize pair iterators (pairs only exist in fastq)
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        Context *ctx = &vb->contexts[did_i];
        if (buf_is_allocated (&ctx->pair))
            ctx->pair_b250_iter = (SnipIterator){ .next_b250 = FIRSTENT (uint8_t, ctx->pair),
                                                  .prev_word_index = -1 };
    }

    ctx_map_aliases (vb);

    return section_i;
}

static void piz_uncompress_one_vb (VBlock *vb)
{
    START_TIMER;

    ASSERT0 (!flag.reference || genome.ref.num_of_bits, "Error in piz_uncompress_one_vb: reference is not loaded correctly");

    // we read the header and ctxs for all data_types
    ARRAY (const uint32_t, section_index, vb->z_section_headers); 

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);      
    vb->lines.len        = BGEN32 (header->num_lines);       
    vb->longest_line_len = BGEN32 (header->longest_line_len);
    vb->digest_so_far  = header->digest_so_far;

    // in case of --unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every vcf component
    if (flag.unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);

    if (flag.show_vblocks) 
        fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u txt_size=%u genozip_size=%u longest_line_len=%u\n",
                    vb->vblock_i, vb->first_line, (uint32_t)vb->lines.len, vb->vb_data_size, BGEN32 (header->z_data_bytes), vb->longest_line_len);

    DtTranslation trans = dt_get_translation(); // in case we're translating from one data type to another

    // note: txt_data is fully allocated in advance and cannot be extended mid-reconstruction (container_reconstruct_do and possibly others rely on this)
    buf_alloc (vb, &vb->txt_data, vb->vb_data_size * trans.factor + 10000, 1.1, "txt_data"); // +10000 as sometimes we pre-read control data (eg container templates) and then roll back

    piz_uncompress_all_ctxs (vb, 0);

    // genocat flags that with which we needn't reconstruct
    if (exe_type == EXE_GENOCAT && flag.genocat_info_only) goto done;

    // reconstruct from top level snip
    piz_reconstruct_from_ctx (vb, trans.toplevel, 0, true);

    // compress txt_data into BGZF blocks (in vb->compressed) if applicable
    if (vb->bgzf_blocks.len) bgzf_compress_vb (vb);

    // calculate the MD5 contribution of this VB to the single file and bound files, and the MD5 snapshot of this VB
    if (!digest_is_zero (vb->digest_so_far) && flag.reconstruct_as_src) 
        digest_one_vb (vb); 

done:
    vb->is_processed = true; /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (compute);
}

static void piz_read_all_ctxs (VBlock *vb, const SectionListEntry **next_sl)
{
    // ctxs that have dictionaries are already initialized, but others (eg local data only) are not
    ctx_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts);

    while ((*next_sl)->section_type == SEC_B250 || (*next_sl)->section_type == SEC_LOCAL) {
        uint32_t section_start = vb->z_data.len;
        *ENT (uint32_t, vb->z_section_headers, vb->z_section_headers.len) = section_start; 

        int32_t ret = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", (*next_sl)->section_type, *next_sl); // returns 0 if section is skipped

        if (ret) vb->z_section_headers.len++;
        (*next_sl)++;                             
    }
}

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static DataType piz_read_global_area (Digest *original_file_digest) // out
{
    bool success = zfile_read_genozip_header (original_file_digest, 0, 0, 0);
    
    if (flag.show_stats) stats_read_and_display();

    if (!success) return DT_NONE;

    dict_id_initialize (z_file->data_type); // must run after zfile_read_genozip_header that sets z_file->data_type
    
    // for FASTA and FASTQ we convert a "header_only" flag to "header_one" as flag.header_only has some additional logic
    // that doesn't work for FASTA / FASTQ
    if (flag.header_only && (z_file->data_type == DT_FASTA || z_file->data_type == DT_FASTQ)) {
        flag.header_only = false;
        flag.header_one  = true;
    }

    // check if the genozip file includes a reference
    bool has_ref_sections = !!sections_get_first_section_of_type (SEC_REFERENCE, true);

    ASSERT (!has_ref_sections || flag.reference != REF_EXTERNAL || flag.reading_reference, 
            "Error: cannot use --reference with %s because it was not compressed with --reference", z_name);

    if (!flag.reading_reference && has_ref_sections) 
        flag.reference = REF_STORED; // possibly override REF_EXTERNAL (it will be restored for the next file in )

    // if the user wants to see only the header, we can skip the dictionaries, regions and random access
    if (!flag.header_only) {
        
        ctx_read_all_dictionaries (DICTREAD_ALL); // read all CHROM/RNAME dictionaries - needed for regions_make_chregs()

        // update chrom node indices using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
        regions_make_chregs();

        // if the regions are negative, transform them to the positive complement instead
        regions_transform_negative_to_positive_complement();

        // if this is a stored reference we load the reference random access that will determined which reference sections
        // should be read & uncompressed in case of --regions.
        // note: in case of a data file with stored reference - SEC_REF_RAND_ACC will contain the random access of the reference
        // and SEC_RANDOM_ACCESS will contain the random access of the data. In case of a .ref.genozip file, both sections exist 
        // and are identical. It made the coding easier and their size is negligible.
        random_access_load_ra_section (SEC_RANDOM_ACCESS, &z_file->ra_buf, "z_file->ra_buf", 
                                        flag.show_index ? "Random-access index contents (result of --show-index)" : NULL);

        random_access_load_ra_section (SEC_REF_RAND_ACC, &ref_stored_ra, "ref_stored_ra", 
                                        flag.show_ref_index && !flag.reading_reference ? "Reference random-access index contents (result of --show-index)" : NULL);

        if ((flag.reference == REF_STORED || flag.reference == REF_EXTERNAL) && 
            !flag.reading_reference &&
            !(flag.show_headers && exe_type == EXE_GENOCAT))
            ref_contigs_sort_chroms(); // create alphabetically sorted index for user file chrom word list

        ref_contigs_load_contigs(); // note: in case of REF_EXTERNAL, reference is already pre-loaded

        // mapping of the file's chroms to the reference chroms (for files originally compressed with REF_EXTERNAL/EXT_STORE and have alternative chroms)
        ref_alt_chroms_load();

        if (has_ref_sections) { // note: in case of REF_EXTERNAL, reference is already pre-loaded
            ref_load_stored_reference();

            // load the refhash, if we are compressing FASTA or FASTQ
            if (flag.reading_reference && primary_command == ZIP && flag.ref_use_aligner) 
                refhash_load();

            // exit now if all we wanted was just to see the reference (we've already shown it)
            if ((flag.show_reference || flag.show_is_set || flag.show_ref_hash) && exe_type == EXE_GENOCAT) exit_ok;

            if (flag.reading_reference) 
                progress_finalize_component ((flag.test && !flag.reading_reference) ? "Success" : "Done");
        }

        // read dict_id aliases, if there are any
        dict_id_read_aliases();
    }
    
    file_seek (z_file, 0, SEEK_SET, false);

    return z_file->data_type;
}

static bool piz_read_one_vb (VBlock *vb)
{
    START_TIMER; 

    const SectionListEntry *sl = sections_vb_first (vb->vblock_i, false); 

    vb->vb_position_txt_file = txt_file->txt_data_so_far_single;
    
    int32_t vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, &vb->z_data, "z_data", SEC_VB_HEADER, sl++); 
    vb->vb_data_size = BGEN32 (((SectionHeaderVbHeader *)vb->z_data.data)->vb_data_size); // needed by bgzf_calculate_blocks_one_vb

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);
    ctx_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ 

    buf_alloc (vb, &vb->z_section_headers, (MAX_DICTS * 2 + 50) * sizeof(uint32_t), 0, "z_section_headers"); // room for section headers  

    NEXTENT (uint32_t, vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    // read all b250 and local of all fields and subfields
    piz_read_all_ctxs (vb, &sl);

    // read additional sections and other logic specific to this data type
    bool ok_to_compute = DTPZ(piz_read_one_vb) ? DTPZ(piz_read_one_vb)(vb, sl) : true; // true if we should go forward with computing this VB (otherwise skip it)

    // calculate the BGZF blocks that the compute thread is expected to compress
    if (flag.bgzf) bgzf_calculate_blocks_one_vb (vb, vb->vb_data_size);

    COPY_TIMER (piz_read_one_vb); 

    return ok_to_compute;
}

static Digest piz_one_file_verify_md5 (Digest original_file_digest)
{
    if (digest_is_zero (original_file_digest) || flag.genocat_info_only) return DIGEST_NONE; // file was not compressed with --md5 or --test

    Digest decompressed_file_digest = digest_finalize (&txt_file->digest_ctx_bound); // z_file might be a bound file - this is the MD5 of the entire bound file
    char s[200]; 

    if (digest_is_zero (original_file_digest)) { 
        sprintf (s, "%s = %s", digest_name(), digest_display (decompressed_file_digest).s);
        progress_finalize_component (s); 
    }

    else if (digest_is_equal (decompressed_file_digest, original_file_digest)) {

        if (flag.test) { 
            sprintf (s, "%s = %s verified as identical to the original %s", 
                        digest_name(), digest_display (decompressed_file_digest).s, dt_name (txt_file->data_type));
            progress_finalize_component (s); 
        }
    }
    
    else if (flag.test) {
        progress_finalize_component ("FAILED!!!");
        ABORT ("Error: MD5 of original file=%s is different than decompressed file=%s\nPlease contact bugs@genozip.com to help fix this bug in genozip\n",
                digest_display (original_file_digest).s, digest_display (decompressed_file_digest).s);
    }

    else ASSERT (digest_is_zero (original_file_digest), // its ok if we decompressed only a partial file
                 "File integrity error: MD5 of decompressed file %s is %s, but the original %s file's was %s", 
                 txt_file->name, digest_display (decompressed_file_digest).s, dt_name (txt_file->data_type), 
                 digest_display (original_file_digest).s);

    return decompressed_file_digest;
}

// returns false if VB was dispatched, and true if vb was skipped
bool piz_dispatch_one_vb (Dispatcher dispatcher, const SectionListEntry *sl_ent)
{
    static Buffer region_ra_intersection_matrix = EMPTY_BUFFER; // we will move the data to the VB when we get it
    if (!random_access_is_vb_included (sl_ent->vblock_i, &region_ra_intersection_matrix)) return true; // skip this VB if not included

    VBlock *next_vb = dispatcher_generate_next_vb (dispatcher, sl_ent->vblock_i);
    
    if (region_ra_intersection_matrix.data) {
        buf_copy (next_vb, &next_vb->region_ra_intersection_matrix, &region_ra_intersection_matrix, 0,0,0, "region_ra_intersection_matrix");
        buf_free (&region_ra_intersection_matrix); // note: copy & free rather than move - so memory blocks are preserved for VB re-use
    }
    
    // read one VB's genozip data
    bool grepped_out = !piz_read_one_vb (next_vb);

    if (grepped_out || (flag.show_headers && exe_type == EXE_GENOCAT)) 
        dispatcher_abandon_next_vb (dispatcher); 
    else
        dispatcher_compute (dispatcher, piz_uncompress_one_vb);

    return false;
}

// called once per components reconstructed txt_file: i.e.. if unbinding there will be multiple calls.
// returns: false if we're unbinding and there are no more components. true otherwise.
bool piz_one_file (uint32_t unbind_component_i /* 0 if not unbinding */, bool is_last_file)
{
    static Dispatcher dispatcher = NULL; // static dispatcher - with flag.unbind, we use the same dispatcher when unzipping components
    static const SectionListEntry *sl_ent = NULL; // preserve for unbinding multiple files

    // read genozip header
    Digest original_file_digest;

    // read genozip header, dictionaries etc and set the data type when reading the first component of in case of --unbind, 
    static DataType data_type = DT_NONE; 
    bool no_more_headers = false;

    if (unbind_component_i == 0) {

        data_type = piz_read_global_area (&original_file_digest);
        if (data_type == DT_NONE || flag.reading_reference) {
            no_more_headers = true; // reference file has no VBs
            goto finish; 
        }

        sl_ent = NULL; // reset

        ASSERT (!flag.test || !digest_is_zero (original_file_digest), 
                "Error testing %s: --test cannot be used with this file, as it was not compressed with --md5 or --test", z_name);

        if (flag.test || flag.md5) 
            ASSERT0 (dt_get_translation().is_src_dt, "Error: --test or --md5 cannot be used when converting a file to another format"); 

        dispatcher = dispatcher_init (global_max_threads, 0, flag.test, is_last_file, z_file->basename, PROGRESS_PERCENT, 0);
    }
    
    else if (flag.unbind) {
        if (!sections_get_next_section_of_type (&sl_ent, SEC_TXT_HEADER, false, true)) return false; // unbinding - no more components
        sl_ent--; // rewind;
    
        dispatcher_set_input_exhausted (dispatcher, false); // accept more input 
    }
  
    if (DTPZ(piz_initialize)) DTPZ(piz_initialize)();

    // this is the dispatcher loop. In each iteration, it can do one of 3 things, in this order of priority:
    // 1. In input is not exhausted, and a compute thread is available - read a variant block and compute it
    // 2. Wait for the first thread (by sequential order) to complete and write data

    bool header_only_file = true; // initialize - true until we encounter a VB header
    unsigned txt_header_i=0;

    while (!dispatcher_is_done (dispatcher)) {

        // PRIORITY 1: In input is not exhausted, and a compute thread is available - read a variant block and compute it
        if (!dispatcher_is_input_exhausted (dispatcher) && dispatcher_has_free_thread (dispatcher)) {

            bool another_header = sections_get_next_section_of_type2 (&sl_ent, SEC_TXT_HEADER, SEC_VB_HEADER, false, true);

            if (another_header && sl_ent->section_type == SEC_VB_HEADER) 
                header_only_file &= piz_dispatch_one_vb (dispatcher, sl_ent);  // function returns true if VB was skipped

            else if (another_header && sl_ent->section_type == SEC_TXT_HEADER) { // 1st component or 2nd+ component and we're concatenating

                if (flag.unbind && txt_header_i++) goto no_more_data; // if unbinding, we processes only one component (piz_one_file will be called again for the next component)

                txtfile_genozip_to_txt_header (sl_ent, unbind_component_i, (txt_header_i>=2 ? NULL : &original_file_digest)); // read header - it will skip 2nd+ txt header if concatenating

                if (flag.unbind) dispatcher_resume (dispatcher); 

                if (flag.header_only) goto finish;
            }

            else { // we're done (concatenating: no more VBs in the entire file ; unbinding: no more VBs in our component)
                no_more_headers = true;
no_more_data:            
                sl_ent--; // re-read in next call to this function if unbinding
                dispatcher_set_input_exhausted (dispatcher, true);

                if (header_only_file)
                    dispatcher_finalize_one_vb (dispatcher);
            }
        }

        // PRIORITY 2: Wait for the first thread (by sequential order) to complete and write data
        else { // if (dispatcher_has_processed_vb (dispatcher, NULL)) {
            VBlock *processed_vb = dispatcher_get_processed_vb (dispatcher, NULL); 

            // read of a normal file - output uncompressed block (unless we're reading a reference - we don't need to output it)
            if (!flag.reading_reference) txtfile_write_one_vblock (processed_vb);

            z_file->num_vbs++;
            z_file->txt_data_so_far_single += processed_vb->vb_data_size; 

            dispatcher_finalize_one_vb (dispatcher);
        }
    }

    // verifies bgzf reconstructed file against codec_args    

    // verifies reconstructed file against MD5 (if compressed with --md5 or --test) and/or codec_args (if bgzf)
    Digest decompressed_file_digest = piz_one_file_verify_md5 (original_file_digest);

    if (flag.unbind) file_close (&txt_file, true); // close this component file

    if (!flag.test) progress_finalize_component_time ("Done", decompressed_file_digest);

finish:
    // case: we're done with reconstructing this z_file - either concatenated, or this was the last component unbound
    if (no_more_headers) {
        if (dispatcher) dispatcher_finish (&dispatcher, NULL);
    } 
    
    // case: we're unbinding and still have more components - we continue with the same dispatcher in the next component.
    else
        dispatcher_pause (dispatcher);

    DT_FUNC (z_file, piz_finalize)();

    return true;
}
