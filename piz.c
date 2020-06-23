
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
#include "fasta.h"

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
    
    if (delta_snip_len == 1 && delta_snip[0] == '-')
        my_ctx->last_delta = -2 * base_ctx->last_value; // negated previous value

    else if (!delta_snip_len)
        my_ctx->last_delta = -my_ctx->last_delta; // negated previous delta

    else 
        my_ctx->last_delta = (int64_t)strtoull (delta_snip, NULL, 10 /* base 10 */); // strtoull can handle negative numbers, despite its name

    int64_t new_value = base_ctx->last_value + my_ctx->last_delta;  
    RECONSTRUCT_INT (new_value);

    return new_value;
}

static uint32_t piz_reconstruct_from_local_text (VBlock *vb, Context *ctx)
{
    uint32_t start = ctx->next_local; 
    ARRAY (char, data, ctx->local);

    while (ctx->next_local < ctx->local.len && data[ctx->next_local] != SNIP_SEP) ctx->next_local++;
    ASSERT (ctx->next_local < ctx->local.len, 
            "Error reconstructing txt_line=%u: unexpected end of ctx->local data in %s (len=%u)", vb->line_i, ctx->name, (uint32_t)ctx->local.len);

    char *snip = &data[start];
    uint32_t snip_len = ctx->next_local - start; 
    ctx->next_local++; /* skip the tab */ 

    piz_reconstruct_one_snip (vb, ctx, snip, snip_len);

    return snip_len;
}

// for signed numbers, we store them in our "interlaced" format rather than standard ISO format 
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// compressing better in an array that contains both positive and negative
// NOT TESTED YET
#define DEINTERLACE(signedtype,unum) (((unum) & 1) ? -(signedtype)(((unum)>>1)+1) : (signedtype)((unum)>>1))

static int64_t piz_reconstruct_from_local_int (VBlock *vb, Context *ctx, char seperator /* 0 if none */)
{
    unsigned width = ctx_lt_sizeof_one[ctx->ltype];
    bool is_signed = ctx_lt_is_signed [ctx->ltype];

    ASSERT (ctx->next_local < ctx->local.len,  // len is in units of width
            "Error in piz_reconstruct_from_local_int while reconstructing txt_line=%u: unexpected end of %s data (ctx->local.len=%u next=%u)", 
            vb->line_i, ctx->name, (uint32_t)ctx->local.len, ctx->next_local); 

    int64_t num=0;
    if (width == 4) { // check 4 first, as its the most popular
        uint32_t num_big_en = *ENT (uint32_t, ctx->local, ctx->next_local++);
        uint32_t unum = BGEN32 (num_big_en); 
        num = (int64_t)(is_signed ? DEINTERLACE(int32_t,unum) : unum);
    }
    else if (width == 2) {
        uint16_t num_big_en = *ENT (uint16_t, ctx->local, ctx->next_local++);
        uint16_t unum = BGEN16 (num_big_en); 
        num = (int64_t)(is_signed ? DEINTERLACE(int16_t,unum) : unum);
    }
    else if (width == 1) {
        uint8_t unum = *ENT (uint8_t, ctx->local, ctx->next_local++);
        num = (int64_t)(is_signed ? DEINTERLACE(int8_t,unum) : unum);
    }
    else if (width == 8) { // note: for uint64_t the function returns the number correctly, it just needs to be casted to uint64_t
        uint64_t num_big_en = *ENT (uint64_t, ctx->local, ctx->next_local++);
        uint64_t unum = BGEN64 (num_big_en); 
        num = (int64_t)(is_signed ? DEINTERLACE(int64_t,unum) : unum);
    }

    // TO DO: RECONSTRUCT_INT won't reconstruct large uint64_t correctly
    RECONSTRUCT_INT (num);
    if (seperator) RECONSTRUCT1 (seperator);

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

    // special case: it is "*" that was written to " " - we reconstruct it
    if (ctx->local.data[ctx->next_local] == ' ') {
        len = 1;
        if (reconstruct) RECONSTRUCT1 ('*');
    }
    else {
        len = vb->seq_len;
        ASSERT (ctx->next_local + len <= ctx->local.len, "Error reading txt_line=%u: unexpected end of %s data", vb->line_i, ctx->name);

        if (reconstruct) RECONSTRUCT (&ctx->local.data[ctx->next_local], len);
    }

    ctx->last_value = ctx->next_local; // for seq_qual, we use last_value for storing the beginning of the sequence
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

void piz_reconstruct_structured_do (VBlock *vb, const Structured *st, const char *prefixes, uint32_t prefixes_len)
{
    ASSERT (prefixes_len <= STRUCTURED_MAX_PREFIXES_LEN, "Error in piz_reconstruct_structured_do: prefixes_len=%u longer than STRUCTURED_MAX_PREFIXES_LEN=%u", 
            prefixes_len, STRUCTURED_MAX_PREFIXES_LEN);

    ASSERT (!prefixes || prefixes[prefixes_len-1] == SNIP_STRUCTURED, "Error in piz_reconstruct_structured_do: prefixes array does end with a SNIP_STRUCTURED: %.*s",
            prefixes_len, prefixes);

    // structured wide prefix - it will be missing if Structured has no prefixes, or empty if it has only items prefixes
    piz_reconstruct_structured_prefix (vb, &prefixes, &prefixes_len); // item prefix (we will have one per item or none at all)

    for (uint32_t rep_i=0; rep_i < st->repeats; rep_i++) {

        const char *item_prefixes = prefixes; // the remaining after extracting the first prefix - either one per item or none at all
        uint32_t item_prefixes_len = prefixes_len;

        for (unsigned i=0; i < st->num_items; i++) {

            piz_reconstruct_structured_prefix (vb, &item_prefixes, &item_prefixes_len); // item prefix

            const StructuredItem *item = &st->items[i];
            if (item->dict_id.num) // not a prefix-only item
                piz_reconstruct_from_ctx (vb, mtf_get_ctx (vb, item->dict_id)->did_i, 0);
            
            if (item->seperator[0]) RECONSTRUCT1 (item->seperator[0]);
            if (item->seperator[1]) RECONSTRUCT1 (item->seperator[1]);
        }

        if (st->repsep[0]) RECONSTRUCT1 (st->repsep[0]);
        if (st->repsep[1]) RECONSTRUCT1 (st->repsep[1]);
    }

    if ((st->flags & STRUCTURED_DROP_LAST_SEP_OF_LAST_ELEMENT))
        vb->txt_data.len -= (st->items[st->num_items-1].seperator[0] != 0) + 
                            (st->items[st->num_items-1].seperator[1] != 0);
}

static void piz_reconstruct_structured (VBlock *vb, const char *snip, unsigned snip_len)
{
    ASSERT (snip_len <= base64_sizeof(Structured), "Error in piz_reconstruct_structured: snip_len=%u exceed base64_sizeof(Structured)=%u",
            snip_len, base64_sizeof(Structured));

    Structured st;

    unsigned b64_len = snip_len;
    base64_decode (snip, &b64_len, (uint8_t*)&st);
    st.repeats = BGEN32 (st.repeats);
    bool has_prefixes = (b64_len < snip_len);

    piz_reconstruct_structured_do (vb, &st, 
                                   has_prefixes ? &snip[b64_len+1] : NULL,
                                   has_prefixes ? snip_len - (b64_len+1) : 0);
}

static Context *piz_get_other_ctx_from_snip (VBlockP vb, const char **snip, unsigned *snip_len)
{
    unsigned b64_len = base64_sizeof (DictId);
    ASSERT (b64_len + 1 <= *snip_len, "Error in piz_get_other_ctx_from_snip: snip_len=%u but expecting it to be >= %u",
            *snip_len, b64_len + 1);

    DictId dict_id;
    base64_decode ((*snip)+1, &b64_len, dict_id.id);

    Context *other_ctx = mtf_get_ctx (vb, dict_id);

    *snip     += b64_len + 1;
    *snip_len -= b64_len + 1;
    
    return other_ctx;
}

void piz_reconstruct_one_snip (VBlock *vb, Context *snip_ctx, const char *snip, unsigned snip_len)
{
    if (!snip_len) return; // nothing to do
    
    int64_t new_value=0;
    bool have_new_value = false;
    Context *base_ctx = snip_ctx; // this will change if the snip refers us to another data source
    bool store = (snip_ctx->flags & CTX_FL_STORE_VALUE);

    switch (snip[0]) {

    // display the rest of the snip first, and then the lookup up text.
    case SNIP_LOOKUP:
    case SNIP_OTHER_LOOKUP: {

        if (snip[0] == SNIP_LOOKUP) 
            { snip++; snip_len--; }
        else 
            // we are request to reconstruct from another ctx
            base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len

        // case 1: LOCAL is not CTX_LT_SEQUENCE - we reconstruct this snip before adding the looked up data
        if (snip_len && !(base_ctx->ltype == CTX_LT_SEQUENCE)) RECONSTRUCT (snip, snip_len);
        
        if (base_ctx->ltype >= CTX_LT_INT8 && base_ctx->ltype <= CTX_LT_UINT64) {
            new_value = piz_reconstruct_from_local_int (vb, base_ctx, 0);
            have_new_value = true;
        }

        // case 2: CTX_LT_SEQUENCE  - the snip is taken to be the length of the sequence (or if missing, the length will be taken from vb->seq_len)
        else if (base_ctx->ltype == CTX_LT_SEQUENCE) 
            piz_reconstruct_from_local_sequence (vb, base_ctx, snip, snip_len);

/*        // case 2: CTX_LT_SEQUENCE_REF
        else if (base_ctx->ltype == CTX_LT_SEQUENCE_REF) 
            ref_reconstruct (vb, base_ctx, snip, snip_len);
*/
        else piz_reconstruct_from_local_text (vb, base_ctx); // this will call us back recursively with the snip retrieved
                
        break;
    }
    case SNIP_SELF_DELTA:
        new_value = piz_reconstruct_from_delta (vb, snip_ctx, base_ctx, snip+1, snip_len-1);
        have_new_value = true;
        break;

    case SNIP_OTHER_DELTA: 
        base_ctx = piz_get_other_ctx_from_snip (vb, &snip, &snip_len); // also updates snip and snip_len
        new_value = piz_reconstruct_from_delta (vb, snip_ctx, base_ctx, snip, snip_len); 
        have_new_value = true;
        break;

    case SNIP_STRUCTURED:
        piz_reconstruct_structured (vb, snip+1, snip_len-1);
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
        store = false; // override CTX_FL_STORE_VALUE and fall through
        snip++; snip_len--;
        
    default: {
        RECONSTRUCT (snip, snip_len); // simple reconstruction

        if (store) {
            char *after;
            new_value = (int64_t)strtoull (snip, &after, 10); // allows negative values

            // if the snip in its entirety is not a valid integer, don't store the value.
            // this can happen for example when seg_pos_field stores a "nonsense" snip.
            have_new_value = (after == snip + snip_len);
        }
        snip_ctx->last_delta = 0; // delta is 0 since we didn't calculate delta
    }
    }

    // update last_value if needed
    if (have_new_value && store) {
        ASSERT (base_ctx == snip_ctx, "Error in piz_reconstruct_one_snip: attempted to udpate value when base_ctx=%s differs from snip_ctx=%s. Are these supposed to be aliases?",
                base_ctx->name, snip_ctx->name);
        base_ctx->last_value = new_value;
    } 

    snip_ctx->last_line_i = vb->line_i;
}

// returns reconstructed length
uint32_t piz_reconstruct_from_ctx_do (VBlock *vb, uint8_t did_i, char sep)
{
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
        piz_reconstruct_one_snip (vb, ctx, snip, snip_len);        

        // handle chrom and pos to determine whether this line should be grepped-out in case of --regions
        if (did_i == DTF(chrom)) // test original did_i, not the alias target
            vb->chrom_node_index = word_index;

        if (flag_regions && did_i == DTF(pos) && !regions_is_site_included (vb->chrom_node_index, ctx->last_value)) 
            vb->dont_show_curr_line = true;
    }
    
    // case: all data is only in local
    else if (ctx->local.len) {
        if (ctx->ltype >= CTX_LT_INT8 && ctx->ltype <= CTX_LT_UINT64)
            piz_reconstruct_from_local_int(vb, ctx, 0);
        
        else if (ctx->ltype == CTX_LT_SEQUENCE) 
            piz_reconstruct_from_local_sequence (vb, ctx, NULL, 0);
        
        else if (ctx->ltype == CTX_LT_SEQUENCE_REF) 
            sam_piz_reconstruct_seq (vb, ctx);
        
        else if (ctx->ltype == CTX_LT_TEXT)
            piz_reconstruct_from_local_text (vb, ctx);

        else ABORT ("Invalid ltype=%u in ctx=%s of vb_i=%u line_i=%u", ctx->ltype, ctx->name, vb->vblock_i, vb->line_i);
    }

    // case: the entire VB was just \n - so seg dropped the ctx
    else if (ctx->did_i == DTF(eol))
        RECONSTRUCT1('\n');

    else ABORT("Error in piz_reconstruct_from_ctx_do: ctx %s has no data (b250 or local) in vb_i=%u line_i=%u", 
                ctx->name, vb->vblock_i, vb->line_i);

    if (sep) RECONSTRUCT1 (sep); 

    return (uint32_t)(vb->txt_data.len - start);
}

void piz_map_compound_field (VBlock *vb, bool (*predicate)(DictId), SubfieldMapper *mapper)
{
    mapper->num_subfields = 0;

    for (uint8_t did_i=0; did_i < vb->num_dict_ids; did_i++)
        if (predicate (vb->contexts[did_i].dict_id)) {
         
            char index_char = vb->contexts[did_i].dict_id.id[1];
            unsigned index = IS_DIGIT(index_char) ? index_char - '0' : 10 + index_char - 'a';
         
            mapper->did_i[index] = did_i;
            mapper->num_subfields++;
        }
}

uint32_t piz_uncompress_all_ctxs (VBlock *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);
    
    uint32_t section_i = 1;
    while (section_i < vb->z_section_headers.len) {

        if (section_i == vb->z_section_headers.len) break; // no more sections left

        SectionHeaderCtx *header = (SectionHeaderCtx *)ENT (char, vb->z_data, section_index[section_i]);

        bool is_local = header->h.section_type == SEC_LOCAL;
        if (header->h.section_type == SEC_B250 || is_local) {
            Context *ctx = mtf_get_ctx (vb, header->dict_id);
            ctx->flags = header->flags;
            ctx->ltype = header->ltype;

            zfile_uncompress_section (vb, header, is_local ? &ctx->local : &ctx->b250, 
                                      is_local ? "contexts.local" : "contexts.b250", header->h.section_type); 
            section_i++;
        }    

        else break;
    }

    return section_i;
}

static void piz_uncompress_one_vb (VBlock *vb)
{
    START_TIMER;

    // we read the header and ctxs for all data_types, except VCF that does its own special handling
    if (vb->data_type != DT_VCF) {
        ARRAY (const unsigned, section_index, vb->z_section_headers); 

        SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
        vb->first_line       = BGEN32 (header->first_line);      
        vb->lines.len        = BGEN32 (header->num_lines);       
        vb->vb_data_size     = BGEN32 (header->vb_data_size);    
        vb->longest_line_len = BGEN32 (header->longest_line_len);
        if (flag_split) vb->vblock_i = BGEN32 (header->h.vblock_i); /* in case of --split, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher because the dispatcher is re-initialized for every component */ 

        if (flag_show_vblocks) 
            fprintf (stderr, "vb_i=%u first_line=%u num_lines=%u txt_size=%u genozip_size=%u longest_line_len=%u\n",
                     vb->vblock_i, vb->first_line, (uint32_t)vb->lines.len, vb->vb_data_size, BGEN32 (header->z_data_bytes), vb->longest_line_len);

        buf_alloc (vb, &vb->txt_data, vb->vb_data_size + 10000, 1.1, "txt_data", vb->vblock_i); // +10000 as sometimes we pre-read control data (eg structured templates) and then roll back

        piz_uncompress_all_ctxs (vb);
    }

    DTP(uncompress)(vb);

    vb->is_processed = true; /* tell dispatcher this thread is done and can be joined. this operation needn't be atomic, but it likely is anyway */ 
    COPY_TIMER (vb->profile.compute);
}

static void piz_read_all_ctxs (VBlock *vb, SectionListEntry **next_sl)
{
    // ctxs that have dictionaries are already initialized, but others (eg local data only) are not
    mtf_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_dict_ids);

    while ((*next_sl)->section_type == SEC_B250 || (*next_sl)->section_type == SEC_LOCAL) {
        *ENT (unsigned, vb->z_section_headers, vb->z_section_headers.len) = vb->z_data.len; 

        int32_t ret = zfile_read_section (z_file, vb, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data", sizeof(SectionHeaderCtx), 
                                          (*next_sl)->section_type, *next_sl); // returns 0 if section is skipped

        if (ret) vb->z_section_headers.len++;

        (*next_sl)++;                             
    }
}

// Called by PIZ I/O thread: read all the sections at the end of the file, before starting to process VBs
static DataType piz_read_global_area (Md5Hash *original_file_digest) // out
{
    DataType data_type = zfile_read_genozip_header (original_file_digest);
    
    dict_id_initialize (data_type); // must run after zfile_read_genozip_header that sets z_file->data_type

    if (data_type == DT_NONE) return data_type;
    
    // for FASTA and FASTQ we convert a "header_only" flag to "header_one" as flag_header_only has some additional logic
    // that doesn't work for FAST
    if (flag_header_only && (data_type == DT_FASTA || data_type == DT_FASTQ)) {
        flag_header_only = false;
        flag_header_one  = true;
    }

    // check if the genozip file includes a reference
    if (sections_has_reference()) {
        ASSERTW (flag_reference == REF_NONE || flag_reference == REF_STORED, // it will be REF_STORED in the 2nd+ file
        "Note: --reference option is ignored - %s does not require a reference to be uncompressed", z_name);
        flag_reference = REF_STORED;
    }

    // if the user wants to see only the header, we can skip the dictionaries, regions and random access
    if (!flag_header_only) {
        
        // read random access, but only if we are going to need it
        bool need_random_access = flag_regions || flag_show_index || flag_fasta_sequential || flag_reference != REF_NONE;
        if (need_random_access) {
            zfile_read_all_dictionaries (0, DICTREAD_CHROM_ONLY); // read all CHROM/RNAME dictionaries - needed for regions_make_chregs()

            // update chrom node indeces using the CHROM dictionary, for the user-specified regions (in case -r/-R were specified)
            regions_make_chregs (dt_fields[data_type].chrom);

            // if the regions are negative, transform them to the positive complement instead
            regions_transform_negative_to_positive_complement();

            SectionListEntry *ra_sl = sections_get_offset_first_section_of_type (SEC_RANDOM_ACCESS);
            zfile_read_section (z_file, evb, 0, NO_SB_I, &evb->z_data, "z_data", sizeof (SectionHeader), SEC_RANDOM_ACCESS, ra_sl);

            zfile_uncompress_section (evb, evb->z_data.data, &z_file->ra_buf, "z_file->ra_buf", SEC_RANDOM_ACCESS);

            z_file->ra_buf.len /= random_access_sizeof_entry();
            BGEN_random_access();

            if (flag_show_index) {
                random_access_show_index(false);
                if (exe_type == EXE_GENOCAT) exit(0); // in genocat --show-index, we only show the index, not the data
            }

            buf_free (&evb->z_data);
        }

        // get the last vb_i that included in the regions - returns -1 if no vb has the requested regions
        int32_t last_vb_i = flag_regions ? random_access_get_last_included_vb_i() : 0;

        // read dictionaries
        if (last_vb_i >= 0)
            zfile_read_all_dictionaries (last_vb_i, need_random_access ? DICTREAD_EXCEPT_CHROM : DICTREAD_ALL);

        // give the reference FASTA CONTIG dictionary to the reference module
        if (ref_flag_reading_reference)
            ref_consume_ref_fasta_global_area();

        if (flag_reference == REF_STORED) { // note: in case of REF_EXTERNAL, reference is already pre-loaded
            ref_uncompress_all_ranges();
            if (!flag_quiet) fprintf (stderr, (flag_test && !ref_flag_reading_reference) ? "Success\n" : "Done\n");
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

    SectionListEntry *sl = sections_vb_first (vb->vblock_i); 

    int vb_header_offset = zfile_read_section (z_file, vb, vb->vblock_i, NO_SB_I, &vb->z_data, "z_data", 
                                               z_file->data_type == DT_VCF ? sizeof (SectionHeaderVbHeaderVCF) : sizeof (SectionHeaderVbHeader), 
                                               SEC_VB_HEADER, sl++); 

    ASSERT (vb_header_offset != EOF, "Error: unexpected end-of-file while reading vblock_i=%u", vb->vblock_i);
    mtf_overlay_dictionaries_to_vb ((VBlockP)vb); /* overlay all dictionaries (not just those that have fragments in this vblock) to the vb */ 

    buf_alloc (vb, &vb->z_section_headers, (MAX_DICTS * 2 + 50) * sizeof(char*), 0, "z_section_headers", 1); // room for section headers  
    NEXTENT (unsigned, vb->z_section_headers) = vb_header_offset; // vb_header_offset is always 0 for VB header

    // read all b250 and local of all fields and subfields
    piz_read_all_ctxs (vb, &sl);

    // read additional sections and other logic specific to this data type
    bool ok_to_compute = DTPZ(read_one_vb) ? DTPZ(read_one_vb)(vb, sl) : true; // true if we should go forward with computing this VB (otherwise skip it)

    COPY_TIMER (vb->profile.piz_read_one_vb); 

    return ok_to_compute;
}

// returns true is successfully outputted a txt file
bool piz_dispatcher (const char *z_basename, bool is_first_component, bool is_last_file)
{
    // static dispatcher - with flag_split, we use the same dispatcher when unzipping components
    static Dispatcher dispatcher = NULL;
    SectionListEntry *sl_ent = NULL;
    bool piz_successful = true;

    if (flag_split && !sections_has_more_components()) return false; // no more components
    
    // read genozip header
    Md5Hash original_file_digest;

    // read genozip header, dictionaries etc and set the data type when reading the first component of in case of --split, 
    static DataType data_type = DT_NONE; 
    if (is_first_component) {
        data_type = piz_read_global_area (&original_file_digest);

        ASSERT (sections_get_next_header_type(&sl_ent, NULL, NULL) == SEC_TXT_HEADER, "Error: unable to find TXT Header data in %s", z_name);

        ASSERT (!flag_test || !md5_is_zero (original_file_digest), 
                "Error testing %s: --test cannot be used with this file, as it was not compressed with --md5 or --test", z_name);
    }

    if (!dispatcher) 
        dispatcher = dispatcher_init (global_max_threads, 0, flag_test, is_last_file, z_basename);

    // case: we couldn't open the file because we didn't know what type it is - open it now
    if (!flag_split && !ref_flag_reading_reference && !txt_file->file) file_open_txt (txt_file);

    if (data_type == DT_NONE) goto finish;

    // read and write txt header. in split mode this also opens txt_file
    piz_successful = txtfile_genozip_to_txt_header (&original_file_digest);
    
    ASSERT (piz_successful || !is_first_component, "Error: failed to read %s header in %s", dt_name (z_file->data_type), z_name);

    if (!piz_successful || flag_header_only) goto finish;

    if (flag_split) dispatcher_resume (dispatcher); // accept more input 

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

                case SEC_TXT_HEADER: // 2nd+ txt header of a concatenated file
                    if (!flag_split) {
                        txtfile_genozip_to_txt_header (NULL); // skip 2nd+ txt header if concatenating
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
            if (!ref_flag_reading_reference) txtfile_write_one_vblock (processed_vb);

            z_file->num_vbs++;
            
            z_file->txt_data_so_far_single += processed_vb->vb_data_size; 

            dispatcher_finalize_one_vb (dispatcher);
        }

    } while (!dispatcher_is_done (dispatcher));

    // verify file integrity, if the genounzip compress was run with --md5 or --test
    if (flag_md5) {
        Md5Hash decompressed_file_digest;
        md5_finalize (&txt_file->md5_ctx_concat, &decompressed_file_digest); // z_file might be a concatenation - this is the MD5 of the entire concatenation

        if (md5_is_zero (original_file_digest) && !flag_quiet) 
            fprintf (stderr, "MD5 = %s Note: unable to compare this to the original file as file was not originally compressed with --md5\n", md5_display (&decompressed_file_digest, false));
        
        else if (md5_is_equal (decompressed_file_digest, original_file_digest)) {

            if (flag_test && !flag_quiet) fprintf (stderr, "Success          \b\b\b\b\b\b\b\b\b\b\n\n");

            if (flag_md5 && !flag_quiet) 
                fprintf (stderr, "MD5 = %s verified as identical to the original %s\n", 
                         md5_display (&decompressed_file_digest, false), dt_name (txt_file->data_type));
        }
        else if (flag_test) {
            fprintf (stderr, "FAILED!!!          \b\b\b\b\b\b\b\b\b\b\nError: MD5 of original file=%s is different than decompressed file=%s\nPlease contact bugs@genozip.com to help fix this bug in genozip",
                    md5_display (&original_file_digest, false), md5_display (&decompressed_file_digest, false));
            exit (1);
        }
        else ASSERT (md5_is_zero (original_file_digest), // its ok if we decompressed only a partial file
                     "File integrity error: MD5 of decompressed file %s is %s, but the original %s file's was %s", 
                     txt_file->name, md5_display (&decompressed_file_digest, false), dt_name (txt_file->data_type), 
                     md5_display (&original_file_digest, false));
    }

    if (flag_split) file_close (&txt_file, true); // close this component file

    if (!flag_test && !flag_quiet) 
        fprintf (stderr, "Done (%s)           \n", dispatcher_ellapsed_time (dispatcher, false));

finish:
    // in split mode - we continue with the same dispatcher in the next component. otherwise, we finish with it here
    if (!flag_split || !piz_successful) 
        dispatcher_finish (&dispatcher, NULL);
    else
        dispatcher_pause (dispatcher);

    return piz_successful;
}
