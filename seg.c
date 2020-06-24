// ------------------------------------------------------------------
//   seg.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

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

void seg_init_mapper (VBlock *vb, int field_i, Buffer *mapper_buf, const char *name)
{
    if (!buf_is_allocated (&vb->contexts[field_i].ol_mtf)) return;
        
    mapper_buf->len = vb->contexts[field_i].ol_mtf.len;
    
    buf_alloc (vb, mapper_buf, mapper_buf->len * sizeof (SubfieldMapper), 2, name, 0);
    
    for (unsigned i=0; i < mapper_buf->len; i++) 
        ((SubfieldMapper *)mapper_buf->data)[i].num_subfields = (uint8_t)NIL;
}

uint32_t seg_by_ctx (VBlock *vb, const char *snip, unsigned snip_len, Context *ctx, uint32_t add_bytes,
                     bool *is_new) // optional out
{
    buf_alloc (vb, &ctx->mtf_i, MAX (vb->lines.len, ctx->mtf_i.len + 1) * sizeof (uint32_t),
               CTX_GROWTH, "contexts->mtf_i", ctx->did_i);
    
    uint32_t node_index = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, snip, snip_len, is_new);

    ASSERT (node_index < ctx->mtf.len + ctx->ol_mtf.len || node_index == WORD_INDEX_EMPTY_SF, 
            "Error in seg_by_did_i: out of range: dict=%s mtf_i=%d mtf.len=%u ol_mtf.len=%u",  
            ctx->name, node_index, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len);
    
    NEXTENT (uint32_t, ctx->mtf_i) = node_index;
    ctx->txt_len += add_bytes;

    // a snip who is stored in its entirety in local, with just a LOOKUP in the dictionary, is counted as a singleton
    if (snip_len==1 && (snip[0] == SNIP_LOOKUP))
        ctx->num_singletons++;

    return node_index;
} 

const char *seg_get_next_item (void *vb_, const char *str, int *str_len, bool allow_newline, bool allow_tab, bool allow_colon, 
                               unsigned *len, char *separator, 
                               bool *has_13, // out - only needed if allow_newline=true
                               const char *item_name)
{
    VBlockP vb = (VBlockP)vb_;

    unsigned i=0; for (; i < *str_len; i++)
        if ((allow_tab     && str[i] == '\t') ||
            (allow_colon   && str[i] == ':')  ||
            (allow_newline && str[i] == '\n')) {
                *len = i;
                *separator = str[i];
                *str_len -= i+1;

                // check for Windows-style '\r\n' end of line 
                if (i && str[i] == '\n' && str[i-1] == '\r') {
                    (*len)--;
                    ASSERT0 (has_13, "Error in seg_get_next_item: has_13==NULL but expecting it because allow_newline=true");
                    *has_13 = true;
                }

                return str + i+1; // beyond the separator
        }
        else if ((!allow_tab     && str[i] == '\t') ||  // note: a colon with allow_colon=false is not an error, its just part of the string rather than being a separator
                 (!allow_newline && str[i] == '\n')) break;
            
    ASSSEG (*str_len, str, "Error: missing %s field", item_name);

    ASSSEG (str[i] != '\t' || vb->data_type != DT_VCF || (item_name == DTF(names)[VCF_INFO] /* pointer, not string, comparison */), 
            str, "Error: while segmenting %s: expecting a NEWLINE after the INFO field, because this VCF file has no samples (individuals) declared in the header line",
            item_name);

    ABOSEG (str, "Error: while segmenting %s: expecting a %s %s %s after \"%.*s\"", 
            item_name,
            allow_newline ? "NEWLINE" : "", allow_tab ? "TAB" : "", allow_colon ? "\":\"" : "", 
            MIN (i-1, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

const char *seg_get_next_line (void *vb_, const char *str, int *str_len, unsigned *len, bool *has_13 /* out */, const char *item_name)
{
    VBlockP vb = (VBlockP)vb_;

    unsigned i=0; for (; i < *str_len; i++)
        if (str[i] == '\n') {
                *len = i;
                *str_len -= i+1;

                // check for Windows-style '\r\n' end of line 
                if (i && str[i-1] == '\r') {
                    (*len)--;
                    *has_13 = true;
                }

                return str + i+1; // beyond the separator
        }
    
    ASSSEG (*str_len, str, "Error: missing %s field", item_name);

    ABOSEG (str, "Error: while segmenting %s: expecting a NEWLINE after (showing at most 1000 characters): \"%.*s\"", 
            item_name, MIN (i-1, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

void seg_prepare_snip_other (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int32_t parameter, 
                             char *snip, unsigned *snip_len) // out
{
    snip[0] = snip_code;
    *snip_len = 1 + base64_encode (other_dict_id.id, DICT_ID_LEN, &snip[1]);

    if (has_parameter)
        *snip_len += str_int (parameter, &snip[*snip_len]);
}

uint32_t seg_chrom_field (VBlock *vb, const char *chrom_str, unsigned chrom_str_len)
{
    ASSERT0 (chrom_str_len, "Error in seg_chrom_field: chrom_str_len=0");

    uint8_t chrom_did_i = DTF(chrom);
    uint32_t chrom_node_index = seg_by_did_i (vb, chrom_str, chrom_str_len, chrom_did_i, chrom_str_len+1);

    random_access_update_chrom ((VBlockP)vb, chrom_node_index, chrom_str, chrom_str_len);

    return chrom_node_index;
}

// scans a pos field - in case of non-digit or not in the range [0,MAX_POS], either returns -1
// (if allow_nonsense) or errors
int64_t seg_scan_pos_snip (VBlock *vb, const char *snip, unsigned snip_len, bool allow_nonsense)
{
    char *after;
    int64_t value = (int64_t)strtoull (snip, &after, 10);

    if (value >= 0 && value <= MAX_POS && ((unsigned)(after - snip) == snip_len))
        return value; // all good

    ASSSEG (allow_nonsense, snip, "Error: position field must be an integer number between 0 and %"PRId64", seeing: %.*s", 
            MAX_POS, snip_len, snip);

    return -1; // bad number
}

// returns POS value if a valid pos, or 0 if not
int64_t seg_pos_field (VBlock *vb, 
                       uint8_t snip_did_i,    // mandatory: the ctx the snip belongs to
                       uint8_t base_did_i,    // mandatory: base for delta
                       bool allow_non_number, // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
                       const char *pos_str, unsigned pos_len, 
                       bool account_for_separator)
{
    Context *snip_ctx = &vb->contexts[snip_did_i];
    Context *base_ctx = &vb->contexts[base_did_i];

    snip_ctx->flags |= CTX_FL_LOCAL_LZMA  | CTX_FL_NO_STONS; 
    base_ctx->flags |= CTX_FL_STORE_VALUE | CTX_FL_NO_STONS;

    int64_t this_pos = seg_scan_pos_snip (vb, pos_str, pos_len, allow_non_number);

    // < 0  -  caller allows a non-valid-number and this is indeed a non-valid-number, just store the string
    // == 0 - special case where pos=0, e.g. "not available" in SAM_PNEXT. we just store it verbatim
    // In both cases, we store as SNIP_DONT_STORE so that piz doesn't update last_value after reading this value
    if (this_pos <= 0) { 
        SAFE_ASSIGN (1, pos_str-1, SNIP_DONT_STORE);
        seg_by_ctx (vb, pos_str-1, pos_len+1, snip_ctx, pos_len + account_for_separator, NULL); 
        SAFE_RESTORE (1);

        snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
        return 0; // invalid pos
    }

    int64_t pos_delta = this_pos - base_ctx->last_value;
    
    // if we're self-delta'ing - we store the value
    if (snip_ctx == base_ctx) base_ctx->last_value = this_pos; 

    // if the delta is too big, add this_pos (not delta) to local and put SNIP_LOOKUP in the b250
    // EXCEPT if it is the first vb (ie last_pos==0) because we want to avoid creating a whole RANDOM_POS
    // section in every VB just for a single entry in case of a nicely sorted file
    if ((pos_delta > MAX_POS_DELTA || pos_delta < -MAX_POS_DELTA) && base_ctx->last_value) {
        
        // store the value in store it in local - uint32
        buf_alloc (vb, &snip_ctx->local, MAX (snip_ctx->local.len + 1, vb->lines.len) * sizeof (uint32_t), CTX_GROWTH, snip_ctx->name, snip_ctx->did_i);
        NEXTENT (uint32_t, snip_ctx->local) = BGEN32 (this_pos);
        snip_ctx->txt_len += pos_len + account_for_separator;

        snip_ctx->ltype  = CTX_LT_UINT32;

        // add a LOOKUP to b250
        static const char lookup[1] = { SNIP_LOOKUP };
        seg_by_ctx (vb, lookup, 1, snip_ctx, 0, NULL);

        snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
        return this_pos;
    }

    // store the delta in last_delta only if we're also putting in the b250
    snip_ctx->last_delta = pos_delta;
    
    // if the delta is the negative of the previous delta (as happens in unsorted BAM files with the second line in
    // each pair of lines) - we just store an empty snippet
    bool is_negated_last = (snip_ctx->last_delta && snip_ctx->last_delta == -pos_delta);

    // case: add a delta
    if (!is_negated_last) {
        
        char pos_delta_str[100]; // more than enough for the base64
        unsigned total_len;

        if (base_ctx == snip_ctx) {
            pos_delta_str[0] = SNIP_SELF_DELTA;
            total_len = 1;
        }
        else 
            seg_prepare_snip_other (SNIP_OTHER_DELTA, base_ctx->dict_id, false, 0, pos_delta_str, &total_len);

        unsigned delta_len = str_int (pos_delta, &pos_delta_str[total_len]);
        total_len += delta_len;

        seg_by_ctx (vb, pos_delta_str, total_len, snip_ctx, pos_len + account_for_separator, NULL);
    }
    // case: the delta is the negative of the previous delta - add a SNIP_SELF_DELTA with no payload - meaning negated delta
    else {
        char negated_delta = SNIP_SELF_DELTA; // no delta means negate the previous delta
        seg_by_did_i (vb, &negated_delta, 1, snip_did_i, pos_len + account_for_separator);
        snip_ctx->last_delta = 0; // no negated delta next time
    }

    return this_pos;
}

// Commonly (but not always), IDs are SNPid identifiers like "rs17030902". We store the ID divided to 2:
// - We store the final digits, if any exist, and up to 9 digits, as an integer in SEC_NUMERIC_ID_DATA, which is
//   compressed with LZMA
// - In the dictionary we store the prefix up to this number, and \1 if there is a number and a \2 
//   if the caller (as in seg_me23) wants us to store an extra bit.
// example: rs17030902 : in the dictionary we store "rs\1" or "rs\1\2" and in SEC_NUMERIC_ID_DATA we store 17030902.
//          1423       : in the dictionary we store "\1" and 1423 SEC_NUMERIC_ID_DATA
//          abcd       : in the dictionary we store "abcd" and nothing is stored SEC_NUMERIC_ID_DATA
void seg_id_field (VBlock *vb, DictId dict_id, const char *id_snip, unsigned id_snip_len, bool account_for_separator)
{
    int i=id_snip_len-1; for (; i >= 0; i--) 
        if (!IS_DIGIT (id_snip[i])) break;
    
    unsigned num_digits = MIN (id_snip_len - (i+1), 9);

    // leading zeros will be part of the dictionary data, not the number
    for (unsigned i = id_snip_len - num_digits; i < id_snip_len; i++) 
        if (id_snip[i] == '0') 
            num_digits--;
        else 
            break;

    // added to local if we have a trailing number
    if (num_digits) {
        uint32_t id_num = atoi (&id_snip[id_snip_len - num_digits]);
        seg_add_to_local_uint32 (vb, mtf_get_ctx (vb, dict_id), id_num, 0);
    }

    Context *ctx = mtf_get_ctx (vb, dict_id);
    ctx->flags |= CTX_FL_ID;
    ctx->ltype  = CTX_LT_UINT32;

    // prefix the textual part with SNIP_LOOKUP_UINT32 if needed (we temporarily overwrite the previous separator or the buffer underflow area)
    unsigned new_len = id_snip_len - num_digits;
    SAFE_ASSIGN (1, &id_snip[-1], SNIP_LOOKUP); // we assign it anyway bc of the macro convenience, but we included it only if num_digits>0
    seg_by_ctx (vb, id_snip-(num_digits > 0), new_len + (num_digits > 0), ctx, id_snip_len + !!account_for_separator, NULL); // account for the entire length, and sometimes with \t
    SAFE_RESTORE (1);
}

typedef struct { const char *start; unsigned len; DictId dict_id; } InfoItem;

static int sort_by_subfield_name (const void *a, const void *b)  
{ 
    InfoItem *ina = (InfoItem *)a;
    InfoItem *inb = (InfoItem *)b;
    
    return strncmp (ina->start, inb->start, MIN (ina->len, inb->len));
}

// used to seg INFO fields in VCF and ATTRS in GFF3
void seg_info_field (VBlock *vb, SegSpecialInfoSubfields seg_special_subfields, const char *info_str, unsigned info_len)
{
    const int info_field   = DTF(info);
    const char *field_name = DTF(names)[info_field];

    Structured st = { .repeats=1, .num_items=0, .repsep={0,0}, .flags=STRUCTURED_DROP_LAST_SEP_OF_LAST_ELEMENT };

    const char *this_name = info_str, *this_value = NULL;
    int this_name_len = 0, this_value_len=0; // int and not unsigned as it can go negative

    InfoItem info_items[MAX_SUBFIELDS];

    // get name / value pairs - and insert values to the "name" dictionary
    bool reading_name = true;
    for (unsigned i=0; i < info_len + 1; i++) {
        char c = (i==info_len) ? ';' : info_str[i]; // add an artificial ; at the end of the INFO data

        if (reading_name) {

            if (c == '=' || c == ';') {  // end of valueful or valueless name

                bool valueful = (c == '=');

                ASSSEG (this_name_len > 0, info_str, "Error: %s field contains a = or ; without a preceding subfield name", field_name);

                if (this_name_len > 0) { 
                    ASSSEG ((this_name[0] >= 64 && this_name[0] <= 127) || this_name[0] == '.', info_str,
                            "Error: %s field contains a name %.*s starting with an illegal character '%c' (ASCII %u)", 
                            field_name, this_name_len, this_name, this_name[0], this_name[0]);

                    InfoItem *ii = &info_items[st.num_items];
                    ii->start    = this_name; 
                    ii->len      = this_name_len + valueful; // include the '=' if there is one 
                    ii->dict_id  = valueful ? dict_id_type_1 (dict_id_make (this_name, this_name_len)) 
                                            : DICT_ID_NONE;

                    this_value = &info_str[i+1]; 
                    this_value_len = -valueful; // if there is a '=' to be skipped, start from -1
                    reading_name = false; 
                }
            }
            else this_name_len++; // don't count the = or ; in the len
        }
        
        if (!reading_name) {

            if (c == ';') { // end of value
                // If its a valueful item, seg it (either special or regular)
                DictId dict_id = info_items[st.num_items].dict_id;
                if (dict_id.num) { 
                    char optimized_snip[OPTIMIZE_MAX_SNIP_LEN];                
                    bool not_yet_segged = seg_special_subfields (vb, dict_id, &this_value, (unsigned *)&this_value_len, optimized_snip);
                        
                    if (not_yet_segged) seg_by_dict_id (vb, this_value, this_value_len, dict_id, this_value_len);
                }

                reading_name = true;  // end of value - move to the next item
                this_name = &info_str[i+1]; // move to next field in info string
                this_name_len = 0;
                st.num_items++;

                ASSSEG (st.num_items <= MAX_SUBFIELDS, info_str, "A line has too many subfields (tags) in the %s field - the maximum supported is %u",
                        field_name, MAX_SUBFIELDS);
            }
            else this_value_len++;
        }
    }

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag_optimize_sort && st.num_items > 1) 
        qsort (info_items, st.num_items, sizeof(InfoItem), sort_by_subfield_name);

    char prefixes[STRUCTURED_MAX_PREFIXES_LEN]; // these are the Structured prefixes
    prefixes[0] = prefixes[1] = SNIP_STRUCTURED; // initial SNIP_STRUCTURED follow by seperator of empty Structured-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Structured 
    uint32_t total_names_len=0;
    for (unsigned i=0; i < st.num_items; i++) {
        // Set the Structured item and find (or create) a context for this name
        StructuredItem *si = &st.items[i];
        InfoItem *ii  = &info_items[i];
        si->dict_id   = ii->dict_id;
        si->seperator[0] = ';'; 
        si->seperator[1] = 0; 
        si->did_i     = DID_I_NONE; // this must be NONE, it is used only by PIZ

        // add to the prefixes
        ASSSEG (prefixes_len + ii->len + 1 <= STRUCTURED_MAX_PREFIXES_LEN, info_str, 
                "%s contains tag names that, combined (including the '='), exceed the maximum of %u characters", field_name, STRUCTURED_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii->start, ii->len);
        prefixes_len += ii->len;
        prefixes[prefixes_len++] = SNIP_STRUCTURED;

        total_names_len += ii->len;
    }

    seg_structured_by_ctx (vb, &vb->contexts[info_field], &st, prefixes, prefixes_len, 
                            total_names_len /* names inc. = */ + (st.num_items-1) /* the ;s */ + 1 /* \t or \n */);
}

void seg_structured_by_ctx (VBlock *vb, Context *ctx, Structured *st, 
                            // prefixes can be one of 3 options:
                            // 1. NULL
                            // 2. a "structured-wide prefix" that will be reconstructed once, at the beginning of the Structured
                            // 3. a "structured-wide prefix" followed by exactly one prefix per item. the per-item prefixes will be
                            //    displayed once per repeat, before their respective items. in this case, the structured-wide prefix
                            //    may be empty. 
                            // Each prefix is terminated by a SNIP_STRUCTURED character
                            const char *prefixes, unsigned prefixes_len, 
                            unsigned add_bytes)
{
    st->repeats = BGEN32 (st->repeats);
    char snip[1 + base64_sizeof(Structured) + STRUCTURED_MAX_PREFIXES_LEN]; // maximal size
    snip[0] = SNIP_STRUCTURED;
    unsigned b64_len = base64_encode ((uint8_t*)st, sizeof_structured (*st), &snip[1]);
    st->repeats = BGEN32 (st->repeats); // restore

    if (prefixes_len) memcpy (&snip[1+b64_len], prefixes, prefixes_len);

    ctx->flags |= CTX_FL_STRUCTURED;

    seg_by_ctx (vb, snip, 1 + b64_len + prefixes_len, ctx, add_bytes, NULL); 
}

#define MAX_COMPOUND_COMPONENTS 36

void seg_initialize_compound_structured (VBlockP vb, char *name_template, Structured *st)
{
    memset (st, 0, sizeof (Structured));
    unsigned name_len = MIN (strlen (name_template), DICT_ID_LEN);

    char name[DICT_ID_LEN] = {}; // modifiable string
    memcpy (name, name_template, name_len);

    for (unsigned i=0; i < MAX_COMPOUND_COMPONENTS; i++) {
        name[1] = (i < 10) ? ('0' + i) : ('a' + (i-10));
        st->items[i].dict_id = dict_id_type_1 (dict_id_make (name, name_len)); // both FASTQ and SAM use type1 for their compound field (QNAME and DESC respectively)
        st->items[i].did_i   = DID_I_NONE;
        mtf_get_ctx (vb, st->items[i].dict_id); // create ctx
    }

    st->repeats = 1;
}

// We break down the field (eg QNAME in SAM or Description in FASTA/FASTQ) into subfields separated by / and/or : - and/or whitespace 
// these are vendor-defined strings.
// Up to MAX_COMPOUND_COMPONENTS subfields are permitted - if there are more, then all the trailing part is just
// consider part of the last component.
// each subfield is stored in its own dictionary- the second character of the dict_id  the subfield number starting
// from 0 (0->9,a->z)
// The separators are made into a string we call "template" that is stored in the main field dictionary - we
// anticipate that usually all lines have the same format, but we allow lines to have different formats.
void seg_compound_field (VBlock *vb, 
                         Context *field_ctx, const char *field, unsigned field_len, 
                         SubfieldMapper *mapper, Structured st,
                         bool ws_is_sep, // whitespace is separator - separate by ' ' at '\t'
                         unsigned add_for_eol) // account for characters beyond the component seperators
{
    
    const char *snip = field;
    unsigned snip_len = 0;
    unsigned sf_i = 0;
        
    // add each subfield to its dictionary - 2nd char is 0-9,a-z
    for (unsigned i=0; i <= field_len; i++) { // one more than field_len - to finalize the last subfield
    
        char sep = (i==field_len) ? 0 : field[i];

        if (!sep || 
            ((sf_i < MAX_COMPOUND_COMPONENTS-1) && (sep==':' || sep=='/' || sep=='|' || (ws_is_sep && (sep==' ' || sep==1))))) {
        
            // process the subfield that just ended
            Context *sf_ctx;

            if (mapper->num_subfields == sf_i) { // new subfield in this VB (sf_ctx might exist from previous VBs)
                sf_ctx = mtf_get_ctx (vb, st.items[sf_i].dict_id);
                mapper->did_i[sf_i] = sf_ctx->did_i;
                mapper->num_subfields++;
            }
            else 
                sf_ctx = MAPPER_CTX (mapper, sf_i);

            ASSERT0 (sf_ctx, "Error in seg_compound_field: sf_ctx is NULL");

            // allocate memory if needed
            buf_alloc (vb, &sf_ctx->mtf_i, MAX (vb->lines.len, sf_ctx->mtf_i.len + 1) * sizeof (uint32_t),
                       CTX_GROWTH, "contexts->mtf_i", sf_ctx->did_i);

            NEXTENT (uint32_t, sf_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, sf_ctx, snip, snip_len, NULL);
            sf_ctx->txt_len += snip_len;

            // finalize this subfield and get ready for reading the next one
            if (i < field_len) {    
                st.items[sf_i].seperator[0] = field[i];
                st.items[sf_i].seperator[1] = 0;
                snip = &field[i+1];
                snip_len = 0;
            }
            sf_i++;
        }
        else snip_len++;
    }

    st.num_items = sf_i;

    seg_structured_by_ctx (vb, field_ctx, &st, NULL, 0, sf_i-1 + add_for_eol);
}

void seg_array_field (VBlock *vb, DictId dict_id, const char *value, unsigned value_len, 
                      SegOptimize optimize) // optional optimization function
{   
    const char *str = value; 
    int str_len = (int)value_len; // must be int, not unsigned, for the for loop
    
    Structured st = { .num_items = 1, .flags = STRUCTURED_DROP_LAST_SEP_OF_LAST_ELEMENT, 
                      .repsep = {0,0}, .items = { { .seperator = {','}, .did_i = DID_I_NONE } } };
    DictId arr_dict_id = dict_id_make ("XX_ARRAY", 8);
    arr_dict_id.id[0]      = FLIP_CASE (dict_id.id[0]);
    arr_dict_id.id[1]      = FLIP_CASE (dict_id.id[1]);
    st.items[0].dict_id    = sam_dict_id_optnl_sf (arr_dict_id);
    
    Context *arr_ctx = mtf_get_ctx (vb, st.items[0].dict_id);

    for (st.repeats=0; st.repeats < STRUCTURED_MAX_REPEATS && str_len > 0; st.repeats++) { // str_len will be -1 after last number

        const char *snip = str;
        for (; str_len && *str != ','; str++, str_len--) {};

        unsigned number_len = (unsigned)(str - snip);

        if (st.repeats == STRUCTURED_MAX_REPEATS-1) // final permitted repeat - take entire remaining string
            number_len += str_len;

        unsigned snip_len = number_len; 
             
        char new_number_str[30];
        if (optimize && st.repeats < STRUCTURED_MAX_REPEATS-1 && snip_len < 25)
            optimize (&snip, &snip_len, new_number_str);

        seg_by_ctx (vb, snip, snip_len, arr_ctx, number_len+1, NULL);
        
        str_len--; // skip comma
        str++;
    }

    seg_structured_by_dict_id (vb, dict_id, &st, 0);
}

void seg_add_to_local_text (VBlock *vb, Context *ctx, 
                            const char *snip, unsigned snip_len, 
                            unsigned add_bytes)  // bytes in the original text file accounted for by this snip
{
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + snip_len + 1, vb->lines.len * (snip_len+1)), 2, ctx->name, ctx->did_i); 
    if (snip_len) buf_add (&ctx->local, snip, snip_len); 
    
    static const char sep[1] = { SNIP_SEP };
    buf_add (&ctx->local, sep, 1); 

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_fixed (VBlock *vb, Context *ctx, const void *data, unsigned data_len)  // bytes in the original text file accounted for by this snip
{
    if (data_len) {
        buf_alloc (vb, &ctx->local, MAX (ctx->local.len + data_len, vb->lines.len * data_len), CTX_GROWTH, ctx->name, ctx->did_i); 
        buf_add (&ctx->local, data, data_len); 
    }
}

// SIGNED NUMBERS ARE NOT UNTEST YET! NOT USE YET BY ANY SEG
// for signed numbers, we store them in our "interlaced" format rather than standard ISO format 
// example signed: 2, -5 <--> interlaced: 4, 9. Why? for example, a int32 -1 will be 0x00000001 rather than 0xfffffffe - 
// compressing better in an array that contains both positive and negative
#define SAFE_NEGATE(type,n) ((u##type)(-((int64_t)n))) // careful negation to avoid overflow eg -(-128)==0 in int8_t
#define INTERLACE(type,n) ((((type)n) < 0) ? ((SAFE_NEGATE(type,n) << 1) - 1) : (((u##type)n) << 1))

void seg_add_to_local_uint8 (VBlockP vb, ContextP ctx, uint8_t value, unsigned add_bytes)
{
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + 1, vb->lines.len) * sizeof (uint8_t), CTX_GROWTH, ctx->name, ctx->did_i);

    if (ctx_lt_is_signed[ctx->ltype]) value = INTERLACE (int8_t, value);
    NEXTENT (uint8_t, ctx->local) = value;

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_uint16 (VBlockP vb, ContextP ctx, uint16_t value, unsigned add_bytes)
{
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + 1, vb->lines.len) * sizeof (uint16_t), CTX_GROWTH, ctx->name, ctx->did_i);

    if (ctx_lt_is_signed[ctx->ltype]) value = INTERLACE (int16_t, value);
    NEXTENT (uint16_t, ctx->local) = BGEN16 (value);

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_uint32 (VBlockP vb, ContextP ctx, uint32_t value, unsigned add_bytes)
{
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + 1, vb->lines.len) * sizeof (uint32_t), CTX_GROWTH, ctx->name, ctx->did_i);

    if (ctx_lt_is_signed[ctx->ltype]) value = INTERLACE (int32_t, value);
    NEXTENT (uint32_t, ctx->local) = BGEN32 (value);

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_uint64 (VBlockP vb, ContextP ctx, uint64_t value, unsigned add_bytes)
{
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + 1, vb->lines.len) * sizeof (uint64_t), CTX_GROWTH, ctx->name, ctx->did_i);

    if (ctx_lt_is_signed[ctx->ltype]) value = INTERLACE (int64_t, value);
    NEXTENT (uint64_t, ctx->local) = BGEN64 (value);

    if (add_bytes) ctx->txt_len += add_bytes;
}

static void seg_set_hash_hints (VBlock *vb, int third_num)
{
    if (third_num == 1) 
        vb->num_lines_at_1_3 = vb->line_i + 1;
    else 
        vb->num_lines_at_2_3 = vb->line_i + 1;

    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {

        Context *ctx = &vb->contexts[did_i];
        if (ctx->global_hash_prime) continue; // our service is not needed - global_cache for this dict already exists

        if (third_num == 1) 
            ctx->mtf_len_at_1_3 = ctx->mtf.len;

        else 
            ctx->mtf_len_at_2_3 = ctx->mtf.len;
    }
}

static uint32_t seg_estimate_num_lines (VBlock *vb)
{
    // get first line length
    uint32_t len=0; for (; len < vb->txt_data.len && vb->txt_data.data[len] != '\n'; len++) {};

    char s[30];
    ASSSEG (len < vb->txt_data.len, vb->txt_data.data, "Error: a line in the file is longer than %s characters (a maximum defined by vblock). If this is intentional, use --vblock to increase the vblock size", 
            str_uint_commas (global_max_memory_per_vb, s));

    return MAX (100, (uint32_t)(((double)vb->txt_data.len / (double)len) * 1.2));
}

static void seg_more_lines (VBlock *vb, unsigned sizeof_line)
{
    uint32_t num_old_lines = vb->lines.len;
    
    // note: sadly, we cannot use the normal Buffer macros here because each data_type has its own line type
    buf_alloc (vb, &vb->lines, vb->lines.size + sizeof_line, 2, "lines", vb->vblock_i);
    memset (&vb->lines.data[num_old_lines * sizeof_line], 0, vb->lines.size - num_old_lines * sizeof_line);
    
    vb->lines.len = vb->lines.size / sizeof_line;

    // allocate more to the mtf_i buffer of the fields, which each have num_lines entries
    for (int f=0; f < DTF(num_fields); f++) 
        buf_alloc_more_zero (vb, &vb->contexts[f].mtf_i, vb->lines.len - num_old_lines, 0, uint32_t, 1);
}

static void seg_verify_file_size (VBlock *vb)
{
    uint32_t reconstructed_vb_size = 0;

    for (unsigned sf_i=0; sf_i < vb->num_dict_ids; sf_i++) 
        reconstructed_vb_size += vb->contexts[sf_i].txt_len;
        
    if (vb->vb_data_size != reconstructed_vb_size && !flag_optimize) {

        fprintf (stderr, "Txt lengths:\n");
        for (unsigned sf_i=0; sf_i < vb->num_dict_ids; sf_i++) {
            Context *ctx = &vb->contexts[sf_i];
            fprintf (stderr, "%s: %u\n", ctx->name, (uint32_t)ctx->txt_len);
        }
        
        char s1[30], s2[30];
        ABOSEG (vb->txt_data.data, "Error while verifying reconstructed vblock size: "
                "reconstructed_vb_size=%s (calculated bottoms-up) but vb->vb_data_size=%s (calculated tops-down) (diff=%d)", 
                str_uint_commas (reconstructed_vb_size, s1), str_uint_commas (vb->vb_data_size, s2), 
                (int32_t)reconstructed_vb_size - (int32_t)vb->vb_data_size);
    }
}

// split each lines in this variant block to its components
void seg_all_data_lines (VBlock *vb)
{
    START_TIMER;

    mtf_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_dict_ids); // Create ctx for the fields in the correct order 

    mtf_verify_field_ctxs (vb);
    
    uint32_t sizeof_line = DTP(sizeof_zip_dataline) ? DTP(sizeof_zip_dataline)() : 0;

    if (!sizeof_line) sizeof_line=1; // we waste a little bit of memory to avoid making exceptions throughout the code logic
 
    // allocate lines
    buf_alloc (vb, &vb->lines, seg_estimate_num_lines(vb) * sizeof_line, 1.2, "lines", vb->vblock_i);
    buf_zero (&vb->lines);
    vb->lines.len = vb->lines.size / sizeof_line;

    // allocate the mtf_i for the fields which each have num_lines entries
    for (int f=0; f < DTF(num_fields); f++) 
        buf_alloc (vb, &vb->contexts[f].mtf_i, vb->lines.len * sizeof (uint32_t), 1, "contexts->mtf_i", f);
    
    if (DTP(seg_initialize)) DTP(seg_initialize) (vb); // data-type specific initialization

    const char *field_start = vb->txt_data.data;
    bool hash_hints_set_1_3 = false, hash_hints_set_2_3 = false;
    bool does_any_line_have_13 = false;
    for (vb->line_i=0; vb->line_i < vb->lines.len; vb->line_i++) {

        if (field_start - vb->txt_data.data == vb->txt_data.len) { // we're done
            vb->lines.len = vb->line_i; // update to actual number of lines
            break;
        }

        //fprintf (stderr, "vb->line_i=%u\n", vb->line_i);
        bool has_13 = false;
        const char *next_field = DTP(seg_txt_line) (vb, field_start, &has_13);
        if (has_13) does_any_line_have_13 = true;

        vb->longest_line_len = MAX (vb->longest_line_len, (next_field - field_start));
        field_start = next_field;

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
    }

    // if no line has special EOL, we can get rid of the EOL ctx
    if (!does_any_line_have_13) {
        Context *eol_ctx = &vb->contexts[DTF(eol)];
        buf_free (&eol_ctx->dict);
        buf_free (&eol_ctx->mtf);
        buf_free (&eol_ctx->mtf_i);
        buf_free (&eol_ctx->local);
    }

    seg_verify_file_size (vb);

    COPY_TIMER(vb->profile.seg_all_data_lines);
}