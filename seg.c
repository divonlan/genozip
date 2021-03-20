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
#include "data_types.h"
#include "container.h"
#include "codec.h"
#include "reference.h"

WordIndex seg_by_ctx_do (VBlock *vb, const char *snip, unsigned snip_len, Context *ctx, uint32_t add_bytes,
                         bool *is_new) // optional out
{
    ASSERTE0 (ctx, "ctx is NULL");

    buf_alloc_more (vb, &ctx->b250, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->b250");
    
    WordIndex node_index = ctx_evaluate_snip_seg ((VBlockP)vb, ctx, snip, snip_len, is_new);

    ASSERTE (node_index < ctx->nodes.len + ctx->ol_nodes.len || node_index == WORD_INDEX_EMPTY_SF || node_index == WORD_INDEX_MISSING_SF, 
             "out of range: dict=%s node_index=%d nodes.len=%u ol_nodes.len=%u",  
             ctx->name, node_index, (uint32_t)ctx->nodes.len, (uint32_t)ctx->ol_nodes.len);
    
    NEXTENT (uint32_t, ctx->b250) = node_index;
    ctx->txt_len += add_bytes;

    // a snip that is stored in its entirety in local, with just a LOOKUP in the dictionary, is counted as a singleton
    if (snip_len==1 && (snip[0] == SNIP_LOOKUP))
        ctx->num_singletons++;

    return node_index;
} 

void seg_simple_lookup (VBlockP vb, ContextP ctx, unsigned add_bytes)
{
    static const char lookup[1] = { SNIP_LOOKUP };
    seg_by_ctx (vb, lookup, 1, ctx, add_bytes);
}

const char *seg_get_next_item (void *vb_, const char *str, int *str_len, 
                               GetNextAllow newline, GetNextAllow tab, GetNextAllow space,
                               unsigned *len, char *separator, 
                               bool *has_13, // out - only modified if '\r' detected ; only needed if newline=GN_SEP
                               const char *item_name)
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
                    ASSERTE0 (has_13, "has_13=NULL but expecting it because newline=GN_SEP");
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
            MIN (i, 1000), str);

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
    
    ASSSEG (*str_len, str, "missing %s field", item_name);

    ABOSEG (str, "while segmenting %s: expecting a NEWLINE after (showing at most 1000 characters): \"%.*s\"", 
            item_name, MIN (i-1, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

void seg_integer_do (VBlockP vb, DidIType did_i, int64_t n, unsigned add_bytes)
{
    char snip[40];
    unsigned snip_len = str_int (n, snip);
    seg_by_did_i (vb, snip, snip_len, did_i, add_bytes);
}

// prepare snips that contain code + dict_id + optional parameter (SNIP_LOOKUP_OTHER, SNIP_OTHER_DELTA, SNIP_REDIRECTION...)
void seg_prepare_snip_other (uint8_t snip_code, DictId other_dict_id, bool has_parameter, int32_t parameter, 
                             char *snip, unsigned *snip_len) //  in / out
{
    // make sure we have enough memory
    unsigned required_len = 1 + base64_size (DICT_ID_LEN) + 11 /* max length of a int32 -1000000000 */ + 1 /* \0 */ ;
    ASSERTE (*snip_len >= required_len, "*snip_len=%u, but it needs to be at least %u", *snip_len, required_len);

    snip[0] = snip_code;
    *snip_len = 1 + base64_encode (other_dict_id.id, DICT_ID_LEN, &snip[1]);

    if (has_parameter)
        *snip_len += str_int (parameter, &snip[*snip_len]);
}

WordIndex seg_chrom_field (VBlock *vb, const char *chrom_str, unsigned chrom_str_len)
{
    ASSERTE0 (chrom_str_len, "chrom_str_len=0");

    bool is_new;
    WordIndex chrom_node_index = seg_by_did_i_ex (vb, chrom_str, chrom_str_len, CHROM, chrom_str_len+1, &is_new);

    // don't allow adding chroms if const_chroms (except "*")
    ASSSEG (!is_new || !flag.const_chroms || (chrom_str_len == 1 && chrom_str[0] == '*'), 
            chrom_str, "contig '%.*s' appears in file, but is missing in the %s header", chrom_str_len, chrom_str, dt_name (txt_file->data_type));
            
    random_access_update_chrom ((VBlockP)vb, chrom_node_index, chrom_str, chrom_str_len);

    return chrom_node_index;
}

// scans a pos field - in case of non-digit or not in the range [0,MAX_POS], either returns -1
// (if allow_nonsense) or errors
PosType seg_scan_pos_snip (VBlock *vb, const char *snip, unsigned snip_len, 
                           SegError *err) // out - if NULL, exception if error
{
    PosType value=0;
    bool is_int = str_get_int (snip, snip_len, &value); // value unchanged if not integer

    if (is_int && value >= 0 && value <= MAX_POS) {
        *err = ERR_SEG_NO_ERROR;
        return value; // all good
    }

    ASSSEG (err, snip, "position field must be an integer number between 0 and %"PRId64", seeing: %.*s", 
            MAX_POS, snip_len, snip);

    *err = is_int ? ERR_SEG_OUT_OF_RANGE : ERR_SEG_NOT_INTEGER;
    return value; // in- or out-of- range integer, or 0 if value is not integer
}

// returns POS value if a valid pos, or 0 if not
PosType seg_pos_field (VBlock *vb, 
                       DidIType snip_did_i,    // mandatory: the ctx the snip belongs to
                       DidIType base_did_i,    // mandatory: base for delta
                       bool seg_bad_snips_too, // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
                       const char *pos_str, unsigned pos_len, // option 1
                       PosType this_pos,       // option 2
                       unsigned add_bytes)
{
    Context *snip_ctx = &vb->contexts[snip_did_i];
    Context *base_ctx = &vb->contexts[base_did_i];

    snip_ctx->no_stons = true;
    base_ctx->flags.store = STORE_INT;

    base_ctx->no_stons = true;

    SegError err = ERR_SEG_NO_ERROR;
    if (pos_str) {
        this_pos = seg_scan_pos_snip (vb, pos_str, pos_len, &err);
        ASSERTE (seg_bad_snips_too || !err, "invalid value %.*s in %s vb=%u line_i=%u", 
                 pos_len, pos_str, vb->contexts[snip_did_i].name, vb->vblock_i, vb->line_i);

        // we accept out-of-range integer values for non-self-delta
        if (snip_did_i != base_did_i && err == ERR_SEG_OUT_OF_RANGE) err = ERR_SEG_NO_ERROR;

        if (err) {
            SAFE_ASSIGN (pos_str-1, SNIP_DONT_STORE);
            seg_by_ctx (vb, pos_str-1, pos_len+1, snip_ctx, add_bytes); 
            SAFE_RESTORE;
            snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
            return 0; // invalid pos
        }
    }

    else {
        // check out-of-range for self-delta
        if (snip_did_i == base_did_i && (this_pos < 0 || this_pos > MAX_POS)) {
            err = ERR_SEG_OUT_OF_RANGE;
            char snip[15] = { SNIP_DONT_STORE };
            unsigned snip_len = 1 + str_int (this_pos, &snip[1]);
            seg_by_ctx (vb, snip, snip_len, snip_ctx, add_bytes); 
            snip_ctx->last_delta = 0;  // on last_delta as we're PIZ won't have access to it - since we're not storing it in b250 
            return 0; // invalid pos
        }
    }

    PosType pos_delta = this_pos - base_ctx->last_value.i;
    
    snip_ctx->last_value.i = this_pos;

    // if the delta is too big, add this_pos (not delta) to local and put SNIP_LOOKUP in the b250
    // EXCEPT if it is the first vb (ie last_pos==0) because we want to avoid creating a whole RANDOM_POS
    // section in every VB just for a single entry in case of a nicely sorted file
    if ((pos_delta > MAX_POS_DELTA || pos_delta < -MAX_POS_DELTA) && base_ctx->last_value.i) {

        // store the value in store it in local - uint32
        buf_alloc (vb, &snip_ctx->local, MAX (snip_ctx->local.len + 1, vb->lines.len) * sizeof (uint32_t), CTX_GROWTH, "contexts->local");
        NEXTENT (uint32_t, snip_ctx->local) = BGEN32 (this_pos);
        snip_ctx->txt_len += add_bytes;

        snip_ctx->ltype  = LT_UINT32;

        // add a LOOKUP to b250
        static const char lookup[1] = { SNIP_LOOKUP };
        seg_by_ctx (vb, lookup, 1, snip_ctx, 0);

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
        unsigned total_len = sizeof (pos_delta_str);

        if (base_ctx == snip_ctx) {
            pos_delta_str[0] = SNIP_SELF_DELTA;
            total_len = 1;
        }
        else 
            seg_prepare_snip_other (SNIP_OTHER_DELTA, base_ctx->dict_id, false, 0, pos_delta_str, &total_len);

        unsigned delta_len = str_int (pos_delta, &pos_delta_str[total_len]);
        total_len += delta_len;

        seg_by_ctx (vb, pos_delta_str, total_len, snip_ctx, add_bytes);
    }
    // case: the delta is the negative of the previous delta - add a SNIP_SELF_DELTA with no payload - meaning negated delta
    else {
        char negated_delta = SNIP_SELF_DELTA; // no delta means negate the previous delta
        seg_by_did_i (vb, &negated_delta, 1, snip_did_i, add_bytes);
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
        seg_add_to_local_uint32 (vb, ctx_get_ctx (vb, dict_id), id_num, 0);
    }

    Context *ctx = ctx_get_ctx (vb, dict_id);
    ctx->no_stons = true;
    ctx->ltype  = LT_UINT32;

    // prefix the textual part with SNIP_LOOKUP_UINT32 if needed (we temporarily overwrite the previous separator or the buffer underflow area)
    unsigned new_len = id_snip_len - num_digits;
    SAFE_ASSIGN (&id_snip[-1], SNIP_LOOKUP); // we assign it anyway bc of the macro convenience, but we included it only if num_digits>0
    seg_by_ctx (vb, id_snip-(num_digits > 0), new_len + (num_digits > 0), ctx, id_snip_len + !!account_for_separator); // account for the entire length, and sometimes with \t
    SAFE_RESTORE;
}

// returns true if it was an integer
bool seg_integer_or_not (VBlockP vb, ContextP ctx, 
                         const char *this_value, unsigned this_value_len, unsigned add_bytes)
{
    // case: its an integer
    if (!ctx->no_stons && // we interpret no_stons as means also no moving ints to local (one of the reasons is that an int might actually be a float)
        str_get_int (this_value, this_value_len, &ctx->last_value.i) &&
        ctx->last_value.i >= 0 && ctx->last_value.i <= 0xffffffffULL) {

        // if this is the first snip in the ctx, 
        if (!ctx->local.len) { // first number
            ctx->dynamic_size_local = true;
            ctx->ltype = LT_UINT32; // set only upon storing the first number - if there are no numbers, leave it as LT_TEXT so it can be used for singletons
            
            if (!ctx->b250.len) // no non-numbers yet
                ctx->numeric_only = true; // until proven otherwise
        }

        // add to local
        buf_alloc_more (vb, &ctx->local, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->local");
        NEXTENT (uint32_t, ctx->local) = ctx->last_value.i;
        
        // add to b250
        seg_simple_lookup ((VBlockP)vb, ctx, add_bytes);

        return true;
    }

    // case: non-numeric snip
    else { 
        seg_by_ctx (vb, this_value, this_value_len, ctx, add_bytes);
        ctx->numeric_only = false;
        return false;
    }
}

typedef struct { const char *start; unsigned len; DictId dict_id; } InfoItem;

static int sort_by_subfield_name (const void *a, const void *b)  
{ 
    InfoItem *ina = (InfoItem *)a;
    InfoItem *inb = (InfoItem *)b;
    
    return strncmp (ina->start, inb->start, MIN (ina->len, inb->len));
}

// if ostatus is not known yet (seg of a Primary line), check for LIFTREJD subfield
static inline LiftOverStatus seg_info_field_get_ostatus_by_con (unsigned nitems, InfoItem *info_items)
{
    for (unsigned i=0; i < nitems; i++)
        if (info_items[i].dict_id.num == dict_id_INFO_LIFTREJD) 
            return LO_UNSUPPORTED;

    return LO_OK; // no LIFTREJD - Primary line is not rejected
}

static void seg_info_field_correct_for_dual_coordinates (VBlock *vb, Container *con, InfoItem *info_items, LiftOverStatus ostatus)
{
    // case: --chain and INFO is '.' - remove the '.' as we are adding INFO/LIFTOVER or LIFTREJD
    if (ostatus >= 0 && (con_nitems(*con) == 1 && info_items->len == 1 && *info_items->start == '.')) {
        con_dec_nitems (*con);
        vb->vb_data_size--;
    }

    // if ostatus is not known yet (seg of a Primary line), check for LIFTREJD subfield
    if (ostatus == LO_UNKNOWN) 
        ostatus = seg_info_field_get_ostatus_by_con (con_nitems(*con), info_items);

    // cases we add INFO_LIFTOVER: liftover is possible (LO_OK), but not if this is a Primary file that already contains LIFTOVER)
    if (ostatus == LO_OK && txt_file->dual_coords != DC_PRIMARY) {

        info_items[con_nitems(*con)] = (InfoItem){ .start   = INFO_LIFTOVER"=", 
                                                   .len     = INFO_LIFTOVER_LEN + 1, // +1 for the '='
                                                   .dict_id = (DictId)dict_id_INFO_LIFTOVER };  
        
        // case Laft: we're replacing LIFTBACK with LIFTOVER in the default reconstruction
        if (txt_file->dual_coords == DC_LAFT)
            vb->vb_data_size += (INFO_LIFTOVER_LEN - INFO_LIFTBACK_LEN);

        // case: --chain - we're adding this subfield to the default reconstruction
        else
            vb->vb_data_size += INFO_LIFTOVER_LEN + 1 + (con_nitems(*con) > 0); // +1 for '=', +1 for ';' if we already have item(s).

        con_inc_nitems (*con);
        con->filter_items = true; // filter will select which of LIFTOVER and LIFTBACK is reconstructed
    }

    // cases we add INFO_LIFTBACK: lift is possible (LO_OK), but not if this is a Laft file that already contains LIFTBACK
    if (ostatus == LO_OK && txt_file->dual_coords != DC_LAFT) {

        if (ostatus == LO_OK) {
            info_items[con_nitems(*con)] = (InfoItem){ .start   = INFO_LIFTBACK"=", 
                                                       .len     = INFO_LIFTBACK_LEN + 1, // +1 for the '='
                                                       .dict_id = (DictId)dict_id_INFO_LIFTBACK };  
            con_inc_nitems (*con);
            con->filter_items = true; // filter will select which of LIFTOVER and LIFTBACK is reconstructed
        }
        // note: we don't increase vb_data_size for LIFTBACK because its not displayed in the default reconstruction of the file
    }

    // case we add LIFTREJD: --chain with a liftover error (not Primary or Laft files - they already contain these)
    if (ostatus >= 1 && chain_is_loaded) {
        info_items[con_nitems(*con)] = (InfoItem){ .start   = INFO_LIFTREJD"=", 
                                                   .len     = INFO_LIFTREJD_LEN + 1, 
                                                   .dict_id = (DictId)dict_id_INFO_LIFTREJD };
        
        vb->vb_data_size += INFO_LIFTREJD_LEN + 1 + (con_nitems(*con) > 0); // +1 for '=', +1 for ';' if we already have item(s)

        con_inc_nitems (*con);
    }
}

// used to seg INFO fields in VCF and ATTRS in GFF3
void seg_info_field (VBlock *vb, SegSpecialInfoSubfields seg_special_subfields, const char *info_str, unsigned info_len,
                     LiftOverStatus ostatus) // >=0 --chain or dual coordinate file, or -1 if not
{
    const int info_field   = DTF(info);
    const char *field_name = DTF(names)[info_field];

    // max number of Container items consumed by Liftover stuff, reducing the max number of other subfields 
    unsigned max_liftover_items = ostatus == LO_NONE    ? 0  // not a dual-coordinates frile
                                : ostatus == LO_OK      ? 2  // INFO_LIFTOVER and INFO_LIFTBACK
                                : ostatus == LO_UNKNOWN ? 2  // Either INFO_LIFTOVER and INFO_LIFTBACK Or INFO_LIFTREJD
                                :                         1; // INFO_LIFTREJD

    Container con = { .repeats             = 1, 
                      .drop_final_item_sep = true };

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
                            "%s field contains a name %.*s starting with an illegal character '%c' (ASCII %u)", 
                            field_name, this_name_len, this_name, this_name[0], this_name[0]);

                    InfoItem *ii = &info_items[con_nitems(con)];
                    ii->start    = this_name; 
                    ii->len      = this_name_len + valueful; // include the '=' if there is one 
                    ii->dict_id  = valueful ? dict_id_make (this_name, this_name_len, DTYPE_1) 
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
                DictId dict_id = info_items[con_nitems(con)].dict_id;
                if (dict_id.num) { 
                    char optimized_snip[OPTIMIZE_MAX_SNIP_LEN];                
                    if (seg_special_subfields (vb, dict_id, &this_value, (unsigned *)&this_value_len, optimized_snip))
                        seg_integer_or_not (vb, ctx_get_ctx (vb, dict_id), this_value, this_value_len, this_value_len);

                    if (txt_file->dual_coords && (dict_id.num == dict_id_INFO_LIFTOVER || dict_id.num == dict_id_INFO_LIFTBACK))
                        con.filter_items = true; // filter will select which of LIFTOVER and LIFTBACK is reconstructed
                }

                reading_name = true;  // end of value - move to the next item
                this_name = &info_str[i+1]; // move to next field in info string
                this_name_len = 0;
                con_inc_nitems (con);

                ASSSEG (con_nitems(con) <= MAX_SUBFIELDS - max_liftover_items, info_str, 
                        "A line has too many subfields (tags) in the %s field - the maximum supported is %u",
                        field_name, MAX_SUBFIELDS - max_liftover_items);
            }
            else this_value_len++;
        }
    }

    // handle dual coordinates arising from --chain, as well as compressing DC_PRIMARY and DC_LAFT files
    seg_info_field_correct_for_dual_coordinates (vb, &con, info_items, ostatus);

    // finalize special handling
    seg_special_subfields (vb, DICT_ID_NONE, NULL, NULL, NULL);

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag.optimize_sort && con_nitems(con) > 1) 
        qsort (info_items, con_nitems(con), sizeof(InfoItem), sort_by_subfield_name);

    char prefixes[CONTAINER_MAX_PREFIXES_LEN]; // these are the Container prefixes
    prefixes[0] = prefixes[1] = CON_PREFIX_SEP; // initial CON_PREFIX_SEP follow by seperator of empty Container-wide prefix
    unsigned prefixes_len = 2;

    // Populate the Container 
    uint32_t total_names_len=0;
    for (unsigned i=0; i < con_nitems(con); i++) {
        // Set the Container item and find (or create) a context for this name
        InfoItem *ii = &info_items[i];
        con.items[i] = (ContainerItem){ .dict_id   = ii->dict_id,
                                        .seperator = { ';' } }; 
        // add to the prefixes
        ASSSEG (prefixes_len + ii->len + 1 <= CONTAINER_MAX_PREFIXES_LEN, info_str, 
                "%s contains tag names that, combined (including the '='), exceed the maximum of %u characters", field_name, CONTAINER_MAX_PREFIXES_LEN);

        memcpy (&prefixes[prefixes_len], ii->start, ii->len);
        prefixes_len += ii->len;
        prefixes[prefixes_len++] = CON_PREFIX_SEP;

        // don't include LIFTBACK in a dual-coordinates file in add_bytes because its not reconstructed by default
        if (!(ii->dict_id.num == dict_id_INFO_LIFTBACK && 
              (txt_file->dual_coords || ostatus == LO_OK /* --chain */)))  
            total_names_len += ii->len + 1; // +1 for ; \t or \n separator
    }

    container_seg_by_ctx (vb, &vb->contexts[info_field], &con, prefixes, prefixes_len, 
                          total_names_len /* names inc. = and separator */);
}

WordIndex seg_delta_vs_other (VBlock *vb, Context *ctx, Context *other_ctx, const char *value, unsigned value_len,
                              int64_t max_delta /* max abs value of delta - beyond that, seg as is, ignored if < 0 */)
{
    if (!other_ctx) goto fallback;

    if (!str_get_int (value, value_len, &ctx->last_value.i)) goto fallback;

    int64_t delta = ctx->last_value.i - other_ctx->last_value.i; 
    if (max_delta >= 0 && (delta > max_delta || delta < -max_delta)) goto fallback;

    char snip[100]; 
    unsigned snip_len = sizeof (snip);
    seg_prepare_snip_other (SNIP_OTHER_DELTA, other_ctx->dict_id, true, (int32_t)delta, snip, &snip_len);

    other_ctx->no_stons = true;
    other_ctx->flags.store = STORE_INT;

    ctx->numeric_only = false;
    return seg_by_ctx (vb, snip, snip_len, ctx, value_len);

fallback:
    return seg_by_ctx (vb, value, value_len, ctx, value_len);
}

#define MAX_COMPOUND_COMPONENTS 36
Container seg_initialize_container_array (VBlockP vb, DictId dict_id, bool type_1_items)
{
    Container con = (Container){ .repeats = 1 };

    for (unsigned i=0; i < MAX_COMPOUND_COMPONENTS; i++) {
        const uint8_t *id = dict_id.id;
        uint8_t char2 = (i < 10) ? ('0' + i) : ('a' + (i-10));
        
        char dict_id_str[8] = { id[0], char2, id[1], id[2], id[3], id[4], id[5], id[6] };
        
        con.items[i].dict_id = dict_id_make (dict_id_str, 8, type_1_items ? DTYPE_1 : DTYPE_2);
    }

    return con;
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
                         SegCompoundArg arg,   
                         unsigned nonoptimized_len, // if non-zero, we don't account for the string given, instead, only for this amount (+add_for_eol)
                         unsigned add_for_eol)      // account for characters beyond the component seperators
{
    // we use nodes.param in D?ESC contexts to track whether all snips in in this VB are the same
    // defaults to 0 (the same) and we set it to 1 if we encounter a different one
    #define not_all_the_same nodes.param 

    const char *snip = field;
    unsigned snip_len = 0;
    unsigned num_double_sep = 0;

    Container con = seg_initialize_container_array ((VBlockP)vb, field_ctx->dict_id, true); 

    // add each subfield to its dictionary - 2nd char is 0-9,a-z
    for (unsigned i=0; i <= field_len; i++) { // one more than field_len - to finalize the last subfield
    
        char sep = (i==field_len) ? 0 : field[i];

        if (!sep || 
            (con_nitems(con) < MAX_COMPOUND_COMPONENTS-1 && 
             ((sep==':' && arg.colon) || 
              (sep=='/' && arg.slash) ||
              (sep=='|' && arg.pipe)  || 
              (sep=='.' && arg.dot)   || 
              ((sep==' ' || sep=='\t' || sep==1) && arg.whitespace)))) {
        
            // process the subfield that just ended
            Context *sf_ctx = ctx_get_ctx (vb, con.items[con_nitems(con)].dict_id);
            ASSERTE (sf_ctx, "sf_ctx for %s is NULL", dis_dict_id (con.items[con_nitems(con)].dict_id).s);

            sf_ctx->st_did_i = field_ctx->did_i;

            // allocate memory if needed
            buf_alloc_more (vb, &sf_ctx->b250, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->b250");

            // if snip is an integer, we store a delta
            char delta_snip[30];
            unsigned original_snip_len = snip_len;

            PosType this_value;
            if (str_get_int (snip, snip_len, &this_value)) {
                delta_snip[0] = SNIP_SELF_DELTA;

                PosType delta = this_value - sf_ctx->last_value.i;

                // note: if all the snips so far in this VB are the same - store just the snip, so that if the 
                // entire b250 is the same, it can be removed
                if (delta || sf_ctx->not_all_the_same) {
                    snip_len = 1 + str_int (delta, &delta_snip[1]);
                    snip = delta_snip;

                    sf_ctx->flags.store = STORE_INT;
                    sf_ctx->not_all_the_same = true;
                }

                sf_ctx->last_value.i = this_value;
            }
            else if (flag.pair == PAIR_READ_1)
                sf_ctx->no_stons = true; // prevent singletons, so pair-2 can compare to us
            
            // we are evaluating but might throw away this snip and use SNIP_PAIR_LOOKUP instead - however, we throw away if its in the pair file,
            // i.e. its already in the dictionary and hash table - so no resources wasted
            uint32_t word_index = ctx_evaluate_snip_seg ((VBlockP)vb, sf_ctx, snip, snip_len, NULL);

            // case we are compressing fastq pairs - read 1 is the basis and thus must have a b250 node index,
            // and read 2 might have SNIP_PAIR_LOOKUP
            if (flag.pair == PAIR_READ_2) {
 
                // if the number of components in the compound is not exactly the same for every line of
                // pair 1 and pair 2 for this vb, the readings from the b250 will be incorrect, causing missed opportunities 
                // for SNIP_PAIR_LOOKUP and hence worse compression. This conditions makes sure this situation
                // doesn't result in an error (TO DO: overcome this, bug 159)
                // note: this can also happen if the there is no sf_ctx->pair do it being fully singletons and moved to local 
                if (sf_ctx->pair_b250_iter.next_b250 < AFTERENT (uint8_t, sf_ctx->pair) &&
                    // for pairing to word with SNIP_DELTA, if we have SNIP_PAIR_LOOKUP then all previous lines
                    // this VB must have been SNIP_PAIR_LOOKUP as well. Therefore, the first time we encounter an
                    // inequality - we stop the pairing going forward till the end of this VB
                    !sf_ctx->stop_pairing) {
                    
                    WordIndex pair_word_index = base250_decode (&sf_ctx->pair_b250_iter.next_b250, !sf_ctx->pair_flags.all_the_same, sf_ctx->name);  
                    
                    if (pair_word_index == WORD_INDEX_ONE_UP) 
                        pair_word_index = sf_ctx->pair_b250_iter.prev_word_index + 1;
                    
                    sf_ctx->pair_b250_iter.prev_word_index = pair_word_index;
                    
                    // note: if the pair word is a singleton in pair_1 file, then pair_word_index will be the index of {SNIP_LOOKUP}
                    // rather than the snip (as replaced in ctx_evaluate_snip_merge), therefore this condition will fail. This is quite
                    // rare, so not worth handling this case
                    if (word_index == pair_word_index) {
                        sf_ctx->pair_b250 = true;
                        static const char lookup_pair_snip[1] = { SNIP_PAIR_LOOKUP };
                        word_index = ctx_evaluate_snip_seg ((VBlockP)vb, sf_ctx, lookup_pair_snip, 1, NULL);
                    } else
                        // To improve: currently, pairing stops at the first non-match
                        sf_ctx->stop_pairing = true;
                }
                else
                    sf_ctx->stop_pairing = true;
            }

            NEXTENT (uint32_t, sf_ctx->b250) = word_index;

            sf_ctx->txt_len += nonoptimized_len ? 0 : original_snip_len;

            // case double space (common in fasta reference file description)
            bool double_sep = (sep==' ') && (i < field_len-1) && (field[i+1] == sep);
            if (double_sep) {
                i++;
                num_double_sep++;
            }

            // finalize this subfield and get ready for reading the next one
            if (i < field_len) {    
                con.items[con_nitems(con)].seperator[0] = sep;
                con.items[con_nitems(con)].seperator[1] = double_sep ? sep : 0;
                snip = &field[i+1];
                snip_len = 0;
            }
            con_inc_nitems (con);
        }
        else snip_len++;
    }

    container_seg_by_ctx (vb, field_ctx, &con, NULL, 0, (nonoptimized_len ? nonoptimized_len : con_nitems(con) + num_double_sep - 1) + add_for_eol);
}

// an array or array of arrays
// note: if container_ctx->flags.store=STORE_INT, container_ctx->last_value.i will be set to sum of integer
// elements including recursively from sub-arrays (non-integer elements will be ignored)
WordIndex seg_array (VBlock *vb, Context *container_ctx, DidIType stats_conslidation_did_i, 
                     const char *value, int32_t value_len, // must be signed
                     char sep, 
                     char subarray_sep,      // if non-zero, will attempt to find internal arrays
                     bool use_integer_delta) // first item stored as is, subsequent items stored as deltas
{
    MiniContainer *con;
    DictId arr_dict_id;
    Context *arr_ctx;

    // first use in this VB - prepare context where array elements will go in
    if (!container_ctx->con_cache.len) {
        const uint8_t *id = container_ctx->dict_id.id;
        arr_dict_id = (DictId){ .id = { id[0], 
                                        ((id[1]+1) % 256) | 128, // different IDs for top array, subarray and items
                                        id[2], id[3], id[4], id[5], id[6], id[7] } };
        
        buf_alloc (vb, &container_ctx->con_cache, sizeof (MiniContainer), 1, "contexts->con_cache");

        con = FIRSTENT (MiniContainer, container_ctx->con_cache);
        *con = (MiniContainer){ .nitems_lo = 1, 
                                .drop_final_item_sep = true,
                                .items     = { { .dict_id   = arr_dict_id, // only one item
                                                 .seperator = { sep } } } };

        arr_ctx = ctx_get_ctx (vb, arr_dict_id);
        arr_ctx->st_did_i = stats_conslidation_did_i;
    }
    else { 
        con         = FIRSTENT (MiniContainer, container_ctx->con_cache);
        arr_dict_id = con->items[0].dict_id;
        arr_ctx     = ctx_get_ctx (vb, arr_dict_id);
    }

    if (use_integer_delta || container_ctx->flags.store == STORE_INT) 
        arr_ctx->flags.store = STORE_INT;

    uint32_t add_bytes = (uint32_t)value_len; // we will deduct bytes added by sub-arrays

    // count repeats (1 + number of seperators)
    con->repeats=1;
    for (int32_t i=0; i < value_len; i++) 
        if (value[i] == sep) con->repeats++;

    if (container_ctx->flags.store == STORE_INT)
        container_ctx->last_value.i = 0; // initialize for summing elements

    for (uint32_t i=0; i < con->repeats; i++) { // value_len will be -1 after last number

        const char *this_item = value;
        int64_t this_item_value=0;

        bool is_subarray = false;
        for (; value_len && *value != sep; value++, value_len--) 
            if (*value == subarray_sep) 
                is_subarray = true;

        unsigned this_item_len  = (unsigned)(value - this_item);

        // case: its a sub-array
        if (is_subarray) {
            seg_array (vb, arr_ctx, stats_conslidation_did_i, this_item, this_item_len, subarray_sep, 0, use_integer_delta);
            add_bytes -= this_item_len; // sub-array will account for itself
            arr_ctx->numeric_only = false;
            this_item_value = arr_ctx->last_value.i;
        }

        // case: its an scalar (we don't delta arrays that have sub-arrays and we don't delta the first item)
        else if (!use_integer_delta || subarray_sep || i==0) {
            bool is_int = seg_integer_or_not (vb, arr_ctx, this_item, this_item_len, 0);
            if (is_int) this_item_value = arr_ctx->last_value.i;
        }

        // case: delta of 2nd+ item
        else if (str_get_int (this_item, this_item_len, &this_item_value)) {

            char delta_snip[30] = { SNIP_SELF_DELTA };
            unsigned delta_snip_len = 1 + str_int (this_item_value - arr_ctx->last_value.i, &delta_snip[1]);

            seg_by_ctx (vb, delta_snip, delta_snip_len, arr_ctx, 0);

            arr_ctx->last_value.i = this_item_value;
            arr_ctx->numeric_only = false;
        }

        // non-integer that cannot be delta'd - store as-is
        else {
            seg_by_ctx (vb, this_item, this_item_len, arr_ctx, 0);
            arr_ctx->numeric_only = false;
        }

        if (container_ctx->flags.store == STORE_INT)
            container_ctx->last_value.i += this_item_value;

        value_len--; // skip seperator
        value++;
    }

    return container_seg_by_ctx (vb, container_ctx, (ContainerP)con, 0, 0, add_bytes);
}

void seg_add_to_local_text (VBlock *vb, Context *ctx, 
                            const char *snip, unsigned snip_len, 
                            unsigned add_bytes)  // bytes in the original text file accounted for by this snip
{
    buf_alloc_more (vb, &ctx->local, snip_len + 1, vb->lines.len * (snip_len+1), char, CTX_GROWTH, "contexts->local");
    if (snip_len) buf_add (&ctx->local, snip, snip_len); 
    
    static const char sep[1] = { SNIP_SEP };
    buf_add (&ctx->local, sep, 1); 

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_fixed (VBlock *vb, Context *ctx, const void *data, unsigned data_len)  // bytes in the original text file accounted for by this snip
{
    if (data_len) {
        buf_alloc_more (vb, &ctx->local, data_len, vb->lines.len * data_len, char, CTX_GROWTH, "contexts->local");
        buf_add (&ctx->local, data, data_len); 
    }
}

void seg_add_to_local_uint8 (VBlockP vb, ContextP ctx, uint8_t value, unsigned add_bytes)
{
    buf_alloc_more (vb, &ctx->local, 1, vb->lines.len, uint8_t, CTX_GROWTH, "contexts->local");

    if (lt_desc[ctx->ltype].is_signed) value = INTERLACE (int8_t, value);
    NEXTENT (uint8_t, ctx->local) = value;

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_uint16 (VBlockP vb, ContextP ctx, uint16_t value, unsigned add_bytes)
{
    buf_alloc_more (vb, &ctx->local, 1, vb->lines.len, uint16_t, CTX_GROWTH, "contexts->local");

    if (lt_desc[ctx->ltype].is_signed) value = INTERLACE (int16_t, value);
    NEXTENT (uint16_t, ctx->local) = BGEN16 (value);

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_uint32 (VBlockP vb, ContextP ctx, uint32_t value, unsigned add_bytes)
{
    buf_alloc_more (vb, &ctx->local, 1, vb->lines.len, uint32_t, CTX_GROWTH, "contexts->local");

    if (lt_desc[ctx->ltype].is_signed) value = INTERLACE (int32_t, value);
    NEXTENT (uint32_t, ctx->local) = BGEN32 (value);

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_uint64 (VBlockP vb, ContextP ctx, uint64_t value, unsigned add_bytes)
{
    buf_alloc_more (vb, &ctx->local, 1, vb->lines.len, uint64_t, CTX_GROWTH, "contexts->local");

    if (lt_desc[ctx->ltype].is_signed) value = INTERLACE (int64_t, value);
    NEXTENT (uint64_t, ctx->local) = BGEN64 (value);

    if (add_bytes) ctx->txt_len += add_bytes;
}

static void seg_set_hash_hints (VBlock *vb, int third_num)
{
    if (third_num == 1) 
        vb->num_lines_at_1_3 = vb->line_i + 1;
    else 
        vb->num_lines_at_2_3 = vb->line_i + 1;

    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {

        Context *ctx = &vb->contexts[did_i];
        if (ctx->global_hash_prime) continue; // our service is not needed - global_cache for this dict already exists

        if (third_num == 1) 
            ctx->nodes_len_at_1_3 = ctx->nodes.len;

        else 
            ctx->nodes_len_at_2_3 = ctx->nodes.len;
    }
}

static uint32_t seg_estimate_num_lines (VBlock *vb)
{
#   define NUM_LINES_IN_TEST 10 // take an average of this number of lines
    uint32_t len=0; 

    if (!DTP (line_height)) return 0; // this data type doesn't use textual lines

    int newlines=0; for (; newlines < DTP (line_height) * NUM_LINES_IN_TEST; newlines++, len++)
        for (; len < vb->txt_data.len && vb->txt_data.data[len] != '\n'; len++) {};

    len /= NUM_LINES_IN_TEST; // average length of a line

    ASSERTE (vb->txt_data.len, "vb=%u: txt_data is empty", vb->vblock_i); 

    ASSSEG (newlines==DTP (line_height) || len < vb->txt_data.len, vb->txt_data.data, 
            "a line in the file is longer than %s characters (a maximum defined by vblock). If this is intentional, use --vblock to increase the vblock size", 
            str_uint_commas (flag.vblock_memory).s);

    return MAX (100, (uint32_t)(((double)vb->txt_data.len / (double)len) * 1.2));
}

// double the number of lines if we've run out of lines
static void seg_more_lines (VBlock *vb, unsigned sizeof_line)
{
    uint32_t num_old_lines = vb->lines.len;
    
    // note: sadly, we cannot use the normal Buffer macros here because each data_type has its own line type
    buf_alloc (vb, &vb->lines, (vb->lines.len + 1) * sizeof_line, 2, "lines");
    
    memset (&vb->lines.data[num_old_lines * sizeof_line], 0, vb->lines.size - num_old_lines * sizeof_line);
    
    vb->lines.len = vb->lines.size / sizeof_line;

    // allocate more to the b250 buffer of the fields, which each have num_lines entries
    for (int f=0; f < DTF(num_fields); f++) 
        if (buf_is_allocated (&vb->contexts[f].b250))
            buf_alloc_more_zero (vb, &vb->contexts[f].b250, vb->lines.len - num_old_lines, 0, uint32_t, 1, "contexts->b250");
}

static void seg_verify_file_size (VBlock *vb)
{
    uint32_t reconstructed_vb_size = 0;

    for (DidIType sf_i=0; sf_i < vb->num_contexts; sf_i++) 
        reconstructed_vb_size += vb->contexts[sf_i].txt_len;
        
    if (vb->vb_data_size != reconstructed_vb_size/* && !flag.data_modified*/) {

        fprintf (stderr, "Txt lengths:\n");
        for (DidIType sf_i=0; sf_i < vb->num_contexts; sf_i++) {
            Context *ctx = &vb->contexts[sf_i];
            fprintf (stderr, "%s: %u\n", ctx->name, (uint32_t)ctx->txt_len);
        }
        
        ABOSEG (vb->txt_data.data, "Error while verifying reconstructed vblock size: "
                "reconstructed_vb_size=%s (calculated by adding up ctx.txt_data after segging) but vb->vb_data_size=%s (calculated when reading data from the file) (diff=%d)", 
                str_uint_commas (reconstructed_vb_size).s, str_uint_commas (vb->vb_data_size).s, 
                (int32_t)reconstructed_vb_size - (int32_t)vb->vb_data_size);
    }
}

// split each lines in this variant block to its components
void seg_all_data_lines (VBlock *vb)
{
    START_TIMER;

    ctx_initialize_primary_field_ctxs (vb->contexts, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_contexts); // Create ctx for the fields in the correct order 

    ctx_verify_field_ctxs (vb);
 
    // allocate the b250 for the fields which each have num_lines entries
    for (int f=0; f < DTF(num_fields); f++) 
        buf_alloc (vb, &vb->contexts[f].b250, vb->lines.len * sizeof (uint32_t), 1, "contexts->b250");
    
    DT_FUNC (vb, seg_initialize)(vb);  // data-type specific initialization

    // get estimated number of lines, if we haven't already (eg in bam_seg_initialize)
    if (!vb->lines.len)
        vb->lines.len = seg_estimate_num_lines(vb);

    // allocate lines
    uint32_t sizeof_line = DT_FUNC_OPTIONAL (vb, sizeof_zip_dataline, 1)(); // 1 - we waste a little bit of memory to avoid making exceptions throughout the code logic
    buf_alloc (vb, &vb->lines, vb->lines.len * sizeof_line, 1, "lines");
    buf_zero (&vb->lines);

    const char *field_start = vb->txt_data.data;
    bool hash_hints_set_1_3 = false, hash_hints_set_2_3 = false;
    for (vb->line_i=0; vb->line_i < vb->lines.len; vb->line_i++) {

        uint32_t remaining_txt_len = AFTERENT (char, vb->txt_data) - field_start;
        
        if (!remaining_txt_len) { // we're done
            vb->lines.len = vb->line_i; // update to actual number of lines
            break;
        }

        //fprintf (stderr, "vb->line_i=%u\n", vb->line_i);
        bool has_13 = false;
        const char *next_field = DT_FUNC (vb, seg_txt_line) (vb, field_start, remaining_txt_len, &has_13);

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

    DT_FUNC (vb, seg_finalize)(vb); // data-type specific finalization

    if (!flag.make_reference) seg_verify_file_size (vb);

    COPY_TIMER (seg_all_data_lines);
}