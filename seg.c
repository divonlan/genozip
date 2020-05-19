// ------------------------------------------------------------------
//   seg.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "endianness.h"
#include "file.h"
#include "strings.h"
#include "optimize.h"
#include "random_access.h"
#include "dict_id.h"

void seg_init_mapper (VBlock *vb, int field_i, Buffer *mapper_buf, const char *name)
{
    if (!buf_is_allocated (&vb->mtf_ctx[field_i].ol_mtf)) return;
        
    mapper_buf->len = vb->mtf_ctx[field_i].ol_mtf.len;
    
    buf_alloc (vb, mapper_buf, mapper_buf->len * sizeof (SubfieldMapper), 2, name, 0);
    
    for (unsigned i=0; i < mapper_buf->len; i++) 
        ((SubfieldMapper *)mapper_buf->data)[i].num_subfields = (uint8_t)NIL;
}

// store src_buf in dst_buf, and frees src_buf. we attempt to "allocate" dst_buf using memory from txt_data,
// but in the part of txt_data has been already consumed and no longer needed.
// if there's not enough space in txt_data, we allocate on txt_data_spillover
void seg_store (VBlock *vb, 
                bool *dst_is_spillover, uint32_t *dst_start, uint32_t *dst_len, // out
                Buffer *src_buf, uint32_t size, // Either src_buf OR size must be given
                const char *limit_txt_data, // we cannot store in txt starting here. if NULL always allocates in txt_data_spillover
                bool align32) // does start address need to be 32bit aligned to prevent aliasing issues
{
    if (src_buf) size = src_buf->len;

    // align to 32bit (4 bytes) if needed
    if (align32 && (vb->txt_data_next_offset % 4))
        vb->txt_data_next_offset += 4 - (vb->txt_data_next_offset % 4);

    if (limit_txt_data && (vb->txt_data.data + vb->txt_data_next_offset + size < limit_txt_data)) { // we have space in txt_data - we can overlay
        *dst_is_spillover = false;
        *dst_start        = vb->txt_data_next_offset;
        *dst_len          = size;
        if (src_buf) memcpy (&vb->txt_data.data[*dst_start], src_buf->data, size);

        vb->txt_data_next_offset += size;
    }

    else {
        *dst_is_spillover = true;
        *dst_start = vb->txt_data_spillover.len;
        *dst_len = size;
        
        vb->txt_data_spillover.len += size;
        buf_alloc (vb, &vb->txt_data_spillover, MAX (1000, vb->txt_data_spillover.len), 1.5, 
                   "txt_data_spillover", vb->vblock_i);

        if (src_buf) memcpy (&vb->txt_data_spillover.data[*dst_start], src_buf->data, size);
    }

    if (src_buf) buf_free (src_buf);
}

static inline uint32_t seg_one_snip_do (VBlock *vb, const char *str, unsigned len, MtfContext *ctx, uint32_t add_bytes,
                                        bool *is_new) // optional out
{
    buf_alloc (vb, &ctx->mtf_i, MAX (vb->lines.len, ctx->mtf_i.len + 1) * sizeof (uint32_t),
               CTX_GROWTH, "mtf_ctx->mtf_i", ctx->did_i);
    
    uint32_t node_index = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, str, len, is_new);

    ASSERT (node_index < ctx->mtf.len + ctx->ol_mtf.len || node_index == WORD_INDEX_EMPTY_SF, 
            "Error in seg_one_field: out of range: dict=%s mtf_i=%d mtf.len=%u ol_mtf.len=%u",  
            ctx->name, node_index, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len);
    
    NEXTENT (uint32_t, ctx->mtf_i) = node_index;
    ctx->txt_len += add_bytes;

    return node_index;
} 

// returns the node index
uint32_t seg_one_subfield (VBlock *vb, const char *str, unsigned len, DictIdType dict_id, uint32_t add_bytes)
{
    MtfContext *ctx = mtf_get_ctx (vb, dict_id);
    return seg_one_snip_do (vb, str, len, ctx, add_bytes, NULL);
}

// returns the node index
uint32_t seg_one_snip (VBlock *vb, const char *str, unsigned len, int did_i, uint32_t add_bytes,
                       bool *is_new) // optional out
{
    MtfContext *ctx = &vb->mtf_ctx[did_i];
    return seg_one_snip_do (vb, str, len, ctx, add_bytes, is_new);
} 

const char *seg_get_next_item (void *vb_, const char *str, int *str_len, bool allow_newline, bool allow_tab, bool allow_colon, 
                               unsigned *len, char *separator, bool *has_13, // out - only needed if allow_newline=true
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

#define MAX_POS 0x7fffffff // maximum allowed value for POS

// reads a tab-terminated POS string
static int32_t seg_pos_snip_to_int (VBlock *vb, const char *pos_str, const char *field_name,
                                    bool *is_nonsense /* optional out */)
{
    // scan by ourselves - hugely faster the sscanf
    int64_t this_pos_64=0; // int64_t so we can test for overflow
    bool all_digits = true;
    const char *s; for (s=pos_str; *s != '\t' && *s != '\n' && *s != '\r' && *s != ';'; s++) {
        if (!IS_DIGIT (*s)) all_digits=false;

        if (all_digits)
            this_pos_64 = this_pos_64 * 10 + (*s - '0');
    }

    if (s == pos_str || !all_digits || this_pos_64 < 0 || this_pos_64 > 0x7fffffff) {
        ASSSEG (is_nonsense, pos_str, "Error: '%s' field must be an integer number between 0 and %u, seeing: %.*s", 
                field_name, MAX_POS, (int32_t)(s-pos_str), pos_str);

        *is_nonsense = true;
        return (int32_t)(s-pos_str);
    }
    else if (is_nonsense) *is_nonsense = false;

    return (int32_t)this_pos_64;
}

uint32_t seg_chrom_field (VBlock *vb, const char *chrom_str, unsigned chrom_str_len)
{
    ASSERT0 (chrom_str_len, "Error in seg_chrom_field: chrom_str_len=0");

    uint32_t chrom_node_index = seg_one_field (vb, chrom_str, chrom_str_len, DTF(chrom), chrom_str_len+1);

    random_access_update_chrom ((VBlockP)vb, chrom_node_index);

    return chrom_node_index;
}

// get number for storage in RANDOM_POS and check if it is a valid number
uint32_t seg_add_to_random_pos_data (VBlock *vb, const char *snip, unsigned snip_len, unsigned add_bytes, const char *field_name)
{
    bool standard_random_pos_encoding = true;

    // it is eligable for standard encoding only if the length is within the range
    if (!snip_len || snip_len > 10) standard_random_pos_encoding = false; // more than 10 digits is for sure bigger than 4GB=32bits

    // a multi-digit number cannot have a leading zero
    if (snip_len > 1 && snip[0] == '0') standard_random_pos_encoding = false;

    // it is eligable for standard encoding only if all the characters are digits
    uint64_t n64=0; // 64 bit just in case we go above 32 bit
    if (standard_random_pos_encoding)
        for (unsigned i=0; i < snip_len; i++) {
            if (!IS_DIGIT (snip[i])) {
                standard_random_pos_encoding = false;
                break;
            }
            n64 = n64 * 10 + (snip[i] - '0');
        }

    // it is eligable for standard encoding only if it fits in 32 bit (0xffffffff is reserved for a future escape, if needed)
    if (n64 > 0xfffffffe) standard_random_pos_encoding = false;

    ASSSEG (standard_random_pos_encoding, snip, "%s: Error: Bad position data in field %s - expecting an integer between 0 and %u without leading zeros, but instead seeing \"%.*s\"",
            global_cmd, field_name, 0xfffffffe, snip_len, snip);

    // note: random pos data always goes into the primary pos local - no point in anything several local sections with the same type of data
    MtfContext *ctx = &vb->mtf_ctx[DTF(pos)];
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + 1, vb->lines.len) * sizeof (uint32_t), CTX_GROWTH, ctx->name, sizeof (uint32_t));

    // 32 bit number - this compresses better than textual numbers with LZMA
    uint32_t n32 = (uint32_t)n64;
    NEXTENT (uint32_t, ctx->local) = BGEN32 (n32);

    ctx->txt_len += add_bytes;

    return n32;
}

int32_t seg_pos_field (VBlock *vb, int32_t last_pos, int32_t *last_pos_delta /*in /out */, 
                       bool allow_non_number, // should be FALSE if the file format spec expects this field to by a numeric POS, and true if we empirically see it is a POS, but we have no guarantee of it
                       int did_i, const char *pos_str, unsigned pos_len, 
                       const char *field_name, bool account_for_separator)
{
    bool is_nonsense = false;
    
    int32_t this_pos = seg_pos_snip_to_int (vb, pos_str, field_name, allow_non_number ? &is_nonsense : NULL);

    // if caller allows a non-valid-number and this is indeed a non-valid-number, just store the string, prefixed by SNIP_VERBTIM
    if (is_nonsense) { 
        int32_t nonsense_len = this_pos; // in case nonsense, seg_pos_snip_to_int returns the length
        char save = *(pos_str-1);
        *(char*)(pos_str-1) = SNIP_VERBTIM; // note: even if its the very first character in txt_data, we're fine - it will temporarily overwrite the buffer underflow

        seg_one_snip (vb, pos_str-1, nonsense_len+1, did_i, nonsense_len + account_for_separator, NULL); 

        *(char*)(pos_str-1) = save; // restore
        return last_pos; // unchanged last_pos
    }
    int32_t pos_delta = this_pos - last_pos;
    if (last_pos_delta) *last_pos_delta = pos_delta;
    
    // if the delta is too big, add the POS to RANDOM_POS and put \5 in the b250
    // EXCEPT if it is the first vb (ie last_pos==0) because we want to avoid creating a whole RANDOM_POS
    // section in every VB just for a single entry in case of a nicely sorted file
    if ((pos_delta > MAX_POS_DELTA || pos_delta < -MAX_POS_DELTA) && last_pos) {
        seg_add_to_random_pos_data (vb, pos_str, pos_len, pos_len+account_for_separator, field_name);
        static const char pos_lookup[1] = { SNIP_LOOKUP_UINT32 };
        seg_one_snip (vb, pos_lookup, 1, did_i, 0, NULL);

        return this_pos;
    }

    // print our string without expensive sprintf
    char pos_delta_str[50], reverse_pos_delta_str[50];
    
    // if the delta is the negative of the previous delta (as happens in unsorted BAM files with the second line in
    // each pair of lines) - we just store an empty snippet
    bool is_negated_last = (last_pos_delta && *last_pos_delta && (*last_pos_delta == -pos_delta));

    if (!is_negated_last) {

        bool negative = (pos_delta < 0);
        if (negative) pos_delta = -pos_delta;

        // create reverse string
        unsigned delta_len=0; 
        if (pos_delta) { 
            while (pos_delta) {
                reverse_pos_delta_str[delta_len++] = '0' + (pos_delta % 10);
                pos_delta /= 10;
            }
            if (negative) reverse_pos_delta_str[delta_len++] = '-';

            // reverse it
            for (unsigned i=0; i < delta_len; i++) pos_delta_str[i] = reverse_pos_delta_str[delta_len-i-1];
        }
        else { //pos_delta==0
            pos_delta_str[0] = '0';
            delta_len = 1;
        }
        
        seg_one_snip (vb, pos_delta_str, delta_len, did_i, pos_len + account_for_separator, NULL);
    }
    else {
        seg_one_snip (vb, "", 0, did_i, pos_len + account_for_separator, NULL);
        *last_pos_delta = 0; // no negated delta next time
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
void seg_id_field (VBlock *vb, DictIdType dict_id, const char *id_snip, unsigned id_snip_len, bool account_for_separator)
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
        MtfContext *ctx = mtf_get_ctx (vb, dict_id);
        uint32_t id_num = atoi (&id_snip[id_snip_len - num_digits]);
        
        buf_alloc (vb, &ctx->local, MAX (ctx->local.len + 1, vb->lines.len) * sizeof (uint32_t), CTX_GROWTH, ctx->name, sizeof(uint32_t));
        NEXTENT (uint32_t, ctx->local) = BGEN32 (id_num);
    }

    // append the textual part with SNIP_LOOKUP_UINT32 if needed - we have enough space - we have a tab following the field
    unsigned new_len = id_snip_len - num_digits;
    SAFE_ASSIGN (1, &id_snip[new_len], SNIP_LOOKUP_UINT32); // we assign it anyway bc of the macro convenience, but we included it only if num_digits>0
    seg_one_subfield (vb, id_snip, new_len + (num_digits > 0), dict_id, id_snip_len + !!account_for_separator); // account for the entire length, and sometimes with \t
    SAFE_RESTORE (1);
}

typedef struct { const char *start; unsigned len; } InfoNames;

static int sort_by_subfield_name (const void *a, const void *b)  
{ 
    InfoNames *ina = (InfoNames *)a;
    InfoNames *inb = (InfoNames *)b;
    
    return strncmp (ina->start, inb->start, MIN (ina->len, inb->len));
}

#define MAX_INFO_NAMES_LEN 1000 // max len of just the names string, without the data eg "INFO1=INFO2=INFO3="

static void seg_sort_iname (InfoNames *names, unsigned num_names, char *iname, unsigned *iname_len)
{
    if (! (*iname_len) || ((*iname_len) == 1 && iname[0]=='#')) return ;// nothing to sort

    if (iname[(*iname_len) - 1] == '#') num_names--; // we keep a final ";#" in its place

    qsort (names, num_names, sizeof(InfoNames), sort_by_subfield_name);

    char *next = iname;
    for (unsigned i=0; i < num_names; i++) {
        memcpy (next, // pointer into automatic variable iname of seg_info_field
                names[i].start, // pointer into db->txt_data 
                names[i].len);
        next += names[i].len;

        if (next[-1] != '=' && i < num_names-1) 
            *(next++) = ';'; // add ; after value-less name, except in the end
    }

    *iname_len = (unsigned)(next - iname);
}

// segments fields that look like INFO in VCF or ATTRIBUTES in GFF3
void seg_info_field (VBlock *vb, uint32_t *dl_info_mtf_i, Buffer *iname_mapper_buf, uint8_t *num_info_subfields,
                     SegSpecialInfoSubfields seg_special_subfields,
                     const char *info_str, unsigned info_len, 
                     bool this_field_has_13, // this is the last field in the line, and it ends with a Windows-style \r\n - we account for it in txt_len
                     bool this_line_has_13)  // this line ends with \r\n (this field may or may not be the last field) - we store this information as an info subfield for PIZ to recover
{
    // data type de-multiplexors
    #define info_field DTF(info)
    #define field_name DTF(names)[info_field]

    char iname[MAX_INFO_NAMES_LEN];
    unsigned iname_len = 0;
    const char *this_name = info_str;
    unsigned this_name_len = 0;
    const char *this_value = NULL;
    unsigned this_value_len=0;
    unsigned sf_i=0;
    char save_1, save_2=0 /* init to avoid compiler warning */;

    InfoNames names[MAX_SUBFIELDS];
    unsigned num_names=0;

    MtfContext *info_ctx = &vb->mtf_ctx[info_field];

    // if the txt file line ends with \r\n when we add an artificial additional info subfield "#"
    // we know we have space for adding ":#" because the line as at least a "\r\n" appearing somewhere after the INFO field
    if (this_line_has_13) {
        if (info_len) {
            save_2 = info_str[info_len];
            ((char*)info_str)[info_len++] = ';';
        }
        save_1 = info_str[info_len];
        ((char*)info_str)[info_len++] = '#';
    }
    
    // count infos
    SubfieldMapper iname_mapper;
    memset (&iname_mapper, 0, sizeof (iname_mapper));

    for (unsigned i=0; i < info_len; i++) 
        if (info_str[i] == '=') iname_mapper.num_subfields++;
    
    // get name / value pairs - and insert values to the "name" dictionary
    bool reading_name = true;
    for (unsigned i=0; i < info_len + 1; i++) {
        char c = (i==info_len) ? ';' : info_str[i]; // add an artificial ; at the end of the INFO data

        if (reading_name) {
            iname[iname_len++] = c; // info names inc. the =. the = terminats each name
            ASSSEG (iname_len <= MAX_INFO_NAMES_LEN, info_str, "Error: %s field too long, MAX_INFO_NAMES_LEN=%u", field_name, MAX_INFO_NAMES_LEN);

            if (c == '=') {  // end of name

                ASSSEG (this_name_len > 0, info_str, "Error: %s field contains a = without a preceding subfield name", field_name);

                ASSSEG (this_name[0] >= 64 && this_name[0] <= 127, info_str,
                        "Error: %s field contains a name %.*s starting with an illegal character", field_name, this_name_len, this_name);

                reading_name = false; 
                this_value = &info_str[i+1]; 
                this_value_len = 0;

                names[num_names].start = this_name; 
                names[num_names++].len = this_name_len + 1; // +1 for the '='
            }

            // name without value - valid in GFF3/VCF format
            else if (c == ';') {

                names[num_names].start = this_name; 
                names[num_names++].len = this_name_len;

                if (i==info_len) { // our artificial ; terminator
                    iname_len--; // remove ;
                    info_ctx->txt_len++;  // account for the separator (; or \t or \n)
                }
                else {
                    this_name = &info_str[i+1]; // skip the value-less name
                    this_name_len = 0;
                }
                continue;
            }
            
            else this_name_len++; // don't count the = or ; in the len
        }
        else {
            if (c == ';') { // end of value

                ASSSEG (this_value_len > 0, info_str,
                        "Error: %s field subfield %.*s, does not contain a value", field_name, this_name_len, this_name);

                // find (or create) an MTF context (= a dictionary) for this name
                DictIdType dict_id = dict_id_type_1 (dict_id_make (this_name, this_name_len));

                // find which DictId (did_i) this subfield belongs to (+ create a new ctx if this is the first occurance)
                MtfContext *ctx = mtf_get_ctx_sf (vb, num_info_subfields, dict_id);
                iname_mapper.did_i[sf_i] = ctx ? ctx->did_i : (uint8_t)NIL;

                // allocate memory if needed
                Buffer *mtf_i_buf = &ctx->mtf_i;
                buf_alloc (vb, mtf_i_buf, MIN (vb->lines.len, mtf_i_buf->len + 1) * sizeof (uint32_t),
                           CTX_GROWTH, "mtf_ctx->mtf_i", ctx->did_i);

                // Call back to handle special subfields
                char optimized_snip[OPTIMIZE_MAX_SNIP_LEN];                
                bool needs_evaluate = seg_special_subfields (vb, ctx, &this_value, &this_value_len, optimized_snip);
                    
                if (needs_evaluate) {
                    NEXTENT (uint32_t, *mtf_i_buf) = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, this_value, this_value_len, NULL);
                    ctx->txt_len += this_value_len;
                }

                info_ctx->txt_len++; // account for the separator (; or \t or \n) 

                reading_name = true;  // end of value - move to the next time
                this_name = &info_str[i+1]; // move to next field in info string
                this_name_len = 0;
                sf_i++;
            }
            else  
                this_value_len++;
        }
    }

    // if requested, we will re-sort the info fields in alphabetical order. This will result less words in the dictionary
    // thereby both improving compression and improving --regions speed. 
    if (flag_optimize_sort) seg_sort_iname (names, num_names, iname, &iname_len);

    // now insert the info names - a snip is a string that looks like: "INFO1=INFO2=INFO3="
    // 1. find it's mtf_i (and add to dictionary if a new name)
    // 2. place mtf_i in INFO section of this VB
    ARRAY (uint32_t, info_field_mtf_i, info_ctx->mtf_i);
    bool is_new;
    uint32_t node_index = mtf_evaluate_snip_seg ((VBlockP)vb, info_ctx, iname, iname_len, &is_new);
    info_field_mtf_i[vb->mtf_ctx[info_field].mtf_i.len++] = node_index;

    // if this is a totally new iname (first time in this file) - make a new SubfieldMapper for it.
    if (is_new) {   
        ASSSEG (node_index == iname_mapper_buf->len, info_str, "Error: node_index=%u different than iname_mapper_buf->len=%u", 
                node_index, (uint32_t)iname_mapper_buf->len);
    
        iname_mapper_buf->len++;
        buf_alloc (vb, iname_mapper_buf, MAX (100, iname_mapper_buf->len) * sizeof (SubfieldMapper), 1.5, "iname_mapper_buf", 0);
    }

    // it is possible that the iname_mapper is not set yet even though not new - if the node is from a previous VB and
    // we have not yet encountered in node in this VB
    *ENT (SubfieldMapper, *iname_mapper_buf, node_index) = iname_mapper;

    *dl_info_mtf_i = node_index;

    info_ctx->txt_len += iname_len + this_field_has_13; // this includes all the = and, in case INFO is the last field and terminated by \r\n, account for the \r
    
    // recover characters we temporarily changed
    if (this_line_has_13) {
        ((char*)info_str)[info_len-1] = save_1;
        info_ctx->txt_len--; // we accounted for this character, but it doesn't appear in the original txt
        if (info_len > 1) {
            ((char*)info_str)[info_len-2] = save_2;
            info_ctx->txt_len--; // we accounted for this character, but it doesn't appear in the original txt
        }
    }
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
                         MtfContext *field_ctx, const char *field, unsigned field_len, 
                         SubfieldMapper *mapper, DictIdType sf_dict_id,
                         bool ws_is_sep, // whitespace is separator - separate by ' ' at '\t'
                         bool account_for_13)
{
#define MAX_COMPOUND_COMPONENTS (10+26)

    const char *snip = field;
    unsigned snip_len = 0;
    unsigned sf_i = 0;
    char template[MAX_COMPOUND_COMPONENTS-1]; // separators is one less than the subfields

    // add each subfield to its dictionary - 2nd char is 0-9,a-z
    for (unsigned i=0; i <= field_len; i++) { // one more than field_len - to finalize the last subfield
    
        char sep = (i==field_len) ? 0 : field[i];

        if (!sep || 
            ((sf_i < MAX_COMPOUND_COMPONENTS-1) && (sep==':' || sep=='/' || sep=='|' || (ws_is_sep && (sep==' ' || sep==1))))) {
        
            // process the subfield that just ended
            MtfContext *sf_ctx;

            if (mapper->num_subfields == sf_i) { // new subfield in this VB (sf_ctx might exist from previous VBs)
                sf_dict_id.id[1] = (sf_i <= 9) ? (sf_i + '0') : (sf_i-10 + 'a');

                sf_ctx = mtf_get_ctx (vb, sf_dict_id);
                mapper->did_i[sf_i] = sf_ctx->did_i;
                mapper->num_subfields++;
            }
            else 
                sf_ctx = MAPPER_CTX (mapper, sf_i);

            ASSERT0 (sf_ctx, "Error in seg_compound_field: sf_ctx is NULL");

            // allocate memory if needed
            buf_alloc (vb, &sf_ctx->mtf_i, MAX (vb->lines.len, sf_ctx->mtf_i.len + 1) * sizeof (uint32_t),
                       CTX_GROWTH, "mtf_ctx->mtf_i", sf_ctx->did_i);

            NEXTENT (uint32_t, sf_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, sf_ctx, snip, snip_len, NULL);
            sf_ctx->txt_len += snip_len;

            // finalize this subfield and get ready for reading the next one
            if (i < field_len) {    
                template[sf_i] = field[i];
                snip = &field[i+1];
                snip_len = 0;
            }
            sf_i++;
        }
        else snip_len++;
    }

    // if template is empty, make it "*"
    if (sf_i==1) template[0] = '*';

    // add template to the field dictionary (note: template may be of length zero if field has no / or :)
    NEXTENT (uint32_t, field_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, field_ctx, template, MAX (1, sf_i-1), NULL);
    field_ctx->txt_len += sf_i + account_for_13; // sf_i has 1 for each separator including the terminating \t or \n / \r\n
}

void seg_add_to_local_text (VBlock *vb, MtfContext *ctx, 
                            const char *snip, unsigned snip_len, 
                            unsigned add_bytes)  // bytes in the original text file accounted for by this snip
{
    buf_alloc (vb, &ctx->local, MAX (ctx->local.len + snip_len + 1, vb->lines.len * (snip_len+1)), 2, ctx->name, sizeof(char)); // param=quantum of len
    if (snip_len) buf_add (&ctx->local, snip, snip_len); 
    
    static const char sep[1] = { SNIP_SEP };
    buf_add (&ctx->local, sep, 1); 

    if (add_bytes) ctx->txt_len += add_bytes;
}

void seg_add_to_local_fixed (VBlock *vb, MtfContext *ctx, const void *data, unsigned data_len)  // bytes in the original text file accounted for by this snip
{
    if (data_len) {
        buf_alloc (vb, &ctx->local, MAX (ctx->local.len + data_len, vb->lines.len * data_len), CTX_GROWTH, ctx->name, sizeof (char)); 
        buf_add (&ctx->local, data, data_len); 
    }
}

static void seg_set_hash_hints (VBlock *vb, int third_num)
{
    if (third_num == 1) 
        vb->num_lines_at_1_3 = vb->line_i + 1;
    else 
        vb->num_lines_at_2_3 = vb->line_i + 1;

    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {

        MtfContext *ctx = &vb->mtf_ctx[did_i];
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

    ASSSEG0 (len < vb->txt_data.len, vb->txt_data.data, "Error: cannot find a newline in the entire vb");

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
        buf_alloc_more_zero (vb, &vb->mtf_ctx[f].mtf_i, vb->lines.len - num_old_lines, 0, uint32_t, 1);
}

static void seg_verify_file_size (VBlock *vb)
{
    uint32_t reconstructed_vb_size = 0;

    for (unsigned sf_i=0; sf_i < vb->num_dict_ids; sf_i++) 
        reconstructed_vb_size += vb->mtf_ctx[sf_i].txt_len;
        
    if (vb->vb_data_size != reconstructed_vb_size && !flag_optimize) {

        fprintf (stderr, "Txt lengths:\n");
        for (unsigned sf_i=0; sf_i < vb->num_dict_ids; sf_i++) {
            MtfContext *ctx = &vb->mtf_ctx[sf_i];
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

    mtf_initialize_primary_field_ctxs (vb->mtf_ctx, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_dict_ids); // Create ctx for the fields in the correct order 

    mtf_verify_field_ctxs (vb);
    
    uint32_t sizeof_line = (uint32_t)DTP(sizeof_zip_dataline);

    if (!sizeof_line) sizeof_line=1; // we waste a little bit of memory to avoid making exceptions throughout the code logic
 
    // allocate lines
    buf_alloc (vb, &vb->lines, seg_estimate_num_lines(vb) * sizeof_line, 1.2, "lines", vb->vblock_i);
    buf_zero (&vb->lines);
    vb->lines.len = vb->lines.size / sizeof_line;

    // allocate the mtf_i for the fields which each have num_lines entries
    for (int f=0; f < DTF(num_fields); f++) 
        buf_alloc (vb, &vb->mtf_ctx[f].mtf_i, vb->lines.len * sizeof (uint32_t), 1, "mtf_ctx->mtf_i", f);
    
    if (DTP(seg_initialize)) DTP(seg_initialize) (vb); // data-type specific initialization

    const char *field_start = vb->txt_data.data;
    bool hash_hints_set_1_3 = false, hash_hints_set_2_3 = false;
    for (vb->line_i=0; vb->line_i < vb->lines.len; vb->line_i++) {

        if (field_start - vb->txt_data.data == vb->txt_data.len) { // we're done
            vb->lines.len = vb->line_i; // update to actual number of lines
            break;
        }

        //fprintf (stderr, "vb->line_i=%u\n", vb->line_i);

        const char *next_field = DTP(seg_data_line) (vb, field_start);
        
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

    seg_verify_file_size (vb);

    COPY_TIMER(vb->profile.seg_all_data_lines);
}