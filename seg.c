// ------------------------------------------------------------------
//   seg.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "profiler.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "header.h"
#include "endianness.h"
#include "file.h"
#include "strings.h"

// store src_bug in dst_buf, and frees src_buf. we attempt to "allocate" dst_buf using memory from txt_data,
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

// returns the node index
uint32_t seg_one_subfield (VBlock *vb, const char *str, unsigned len,
                           DictIdType dict_id, SectionType sec_b250, int accounts_for_chars)
{
    MtfNode *node;
    MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, dict_id, sec_b250 - 1);

    // allocate memory if needed
    buf_alloc (vb, &ctx->mtf_i, MAX (vb->lines.len, ctx->mtf_i.len + 1) * sizeof (uint32_t),
               CTX_GROWTH, "mtf_ctx->mtf_i", ctx->did_i);

    uint32_t node_index = mtf_evaluate_snip_seg (vb, ctx, str, len, &node, NULL);

    ASSERT (node_index < ctx->mtf.len + ctx->ol_mtf.len || node_index == WORD_INDEX_EMPTY_SF, 
            "Error in seg_one_subfield: out of range: dict=%s %s mtf_i=%d mtf.len=%u ol_mtf.len=%u",  
            err_dict_id (ctx->dict_id), st_name (ctx->dict_section_type),
            node_index, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len);

    NEXTENT (uint32_t, ctx->mtf_i) = node_index;

    vb->txt_section_bytes[sec_b250] += accounts_for_chars; 

    return node_index;
}

// returns the node index
uint32_t seg_one_snip (VBlock *vb, const char *str, unsigned len, int did_i, SectionType sec_b250,
                       bool *is_new) // optional out
{
    MtfContext *ctx = &vb->mtf_ctx[did_i];

    buf_alloc (vb, &ctx->mtf_i, (ctx->mtf_i.len + 1) * sizeof (uint32_t), 2, "mtf_ctx->mtf_i", did_i);
    
    MtfNode *node;
    uint32_t node_index = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, str, len, &node, is_new);

    ASSERT (node_index < ctx->mtf.len + ctx->ol_mtf.len || node_index == WORD_INDEX_EMPTY_SF, 
            "Error in seg_one_field: out of range: dict=%s %s mtf_i=%d mtf.len=%u ol_mtf.len=%u",  
            err_dict_id (ctx->dict_id), st_name (ctx->dict_section_type),
            node_index, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len);
    
    NEXTENT (uint32_t, vb->mtf_ctx[did_i].mtf_i) = node_index;

    vb->txt_section_bytes[sec_b250] += len + 1;
    return node_index;
} 

const char *seg_get_next_item (void *vb_, const char *str, int *str_len, bool allow_newline, bool allow_tab, bool allow_colon, 
                               unsigned *len, char *separator, bool *has_13, // out
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
                    *has_13 = true;
                }

                return str + i+1; // beyond the separator
        }
        else if ((!allow_tab     && str[i] == '\t') ||  // note: a colon with allow_colon=false is not an error, its just part of the string rather than being a separator
                 (!allow_newline && str[i] == '\n')) break;
            
    ASSSEG (*str_len, str, "Error: missing %s field", item_name);

    ASSSEG (str[i] != '\t' || vb->data_type != DT_VCF || strcmp (item_name, field_names[DT_VCF][VCF_INFO]), 
            str, "Error: while segmenting %s: expecting a NEWLINE after the INFO field, because this VCF file has no samples (individuals) declared in the header line",
            item_name);

    ABOSEG (str, "Error: while segmenting %s: expecting a %s %s %s after \"%.*s\"", 
            item_name,
            allow_newline ? "NEWLINE" : "", allow_tab ? "TAB" : "", allow_colon ? "\":\"" : "", 
            MIN (i-1, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

#define MAX_POS 0x7fffffff // maximum allowed value for POS

// reads a tab-terminated POS string
int32_t seg_pos_snip_to_int (VBlock *vb, const char *pos_str, const char *field_name)
{
    // scan by ourselves - hugely faster the sscanf
    int64_t this_pos_64=0; // int64_t so we can test for overflow
    const char *s; for (s=pos_str; *s != '\t' && *s != '\n' && *s != '\r'; s++) {
        ASSSEG (IS_DIGIT (*s), pos_str, "Error: '%s' field must be an integer number between 0 and %u, seeing: %.*s", 
                field_name, MAX_POS, (int)(s-pos_str+1), pos_str);

        this_pos_64 = this_pos_64 * 10 + (*s - '0');
    }
    ASSSEG (this_pos_64 >= 0 && this_pos_64 <= 0x7fffffff, pos_str,
            "Error: Invalid '%s' field - value should be between 0 and %u, but found %"PRIu64, field_name, MAX_POS, this_pos_64);
    
    return (int32_t)this_pos_64;
}

int32_t seg_pos_field (VBlock *vb, int32_t last_pos, int did_i, SectionType sec_pos_b250,
                       const char *pos_str, unsigned pos_len, const char *field_name)
{
    int32_t this_pos = seg_pos_snip_to_int (vb, pos_str, field_name);

    int32_t pos_delta = this_pos - last_pos;
    
    // print our string without expensive sprintf
    char pos_delta_str[50], reverse_pos_delta_str[50];
    
    bool negative = (pos_delta < 0);
    if (negative) pos_delta = -pos_delta;

    // create reverse string
    unsigned len=0; 
    if (pos_delta) { 
        while (pos_delta) {
            reverse_pos_delta_str[len++] = '0' + (pos_delta % 10);
            pos_delta /= 10;
        }
        if (negative) reverse_pos_delta_str[len++] = '-';

        // reverse it
        for (unsigned i=0; i < len; i++) pos_delta_str[i] = reverse_pos_delta_str[len-i-1];
    }
    else { //pos_delta==0
        pos_delta_str[0] = '0';
        len = 1;
    }

    seg_one_snip (vb, pos_delta_str, len, did_i, sec_pos_b250, NULL);

    vb->txt_section_bytes[sec_pos_b250] += pos_len - len; // re-do the calculation - seg_one_field doesn't do it good in our case

    return this_pos;
}

// Commonly (but not always), IDs are SNPid identifiers like "rs17030902". We store the ID divided to 2:
// - We store the final digits, if any exist, and up to 9 digits, as an integer in SEC_ID_DATA, which is
//   compressed with LZMA
// - In the dictionary we store the prefix up to this number, and \1 if there is a number and a \2 
//   if the caller (as in seg_me23) wants us to store an extra bit.
// example: rs17030902 : in the dictionary we store "rs\1" or "rs\1\2" and in SEC_ID_DATA we store 17030902.
//          1423       : in the dictionary we store "\1" and 1423 SEC_ID_DATA
//          abcd       : in the dictionary we store "abcd" and nothing is stored SEC_ID_DATA
void seg_id_field (VBlock *vb, Buffer *id_buf, int id_field, char *id_snip, unsigned id_snip_len, bool extra_bit)
{
    int i=id_snip_len-1; for (; i >= 0; i--) 
        if (!IS_DIGIT (id_snip[i])) break;
    
    unsigned num_digits = MIN (id_snip_len - (i+1), 9);

    // added to SEC_ID_DATA if we have a trailing number
    if (num_digits) {
        uint32_t id_num = BGEN32 (atoi (&id_snip[id_snip_len - num_digits]));
        seg_add_to_data_buf (vb, id_buf, SEC_ID_DATA, (char*)&id_num, sizeof (id_num), 0, num_digits); // account for the digits
    }

    // append the textual part with \1 and \2 as needed - we have enough space - we have a tab following the field
    // and if we need to add \1 we also know that we have at least one digit
    unsigned new_len = id_snip_len - num_digits;
    if (num_digits) id_snip[new_len++] = 1;
    if (extra_bit)  id_snip[new_len++] = 2;

    seg_one_field (vb, id_snip, new_len, id_field); // acounts for the length and \t
    vb->txt_section_bytes[FIELD_TO_B250_SECTION (vb->data_type, id_field)] -= (!!num_digits) + extra_bit; // we don't account for the \1 and \2 that were not in the txt file
}

// We break down the field (eg QNAME in SAM or Description in FASTA/FASTQ) into subfields separated by / and/or : -
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
                         char extra_separator, // a separator other than : / | (or 0 is there isn't one)
                         SectionType field_b250_sec, SectionType sf_b250_sec)
{
#define MAX_COMPOUND_COMPONENTS (10+26)

    const char *snip = field;
    unsigned snip_len = 0;
    unsigned sf_i = 0;
    char template[MAX_COMPOUND_COMPONENTS-1]; // separators is one less than the subfields
    MtfNode *node;

    // add each subfield to its dictionary - 2nd char is 0-9,a-z
    for (unsigned i=0; i <= field_len; i++) { // one more than field_len - to finalize the last subfield

        if (i==field_len || 
            ((field[i]==':' || field[i]=='/' || field[i]=='|' || field[i]==extra_separator) && sf_i < MAX_COMPOUND_COMPONENTS-1)) { // a subfield ended - separator between subfields
            
            // process the subfield that just ended
            MtfContext *sf_ctx;

            if (mapper->num_subfields == sf_i) { // new subfield in this VB (sf_ctx might exist from previous VBs)
                sf_dict_id.id[1] = (sf_i <= 9) ? (sf_i + '0') : (sf_i-10 + 'a');

                sf_ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, 
                                                 sf_dict_id, sf_b250_sec-1);
                mapper->did_i[sf_i] = sf_ctx->did_i;
                mapper->num_subfields++;
            }
            else 
                sf_ctx = MAPPER_CTX (mapper, sf_i);

            ASSERT0 (sf_ctx, "Error in seg_compound_field: sf_ctx is NULL");

            // allocate memory if needed
            buf_alloc (vb, &sf_ctx->mtf_i, MAX (vb->lines.len, sf_ctx->mtf_i.len + 1) * sizeof (uint32_t),
                       CTX_GROWTH, "mtf_ctx->mtf_i", sf_ctx->did_i);

            NEXTENT (uint32_t, sf_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, sf_ctx, snip, snip_len, &node, NULL);

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
    NEXTENT (uint32_t, field_ctx->mtf_i) = mtf_evaluate_snip_seg ((VBlockP)vb, field_ctx, template, MAX (1, sf_i-1), &node, NULL);

    // byte counts for --show-sections 
    vb->txt_section_bytes[field_b250_sec] += sf_i; // sf_i has 1 for each separator including the terminating \t or \n
    vb->txt_section_bytes[sf_b250_sec]    += field_len - (sf_i-1); // the entire field except for the / and : separators
}

void seg_add_to_data_buf (VBlock *vb, Buffer *buf, SectionType sec, 
                          const char *snip, unsigned snip_len, 
                          char add_separator,  // seperator to add to the buffer after the snip. 0 if none.
                          unsigned add_bytes)  // bytes in the original text file accounted for by this snip
{
    ASSERT0 (buf_is_allocated (buf), "Error in seg_add_to_data_buf: buf is not allocated");

    buf_alloc_more (vb, buf, snip_len + !!add_separator, 0, char, 2); // buffer must be pre-allocated before first call to seg_add_to_data_buf
    if (snip_len) buf_add (buf, snip, snip_len); 
    if (add_separator) buf_add (buf, &add_separator, 1); 
    vb->txt_section_bytes[sec] += add_bytes;
}

static void seg_set_hash_hints (VBlock *vb)
{
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {

        MtfContext *ctx = &vb->mtf_ctx[did_i];
        if (ctx->global_hash_prime) continue; // our service is not needed - global_cache for this dict already exists

        ctx->mtf_len_at_half = ctx->mtf.len;
        ctx->num_lines_at_half = vb->line_i + 1;
    }
}

static uint32_t seg_estimate_num_lines (VBlock *vb)
{
    // get first line length
    uint32_t len=0; for (; len < vb->txt_data.len && vb->txt_data.data[len] != '\n'; len++) {};

    ASSSEG0 (len < vb->txt_data.len, vb->txt_data.data, "Error: cannot find a newline in the entire vb");

    return MAX (100, (uint32_t)(((double)vb->txt_data.len / (double)len) * 1.2));
}

// split each lines in this variant block to its components
void seg_all_data_lines (VBlock *vb,
                         SegDataLineFuncType seg_data_line, 
                         SegInitializer seg_initialize, 
                         unsigned sizeof_line)
{
    START_TIMER;

    mtf_initialize_primary_field_ctxs (vb->mtf_ctx, vb->data_type, vb->dict_id_to_did_i_map, &vb->num_dict_ids); // Create ctx for the fields in the correct order 
    
    if (!sizeof_line) sizeof_line=1; // we waste a little bit of memory to avoid making exceptions throughout the code logic

    // allocate lines
    buf_alloc (vb, &vb->lines, seg_estimate_num_lines(vb) * sizeof_line, 1.2, "lines", vb->vblock_i);
    buf_zero (&vb->lines);
    vb->lines.len = vb->lines.size / sizeof_line;

    // allocate the mtf_i for the fields which each have num_lines entries
    for (int f=0; f <= datatype_last_field[vb->data_type]; f++) 
        buf_alloc (vb, &vb->mtf_ctx[f].mtf_i, vb->lines.len * sizeof (uint32_t), 1, "mtf_ctx->mtf_i", f);
    
    if (seg_initialize) seg_initialize (vb); // data-type specific initialization

    const char *field_start = vb->txt_data.data;
    bool hash_hints_set = false;
    for (vb->line_i=0; vb->line_i < vb->lines.len; vb->line_i++) {

        if (field_start - vb->txt_data.data == vb->txt_data.len) { // we're done
            vb->lines.len = vb->line_i; // update to actual number of lines
            break;
        }

        //fprintf (stderr, "vb->line_i=%u\n", vb->line_i);

        const char *next_field = seg_data_line (vb, field_start);
        
        vb->longest_line_len = MAX (vb->longest_line_len, (next_field - field_start));
        field_start = next_field;

        // if our estimate number of lines was too small, increase it
        if (vb->line_i == vb->lines.len-1 && field_start - vb->txt_data.data != vb->txt_data.len) {
            buf_alloc_more_zero (vb, &vb->lines, sizeof_line, 0, char, 1.5);
            vb->lines.len = vb->lines.size / sizeof_line;
        }

        // if there is no global_hash yet, and we've past half of the data,
        // collect stats to help mtf_merge create one when we merge
        if (!hash_hints_set && (field_start - vb->txt_data.data) > vb->txt_data.len / 2) {
            seg_set_hash_hints (vb);
            hash_hints_set = true;
        }
    }

    COPY_TIMER(vb->profile.seg_all_data_lines);
}