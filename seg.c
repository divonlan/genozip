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
uint32_t seg_one_subfield (VBlock *vb, const char *str, unsigned len, unsigned vb_line_i,
                           DictIdType dict_id, SectionType sec_b250, int accounts_for_chars)
{
    MtfNode *node;
    MtfContext *ctx = mtf_get_ctx_by_dict_id (vb->mtf_ctx, &vb->num_dict_ids, NULL, dict_id, sec_b250 - 1);

    // allocate memory if needed
    buf_alloc (vb, &ctx->mtf_i, MAX (vb->num_lines, ctx->mtf_i.len + 1) * sizeof (uint32_t),
               CTX_GROWTH, "mtf_ctx->mtf_i", ctx->did_i);

    uint32_t node_index = mtf_evaluate_snip_seg (vb, ctx, str, len, &node, NULL);

    ASSERT (node_index < ctx->mtf.len + ctx->ol_mtf.len || node_index == WORD_INDEX_EMPTY_SF, 
            "Error in seg_one_subfield: out of range: dict=%.*s %s mtf_i=%d mtf.len=%u ol_mtf.len=%u",  
            DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, st_name (ctx->dict_section_type),
            node_index, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len);

    NEXTENT (uint32_t, ctx->mtf_i) = node_index;

    vb->txt_section_bytes[sec_b250] += accounts_for_chars; 

    return node_index;
}

// returns the node index
uint32_t seg_one_field (VBlock *vb, const char *str, unsigned len, unsigned vb_line_i, int f, SectionType sec_b250,
                        bool *is_new) // optional out
{
    buf_alloc (vb, &vb->mtf_ctx[f].mtf_i, (vb->mtf_ctx[f].mtf_i.len + 1) * sizeof (uint32_t), 2, "mtf_ctx->mtf_i", f);
    
    MtfContext *ctx = &vb->mtf_ctx[f];
    MtfNode *node;
    uint32_t node_index = mtf_evaluate_snip_seg ((VBlockP)vb, ctx, str, len, &node, is_new);

    ASSERT (node_index < ctx->mtf.len + ctx->ol_mtf.len || node_index == WORD_INDEX_EMPTY_SF, "Error in seg_one_field: out of range: dict=%.*s %s mtf_i=%d mtf.len=%u ol_mtf.len=%u",  
            DICT_ID_LEN, dict_id_printable (ctx->dict_id).id, st_name (ctx->dict_section_type),
            node_index, (uint32_t)ctx->mtf.len, (uint32_t)ctx->ol_mtf.len);
    
    NEXTENT (uint32_t, vb->mtf_ctx[f].mtf_i) = node_index;

    vb->txt_section_bytes[sec_b250] += len + 1;
    return node_index;
} 

const char *seg_get_next_item (const char *str, int *str_len, bool allow_newline, bool allow_tab, bool allow_colon, unsigned vb_line_i, // line in vcf file,
                               unsigned *len, char *separator, bool *has_13, // out
                               const char *item_name)
{
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
            
    ASSERT (*str_len, "Error: missing %s field in vb_line_i=%u", item_name, vb_line_i);

    ASSERT (str[i] != '\t' || txt_file->data_type != DATA_TYPE_VCF || strcmp (item_name, vcf_field_names[VCF_INFO]), 
           "Error: while segmenting %s in vb_line_i=%u: expecting a NEWLINE after the INFO field, because this VCF file has no samples (individuals) declared in the header line",
            item_name, vb_line_i);

    ABORT ("Error: while segmenting %s in vb_line_i=%u: expecting a %s %s %s after \"%.*s\"", 
            item_name, vb_line_i, 
            allow_newline ? "NEWLINE" : "", allow_tab ? "TAB" : "", allow_colon ? "\":\"" : "", 
            MIN (i-1, 1000), str);

    return 0; // avoid compiler warning - never reaches here
}

#define MAX_POS 0x7fffffff // maximum allowed value for POS

// reads a tab-terminated POS string
int32_t seg_pos_snip_to_int (const char *pos_str, unsigned vb_line_i)
{
    // scan by ourselves - hugely faster the sscanf
    int64_t this_pos_64=0; // int64_t so we can test for overflow
    const char *s; for (s=pos_str; *s != '\t'; s++) {
        ASSERT (*s >= '0' && *s <= '9', "Error: POS field in vb_line_i=%u must be an integer number between 0 and %u", vb_line_i, MAX_POS);
        this_pos_64 = this_pos_64 * 10 + (*s - '0');
    }
    ASSERT (this_pos_64 >= 0 && this_pos_64 <= 0x7fffffff, 
            "Error: Invalid POS in vb_line_i=%u - value should be between 0 and %u, but found %"PRIu64, vb_line_i, MAX_POS, this_pos_64);
    
    return (int32_t)this_pos_64;
}

int32_t seg_pos_field (VBlock *vb, int32_t last_pos, int pos_field, SectionType sec_pos_b250,
                       const char *pos_str, unsigned pos_len, unsigned vb_line_i)
{
    int32_t this_pos = seg_pos_snip_to_int (pos_str, vb_line_i);

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

    seg_one_field (vb, pos_delta_str, len, vb_line_i, pos_field, sec_pos_b250, NULL);

    vb->txt_section_bytes[sec_pos_b250] += pos_len - len; // re-do the calculation - seg_vcf_one_field doesn't do it good in our case

    return this_pos;
}

static void seg_set_hash_hints (VBlock *vb, uint32_t vb_line_i)
{
    for (unsigned did_i=0; did_i < vb->num_dict_ids; did_i++) {

        MtfContext *ctx = &vb->mtf_ctx[did_i];
        if (ctx->global_hash_prime) continue; // our service is not needed - global_cache for this dict already exists

        ctx->mtf_len_at_half = ctx->mtf.len;
        ctx->num_lines_at_half = vb_line_i + 1;
    }
}

static void seg_realloc_datalines (VBlock *vb, uint32_t new_num_data_lines, unsigned sizeof_line)
{
    vb->data_lines = REALLOC (vb->data_lines, new_num_data_lines * sizeof_line);
    
    memset (&((char*)vb->data_lines)[vb->num_lines_alloced * sizeof_line], 0, 
            (new_num_data_lines - vb->num_lines_alloced) * sizeof_line);
    
    vb->num_lines_alloced = new_num_data_lines;
}

static void seg_allocate_per_line_memory (VBlock *vb, unsigned sizeof_line, int first_field, int last_field)
{
    ASSERT (!!vb->data_lines == !!vb->num_lines_alloced, 
            "Error: expecting vb->data_lines to be nonzero iff vb->num_lines_alloced is nonzero. vb_i=%u", vb->vblock_i);

    ASSERT (vb->num_lines_alloced >= vb->num_lines, 
            "Error: expecting vb->num_lines_alloced >= vb->num_lines. vb_i=%u", vb->vblock_i);

    // first line in this vb.id - we calculate an estimated number of lines
    if (!vb->num_lines) { 
        // get first line length
        uint32_t len=0; for (; len < vb->txt_data.len && vb->txt_data.data[len] != '\n'; len++) {};

        ASSERT (len < vb->txt_data.len, "Error: cannot from a newline in the entire vb. vb_i=%u", vb->vblock_i);

        // set initial number of lines based on an estimate. it might grow during segmentation if it turns out the
        // estimate is too low, and will be set to its actual real value at the end of segmentation
        
        uint32_t lower_end_estimate  = MAX (100, (uint32_t)((double)vb->txt_data.len / (double)len));
        uint32_t higher_end_estimate = MAX (100, (uint32_t)(((double)vb->txt_data.len / (double)len) * 1.2));

        // if we have enough according to the lower end of the estimate, go for it
        if (vb->num_lines_alloced >= lower_end_estimate) 
            vb->num_lines = vb->num_lines_alloced;
        
        else if (!vb->num_lines_alloced) {
            vb->num_lines = vb->num_lines_alloced = higher_end_estimate;
            vb->data_lines = calloc (vb->num_lines, sizeof_line);
        }
        // num_lines_alloced is non-zero, but too low - reallocate 
        else { 
            seg_realloc_datalines (vb, higher_end_estimate, sizeof_line); // uses and updates vb->num_lines_alloced       
            vb->num_lines = vb->num_lines_alloced;
        }
    }

    // first line in the VB, but allocation exists from previous reincarnation of this VB structure - keep
    // what's already allocated
    else if (vb->num_lines < vb->num_lines_alloced) 
        vb->num_lines = vb->num_lines_alloced;

    // we already have lines - but we need more. reallocate.
    else {
        seg_realloc_datalines (vb, vb->num_lines_alloced * 1.5, sizeof_line);        
        vb->num_lines = vb->num_lines_alloced;
    }

    // allocate (or realloc) the mtf_i for the fields which each have num_lines entries
    for (int f=first_field; f <= last_field; f++) 
        buf_alloc (vb, &vb->mtf_ctx[f].mtf_i, vb->num_lines * sizeof (uint32_t), 1, "mtf_ctx->mtf_i", f);
}

// split each lines in this variant block to its components
void seg_all_data_lines (VBlock *vb,
                         SegDataLineFuncType seg_data_line, 
                         unsigned sizeof_line,
                         int first_field, int last_field,
                         const char **field_names,
                         SectionType first_field_dict_section)
{
    START_TIMER;

    vb->num_dict_ids = MAX (last_field+1, vb->num_dict_ids); // first mtf_ctx are reserved for the field (vb->num_dict_ids might be already higher due to previous VBs)

    // Set ctx stuff for QNAME->OPTIONAL fields (note: mtf_i is allocated by seg_allocate_per_line_memory)
    for (int f=first_field; f <= last_field; f++) {
        vb->mtf_ctx[f].dict_id = dict_id_field (dict_id_make (field_names[f], strlen(field_names[f])));
        vb->mtf_ctx[f].b250_section_type = first_field_dict_section + f*2 + 1; // b250 is always one after its dict
        vb->mtf_ctx[f].dict_section_type = first_field_dict_section + f*2;
    }

    seg_allocate_per_line_memory (vb, sizeof_line, first_field, last_field); // set vb->num_lines to an initial estimate

    const char *field_start = vb->txt_data.data;
    bool hash_hints_set = false;
    for (unsigned vb_line_i=0; vb_line_i < vb->num_lines; vb_line_i++) {

        if (field_start - vb->txt_data.data == vb->txt_data.len) { // we're done
            vb->num_lines = vb_line_i; // update to actual number of lines
            break;
        }

        //fprintf (stderr, "vb_line_i=%u\n", vb_line_i);

        field_start = seg_data_line (vb, field_start, vb_line_i);

        // if our estimate number of lines was too small, increase it
        if (vb_line_i == vb->num_lines-1 && field_start - vb->txt_data.data != vb->txt_data.len) 
            seg_allocate_per_line_memory (vb, sizeof_line, first_field, last_field); // increase number of lines as evidently we need more

        // if there is no global_hash yet, and we've past half of the data,
        // collect stats to help mtf_merge create one when we merge
        if (!hash_hints_set && (field_start - vb->txt_data.data) > vb->txt_data.len / 2) {
            seg_set_hash_hints (vb, vb_line_i);
            hash_hints_set = true;
        }
    }

    COPY_TIMER(vb->profile.seg_all_data_lines);
}