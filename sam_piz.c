// ------------------------------------------------------------------
//   sam_zip.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "sam_private.h"
#include "seg.h"
#include "move_to_front.h"
#include "piz.h"
#include "strings.h"
#include "dict_id.h"

// CIGAR - calculate vb->seq_len from the CIGAR string, and if original CIGAR was "*" - recover it
void sam_piz_special_CIGAR (VBlock *vb, MtfContext *ctx, const char *snip, unsigned snip_len)
{
    vb->seq_len = sam_seq_len_from_cigar (snip, snip_len);

    if (snip[snip_len-1] == '*') 
        RECONSTRUCT1 ('*');
    else
        RECONSTRUCT (snip, snip_len);
}   

void sam_piz_special_TLEN (VBlock *vb, MtfContext *ctx, const char *snip, unsigned snip_len)
{
    ASSERT0 (snip_len, "Error in sam_piz_special_TLEN: snip_len=0");

    int32_t tlen_by_calc = atoi (snip);
    int32_t tlen_val = tlen_by_calc + vb->contexts[SAM_PNEXT].last_delta + vb->seq_len;

    ctx->last_value = tlen_val;

    RECONSTRUCT_INT (tlen_val);
}

void sam_piz_special_AS (VBlock *vb, MtfContext *ctx, const char *snip, unsigned snip_len)
{
    RECONSTRUCT_INT (vb->seq_len - atoi (snip));
}

// logic: snip is eg "119C" (possibly also "") - we reconstruct the original, eg "119C31" 
// by concating a number which is (seq_len - partial_seq_len_by_md_field)
void sam_piz_special_MD (VBlock *vb, MtfContext *ctx, const char *snip, unsigned snip_len)
{
    if (snip_len) RECONSTRUCT (snip, snip_len);

    unsigned partial_seq_len_by_md_field = sam_seg_get_seq_len_by_MD_field (snip, snip_len);
    RECONSTRUCT_INT (vb->seq_len - partial_seq_len_by_md_field);
}

// BI is a delta from the BD. Note: if BD doesn't appear on this line, then the snip is LOOKUP and not SPECIAL
// and we won't arrive at this function
void sam_piz_special_BI (VBlock *vb, MtfContext *ctx, const char *snip, unsigned snip_len)
{
    ASSERT (ctx->next_local + snip_len <= ctx->local.len, "Error reading txt_line=%u: unexpected end of BI data", vb->line_i);

    MtfContext *BD_ctx = mtf_get_ctx (vb, (DictIdType)dict_id_OPTION_BD);
    char *bi_dst = AFTERENT (char, vb->txt_data);
    const char *bd_txt = ENT (char, BD_ctx->local, (uint32_t)BD_ctx->last_value);
    const char *bi_txt = ENT (char, ctx->local, ctx->next_local);

    for (uint32_t i=0; i < vb->seq_len; i++) 
        bi_dst[i] = bi_txt[i] + bd_txt[i];

    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;
}

static void sam_piz_map_optional_subfields (VBlockSAM *vb)
{
    // terminology: we call a list of INFO subfield names, an "oname". An oname looks something like
    // this: "MX:Z:ab:i:". Each oname consists of OPTIONAL subfields. These subfields are not unique to this
    // oname and can appear in other onames. Each optional subfield is made of a segment of the template eg "MX:Z:"
    // and a value eg "abcded" which is stored in the b250 of that optional subfield.

    const MtfContext *optional_ctx = &vb->contexts[SAM_OPTIONAL];
    vb->optional_mapper_buf.len = optional_ctx->word_list.len;
    buf_alloc (vb, &vb->optional_mapper_buf, sizeof (SubfieldMapper) * vb->optional_mapper_buf.len,
               1, "optional_mapper_buf", 0);
    buf_zero (&vb->optional_mapper_buf);

    ARRAY (SubfieldMapper, all_optional_mappers, vb->optional_mapper_buf);

    ARRAY (const MtfWord, all_optionals, optional_ctx->word_list);

    for (unsigned optional_i=0; optional_i < vb->optional_mapper_buf.len; optional_i++) {

        char *optional = ENT (char, optional_ctx->dict, all_optionals[optional_i].char_index); // e.g. "MX:Z:ab:i:" - pointer into the OPTIONAL dictionary
        SubfieldMapper *optional_mapper = &all_optional_mappers[optional_i]; // optional_mapper of this specific set of names "I1=I2=I3="
        optional_mapper->num_subfields = all_optionals[optional_i].snip_len / 5; 

        // traverse the subfields of one optional. E.g. if the optional is "MX:Z:ab:i:" then we traverse MX, ab
        // note - optional 
        for (unsigned i=0; i < optional_mapper->num_subfields; i++) {
                        
            DictIdType dict_id = sam_dict_id_optnl_sf (dict_id_make (&optional[i*5], 4)); // get dict_id from first 4 characters eg "MX:i"
            optional_mapper->did_i[i] = mtf_get_ctx (vb, dict_id)->did_i; 
        }
    }
}

static void sam_piz_reconstruct_optional_fields (VBlockSAM *vb, const char *oname, unsigned oname_len, uint32_t opt_word_index)
{

    SubfieldMapper *opt_map = ENT (SubfieldMapper, vb->optional_mapper_buf, opt_word_index);

    ASSERT (opt_map->num_subfields == oname_len / 5, "Error: opt_map->num_subfields=%u but oname=%.*s indicates %u optional fields. sam_line=%u", 
            opt_map->num_subfields, oname_len, oname, oname_len/5, vb->line_i);

    for (unsigned sf_i=0; sf_i < opt_map->num_subfields; sf_i++) {

        RECONSTRUCT (&oname[sf_i*5], 5)

        uint8_t did_i = opt_map->did_i[sf_i];
        piz_reconstruct_from_ctx (vb, did_i, 0);        

        if (sf_i != opt_map->num_subfields-1)
            RECONSTRUCT1 ('\t');
    }
}

void sam_piz_reconstruct_vb (VBlockSAM *vb)
{
    piz_map_compound_field ((VBlockP)vb, sam_dict_id_is_qname_sf, &vb->qname_mapper);
    sam_piz_map_optional_subfields (vb);

    DECLARE_SNIP;
    for (uint32_t vb_line_i=0; vb_line_i < vb->lines.len; vb_line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;
        vb->line_i = vb->first_line + vb_line_i;

        piz_reconstruct_from_ctx (vb, SAM_QNAME, '\t');
        piz_reconstruct_from_ctx (vb, SAM_FLAG,  '\t');
        piz_reconstruct_from_ctx (vb, SAM_RNAME, '\t');
        piz_reconstruct_from_ctx (vb, SAM_POS,   '\t');
        piz_reconstruct_from_ctx (vb, SAM_MAPQ,  '\t'); 
        piz_reconstruct_from_ctx (vb, SAM_CIGAR, '\t');
        piz_reconstruct_from_ctx (vb, SAM_RNEXT, '\t'); 
        piz_reconstruct_from_ctx (vb, SAM_PNEXT, '\t');
        piz_reconstruct_from_ctx (vb, SAM_TLEN,  '\t');
        piz_reconstruct_from_ctx (vb, SAM_SEQ,   '\t');
        piz_reconstruct_from_ctx (vb, SAM_QUAL,     0);

        // OPTIONAL fields, and Windows-style \r if needed
        uint32_t word_index = LOAD_SNIP (SAM_OPTIONAL);
        if (snip_len != 1 || snip[0] != '*') RECONSTRUCT1 ('\t');
        
        sam_piz_reconstruct_optional_fields (vb, snip, snip_len, word_index);

        piz_reconstruct_from_ctx (vb, SAM_EOL, 0);

        // after consuming the line's data, if it is not to be outputted - trim txt_data back to start of line
        if (vb->dont_show_curr_line) vb->txt_data.len = txt_data_start; 
    }
}
