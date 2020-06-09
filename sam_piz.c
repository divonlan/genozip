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
void sam_piz_special_CIGAR (VBlock *vb_, Context *ctx, const char *snip, unsigned snip_len)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    sam_analyze_cigar (snip, snip_len, &vb->seq_len, &vb->ref_consumed);

    if (snip[snip_len-1] == '*') { // change back eg "151*" -> "*"
        snip = &snip[snip_len-1];
        snip_len = 1;
    }

    RECONSTRUCT (snip, snip_len);
    vb->last_cigar = snip;
}   

void sam_piz_special_TLEN (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    ASSERT0 (snip_len, "Error in sam_piz_special_TLEN: snip_len=0");

    int32_t tlen_by_calc = atoi (snip);
    int32_t tlen_val = tlen_by_calc + vb->contexts[SAM_PNEXT].last_delta + vb->seq_len;

    ctx->last_value = tlen_val;

    RECONSTRUCT_INT (tlen_val);
}

void sam_piz_special_AS (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    RECONSTRUCT_INT (vb->seq_len - atoi (snip));
}

// logic: snip is eg "119C" (possibly also "") - we reconstruct the original, eg "119C31" 
// by concating a number which is (seq_len - partial_seq_len_by_md_field)
void sam_piz_special_MD (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    if (snip_len) RECONSTRUCT (snip, snip_len);

    unsigned partial_seq_len_by_md_field = sam_seg_get_seq_len_by_MD_field (snip, snip_len);
    RECONSTRUCT_INT (vb->seq_len - partial_seq_len_by_md_field);
}

// BI is a delta from the BD. Note: if BD doesn't appear on this line, then the snip is LOOKUP and not SPECIAL
// and we won't arrive at this function
void sam_piz_special_BI (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    ASSERT (ctx->next_local + snip_len <= ctx->local.len, "Error reading txt_line=%u: unexpected end of BI data", vb->line_i);

    Context *BD_ctx = mtf_get_ctx (vb, (DictId)dict_id_OPTION_BD);
    char *bi_dst = AFTERENT (char, vb->txt_data);
    const char *bd_txt = ENT (char, BD_ctx->local, (uint32_t)BD_ctx->last_value);
    const char *bi_txt = ENT (char, ctx->local, ctx->next_local);

    for (uint32_t i=0; i < vb->seq_len; i++) 
        bi_dst[i] = bi_txt[i] + bd_txt[i];

    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;
}

void sam_piz_reconstruct_vb (VBlockSAM *vb)
{
    piz_map_compound_field ((VBlockP)vb, sam_dict_id_is_qname_sf, &vb->qname_mapper);

    for (vb->line_i=vb->first_line; vb->line_i < vb->first_line + vb->lines.len; vb->line_i++) {

        uint32_t txt_data_start = vb->txt_data.len;

        piz_reconstruct_from_ctx (vb, SAM_QNAME,    '\t');
        piz_reconstruct_from_ctx (vb, SAM_FLAG,     '\t');
        piz_reconstruct_from_ctx (vb, SAM_RNAME,    '\t');
        piz_reconstruct_from_ctx (vb, SAM_POS,      '\t');
        piz_reconstruct_from_ctx (vb, SAM_MAPQ,     '\t'); 
        piz_reconstruct_from_ctx (vb, SAM_CIGAR,    '\t');
        piz_reconstruct_from_ctx (vb, SAM_RNEXT,    '\t'); 
        piz_reconstruct_from_ctx (vb, SAM_PNEXT,    '\t');
        piz_reconstruct_from_ctx (vb, SAM_TLEN,     '\t');
        piz_reconstruct_from_ctx (vb, SAM_SEQ,      '\t');
        piz_reconstruct_from_ctx (vb, SAM_QUAL,     '\t');
        piz_reconstruct_from_ctx (vb, SAM_OPTIONAL, 0   ); // the optional subfields (if there are any) provide the \t separators

        vb->txt_data.len--; // remove last \t (the line has ended either after QUAL or after the last OPTIONAL subfield)
        piz_reconstruct_from_ctx (vb, SAM_EOL, 0);

        // after consuming the line's data, if it is not to be outputted - trim txt_data back to start of line
        if (vb->dont_show_curr_line) vb->txt_data.len = txt_data_start; 
    }
}
 