// ------------------------------------------------------------------
//   sam_zip.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "sam_private.h"
#include "seg.h"
#include "context.h"
#include "piz.h"
#include "strings.h"
#include "dict_id.h"
#include "reference.h"

// PIZ: SEQ contains : 
// '-' - data should be taken from the reference (-) 
// '.' - data should be taken from SQnonref.local (.)
// other - verbatim (this happens if there is no reference, eg unaligned BAM)
void sam_piz_reconstruct_seq (VBlock *vb_, Context *seq_ctx)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    ASSERT0 (seq_ctx, "Error in ref_reconstruct: seq_ctx is NULL");

    if (piz_is_skip_section (vb, SEC_LOCAL, seq_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx = mtf_get_ctx (vb, (DictId)dict_id_SAM_SQnonref);
    const char *nonref = &nonref_ctx->local.data[nonref_ctx->next_local]; // possibly, this VB has no nonref (i.e. everything is ref), in which cse nonref would be an invalid pointer. That's ok, as it will not be accessed.
    const char *nonref_start = nonref;
    const char *seq = &seq_ctx->local.data[seq_ctx->next_local];
    unsigned subcigar_len = 0;
    char cigar_op=0;
    
    // case where seq is '*' (rewritten as ' ' by the zip callback)
    if (seq[0] == ' ') {
        RECONSTRUCT1 ('*');
        seq_ctx->last_value = seq_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence
        seq_ctx->next_local++; // only 1
        return;
    }

    // get pointer to ref @ last chrom and pos. ref[0] is the base at POS.
    const char *ref = NULL;
    if (vb->last_cigar[0] != '*') 
        ref = ref_get_ref (vb_, (uint32_t)vb->contexts[SAM_POS].last_value, vb->ref_consumed);

    unsigned seq_consumed=0, ref_consumed=0;
    while (seq_consumed < vb->seq_len) {
        
        if (!subcigar_len) {
            subcigar_len = strtod (vb->last_cigar, (char **)&vb->last_cigar); // get number and advance next_cigar
            cigar_op = *(vb->last_cigar++);

            // case: Deletion or Skipping - skip some of the reference
            if (cigar_op == 'D' || cigar_op == 'N') {
                ref_consumed += subcigar_len;
                subcigar_len = 0;
                continue;
            } 

            // case: hard clipping - just ignore this subcigar
            if (cigar_op == 'H' || cigar_op == 'P') {
                subcigar_len = 0;
                continue;
            }
        }

        char c = seq[seq_consumed++];

        if (c == '-') {
            ASSERT0 (ref, "SEQ shows -, but ref is unavailable (eg bc entire chromosome is POS=0 or CIGAR=* or SEQ=*");
            RECONSTRUCT1 (ref[ref_consumed++]); 
        }

        else if (c == '.') {
            RECONSTRUCT1 (*nonref);
            nonref++; // consume non-ref even if not reconstructing this seq

            // advance ref if this is a SNP (but not in case of 'I' or 'S')
            if (cigar_op == 'M' || cigar_op == '=' || cigar_op == 'X') ref_consumed++; 
        }
        
        // in case of SAM having CIGAR='*' or POS=0 (eg unaligned BAM), we just copy verbatim
        else RECONSTRUCT1 (c); 


        subcigar_len--;
    }

    seq_ctx->last_value = seq_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence
    seq_ctx->next_local += vb->seq_len;
    
    nonref_ctx->next_local += (uint32_t)(nonref - nonref_start);
}

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
        vb->dont_show_curr_line = false; 

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
