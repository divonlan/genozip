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
#include "codec.h" // must be included before reference.h
#include "reference.h"
#include "regions.h"
#include "aligner.h"
#include "file.h"

// returns true if section is to be skipped reading / uncompressing
bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT sections

    return false;
}

// PIZ: SEQ reconstruction 
void sam_piz_reconstruct_seq (VBlock *vb_, Context *bitmap_ctx, const char *unused, unsigned unused2)
{
#define ROUNDUP_TO_NEAREST_4(x) ((uint32_t)(x) + 3) & ~((uint32_t)0x3)

    VBlockSAMP vb = (VBlockSAMP)vb_;
    ASSERT0 (bitmap_ctx && bitmap_ctx->did_i == SAM_SEQ_BITMAP, "Error in sam_piz_reconstruct_seq: context is not SAM_SEQ_BITMAP");

    if (piz_is_skip_section (vb, SEC_LOCAL, bitmap_ctx->dict_id)) return; // if case we need to skip the SEQ field (for the entire file)

    Context *nonref_ctx      = &vb->contexts[SAM_NONREF];
    const char *nonref       = ENT (const char, nonref_ctx->local, nonref_ctx->next_local); // possibly, this VB has no nonref (i.e. everything is ref), in which case nonref would be an invalid pointer. That's ok, as it will not be accessed.
    const char *nonref_start = nonref;
    unsigned subcigar_len    = 0;
    char cigar_op            = 0;
    const PosType pos        = vb->contexts[SAM_POS].last_value.i;
    const Range *range       = NULL;
    unsigned seq_consumed=0, ref_consumed=0;

    // case: unaligned sequence - pos is 0 
    if (!pos) {
        // case: compressed with a reference, using our aligner
        if (z_file->flags & GENOZIP_FL_ALIGNER) {
            aligner_reconstruct_seq ((VBlockP)vb, bitmap_ctx, vb->seq_len, false);
            nonref_ctx->next_local = ROUNDUP_TO_NEAREST_4 (nonref_ctx->next_local);
        }
        // case: no reference was used - in this case, the sequence is not encoded in the bitmap at all. we just copy it from SEQ_NOREF
        else {
            RECONSTRUCT (nonref, vb->seq_len); 
            nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (vb->seq_len);
        }
        return;
    }

    // case: missing sequence - sequence is '*' (which zip marked with a '-' in the cigar) - just reconstruct the '*'
    if (*vb->last_cigar == '-') {
        RECONSTRUCT1 ('*');
        vb->last_cigar++; // skip the '-' so it doesn't affect a subsequent E2 on the same line
        return;
    }

    const char *next_cigar = vb->last_cigar; // don't change vb->last_cigar as we may still need it, eg if we have an E2 optional field
    range = ref_piz_get_range (vb_, pos, vb->ref_consumed);
    
    while (seq_consumed < vb->seq_len || ref_consumed < vb->ref_consumed) {
        
        if (!subcigar_len) {
            subcigar_len = strtod (next_cigar, (char **)&next_cigar); // get number and advance next_cigar
        
            cigar_op = cigar_lookup[(uint8_t)*(next_cigar++)];
            ASSERT (cigar_op != CIGAR_INVALID, "Error in sam_piz_reconstruct_seq: Invalid CIGAR op while reconstructing line %u: '%c' (ASCII %u)", vb->line_i, *(next_cigar-1), *(next_cigar-1));
        }

        if (cigar_op & CIGAR_CONSUMES_QUERY) {

            if ((cigar_op & CIGAR_CONSUMES_REFERENCE) && NEXTLOCALBIT (bitmap_ctx)) /* copy from reference */ {

                if (!vb->dont_show_curr_line) { // note: if this line is excluded with --regions, then the reference section covering it might not be loaded
                    uint32_t idx = (pos - range->first_pos) + ref_consumed ;

                    if (!ref_is_nucleotide_set (range, idx)) { 
                        ref_print_is_set (range, pos + ref_consumed);
                        ABORT ("Error in sam_piz_reconstruct_seq: while reconstructing line %u (vb_i=%u: last_txt_line=%u num_lines=%u): reference is not set: chrom=%u \"%.*s\" pos=%"PRId64" range=[%"PRId64"-%"PRId64"]"
                            " (cigar=%s seq_start_pos=%"PRId64" ref_consumed=%u seq_consumed=%u)",
                            vb->line_i, vb->vblock_i, (uint32_t)(vb->first_line + vb->lines.len - 1), (uint32_t)vb->lines.len, range->chrom, range->chrom_name_len, range->chrom_name, pos + ref_consumed, range->first_pos, range->last_pos, vb->last_cigar, pos, ref_consumed, seq_consumed);
                    }

                    char ref = ref_get_nucleotide (range, idx);
                    RECONSTRUCT1 (ref); 
                }
            }
            else 
                RECONSTRUCT1 (*nonref++);

            seq_consumed++;
        }

        if (cigar_op & CIGAR_CONSUMES_REFERENCE) 
            ref_consumed++;

        subcigar_len--;
    }

    ASSERT (seq_consumed == vb->seq_len,      "Error in sam_piz_reconstruct_seq: expecting seq_consumed(%u) == vb->seq_len(%u)", seq_consumed, vb->seq_len);
    ASSERT (ref_consumed == vb->ref_consumed, "Error in sam_piz_reconstruct_seq: expecting ref_consumed(%u) == vb->ref_consumed(%u)", ref_consumed, vb->ref_consumed);

    bitmap_ctx->last_value.i = bitmap_ctx->next_local; // for SEQ, we use last_value for storing the beginning of the sequence
    
    nonref_ctx->next_local += ROUNDUP_TO_NEAREST_4 (nonref - nonref_start);
}

// CIGAR - calculate vb->seq_len from the CIGAR string, and if original CIGAR was "*" - recover it
void sam_piz_special_CIGAR (VBlock *vb_, Context *ctx, const char *snip, unsigned snip_len)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;

    sam_analyze_cigar (snip, snip_len, &vb->seq_len, &vb->ref_consumed, NULL);

    if (snip[snip_len-1] == '*') // eg "151*" - zip added the "151" to indicate seq_len - we don't reconstruct it, just the '*'
        RECONSTRUCT1 ('*');
    
    else if (snip[0] == '-') // eg "-151M" or "-151*" - zip added the "-" to indicate a '*' SEQ field - we don't reconstruct it
        RECONSTRUCT (snip + 1, snip_len - 1)

    else
        RECONSTRUCT (snip, snip_len);    

    vb->last_cigar = snip;

    if (flag_regions && vb->chrom_node_index != WORD_INDEX_NONE && vb->contexts[SAM_POS].last_value.i && 
        !regions_is_range_included (vb->chrom_node_index, vb->contexts[SAM_POS].last_value.i, vb->contexts[SAM_POS].last_value.i + vb->ref_consumed - 1, true))
        vb->dont_show_curr_line = true;
}   

void sam_piz_special_TLEN (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    ASSERT0 (snip_len, "Error in sam_piz_special_TLEN: snip_len=0");

    int32_t tlen_by_calc = atoi (snip);
    int32_t tlen_val = tlen_by_calc + vb->contexts[SAM_PNEXT].last_delta + vb->seq_len;

    ctx->last_value.i = tlen_val;

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

// BD and BI - reconstruct from BD_BI context which contains interlaced BD and BI data. 
void sam_piz_special_BD_BI (VBlock *vb, Context *ctx, const char *snip, unsigned snip_len)
{
    if (!vb->seq_len) return;

    Context *bdbi_ctx = mtf_get_existing_ctx (vb, dict_id_OPTION_BD_BI);

    // note: bd and bi use their own next_local to retrieve data from bdbi_ctx. the actual index
    // in bdbi_ctx.local is calculated given the interlacing
    ASSERT (ctx->next_local + vb->seq_len * 2 <= bdbi_ctx->local.len, "Error reading txt_line=%u: unexpected end of %s data", vb->line_i, err_dict_id (ctx->dict_id));

    char *dst        = AFTERENT (char, vb->txt_data);
    const char *src  = ENT (char, bdbi_ctx->local, ctx->next_local * 2);
    uint32_t seq_len = vb->seq_len; // automatic var for effeciency

    if (ctx->dict_id.num == dict_id_OPTION_BD)
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src;
    else
        for (uint32_t i=0; i < seq_len; i++, src+=2, dst++) *dst = *src + *(src+1);
    
    vb->txt_data.len += vb->seq_len;    
    ctx->next_local  += vb->seq_len;
}
