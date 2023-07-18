// ------------------------------------------------------------------
//   sam_pacbio.c
//   Copyright (C) 2023-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "reconstruct.h"

//-------------------------------
// QUAL - predicted by iq, sq, dq
//-------------------------------

static inline char prediction (char iq, char sq, char dq)
{
    // create sorted array - bubble sort in ascending order
    char s[3] = { iq, sq, dq };
    if (s[0] > s[1]) SWAP (s[0], s[1]);
    if (s[1] > s[2]) SWAP (s[1], s[2]);
    if (s[0] > s[1]) SWAP (s[0], s[1]);    

    char prediction = (s[0] == '~')     ? '~'
                    : (s[1] - s[0] < 4) ? MAX_('!', s[0] - 2)
                    :                     s[0];

    return prediction;
}

bool sam_seg_pacbio_qual (VBlockSAMP vb, STRp(qual)/*textual*/, unsigned add_bytes)
{
    ContextP diff_ctx = CTX(SAM_QUAL_PACBIO_DIFF);

    ASSINP (has(dq_Z) && has(iq_Z) && has(sq_Z), "%s: Expecting line to have dq:Z, iq:Z and sq:Z but some are missing", LN_NAME);

    STR(iq); sam_seg_get_aux_Z (vb, vb->idx_iq_Z, pSTRa(iq), IS_BAM_ZIP);
    STR(sq); sam_seg_get_aux_Z (vb, vb->idx_sq_Z, pSTRa(sq), IS_BAM_ZIP);
    STR(dq); sam_seg_get_aux_Z (vb, vb->idx_dq_Z, pSTRa(dq), IS_BAM_ZIP);

    if (dq_len != qual_len || iq_len != qual_len || sq_len != qual_len) return false;

    buf_alloc (VB, &diff_ctx->local, qual_len, MIN_(vb->lines.len * qual_len, Ltxt / 5/*SEQ+QUAL+dq+iq+sq*/), 
               char, CTX_GROWTH, "local");

    char *next = BAFTc(diff_ctx->local);

    for (int i=0; i < qual_len; i++) 
        *next++ = qual[i] - prediction (iq[i], sq[i], dq[i]);
    
    diff_ctx->local.len32 = BNUM (diff_ctx->local, next);
    diff_ctx->txt_len += add_bytes;

    return true;
}

// callback function for compress to get data of one line
#define PACBIO_QV(f)                                    \
COMPRESSOR_CALLBACK (sam_zip_##f)                       \
{                                                       \
    TxtWord w = DATA_LINE (vb_line_i)->f;               \
                                                        \
    *line_data_len = MIN_(maximum_size, w.len);         \
                                                        \
    if (!line_data) return; /* only lengths requested */\
                                                        \
    *line_data = Btxt (w.index);                        \
}                                                         
PACBIO_QV(dq)
PACBIO_QV(iq)
PACBIO_QV(sq)

void sam_recon_pacbio_qual (VBlockSAMP vb, ContextP ctx, bool reconstruct)
{
    ContextP diff_ctx = CTX(SAM_QUAL_PACBIO_DIFF);

    if (!reconstruct) return;

    ASSPIZ (diff_ctx->next_local + vb->seq_len <= diff_ctx->local.len32, "SAM_QUAL_PACBIO_DIFF.local exhausted: next_local=%u + seq_len=%u > len=%u",
            diff_ctx->next_local, vb->seq_len, diff_ctx->local.len32);

    rom diff = Bc(diff_ctx->local, diff_ctx->next_local);
    rom il   = Bc(CTX(OPTION_iq_sq_dq)->local, CTX(OPTION_iq_Z)->next_local * 3); // interlaced data

    char *qual = BAFTtxt;

    for (int i=0; i < vb->seq_len; i++) { 
        *qual++ = (*diff++) + prediction (il[0], il[1] - il[0], -il[2]);
        il += 3;
    }

    Ltxt += vb->seq_len;
    diff_ctx->next_local += vb->seq_len;
}


// -------------
// iq, sq and dq
// -------------

void sam_seg_pacbio_xq (VBlockSAMP vb, ZipDataLineSAM *dl, Did did_i, TxtWord *dl_word, STRp(value), unsigned add_bytes)    
{                                                          
    decl_ctx (did_i);

    ASSINP (value_len == dl->SEQ.len, "%s: Expecting %s.len=%u == SEQ.len=%u. %s=\"%.*s\"",       
            LN_NAME, ctx->tag_name, value_len, dl->SEQ.len, ctx->tag_name, STRf(value));     

    *dl_word = TXTWORD(value);                                                           
    CTX(OPTION_iq_sq_dq)->txt_len += add_bytes;                      
    CTX(OPTION_iq_sq_dq)->local.len32 += value_len;                  

    seg_by_ctx (VB, ((char[]){ SNIP_SPECIAL, SAM_SPECIAL_iqsqdq }), 2, ctx, 0);
}                                                          

// callback function for compress to get iq_sq_dq data of one line: this is an
// interlaced line containing a character from iq followed by a character from sq followed by dq 
COMPRESSOR_CALLBACK (sam_zip_iq_sq_dq)
{
    ZipDataLineSAM *dl = DATA_LINE (vb_line_i);
    
    rom iq = dl->iq.index ? Btxt (dl->iq.index) : NULL;
    rom sq = dl->sq.index ? Btxt (dl->sq.index) : NULL;
    rom dq = dl->dq.index ? Btxt (dl->dq.index) : NULL;
    
    if (!iq && !sq && !dq) return; // no iq/sq/dq on this line

    ASSERT (iq && sq && dq, "%s/%u: line has one or two of the iq:Z/sq:Z/dq:Z triplet - Genozip can only compress lines that have either all these or none", VB_NAME, vb_line_i); 
    
    // note: maximum_len might be shorter than the data available if we're just sampling data in codec_assign_best_codec
    *line_data_len  = MIN_(maximum_size, dl->SEQ.len * 3);
    if (is_rev) *is_rev = dl->FLAG.rev_comp;

    if (!line_data) return; // only length was requested

    buf_alloc_exact (vb, VB_SAM->interlaced, dl->SEQ.len * 3, uint8_t, "interlaced");

    uint8_t *next = B1ST8 (VB_SAM->interlaced);
    for (uint32_t i=0; i < dl->SEQ.len; i++) {
        *next++ = iq[i];
        *next++ = iq[i] + sq[i];
        *next++ = -dq[i];
    }

    *line_data = B1STc (VB_SAM->interlaced);
}   

// iq, sq, dq - reconstruct from iq_sq_dq context which contains interlaced data. 
SPECIAL_RECONSTRUCTOR (sam_piz_special_iq_sq_dq)
{
    if (!vb->seq_len || !reconstruct) goto done;

    Context *il_ctx = CTX(OPTION_iq_sq_dq); // interlaced data

    // note: iq,sq,dq use their own next_local to retrieve data from il_ctx. the actual index
    // in il_ctx.local is calculated given the interlacing
    ASSPIZ (ctx->next_local + vb->seq_len * 3 <= il_ctx->local.len, "Unexpected end of %s data. Expecting ctx->next_local=%u + vb->seq_len=%u * 3 <= il_ctx->local.len=%u", 
            il_ctx->tag_name, ctx->next_local, vb->seq_len, il_ctx->local.len32);

    char *dst        = BAFTtxt;
    rom src          = Bc (il_ctx->local, ctx->next_local * 3);
    uint32_t seq_len = vb->seq_len; // automatic var for effeciency

    if (ctx->did_i == OPTION_iq_Z)
        for (uint32_t i=0; i < seq_len; i++, src+=3, dst++) 
            *dst = src[0];
    
    else if (ctx->did_i == OPTION_sq_Z) 
        for (uint32_t i=0; i < seq_len; i++, src+=3, dst++) 
            *dst = src[1] - src[0];

    else 
        for (uint32_t i=0; i < seq_len; i++, src+=3, dst++) 
            *dst = -src[2];

    Ltxt += vb->seq_len;    
    ctx->next_local  += vb->seq_len;

done:
    return NO_NEW_VALUE;
}
