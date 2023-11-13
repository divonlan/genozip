// ------------------------------------------------------------------
//   sam_pacbio.c
//   Copyright (C) 2023-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "reconstruct.h"
#include "piz.h"

void sam_pacbio_seg_initialize (VBlockSAMP vb)
{
    ctx_set_no_stons (VB, OPTION_dt_Z, OPTION_mq_Z, OPTION_st_Z, DID_EOL);

    ctx_set_ltype (VB, LT_DYN_INT, OPTION_qs_i, OPTION_qe_i, OPTION_np_i, DID_EOL);
    
    ctx_set_store (VB, STORE_INT, OPTION_qs_i, OPTION_qe_i, OPTION_zm_i, OPTION_np_i, DID_EOL);

    ctx_set_store (VB, STORE_FLOAT, OPTION_ec_f, DID_EOL);

    if (segconf.has[OPTION_ec_f]) 
        CTX(OPTION_np_i)->flags.same_line = true;  // np segged as delta vs ec, and np needs to be peeked for QUAL, before reconstructing ec

    if (segconf.use_pacbio_iqsqdq) {
        ctx_set_ltype (VB, LT_SEQUENCE, OPTION_dq_Z, OPTION_iq_Z, OPTION_sq_Z, DID_EOL);        
        ctx_consolidate_stats (VB, OPTION_iq_sq_dq, OPTION_dq_Z, OPTION_iq_Z, OPTION_sq_Z, DID_EOL);
    }
}

// -------------
// zm:i
// -------------

void sam_seg_pacbio_zm (VBlockSAMP vb, int64_t zm, unsigned add_bytes)    
{                               
    if (zm == CTX(SAM_Q1NAME)->last_value.i) 
        seg_by_ctx (VB, STRa(copy_Q1NAME_int), CTX(OPTION_zm_i), add_bytes);

    else 
        seg_integer (VB, CTX(OPTION_zm_i), zm, true, add_bytes);
}

// -------------
// np:i
// -------------

void sam_seg_pacbio_np (VBlockSAMP vb, ZipDataLineSAM *dl, int64_t np, unsigned add_bytes)    
{                     
    dl->np = np;

    if (segconf.has[OPTION_ec_f] && ctx_has_value_in_line_(vb, CTX(OPTION_ec_f)))
        seg_delta_vs_other_do (VB, CTX(OPTION_np_i), CTX(OPTION_ec_f), 0, 0, np, -1, add_bytes);

    else 
        seg_integer (VB, CTX(OPTION_np_i), np, true, add_bytes);
}

int32_t sam_zip_get_np (VBlockP vb, LineIType line_i)
{
    return DATA_LINE(line_i)->np;
}

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
               char, CTX_GROWTH, CTX_TAG_LOCAL);

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
    ContextP il_ctx   = CTX(OPTION_iq_sq_dq);

    ASSISLOADED(il_ctx);
    ASSISLOADED(diff_ctx);
    
    ASSPIZ (diff_ctx->next_local + vb->seq_len <= diff_ctx->local.len32, "SAM_QUAL_PACBIO_DIFF.local exhausted: next_local=%u + seq_len=%u > len=%u",
            diff_ctx->next_local, vb->seq_len, diff_ctx->local.len32);

    ASSPIZ (il_ctx->next_local + vb->seq_len * 3 <= il_ctx->local.len32, "OPTION_iq_sq_dq.local exhausted: next_local=%u + 3*seq_len=%u > len=%u",
            il_ctx->next_local, 3*vb->seq_len, il_ctx->local.len32);

    if (reconstruct) {
        rom diff = Bc(diff_ctx->local, diff_ctx->next_local);
        rom il   = Bc(il_ctx->local, il_ctx->next_local); // interlaced data

        char *qual = BAFTtxt, *after = BAFTtxt + vb->seq_len;

        while (qual < after) { 
            *qual++ = (*diff++) + prediction (il[0], il[1] - il[0], -il[2]);
            il += 3;
        }

        Ltxt += vb->seq_len;
    }

    diff_ctx->next_local += vb->seq_len;

    // note: OPTION_iq_sq_dq.next_local is used for reconstructing QUAL, while OPTION_[idx]q_Z.next_local are used for iq/dq/sq
    il_ctx->next_local   += vb->seq_len * 3; 
}

// advance il_ctx->next_local in case this line's QUAL was NOT segged with iq_sq_dq (eg it was monochar)
void sam_recon_skip_pacbio_qual (VBlockSAMP vb)
{
    if (container_peek_has_item (VB, CTX(SAM_AUX), _OPTION_iq_Z, false))
        CTX(OPTION_iq_sq_dq)->next_local += vb->seq_len * 3; 
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

    BufferP il = &CTX(OPTION_iq_sq_dq)->interlaced;
    buf_alloc_exact (vb, *il, dl->SEQ.len * 3, uint8_t, "interlaced");

    uint8_t *next = B1ST8 (*il);
    for (uint32_t i=0; i < dl->SEQ.len; i++) {
        *next++ = iq[i];
        *next++ = iq[i] + sq[i];
        *next++ = -dq[i];
    }

    *line_data = B1STc (*il);
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

// SAM + FASTQ: callback defined in qname_flavor for QF_PACBIO_rng 
void seg_qname_rng2seq_len_cb (VBlockP vb, ContextP ctx, STRp(value))
{
    seg_by_ctx (vb, (char[]){ SNIP_SPECIAL, (VB_DT(FASTQ) ? FASTQ_SPECIAL_qname_rng2seq_len : SAM_SPECIAL_qname_rng2seq_len) }, 2, ctx, 0);

    ctx_set_last_value (vb, ctx, (ctx-1)->last_value.i - (ctx-2)->last_value.i);
}

// SAM + FASTQ: since item has CI0_INVISIBLE, this function is always called with reconstruct=false
SPECIAL_RECONSTRUCTOR (special_qname_rng2seq_len)
{
    new_value->i = (ctx-1)->last_value.i - (ctx-2)->last_value.i;

    return HAS_NEW_VALUE;
}