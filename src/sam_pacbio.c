// ------------------------------------------------------------------
//   sam_pacbio.c
//   Copyright (C) 2023-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

#include "sam_private.h"
#include "reconstruct.h"
#include "piz.h"
#include "qname.h"

void sam_pacbio_seg_initialize (VBlockSAMP vb)
{
    ctx_set_no_stons (VB, OPTION_dt_Z, OPTION_mq_Z, OPTION_st_Z, DID_EOL);

    ctx_set_dyn_int (VB, OPTION_qs_i, OPTION_qe_i, OPTION_np_i, OPTION_zm_i, DID_EOL);
    
    ctx_set_store (VB, STORE_FLOAT, OPTION_ec_f, DID_EOL);

    ctx_set_store (VB, STORE_INT, OPTION_ws_i,
                   T(flag.best && segconf.has[OPTION_we_i], OPTION_pw_B_C), // store sum of array elements
                   T(flag.best && segconf.has[OPTION_we_i], OPTION_ip_B_C), DID_EOL);

    ctx_set_same_line (VB, T(segconf.has[OPTION_ec_f], OPTION_np_i), DID_EOL);  // np segged as delta vs ec, and np needs to be peeked for QUAL, before reconstructing ec

    if (segconf.use_pacbio_iqsqdq) {
        ctx_set_ltype (VB, LT_BLOB, OPTION_dq_Z, OPTION_iq_Z, OPTION_sq_Z, DID_EOL);        
        ctx_consolidate_stats (VB, OPTION_iq_sq_dq, OPTION_dq_Z, OPTION_iq_Z, OPTION_sq_Z, DID_EOL);
    }

    if (segconf.sam_use_sn_mux)
        seg_mux_init (vb, OPTION_sn_B_f, SAM_SPECIAL_DEMUX_sn, false, sn);
}

// -------------
// zm:i
// -------------

void sam_seg_pacbio_zm (VBlockSAMP vb, int64_t zm, unsigned add_bytes)    
{                  
    // note: we can't delta vs Q1NAME in PRIM/DEPN bc QNAME is copied as a whole from SAG             
    if (IS_MAIN(vb) && zm == CTX(SAM_Q1NAME)->last_value.i) 
        seg_by_ctx (VB, STRa(copy_Q1NAME_int), CTX(OPTION_zm_i), add_bytes);

    else 
        seg_integer (VB, CTX(OPTION_zm_i), zm, true, add_bytes);
}

// -------------
// np:i
// -------------

void sam_seg_pacbio_np (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t np, unsigned add_bytes)    
{                     
    dl->np = np;

    if (segconf.pacbio_subreads && np == 1)
        seg_by_did (VB, "1", 1, OPTION_np_i, add_bytes); // in subreads, we expect np=1 all-the-same

    else if (segconf.has[OPTION_ec_f] && ctx_has_value_in_line_(vb, CTX(OPTION_ec_f))) 
        seg_delta_vs_other_localN (VB, CTX(OPTION_np_i), CTX(OPTION_ec_f), np, -1, add_bytes);

    else 
        seg_integer (VB, CTX(OPTION_np_i), np, true, add_bytes);
}

int32_t sam_zip_get_np (VBlockP vb, LineIType line_i)
{
    return DATA_LINE(line_i)->np;
}

// ----------------------------------------------
// qs:i - 0-based start of query in the ZMW read
// qe:i - 0-based end   of query in the ZMW read
// ----------------------------------------------

void sam_seg_pacbio_qs (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t qs, unsigned add_bytes)    
{
    decl_ctx (OPTION_qs_i);

    // in subreads: e.g. QNAME=m54284U_210913_013042/0/6137_11205 qe=11205 qs=6137  
    if (segconf.pacbio_subreads) {
        if (segconf_qf_id (QNAME1) == QF_PACBIO_rng && ctx_has_value_in_line_(VB, CTX(SAM_Q2NAME)) && qs == CTX(SAM_Q2NAME)->last_value.i)
            seg_by_ctx (VB, STRa(copy_Q2NAME_int), ctx, add_bytes);

        else
            seg_integer (VB, ctx, qs, true, add_bytes); 
    }

    // in ccs: qs has a small number of possible low values, often mono-value
    else {
        if (segconf.running) {
            if (!vb->line_i) 
                segconf.sam_first_qs = qs;
            else             
                segconf.sam_diverse_qs |= (qs != segconf.sam_first_qs);
        }

        if (segconf.sam_diverse_qs)
            seg_integer (VB, ctx, qs, true, add_bytes); 

        else 
            seg_integer_as_snip_do (VB, ctx, qs, add_bytes); // this is usually all-the-same
    }
}

void sam_seg_pacbio_qe (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t qe, unsigned add_bytes)    
{
    decl_ctx (OPTION_qe_i);

    // in subreads: e.g. QNAME=m54284U_210913_013042/0/6137_11205 qe=11205 qs=6137  
    if (segconf.pacbio_subreads) {
        if (segconf_qf_id (QNAME1) == QF_PACBIO_rng && ctx_has_value_in_line_(VB, CTX(SAM_Q3NAME)) && qe == CTX(SAM_Q3NAME)->last_value.i)
            seg_by_ctx (VB, STRa(copy_Q3NAME_int), ctx, add_bytes);

        else
            seg_integer (VB, ctx, qe, true, add_bytes); 
    }

    // in ccs: qe = qs + seq_len
    else {
        int32_t qs;
        if (sam_seg_peek_int_field (vb, OPTION_qs_i, vb->idx_qs_i, 0, 0x7ffffff, true, &qs) && qs + dl->SEQ.len == qe) 
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PACBIO_qe }, 2, ctx, add_bytes);
        
        else
            seg_integer (VB, ctx, qe, true, add_bytes); 
    }
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_PACBIO_qe)
{
    int64_t qs = reconstruct_peek_(vb, _OPTION_qs_i, 0, 0).i;
    new_value->i = qs + vb->seq_len;

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

// ----------------------------------------------
// ws:i - Start of first base of the query (‘qs’) in approximate raw frame count since start of movie. For a CCS read, the start of the first base of the first incorporated subread.
// we:i - Start of last base of the query (‘qe - 1’) in approximate raw frame count since start of movie. For a CCS read, the start of the last base of the last incorporated subread.
// ----------------------------------------------

void sam_seg_pacbio_we (VBlockSAMP vb, ZipDataLineSAMP dl, int64_t we, unsigned add_bytes)    
{
    decl_ctx (OPTION_we_i);

    int32_t ws;

    // prediction (subreads data): we ~= ws + (sum(ip) + sum (pw) + seq_len)
    if (flag.best && // summing ip/pw add time - not worth the tiny benefit of this method unless --best
        segconf.has[OPTION_we_i] && // condition for ip and pw sum to be calculated
        sam_seg_peek_int_field (vb, OPTION_ws_i, vb->idx_ws_i, 0, 0x7ffffff, true, &ws)  &&
        ctx_has_value_in_line_(vb, CTX(OPTION_ip_B_C)) && 
        ctx_has_value_in_line_(vb, CTX(OPTION_pw_B_C))) {
     
        int64_t sum_ip = CTX(OPTION_ip_B_C)->last_value.i;
        int64_t sum_pw = CTX(OPTION_pw_B_C)->last_value.i;
        int64_t delta = we - (ws + sum_ip + sum_pw + dl->SEQ.len);

        seg_integer (VB, ctx, delta, false, 0);
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_PACBIO_we }, 2, ctx, add_bytes);
    }
    
    else 
        seg_integer (VB, ctx, we, true, add_bytes); // delta vs ws is worse 
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_PACBIO_we)
{
    int64_t sum_ip = CTX(OPTION_ip_B_C)->last_value.i;
    int64_t sum_pw = CTX(OPTION_pw_B_C)->last_value.i;
    int64_t ws = reconstruct_peek (vb, CTX(OPTION_ws_i), 0, 0).i;
    int64_t delta = reconstruct_from_local_int (vb, ctx, 0, RECON_OFF);

    new_value->i = delta + (ws + sum_ip + sum_pw + vb->seq_len);

    if (reconstruct) RECONSTRUCT_INT (new_value->i);

    return HAS_NEW_VALUE;
}

//-------------------------------
// QUAL - predicted by iq, sq, dq
//-------------------------------

static inline char qual_prediction (char iq, char sq, char dq)
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
        *next++ = qual[i] - qual_prediction (iq[i], sq[i], dq[i]);
    
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
            *qual++ = (*diff++) + qual_prediction (il[0], il[1] - il[0], -il[2]);
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

void sam_seg_pacbio_xq (VBlockSAMP vb, ZipDataLineSAMP dl, Did did_i, TxtWord *dl_word, STRp(value), unsigned add_bytes)    
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
    ZipDataLineSAMP dl = DATA_LINE (vb_line_i);
    
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

// ---------------------------------------------------------------------------------------------------------
// sn:B:f 4 floats for the average signal-to-noise ratio of A, C, G, and T (in that order) over the HQRegion
// ---------------------------------------------------------------------------------------------------------

static int get_sn_channel_i (VBlockSAMP vb)
{
    return ctx_encountered_in_prev_line (VB, OPTION_sn_B_f) && CTX(SAM_Q1NAME)->last_delta == 0;
}

void sam_seg_pacbio_sn (VBlockSAMP vb, ZipDataLineSAMP dl, rom sn, int/*signed*/ sn_len)
{
    decl_ctx (OPTION_sn_B_f);
    ContextP channel_ctx = ctx; // fallback

    // sn:B is identical for all subreads of a molecule
    if (segconf.sam_use_sn_mux) {
        channel_ctx = seg_mux_get_channel_ctx (VB, OPTION_sn_B_f, (MultiplexerP)&vb->mux_sn, get_sn_channel_i (vb));

        STRlast (prev_sn, OPTION_sn_B_f);
        if (str_issame (sn, prev_sn)) {
            unsigned add_bytes =  (IS_BAM_ZIP ? (sn_len * 4/*width of type 'f'*/ + 4/*count*/ + 1/*type*/) 
                                              : (sn_len + 2/*type - eg "i,"*/ + 1/*\t or \n*/));

            seg_by_ctx (VB, STRa(copy_sn_snip), channel_ctx, add_bytes);
        }

        else 
            sam_seg_array_one_ctx (vb, dl, channel_ctx->dict_id, 'f', STRa(sn), 0, 0, NULL);

        seg_by_did (VB, STRa(vb->mux_sn.snip), OPTION_sn_B_f, 0); // de-multiplexer
    }

    else 
        sam_seg_array_one_ctx (vb, dl, _OPTION_sn_B_f, 'f', STRa(sn), 0, 0, NULL);

    seg_set_last_txt (VB, ctx, STRa(sn));
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_DEMUX_sn)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), get_sn_channel_i (VB_SAM), new_value, reconstruct);
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