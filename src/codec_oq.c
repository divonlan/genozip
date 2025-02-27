// ------------------------------------------------------------------
//   oq.c
//   Copyright (C) 2024-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// WARNING: THIS FILE CONTAINS A METHOD THAT IS PATENT PENDING.

#include "compressor.h"
#include "sam_private.h"

#define NUM_OQ_CTXS 94

#define OQ_DICT_ID_Q(q) DICT_ID_MAKE2_5(((char[]){'O', q+33, 'Q',':','Z'}))

#define decl_oq_ctxs_zip        \
    ContextP ctxs[NUM_OQ_CTXS]; \
    for (int q=0; q < NUM_OQ_CTXS; q++) ctxs[q] = ctx_get_ctx (vb, OQ_DICT_ID_Q(q));  

#define decl_oq_ctxs_piz        \
    ContextP ctxs[NUM_OQ_CTXS]; \
    for (int q=0; q < NUM_OQ_CTXS; q++) ctxs[q] = ECTX (OQ_DICT_ID_Q(q));  

//--------------
// ZIP side
//--------------

bool codec_oq_comp_init (VBlockP vb)
{
    // verify that OQ can be segged against QUAL (i.e. QUAL exists iff OQ exists)
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {
        ZipDataLineSAMP dl = DATA_LINE (line_i);
        if (dl->OQ && dl->SEQ.len != dl->QUAL.len) return false;
    }

    decl_ctx (OPTION_OQ_Z);
    
    ctx->ltype  = LT_CODEC;
    ctx->lcodec = CODEC_OQ;

    decl_oq_ctxs_zip;
    ctx_consolidate_statsA (vb, OPTION_OQ_Z, ctxs, NUM_OQ_CTXS);

    for (int q=0; q < NUM_OQ_CTXS; q++) { 
        ctxs[q]->local_dep = DEP_L1;  
        ctxs[q]->ltype = LT_SUPP; // data used by codec - not reconstructed directly
    }

    return true;
}

COMPRESS (codec_oq_compress)
{
    START_TIMER;

    decl_oq_ctxs_zip;

    // first pass - count qual scores to allocate memory and allocate txt_len
    int32_t count_q[NUM_OQ_CTXS] = {};
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        ZipDataLineSAMP dl = DATA_LINE (line_i);
        txtSTR (qual, dl->QUAL);

        for (uint32_t i=0; i < qual_len; i++) {
            #ifdef DEBUG // note: codecs that are destructive to QUAL data (DOMQ, NORMQ...) must set QUAL context to local_dep >= DEP_L1
            ASSERT (qual[i] >= 33 && qual[i] <= 126, "%s/%u: Invalid QUAL[%u]=%d", VB_NAME, line_i, i, qual[i]);
            #endif
            count_q[(int)qual[i] - 33]++; // BAM note: qual in txt_data has been already converted to SAM values in bam_rewrite_qual
        }
    }

    char *next[NUM_OQ_CTXS] = {};
    for (int q=0; q < NUM_OQ_CTXS; q++) {
        if (!count_q[q]) continue;

        buf_alloc_exact (vb, ctxs[q]->local, count_q[q], char, CTX_TAG_LOCAL); 
        next[q] = B1STc (ctxs[q]->local);

        // move stats allocation to the individual contexts
        // note: QUAL context might be already merged, or not yet. if not yet, zctx->txt_len will go negative, and become positive again after the merge.
        ctx_update_zctx_txt_len (vb, ctx,    -count_q[q]);
        ctx_update_zctx_txt_len (vb, ctxs[q], count_q[q]);
    }

    // second pass - seg OQ into channels - mux by QUAL
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        ZipDataLineSAMP dl = DATA_LINE (line_i);
        txtSTR (qual, dl->QUAL);
        rom oq = Btxt(dl->OQ);
        uint32_t oq_len = dl->SEQ.len;

        if (!oq_len) continue; // note: oq_len=0 can either mean no OQ field in this line, or empty OQ field exists (in which case SEQ=*, enforced by sam_seg_other_qual)

        for (uint32_t i=0; i < qual_len; i++) 
            *next[(uint8_t)qual[i] - 33]++ = oq[i]; // qual[i] is always in SAM terms (even in BAM), but the channels start from 0 (i.e. in BAM terms)
    }

    // for monochar channels, we just store the mono character in OQ, and cancel the separate context
    char monochars[NUM_OQ_CTXS] = {};
    for (int q=0; q < NUM_OQ_CTXS; q++) 
        if (count_q[q] && str_is_monochar (B1STc(ctxs[q]->local), count_q[q])) {
            monochars[q] = *B1STc(ctxs[q]->local);
            ctxs[q]->local.len32 = 0;
        }

    // compress the tiny monochar array
    header->sub_codec = CODEC_RANB; // note: if changing this codec, also update the est_size function in CODEC_ARGS.OQ
    CodecCompress *compress = codec_args[header->sub_codec].compress;
    *uncompressed_len = NUM_OQ_CTXS;

    COPY_TIMER_COMPRESS (compressor_oq); 
    return compress (vb, ctx, header, monochars, uncompressed_len, NULL, STRa(compressed), false, name);

    return true;
}

//--------------
// PIZ side
//--------------

CODEC_RECONSTRUCT (codec_oq_reconstruct)
{
    START_TIMER;
    decl_oq_ctxs_piz;
    char sam_diff = (OUT_DT(BAM) || OUT_DT(CRAM)) ? 0 : 33;

    STRlast (qual, IS_PRIM(vb) ? SAM_QUALSA : SAM_QUAL);
    ASSPIZ (len == qual_len, "expecting len=%u == qual_len=%u", len, qual_len);

    rom monochars = B1STc(ctx->local); // 93 values, non-zero if channel_i is monochar

    char *recon = BAFTtxt;
    rom next[NUM_OQ_CTXS] = {}, after[NUM_OQ_CTXS] = {};

    for (int q=0; q < NUM_OQ_CTXS; q++) 
        if (ctxs[q]) {
            next[q]  = Bc(ctxs[q]->local, ctxs[q]->next_local);
            after[q] = BAFTc(ctxs[q]->local);
        }

    for (int32_t i=0; i < len; i++) {
        int channel_i = qual[i] - sam_diff;

        if (monochars[channel_i]) 
            recon[i] = monochars[channel_i];

        else {
            ASSPIZ (next[channel_i] < after[channel_i], "channel='%c' is out of data: len=%u",
                    channel_i + 33, ctxs[channel_i] ? ctxs[channel_i]->local.len32 : 0);

            recon[i] = *next[channel_i]++;
        }
    }

    if (reconstruct) Ltxt += len;

    for (int q=0; q < NUM_OQ_CTXS; q++) 
        if (ctxs[q]) 
            ctxs[q]->next_local = BNUM (ctxs[q]->local, next[q]);

    COPY_TIMER(codec_oq_reconstruct);
}
