// ------------------------------------------------------------------
//   codec_smux.c
//   Copyright (C) 2024-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// WARNING: THIS FILE CONTAINS A METHOD THAT IS PATENT PENDING.

#include <math.h>
#include "reconstruct.h"
#include "compressor.h"

// note: dict_id compatible with codec_pacb_smux_is_qual()
#define decl_smux_ctxs_zip(id)                                                                                                       \
    ContextP ctxs[5] = { ctx_get_ctx (vb, DICT_ID_MAKEF_8(((char[]){'A','A','A','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),  \
                         ctx_get_ctx (vb, DICT_ID_MAKEF_8(((char[]){'C','C','C','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),  \
                         ctx_get_ctx (vb, DICT_ID_MAKEF_8(((char[]){'G','G','G','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),  \
                         ctx_get_ctx (vb, DICT_ID_MAKEF_8(((char[]){'T','T','T','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),  \
                         ctx_get_ctx (vb, DICT_ID_MAKEF_8(((char[]){'N','N','N','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))) }; // N and other IUPACs

#define decl_smux_ctxs_piz(id)                                                                                                       \
    ContextP ctxs[5] = { ECTX(DICT_ID_MAKEF_8(((char[]){'A','A','A','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),              \
                         ECTX(DICT_ID_MAKEF_8(((char[]){'C','C','C','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),              \
                         ECTX(DICT_ID_MAKEF_8(((char[]){'G','G','G','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),              \
                         ECTX(DICT_ID_MAKEF_8(((char[]){'T','T','T','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))),              \
                         ECTX(DICT_ID_MAKEF_8(((char[]){'N','N','N','-',(id[0] & 0x7f) | 0x40, id[1], id[2], id[3]}))) }; // NULL if not data for this base


//--------------
// ZIP side
//--------------

static void normalize_histo (double histo[94])
{
    double total = 0;
    for (int i=0; i < 93; i++) total += histo[i];
    
    if (total)
        for (int i=0; i < 93; i++) histo[i] /= total;
} 

static void set_stdv (double *histo, int histo_len)
{
    #define HISTO(x,q) histo[(x)*94 + (q)]

    double max_var = -1;
    int8_t max_var_q = -1;
    for (int q=0; q < 94; q++) {
        double av_histo = 0;
        for (int x=0; x < histo_len; x++) 
            av_histo += HISTO(x,q) / (double)histo_len;

        double var = 0;
        for (int x=0; x < histo_len; x++) 
            var += SQR (HISTO(x,q) - av_histo) / (double)histo_len;

        if (var > max_var) {
            max_var = var;
            max_var_q = q;
        }
    }

    segconf.smux_max_stdv = sqrt (max_var);
    segconf.smux_max_stdv_q = max_var_q + '!';

    if (flag.show_qual)      
        iprintf ("SMUX threadshold: '%c' has the highest standard deviation: %2.1f%%\n", 
                 segconf.smux_max_stdv_q, 100.0 * segconf.smux_max_stdv);
}

// calculate stats for smux AND qmux as well as stats
void codec_smux_calc_stats (VBlockP vb)
{
    decl_ctx (SAM_QUAL); // == FASTQ_QUAL
    LocalGetLineCB *seq_callback  = (VB_DT(FASTQ) ? fastq_zip_seq  : sam_zip_seq);
    LocalGetLineCB *qual_callback = (VB_DT(FASTQ) ? fastq_zip_qual : sam_zip_qual);

    #define BP_IN_SAMPLE 30000

    double histo[94]              = {}; // independent histogram
    double histo_by_seq[5][94]    = {}; // histogram, dependent on SEQ
    
    double count_base[5]= {}; // in FASTQ terms - i.e. un-revcomped SAM sequences if needed
    int32_t num_bp = 0;       // only for lines included in QUAL

    // calculate histograms
    for (LineIType line_i=0; line_i < vb->lines.len32 && num_bp < BP_IN_SAMPLE; line_i++) {   
        bytes qual; 
        uint32_t qual_len;
        qual_callback (vb, ctx, line_i, (char **)pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        
        if (!qual_len) continue; // case in SAM: dl->dont_compress_QUAL

        STRw(seq);
        bool is_rev;
        seq_callback (vb, NULL, line_i, pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, &is_rev);

        if (qual[0] != ' ') { // ignore missing qual
            if (!is_rev)
                for (unsigned i=0; i < qual_len; i++) { // note: qual_len==1 if missing qual
                    count_base[nuke_encode(seq[i])]++;
                    histo[qual[i]-'!']++;
                    histo_by_seq[nuke_encode(seq[i])][qual[i]-33]++;
                }
            else
                for (unsigned i=0; i < qual_len; i++) {
                    count_base[nuke_encode_comp(seq[i])]++;
                    histo[qual[i]-'!']++;
                    histo_by_seq[nuke_encode_comp(seq[i])][qual[i]-33]++;
                }
        }
        
        num_bp += qual_len;
    }

    // normalize each histogram so that Σ(q=0...93)histo = 1 
    normalize_histo (histo);
    for (int b=0; b < 5; b++) normalize_histo (histo_by_seq[b]);

    #define SHOW_THREASHOLD 0.03

    if (flag.show_qual) { 
        decl_acgt_decode;
        iprintf ("%s: For each line, showing only qual scores which have > %2.0f%% prevalence in that line, and only lines with size>%2.0f%%:\n", 
                 VB_NAME, SHOW_THREASHOLD * 100.0, SHOW_THREASHOLD * 100.0);
        
        iprint0 ("Unmuxed:\nsize=100%: ");
        for (int q=0; q < 94; q++) 
            if (flag.show_qual && histo[q] > SHOW_THREASHOLD) 
                iprintf ("'%c'=%2.1f%% ", q + '!', 100.0 * histo[q]);
                
        iprint0 ("\nMux by SEQ:\n");
        for (int b=0; b < 4; b++) {
            iprintf ("%c: size=%2.1f%%: ", acgt_decode(b), 100.0 * count_base[b] / num_bp);
            for (int q=0; q < 94; q++) {
                if (histo_by_seq[b][q] > SHOW_THREASHOLD) 
                    iprintf ("'%c'=%2.1f%% ", q + '!', 100.0 * histo_by_seq[b][q]);
            }
            iprint_newline();
        }

        iprint_newline();
    }

    set_stdv ((double *)histo_by_seq, 4);
}

// ZIP: called for QUAL-like dids
bool codec_smux_maybe_used (Did did_i)
{
    #define SMUX_STDV_THREADHOLD 0.08 // between 0 and 1. The higher this is, the more drastic differences between the quality histogram associated with each base is needed to qualify for SMUX

    return did_i == SAM_QUAL/*==FASTQ_QUAL*/ && // we only calculated stats for SAM_QUAL (not OQ)
               (flag.force_qual_codec == CODEC_SMUX || 
                (TECH(MGI) && segconf.nontrivial_qual && !flag.no_smux && segconf.smux_max_stdv > SMUX_STDV_THREADHOLD)); // not yet seen benefit for non-MGI files);
}

bool codec_smux_comp_init (VBlockP vb, Did qual_did_i, LocalGetLineCB get_line_cb, bool force)
{
    if (!force && !codec_smux_maybe_used (qual_did_i)) return false;

    decl_ctx (qual_did_i);
    
    ctx->ltype     = LT_CODEC;
    ctx->lcodec    = CODEC_SMUX;
    ctx->local_dep = DEP_L1; // yield to other codecs (eg CODEC_OQ) that need to query QUAL before we destroy it

    decl_smux_ctxs_zip (ctx->dict_id.id);
    ctx_consolidate_statsA (vb, qual_did_i, ctxs, 5);

    for (int b=0; b < 5; b++) { 
        ctxs[b]->local_dep = DEP_L2;  
        ctxs[b]->ltype = LT_SUPP; // data used by codec - not reconstructed directly
    }

    return true;
}

COMPRESS (codec_smux_compress)
{
    START_TIMER;

    LocalGetLineCB *seq_callback = (VB_DT(FASTQ) ? fastq_zip_seq : sam_zip_seq);

    decl_smux_ctxs_zip (ctx->dict_id.id);

    int32_t count_base[5] = {};
    STRw(seq); STRw(qual);
    char *next[5];
    bool is_rev;

    // first pass - count bases to allocate memory and allocate txt_len
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        get_line_cb (vb, ctx, line_i, pSTRa(qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        if (!qual_len) continue;

        seq_callback (vb, NULL, line_i, pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, &is_rev);
        
        if (!is_rev)
            for (uint32_t i=0; i < qual_len; i++) count_base[nuke_encode(seq[i])]++;

        else if (IS_SPACE(qual))
            count_base[nuke_encode_comp(seq[seq_len-1])]++; 

        else
            for (uint32_t i=0; i < qual_len; i++) count_base[nuke_encode_comp(seq[i])]++;
    }

    for (int b=0; b < 5; b++) {
        buf_alloc_exact (vb, ctxs[b]->local, count_base[b], char, CTX_TAG_LOCAL); 
        next[b] = B1STc(ctxs[b]->local);

        // move stats allocation to the individual contexts
        // note: QUAL context might be already merged, or not yet. if not yet, zctx->txt_len will go negative, and become positive again after the merge.
        ctx_update_zctx_txt_len (vb, ctx, -count_base[b]);
        ctx_update_zctx_txt_len (vb, ctxs[b], count_base[b]);
    }

    // second pass - move each qual base to its correct A,C,G,T or N context.
    for (LineIType line_i=0; line_i < vb->lines.len32; line_i++) {   
        
        get_line_cb (vb, ctx, line_i, pSTRa (qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        if (!qual_len) continue;

        STRw(seq);
        bool is_rev;
        seq_callback (vb, NULL, line_i,  pSTRa(seq), CALLBACK_NO_SIZE_LIMIT, &is_rev);

        // note: if missing qual, qual_len==1 and qual==' ' and we seg to the channel of the first base. 
        // if also missing seq, then seq=='*' and we seg to channel 4
        if (!is_rev)
            for (int32_t i=0; i < qual_len; i++) 
                // note: if qual==" " (missing) it goes into channel by seq[0], and loop ends here bc qual_len==1
                *next[nuke_encode(seq[i])]++ = qual[i]; 
        
        // note: if reversed qual==" " (missing), it is put into channel by seq[seq_len-1] and not seq[0] - 
        // this is because recon won't know there's a missing quality until it reconstructs
        else if (IS_SPACE(qual))
            *next[nuke_encode_comp(seq[seq_len-1])]++ = ' '; 

        // case: reversed seq, qual NOT missing
        else 
            for (int32_t i=qual_len-1; i >= 0; i--) 
                *next[nuke_encode_comp(seq[i])]++ = qual[i]; 
    }

    for (int b=0; b < 5; b++) 
        ASSERT (next[b] == BAFTc(ctxs[b]->local), "%s: error writing into %s.local (b=%d): expecting bytes=%u == len=%u", 
                VB_NAME, ctxs[b]->tag_name, b, BNUM(ctxs[b]->local, next[b]), ctxs[b]->local.len32);

    // usually all non-ACGT bases are 'N' and they all of the same score - we can drop the section
    if (ctxs[4]->local.len && str_is_monochar (STRb(ctxs[4]->local))) {
        ((SectionHeaderCtxP)header)->param = *B1ST8(ctxs[4]->local);

        ctx_update_zctx_txt_len (vb, ctx, ctxs[4]->local.len); // give it back to ctx, as ctxs[4] won't exist
        ctx_update_zctx_txt_len (vb, ctxs[4], -(int64_t)ctxs[4]->local.len);
        ctxs[4]->local.len = 0;
    }

    // write one byte to the SEC_LOCAL section, just to trigger reconstruction
    header->sub_codec = CODEC_NONE;
    *compressed = 0;
    *uncompressed_len = *compressed_len = 1;
    
    COPY_TIMER_COMPRESS (compressor_smux); 
    return true;
}

//--------------
// PIZ side
//--------------

CODEC_RECONSTRUCT (codec_smux_reconstruct)
{
    START_TIMER;
    decl_smux_ctxs_piz (ctx->dict_id.id);

    rom seq = VB_DT(SAM) ? sam_piz_get_textual_seq(vb) : last_txtx (vb, CTX(FASTQ_SQBITMAP)); 
    if (*seq == '*') 
        len = 1;
    
    // case: Deep, and len is only the trimmed suffix as the rest if copied from SAM (see fastq_special_deep_copy_QUAL)
    else if (flag.deep && len < vb->seq_len)
        seq += (vb->seq_len - len); // advance seq to the trimmed part too
    else
        ASSPIZ (len == vb->seq_len, "expecting len=%u == vb->seq_len=%u", len, vb->seq_len);

    bool is_src_rev   = VB_DT(SAM) ? sam_is_last_flags_rev_comp (vb) : false;
    bool is_recon_rev = OUT_DT(FASTQ) ? false : is_src_rev;

    char *recon = BAFTtxt;
    rom next_b[5], after_b[5];

    for (int b=0; b < 5; b++) {
        next_b[b]  = ctxs[b] ? Bc(ctxs[b]->local, ctxs[b]->next_local) : NULL;
        after_b[b] = ctxs[b] ? BAFTc(ctxs[b]->local) : NULL;
    }

    if (!is_recon_rev) {
        for (int32_t i=0; i < len; i++) {
            int channel_i = nuke_encode(seq[i]);

            if (channel_i == 4 && ctx->local.param) // case: N quality score is monochar and sent through param 
                recon[i] = ctx->local.param;

            else {
                ASSPIZ (next_b[channel_i] < after_b[channel_i], "%s is out of data: len=%u (channel_i=%u)",
                        ctxs[channel_i] ? ctxs[channel_i]->tag_name    : "", 
                        ctxs[channel_i] ? ctxs[channel_i]->local.len32 : 0, 
                        channel_i);

                recon[i] = *next_b[channel_i]++;
            }

            if (recon[i] == ' ') {
                sam_reconstruct_missing_quality (vb, reconstruct);
                goto done;
            }
        }

        // case: translating to FASTQ, and alignment is revcomped in SAM - reverse it
        if (is_src_rev)
            str_reverse_in_place (recon, len);
    }
    else
        for (int32_t i=len-1; i >= 0; i--) {
            int channel_i = nuke_encode_comp(seq[i]);

            if (channel_i == 4 && ctx->local.param) // case: N quality score is monochar and sent through param 
                recon[i] = ctx->local.param;

            else {
                ASSPIZ (!ctxs[channel_i] || next_b[channel_i] < after_b[channel_i], "%s is out of data: len=%u",
                        ctxs[channel_i]->tag_name, ctxs[channel_i]->local.len32);

                ASSPIZ (next_b[channel_i], "next_b[%u] is NULL", channel_i);

                recon[i] = *next_b[channel_i]++;
            }

            if (recon[i] == ' ') {
                sam_reconstruct_missing_quality (vb, reconstruct);
                goto done;
            }
        }

    if (reconstruct) Ltxt += len;

done:
    for (int b=0; b < 5; b++)
        if (ctxs[b]) 
            ctxs[b]->next_local = BNUM (ctxs[b]->local, next_b[b]);

    COPY_TIMER(codec_smux_reconstruct);
}
