// ------------------------------------------------------------------
//   codec_tmpl.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// WARNING: THIS FILE CONTAINS A METHOD THAT IS PATENT PENDING.

#include <math.h>
#include "reconstruct.h"
#include "compressor.h"

//--------------
// ZIP side
//--------------

static int count_q_in_tmpl[94]; // how many times each value q appears in the template

// need to do for for each component as template may differ
void codec_tmpl_segconf_finalize (VBlockP vb, Did did_i, LocalGetLineCB get_line_cb)
{
    ASSERTNOTINUSE (vb->scratch);
    decl_ctx (did_i);
    decl_zctx (did_i); // note: ZCTX(did_i) permitted, because did_i is predefined
    
    if (flag.force_qual_codec && flag.force_qual_codec != CODEC_TMPL)
        return;

    uint32_t template_len = IS_R2 ? segconf.std_seq_lR2 : segconf.std_seq_len;

    buf_alloc_exact_zero (vb, vb->scratch, template_len * 94, uint32_t, "scratch");
    uint32_t (*histogram)[94] = (void*)B1STc(vb->scratch);
    
    for_line {
        STRw(qual);
        get_line_cb (vb, ctx, line_i, pSTRa(qual), CALLBACK_NO_SIZE_LIMIT, NULL);

        MINIMIZE(qual_len, template_len); // note: template_len is the longest seq_len in the segconf data

        for (uint32_t i=0; i < qual_len; i++)
            histogram[i][(int)qual[i] - 33]++;
    }

    int n_have_majority_score = 0;

    buf_alloc_exact (evb, (zctx+1)->template, template_len, char, "template"); // note: possibly override previous component's template
    char *template = B1STc((zctx+1)->template);
    
    memset (count_q_in_tmpl, 0, sizeof(count_q_in_tmpl));

    for (uint32_t i=0; i < template_len; i++) {
        char max_q=0;
        uint32_t max_q_count=0;

        for (int q=0; q < 94; q++)
            if (histogram[i][q] > max_q_count) {
                max_q_count = histogram[i][q];
                max_q = q + 33;
            }
        
        template[i] = max_q;
        count_q_in_tmpl[max_q - 33]++;

        if (percent (max_q_count, vb->lines.len32) >= 40)
            n_have_majority_score++;
    }

    buf_free (vb->scratch);

    // We can use TMPL if over 30% of positions have a dominant score (=a score appearing in over 40% of lines)
    if (percent (n_have_majority_score, template_len) < 30) {
        // fallback to BSC
        zctx->lcodec = zctx->qual_codec = (flag.fast ? CODEC_ARTB : CODEC_BSC);
        zctx->lcodec_hard_coded = true;
        zctx->tmpl_calculated   = false;
        return;
    }

    // create subdicts array (transmitted to PIZ)
    DictId tmpl_dict_id = (DictId)DICT_ID_MAKEF_4("tmpl"); // note: don't use ctx->dict_id==QUAL because subdicts Q?UAL will conflict with Q?NAME subdicts in d2d_map

    // subdicts contains multiple concatenated arrays of DictId, one for each FASTQ component
    buf_alloc_zero (evb, &zctx->subdicts, template_len, 0, DictId, 0, "subdicts");
    DictId *subdicts = BAFT(DictId, zctx->subdicts);
    zctx->subdicts.len += template_len;

    // populate zctx->subdicts and create a zctx for each subdict
    for (uint32_t i=0; i < template_len; i++) 
        if (!subdicts[i].num) {
            subdicts[i] = sub_dict_id (tmpl_dict_id, template[i]);
            
            // create zctx
            ContextP vctx = ctx_get_ctx (vb, subdicts[i]);
            ctx_get_zctx_from_vctx (vctx, true, false);
        }
    
    zctx->tmpl_calculated  = true;
    zctx->subdicts_section = true; // we store the dict_id array in SEC_SUBDICTS (global section)

    if (flag.show_qual)
        iprintf ("TMPL: %s template: %.*s\n", ctx->tag_name, STRf(template));
}

// ZIP: called for QUAL-like dids
bool codec_tmpl_maybe_used (Did did_i)
{
    return TXT_DT(FASTQ) && // TO DO: support for SAM, requires handling reversing and trimming
           (did_i == SAM_QUAL/*==FASTQ_QUAL*/ || did_i == OPTION_OQ_Z) &&
           !flag.deep &&    // Can't use in Deep, because std_seq_len relates to SAM 
           (segconf.running || ZCTX(did_i)->tmpl_calculated) && 
               (flag.force_qual_codec == CODEC_TMPL || 
                (TECH(ELEMENT) && segconf.nontrivial_qual && !flag.no_tmpl)); 
}

bool codec_tmpl_comp_init (VBlockP vb, Did qual_did_i, bool force)
{
    if (!force && !codec_tmpl_maybe_used (qual_did_i)) return false;

    decl_ctx (qual_did_i);
    decl_zctx (qual_did_i);
    
    ctx->ltype     = LT_CODEC;
    ctx->lcodec    = CODEC_TMPL;
    ctx->local_dep = DEP_L1; 

    for_buf (DictId, dict_id_p, zctx->subdicts) { // a bit redundant, as each value may appear many times, but no harm
        ContextP subctx   = ctx_get_ctx (vb, *dict_id_p);
        subctx->ltype     = LT_SUPP; // data used by codec - not reconstructed directly
        subctx->local_dep = DEP_L2;
        subctx->lcodec    = CODEC_BSC;
        subctx->lcodec_hard_coded = true;

        ctx_consolidate_stats (vb, qual_did_i, subctx->did_i, DID_EOL);
    }

    // excess values beyond template_len go into (ctx+1)->local
    (ctx+1)->ltype     = LT_SUPP;
    (ctx+1)->local_dep = DEP_L2;
        
    return true;
}

COMPRESS (codec_tmpl_compress)
{
    START_TIMER;
    decl_zctx (ctx->did_i);

    rom template = B1STc((zctx+1)->template);
    uint32_t template_len = (zctx+1)->template.len32;

    ASSERTNOTNULL(template);

    char *chan_p[94] = {};
    ContextP subctxs[94] = {};

    // allocate channel buffers (maximum - assuming all reads are template_len)
    const DictId *subdicts = B(DictId, zctx->subdicts, zctx->subdicts.len32 - template_len);

    for (uint32_t i=0; i < template_len; i++) {
        int q = template[i] - 33;
        if (!subctxs[q]) {
            subctxs[q] = ctx_get_ctx (vb, subdicts[i]);
            buf_alloc (vb, &subctxs[q]->local, 0, count_q_in_tmpl[q] * vb->lines.len, char, 0, C_LOCAL);

            chan_p[q] = B1STc(subctxs[q]->local);
        }
    }

    // copy qual scores, each to its channel determined by the template
    for_line {
        STRw(qual);
        get_line_cb (vb, ctx, line_i, pSTRa(qual), CALLBACK_NO_SIZE_LIMIT, NULL);
        
        uint32_t limited_qual_len = MIN_(qual_len, template_len);
        for (uint32_t i=0; i < limited_qual_len; i++) {
            int q = template[i] - 33; // channel
            *chan_p[q]++ = qual[i];
        }

        // excess qual data beyond template_len to ctx+1
        if (qual_len > limited_qual_len)
            buf_add_more (vb, &(ctx+1)->local, &qual[limited_qual_len], qual_len - limited_qual_len, C_LOCAL);
    }

    for (int q=0; q < 94; q++)
        if (count_q_in_tmpl[q]) {
            ASSERT (IN_RANGX(chan_p[q], subctxs[q]->local.data, subctxs[q]->local.data + subctxs[q]->local.size),
                    "chan_p[%u='%c'] out of range", q, q+33); // sanity

            // update lengths
            subctxs[q]->local.len32 = BNUM (subctxs[q]->local, chan_p[q]);

            // update accounting (doing it in Z, because already merged)
            ContextP sub_zctx = ctx_get_zctx_from_vctx (subctxs[q], false, false);
            ASSERTNOTNULL (sub_zctx);

            add_relaxed (sub_zctx->txt_len, subctxs[q]->local.len32);
            sub_relaxed (zctx->txt_len,     subctxs[q]->local.len32);
        }

    // update accounting for excess
    add_relaxed ((zctx+1)->txt_len, (ctx+1)->local.len32);
    sub_relaxed (zctx->txt_len,     (ctx+1)->local.len32); // all transferred to subcontexts already
    
    // write one byte to the SEC_LOCAL section, just to trigger reconstruction
    header->sub_codec = CODEC_NONE;
    *compressed = 0;
    *uncompressed_len = *compressed_len = 1;

    COPY_TIMER_COMPRESS_BY_CODEC (compressor_tmpl); 
    return true;
}

//--------------
// PIZ side
//--------------

CODEC_RECONSTRUCT (codec_tmpl_reconstruct)
{
    START_TIMER;
    decl_zctx (ctx->did_i);

    uint32_t template_len = (vb->comp_i == FQ_COMP_R2 && flag.pair) ? segconf.std_seq_lR2 : segconf.std_seq_len;

    if (!ctx->is_initialized) {
        buf_alloc_exact (vb, ctx->channel_data, template_len, ContextP, "channel_data"); 
        
        const DictId *subdicts = B(DictId, zctx->subdicts, vb->comp_i * segconf.std_seq_len); // note: whether or not we are R2, all previous component's template length is std_seq_len

        for_buf (ContextP, subctx, ctx->channel_data)
            *subctx = ECTX(*subdicts++);

        ctx->is_initialized = true;
    }

    ARRAY (ContextP, subctx, ctx->channel_data);

    char *next_recon = BAFTtxt;

    uint32_t limited_qual_len = MIN_(len, zctx->subdicts.len32);

    if (reconstruct) {
        for (uint32_t i=0; i < limited_qual_len; i++) 
            *next_recon++ = *Bc(subctx[i]->local, subctx[i]->next_local++);

        Ltxt += limited_qual_len;
    }

    else // consume without reconstructing
        for (uint32_t i=0; i < limited_qual_len; i++) 
            subctx[i]->next_local++;

    if (limited_qual_len < len) 
        RECONSTRUCT_NEXT (ctx+1, len - limited_qual_len); // works whether reconstruct or not

    COPY_TIMER(codec_tmpl_reconstruct);
}
