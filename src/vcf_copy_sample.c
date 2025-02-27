// ------------------------------------------------------------------
//   vcf_copy.c
//   Copyright (C) 2024-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include "vcf_private.h"
#include "zip_dyn_int.h"
#include "lookback.h"

// TxtWord of last sample copied for this sample_i and FORMAT
#define LAST_SAMPLE_SAME_FMT_ZIP *B(TxtWord, CTX(VCF_SAMPLES)->last_samples, dl->format_node_i * vcf_num_samples + vb->sample_i)
#define LAST_SAMPLE_SAME_FMT_PIZ *B(TxtWord, CTX(VCF_SAMPLES)->last_samples, CTX(VCF_FORMAT)->last_value.i * vcf_num_samples + VB_VCF->sample_i)

// true if previous sample of this sample_i and FORMAT was copied
#define SAMPLE_COPIED_SAME_FMT_ZIP *B(bool, CTX(VCF_COPY_SAMPLE)->sample_copied, dl->format_node_i * vcf_num_samples + vb->sample_i)
#define SAMPLE_COPIED_SAME_FMT_PIZ *B(bool, CTX(VCF_COPY_SAMPLE)->sample_copied, CTX(VCF_FORMAT)->last_value.i * vcf_num_samples + VB_VCF->sample_i)

//------------
// ZIP
//------------

void vcf_copy_sample_seg_initialize (VBlockVCFP vb)
{
    decl_ctx (VCF_COPY_SAMPLE);

    if (segconf_running) {
        if ((segconf.vcf_is_gvcf && vcf_num_samples >= 1) || vcf_num_samples >= 5) // gvcf or any file with 5 or more samples
            segconf.vcf_sample_copy = true; // initialize optimistically
    }

    if (segconf.vcf_sample_copy) {
        uint32_t n_fmts = CTX(VCF_FORMAT)->ol_nodes.len;
        buf_alloc_exact_zero (vb, CTX(VCF_SAMPLES)->last_samples,  n_fmts * vcf_num_samples, TxtWord, "contexts->last_samples");
        buf_alloc_exact_zero (vb, CTX(VCF_COPY_SAMPLE)->sample_copied, n_fmts * vcf_num_samples, bool, "contexts->sample_copied");
   
        seg_mux_init (vb, FORMAT_GT, VCF_SPECIAL_MUX_BY_PREV_COPIED, false, GT);

        seg_init_all_the_same (VB, VCF_COPY_SAMPLE, (char[]){ SNIP_SPECIAL, VCF_SPECIAL_COPY_SAMPLE }, 2);
        
        ctx_set_dyn_int (VB, VCF_COPY_SAMPLE, DID_EOL);
        ctx->dyn_transposed = true;
        ctx->local_dep = DEP_L1; // transpose this matrix itself only after using it to transpose all other fields 

        buf_alloc (vb, &ctx->local, 0, vb->lines.len * vcf_num_samples, bool, 0, CTX_TAG_LOCAL); // initial allocation
    }
}

void vcf_copy_samples_segconf_finalize (VBlockVCFP vb)
{
    // dis-vcf_sample_copy if we have evidence that too few samples were copied 
    // note: risk that copied percent is non-characteristically high in segconf bc variants are*near the telomere
    double pc_samples_copied = (segconf.vcf_sample_copy && vb->lines.len > 1) 
        ? ((double)CTX(VCF_COPY_SAMPLE)->num_samples_copied / (double)((vb->lines.len-1/*line=0 no copies*/) * vcf_num_samples)) : 0;
    
    if (vb->lines.len > 1 && pc_samples_copied < 0.50) // note: if very-high-sample-count files, segconf will be only 1 line, in which case we will default to vcf_sample_copy
        segconf.vcf_sample_copy = false;
}

void vcf_copy_sample_seg_finalize (VBlockVCFP vb)
{
    decl_ctx (VCF_COPY_SAMPLE);

    // remove VCF_COPY_SAMPLE if no samples were copied
    if (*B1ST8(ctx->local) == false && str_is_monochar (STRb(ctx->local))) 
        ctx->local.len = 0;
}

unsigned vcf_seg_copy_one_sample (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP *ctxs, ContainerP format, STRp(sample))
{
    START_TIMER;

    decl_ctx (VCF_COPY_SAMPLE);
    bool success = false;

    seg_all_the_same (VB, ctx, 0);

    // case: sample is identical to previous line - just copy
    if (str_issame_(STRa(sample), STRtxt (LAST_SAMPLE_SAME_FMT_ZIP))) {        
        success = true;

        // update txt_len for all subfields
        str_split (sample, sample_len, con_nitems (*format) - 1, ':', sf, false);

        for (uint32_t i=0; i < n_sfs; i++) {
            // account for subfields
            ctxs[i]->txt_len += sf_lens[i];

            set_last_txt_(ctxs[i]->did_i, sfs[i], sf_lens[i]);

            // save values needed for future FORMAT fields, or deferred INFO fields
            switch (ctxs[i]->did_i) {
                case FORMAT_GT:
                    vcf_seg_analyze_copied_GT (vb, STRi(sf,i));
                    break;

                case FORMAT_DP: { // for vcf_seg_INFO_DP_by_FORMAT_DP
                    int64_t dp;
                    if (!str_is_1chari(sf,i,'.') && str_get_int (STRi(sf,i), &dp)) {
                        if (segconf.INFO_DP_method == BY_FORMAT_DP) 
                            CTX(INFO_DP)->dp.sum_format_dp += dp;

                        vcf_seg_sum_DP_for_QD (vb, dp);
                    }
                    break;
                }
                case FORMAT_SB: // for vcf_seg_INFO_AS_SB_TABLE
                    vcf_sum_SB_for_AS_SB_TABLE (vb, STRi(sf,i));
                    break;

                case FORMAT_PS:
                case FORMAT_PID:
                    lookback_insert (VB, VCF_LOOKBACK, ctxs[i]->did_i, false, TXTWORDi(sf,i));
                    break;

                // note: vcf_seg_analyze_GT takes care of analyzing GT for vcf_seg_INFO_SF_seg
                default : {}
            }
        }

        ctx->num_samples_copied++;
        dyn_int_append (VB, ctx, true, 0);
    }

    // case: not the same as in previous line - return false, to seg normally
    else  
        dyn_int_append (VB, ctx, false, 0);

    if (!success) // no need to change the word if successful - it is identical    
        LAST_SAMPLE_SAME_FMT_ZIP = TXTWORD(sample);

    COPY_TIMER (vcf_seg_copy_one_sample);
    return success;
}

void vcf_copy_sample_seg_set_copied (VBlockVCFP vb, ZipDataLineVCFP dl, bool is_copied)
{
    SAMPLE_COPIED_SAME_FMT_ZIP = is_copied;
}

//------------
// PIZ
//------------

void vcf_sample_copy_piz_init_vb (VBlockVCFP vb)
{
    // hoist VCF_COPY_SAMPLE.local as it needs to be prepared (untranposed etc) before other transposed sections (AaD, DP...) are untransposed
    for_buf (uint32_t, header_offset, vb->z_section_headers) {
        SectionHeaderCtxP ctx_header = (SectionHeaderCtxP)Bc(vb->z_data, *header_offset);
        if (ctx_header->section_type == SEC_LOCAL && ctx_header->dict_id.num == _VCF_COPY_SAMPLE) {
            SWAP (*header_offset, *B1ST32(vb->z_section_headers));
            break;
        }
    }

    buf_alloc_exact_zero (vb, CTX(VCF_SAMPLES)->last_samples,      vcf_num_samples * ZCTX(VCF_FORMAT)->word_list.len, TxtWord, "contexts->last_samples");
    buf_alloc_exact_zero (vb, CTX(VCF_COPY_SAMPLE)->sample_copied, vcf_num_samples * ZCTX(VCF_FORMAT)->word_list.len, bool, "contexts->sample_copied");
}

void vcf_copy_sample_piz_store (VBlockVCFP vb, STRp(recon_sample))
{
    LAST_SAMPLE_SAME_FMT_PIZ   = TXTWORD(recon_sample);
    SAMPLE_COPIED_SAME_FMT_PIZ = CTX(VCF_COPY_SAMPLE)->last_value.i;
}

SPECIAL_RECONSTRUCTOR_DT (vcf_piz_special_COPY_SAMPLE)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    new_value->i = ctx->local.len ? reconstruct_from_local_int (VB, ctx, 0, false) 
                                  : 0; // local was dropped because no sample is copied (see vcf_copy_sample_seg_finalize) 

    rom recon = BAFTtxt;

    if (new_value->i) {
        TxtWord tw = LAST_SAMPLE_SAME_FMT_PIZ;
        RECONSTRUCT (Btxt(tw.index), tw.len);
        if (*BLSTtxt == '\t') Ltxt--; // \t will be reconstructed by the container

        ConstContainerP con = vb->con_stack[vb->con_stack_len-1].con;

        str_split (recon, BAFTtxt - recon, con_nitems(*con)-1, ':', sf, false);

        for (int i=0; i  < n_sfs; i++) {
            ContextP item_ctx = !con->items[i+1].dict_id.num      ? NULL // note: +1 bc first con element in VCF_COPY_SAMPLE
                              : con->items[i+1].did_i_small < 255 ? CTX(con->items[i+1].did_i_small)
                              :                                     ECTX(con->items[i+1].dict_id);

            set_last_txt_(item_ctx->did_i, sfs[i], sf_lens[i]);

            switch (item_ctx->did_i) {
                case FORMAT_GT: // handle GT for INFO_SF: logic simlar to vcf_piz_container_cb
                    vcf_piz_GT_update_other_fields (vb, sfs[i]);
                    break;

                case FORMAT_DP: 
                    if (!str_is_1chari (sf,i,'.')) {
                        int64_t dp;
                        ASSPIZ (str_get_int (STRi(sf,i), &dp), "Expecting FORMAT/DP to be an integer: %.*s", STRfi(sf,i));
                        ctx_set_last_value (VB, item_ctx, dp);
                    }
                    // fallthrough

                case FORMAT_SB:
                case FORMAT_PS:
                case FORMAT_PID:
                    vcf_piz_con_item_cb (VB, &con->items[i+1], STRi(sf, i));
                    break;

                default: {}
            }
        }
    }

    return HAS_NEW_VALUE;
}

void seg_mux_by_is_prev_sample_copied (VBlockVCFP vb, ZipDataLineVCF *dl, ContextP ctx, Multiplexer2P mux, STRp(value))
{
    ASSERTINRANGE (vb->sample_i, 0, vcf_num_samples);
    int channel_i = SAMPLE_COPIED_SAME_FMT_ZIP;

    ContextP channel_ctx = seg_mux_get_channel_ctx (VB, ctx->did_i, (MultiplexerP)mux, channel_i);

    vcf_seg_field_fallback (vb, channel_ctx, STRa(value));
    seg_by_ctx (VB, STRa(mux->snip), ctx, 0);
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_PREV_COPIED)
{    
    ASSERTINRANGE (vb->sample_i, 0, vcf_num_samples);
    int channel_i = SAMPLE_COPIED_SAME_FMT_PIZ;

    return reconstruct_demultiplex (vb, ctx, STRa(snip), channel_i, new_value, reconstruct);
}  
