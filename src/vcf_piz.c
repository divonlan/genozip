// ------------------------------------------------------------------
//   vcf_piz.c
//   Copyright (C) 2019-2025 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#include <math.h>
#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"

void vcf_piz_finalize (bool is_last_z_file)
{
    vcf_header_finalize();
}

void vcf_piz_genozip_header (ConstSectionHeaderGenozipHeaderP header)
{
    if (VER(14)) 
        segconf.has[FORMAT_RGQ] = header->vcf.segconf_has_RGQ;

    if (VER(15)) {
        z_file->max_ploidy_for_mux    = header->vcf.max_ploidy_for_mux; // since 15.0.36
        segconf.FMT_GQ_method         = header->vcf.segconf_GQ_method;
        segconf.FMT_DP_method         = header->vcf.segconf_FMT_DP_method;
        segconf.INFO_DP_method        = header->vcf.segconf_INF_DP_method;
        segconf.MATEID_method         = header->vcf.segconf_MATEID_method;
        segconf.vcf_del_svlen_is_neg  = header->vcf.segconf_del_svlen_is_neg;
        segconf.vcf_sample_copy          = header->vcf.segconf_sample_copy;
        segconf.Q_to_O                = BGEN32F (header->vcf.segconf_Q_to_O);
        segconf.wid[INFO_AC].width    = header->vcf.width.AC;
        segconf.wid[INFO_AF].width    = header->vcf.width.AF;
        segconf.wid[INFO_AN].width    = header->vcf.width.AN;
        segconf.wid[INFO_DP].width    = header->vcf.width.DP;
        segconf.wid[INFO_QD].width    = header->vcf.width.QD;
        segconf.wid[INFO_SF].width    = header->vcf.width.SF;
        segconf.wid[INFO_MLEAC].width = header->vcf.width.MLEAC;
        segconf.wid[INFO_AS_SB_TABLE].width = header->vcf.width.AS_SB_TABLE;
        segconf.wid[VCF_ID].width     = header->vcf.width.ID;
        segconf.wid[VCF_QUAL].width   = header->vcf.width.QUAL;
        segconf.wid[INFO_BaseCounts].width  = header->vcf.width.BaseCounts;
        segconf.wid[INFO_DPB].width   = header->vcf.width.DPB;
    }
}

bool vcf_piz_init_vb (VBlockP vb_, ConstSectionHeaderVbHeaderP header)
{ 
    VBlockVCFP vb = (VBlockVCFP)vb_;

    vb->recon_size = BGEN32 (header->recon_size); 
    
    // HT_n_lines: number of lines that are included in HT_MATRIX. Since 15.0.48, we don't
    // include lines with no FORMAT/GT or lines copied in vcf_seg_sv_SAMPLES. Prior to that, 
    // all lines were always included.
    decl_ctx (FORMAT_GT_HT); 
    if (VER(15)) {
        ctx->HT_n_lines = BGEN32 (header->vcf_HT_n_lines); // since 15.0.48.
        
        // case: vcf_zip_set_vb_header_specific set header->vcf_HT_n_lines, but truly no lines have GT
        if (ctx->HT_n_lines == 0xffffffff) 
            ctx->HT_n_lines = 0;
        
        // case: file is 15.0.0 to 15.0.46 where this field was not set
        else if (ctx->HT_n_lines == 0)  
            goto fallback; 
    }
    
    else fallback: {
        // in older files, we always included all lines in the HT matrix if any line was included
        if (section_get_section (vb->vblock_i, SEC_LOCAL, _FORMAT_PBWT_RUNS)) 
            ctx->HT_n_lines = vb->lines.len32; 
    }

    if (segconf.vcf_sample_copy) 
        vcf_sample_copy_piz_init_vb (vb);

    CTX(INFO_END)->last_end_line_i = LAST_LINE_I_INIT;

    return true; // all good*
}

void vcf_piz_vb_recon_init (VBlockP vb)
{
}

// returns true if section is to be skipped reading / uncompressing
IS_SKIP (vcf_piz_is_skip_section)
{
    if (flag.drop_genotypes && // note: if all samples are filtered out with --samples then flag.drop_genotypes=true (set in vcf_samples_analyze_field_name_line)
        (dict_id.num == _VCF_FORMAT || 
         (dict_id.num == _VCF_SAMPLES && !segconf.has[FORMAT_RGQ]) || // note: if has[RGQ], vcf_piz_special_REFALT peeks SAMPLES so we need it 
         dict_id_is_vcf_format_sf (dict_id)))
        return true;

    if (flag.gt_only && IS_DICTED_SEC (st) && dict_id_is_vcf_format_sf (dict_id) 
        && dict_id.num != _FORMAT_GT
        && dict_id.num != _FORMAT_GT_HT
        && dict_id.num != _FORMAT_PBWT_RUNS
        && dict_id.num != _FORMAT_PBWT_FGRC
        && dict_id.num != _FORMAT_GT_HT_BIG)
        return true;

    // if --count, we only need TOPLEVEL and the fields needed for the available filters (--regions)
    if (flag.count && IS_DICTED_SEC (st) &&
        (     dict_id.num != _VCF_TOPLEVEL && 
              dict_id.num != _VCF_CHROM    && // easier to always have CHROM
             (dict_id.num != _VCF_POS || !flag.regions))) return true;

    return false;
}

// insert a field after following fields have already been reconstructed
// IMPORTANT: if calling this function on a new ctx, add to dids array below 
void vcf_piz_insert_field (VBlockVCFP vb, ContextP ctx, STRp(value))
{
#ifdef DEBUG
    DO_ONCE WARN_IF (segconf.wid[ctx->did_i].width==0, "segconf.wid[%s].width=0. This can legimately happen if this field is not encountered in segconf, or all segconf occurances are longer than %u, or this field didn't exist in version %s, or possibly a bug", ctx->tag_name, SEGCONF_MAX_WIDTH, file_version().s);
#endif

    if (!IS_RECON_INSERTION(ctx)) return;
    
    char *addr = last_txtx (VB, ctx);
    int move_by = (int)value_len - segconf.wid[ctx->did_i].width;

    // actually insert or remove space in txt_data if different than what was reserved
    if (move_by) {
        if (move_by > 0) memmove (addr + move_by, addr, BAFTtxt - addr);             // make room if chars_reserved is not enough
        else             memmove (addr, addr - move_by, BAFTtxt - (addr - move_by)); // shrink room added by vcf_piz_defer if it was too much

        Ltxt += move_by;

        // note: keep txt_data.len 64b to detect bugs
        ASSPIZ (Rtxt, "txt_data overflow: len=%"PRIu64" > size=%"PRIu64". vb->txt_data dumped to %s", 
                vb->txt_data.len, (uint64_t)vb->txt_data.size, txtfile_dump_vb (VB, z_name, NULL).s);
        
        // adjust last_txt of other INFO contexts that might need insertion (and hence last_txt)
        if (ctx->did_i != VCF_ID) { // no need to adjust after inserting ID, as it is inserted during REFALT reconstruction (not at end of TOPLEVEL like the rest)
            Did dids[] = { VCF_QUAL, INFO_QD, INFO_SF, INFO_DP, INFO_AN, INFO_AS_SB_TABLE, INFO_BaseCounts, INFO_DPB };
            uint32_t last_txt_index = ctx->last_txt.index; // note: smaller than index of samples

            bool found_me = false;
            for (int i=0; i < ARRAY_LEN(dids); i++) {
                if (CTX(dids[i])->last_txt.index > last_txt_index)
                    CTX(dids[i])->last_txt.index += move_by;
                
                if (dids[i] == ctx->did_i) found_me = true;
            }

            ASSPIZ (found_me, "Missing ctx=%s in dids[]", ctx->tag_name);   

            // move entries of vb->last_samples set in this line (i.e. for this line's FORMAT)
            if (segconf.vcf_sample_copy) {
                ARRAY (TxtWord, tw, CTX(VCF_SAMPLES)->last_samples);
                uint32_t start = CTX(VCF_FORMAT)->last_value.i * vcf_num_samples;

                for (uint32_t i=start; i < start + vcf_num_samples; i++)
                    if (tw[i].index > last_txt_index)
                        tw[i].index += move_by;
            }
        }

        // ID: adjust REF ALT
        else {
            vb->REF += move_by;
            vb->ALT += move_by;
            CTX(VCF_REFALT)->last_txt.index += move_by;

            for_alt
                alt->alt += move_by;
        }
        
        // adjust lookback addresses that might be affected by this insertion
        vcf_piz_ps_pid_lookback_shift (VB, addr, move_by);
    }
    
    memcpy (addr, value, value_len); // copy
    ctx->last_txt.len = value_len;
}

bool vcf_piz_line_has_RGQ (VBlockVCFP vb)
{
    decl_ctx (FORMAT_RGQ);

    if (ctx->line_has_RGQ == unknown)
        ctx->line_has_RGQ = segconf.has[FORMAT_RGQ] && 
                            container_peek_has_item (VB, CTX(VCF_SAMPLES), (DictId)_FORMAT_RGQ, false); // note: segconf.has[FORMAT_RGQ] is never set in PIZ prior to v14

    return ctx->line_has_RGQ;
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_MUX_BY_HAS_RGQ)
{
    return reconstruct_demultiplex (vb, ctx, STRa(snip), vcf_piz_line_has_RGQ (VB_VCF), new_value, reconstruct);
}

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be reconstructed. contexts are not consumed.
CONTAINER_FILTER_FUNC (vcf_piz_filter)
{
    uint64_t dnum = (item >= 0) ? con->items[item].dict_id.num : 0;

    switch (dict_id.num) {

        case _VCF_TOPLEVEL:
            // After reconstructing POS (next field is ID) - save POS.last_value, as it might be changed by an INFO/END. 
            if (dnum == _VCF_ID) {
                CTX(VCF_POS)->pos_last_value = CTX(VCF_POS)->last_value.i; // consumed by regions_is_site_included
                break;
            }

            else if (dnum == _VCF_QUAL) {
                if (!VER(15)) vcf_piz_refalt_parse (VB_VCF); // starting 14.0.12, called from vcf_piz_con_item_cb (no harm if called twice)
            }

            // --drop-genotypes: remove the two tabs at the end of the line
            else if (dnum == _VCF_EOL) { 
                if (flag.drop_genotypes) Ltxt -= 2;
            }

            // --drop-genotypes: drop SAMPLES (might be loaded if has[RGQ], as needed for vcf_piz_special_REFALT)
            else if (dnum == _VCF_SAMPLES) {
                if (flag.drop_genotypes) return false;
            }

            // DVCF-related fields that exist if files compressed 12.0.0 - 15.0.41
            else if (dnum == _VCF_oSTATUS || dnum == _VCF_COORDS) 
                return false; // filter out entirely without consuming 
        
            break;
            
        case _VCF_SAMPLES:
        case _VCF_SAMPLES_0: // "no mate" channel of SAMPLES multiplexor
            // Set sample_i before reconstructing each sample
            if (item == -1) {
                vb->sample_i = rep;
                
                if (segconf.vcf_sample_copy)  
                    CTX(VCF_COPY_SAMPLE)->last_value.i = 0;
            }

            // case: sample as copied from previous line - filter out all items following the VCF_COPY_SAMPLE
            else if (CTX(VCF_COPY_SAMPLE)->last_value.i)
                return false;

            // --GT-only - don't reconstruct non-GT sample subfields
            else if (flag.gt_only && dnum != _FORMAT_GT) 
                return false; 

            break;

        default: break;
    }

    return true;    

}

// called from vcf_piz_refalt_parse
void vcf_piz_insert_VCF_ID (VBlockVCFP vb)
{
    // case: ID_is_variant
    if (CTX(VCF_ID)->deferred_snip.len32) { 
        decl_ctx (VCF_ID);
        
        if (IS_RECON_INSERTION(ctx)) {
            rom id_recon = BAFTtxt;
            char sep = *B1STc(CTX(VCF_ID)->deferred_snip);

            RECONSTRUCT_str (vb->chrom_name);
            RECONSTRUCT1 (sep);
            RECONSTRUCT_INT (CTX(VCF_POS)->last_value.i);
            RECONSTRUCT1 (sep);
            RECONSTRUCT_str (vb->REF);
            RECONSTRUCT1 (sep);
            RECONSTRUCT_str (vb->ALT);
            
            uint32_t id_len = BAFTtxt - id_recon;
            Ltxt -= id_len; 

            char id[id_len];
            memcpy (id, id_recon, id_len);
            
            vcf_piz_insert_field (vb, ctx, STRa(id));
        }
    }

    // case: deferred with no snip = PBSV
    else
        vcf_piz_insert_pbsv_ID (vb);    
}

// ------------------------------------
// callback called after reconstruction
// ------------------------------------

static void inline vcf_piz_SAMPLES_subset_samples (VBlockVCFP vb, unsigned rep, unsigned num_reps, int32_t recon_len)
{
    if (!samples_am_i_included (rep))
        Ltxt -= recon_len + (rep == num_reps - 1); // if last sample, we also remove the preceeding \t (recon_len includes the sample's separator \t, except for the last sample that doesn't have a separator)
}

CONTAINER_ITEM_CALLBACK (vcf_piz_con_item_cb)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    switch (con_item->dict_id.num) {
        // IMPORTANT: when adding a "case", also update set CI1_ITEM_CB in vcf_seg_FORMAT
        
        // occurs if GT was segged as a snip (see vcf_seg_FORMAT) (if pbwt, it is done in vcf_piz_container_cb)
        case _FORMAT_GT:
            vcf_piz_GT_update_other_fields (vb, recon);
            break;

        case _FORMAT_DP:
            if (ctx_has_value (VB, FORMAT_DP)) { // not '.' or missing
                if (segconf.INFO_DP_method == BY_FORMAT_DP) 
                    CTX(INFO_DP)->dp.sum_format_dp += CTX(FORMAT_DP)->last_value.i;

                // add up DP's of samples with GT!=0/0, for consumption by INFO/QD predictor
                QdPredType pd = CTX(INFO_QD)->qd.pred_type;
                if (pd == QD_PRED_SUM_DP || pd == QD_PRED_SUM_DP_P001 || pd == QD_PRED_SUM_DP_M001)
                    vcf_piz_sum_DP_for_QD (vb, STRa(recon));
            }
            break;
            
        case _FORMAT_SB:
            vcf_sum_SB_for_AS_SB_TABLE (vb, STRa(recon));
            break;
            
        case _FORMAT_PS: // since v13: PS has item_cb
            vcf_piz_ps_pid_lookback_insert (vb, FORMAT_PS, STRa(recon));
            break;

        case _FORMAT_PID: // since v13: PID has item_cb
            vcf_piz_ps_pid_lookback_insert (vb, FORMAT_PID, STRa(recon));
            break;
        
        case _VCF_REFALT: // files compressed starting v14.0.12
            vcf_piz_refalt_parse (vb);
            break;

        default:
            ABORT_PIZ ("vcf_piz_con_item_cb doesn't know how to handle dict_id=%s. %s", 
                       dis_dict_id (con_item->dict_id).s, genozip_update_msg());
    }
}

// called from toplevel callback
static void vcf_piz_insert_by_snip (VBlockVCFP vb, ContextP ctx)
{
    rom recon = BAFTtxt;

    reconstruct_one_snip (VB, ctx, WORD_INDEX_NONE, STRb(ctx->deferred_snip), true, __FUNCLINE);
    uint32_t recon_len = BAFTtxt - recon;

    // we can't send recon in txt_data since its going to memmove, so we copy it to a buffer
    ctx->deferred_snip.len32 = 0;
    buf_add_moreS (vb, &ctx->deferred_snip, recon, "ctx->deferred_snip");
    
    vcf_piz_insert_field (vb, ctx, STRb(ctx->deferred_snip));
    Ltxt -= recon_len;
}

// defer reconstruction after reconstruction of FORMAT/SB - happens in vcf_piz_insert_INFO_AS_SB_TABLE
SPECIAL_RECONSTRUCTOR (vcf_piz_special_DEFER)
{
    // save snip for later (note: the SNIP_SPECIAL+code are already removed)
    if (snip_len) {
        ctx->deferred_snip.len32 = 0;
        buf_add_moreS (vb, &ctx->deferred_snip, snip, "contexts->deferred_snip");
        *BAFTc (ctx->deferred_snip) = 0; // buf_add_moreS allocates room for this nul-terminator
    }
    
    vcf_piz_defer (ctx);

    return NO_NEW_VALUE;
}

CONTAINER_CALLBACK (vcf_piz_container_cb)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    // note: occurs if GT is a container, not if it is segged as a snip (in which case it happens in vcf_piz_con_item_cb)
    if (dict_id.num == _FORMAT_GT) {
        // after first HT is reconstructed, re-write the predictable phase (note: predictable phases are only used for ploidy=2)
        if (rep==0 && con->repsep[0] == '&') 
            vcf_piz_FORMAT_GT_rewrite_predicted_phase (vb, STRa (recon));
   
        if (rep==1 && vb->flags.vcf.use_null_DP_method && con->repeats==2 && 
            con->items[0].separator[1] != CI1_ITEM_PRIVATE/*override flag*/) 
            vcf_piz_GT_cb_null_GT_if_null_DP (vb, recon);

        if (rep==0/*first HT*/) 
            vcf_piz_GT_update_other_fields (vb, recon);
    }

    else if (is_top_level) {
        // insert fields whose value is determined by the sample fields that we just finished reconstructing    
        if (is_deferred(VCF_QUAL))         vcf_piz_insert_QUAL_by_GP (vb); 
        if (is_deferred(INFO_BaseCounts))  vcf_piz_insert_INFO_BaseCounts_by_AD (vb); 
        if (is_deferred(INFO_DP))          vcf_piz_insert_INFO_DP (vb);    // constraint: must be after BaseCounts (might depends on BaseCounts.last_value)
        if (is_deferred(INFO_DPB))         vcf_piz_insert_by_snip (vb, CTX(INFO_DPB)); // constraint: must be after DP (depends on DP.last_value)
        if (is_deferred(INFO_SF))          vcf_piz_insert_INFO_SF (vb);    
        if (is_deferred(INFO_QD))          vcf_piz_insert_INFO_QD (vb);    // constraint: must be inserted after QUAL (it might depend on it)
        if (is_deferred(INFO_AS_SB_TABLE)) vcf_piz_insert_INFO_AS_SB_TABLE (vb);

        if (flag.snps_only && !vcf_refalt_piz_is_variant_snp (vb))
            vb->drop_curr_line = "snps_only";

        if (flag.indels_only && !vcf_refalt_piz_is_variant_indel (vb))
            vb->drop_curr_line = "indels_only";
    }

    // store reference to sample, in case the same sample in the next line needs to copy it
    else if (dict_id.num == _VCF_SAMPLES) {
        ctx_set_last_value (VB, CTX(VCF_LOOKBACK), (ValueType){ .i = con->repeats });

        if (segconf.vcf_sample_copy)  // since 15.0.69
            vcf_copy_sample_piz_store (vb, STRa(recon));

        if (flag.samples) 
            vcf_piz_SAMPLES_subset_samples (vb, rep, con->repeats, recon_len);
    }
}
