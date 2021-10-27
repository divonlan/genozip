// ------------------------------------------------------------------
//   vcf_piz.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <math.h>
#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "dict_id.h"
#include "strings.h"
#include "codec.h"
#include "reconstruct.h"
#include "piz.h"

bool vcf_piz_read_one_vb (VBlock *vb, Section sl)
{ 
    VB_VCF->last_end_line_i = LAST_LINE_I_INIT;
    return true;
}

// determine whether we should reconstruct this VB in LUFT coordinates. this is the is_translation callback defined in TRANSLATIONS
bool vcf_vb_is_luft (VBlockP vb)
{
    if (!vb) return false; // txt header translation is handled by vcf_inspect_txt_header_piz, not a header translator

    return vb->vb_coords == DC_LUFT;
}

// returns true if section is to be skipped reading / uncompressing
bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (flag.drop_genotypes && // note: if all samples are filtered out with --samples then flag.drop_genotypes=true (set in vcf_samples_analyze_field_name_line)
        (dict_id.num == _VCF_FORMAT || dict_id.num == _VCF_SAMPLES || dict_id_is_vcf_format_sf (dict_id)))
        return true;

    if (flag.gt_only && sections_has_dict_id (st) && dict_id_is_vcf_format_sf (dict_id) 
        && dict_id.num != _FORMAT_GT
        && dict_id.num != _FORMAT_GT_HT
        && dict_id.num != _FORMAT_PBWT_RUNS
        && dict_id.num != _FORMAT_PBWT_FGRC
        && dict_id.num != _FORMAT_GT_HT_INDEX)
        return true;

    // if --count, we only need TOPLEVEL and the fields needed for the available filters (--regions)
    if (flag.count && sections_has_dict_id (st) &&
        (     dict_id.num != _VCF_TOPLEVEL && 
              dict_id.num != _VCF_CHROM    && // easier to always have CHROM
             (dict_id.num != _VCF_POS || !flag.regions))) return true;

    return false;
}

static void vcf_piz_replace_pos_with_gpos (VBlockP vb)
{
    DidIType chrom_did_i = (vb->vb_coords == DC_LUFT ? VCF_oCHROM : VCF_CHROM);
    Context *pos_ctx = CTX (vb->vb_coords == DC_LUFT ? VCF_oPOS : VCF_POS);

    WordIndex ref_chrom_index = ref_contigs_get_by_name (gref, last_txt (vb, chrom_did_i), vb->last_txt_len (chrom_did_i), true, true);
    Range *r = ref_get_range_by_ref_index (vb, gref, ref_chrom_index); // possibly NULL

    // remove CHROM and POS and two \t
    vb->txt_data.len -= pos_ctx->last_txt_len + vb->last_txt_len (chrom_did_i) + 2; // remove main CHROM\tPOS\t

    PosType pos = pos_ctx->last_value.i;
    PosType gpos = (r && pos >= r->first_pos && pos <= r->last_pos) ? r->gpos + (pos - r->first_pos) : 0; 

    // if CHROM exists in the reference, POS must fit withing it
    ASSPIZ (!r || (pos >= r->first_pos && pos <= r->last_pos), "POS=%"PRId64" is out of bounds for reference of CHROM: [%"PRId64",%"PRId64"]", 
            pos, r->first_pos, r->last_pos);

    RECONSTRUCT ("GPOS\t", 5); 
    RECONSTRUCT_INT (gpos);
    RECONSTRUCT1('\t');
}

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be reconstructed. contexts are not consumed.
CONTAINER_FILTER_FUNC (vcf_piz_filter)
{
    // before INFO is reconstructed, we save POS.last_value, as it might be changed by an INFO/END. We restore it in the INFO container callback.
    if (item < 0) {
        CTX(VCF_POS)->rback_last_value = CTX(VCF_POS)->last_value; // consumed by regions_is_site_included
        return true; // show repeat as normal
    }
    
    uint64_t dnum = con->items[item].dict_id.num;

    switch (dict_id.num) {

        case _VCF_TOPLEVEL:
        case _VCF_TOPLUFT:
            // --drop-genotypes: remove the two tabs at the end of the line
            if (flag.drop_genotypes && dnum == _VCF_EOL)
                vb->txt_data.len -= 2;

            // in dual-coordinates files - get the COORDS and oSTATUS at the beginning of each line
            else if (dnum == _VCF_oSTATUS || dnum == _VCF_COORDS) {
                if (flag.show_dvcf)
                    return true; // show
                else if (z_dual_coords)
                    *reconstruct = false; // set last_index, without reconstructing
                else
                    return false; // filter out entirely without consuming (non-DC files have no oStatus data)
            }

            // --gpos: after main POS or oPOS reconstructed - re-write POS as GPOS
            else if (item == 4 && flag.gpos) 
                vcf_piz_replace_pos_with_gpos (vb);

            break;

        case _VCF_INFO:
            // for dual-coordinates genozip files - select which DVCF item to show based on flag.luft
            if ((vb->vb_coords == DC_LUFT    && (dnum == _INFO_LUFT || dnum == _INFO_LREJ)) ||
                (vb->vb_coords == DC_PRIMARY && (dnum == _INFO_PRIM || dnum == _INFO_PREJ)))
                return false;
        
            // case: --single-coord - get rid of the DVCF INFO item that would be displayed (this kills the prefixes, but not the values)
            if (flag.single_coord && (dnum == _INFO_PRIM || dnum == _INFO_LUFT || dnum == _INFO_LREJ || dnum == _INFO_PREJ)) {
                if (con_nitems (*con) == 2 && z_dual_coords) RECONSTRUCT1 ('.'); // if this variant's only INFO fields are DVCF - replace with '.'
                *reconstruct = false; // consume but don't reconstruct
            }
            break;

        case _INFO_PRIM: 
        case _INFO_LUFT:
        case _INFO_LREJ:
        case _INFO_PREJ:
            // case: --single-coord - get rid of all of the DVCF INFO item's subitems (needed, since the "reconstruct" in _VCF_INFO for the _INFO_PRIM item, doesn't propagate to the _INFO_PRIM container)
            if (flag.single_coord)
                *reconstruct = false; // consume but don't reconstruct
            break;

        case _VCF_SAMPLES:
            // --GT-only - don't reconstruct non-GT sample subfields
            if (flag.gt_only && dnum != _FORMAT_GT) 
                return false; 
            break;

        default: break;
    }

    return true;    
}

// FORMAT - obey GT-only ; load haplotype line (note: if drop_genotypes we don't reach here, as FORMAT is eliminated by vcf_piz_is_skip_section)
SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT)
{
    bool has_GT = (snip_len>=2 && snip[0]=='G' && snip[1] == 'T' && (snip_len==2 || snip[2] == ':'));
    
    if (reconstruct) {
        if (flag.gt_only) {
            if (has_GT)
                RECONSTRUCT ("GT", 2);
        }
        else 
            RECONSTRUCT (snip, snip_len);
    }

    // initialize haplotype stuff
    if (has_GT && !vb->ht_matrix_ctx) 
        codec_hapmat_piz_calculate_columns (vb);

    return false; // no new value
}

// ------------------------------------
// callback called after reconstruction
// ------------------------------------

static void inline vcf_piz_SAMPLES_subset_samples (VBlockVCFP vb, unsigned rep, unsigned num_reps, int32_t recon_len)
{
    if (!samples_am_i_included (rep))
        vb->txt_data.len -= recon_len + (rep == num_reps - 1); // if last sample, we also remove the preceeding \t (recon_len includes the sample's separator \t, except for the last sample that doesn't have a separator)
}

static void inline vcf_piz_append_ostatus_to_INFO (VBlockP vb)
{
    STR0(snip);
    ContextP ostatus_ctx = CTX (VCF_oSTATUS);
    ctx_get_snip_by_word_index (ostatus_ctx, ostatus_ctx->last_value.i, &snip, &snip_len);
    RECONSTRUCT (";oSTATUS=", 9);
    RECONSTRUCT (snip, snip_len);
}

CONTAINER_CALLBACK (vcf_piz_container_cb)
{
    VBlockVCFP vcf_vb = (VBlockVCFP)vb;
    #define have_INFO_SF  (vcf_vb->sf_snip.len > 0)

    // case: we have an INFO/SF field and we reconstructed and we reconstructed a repeat (i.e. one ht) GT field of a sample 
    if (dict_id.num == _FORMAT_GT && have_INFO_SF) 
        vcf_piz_GT_cb_calc_INFO_SF (vcf_vb, rep, recon, recon_len);

    else if (is_top_level) {

        // case: we have an INFO/SF field and we reconstructed one VCF line
        if (have_INFO_SF) 
            vcf_piz_TOPLEVEL_cb_insert_INFO_SF (vcf_vb); // cleans up allocations - call even if line will be dropped due oSTATUS

        // case: we are reconstructing with --luft and we reconstructed one VCF line
        if (z_dual_coords) 
            vcf_lo_piz_TOPLEVEL_cb_filter_line (vb);

        if (flag.snps_only && !vcf_refalt_piz_is_variant_snp (vb))
            vb->drop_curr_line = "snps_only";

        if (flag.indels_only && !vcf_refalt_piz_is_variant_indel (vb))
            vb->drop_curr_line = "indels_only";
    }

    // if requested, add oSTATUS at thend of INFO
    else if (flag.show_ostatus && dict_id.num == _VCF_INFO) 
        vcf_piz_append_ostatus_to_INFO (vb);

    else if (dict_id.num == _VCF_SAMPLES && flag.samples) 
        vcf_piz_SAMPLES_subset_samples (vcf_vb, rep, num_reps, recon_len);
}
