// ------------------------------------------------------------------
//   vcf_vblock.c
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "vcf_private.h"
#include "strings.h"

unsigned vcf_vb_size (DataType dt) { return sizeof (VBlockVCF); }
unsigned vcf_vb_zip_dl_size (void) { return sizeof (ZipDataLineVCF); }
bool vcf_vb_has_haplotype_data (VBlockP vb) { return !!VB_VCF->ht_matrix_ctx; }

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
void vcf_vb_release_vb (VBlockVCFP vb) 
{
    vb->ploidy = 0;
    vb->use_special_sf = 0;
    vb->main_ref = vb->main_alt = NULL;
    vb->main_ref_len = vb->main_alt_len = 0;
    memset (vb->ad_values, 0, sizeof (vb->ad_values));
    vb->new_ref = 0;
    vb->is_del_sv = 0;
    vb->vcf_version = 0;
    vb->PL_mux_by_DP = 0;
    vb->pos_aln_i = 0;
    vb->vb_coords = vb->line_coords = 0;
    vb->is_rejects_vb = vb->is_unsorted[0] = vb->is_unsorted[1] = false;    
    vb->recon_size_luft = vb->reject_bytes = 0;
    vb->sort = false;
    vb->first_line = 0;
    vb->line_has_RGQ = 0;

    memset (&vb->first_mux, 0, (char*)&vb->after_mux - (char*)&vb->first_mux); // all mux's 
    
    buf_free (vb->sf_txt);
    buf_free (vb->sf_snip);
    buf_free (vb->hapmat_helper_index_buf);
    buf_free (vb->hapmat_columns_data);
    buf_free (vb->hapmat_one_array);
    buf_free (vb->hapmat_column_of_zeros);
    buf_free (vb->format_mapper_buf);
    buf_free (vb->format_contexts);
    buf_free (vb->info_items);
    buf_free (vb->tags);
    buf_free (vb->rejects_report);
    buf_free (vb->last_format);
}

void vcf_vb_destroy_vb (VBlockVCFP vb)
{
    buf_destroy (vb->sf_txt);
    buf_destroy (vb->sf_snip);
    buf_destroy (vb->hapmat_helper_index_buf);
    buf_destroy (vb->hapmat_columns_data);
    buf_destroy (vb->hapmat_one_array);
    buf_destroy (vb->hapmat_column_of_zeros);
    buf_destroy (vb->format_mapper_buf);
    buf_destroy (vb->format_contexts);
    buf_destroy (vb->info_items);
    buf_destroy (vb->tags);
    buf_destroy (vb->rejects_report);
    buf_destroy (vb->last_format);
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vcf_vb_cleanup_memory (VBlockP vb_)
{
    vcf_num_samples = 0;
}

// this is also used to determine whether we should reconstruct this VB in LUFT coordinates - the is_translation callback defined in TRANSLATIONS
bool vcf_vb_is_luft    (VBlockP vb) { return vb && VB_VCF->vb_coords == DC_LUFT; }
bool vcf_vb_is_primary (VBlockP vb) { return vb && VB_VCF->vb_coords == DC_PRIMARY; }

int32_t vcf_vb_get_reject_bytes (VBlockP vb) { return VB_VCF->reject_bytes; }

rom vcf_coords_name (int coord)
{
    static rom coords_names[4] = { "NONE", "PRIM", "LUFT", "BOTH" };
    
    return (coord < 0 || coord >= NUM_COORDS) ? "(invalid coord)" : coords_names[coord];
}

// ZIP/PIZ: called before segging / reconstructing each line
void vcf_reset_line (VBlockP vb_)
{
    VBlockVCFP vb = (VBlockVCFP)vb_;

    vb->sample_i = 0;
    vb->line_has_RGQ = RGQ_UNKNOWN;

    CTX(INFO_DP)->sum_dp_this_line = 0;
    CTX(INFO_DP)->is_initialized = false;        

    CTX(INFO_QD)->qd.pred_type = QD_PRED_NONE; // sum_dp_with_dosage and pred_type
    CTX(INFO_QD)->qd.sum_dp_with_dosage = 0;

    if (IS_ZIP) {
        for (Did did_i=0; did_i < NUM_VCF_FIELDS; did_i++)
            CTX(did_i)->sf_i = -1; // initialize
    }
}
