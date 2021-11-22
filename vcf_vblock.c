// ------------------------------------------------------------------
//   vcf_vblock.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "vcf_private.h"
#include "strings.h"

unsigned vcf_vb_size (DataType dt) { return sizeof (VBlockVCF); }
unsigned vcf_vb_zip_dl_size (void) { return sizeof (ZipDataLineVCF); }
bool vcf_vb_has_haplotype_data (VBlockP vb) { return !!VB_VCF->ht_matrix_ctx; }

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
void vcf_vb_release_vb (VBlockVCF *vb) 
{
    vb->ploidy = 0;
    vb->use_special_sf = 0;
    vb->gt_prev_ploidy = vb->num_dps_this_line = 0;
    vb->gt_prev_phase = 0;
    vb->main_refalt = NULL;
    vb->main_ref_len = vb->main_alt_len = 0;
    vb->last_end_line_i = 0;
    memset (vb->ad_values, 0, sizeof (vb->ad_values));
    vb->new_ref = 0;
    vb->is_del_sv = 0;
    vb->vcf_version = 0;
    vb->PL_mux_by_DP = 0;
    vb->PS_encountered_last_line = 0;
    vb->sum_dp_this_line = 0;

    memset (&vb->mux_PLn,    0, sizeof(vb->mux_PLn));
    memset (&vb->mux_GL,     0, sizeof(vb->mux_GL));
    memset (&vb->mux_GP,     0, sizeof(vb->mux_GP));
    memset (&vb->mux_PRI,    0, sizeof(vb->mux_PRI));
    memset (&vb->mux_DS,     0, sizeof(vb->mux_DS));
    memset (&vb->mux_PP,     0, sizeof(vb->mux_PP));
    memset (&vb->mux_PVAL,   0, sizeof(vb->mux_PVAL));
    memset (&vb->mux_FREQ,   0, sizeof(vb->mux_FREQ));
    memset (&vb->mux_RD,     0, sizeof(vb->mux_RD));
    memset (&vb->mux_GQ,     0, sizeof(vb->mux_GQ));
    memset (&vb->mux_AD,     0, sizeof(vb->mux_AD));
    memset (&vb->mux_ADALL,  0, sizeof(vb->mux_ADALL));
    memset (&vb->mux_PLy,    0, sizeof(vb->mux_PLy));
    
    buf_free (&vb->sf_txt);
    buf_free (&vb->sf_snip);
    buf_free (&vb->hapmat_helper_index_buf);
    buf_free (&vb->hapmat_columns_data);
    buf_free (&vb->hapmat_one_array);
    buf_free (&vb->hapmat_column_of_zeros);
    buf_free (&vb->format_mapper_buf);
    buf_free (&vb->format_contexts);
    buf_free (&vb->info_items);
    buf_free (&vb->tags);
    buf_free (&vb->rejects_report);
    buf_free (&vb->last_format);
}

void vcf_vb_destroy_vb (VBlockVCF *vb)
{
    buf_destroy (&vb->sf_txt);
    buf_destroy (&vb->sf_snip);
    buf_destroy (&vb->hapmat_helper_index_buf);
    buf_destroy (&vb->hapmat_columns_data);
    buf_destroy (&vb->hapmat_one_array);
    buf_destroy (&vb->hapmat_column_of_zeros);
    buf_destroy (&vb->format_mapper_buf);
    buf_destroy (&vb->format_contexts);
    buf_destroy (&vb->info_items);
    buf_destroy (&vb->tags);
    buf_destroy (&vb->rejects_report);
    buf_destroy (&vb->last_format);
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vcf_vb_cleanup_memory (VBlock *vb_)
{
    vcf_num_samples = 0;
}

