// ------------------------------------------------------------------
//   vcf_vblock.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "vcf_private.h"
#include "strings.h"

unsigned vcf_vb_size (void) { return sizeof (VBlockVCF); }
unsigned vcf_vb_zip_dl_size (void) { return sizeof (ZipDataLineVCF); }
bool vcf_vb_has_haplotype_data (VBlockP vb) { return !!((VBlockVCFP)vb)->ht_matrix_ctx; }

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
void vcf_vb_release_vb (VBlockVCF *vb) 
{
    vb->ploidy = 0;
    vb->use_special_sf = 0;
    vb->gt_prev_ploidy = 0;
    vb->gt_prev_phase = 0;
    vb->sf_ctx = NULL;
    vb->main_refalt = NULL;
    vb->main_ref_len = vb->main_alt_len = 0;
    vb->last_end_line_i = 0;
    
    buf_free (&vb->sf_txt);
    buf_free (&vb->sf_snip);
    buf_free (&vb->hapmat_helper_index_buf);
    buf_free (&vb->hapmat_columns_data);
    buf_free (&vb->hapmat_one_array);
    buf_free (&vb->hapmat_column_of_zeros);
    buf_free (&vb->format_mapper_buf);
    buf_free (&vb->info_items);
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
    buf_destroy (&vb->info_items);
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vcf_vb_cleanup_memory (VBlock *vb_)
{
    vcf_num_samples = 0;
}

