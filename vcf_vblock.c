// ------------------------------------------------------------------
//   vcf_vblock.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "vcf_private.h"

unsigned vcf_vb_size (void) { return sizeof (VBlockVCF); }
unsigned vcf_vb_zip_dl_size (void) { return sizeof (ZipDataLineVCF); }
bool vcf_vb_has_haplotype_data (VBlockP vb) { return ((VBlockVCFP)vb)->has_haplotype_data; }

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
void vcf_vb_release_vb (VBlockVCF *vb) 
{
    vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->max_gt_line_len = vb->max_genotype_section_len = 0;
    vb->ac = vb->an = vb->af = NULL;
    vb->ac_len = vb->an_len = vb->af_len = 0;
    vb->is_af_before_ac = vb->is_an_before_ac = false;

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }

    buf_free(&vb->line_gt_data);
    buf_free(&vb->line_ht_data);
    buf_free(&vb->line_phase_data);
    buf_free(&vb->sample_iterator);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->is_sb_included);
    buf_free(&vb->genotype_section_lens_buf);
    buf_free(&vb->format_mapper_buf);
    buf_free(&vb->optimized_gl_dict);
    buf_free(&vb->haplotype_permutation_index);
    buf_free(&vb->haplotype_permutation_index_squeezed);
    buf_free(&vb->gt_sb_line_starts_buf);
    buf_free(&vb->gt_sb_line_lengths_buf);
    buf_free(&vb->helper_index_buf);
    buf_free(&vb->ht_columns_data);
    buf_free(&vb->column_of_zeros);
    buf_free(&vb->gtshark_db_db_data);
    buf_free(&vb->gtshark_db_gt_data);
    buf_free(&vb->gtshark_exceptions_line_i);
    buf_free(&vb->gtshark_exceptions_ht_i);
    buf_free(&vb->gtshark_exceptions_allele);
    buf_free(&vb->gtshark_vcf_data);
}

void vcf_vb_destroy_vb (VBlockVCF *vb)
{
    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_destroy (&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_destroy (&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_destroy (&vb->phase_sections_data[i]);
    }

    buf_destroy (&vb->line_gt_data);
    buf_destroy (&vb->line_ht_data);
    buf_destroy (&vb->line_phase_data);
    buf_destroy (&vb->sample_iterator);
    buf_destroy (&vb->genotype_one_section_data);
    buf_destroy (&vb->is_sb_included);
    buf_destroy (&vb->genotype_section_lens_buf);
    buf_destroy (&vb->format_mapper_buf);
    buf_destroy (&vb->optimized_gl_dict);
    buf_destroy (&vb->haplotype_permutation_index);
    buf_destroy (&vb->haplotype_permutation_index_squeezed);
    buf_destroy (&vb->gt_sb_line_starts_buf);
    buf_destroy (&vb->gt_sb_line_lengths_buf);
    buf_destroy (&vb->helper_index_buf);
    buf_destroy (&vb->ht_columns_data);
    buf_destroy (&vb->format_mapper_buf);
    buf_destroy (&vb->column_of_zeros);
    buf_destroy (&vb->gtshark_db_db_data);
    buf_destroy (&vb->gtshark_db_gt_data);
    buf_destroy (&vb->gtshark_exceptions_line_i);
    buf_destroy (&vb->gtshark_exceptions_ht_i);
    buf_destroy (&vb->gtshark_exceptions_allele);
    buf_destroy (&vb->gtshark_vcf_data);
}

// freeing arrays of buffers allocated by calloc 
static void vcf_vb_free_buffer_array (Buffer **buf_array, unsigned buf_array_len)
{
    if (! (*buf_array)) return; // array not allocated - nothing to do

    for (unsigned i=0; i < buf_array_len; i++) 
        buf_destroy (&(*buf_array)[i]);

    FREE (*buf_array);
    *buf_array = NULL;
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vcf_vb_cleanup_memory (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    if (vb_ == evb) return; // nothing to cleanup in evb with VCF

    vcf_vb_free_buffer_array (&vb->genotype_sections_data, vb->num_sample_blocks);
    vcf_vb_free_buffer_array (&vb->haplotype_sections_data, vb->num_sample_blocks);
    vcf_vb_free_buffer_array (&vb->phase_sections_data, vb->num_sample_blocks);
    vb->num_sample_blocks = 0;

    global_vcf_num_samples = 0;
}


unsigned vcf_vb_num_samples_in_sb (const VBlockVCF *vb, unsigned sb_i)
{
    // case: last block has less than a full block of samples
    if (sb_i == vb->num_sample_blocks-1 && global_vcf_num_samples % vb->num_samples_per_block)
        return global_vcf_num_samples % vb->num_samples_per_block;

    else
        return vb->num_samples_per_block;
} 
