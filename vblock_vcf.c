// ------------------------------------------------------------------
//   vblock_vcf.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "vblock.h"

unsigned vb_vcf_num_samples_in_sb (const VBlockVCF *vb, unsigned sb_i)
{
    // case: last block has less than a full block of samples
    if (sb_i == vb->num_sample_blocks-1 && global_num_samples % vb->num_samples_per_block)
        return global_num_samples % vb->num_samples_per_block;

    else
        return vb->num_samples_per_block;
} 

unsigned vb_vcf_num_sections(VBlockVCF *vb) 
{
    return 1 + vb->has_genotype_data + (vb->phase_type == PHASE_MIXED_PHASED) + (vb->num_haplotypes_per_line > 0);
}


// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
void vb_vcf_release_vb (VBlockVCF *vb) 
{
    // note: vb->data_line is not freed but rather used by subsequent vbs
    if (command == ZIP && vb->data_lines.zip)
        memset (vb->data_lines.zip, 0, sizeof(ZipDataLine) * vb->num_data_lines_allocated);

    else if (command != ZIP && vb->data_lines.piz) {
        for (unsigned i=0; i < vb->num_data_lines_allocated; i++) {
            PizDataLine *dl = &vb->data_lines.piz[i];
            
            dl->has_haplotype_data = dl->has_genotype_data = 0;
            dl->format_mtf_i = 0;

            buf_free(&dl->line);
            buf_free(&dl->v1_variant_data);
        }
    }

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }

    vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->last_pos = 0;
    vb->curr_ra_ent = NULL; 
    vb->curr_ra_ent_is_initialized = false;

    buf_free(&vb->line_variant_data);
    buf_free(&vb->line_gt_data);
    buf_free(&vb->line_ht_data);
    buf_free(&vb->line_phase_data);

    buf_free(&vb->sample_iterator);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->is_sb_included);
    buf_free(&vb->genotype_section_lens_buf);

    buf_free (&vb->format_mapper_buf);
    buf_free (&vb->iname_mapper_buf);

    vb->num_info_subfields = vb->num_format_subfields = 0;

    buf_free(&vb->optimized_gl_dict);
    buf_free(&vb->haplotype_permutation_index);
    buf_free(&vb->haplotype_permutation_index_squeezed);
    
    buf_free(&vb->gt_sb_line_starts_buf);
    buf_free(&vb->gt_sb_line_lengths_buf);
    buf_free(&vb->helper_index_buf);
    buf_free(&vb->ht_columns_data);
    buf_free(&vb->format_info_buf);
    buf_free(&vb->ra_buf);
    buf_free(&vb->region_ra_intersection_matrix);
    buf_free(&vb->column_of_zeros);

    buf_free(&vb->gtshark_db_db_data);
    buf_free(&vb->gtshark_db_gt_data);
    buf_free(&vb->gtshark_exceptions_line_i);
    buf_free(&vb->gtshark_exceptions_ht_i);
    buf_free(&vb->gtshark_exceptions_allele);
    buf_free(&vb->gtshark_vcf_data);

    // backward compatibility with genozip v1
    buf_free(&vb->v1_subfields_start_buf);        
    buf_free(&vb->v1_subfields_len_buf);
    buf_free(&vb->v1_num_subfields_buf);
    buf_free(&vb->v1_variant_data_section_data);
}

// freeing arrays of buffers allocated by calloc 
static void vb_vcf_free_buffer_array (Buffer **buf_array, unsigned buf_array_len)
{
    if (! (*buf_array)) return; // array not allocated - nothing to do

    for (unsigned i=0; i < buf_array_len; i++) 
        buf_destroy (&(*buf_array)[i]);

    FREE (*buf_array);

    *buf_array = NULL;
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vb_vcf_cleanup_memory (VBlockVCF *vb)
{
    vb_vcf_free_buffer_array (&vb->genotype_sections_data, vb->num_sample_blocks);
    vb_vcf_free_buffer_array (&vb->haplotype_sections_data, vb->num_sample_blocks);
    vb_vcf_free_buffer_array (&vb->phase_sections_data, vb->num_sample_blocks);
    vb->num_sample_blocks = 0;
}
