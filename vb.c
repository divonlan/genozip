// ------------------------------------------------------------------
//   vb.c
//   Copyright (C) 2019-2020 Divon Lan <genozip@blackpawventures.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VariantBlock - i.e. one block of VARIANTS_PER_BLOCK data lines in the vcf file

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "genozip.h"

unsigned vb_num_samples_in_sb(const VariantBlock *vb, unsigned sb_i)
{
    return sb_i < vb->num_sample_blocks-1 ? vb->num_samples_per_block 
                                          : global_num_samples % vb->num_samples_per_block;
} 

unsigned vb_num_sections(VariantBlock *vb) 
{
    return 1 + vb->has_genotype_data + (vb->phase_type == PHASE_MIXED_PHASED) + (vb->num_haplotypes_per_line > 0);
}

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb (VariantBlock **vb_p) 
{
    if (! *vb_p) return; // nothing to release

    VariantBlock *vb = *vb_p;
    *vb_p = NULL;

    for (unsigned i=0; i < VARIANTS_PER_BLOCK; i++) {
        DataLine *dl = &vb->data_lines[i];
         
        dl->line_i = dl->num_subfields = dl->genotype_data.len = 0;
        dl->phase_type = PHASE_UNKNOWN;
        dl->has_haplotype_data = dl->has_genotype_data = 0;

        memset (dl->sf_i, 0, MAX_SUBFIELDS * sizeof(dl->sf_i[0]));

        buf_free(&dl->line);
        buf_free(&dl->variant_data);
        buf_free(&dl->genotype_data);
        buf_free(&dl->haplotype_data);
        buf_free(&dl->phase_data);
    }

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }
            

    vb->num_lines = vb->first_line = vb->variant_block_i = 0;
    vb->vcf_data_size = vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = vb->ready_to_dispatch = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->vcf_file = vb->z_file = NULL;

    vb->min_pos = vb->max_pos = vb->last_pos = vb->chrom[0] = vb->is_sorted_by_pos = 0;

    memset(vb->add_bytes, 0, sizeof(vb->add_bytes));
    memset(vb->vcf_section_bytes, 0, sizeof(vb->vcf_section_bytes));
    memset(vb->z_section_bytes, 0, sizeof(vb->z_section_bytes));

    memset (&vb->profile, 0, sizeof (vb->profile));

    buf_free(&vb->line_variant_data);
    buf_free(&vb->line_gt_data);
    buf_free(&vb->line_ht_data);
    buf_free(&vb->line_phase_data);

    buf_free(&vb->next_gt_in_sample);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->genotype_section_lens_buf);

    buf_free(&vb->optimized_gl_dict);
    buf_free(&vb->variant_data_section_data);
    buf_free(&vb->haplotype_permutation_index);
    buf_free(&vb->gt_sb_line_starts_buf);
    buf_free(&vb->gt_sb_line_lengths_buf);
    buf_free(&vb->helper_index_buf);
    buf_free(&vb->vardata_header_buf);
    buf_free(&vb->compressed);
    buf_free(&vb->z_data);
    buf_free(&vb->z_section_headers);
    buf_free(&vb->ht_columns_data);

    for (unsigned i=0; i < MAX_SUBFIELDS; i++) 
        mtf_free_context (&vb->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->z_file = vb->vcf_file = NULL; // we're not freeing them, must setting the point to null

    vb->in_use = false; // released the VB back into the pool - it may now be reused

    // things that persist:
    
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used
    
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    // between VBs of a file or concatenated files.

}

VariantBlockPool *vb_construct_pool (unsigned num_vbs, unsigned vb_id_prefix)
{
    VariantBlockPool *pool = (VariantBlockPool *)calloc (1, sizeof (VariantBlockPool) + num_vbs * sizeof (VariantBlock)); // note we can't use Buffer yet, because we don't have VBs yet...
    ASSERT0 (pool, "Error: failed to calloc pool");

    pool->num_vbs      = num_vbs;
    pool->vb_id_prefix = vb_id_prefix;

    return pool;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VariantBlock *vb_get_vb(VariantBlockPool *pool, 
                        File *vcf_file, File *z_file,
                        unsigned variant_block_i)
{
    VariantBlock *vb;

    if (!pool) { // should only be used for unit testing - memory-leaks a VB
        vb = calloc (1, sizeof(VariantBlock));
        ASSERT0 (vb, "Error: failed to calloc vb");
    }
    else {
        // see if there's a VB avaiable for recycling
        unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        
            vb = &pool->vb[vb_i];
        
            if (!vb->in_use) {
                vb->id = pool->vb_id_prefix + vb_i;
                break;
            }
        }

        ASSERT (vb_i < pool->num_vbs, "Error: VB pool vb_id_prefix=%u maxed out", pool->vb_id_prefix)
    }

    vb->in_use           = true;
    vb->variant_block_i  = variant_block_i;
    vb->vcf_file         = vcf_file;
    vb->z_file           = z_file;
    vb->is_sorted_by_pos = true;
    vb->min_pos = vb->max_pos = -1;

    return vb;
}

void vb_free_buffer_array (Buffer **buf_array, unsigned buf_array_len)
{
    if (! (*buf_array)) return; // array not allocated - nothing to do

    for (unsigned i=0; i < buf_array_len; i++) 
        if ((*buf_array)[i].memory) free ((*buf_array)[i].memory);

    free (*buf_array);

    *buf_array = NULL;
}

void vb_cleanup_memory(VariantBlockPool *pool)
{
    // see if there's a VB avaiable for recycling
    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VariantBlock *vb = &pool->vb[vb_i];
     
        vb_free_buffer_array (&vb->genotype_sections_data, vb->num_sample_blocks);
        vb_free_buffer_array (&vb->haplotype_sections_data, vb->num_sample_blocks);
        vb_free_buffer_array (&vb->phase_sections_data, vb->num_sample_blocks);
     
        vb->num_sample_blocks = 0;
    }
}
