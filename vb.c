// ------------------------------------------------------------------
//   vb.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VariantBlock - i.e. one block of VARIANTS_PER_BLOCK data lines in the vcf file

#include "genozip.h"
#include "vb.h"
#include "move_to_front.h"

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

    if (vb->data_lines) {
        for (unsigned i=0; i < global_max_lines_per_vb; i++) {
            DataLine *dl = &vb->data_lines[i];
            
            dl->line_i = dl->num_subfields = dl->genotype_data.len = 0;
            dl->phase_type = PHASE_UNKNOWN;
            dl->has_haplotype_data = dl->has_genotype_data = 0;

            memset (dl->sf_i, 0, MAX_DICTS * sizeof(dl->sf_i[0]));

            buf_free(&dl->line);
            buf_free(&dl->variant_data);
            buf_free(&dl->genotype_data);
            buf_free(&dl->haplotype_data);
            buf_free(&dl->phase_data);
        }

        // note: vb->data_line is not freed but rather used by subsequent vbs
        // all VBs for a specific file (including concatenated files) have the same number of data_lines which is global_max_lines_per_vb
    }

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }
            

    vb->num_lines = vb->first_line = vb->variant_block_i = 0;
    vb->vb_data_size = vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = vb->ready_to_dispatch = vb->is_processed = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->vcf_file = vb->z_file = NULL;
    vb->z_next_header_i = 0;
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
    buf_free(&vb->spiced_pw);
    buf_free(&vb->subfields_start_buf);        
    buf_free(&vb->subfields_len_buf);
    buf_free(&vb->num_subfields_buf);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->mtf_ctx[i].dict_id.num)
            mtf_free_context (&vb->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->z_file = vb->vcf_file = NULL; // we're not freeing them, must setting the point to null

    vb->in_use = false; // released the VB back into the pool - it may now be reused
    // things that persist:
    
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used
    
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    // between VBs of a file or concatenated files.

    // vb->column_of_zeros : we don't free this as its a constant array of zeros, of size global_max_lines_per_vb

}

static VariantBlockPool *pools[NUM_POOLS] = {NULL, NULL}; // the pools remains even between vcf files

void vb_create_pool (PoolId pool_id, unsigned num_vbs)
{
    ASSERT (!pools[pool_id] || num_vbs==pools[pool_id]->num_vbs, 
            "Error: pool %u already exists, but with the wrong number of vbs - expected %u but it has %u", pool_id, num_vbs, pools[pool_id]->num_vbs);

    if (!pools[pool_id])  {
        pools[pool_id] = (VariantBlockPool *)calloc (1, sizeof (VariantBlockPool) + num_vbs * sizeof (VariantBlock)); // note we can't use Buffer yet, because we don't have VBs yet...
        ASSERT0 (pools[pool_id], "Error: failed to calloc pool");

        pools[pool_id]->num_vbs = num_vbs; 
    }
}

VariantBlockPool *vb_get_pool (PoolId pool_id)
{
    ASSERT ((int)pool_id < NUM_POOLS && pools[pool_id], "Error: pool %u doesn't exists", pool_id);
    return pools[pool_id];
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VariantBlock *vb_get_vb (PoolId pool_id, 
                         FileP vcf_file, FileP z_file,
                         unsigned variant_block_i)
{
    VariantBlock *vb=NULL;

    if (pool_id == POOL_ID_UNIT_TEST) { // should only be used for unit testing - memory-leaks a VB
        vb = (VariantBlock *)calloc (1, sizeof(VariantBlock));
        ASSERT0 (vb, "Error: failed to calloc vb");
    }
    else {
        // see if there's a VB avaiable for recycling
        unsigned vb_i; for (vb_i=0; vb_i < pools[pool_id]->num_vbs; vb_i++) {
        
            vb = &pools[pool_id]->vb[vb_i];
        
            if (!vb->in_use) {
                vb->id = vb_i;
                vb->pool_id = pool_id;
                break;
            }
        }

        ASSERT (vb_i < pools[pool_id]->num_vbs, "Error: VB pool pool_id=%u is full - it already has %u VBs", pool_id, pools[pool_id]->num_vbs)
    }

    vb->in_use           = true;
    vb->variant_block_i  = variant_block_i;
    vb->vcf_file         = vcf_file;
    vb->z_file           = z_file;
    vb->is_sorted_by_pos = true;
    vb->min_pos = vb->max_pos = -1;

    return vb;
}

// freeing arrays of buffers allocated by calloc 
static void vb_free_buffer_array (VariantBlock *vb, Buffer **buf_array, unsigned buf_array_len)
{
    if (! (*buf_array)) return; // array not allocated - nothing to do

    for (unsigned i=0; i < buf_array_len; i++) 
        buf_destroy (vb, &(*buf_array)[i]);

    free (*buf_array);

    *buf_array = NULL;
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vb_cleanup_memory (PoolId pool_id)
{
    // see if there's a VB avaiable for recycling
    for (unsigned vb_i=0; vb_i < pools[pool_id]->num_vbs; vb_i++) {
        VariantBlock *vb = &pools[pool_id]->vb[vb_i];
     
        vb_free_buffer_array (vb, &vb->genotype_sections_data, vb->num_sample_blocks);
        vb_free_buffer_array (vb, &vb->haplotype_sections_data, vb->num_sample_blocks);
        vb_free_buffer_array (vb, &vb->phase_sections_data, vb->num_sample_blocks);

        if (vb->data_lines) {
            for (unsigned i=0; i < global_max_lines_per_vb; i++) {
                DataLine *dl = &vb->data_lines[i];
                
                buf_destroy(vb, &dl->line);
                buf_destroy(vb, &dl->variant_data);
                buf_destroy(vb, &dl->genotype_data);
                buf_destroy(vb, &dl->haplotype_data);
                buf_destroy(vb, &dl->phase_data);
            }
            free (vb->data_lines);
            vb->data_lines = NULL;
        }

        vb->num_sample_blocks = 0;
        buf_free (&vb->column_of_zeros); // a constant array of 0 of length global_max_lines_per_vb - free only when global_max_lines_per_vb changes
    }

    // note: this function is never called in test mode because it is always is_last_file. this is lucky, bc both sides use these globals
    global_max_lines_per_vb = 0;
    global_num_samples = 0;
}
