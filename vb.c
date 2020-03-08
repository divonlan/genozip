// ------------------------------------------------------------------
//   vb.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VariantBlock - i.e. one block of VARIANTS_PER_BLOCK data lines in the vcf file

#include "genozip.h"
#include "vb.h"
#include "move_to_front.h"

// pool of VBs allocated based on number of threads
static VariantBlockPool *pool = NULL;

// one VB outside of pool
VariantBlock *evb = NULL;

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
            
            dl->vcf_line_i = dl->genotype_data.len = 0;
            dl->phase_type = PHASE_UNKNOWN;
            dl->has_haplotype_data = dl->has_genotype_data = 0;
            dl->format_mtf_i = dl->info_mtf_i = 0;

            buf_free(&dl->line);
            buf_free(&dl->v1_variant_data);
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

    vb->num_lines = vb->first_line = vb->variant_block_i = vb->vcf_data_next_offset = 0;
    vb->vb_data_size = vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = vb->ready_to_dispatch = vb->is_processed = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->vcf_file = vb->z_file = NULL;
    vb->z_next_header_i = 0;
    vb->num_dict_ids = vb->num_format_subfields = 0;
    vb->last_pos = 0;
    vb->curr_ra_ent = NULL; 
    vb->curr_ra_ent_is_initialized = false;

    memset(vb->vcf_section_bytes, 0, sizeof(vb->vcf_section_bytes));
    memset(vb->z_section_bytes, 0, sizeof(vb->z_section_bytes));
    memset(vb->z_num_sections, 0, sizeof(vb->z_num_sections));
    memset(vb->z_section_entries, 0, sizeof(vb->z_section_entries));

    memset (&vb->profile, 0, sizeof (vb->profile));

    buf_free(&vb->line_variant_data);
    buf_free(&vb->line_gt_data);
    buf_free(&vb->line_ht_data);
    buf_free(&vb->line_phase_data);

    buf_free(&vb->sample_iterator);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->genotype_section_lens_buf);

    buf_free (&vb->format_mapper_buf);
    buf_free (&vb->iname_mapper_buf);

    vb->num_info_subfields=0;

    buf_free(&vb->optimized_gl_dict);
    buf_free(&vb->haplotype_permutation_index);
    buf_free(&vb->haplotype_permutation_index_squeezed);
    
    buf_free(&vb->gt_sb_line_starts_buf);
    buf_free(&vb->gt_sb_line_lengths_buf);
    buf_free(&vb->helper_index_buf);
    buf_free(&vb->compressed);
    buf_free(&vb->vcf_data);
    buf_free(&vb->z_data);
    buf_free(&vb->z_section_headers);
    buf_free(&vb->ht_columns_data);
    buf_free(&vb->spiced_pw);
    buf_free(&vb->format_info_buf);
    buf_free(&vb->ra_buf);
    buf_free(&vb->show_headers_buf);
    buf_free(&vb->show_b250_buf);
    buf_free(&vb->section_list_buf);
    buf_free(&vb->region_ra_intersection_matrix);
    
    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->mtf_ctx[i].dict_id.num)
            mtf_free_context (&vb->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->z_file = vb->vcf_file = NULL; // we're not freeing them, must setting the point to null

    vb->in_use = false; // released the VB back into the pool - it may now be reused

    // backward compatability with genozip v1
    buf_free(&vb->v1_subfields_start_buf);        
    buf_free(&vb->v1_subfields_len_buf);
    buf_free(&vb->v1_num_subfields_buf);
    buf_free(&vb->v1_variant_data_section_data);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    //                         between VBs of a file or concatenated files.
    // vb->column_of_zeros : we don't free this as its a constant array of zeros, of size global_max_lines_per_vb
}

void vb_create_pool (unsigned num_vbs)
{
    ASSERT (!pool || num_vbs==pool->num_vbs, 
            "Error: vb pool already exists, but with the wrong number of vbs - expected %u but it has %u", num_vbs, pool->num_vbs);

    if (!pool)  {
        pool = (VariantBlockPool *)calloc (1, sizeof (VariantBlockPool) + num_vbs * sizeof (VariantBlock)); // note we can't use Buffer yet, because we don't have VBs yet...
        ASSERT0 (pool, "Error: failed to calloc pool");

        pool->num_vbs = num_vbs; 
    }
}

VariantBlockPool *vb_get_pool ()
{
    return pool;
}

void vb_external_vb_initialize()
{
    ASSERT0 (!evb, "Error: evb already initialized");

    evb = calloc (1, sizeof (VariantBlock));
    ASSERT0 (evb, "Error: failed to calloc evb");
    evb->id = -1;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VariantBlock *vb_get_vb (FileP vcf_file, FileP z_file, unsigned variant_block_i)
{
    VariantBlock *vb=NULL;
    
    // see if there's a VB avaiable for recycling
    unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
    
        vb = &pool->vb[vb_i];
    
        if (!vb->in_use) {
            vb->id = vb_i;
            break;
        }
    }

    ASSERT (vb_i < pool->num_vbs, "Error: VB pool is full - it already has %u VBs", pool->num_vbs)

    vb->in_use           = true;
    vb->variant_block_i  = variant_block_i;
    vb->vcf_file         = vcf_file;
    vb->z_file           = z_file;
    
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
void vb_cleanup_memory ()
{
    // see if there's a VB avaiable for recycling
    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VariantBlock *vb = &pool->vb[vb_i];
     
        vb_free_buffer_array (vb, &vb->genotype_sections_data, vb->num_sample_blocks);
        vb_free_buffer_array (vb, &vb->haplotype_sections_data, vb->num_sample_blocks);
        vb_free_buffer_array (vb, &vb->phase_sections_data, vb->num_sample_blocks);

        if (vb->data_lines) {
            for (unsigned i=0; i < global_max_lines_per_vb; i++) {
                DataLine *dl = &vb->data_lines[i];
                
                buf_destroy(vb, &dl->line);
                buf_destroy(vb, &dl->v1_variant_data);
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
