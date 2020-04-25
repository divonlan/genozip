// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "genozip.h"
#include "move_to_front.h"
#include "vblock.h"
#include "file.h"

// pool of VBs allocated based on number of threads
static VBlockPool *pool = NULL;

// one VBlock outside of pool
VBlock *evb = NULL;

//--------------------------------
// VCF stuff
//--------------------------------

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
static void vb_vcf_release_vb (VBlock *vb_) 
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }

    vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->max_gt_line_len = vb->next_numeric_id = 0;

    buf_free(&vb->line_gt_data);
    buf_free(&vb->line_ht_data);
    buf_free(&vb->line_phase_data);

    buf_free(&vb->sample_iterator);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->is_sb_included);
    buf_free(&vb->genotype_section_lens_buf);
    buf_free(&vb->id_numeric_data);

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

static void vb_vcf_destroy_vb (VBlock **vb_)
{
    VBlockVCF **vb = (VBlockVCF **)vb_;

    for (unsigned i=0; i < (*vb)->num_sample_blocks; i++) {
        if ((*vb)->haplotype_sections_data) buf_destroy (&(*vb)->haplotype_sections_data[i]);
        if ((*vb)->genotype_sections_data)  buf_destroy (&(*vb)->genotype_sections_data[i]);
        if ((*vb)->phase_sections_data)     buf_destroy (&(*vb)->phase_sections_data[i]);
    }

    buf_destroy (&(*vb)->line_gt_data);
    buf_destroy (&(*vb)->line_ht_data);
    buf_destroy (&(*vb)->line_phase_data);
    buf_destroy (&(*vb)->sample_iterator);
    buf_destroy (&(*vb)->genotype_one_section_data);
    buf_destroy (&(*vb)->is_sb_included);
    buf_destroy (&(*vb)->genotype_section_lens_buf);
    buf_destroy (&(*vb)->format_mapper_buf);
    buf_destroy (&(*vb)->iname_mapper_buf);
    buf_destroy (&(*vb)->optimized_gl_dict);
    buf_destroy (&(*vb)->haplotype_permutation_index);
    buf_destroy (&(*vb)->haplotype_permutation_index_squeezed);
    buf_destroy (&(*vb)->gt_sb_line_starts_buf);
    buf_destroy (&(*vb)->gt_sb_line_lengths_buf);
    buf_destroy (&(*vb)->helper_index_buf);
    buf_destroy (&(*vb)->ht_columns_data);
    buf_destroy (&(*vb)->format_info_buf);
    buf_destroy (&(*vb)->id_numeric_data);
    buf_destroy (&(*vb)->column_of_zeros);
    buf_destroy (&(*vb)->gtshark_db_db_data);
    buf_destroy (&(*vb)->gtshark_db_gt_data);
    buf_destroy (&(*vb)->gtshark_exceptions_line_i);
    buf_destroy (&(*vb)->gtshark_exceptions_ht_i);
    buf_destroy (&(*vb)->gtshark_exceptions_allele);
    buf_destroy (&(*vb)->gtshark_vcf_data);
    buf_destroy (&(*vb)->v1_subfields_start_buf);        
    buf_destroy (&(*vb)->v1_subfields_len_buf);
    buf_destroy (&(*vb)->v1_num_subfields_buf);
    buf_destroy (&(*vb)->v1_variant_data_section_data);
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
static void vb_vcf_cleanup_memory (VBlock *vb_)
{
    VBlockVCF *vb = (VBlockVCF *)vb_;

    vb_vcf_free_buffer_array (&vb->genotype_sections_data, vb->num_sample_blocks);
    vb_vcf_free_buffer_array (&vb->haplotype_sections_data, vb->num_sample_blocks);
    vb_vcf_free_buffer_array (&vb->phase_sections_data, vb->num_sample_blocks);
    vb->num_sample_blocks = 0;
}

//--------------------------------
// SAM stuff
//--------------------------------

static void vb_sam_initialize_vb (VBlock *vb_)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;

    vb->last_rname_node_index = (uint32_t)-1;
}

static void vb_sam_release_vb (VBlock *vb_)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;

    vb->num_optional_subfield_b250s = vb->next_seq = vb->next_qual = vb->next_random_pos = vb->next_md = 0;
    vb->nm_did_i = vb->strand_did_i = vb->last_tlen_abs_len = 0;
    vb->last_rname_node_index = 0;
    vb->last_tlen_abs = 0;
    vb->last_tlen_is_positive = 0;

    memset (&vb->qname_mapper, 0, sizeof (vb->qname_mapper));
    
    buf_free (&vb->random_pos_data);
    buf_free (&vb->optional_mapper_buf);
    buf_free (&vb->seq_data);
    buf_free (&vb->qual_data);
    buf_free (&vb->md_data);
}

// free all memory of a VB
static void vb_sam_destroy_vb (VBlock **vb_)
{
    VBlockSAM **vb = (VBlockSAM **)vb_;

    buf_destroy (&(*vb)->random_pos_data);
    buf_destroy (&(*vb)->optional_mapper_buf);
    buf_destroy (&(*vb)->seq_data);
    buf_destroy (&(*vb)->qual_data);    
    buf_destroy (&(*vb)->md_data);    
}

//--------------------------------
// FASTQ stuff
//--------------------------------

static void vb_fast_release_vb (VBlock *vb_)
{
    VBlockFAST *vb = (VBlockFAST *)vb_;

    vb->next_seq = vb->next_qual = vb->next_comment = vb->last_line = 0;

    memset (&vb->desc_mapper, 0, sizeof (vb->desc_mapper));
    
    buf_free (&vb->seq_data);
    buf_free (&vb->qual_data);
    buf_free (&vb->comment_data);
}

// free all memory of a VB
static void vb_fast_destroy_vb (VBlock **vb_)
{
    VBlockFAST **vb = (VBlockFAST **)vb_;

    buf_destroy (&(*vb)->seq_data);
    buf_destroy (&(*vb)->qual_data);    
    buf_destroy (&(*vb)->comment_data);
}

//--------------------------------
// ME23 stuff
//--------------------------------

static void vb_me23_release_vb (VBlock *vb_)
{
    VBlockME23 *vb = (VBlockME23 *)vb_;

    vb->next_numeric_id = vb->next_genotype = 0;

    buf_free (&vb->id_numeric_data);
    buf_free (&vb->genotype_data);
}

// free all memory of a VB
static void vb_me23_destroy_vb (VBlock **vb_)
{
    VBlockME23 **vb = (VBlockME23 **)vb_;

    buf_destroy (&(*vb)->id_numeric_data);
    buf_destroy (&(*vb)->genotype_data);
}

//--------------------------------
// Non-type-specific stuff
//--------------------------------

static void donothing() {}

// data type multiplexors
typedef void (*VbFunc) (VBlock *vb);
typedef void (*VbFuncP) (VBlock **vb);

static size_t vb_size_by_dt[NUM_DATATYPES] = { sizeof (VBlockVCF), sizeof (VBlockSAM), 
                                               sizeof (VBlockFAST), sizeof (VBlockFAST), sizeof (VBlockME23) };

static VbFunc vb_release_by_dt[NUM_DATATYPES] = { vb_vcf_release_vb, vb_sam_release_vb, 
                                                  vb_fast_release_vb, vb_fast_release_vb, vb_me23_release_vb };

static VbFuncP vb_destroy_by_dt[NUM_DATATYPES] = { vb_vcf_destroy_vb, vb_sam_destroy_vb, 
                                                   vb_fast_destroy_vb, vb_fast_destroy_vb, vb_me23_destroy_vb };

static VbFunc vb_init_by_dt[NUM_DATATYPES] = { donothing, vb_sam_initialize_vb, donothing, donothing, donothing };

static VbFunc vb_cleanup_by_dt[NUM_DATATYPES] = { vb_vcf_cleanup_memory, donothing, donothing, donothing, donothing };

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb (VBlock *vb) 
{
    if (!vb) return; // nothing to release

    vb->first_line = vb->vblock_i = vb->txt_data_next_offset = 0;
    vb->vb_data_size = vb->vb_data_read_size = vb->last_pos = vb->longest_line_len = vb->line_i = 0;
    vb->ready_to_dispatch = vb->is_processed = false;
    vb->z_next_header_i = 0;
    vb->num_dict_ids = 0;
    vb->chrom_node_index = 0; 
    vb->vb_position_txt_file = 0;

    memset(vb->txt_section_bytes, 0, sizeof(vb->txt_section_bytes));
    memset(vb->z_section_bytes, 0, sizeof(vb->z_section_bytes));
    memset(vb->z_num_sections, 0, sizeof(vb->z_num_sections));
    memset(vb->z_section_entries, 0, sizeof(vb->z_section_entries));
    memset(&vb->profile, 0, sizeof (vb->profile));
    memset(vb->dict_id_to_did_i_map, 0, sizeof(vb->dict_id_to_did_i_map));

    buf_free(&vb->lines);
    buf_free(&vb->ra_buf);
    buf_free(&vb->compressed);
    buf_free(&vb->txt_data);
    buf_free(&vb->txt_data_spillover);
    buf_free(&vb->z_data);
    buf_free(&vb->z_section_headers);
    buf_free(&vb->spiced_pw);
    buf_free(&vb->show_headers_buf);
    buf_free(&vb->show_b250_buf);
    buf_free(&vb->section_list_buf);
    buf_free(&vb->region_ra_intersection_matrix);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->mtf_ctx[i].dict_id.num)
            mtf_free_context (&vb->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->in_use = false; // released the VB back into the pool - it may now be reused

    // release data_type -specific fields
    if (vb->data_type != DT_NONE)
        vb_release_by_dt[vb->data_type](vb);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->num_lines_alloced
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    //                         between VBs of a file or concatenated files.
    // vb->data_type : type of this vb 
}

void vb_destroy_vb (VBlockP *vb)
{
    buf_destroy (&(*vb)->ra_buf);
    buf_destroy (&(*vb)->compressed);
    buf_destroy (&(*vb)->txt_data);
    buf_destroy (&(*vb)->txt_data_spillover);
    buf_destroy (&(*vb)->z_data);
    buf_destroy (&(*vb)->z_section_headers);
    buf_destroy (&(*vb)->spiced_pw);
    buf_destroy (&(*vb)->show_headers_buf);
    buf_destroy (&(*vb)->show_b250_buf);
    buf_destroy (&(*vb)->section_list_buf);
    buf_destroy (&(*vb)->region_ra_intersection_matrix);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if ((*vb)->mtf_ctx[i].dict_id.num)
            mtf_destroy_context (&(*vb)->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_destroy (&(*vb)->compress_bufs[i]);

    // destory data_type -specific buffers
    if ((*vb)->data_type != DT_NONE)
        vb_destroy_by_dt[(*vb)->data_type](vb);

    FREE (*vb);
    *vb = NULL;
}

void vb_create_pool (unsigned num_vbs)
{
    ASSERT (!pool || num_vbs==pool->num_vbs, 
            "Error: vb pool already exists, but with the wrong number of vbs - expected %u but it has %u", num_vbs, pool->num_vbs);

    if (!pool)  {
        // allocation includes array of pointers (initialized to NULL)
        pool = (VBlockPool *)calloc (1, sizeof (VBlockPool) + num_vbs * sizeof (VBlock *)); // note we can't use Buffer yet, because we don't have VBs yet...
        ASSERT0 (pool, "Error: failed to calloc pool");

        pool->num_vbs = num_vbs; 
    }
}

VBlockPool *vb_get_pool (void)
{
    return pool;
}

void vb_external_vb_initialize(void)
{
    ASSERT0 (!evb, "Error: evb already initialized");

    evb = calloc (1, sizeof (VBlock));
    ASSERT0 (evb, "Error: failed to calloc evb");
    evb->data_type = DT_NONE;
    evb->id = -1;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlock *vb_get_vb (unsigned vblock_i)
{
    // see if there's a VB avaiable for recycling
    unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
    
        // free if this is a VB allocated by a previous file, with a different data type
        if (pool->vb[vb_i] && pool->vb[vb_i]->data_type != z_file->data_type) {
            vb_destroy_vb (&pool->vb[vb_i]);
            pool->num_allocated_vbs--; // we will immediately allocate and increase this back
        }
        
        if (!pool->vb[vb_i]) { // VB is not allocated - allocate it
            pool->vb[vb_i] = calloc (vb_size_by_dt[z_file->data_type], 1); 
            pool->num_allocated_vbs++;
            pool->vb[vb_i]->data_type = z_file->data_type;
        }

        if (!pool->vb[vb_i]->in_use) break;
    }

    ASSERT (vb_i < pool->num_vbs, "Error: VB pool is full - it already has %u VBs", pool->num_vbs)

    // initialize VB fields that need to be a value other than 0
    VBlock *vb = pool->vb[vb_i];
    vb->id             = vb_i;
    vb->in_use         = true;
    vb->vblock_i       = vblock_i;
    vb->buffer_list.vb = vb;
    memset (vb->dict_id_to_did_i_map, DID_I_NONE, sizeof(vb->dict_id_to_did_i_map));

    // initialize data-type-specific fields
    if (vb->data_type != DT_NONE)    
        vb_init_by_dt[z_file->data_type](vb);

    return vb;
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vb_cleanup_memory (void)
{
    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VBlock *vb = pool->vb[vb_i];
        if (vb && vb->data_type != DT_NONE) vb_cleanup_by_dt[z_file->data_type](vb);
    }

    global_vcf_num_samples = 0;
}

unsigned vb_vcf_num_samples_in_sb (const VBlockVCF *vb, unsigned sb_i)
{
    // case: last block has less than a full block of samples
    if (sb_i == vb->num_sample_blocks-1 && global_vcf_num_samples % vb->num_samples_per_block)
        return global_vcf_num_samples % vb->num_samples_per_block;

    else
        return vb->num_samples_per_block;
} 

unsigned vb_vcf_num_sections(VBlockVCF *vb) 
{
    return 1 + vb->has_genotype_data + (vb->phase_type == PHASE_MIXED_PHASED) + (vb->num_haplotypes_per_line > 0);
}

// NOT thread safe, use only in execution-terminating messages
const char *err_vb_pos (void *vb)
{
    static char s[80];
    sprintf (s, "vb i=%u position in %s file=%"PRIu64, 
             ((VBlockP)vb)->vblock_i, dt_name (txt_file->data_type), ((VBlockP)vb)->vb_position_txt_file);
    return s;
}
