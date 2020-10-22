// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "genozip.h"
#include "context.h"
#include "vblock.h"
#include "file.h"
#include "reference.h"
#include "md5.h"

// pool of VBs allocated based on number of threads
static VBlockPool *pool = NULL;

// one VBlock outside of pool
VBlock *evb = NULL;

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb (VBlock *vb) 
{
    if (!vb) return; // nothing to release

    vb->first_line = vb->vblock_i = vb->txt_data_next_offset = vb->num_haplotypes_per_line = 0;
    vb->vb_data_size = vb->vb_data_read_size = vb->longest_line_len = vb->line_i = vb->grep_stages = 0;
    vb->ready_to_dispatch = vb->is_processed = vb->dont_show_curr_line = false;
    vb->z_next_header_i = 0;
    vb->num_contexts = 0;
    vb->chrom_node_index = vb->chrom_name_len = vb->seq_len = 0; 
    vb->vb_position_txt_file = vb->line_start = 0;
    vb->num_lines_at_1_3 = vb->num_lines_at_2_3 = 0;
    vb->dont_show_curr_line = vb->has_non_agct = false;    
    vb->num_type1_subfields = vb->num_type2_subfields = 0;
    vb->range = NULL;
    vb->chrom_name = NULL;
    vb->prev_range = NULL;
    vb->ht_matrix_ctx = NULL;
    vb->gtshark_gt = vb->gtshark_db = vb->gtshark_x_line = vb->gtshark_x_ht = vb->gtshark_x_allele = NULL;
    vb->prev_range_chrom_node_index = vb->prev_range_range_i = vb->range_num_set_bits = 0;
    vb->md5_hash_so_far = MD5HASH_NONE;
    vb->refhash_layer = vb->refhash_start_in_layer = 0;
    
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
    buf_free(&vb->hapmat_helper_index_buf);
    buf_free(&vb->hapmat_columns_data);
    buf_free(&vb->hapmat_one_array);
    buf_free(&vb->hapmat_column_of_zeros);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->contexts[i].dict_id.num)
            mtf_free_context (&vb->contexts[i]);

    for (unsigned i=0; i < NUM_CODEC_BUFS; i++)
        buf_free (&vb->codec_bufs[i]);
        
    vb->in_use = false; // released the VB back into the pool - it may now be reused

    // release data_type -specific fields
    if (vb->data_type != DT_NONE && DTP(release_vb)) 
        DTP(release_vb)(vb);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->num_lines_alloced
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    //                         between VBs of a file or bound files.
    // vb->data_type : type of this vb 
}

void vb_destroy_vb (VBlockP *vb_p)
{
    VBlockP vb = *vb_p;

    buf_destroy (&vb->ra_buf);
    buf_destroy (&vb->compressed);
    buf_destroy (&vb->txt_data);
    buf_destroy (&vb->txt_data_spillover);
    buf_destroy (&vb->z_data);
    buf_destroy (&vb->z_section_headers);
    buf_destroy (&vb->spiced_pw);
    buf_destroy (&vb->show_headers_buf);
    buf_destroy (&vb->show_b250_buf);
    buf_destroy (&vb->section_list_buf);
    buf_destroy (&vb->region_ra_intersection_matrix);
    buf_destroy (&vb->hapmat_helper_index_buf);
    buf_destroy (&vb->hapmat_columns_data);
    buf_destroy (&vb->hapmat_one_array);
    buf_destroy (&vb->hapmat_column_of_zeros);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->contexts[i].dict_id.num)
            mtf_destroy_context (&vb->contexts[i]);

    for (unsigned i=0; i < NUM_CODEC_BUFS; i++)
        buf_destroy (&vb->codec_bufs[i]);

    // destory data_type -specific buffers
    if (vb->data_type != DT_NONE && DTP(destroy_vb))
        DTP(destroy_vb)(vb);

    FREE (*vb_p);
}

void vb_create_pool (unsigned num_vbs)
{
    ASSERT (!pool || num_vbs==pool->num_vbs, 
            "Error: vb pool already exists, but with the wrong number of vbs - expected %u but it has %u", num_vbs, pool->num_vbs);

    if (!pool)  {
        // allocation includes array of pointers (initialized to NULL)
        pool = (VBlockPool *)CALLOC (sizeof (VBlockPool) + num_vbs * sizeof (VBlock *)); // note we can't use Buffer yet, because we don't have VBs yet...
        pool->num_vbs = num_vbs; 
    }
}

VBlockPool *vb_get_pool (void)
{
    return pool;
}

void vb_initialize_evb(void)
{
    ASSERT0 (!evb, "Error: evb already initialized");

    evb = CALLOC (sizeof (VBlock));
    evb->data_type = DT_NONE;
    evb->id = -1;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlock *vb_get_vb (unsigned vblock_i)
{
    // see if there's a VB avaiable for recycling
    unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
    
        // free if this is a VB allocated by a previous file, with a different data type
        // note: if z_file is DT_NONE, we're performing a generic fan_out task, and  data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozipwe can use what ever VB already exists
        if (pool->vb[vb_i] && z_file && pool->vb[vb_i]->data_type != z_file->data_type) {
            vb_destroy_vb (&pool->vb[vb_i]);
            pool->num_allocated_vbs--; // we will immediately allocate and increase this back
        }
        
        if (!pool->vb[vb_i]) { // VB is not allocated - allocate it
            unsigned sizeof_vb = z_file && DTPZ(sizeof_vb) ? DTPZ(sizeof_vb)() : sizeof (VBlock);
            pool->vb[vb_i] = CALLOC (sizeof_vb); 
            pool->num_allocated_vbs++;
            pool->vb[vb_i]->data_type = z_file ? z_file->data_type : DT_NONE;
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

    return vb;
}

// free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files
void vb_cleanup_memory (void)
{
    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VBlock *vb = pool->vb[vb_i];
        if (vb && 
            vb->data_type == z_file->data_type && // skip VBs that were initialized by a previous file of a different data type and not used by this file
            DTPZ(cleanup_memory)) 
            DTPZ(cleanup_memory)(vb);
    }

    if (z_file->data_type != DT_NONE && DTPZ(cleanup_memory))
        DTPZ(cleanup_memory)(evb);

    ref_unload_reference (false);
}

// NOT thread safe, use only in execution-terminating messages
const char *err_vb_pos (void *vb)
{
    static char s[80];
    sprintf (s, "vb i=%u position in %s file=%"PRIu64, 
             ((VBlockP)vb)->vblock_i, dt_name (txt_file->data_type), ((VBlockP)vb)->vb_position_txt_file);
    return s;
}

void vb_set_global_max_memory_per_vb (const char *mem_size_mb_str)
{
    const char *err_msg = "Error: invalid argument of --vblock: %s. Expecting an integer between 1 and 2048. The file will be read and processed in blocks of this number of megabytes.";

    unsigned len = strlen (mem_size_mb_str);
    ASSERT (len <= 4 || (len==1 && mem_size_mb_str[0]=='0'), err_msg, mem_size_mb_str);
    ASSERT (strspn (mem_size_mb_str, "0123456789") == len, err_msg, mem_size_mb_str);

    unsigned mem_size_mb = atoi (mem_size_mb_str);
    ASSERT (mem_size_mb <= 2048, err_msg, mem_size_mb_str);

    global_max_memory_per_vb = mem_size_mb * 1024 * 1024;
}
