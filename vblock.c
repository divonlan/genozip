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

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb (VBlock **vb_p) 
{
    if (! *vb_p) return; // nothing to release

    VBlock *vb = *vb_p;
    *vb_p = NULL;

    vb->num_lines = vb->first_line = vb->vblock_i = vb->txt_data_next_offset = 0;
    vb->vb_data_size = vb->vb_data_read_size = 0;
    vb->ready_to_dispatch = vb->is_processed = false;
    vb->z_next_header_i = 0;
    vb->num_dict_ids = 0;

    memset(vb->txt_section_bytes, 0, sizeof(vb->txt_section_bytes));
    memset(vb->z_section_bytes, 0, sizeof(vb->z_section_bytes));
    memset(vb->z_num_sections, 0, sizeof(vb->z_num_sections));
    memset(vb->z_section_entries, 0, sizeof(vb->z_section_entries));

    memset (&vb->profile, 0, sizeof (vb->profile));

    buf_free(&vb->compressed);
    buf_free(&vb->txt_data);
    buf_free(&vb->txt_data_spillover);
    buf_free(&vb->z_data);
    buf_free(&vb->z_section_headers);
    buf_free(&vb->spiced_pw);
    buf_free(&vb->show_headers_buf);
    buf_free(&vb->show_b250_buf);
    buf_free(&vb->section_list_buf);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->mtf_ctx[i].dict_id.num)
            mtf_free_context (&vb->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->in_use = false; // released the VB back into the pool - it may now be reused

    if      (vb->data_type == DATA_TYPE_VCF) vb_vcf_release_vb ((VBlockVCF*)vb);
    else if (vb->data_type == DATA_TYPE_SAM) vb_sam_release_vb ((VBlockSAM*)vb);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->num_data_lines_allocated
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    //                         between VBs of a file or concatenated files.
    // vb->data_type : type of this vb 
}

void vb_create_pool (unsigned num_vbs)
{
    ASSERT (!pool || num_vbs==pool->num_vbs, 
            "Error: vb pool already exists, but with the wrong number of vbs - expected %u but it has %u", num_vbs, pool->num_vbs);

    if (!pool)  {
        // allocation includes array of pointers (initialized to NULL)
        pool = (VBlockPool *)calloc (1, sizeof (VBlockPool) + num_vbs * sizeof (VBlockVCF *)); // note we can't use Buffer yet, because we don't have VBs yet...
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
    evb->data_type = DATA_TYPE_NONE;
    evb->id = -1;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlock *vb_get_vb (unsigned vblock_i)
{
    // see if there's a VB avaiable for recycling
    unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
    
        // free if this is a VB allocated by a previous file, with a different data type
        if (pool->vb[vb_i] && pool->vb[vb_i]->data_type != z_file->data_type) {
            FREE (pool->vb[vb_i]);
            pool->vb[vb_i] = NULL; 
            pool->num_allocated_vbs--; // we will immediately allocate and increase this back
        }
        
        if (!pool->vb[vb_i]) { // VB is not allocated - allocate it
            switch (z_file->data_type) {
                case DATA_TYPE_VCF : pool->vb[vb_i] = calloc (sizeof (VBlockVCF), 1); break;
                case DATA_TYPE_SAM : pool->vb[vb_i] = calloc (sizeof (VBlockSAM), 1); break;
                default            : ABORT ("Error in vb_get_vb: Invalid data_type=%d", z_file->data_type);
            }
            pool->num_allocated_vbs++;
            pool->vb[vb_i]->data_type = z_file->data_type;
        }

        if (!pool->vb[vb_i]->in_use) {
            pool->vb[vb_i]->id = vb_i;
            break;
        }
    }

    ASSERT (vb_i < pool->num_vbs, "Error: VB pool is full - it already has %u VBs", pool->num_vbs)

    pool->vb[vb_i]->in_use           = true;
    pool->vb[vb_i]->vblock_i         = vblock_i;
    pool->vb[vb_i]->buffer_list.vb_i = vblock_i;

    return pool->vb[vb_i];
}

// free memory allocations that assume subsequent files will have the same number of samples.
void vb_cleanup_memory (void)
{
    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VBlock *vb = pool->vb[vb_i];
        if      (vb && vb->data_type == DATA_TYPE_VCF) vb_vcf_cleanup_memory ((VBlockVCF *)vb);
        else if (vb && vb->data_type == DATA_TYPE_SAM) vb_sam_cleanup_memory ((VBlockSAM *)vb);
    }

    global_vcf_num_samples = 0;
}
