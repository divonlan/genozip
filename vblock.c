// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "libdeflate/libdeflate.h"
#include "genozip.h"
#include "context.h"
#include "vblock.h"
#include "file.h"
#include "reference.h"
#include "digest.h"
#include "bgzf.h"
#include "strings.h"
#include "threads.h"

// pool of VBs allocated based on number of threads
static VBlockPool *pool = NULL;

// one VBlock outside of pool
VBlock *evb = NULL;

#define FINALIZE_VB_BUFS(func, ctx_func, vb_func) \
    func (&vb->lines);               \
    func (&vb->ra_buf);              \
    func (&vb->compressed);          \
    func (&vb->txt_data);            \
    func (&vb->z_data);              \
    func (&vb->z_section_headers);   \
    func (&vb->spiced_pw);           \
    func (&vb->show_headers_buf);    \
    func (&vb->show_b250_buf);       \
    func (&vb->section_list_buf);    \
    func (&vb->bgzf_blocks);         \
    func (&vb->coverage);            \
    func (&vb->read_count);          \
    func (&vb->unmapped_read_count); \
    func (&vb->liftover);            \
    func (&vb->liftover_rejects);    \
    for (unsigned i=0; i < NUM_CODEC_BUFS; i++) func (&vb->codec_bufs[i]); \
    for (unsigned i=0; i < MAX_DICTS; i++) if (vb->contexts[i].dict_id.num) ctx_func (&vb->contexts[i]); \
    if (vb->data_type != DT_NONE) DT_FUNC (vb, vb_func)(vb);    

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb_do (VBlock **vb_p, const char *func) 
{
    VBlock *vb = *vb_p;

    if (!vb) return; // nothing to release

    ASSERT (vb->in_use || vb==evb, "Cannot release VB because it is not in_use (called from %s): vb->id=%u vb->vblock_id=%u", 
            func, vb->id, vb->vblock_i);

    threads_log_by_vb (vb, vb->compute_task ? vb->compute_task : func, "RELEASING VB", 0);

    if (flag.show_vblocks) 
        iprintf ("VB_RELEASE(id=%u) vb_i=%d\n", vb->id, vb->vblock_i);

    if (flag.show_time) 
        profiler_add (vb);

    buf_test_overflows(vb, func); 

    // verify that gzip_compressor was released after use
    ASSERT (!vb->gzip_compressor, "vb=%u: expecting gzip_compressor=NULL", vb->vblock_i);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->data_type   : type of this vb 

    vb->first_line = vb->vblock_i = vb->fragment_len = vb->fragment_num_words = 0;
    vb->vb_data_size = vb->vb_data_size_0 = vb->luft_reject_bytes = vb->longest_line_len = vb->line_i = vb->component_i = vb->grep_stages = 0;
    vb->ready_to_dispatch = vb->is_processed = vb->is_unsorted[0] = vb->is_unsorted[1] = vb->drop_curr_line = false;
    vb->z_next_header_i = 0;
    vb->num_contexts = 0;
    vb->chrom_node_index = vb->chrom_name_len = vb->seq_len = 0; 
    vb->vb_position_txt_file = vb->line_start = 0;
    vb->num_lines_at_1_3 = vb->num_lines_at_2_3 = vb->num_nondrop_lines = 0;
    vb->drop_curr_line = vb->has_non_agct = false;    
    vb->num_type1_subfields = vb->num_type2_subfields = 0;
    vb->range = NULL;
    vb->chrom_name = vb->fragment_start = NULL;
    vb->prev_range = NULL;
    vb->prev_range_chrom_node_index = vb->prev_range_range_i = vb->range_num_set_bits = 0;
    vb->digest_so_far = DIGEST_NONE;
    vb->refhash_layer = vb->refhash_start_in_layer = 0;
    vb->fragment_ctx = vb->ht_matrix_ctx = vb->runs_ctx = vb->fgrc_ctx = NULL;
    vb->fragment_codec = vb->codec_using_codec_bufs = 0;
    vb->ht_per_line = 0;
    vb->is_rejects_vb = 0;
    vb->vb_header_flags = (struct FlagsVbHeader){};
    vb->compute_thread_id = 0;
    vb->compute_task = NULL;
    vb->compute_func = NULL;
    memset(&vb->profile, 0, sizeof (vb->profile));
    memset(vb->dict_id_to_did_i_map, 0, sizeof(vb->dict_id_to_did_i_map));
    mutex_destroy (vb->vb_ready_for_compute_thread);

    FINALIZE_VB_BUFS (buf_free, ctx_free_context, release_vb);

    // this release can be run by either the main or writer thread. we make sure to update in_use as the very
    // last change, and do so atomically

    // case: this VB is from the pool (i.e. not evb)
    if (vb != evb) {
        __atomic_store_n (&vb->in_use, (bool)0, __ATOMIC_RELAXED); // released the VB back into the pool - it may now be reused 
        *vb_p = NULL;
    }
}

void vb_destroy_vb (VBlockP *vb_p)
{
    ASSERTMAINTHREAD;

    VBlockP vb = *vb_p;

    if (!vb) return;

    FINALIZE_VB_BUFS (buf_destroy, ctx_destroy_context, destroy_vb);

    FREE (*vb_p);
}

void vb_create_pool (unsigned num_vbs)
{
    ASSERTMAINTHREAD;

    if (!pool)  {
        // allocation includes array of pointers (initialized to NULL)
        pool = (VBlockPool *)CALLOC (sizeof (VBlockPool) + num_vbs * sizeof (VBlock *)); // note we can't use Buffer yet, because we don't have VBs yet...
        pool->num_vbs = num_vbs; 
    }

    // case: old pool is too small - realloc it (the pool contains only pointers to VBs, so the VBs themselves are not realloced)
    else if (pool->num_vbs < num_vbs) {
        pool = (VBlockPool *)REALLOC (pool, sizeof (VBlockPool) + num_vbs * sizeof (VBlock *), "VBlockPool"); 
        memset (&pool->vb[pool->num_vbs], 0, (num_vbs - pool->num_vbs) * sizeof (VBlockP)); // initialize new entries
        pool->num_vbs = num_vbs; 
    }
}

VBlockPool *vb_get_pool (void)
{
    return pool;
}

VBlockP vb_initialize_nonpool_vb (int vb_id)
{
    VBlockP vb = CALLOC (sizeof (VBlock));
    vb->data_type = DT_NONE;
    vb->id = vb_id;
    vb->compute_task = "evb";
    return vb;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlock *vb_get_vb (const char *task_name, uint32_t vblock_i)
{
    ASSERTMAINTHREAD;

    // circle around until a VB becomes available (busy wait)
    unsigned vb_id; for (vb_id=0; ; vb_id = (vb_id+1) % pool->num_vbs) {
    
        // free if this is a VB allocated by a previous file, with a different data type
        // note: if z_file is DT_NONE, we're performing a generic fan_out task, and  data/GRCh38_full_analysis_set_plus_decoy_hla.ref.genozipwe can use what ever VB already exists
        if (pool->vb[vb_id] && z_file && pool->vb[vb_id]->data_type != z_file->data_type) {
            vb_destroy_vb (&pool->vb[vb_id]);
            pool->num_allocated_vbs--; // we will immediately allocate and increase this back
        }
        
        if (!pool->vb[vb_id]) { // VB is not allocated - allocate it
            unsigned sizeof_vb = command==ZIP ? (txt_file && DTPT(sizeof_vb) ? DTPT(sizeof_vb)() : sizeof (VBlock))
                                              : (z_file   && DTPZ(sizeof_vb) ? DTPZ(sizeof_vb)() : sizeof (VBlock));
            pool->vb[vb_id] = CALLOC (sizeof_vb); 
            pool->num_allocated_vbs++;
            pool->vb[vb_id]->data_type = command==ZIP ? (txt_file ? txt_file->data_type : DT_NONE)
                                                     : (z_file   ? z_file->data_type   : DT_NONE);
        }

        bool in_use = __atomic_load_n (&pool->vb[vb_id]->in_use, __ATOMIC_RELAXED);
        if (!in_use) break;

        // case: we've checked all the VBs and none is available - wait a bit and check again
        if (vb_id == pool->num_vbs-1) {
            // this happens when a lot VBs are handed over to the writer thread which has not processed them yet.
            // for example, if writer is blocking on write(), waiting for a pipe counterpart to read.
            // it will be released when the writer thread completes one VB.
            usleep (1000); // 1 ms
        }
    }

    // initialize VB fields that need to be a value other than 0
    VBlock *vb = pool->vb[vb_id];
    vb->id                = vb_id;
    vb->in_use            = true;
    vb->vblock_i          = vblock_i;
    vb->buffer_list.vb    = vb;
    vb->compute_thread_id = THREAD_ID_NONE;
    vb->compute_task      = task_name;
    memset (vb->dict_id_to_did_i_map, 0xff, sizeof(vb->dict_id_to_did_i_map)); // DID_I_NONE

    if (flag.show_vblocks) 
        iprintf ("VB_GET_VB(task=%s id=%u) vb_i=%d\n", task_name, vb->id, vb->vblock_i);

    threads_log_by_vb (vb, task_name, "GET VB", 0);

    return vb;
}

bool vb_has_free_vb (void)
{
    for (unsigned vb_id=0; vb_id < pool->num_vbs ; vb_id++) 
        if (!pool->vb[vb_id] || ! __atomic_load_n (&pool->vb[vb_id]->in_use, __ATOMIC_RELAXED)) 
            return true;

    return false;
}

// free memory allocations between files, when compressing multiple non-bound files or decompressing multiple files
void vb_cleanup_memory (void)
{
    ASSERTMAINTHREAD;

    if (!pool) return;

    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VBlock *vb = pool->vb[vb_i];
        if (vb && 
            vb->data_type == z_file->data_type && // skip VBs that were initialized by a previous file of a different data type and not used by this file
            DTPZ(cleanup_memory)) 
            DTPZ(cleanup_memory)(vb);
    }

    if (z_file->data_type != DT_NONE && DTPZ(cleanup_memory))
        DTPZ(cleanup_memory)(evb);

    ref_unload_reference();
}

// frees memory of all VBs, except for non-pool VBs (evb)
void vb_destroy_pool_vbs (void)
{
    ASSERTMAINTHREAD;

    if (!pool) return;

    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) 
        vb_destroy_vb (&pool->vb[vb_i]);

    FREE (pool);
}

// NOT thread safe, use only in execution-terminating messages
const char *err_vb_pos (void *vb)
{
    static char s[80];
    sprintf (s, "vb i=%u position in %s file=%"PRIu64, 
             ((VBlockP)vb)->vblock_i, dt_name (txt_file->data_type), ((VBlockP)vb)->vb_position_txt_file);
    return s;
}

