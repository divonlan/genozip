// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

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
    func (&vb->is_dropped);          \
    func (&vb->ra_buf[0]);           \
    func (&vb->ra_buf[1]);           \
    func (&vb->compressed);          \
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
    func (&vb->chrom2ref_map);       \
    func (&vb->ol_chrom2ref_map);    \
    func (&vb->lo_rejects[0]);       \
    func (&vb->lo_rejects[1]);       \
    func (&vb->frozen_state);        \
    for (unsigned i=0; i < NUM_CODEC_BUFS; i++) func (&vb->codec_bufs[i]); \
    if (vb->data_type != DT_NONE)    \
        for (unsigned i=0; i < MAX_DICTS; i++) if (CTX(i)->dict_id.num || i < DTF(num_fields)) ctx_func (CTX(i), i); /* note: always erase num_fields as they may be set in *_seg_initialize even if not used */\
    func (&vb->txt_data); /* handle contexts first, as vcf_seg_FORMAT_GT overlays FORMAT_GT_HT on txt_data */ \
    if (vb->data_type_alloced != DT_NONE && dt_props[vb->data_type_alloced].vb_func) dt_props[vb->data_type_alloced].vb_func(vb);    

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb_do (VBlockP *vb_p, const char *func) 
{
    START_TIMER;

    VBlock *vb = *vb_p;

    if (!vb) return; // nothing to release

    ASSERT (vb->in_use || vb==evb, "Cannot release VB because it is not in_use (called from %s): vb->id=%d vb->vblock_id=%u", 
            func, vb->id, vb->vblock_i);

    threads_log_by_vb (vb, vb->compute_task ? vb->compute_task : func, "RELEASING VB", 0);

    if (flag.show_vblocks) 
        iprintf ("VB_RELEASE(id=%d) vb_i=%d caller=%s\n", vb->id, vb->vblock_i, func);

    if (flag.show_time) 
        profiler_add (vb);

    if (vb->id != VB_ID_EVB) // cannot test evb, see comment in buf_test_overflows_do
        buf_test_overflows(vb, func); 

    // verify that gzip_compressor was released after use
    ASSERT (!vb->gzip_compressor, "vb=%u: expecting gzip_compressor=NULL", vb->vblock_i);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->data_type   : type of this vb (maybe changed by vb_get_vb)
    // vb->data_type_alloced (maybe changed by vb_get_vb)
    // vb->id

    vb->first_line = vb->vblock_i = vb->fragment_len = vb->fragment_num_words = vb->pos_aln_i = 0;
    vb->recon_size = vb->recon_size_luft = vb->txt_size = vb->txt_size_source_comp = vb->reject_bytes = vb->longest_line_len = vb->line_i = vb->component_i = vb->grep_stages = vb->sample_i = 0;
    vb->recon_num_lines = vb->recon_num_lines_luft = 0;
    vb->ready_to_dispatch = vb->is_processed = vb->is_unsorted[0] = vb->is_unsorted[1] = false;
    vb->z_next_header_i = 0;
    vb->num_contexts = 0;
    vb->line_coords = DC_NONE;
    vb->chrom_node_index = vb->chrom_name_len = vb->seq_len = 0; 
    vb->vb_position_txt_file = vb->line_start = 0;
    vb->num_lines_at_1_3 = vb->num_lines_at_2_3 = vb->num_nondrop_lines = 0;
    vb->is_rejects_vb = false;    
    vb->num_type1_subfields = vb->num_type2_subfields = 0;
    vb->range = NULL;
    vb->drop_curr_line = vb->chrom_name = vb->fragment_start = NULL;
    vb->prev_range[0] = vb->prev_range[1] = NULL;
    vb->prev_range_chrom_node_index[0] = vb->prev_range_chrom_node_index[1] = vb->prev_range_range_i = 0;
    vb->digest_so_far = DIGEST_NONE;
    vb->refhash_layer = vb->refhash_start_in_layer = 0;
    vb->fragment_ctx = vb->ht_matrix_ctx = vb->runs_ctx = vb->fgrc_ctx = NULL;
    vb->fragment_codec = vb->codec_using_codec_bufs = 0;
    vb->ht_per_line = 0;
    vb->vb_coords = DC_NONE;
    vb->compute_thread_id = 0;
    vb->compute_task = NULL;
    vb->compute_func = NULL;
    vb->translation = (DtTranslation){};
    vb->ref = NULL;
    vb->buddy_line_i = 0;

    memset(&vb->profile, 0, sizeof (vb->profile));
    memset(vb->dict_id_to_did_i_map, 0, sizeof(vb->dict_id_to_did_i_map));
    vb->iupacs_last_range[0] = vb->iupacs_last_range[1] = NULL;
    vb->iupacs_last_pos[0] = vb->iupacs_last_pos[1] = vb->iupacs_next_pos[0] = vb->iupacs_next_pos[1] = 0;
    vb->num_rollback_ctxs=0;
    memset (vb->rollback_ctxs, 0, sizeof(vb->rollback_ctxs));
    mutex_destroy (vb->vb_ready_for_compute_thread);

    FINALIZE_VB_BUFS (buf_free, ctx_free_context, release_vb);

    // this release can be run by either the main or writer thread. we make sure to update in_use as the very
    // last change, and do so atomically

    // case: this VB is from the pool (i.e. not evb)
    if (vb != evb) {
        __atomic_store_n (&vb->in_use, (bool)0, __ATOMIC_RELAXED); // released the VB back into the pool - it may now be reused 
        *vb_p = NULL;
    }

    if (vb->id >= 0) COPY_TIMER_VB (evb, vb_release_vb_do)
}

void vb_destroy_vb_do (VBlockP *vb_p, const char *func)
{
    ASSERTMAINTHREAD;
    START_TIMER;

    VBlockP vb = *vb_p;
    if (!vb) return;

    if (flag.show_vblocks) 
        iprintf ("VB_DESTROY(id=%d) vb_i=%d caller=%s\n", vb->id, vb->vblock_i, func);

    bool is_evb = vb->id == VB_ID_EVB;

    FINALIZE_VB_BUFS (buf_destroy, ctx_destroy_context, destroy_vb);

    // test that vb_release_vb_do indeed frees everything
    if (flag.debug && vb->in_use) {
        unsigned sizeof_vb = DT_FUNC(vb, sizeof_vb)(vb->data_type_alloced);
        DataType dt = vb->data_type_alloced;

        vb_release_vb_do (&vb, "vb_destroy_vb"); 
        vb = *vb_p; // re-initialize after release zeroed it

        // clear the few fields purposely not cleared by vb_release_vb_do
        buf_destroy (&(*vb_p)->buffer_list); // used by vb_release_vb_do
        vb->data_type = vb->data_type_alloced = 0; 
        vb->id = 0;
        
        for (char *c=(char*)vb; c  < (char*)(vb) + sizeof_vb; c++)
            if (*c) {
                #define REL_LOC(field) (((char*)(&((VBlock *)vb)->field)) - (char*)vb) // <-- to use private datatype VBlocks, temporarily include the private .h for debugging
                fprintf (stderr, "sizeof_vb=%u sizeof(VBlock)=%u. Bad byte=%u Your field: %u\n", 
                         sizeof_vb, (int)sizeof (VBlock), (unsigned)(c - (char*)vb), (int)REL_LOC(contexts[20])); // <-- to find the offending field, plug in field names and run iteratively
                
                ABORT ("vb_release_vb_do of %s didn't fully clear the VB, byte %u != 0", dt_name(dt), (unsigned)(c - (char*)vb));
            }
    }
    else {
        // clear the few fields purposely not cleared by vb_release_vb_do
        buf_destroy (&(*vb_p)->buffer_list);
        vb->data_type = vb->data_type_alloced = 0; 
        vb->id = 0;
    }

    FREE (*vb_p);

    if (!is_evb) COPY_TIMER_VB (evb, vb_destroy_vb)
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
        REALLOC (&pool, sizeof (VBlockPool) + num_vbs * sizeof (VBlock *), "VBlockPool"); 
        memset (&pool->vb[pool->num_vbs], 0, (num_vbs - pool->num_vbs) * sizeof (VBlockP)); // initialize new entries
        pool->num_vbs = num_vbs; 
    }
}

VBlockPool *vb_get_pool (void)
{
    return pool;
}

VBlockP vb_initialize_nonpool_vb (int vb_id, DataType dt, const char *task)
{
    uint64_t sizeof_vb = (dt != DT_NONE && dt_props[dt].sizeof_vb) ? dt_props[dt].sizeof_vb(dt) : sizeof (VBlock);

    VBlockP vb            = CALLOC (sizeof_vb);
    vb->data_type         = DT_NONE;
    vb->id                = vb_id;
    vb->compute_task      = task;
    vb->data_type         = dt;
    vb->in_use            = true;
    vb->data_type_alloced = dt;
    memset (vb->dict_id_to_did_i_map, 0xff, sizeof(vb->dict_id_to_did_i_map)); // DID_I_NONE
    return vb;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlock *vb_get_vb (const char *task_name, uint32_t vblock_i)
{
    ASSERTMAINTHREAD;
    START_TIMER;

    DataType dt = command==ZIP ? (txt_file ? txt_file->data_type : DT_NONE)
                               : (z_file   ? z_file->data_type   : DT_NONE);
    
    uint64_t sizeof_vb = (dt != DT_NONE && dt_props[dt].sizeof_vb) ? dt_props[dt].sizeof_vb(dt) : sizeof (VBlock);

    DataType alloc_dt = sizeof_vb == sizeof (VBlock)   ? DT_NONE
                      : dt == DT_REF && command == PIZ ? DT_NONE
                      : dt == DT_BAM                   ? DT_SAM
                      : dt == DT_BCF                   ? DT_VCF 
                      :                                  dt;

    // circle around until a VB becomes available (busy wait)
    unsigned vb_id; for (vb_id=0; ; vb_id = (vb_id+1) % pool->num_vbs) {

        // case: the VB is allocated to a different data type - change its data type
        if (pool->vb[vb_id] && z_file && pool->vb[vb_id]->data_type != dt) {

            DataType old_alloced_dt = pool->vb[vb_id]->data_type_alloced;

            // the new data type has a private section in its VB, that is different that the one of alloc_dt - realloc private section
            if (old_alloced_dt != alloc_dt && alloc_dt != DT_NONE) {
                
                // destroy private part of previous data_type
                if (old_alloced_dt != DT_NONE && dt_props[old_alloced_dt].destroy_vb) 
                    dt_props[old_alloced_dt].destroy_vb (pool->vb[vb_id]);    

                // calloc private part of new data_type
                VBlockP old_vb = pool->vb[vb_id];
                REALLOC (&pool->vb[vb_id], sizeof_vb, "VBlock");
                VBlockP new_vb = pool->vb[vb_id]; // new address

                memset ((char*)new_vb + sizeof (VBlock), 0, sizeof_vb - sizeof (VBlock));

                // update buf->vb in all buffers of this VB to new VB address
                buf_update_buf_list_vb_addr_change (new_vb, old_vb);
                new_vb->data_type_alloced = alloc_dt;
            }

            pool->vb[vb_id]->data_type = dt;
        }
        
        if (!pool->vb[vb_id]) { // VB is not allocated - allocate it
            pool->vb[vb_id] = CALLOC (sizeof_vb); 
            pool->num_allocated_vbs++;
            pool->vb[vb_id]->data_type = dt;
            pool->vb[vb_id]->data_type_alloced = alloc_dt;
        }

        bool in_use = __atomic_load_n (&pool->vb[vb_id]->in_use, __ATOMIC_RELAXED);
        if (!in_use) break;

        // case: we've checked all the VBs and none is available - wait a bit and check again
        // this happens when a lot VBs are handed over to the writer thread which has not processed them yet.
        // for example, if writer is blocking on write(), waiting for a pipe counterpart to read.
        // it will be released when the writer thread completes one VB.
        if (vb_id == pool->num_vbs-1) usleep (1000); // 1 ms
    }

    // initialize VB fields that need to be a value other than 0
    VBlock *vb            = pool->vb[vb_id];
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

    COPY_TIMER_VB (evb, vb_get_vb);

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

    if (command == ZIP && (flag.reference == REF_INTERNAL || flag.reference == REF_EXT_STORE))
        ref_unload_reference (gref);
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
             (VB)->vblock_i, dt_name (txt_file->data_type), (VB)->vb_position_txt_file);
    return s;
}

unsigned def_vb_size (DataType dt) { return sizeof (VBlock); }
