// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

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
#include "writer.h"

// pool of VBs allocated based on number of threads
static VBlockPool *pools[NUM_POOL_TYPES] = {};

VBlockP evb = NULL; // outside a pool

rom pool_name (VBlockPoolType pool_type)
{
    static const rom pool_names[] = POOL_NAMES;
    if (pool_type < -1 || pool_type >= NUM_POOL_TYPES) 
        return "INVALID_POOL";
    else if (pool_type == -1)
        return "NO_POOL";
    else
        return pool_names[pool_type];
}

VBlockPool *vb_get_pool (VBlockPoolType type, FailType soft_fail)
{
    ASSERT (pools[type] || soft_fail, "VB Pool %s is not allocated", pool_name (type));
    return pools[type];
}

VBlockP vb_get_from_pool (VBlockPoolP pool, int32_t vb_id) 
{
    ASSERTNOTNULL (pool);
    ASSERT (vb_id == VB_ID_EVB || vb_id < pool->num_vbs, "vb_id=%d out of range for pool %s [0,%d]", vb_id, pool_name (&pool - pools), pool->num_vbs-1);

    return (vb_id == VB_ID_EVB) ? evb : pool->vb[vb_id];
}

static VBlockP nonpool_vbs[NUM_NONPOOL_VBs] = {}; 

#define FINALIZE_VB_BUFS(func, ctx_func, vb_func, all) \
    func (vb->lines);               \
    func (vb->ra_buf[0]);           \
    func (vb->ra_buf[1]);           \
    func (vb->scratch);             \
    func (vb->z_data);              \
    func (vb->z_data_test);         \
    func (vb->z_section_headers);   \
    func (vb->spiced_pw);           \
    func (vb->show_headers_buf);    \
    func (vb->show_b250_buf);       \
    func (vb->section_list_buf);    \
    func (vb->bgzf_blocks);         \
    func (vb->coverage);            \
    func (vb->read_count);          \
    func (vb->unmapped_read_count); \
    func (vb->ol_chrom2ref_map);    \
    func (vb->gencomp_lines);       \
    func (vb->frozen_state);        \
    func (vb->txt_data);            \
    for (unsigned i=0; i < NUM_CODEC_BUFS; i++) func (vb->codec_bufs[i]); \
    if (vb->data_type != DT_NONE)   \
        for (unsigned i=0; i < MAX_(vb->num_contexts, DTF(num_fields)); i++) \
            if (all || CTX(i)->dict_id.num || i < DTF(num_fields)) \
                ctx_func (CTX(i), i); /* note: always erase num_fields as they may be set in *_seg_initialize even if not used */\
    if (vb->data_type_alloced != DT_NONE && dt_props[vb->data_type_alloced].vb_func) dt_props[vb->data_type_alloced].vb_func(vb);    

// test that vb_release_vb_do indeed frees everything
static void vb_release_verify_completeness (VBlockP vb, uint32_t vblock_i, uint32_t num_contexts)
{
    unsigned sizeof_vb = MAX_(DT_FUNC(vb, sizeof_vb)(vb->data_type_alloced), sizeof (VBlock));
    DataType dt = vb->data_type_alloced;

    VBlockP save = MALLOC (sizeof_vb);
    memcpy (save, vb, sizeof_vb);

    // erase all buffers
    vb->num_contexts = MAX_DICTS; 
    FINALIZE_VB_BUFS (buf_erase, ctx_erase_context_bufs, erase_vb_bufs, true);
    memset (&vb->buffer_list, 0, sizeof (Buffer));

    vb->num_contexts = 0;

    // clear the fields purposely not cleared by vb_release_vb_do
    vb->data_type = vb->data_type_alloced = 0; 
    vb->id = 0;
    vb->pool = 0;

    for (char *c=(char*)vb; c < (char*)vb + sizeof_vb; c++)
        if (*c && c != (char *)&vb->in_use) {  // we can't zero in_use yet, bc another thread may grab the vb
            #define field d2d_map // <-- to find the offending field, plug in field names and run iteratively
            #define REL_LOC (((char*)(&((VBlockP)vb)->field) - (char*)vb)) // <-- to use private datatype VBlocks, temporarily include the private .h for debugging
            fprintf (stderr, "sizeof_vb=%u sizeof(VBlock)=%u. Bad byte=%u Your field: %u\n", 
                        sizeof_vb, (int)sizeof (VBlock), (unsigned)(c - (char*)vb), (int)REL_LOC);
            
            ABORT ("vb_release_vb_do of %s vb_i=%u didn't fully clear the VB, byte %u != 0", dt_name(dt), vblock_i, (unsigned)(c - (char*)vb));
        }

    memcpy (vb, save, sizeof_vb);
    FREE (save);
}

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
void vb_release_vb_do (VBlockP *vb_p, rom task_name, rom func) 
{
    START_TIMER;

    VBlockP vb = *vb_p;

    if (!vb) return; // nothing to release

    ASSERT (vb->in_use || vb==evb, "Cannot release VB because it is not in_use (called from %s): vb->id=%d vb->vblock_id=%u", 
            func, vb->id, vb->vblock_i);

    threads_log_by_vb (vb, vb->compute_task ? vb->compute_task : func, "RELEASING VB", 0);

    if (flag.show_time) 
        profiler_add (vb);

    if (vb->id != VB_ID_EVB) // cannot test evb, see comment in buf_test_overflows_do
        buf_test_overflows(vb, func); 

    // verify that gzip_compressor was released after use
    ASSERT (!vb->gzip_compressor, "vb=%s: expecting gzip_compressor=NULL", VB_NAME);

    // note: we need to free the contexts before resetting vb->num_contexts
    uint32_t save_num_contexts = vb->num_contexts;
    VBIType save_vb_i = save_vb_i;
    FINALIZE_VB_BUFS (buf_free, ctx_free_context, release_vb, false);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->data_type   : type of this vb (maybe changed by vb_get_vb)
    // vb->data_type_alloced (maybe changed by vb_get_vb)
    // vb->id
    // vb->pool

    vb->vblock_i = vb->fragment_len = vb->fragment_num_words = 0;
    vb->line_i = 0;
    vb->recon_size = vb->txt_size = vb->longest_line_len = vb->sample_i = 0;
    vb->comp_i = 0;
    vb->dispatch = vb->is_processed = vb->preprocessing = vb->has_ctx_index = vb->show_containers = vb->is_eof = 
    vb->debug_line_hash_skip = false;
    vb->z_next_header_i = 0;
    vb->num_contexts = vb->curr_field = 0;
    vb->chrom_node_index = vb->chrom_name_len = vb->seq_len = 0; 
    vb->vb_position_txt_file = vb->line_start = 0;
    vb->num_lines_at_1_3 = vb->num_lines_at_2_3 = vb->num_nondrop_lines = vb->debug_line_hash = 0;
    vb->num_type1_subfields = vb->num_type2_subfields = 0;
    vb->range = NULL;
    vb->flags = (union FlagsVbHeader){};
    vb->drop_curr_line = vb->chrom_name = vb->fragment_start = NULL;
    vb->prev_range[0] = vb->prev_range[1] = NULL;
    vb->prev_range_chrom_node_index[0] = vb->prev_range_chrom_node_index[1] = vb->prev_range_range_i = 0;
    vb->digest = DIGEST_NONE;
    vb->expected_digest = DIGEST_NONE;
    vb->refhash_layer = vb->refhash_start_in_layer = 0;
    vb->fragment_ctx = vb->ht_matrix_ctx = vb->runs_ctx = vb->fgrc_ctx = NULL;
    vb->codec_using_codec_bufs = 0;
    vb->ht_per_line = 0;
    vb->compute_thread_id = 0;
    vb->compute_task = NULL;
    vb->compute_func = NULL;
    vb->translation = (DtTranslation){};
    vb->ref = NULL;
    vb->rback_id = 0;
    vb->is_dropped = 0;
    vb->vb_bgzf_i = 0;
    vb->line_bgzf_uoffset = 0;
    vb->start_compute_timestamp = 0;

    memset(&vb->profile,  0, sizeof (vb->profile));
    memset(vb->d2d_map,   0, sizeof(vb->d2d_map));
    memset(vb->ctx_index, 0, sizeof(vb->ctx_index));
    vb->iupacs_last_range[0] = vb->iupacs_last_range[1] = NULL;
    vb->iupacs_last_pos[0] = vb->iupacs_last_pos[1] = vb->iupacs_next_pos[0] = vb->iupacs_next_pos[1] = 0;
    vb->num_rollback_ctxs = 0;
    memset (vb->rollback_ctxs, 0, sizeof(vb->rollback_ctxs));
    mutex_destroy (vb->ready_for_compute);

    if (flag.debug && vb->in_use) 
        vb_release_verify_completeness (vb, save_vb_i, save_num_contexts);

    // this release can be run by either the main or writer thread. we make sure to update in_use as the very
    // last change, and do so atomically

    // case: this VB is from the pool (i.e. not evb)
    int32_t num_in_use = -1;
    if (vb != evb) {
        // Logic: num_in_use is always AT LEAST sum(vb)->in_use. i.e. pessimistic. (it can be mometarily less between these too updates)
        __atomic_store_n (&vb->in_use, (bool)0, __ATOMIC_SEQ_CST);   // released the VB back into the pool - it may now be reused 
        if (vb->id >= 0) 
            num_in_use = __atomic_sub_fetch (&pools[vb->pool]->num_in_use, 1, __ATOMIC_SEQ_CST); // atomic to prevent concurrent update by writer thread and main thread (must be after update of in_use)
        *vb_p = NULL;
    }

    if (flag_is_show_vblocks (task_name)) 
        iprintf (vb->id >= 0 ? "VB_RELEASE(task=%s id=%d) vb=%s caller=%s in_use=%d/%d\n"
                             : "VB_RELEASE(task=%s id=%d) vb=%s caller=%s\n", 
                 task_name, vb->id, VB_NAME, func, num_in_use, (vb->id >= 0 ? pools[vb->pool]->num_vbs : -1));

    buf_compact_buf_list (vb);

    if (vb->id >= 0) COPY_TIMER_EVB (vb_release_vb_do);
}

void vb_destroy_vb_do (VBlockP *vb_p, rom func)
{
    ASSERTMAINTHREAD;
    START_TIMER;

    VBlockP vb = *vb_p;
    if (!vb) return;

    if (flag_is_show_vblocks (vb->compute_task)) 
        iprintf ("VB_DESTROY(id=%d) vb_i=%d caller=%s\n", vb->id, vb->vblock_i, func);

    bool is_evb = vb->id == VB_ID_EVB;

    buf_destroy_vb_bufs (vb);

    // case: this is a nonpool VB
    for (int i=0; i < NUM_NONPOOL_VBs; i++)
        if (nonpool_vbs[i] == vb) 
            nonpool_vbs[i] = NULL;

    FREE (*vb_p);

    if (!is_evb) COPY_TIMER_EVB (vb_destroy_vb);
}

void vb_create_pool (VBlockPoolType type)
{
    // only main-thread dispatcher can create a pool. other dispatcher (eg writer's bgzf compression) can must existing pool
    uint32_t num_vbs = (type == POOL_MAIN) ? MAX_(1, global_max_threads) +               // compute thread VBs
                                             (IS_PIZ ? 2 : 0)    +                       // txt header VB and wvb (for PIZ)
                                             (IS_PIZ ? z_file->max_conc_writing_vbs : 0) // thread-less VBs handed over to the writer thread
                                           : writer_get_max_bgzf_threads();
    if (flag_is_show_vblocks (NULL)) 
        iprintf ("CREATING_VB_POOL: type=%s global_max_threads=%u max_conc_writing_vbs=%u num_vbs=%u\n", 
                 pool_name (type), global_max_threads, z_file->max_conc_writing_vbs, num_vbs); 

    uint32_t size = sizeof (VBlockPool) + num_vbs * sizeof (VBlockP);

    if (!pools[type])  
        // allocation includes array of pointers (initialized to NULL)
        pools[type] = (VBlockPool *)CALLOC (size); // note we can't use Buffer yet, because we don't have VBs yet...

    // case: old pool is too small - realloc it (the pool contains only pointers to VBs, so the VBs themselves are not realloced)
    else if (pools[type]->num_vbs < num_vbs) {
        REALLOC (&pools[type], size, "VBlockPool"); 
        memset (&pools[type]->vb[pools[type]->num_vbs], 0, (num_vbs - pools[type]->num_vbs) * sizeof (VBlockP)); // initialize new entries
    }

    pools[type]->size    = size; 
    pools[type]->num_vbs = num_vbs; 
}

VBlockP vb_initialize_nonpool_vb (int vb_id, DataType dt, rom task)
{
    uint64_t sizeof_vb = (dt != DT_NONE && dt_props[dt].sizeof_vb) ? dt_props[dt].sizeof_vb(dt) : sizeof (VBlock);

    VBlockP vb            = CALLOC (sizeof_vb);
    vb->data_type         = DT_NONE;
    vb->id                = vb_id;
    vb->compute_task      = task;
    vb->data_type         = dt;
    vb->in_use            = true;
    vb->data_type_alloced = dt;
    vb->comp_i            = COMP_NONE;
    vb->pool              = NO_POOL;
    init_dict_id_to_did_map (vb->d2d_map); 
    
    nonpool_vbs[NUM_NONPOOL_VBs + vb_id] = vb; // vb_id is a negative integer

    return vb;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlockP vb_get_vb (VBlockPoolType type, rom task_name, VBIType vblock_i, CompIType comp_i)
{
    START_TIMER;
    
    VBlockPool *pool = vb_get_pool (type, HARD_FAIL);

    DataType dt = (type == POOL_BGZF)            ? DT_NONE
                : (flag.deep && flag.zip_comp_i) ? DT_FASTQ
                : (IS_ZIP && segconf.has_embdedded_fasta) ? DT_FASTA // GFF3 with embedded FASTA
                : (IS_ZIP && txt_file)           ? txt_file->data_type
                : (IS_PIZ && z_file && flag.deep && comp_i >= SAM_COMP_FQ00) ? DT_FASTQ
                : (IS_PIZ && z_file)             ? z_file->data_type  
                :                                  DT_NONE;
    
    uint64_t sizeof_vb = (dt != DT_NONE && dt_props[dt].sizeof_vb) ? dt_props[dt].sizeof_vb(dt) : sizeof (VBlock);

    DataType alloc_dt = sizeof_vb == sizeof (VBlock) ? DT_NONE
                      : dt == DT_REF && IS_PIZ       ? DT_NONE
                      : dt == DT_BAM                 ? DT_SAM
                      : dt == DT_BCF                 ? DT_VCF 
                      :                                dt;

    // circle around until a VB becomes available (busy wait)
    uint32_t vb_id; for (vb_id=0; ; vb_id = (vb_id+1) % pool->num_vbs) {

        // case: the VB is allocated to a different data type - change its data type
        if (pool->vb[vb_id] && z_file && pool->vb[vb_id]->data_type != dt && !pool->vb[vb_id]->in_use) {

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
            pool->vb[vb_id]->data_type_alloced = alloc_dt;
        }

        bool in_use = __atomic_load_n (&pool->vb[vb_id]->in_use, __ATOMIC_SEQ_CST);

        if (!in_use) break;

        // case: we've checked all the VBs and none is available - wait a bit and check again
        // this happens when a lot VBs are handed over to the writer thread which has not processed them yet.
        // for example, if writer is blocking on write(), waiting for a pipe counterpart to read.
        // it will be released when the writer thread completes one VB.
        if (vb_id == pool->num_vbs-1) usleep (10000); // 10 ms
    }

    VBlockP vb            = pool->vb[vb_id];

    // Logic: num_in_use is always AT LEAST sum(vb)->in_use. i.e. pessimistic. (it can be mometarily more between these two updates)
    uint32_t num_in_use = __atomic_add_fetch (&pool->num_in_use, 1, __ATOMIC_SEQ_CST); // atomic to prevent concurrent update by writer thread and main thread (must be before update of in_use)
    __atomic_store_n (&vb->in_use, (bool)1, __ATOMIC_SEQ_CST);   

    // initialize VB fields that need to be a value other than 0
    vb->id                = vb_id;
    vb->data_type         = dt;
    vb->vblock_i          = vblock_i;
    vb->comp_i            = comp_i;
    vb->buffer_list.vb    = vb;
    vb->compute_thread_id = THREAD_ID_NONE;
    vb->compute_task      = task_name;
    vb->pool              = type;
    init_dict_id_to_did_map (vb->d2d_map);

    if (flag_is_show_vblocks (task_name)) 
        iprintf ("VB_GET_VB(task=%s id=%u) vb_i=%s/%d in_use=%u/%u\n", 
                  task_name, vb->id, comp_name (vb->comp_i), vb->vblock_i, num_in_use, pool->num_vbs);

    threads_log_by_vb (vb, task_name, "GET VB", 0);

    if (flag.debug_memory)
        iprintf ("vb_get_vb: got vb_i=%d id=%d task=%s dt=%s address=[%p - %p]\n", 
                 vblock_i, vb_id, task_name, dt_name(alloc_dt), vb, (char*)vb + sizeof_vb);

    COPY_TIMER_EVB (vb_get_vb);
    return vb;
}

uint32_t vb_pool_get_num_in_use (VBlockPoolType type, int32_t *id/*optional out*/)
{
    VBlockPool *pool = vb_get_pool (type, HARD_FAIL);
    int num_in_use = __atomic_load_n (&pool->num_in_use, __ATOMIC_SEQ_CST); // atomic, be for POOL_MAIN, writer thread might update concurrently.

    if (id) {
        *id = -999;
        for (uint32_t vb_id=0; vb_id < pool->num_vbs; vb_id++)
            if (pool->vb[vb_id] && pool->vb[vb_id]->in_use) return *id = vb_id;
    }

    return num_in_use;
}

// Note: num_in_use is always AT LEAST sum(vb)->in_use (it can be mometarily less than sum(vb)->in_use as they are getting updated)
// therefore, this function may return true when pool is actually no longer full.
bool vb_pool_is_full (VBlockPoolType type)
{
    return vb_pool_get_num_in_use (type, NULL) == vb_get_pool(type, HARD_FAIL)->num_vbs;
}
 
// Note: As in vb_pool_is_full, if the function returns false (not empty) the pool might in fact already be empty
bool vb_pool_is_empty (VBlockPoolType type)
{
    return vb_pool_get_num_in_use (type, NULL) == 0;
}

bool vb_is_valid (VBlockP vb)
{
    if (vb->pool != NO_POOL)
        for (uint32_t vb_id=0; vb_id < pools[vb->pool]->num_vbs; vb_id++)
            if (vb == pools[vb->pool]->vb[vb_id]) return true;

    for (int i=0; i < NUM_NONPOOL_VBs; i++)
        if (nonpool_vbs[i] == vb) return true;

    return false;
}

// frees memory of all VBs, except for non-pool VBs (evb, segconf...)
void vb_destroy_pool_vbs (VBlockPoolType type)
{
    if (!pools[type]) return;

    for (uint32_t vb_id=0; vb_id < pools[type]->num_vbs; vb_id++) 
        vb_destroy_vb (&pools[type]->vb[vb_id]);

    FREE (pools[type]);
}

// NOT thread safe, use only in execution-terminating messages
rom err_vb_pos (void *vb)
{
    static char s[80];
    sprintf (s, "vb i=%u position in %s file=%"PRIu64, 
             (VB)->vblock_i, dt_name (txt_file->data_type), (VB)->vb_position_txt_file);
    return s;
}

unsigned def_vb_size (DataType dt) { return sizeof (VBlock); }

bool vb_is_processed (VBlockP vb)
{
    return __atomic_load_n (&vb->is_processed, __ATOMIC_RELAXED);
}

void vb_set_is_processed (VBlockP vb)
{
    __atomic_store_n (&vb->is_processed, (bool)true, __ATOMIC_RELAXED); 
}

bool vb_buf_locate (VBlockP vb, ConstBufferP buf)
{
    if (!vb) return false;

    unsigned sizeof_vb = MAX_(DT_FUNC(vb, sizeof_vb)(vb->data_type_alloced), sizeof (VBlock));

    return vb && (char*)buf >= (char*)vb && 
                 (char*)buf <  (char*)vb + sizeof_vb;
}
